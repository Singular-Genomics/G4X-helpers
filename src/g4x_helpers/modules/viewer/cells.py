import logging

import numpy as np
import pandas as pd
import polars as pl
import zarr
from numcodecs import Blosc
from scipy.sparse import csr_matrix
from shapely import to_ragged_array

from ... import c, io
from . import zarr_utils as utils

log = logging.getLogger(__name__)
UNASSIGNED_CELL = 'unassigned'

COMPRESSOR = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)


def write_cells(
    zarr_path: str,
    segmentation_mask: np.ndarray,
    cell_metadata: pl.DataFrame,
    cell_x_gene: pl.DataFrame,
    clustering_umap: pl.DataFrame,
    cell_x_protein: pl.DataFrame | None = None,
    seg_name: str = 'g4x-default',
    overwrite: bool = True,
    show_progress: bool = False,
    **kwargs,
):

    zarr_path = io.pathval.validate_dir_path(zarr_path)

    cell_group = zarr.open_group(zarr_path / 'cells', mode='r+')

    log.debug('Setting up cell data group')
    seg_path = _add_segmentation_attrs(cell_group, seg_name)
    seg_group = cell_group.create_group(seg_path, overwrite=overwrite)

    polygons = extract_polygons(segmentation_mask, show_progress=show_progress, **kwargs)
    cell_x_gene = cell_x_gene.sort(c.CELL_ID_NAME)
    cell_metadata = cell_metadata.sort(c.CELL_ID_NAME)
    clustering_umap = clustering_umap.sort(c.CELL_ID_NAME)

    if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_gene[c.CELL_ID_NAME]):
        raise ValueError('Cell IDs in metadata and gene expression matrix do not match')

    if not all(cell_metadata[c.CELL_ID_NAME].to_numpy() == polygons[c.CELL_ID_NAME].array):
        raise ValueError('Cell IDs in metadata and polygons do not match')

    if cell_x_protein is not None:
        cell_x_protein = cell_x_protein.sort(c.CELL_ID_NAME)
        if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_protein[c.CELL_ID_NAME]):
            raise ValueError('Cell IDs in metadata and protein table do not match')

        cell_metadata = cell_metadata.join(cell_x_protein, on=c.CELL_ID_NAME)

    gex, gene_names = _cellxgene_to_csr(cell_x_gene)

    cell_metadata = cell_metadata.cast({c.CELL_ID_NAME: pl.UInt32})
    clustering_umap = clustering_umap.cast({c.CELL_ID_NAME: pl.UInt32})
    metadata = cell_metadata.join(clustering_umap, on=c.CELL_ID_NAME, how='left')
    metadata = metadata.with_columns(pl.col('^leiden.*$').fill_null(UNASSIGNED_CELL))

    ################################
    clusterings = [c for c in metadata.columns if c.startswith('leiden')]
    clusterings_order = _get_sorted_clusterings(metadata, clusterings)

    protein_columns = [col for col in metadata.columns if c.IMG_INTENSITY_HANDLE in col]
    protein_names = [s.removesuffix(c.IMG_INTENSITY_HANDLE) for s in protein_columns]

    cluster_labels_meta = {}
    for i, key in enumerate(clusterings_order):
        sorted_cluster_ids = _get_sorted_cluster_ids(metadata, cluster_key=key)
        cluster_color_map = _generate_cluster_palette(sorted_cluster_ids)

        cluster_labels_meta[key] = {
            'index': i,
            'clusterID_order': list(cluster_color_map.keys()),
            'clusterID_colors': cluster_color_map,
        }

    seg_group.attrs['cluster_labels'] = cluster_labels_meta
    seg_group.attrs['cluster_labels_order'] = clusterings_order
    seg_group.attrs['genes_shape'] = gex.shape

    gene_name_array = np.array(gene_names).astype('U')
    prot_name_array = np.array(protein_names).astype('U')
    utils.create_array(
        seg_group,
        'gene_names',
        data=gene_name_array,
        compressor=COMPRESSOR,
        chunks=gene_name_array.shape,
    )
    utils.create_array(
        seg_group,
        'protein_names',
        data=prot_name_array,
        compressor=COMPRESSOR,
        chunks=prot_name_array.shape,
    )

    for column in ['total_counts', 'n_genes_by_counts']:
        if column not in metadata.columns:
            metadata = metadata.with_columns(pl.lit(None).cast(pl.UInt16).alias(column))

    ################## preparing to write arrays
    meta_columns = {
        'cell_id': (metadata[c.CELL_ID_NAME], 'uint32'),
        'area': (metadata[c.CELL_AREA_NAME], 'uint16'),
        'position': (metadata.select([c.CELL_COORD_X, c.CELL_COORD_Y]), 'float16'),
        'cluster_id': (metadata.select(clusterings_order), 'U'),
        'total_counts': (metadata['total_counts'], 'uint16'),
        'total_genes': (metadata['n_genes_by_counts'], 'uint16'),
        'protein_values': (metadata.select(protein_columns).fill_null(np.nan), 'float16'),
        'umap': (metadata.select(['UMAP1', 'UMAP2']).fill_null(np.nan), 'float16'),
        'gene_counts': (gex.data, 'uint16'),
        'gene_indices': (gex.indices, 'int32'),
        'gene_indptr': (gex.indptr, 'int32'),
    }

    log.info('Writing cell metadata arrays')
    _write_metadata_arrays(seg_group, meta_columns)

    log.info('Writing cell polygon arrays')
    ragged = to_ragged_array(polygons.geometry, include_z=False, include_m=False)
    verts_xy, offsets = ragged[1], ragged[2][0]
    for array, name in zip([offsets, verts_xy], ['polygon_offsets', 'polygon_vertices_xy']):
        chunks = utils.calculate_chunks(array, target_mb=4)
        utils.create_array(seg_group, name, data=array, compressor=COMPRESSOR, chunks=chunks)


def _write_metadata_arrays(seg_group, meta_columns):
    for key, (arr, dtype) in meta_columns.items():
        if key in seg_group:
            del seg_group[key]

        array = arr.to_numpy() if isinstance(arr, (pl.DataFrame, pl.Series)) else arr
        array = array.astype(dtype)
        chunks = utils.calculate_chunks(array, target_mb=4)

        utils.create_array(seg_group, key, data=array, compressor=COMPRESSOR, chunks=chunks)


def _get_sorted_clusterings(df, cluster_keys: list[str]):
    clusterings_order = (
        df.select(cluster_keys)
        .unpivot(cluster_keys)
        .group_by('variable')
        .agg(pl.col('value').n_unique())
        .sort('value', descending=False)['variable']
        .to_list()
    )
    return clusterings_order


def _get_sorted_cluster_ids(df, cluster_key: str):
    cluster_ids_order = (
        df.select(cluster_key).group_by(cluster_key).agg(pl.len()).sort('len', descending=True)[cluster_key].to_list()
    )

    if UNASSIGNED_CELL in cluster_ids_order:
        cluster_ids_order.remove(UNASSIGNED_CELL)
        cluster_ids_order.append(UNASSIGNED_CELL)

    return cluster_ids_order


def _generate_cluster_palette(ordered_unique_clusters: list, max_colors: int = 256) -> dict:

    n_clusters = len(ordered_unique_clusters)

    if n_clusters <= 20:
        hex_list = c.SG_PALETTE

    else:
        from matplotlib.pyplot import get_cmap

        hex_list = get_cmap('hsv', min(max_colors, n_clusters)).colors

    cluster_palette = {}
    for i, cluster in enumerate(ordered_unique_clusters):
        cluster_palette[str(cluster)] = utils.hex_to_rgb(hex_list[i])

    cluster_palette[UNASSIGNED_CELL] = utils.hex_to_rgb(c.UNASSIGNED_COLOR)

    return cluster_palette


def _cellxgene_to_csr(cell_x_gene: pl.DataFrame) -> tuple[csr_matrix, np.ndarray]:

    cell_x_gene = cell_x_gene.drop(c.CELL_ID_NAME)
    gene_names = np.array(cell_x_gene.columns)

    gex = cell_x_gene.to_numpy().astype(np.uint16)
    gex = csr_matrix(gex)
    return gex, gene_names


def extract_polygons(mask, buffer: float = 4.0, simplify_tolerance: float = 1.0, show_progress: bool = False):
    polygons = io.convert.ndarray_to_gdf(mask, show_progress=show_progress)

    # Keep only largest polygon per label
    polygons['_area'] = polygons.geometry.area
    polygons = (
        polygons.sort_values('_area', ascending=False)
        .drop_duplicates(subset=c.CELL_ID_NAME, keep='first')
        .drop(columns='_area')
        .reset_index(drop=True)
    )

    # Alternative way to keep largest polygon per label
    # idx = gdf.groupby(CELL_ID_NAME)["_area"].idxmax()
    # gdf = gdf.loc[idx].drop(columns="_area").reset_index(drop=True)

    # Simplify geometries
    polygons['geometry_simplified'] = polygons.geometry.buffer(buffer).buffer(-buffer)
    if simplify_tolerance > 0:
        polygons['geometry_simplified'] = polygons.geometry_simplified.simplify(
            tolerance=simplify_tolerance, preserve_topology=True
        )

    polygons['geometry'] = polygons.geometry_simplified
    polygons = polygons.drop(columns='geometry_simplified')

    return polygons.sort_values(c.CELL_ID_NAME)


def _sanitize_path_component(s, replacement='_'):
    invalid = r'<>:"/\\|?*- '
    for ch in invalid:
        s = s.replace(ch, replacement)
    return s.strip(' .')  # Windows disallows trailing space/dot


def _add_segmentation_attrs(cell_group, seg_name):
    seg_path = _sanitize_path_component(seg_name) + '_segmentation'

    seg_sources = cell_group.attrs['segmentation_sources']
    seg_order = cell_group.attrs['segmentation_order']

    seg_sources.update({seg_name: seg_path})
    seg_order.append(seg_name)

    cell_group.attrs['segmentation_sources'] = seg_sources
    cell_group.attrs['segmentation_order'] = list(set(seg_order))
    return seg_path


# NOTE unused
def _prepare_metadata_for_tiling(metadata, tile_size, img_res):
    metadata = metadata.with_row_index()

    image_resolution_hw = img_res
    n_tiles_w = image_resolution_hw[1] // tile_size
    n_tiles_h = image_resolution_hw[0] // tile_size
    n_tiles_w, n_tiles_h

    metadata = metadata.with_columns(
        (pl.col('cell_x') / tile_size).cast(pl.Int32).alias('tile_x'),
        (pl.col('cell_y') / tile_size).cast(pl.Int32).alias('tile_y'),
    ).sort('tile_y', 'tile_x')

    return metadata


def _map_categories(mask: np.ndarray, labels: np.ndarray, categories: np.ndarray, missing_val=-1):
    flat_mask = mask.ravel()

    idx = pd.Index(labels)
    pos = idx.get_indexer(flat_mask)  # -1 where not found

    out_flat = np.full(flat_mask.shape, missing_val, dtype=int)

    valid = pos != -1
    if valid.any():
        out_flat[valid] = categories[pos[valid]]

    return out_flat.reshape(mask.shape)


def _map_clusters_to_mask(meta: pl.DataFrame, cluster_key: str, mask: np.ndarray):
    cluster_cat = pd.Categorical(meta[cluster_key])

    if 'unassigned' not in cluster_cat.categories:
        cluster_cat = cluster_cat.add_categories('unassigned')

    cats = list(cluster_cat.categories)
    cats = ['unassigned'] + [c for c in cats if c != 'unassigned']
    cluster_cat = cluster_cat.reorder_categories(cats)

    cluster_codes = cluster_cat.codes  # integers 0..n-1
    p_mask = _map_categories(mask=mask, labels=meta['cell_id'].to_numpy(), categories=cluster_codes)

    pal = [c.UNASSIGNED_COLOR] + c.SG_PALETTE
    pal = np.array([utils.hex_to_rgb(c, normalized=False) for c in pal])

    rgb = pal[p_mask]
    rgb[p_mask == -1] = [0, 0, 0]

    return rgb


def get_user_cluster_metadata(new_data: pl.DataFrame) -> dict[str, dict]:
    cols = new_data.columns
    cluster_cols = [col for col in cols if col not in (c.CELL_ID_NAME, 'UMAP1', 'UMAP2')]

    clusterings_dict = {
        col: f'{col}_color' if f'{col}_color' in cluster_cols else None
        for col in cluster_cols
        if not col.endswith('_color')
    }

    clusterings_order = _get_sorted_clusterings(new_data, clusterings_dict.keys())

    cluster_labels_meta = {}
    for i, key in enumerate(clusterings_order):
        sorted_cluster_ids = _get_sorted_cluster_ids(new_data, key)
        if clusterings_dict[key] is None:
            cluster_color_map = _generate_cluster_palette(sorted_cluster_ids)
        else:
            color_map = (
                new_data.select([key, clusterings_dict[key]])
                .unique()
                .cast({key: pl.Enum(sorted_cluster_ids)})
                .sort(key)
            )
            cluster_color_map = {}
            for row in color_map.iter_rows():
                cluster_color_map[row[0]] = utils.hex_to_rgb(row[1])

        cluster_labels_meta[key] = {
            'index': i,
            'clusterID_order': list(cluster_color_map.keys()),
            'clusterID_colors': cluster_color_map,
        }

    return cluster_labels_meta, clusterings_order


def get_seg_group(viewer_dir, segmentation_source='g4x_default_segmentation'):
    cell_group = zarr.open(viewer_dir / 'cells', mode='r+')

    if segmentation_source not in cell_group:
        raise ValueError(
            f'Segmentation source "{segmentation_source}" not found in the data. Available sources: {list(cell_group.keys())}'
        )

    return cell_group[segmentation_source]


def get_cell_metadata(seg_group):
    cell_ids = pl.Series(name=c.CELL_ID_NAME, values=seg_group['cell_id'][:])

    seg_meta_dict = dict(seg_group.attrs)
    df = pl.DataFrame(cell_ids)

    umap_dims = ['UMAP1', 'UMAP2']
    for i, label in enumerate(umap_dims):
        s = pl.Series(name=label, values=seg_group['umap'][:, i])
        df = df.with_columns(s)

    cluster_names = seg_meta_dict['cluster_labels_order']
    for i, label in enumerate(cluster_names):
        s = pl.Series(name=label, values=seg_group['cluster_id'][:, i])
        df = df.with_columns(s)

        cmap = seg_meta_dict['cluster_labels'][label]['clusterID_colors']
        cmap = {k: utils.rgb_to_hex(v) for k, v in cmap.items()}
        df = df.with_columns(pl.col(label).replace(cmap).alias(f'{label}_color'))

    return df


def apply_viewer_metadata(seg_group, new_data):
    new_data_df = pl.read_csv(new_data)
    existing_data = get_cell_metadata(seg_group)

    if c.CELL_ID_NAME not in new_data_df.columns:
        raise ValueError(f"Expected column '{c.CELL_ID_NAME}' not found in new data")

    if not new_data_df[c.CELL_ID_NAME].equals(existing_data[c.CELL_ID_NAME]):
        raise ValueError('Cell ID columns do not match between existing and new data')

    cluster_labels_meta, clusterings_order = get_user_cluster_metadata(new_data_df)
    seg_group.attrs['cluster_labels'] = cluster_labels_meta
    seg_group.attrs['cluster_labels_order'] = clusterings_order

    meta_columns = {
        'cluster_id': (new_data_df.select(clusterings_order), 'U'),
        'umap': (new_data_df.select(['UMAP1', 'UMAP2']).fill_null(np.nan), 'float16'),
    }

    _write_metadata_arrays(seg_group, meta_columns)
