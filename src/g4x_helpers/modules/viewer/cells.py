import logging
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
import zarr
from numcodecs import Blosc
from shapely import to_ragged_array

from ... import c, io
from ...schema.definition import CellMetadata, CellxGene, CellxProt, ClusteringUmap, Segmentation
from ..workflow import PRESET_SOURCE, collect_input
from .utils import calculate_chunks, create_array

if TYPE_CHECKING:
    from ...g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)
UNASSIGNED_CELL = 'unassigned'

COMPRESSOR = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)


def write_cells(
    smp: 'G4Xoutput',
    seg_name: str = 'g4x-default',
    components: tuple | None = None,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    log.info('Preparing cell data')

    # mode = 'w' if overwrite else 'a'
    cell_group = zarr.open_group(smp.out.ViewerZarr.p / 'cells', mode='r+')

    log.debug('Setting up cell data group')
    seg_path = _add_segmentation_attrs(cell_group, seg_name)
    seg_group = cell_group.create_group(seg_path, overwrite=overwrite)

    if components is None:
        log.debug('No components provided, processing cell data from source')
        components = process_cell_data(smp, logger=log)
    else:
        log.debug('Using provided components to select data')

    metadata, gex, gene_names, verts_xy, offsets = components

    clusterings = [c for c in metadata.columns if c.startswith('leiden')]
    clusterings_order = get_sorted_clusterings(metadata, clusterings)

    protein_columns = [col for col in metadata.columns if c.IMG_INTENSITY_HANDLE in col]
    protein_names = [s.removesuffix(c.IMG_INTENSITY_HANDLE) for s in protein_columns]

    cluster_labels_meta = {}
    for i, key in enumerate(clusterings_order):
        sorted_cluster_ids = get_sorted_cluster_ids(metadata, cluster_key=key)
        cluster_color_map = generate_cluster_palette(sorted_cluster_ids)

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
    create_array(seg_group, 'gene_names', data=gene_name_array, compressor=COMPRESSOR, chunks=gene_name_array.shape)
    create_array(seg_group, 'protein_names', data=prot_name_array, compressor=COMPRESSOR, chunks=prot_name_array.shape)

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
    for key, (arr, dtype) in meta_columns.items():
        array = arr.to_numpy() if isinstance(arr, (pl.DataFrame, pl.Series)) else arr

        array = array.astype(dtype)
        chunks = calculate_chunks(array, target_mb=4)

        create_array(seg_group, key, data=array, compressor=COMPRESSOR, chunks=chunks)

    log.info('Writing cell polygon arrays')
    for array, name in zip([offsets, verts_xy], ['polygon_offsets', 'polygon_vertices_xy']):
        chunks = calculate_chunks(array, target_mb=4)
        create_array(seg_group, name, data=array, compressor=COMPRESSOR, chunks=chunks)


def prepare_metadata_for_tiling(metadata, tile_size, img_res):
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


def process_cell_data(
    smp,
    segmentation_mask: str = PRESET_SOURCE,
    cell_metadata: str = PRESET_SOURCE,
    cell_x_gene: str = PRESET_SOURCE,
    cell_x_protein: str = PRESET_SOURCE,
    clustering_umap: str = PRESET_SOURCE,
    logger: logging.Logger | None = None,
):
    from scipy.sparse import csr_matrix

    log = logger or LOGGER

    collect_in_partial = partial(collect_input, smp, validate=True, logger=log)
    segment_in = collect_in_partial(segmentation_mask, Segmentation, validate=False)
    cellmet_in = collect_in_partial(cell_metadata, CellMetadata)
    cellxgene_in = collect_in_partial(cell_x_gene, CellxGene)
    clustumap_in = collect_in_partial(clustering_umap, ClusteringUmap)

    cell_metadata = cellmet_in.load()

    # process cell by gene
    cell_x_gene = cellxgene_in.load()

    if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_gene[c.CELL_ID_NAME]):
        raise ValueError('Cell IDs in metadata and gene expression matrix do not match')

    cell_x_gene = cell_x_gene.drop(c.CELL_ID_NAME)
    gene_names = np.array(cell_x_gene.columns)

    gex = cell_x_gene.to_numpy().astype(np.uint16)
    gex = csr_matrix(gex)

    # in case of protein only runs, this will be the case and we need to create empty arrays for the gex to avoid issues downstream
    if len(gex.data) == 0:
        log.warning('Gene expression matrix is empty, creating empty arrays for gex')
        fill = len(gex.indptr)
        gex.data = np.zeros(fill, dtype='uint16')
        gex.indices = np.zeros(fill, dtype='int32')

    del cell_x_gene

    # process cell by protein (optional)
    if smp.src.pr_detected:
        cellxprot_in = collect_in_partial(cell_x_protein, CellxProt)
        cell_x_protein = cellxprot_in.load()

        if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_protein[c.CELL_ID_NAME]):
            raise ValueError('Cell IDs in metadata and protein table do not match')

        cell_metadata = cell_metadata.join(cell_x_protein, on=c.CELL_ID_NAME)
        del cell_x_protein

    # 2: extract cell vertices
    log.info('Extracting cell vertices from segmentation')
    gdf = extract_vertices(segment_in.load(), show_progress=False)

    assert all(cell_metadata[c.CELL_ID_NAME].to_pandas() == gdf[c.CELL_ID_NAME].array)

    ragged = to_ragged_array(gdf.geometry_simplified, include_z=False, include_m=False)
    verts_xy, offsets = ragged[1], ragged[2][0]

    clust_umap = clustumap_in.load()

    cell_metadata = cell_metadata.cast({c.CELL_ID_NAME: pl.UInt64})
    clust_umap = clust_umap.cast({c.CELL_ID_NAME: pl.UInt64})

    cell_metadata = cell_metadata.join(clust_umap, on=c.CELL_ID_NAME, how='left')
    cell_metadata = cell_metadata.with_columns(pl.col('^leiden.*$').fill_null(UNASSIGNED_CELL))

    return cell_metadata, gex, gene_names, verts_xy, offsets


def get_sorted_clusterings(df, cluster_keys: list[str]):
    clusterings_order = (
        df.select(cluster_keys)
        .unpivot(cluster_keys)
        .group_by('variable')
        .agg(pl.col('value').n_unique())
        .sort('value', descending=False)['variable']
        .to_list()
    )
    return clusterings_order


def get_sorted_cluster_ids(df, cluster_key: str):
    cluster_ids_order = (
        df.select(cluster_key).group_by(cluster_key).agg(pl.len()).sort('len', descending=True)[cluster_key].to_list()
    )

    if UNASSIGNED_CELL in cluster_ids_order:
        cluster_ids_order.remove(UNASSIGNED_CELL)
        cluster_ids_order.append(UNASSIGNED_CELL)

    return cluster_ids_order


def generate_cluster_palette(ordered_unique_clusters: list, max_colors: int = 256) -> dict:
    import matplotlib.colors as mcolors

    def hex2rgb(hex: str) -> list[int, int, int]:
        return [int(x * 255) for x in mcolors.to_rgb(hex)]

    n_clusters = len(ordered_unique_clusters)

    if n_clusters <= 20:
        hex_list = c.SG_PALETTE

    else:
        from matplotlib.pyplot import get_cmap

        hex_list = get_cmap('hsv', min(max_colors, n_clusters)).colors

    cluster_palette = {}
    for i, cluster in enumerate(ordered_unique_clusters):
        cluster_palette[str(cluster)] = hex2rgb(hex_list[i])

    cluster_palette[UNASSIGNED_CELL] = hex2rgb(c.UNASSIGNED_COLOR)

    return cluster_palette


def extract_vertices(mask, show_progress: bool = False):
    gdf = io.convert.ndarray_to_gdf(mask, show_progress=show_progress)

    # Keep only largest polygon per label
    gdf['_area'] = gdf.geometry.area
    gdf = (
        gdf.sort_values('_area', ascending=False)
        .drop_duplicates(subset=c.CELL_ID_NAME, keep='first')
        .drop(columns='_area')
        .reset_index(drop=True)
    ).sort_values(c.CELL_ID_NAME)

    # Alternative way to keep largest polygon per label
    # idx = gdf.groupby(CELL_ID_NAME)["_area"].idxmax()
    # gdf = gdf.loc[idx].drop(columns="_area").reset_index(drop=True)

    # Simplify geometries
    gdf['geometry_simplified'] = gdf.geometry.buffer(4).buffer(-4)
    gdf['geometry_simplified'] = gdf.geometry_simplified.simplify(tolerance=1, preserve_topology=True)
    # gdf['geometry_simplified'] = gdf.geometry.simplify(tolerance=1.5, preserve_topology=True)
    # gdf['geometry_simplified'] = gdf['geometry_simplified'].buffer(0)

    return gdf


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


def map_categories(mask: np.ndarray, labels: np.ndarray, categories: np.ndarray, missing_val=-1):
    flat_mask = mask.ravel()

    idx = pd.Index(labels)
    pos = idx.get_indexer(flat_mask)  # -1 where not found

    out_flat = np.full(flat_mask.shape, missing_val, dtype=int)

    valid = pos != -1
    if valid.any():
        out_flat[valid] = categories[pos[valid]]

    return out_flat.reshape(mask.shape)


def hex_to_rgb(hex_color, normalized=False):
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    if normalized:
        return tuple(v / 255 for v in rgb)
    return rgb


def map_clusters_to_mask(meta: pl.DataFrame, cluster_key: str, mask: np.ndarray):
    cluster_cat = pd.Categorical(meta[cluster_key])

    if 'unassigned' not in cluster_cat.categories:
        cluster_cat = cluster_cat.add_categories('unassigned')

    cats = list(cluster_cat.categories)
    cats = ['unassigned'] + [c for c in cats if c != 'unassigned']
    cluster_cat = cluster_cat.reorder_categories(cats)

    cluster_codes = cluster_cat.codes  # integers 0..n-1
    p_mask = map_categories(mask=mask, labels=meta['cell_id'].to_numpy(), categories=cluster_codes)

    pal = [c.UNASSIGNED_COLOR] + c.SG_PALETTE
    pal = np.array([hex_to_rgb(c, normalized=False) for c in pal])

    rgb = pal[p_mask]
    rgb[p_mask == -1] = [0, 0, 0]

    return rgb
