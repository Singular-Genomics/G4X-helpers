import logging
from functools import lru_cache, partial
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
import zarr
from numcodecs import Blosc

from ... import c, io
from ...schema.definition import CellMetadata, CellxGene, CellxProt, ClusteringUmap, Segmentation
from ..workflow import PRESET_SOURCE, collect_input
from .utils import create_array

if TYPE_CHECKING:
    from ...g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)
UNASSIGNED_CELL = 'unassigned'


def write_cells(
    smp: 'G4Xoutput',
    seg_name: str,
    components: tuple | None = None,
    cell_group: zarr.Group | None = None,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
):
    log = LOGGER or logger

    if cell_group is None:
        cell_group = zarr.open_group(smp.out.ViewerZarr.p / 'cells', mode='a')

    seg_path = _sanitize_path_component(seg_name) + '_segmentation'

    seg_sources = cell_group.attrs['segmentation_sources']
    seg_order = cell_group.attrs['segmentation_order']

    seg_sources.update({seg_name: seg_path})
    seg_order.append(seg_name)

    cell_group.attrs['segmentation_sources'] = seg_sources
    cell_group.attrs['segmentation_order'] = list(set(seg_order))

    log.info('Preparing cell data')
    seg_group = cell_group.create_group(seg_path, overwrite=overwrite)
    metadata_group = seg_group.create_group('metadata', overwrite=overwrite)
    protein_group = seg_group.create_group('protein', overwrite=overwrite)
    polygon_group = seg_group.create_group('polygons', overwrite=overwrite)
    genes_group = seg_group.create_group('genes', overwrite=overwrite)

    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    if components is None:
        log.info('Preparing cell group input')
        metadata, gex, gene_names = process_cell_data(smp)
    else:
        log.info('Using provided cell group input')
        metadata, gex, gene_names = components

    clusterings = [c for c in metadata.columns if c.startswith('leiden_')]
    clusterings_order = get_sorted_clusterings(metadata, clusterings)

    meta_columns = {
        'cell_id': metadata[c.CELL_ID_NAME].to_numpy().astype(np.uint32),
        'area': metadata[c.CELL_AREA_NAME].to_numpy().astype(np.uint16),
        'cluster_id': metadata.select(clusterings_order).to_numpy().astype('U20'),
        'total_counts': metadata['total_counts'].to_numpy().astype(np.uint16),
        'total_genes': metadata['n_genes_by_counts'].to_numpy().astype(np.uint16),
        'umap': metadata.select(['UMAP1', 'UMAP2'])
        .rename({'UMAP1': 'umap_1', 'UMAP2': 'umap_2'})
        .to_numpy()
        .astype(np.float16),
    }

    log.info('Writing cell metadata arrays')
    for key, arr in meta_columns.items():
        if overwrite and key in metadata_group:
            del metadata_group[key]
        create_array(metadata_group, key, data=arr, compressor=compressor)

    log.info('Writing cell protein arrays')
    for key in ['protein_values', 'protein_names']:
        if overwrite and key in protein_group:
            del protein_group[key]

    protein_df = metadata.select(*[c for c in metadata.columns if c.IMG_INTENSITY_HANDLE in c])
    arr = protein_df.to_numpy().astype(np.uint16)
    create_array(protein_group, 'protein_values', data=arr, compressor=compressor)
    protein_names = np.array([s.removesuffix(c.IMG_INTENSITY_HANDLE) for s in protein_df.columns]).astype('U50')
    create_array(protein_group, 'protein_names', data=protein_names, compressor=compressor)

    log.info('Writing cell polygon arrays')
    for key in ['polygon_offsets', 'polygon_vertices_xy']:
        if overwrite and key in polygon_group:
            del polygon_group[key]

    # write ragged array of polygon vertices
    x_vert = metadata['vert_x'].to_numpy()
    y_vert = metadata['vert_y'].to_numpy()

    n_cells = len(x_vert)
    lengths = np.array([p.shape[0] for p in x_vert], dtype=np.uint8)
    offsets = np.empty(n_cells + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(lengths, out=offsets[1:])

    x_long = np.concatenate(x_vert, axis=0)
    y_long = np.concatenate(y_vert, axis=0)
    verts_xy = np.stack([x_long, y_long], axis=1)

    create_array(polygon_group, 'polygon_offsets', data=offsets, compressor=compressor)
    create_array(polygon_group, 'polygon_vertices_xy', data=verts_xy, compressor=compressor)

    log.info('Writing cell gene expression arrays')
    for key in ['data', 'indices', 'indptr', 'gene_names']:
        if overwrite and key in genes_group:
            del genes_group[key]

    cluster_labels_meta = {}
    for i, key in enumerate(clusterings_order):
        sorted_cluster_ids = get_sorted_cluster_ids(metadata, cluster_key=key)
        cluster_color_map = generate_cluster_palette(sorted_cluster_ids)

        cluster_labels_meta[key] = {
            'index': i,
            'clusterID_order': list(cluster_color_map.keys()),
            'clusterID_colors': cluster_color_map,
        }

    seg_group['metadata'].attrs['cluster_labels'] = cluster_labels_meta
    seg_group['metadata'].attrs['cluster_labels_order'] = clusterings_order

    seg_group['genes'].attrs['shape'] = gex.shape

    write_csr(genes_group, csr=gex, gene_names=gene_names, compressor=compressor, chunks='auto')


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
    # mask = smp.load_segmentation(expanded=True)
    vertices = extract_vertices_cached(segment_in.p, key=segment_in.main_key, shape=smp.shape)
    cell_metadata = cell_metadata.join(vertices, on=c.CELL_ID_NAME)
    del vertices  # , mask

    # cluster_cols = [c for c in clust_umap.columns if 'leiden_' in c]
    clust_umap = clustumap_in.load()
    cell_metadata = cell_metadata.join(clust_umap, on=c.CELL_ID_NAME, how='left')
    cell_metadata = cell_metadata.with_columns(pl.col('^leiden_.*$').fill_null(UNASSIGNED_CELL))

    return cell_metadata, gex, gene_names


def write_csr(group, csr, gene_names, compressor=None, chunks=None):
    create_array(group, 'data', data=csr.data.astype('int16'), compressor=compressor, chunks=chunks)
    create_array(group, 'indices', data=csr.indices.astype('int32'), compressor=compressor, chunks=chunks)
    create_array(group, 'indptr', data=csr.indptr.astype('int64'), compressor=compressor, chunks=chunks)

    # TODO find maybe better spot for assigning this dtype
    create_array(group, 'gene_names', data=np.array(gene_names).astype('U12'), compressor=compressor, chunks='auto')


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


@lru_cache(maxsize=None)
def extract_vertices_cached(mask, shape, key, show_progress: bool = False):
    mask = io.import_segmentation(mask, expected_shape=shape, labels_key=key, use_cache=True)
    return extract_vertices(mask, show_progress=show_progress)


def extract_vertices(mask, show_progress: bool = False):
    gdf = io.convert.ndarray_to_gdf(mask, show_progress=show_progress)

    # Keep only largest polygon per label
    gdf['_area'] = gdf.geometry.area
    gdf = (
        gdf.sort_values('_area', ascending=False)
        .drop_duplicates(subset=c.CELL_ID_NAME, keep='first')
        .drop(columns='_area')
        .reset_index(drop=True)
    )

    # Alternative way to keep largest polygon per label
    # idx = gdf.groupby(CELL_ID_NAME)["_area"].idxmax()
    # gdf = gdf.loc[idx].drop(columns="_area").reset_index(drop=True)

    # Simplify geometries
    gdf['geometry_simplified'] = gdf.geometry.simplify(tolerance=1.5, preserve_topology=True)
    gdf['geometry_simplified'] = gdf['geometry_simplified'].buffer(0)

    # TODO this section has a lot of room for optimization
    # for now it ensures that the vertices are in the same order as the metadata rows
    res = {
        'x': [poly.exterior.xy[0].tolist() for poly in gdf.geometry_simplified],
        'y': [poly.exterior.xy[1].tolist() for poly in gdf.geometry_simplified],
    }

    gdf['vert_x'] = res['x']
    gdf['vert_y'] = res['y']

    vertices = pl.from_pandas(gdf[[c.CELL_ID_NAME, 'vert_x', 'vert_y']]).sort(c.CELL_ID_NAME)
    return vertices


def _sanitize_path_component(s, replacement='_'):
    invalid = r'<>:"/\\|?*- '
    for ch in invalid:
        s = s.replace(ch, replacement)
    return s.strip(' .')  # Windows disallows trailing space/dot
