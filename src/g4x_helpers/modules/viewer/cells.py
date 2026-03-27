import logging
from functools import lru_cache

import numpy as np
import polars as pl
from numcodecs import Blosc

from ... import c, io
from .setup import create_array, populate_zarr_metadata

LOGGER = logging.getLogger(__name__)


def write_cells(smp, root_group, prechew: None = None, logger: logging.Logger | None = None):
    log = LOGGER or logger

    log.info('Preparing cell data')
    cell_group = root_group['cells']
    metadata_group = cell_group['metadata']
    protein_group = cell_group['protein']
    polygon_group = cell_group['polygons']
    genes_group = cell_group['genes']

    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    log.info('Preparing cell group input')
    if prechew is None:
        metadata, gex, gene_names = prepare_cell_group_input(smp)
    else:
        metadata, gex, gene_names = prechew

    log.info('Writing cell metadata arrays')
    # write arrays for each attribute
    arr = metadata[c.CELL_ID_NAME].to_numpy().astype(np.uint32)
    create_array(metadata_group, 'cell_id', data=arr, compressor=compressor)

    arr = metadata['wholecell_area_um'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'area', data=arr, compressor=compressor)

    arr = metadata['leiden_fine'].to_numpy().astype('U10')
    create_array(metadata_group, 'cluster_id', data=arr, compressor=compressor)

    arr = metadata['total_counts'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'total_counts', data=arr, compressor=compressor)

    arr = metadata['n_genes_by_counts'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'total_genes', data=arr, compressor=compressor)

    protein_df = metadata.select(*[c for c in metadata.columns if '_intensity_mean' in c])
    arr = protein_df.to_numpy().astype(np.uint16)
    create_array(protein_group, 'protein_values', data=arr, compressor=compressor)
    protein_names = np.array([s.removesuffix('_intensity_mean') for s in protein_df.columns]).astype('U50')
    create_array(protein_group, 'protein_names', data=protein_names, compressor=compressor)

    arr = (
        metadata.select(['UMAP1', 'UMAP2']).rename({'UMAP1': 'umap_1', 'UMAP2': 'umap_2'}).to_numpy().astype(np.float16)
    )
    create_array(metadata_group, 'umap', data=arr, compressor=compressor)

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

    uids = sort_clusters(metadata, key='leiden_fine')
    populate_zarr_metadata(root_group, gene_mtx_shape=gex.shape, cluster_ids=generate_cluster_palette(uids))
    write_csr(genes_group, csr=gex, gene_names=gene_names, compressor=compressor, chunks='auto')


def write_csr(group, csr, gene_names, compressor=None, chunks=None):
    create_array(group, 'data', data=csr.data.astype('int16'), compressor=compressor, chunks=chunks)
    create_array(group, 'indices', data=csr.indices.astype('int32'), compressor=compressor, chunks=chunks)
    create_array(group, 'indptr', data=csr.indptr.astype('int64'), compressor=compressor, chunks=chunks)

    # TODO find maybe better spot for dtype conversion
    create_array(group, 'gene_names', data=np.array(gene_names).astype('U12'), compressor=compressor, chunks='auto')


def sort_clusters(metadata, key='cluster_id'):
    uids = metadata.group_by(key).agg(pl.len())
    uids = uids.sort('len', descending=True)[key].to_list()
    if 'unassigned' in uids:
        uids.remove('unassigned')

    uids.append('unassigned')
    return uids


def prepare_cell_group_input(smp):
    from scipy.sparse import csr_matrix

    cell_metadata = smp.tree.CellMetadata.load()
    cell_x_gene = smp.tree.CellxGene.load()
    cell_x_protein = smp.tree.CellxProtein.load()

    if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_protein[c.CELL_ID_NAME]):
        raise ValueError('Cell IDs in metadata and protein table do not match')

    if not cell_metadata[c.CELL_ID_NAME].equals(cell_x_gene[c.CELL_ID_NAME]):
        raise ValueError('Cell IDs in metadata and gene expression matrix do not match')

    cell_x_gene = cell_x_gene.drop(c.CELL_ID_NAME)
    gene_names = np.array(cell_x_gene.columns)

    gex = cell_x_gene.to_numpy().astype(np.uint16)
    gex = csr_matrix(gex)
    del cell_x_gene

    # 2: extract cell vertices
    # mask = smp.load_segmentation(expanded=True)
    vertices = extract_vertices_cached(smp.tree.Segmentation.p, shape=smp.shape)
    # del mask

    cell_metadata = cell_metadata.join(cell_x_protein, on=c.CELL_ID_NAME)
    del cell_x_protein

    cell_metadata = cell_metadata.join(vertices, on=c.CELL_ID_NAME)
    del vertices

    clust_umap = pl.read_csv(smp.tree.ClusteringUmap.p)
    # cluster_cols = [c for c in clust_umap.columns if 'leiden_' in c]

    cell_metadata = cell_metadata.join(clust_umap, on=c.CELL_ID_NAME, how='left')
    cell_metadata = cell_metadata.with_columns(pl.col('^leiden_.*$').fill_null('unassigned'))
    return cell_metadata, gex, gene_names


# def prepare_cell_group_input(smp):
#     # 1: create fresh adata with qc-values
#     tx_panel = smp.tree.TranscriptPanel.p
#     cell_metadata = smp.tree.CellMetadata.p
#     cell_x_gene = smp.tree.CellxGene.p

#     adata = single_cell.init_adata(smp, tx_panel=tx_panel, cell_metadata=cell_metadata, cell_x_gene=cell_x_gene)

#     gex = adata.X
#     del cell_x_gene
#     gene_names = adata.var_names

#     # 2: extract cell vertices
#     mask = smp.load_segmentation(expanded=True)
#     vertices = extract_vertices(mask)
#     del mask

#     # 3: create metadata dataframe
#     obs_df = (
#         pl.from_pandas(adata.obs, include_index=True)
#         .drop('tissue_type', 'block', 'seg_source')
#         .cast({c.CELL_ID_NAME: pl.UInt32})
#     )
#     del adata

#     metadata = obs_df.join(vertices, on=c.CELL_ID_NAME)
#     metadata = metadata.sort(c.CELL_ID_NAME)
#     del obs_df, vertices

#     # 4: load clustering and UMAP coordinates
#     clust_umap = pl.read_csv(
#         smp.tree.ClusteringUmap.p,
#         schema={
#             'idx': pl.UInt16,
#             c.CELL_ID_NAME: pl.UInt32,
#             'leiden_1.00': pl.Utf8,
#             'UMAP1': pl.Float32,
#             'UMAP2': pl.Float32,
#         },
#     )
#     clust_umap = clust_umap.rename({'leiden_1.00': 'cluster_id'})  # TODO this is hard coded right now

#     metadata = metadata.join(clust_umap, on=c.CELL_ID_NAME, how='left')
#     metadata = metadata.with_columns(pl.col('cluster_id').fill_null('unassigned'))

#     # 5: load protein data if available
#     if smp.tree.pr_detected:
#         protein_data = pl.read_csv(smp.data_dir / 'single_cell_data' / 'cell_by_protein.csv.gz')
#         metadata = metadata.join(protein_data, on=c.CELL_ID_NAME, how='left')

#     return metadata, gex, gene_names


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

    cluster_palette['unassigned'] = hex2rgb(c.UNASSIGNED_COLOR)

    return cluster_palette


@lru_cache(maxsize=None)
def extract_vertices_cached(mask, shape):
    mask = io.import_segmentation(mask, expected_shape=shape, labels_key='nuclei_exp', use_cache=True)
    return extract_vertices(mask)


def extract_vertices(mask):
    gdf = io.convert.ndarray_to_gdf(mask)

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
