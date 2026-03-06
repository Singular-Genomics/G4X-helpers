import numpy as np
import polars as pl
from numcodecs import Blosc

from .. import constants as c
from .. import io
from ..modules import single_cell as g4xsc
from .setup import populate_zarr_metadata


def write_csr(group, csr, gene_names, compressor=None, chunks=None):
    group.create_array('data', data=csr.data.astype('int16'), compressor=compressor, chunks=chunks)
    group.create_array('indices', data=csr.indices.astype('int32'), compressor=compressor, chunks=chunks)
    group.create_array('indptr', data=csr.indptr.astype('int64'), compressor=compressor, chunks=chunks)

    # TODO find maybe better spot for dtype conversion
    group.create_array('gene_names', data=np.array(gene_names).astype('U12'), compressor=compressor, chunks='auto')


def write_cells(smp, root_group):
    cell_group = root_group['cells']
    metadata_group = cell_group['metadata']
    protein_group = cell_group['protein']
    polygon_group = cell_group['polygons']
    genes_group = cell_group['genes']

    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    metadata, gex, gene_names = prepare_cell_group_input(smp)

    # write arrays for each attribute
    arr = metadata[c.CELL_ID_NAME].to_numpy().astype(np.uint32)
    metadata_group.create_array('cell_id', data=arr, compressor=compressor)

    arr = metadata['area_um'].to_numpy().astype(np.uint16)
    metadata_group.create_array('area', data=arr, compressor=compressor)

    arr = metadata['cluster_id'].to_numpy().astype('U10')
    metadata_group.create_array('cluster_id', data=arr, compressor=compressor)

    arr = metadata['total_counts'].to_numpy().astype(np.uint16)
    metadata_group.create_array('total_counts', data=arr, compressor=compressor)

    arr = metadata['n_genes_by_counts'].to_numpy().astype(np.uint16)
    metadata_group.create_array('total_genes', data=arr, compressor=compressor)

    protein_df = metadata.select(*[c for c in metadata.columns if '_intensity_mean' in c])
    arr = protein_df.to_numpy().astype(np.uint16)
    protein_group.create_array('protein_values', data=arr, compressor=compressor)
    protein_names = np.array([s.removesuffix('_intensity_mean') for s in protein_df.columns]).astype('U50')
    protein_group.create_array('protein_names', data=protein_names, compressor=compressor)

    umap = ['UMAP1', 'UMAP2']  # TODO <- remove hard code
    arr = metadata.select(umap).rename({umap[0]: 'umap_1', umap[1]: 'umap_2'}).to_numpy().astype(np.float16)
    metadata_group.create_array('umap', data=arr, compressor=compressor)

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

    polygon_group.create_array('polygon_offsets', data=offsets, compressor=compressor)
    polygon_group.create_array('polygon_vertices_xy', data=verts_xy, compressor=compressor)

    uids = sort_clusters(metadata)
    populate_zarr_metadata(root_group, gene_mtx_shape=gex.shape, cluster_ids=generate_cluster_palette(uids))
    write_csr(genes_group, csr=gex, gene_names=gene_names, compressor=compressor, chunks='auto')


def sort_clusters(metadata):
    uids = metadata.group_by('cluster_id').agg(pl.len())
    uids = uids.sort('len', descending=True)['cluster_id'].to_list()
    if 'unassigned' in uids:
        uids.remove('unassigned')

    uids.append('unassigned')
    return uids


def prepare_cell_group_input(smp):
    # 1: create fresh adata with qc-values
    cell_x_gene = pl.read_csv(smp.data_dir / 'single_cell_data' / 'cell_by_gene.csv.gz')

    mask = smp.load_segmentation(expanded=True)
    cell_metadata = g4xsc.create_cell_metadata(smp, mask=mask)
    adata = g4xsc.create_adata(smp, cell_metadata=cell_metadata, cell_x_gene=cell_x_gene)

    gex = adata.X
    del cell_x_gene
    gene_names = adata.var_names

    # 2: extract cell vertices
    vertices = extract_vertices(mask)
    del mask

    # 3: create metadata dataframe
    obs_df = (
        pl.from_pandas(adata.obs, include_index=True)
        .drop('tissue_type', 'block', 'source')
        .cast({c.CELL_ID_NAME: pl.UInt32})
    )
    del adata

    metadata = obs_df.join(vertices, on=c.CELL_ID_NAME)
    metadata = metadata.sort(c.CELL_ID_NAME)
    del obs_df, vertices

    # 4: load clustering and UMAP coordinates
    clust_umap = pl.read_csv(
        smp.data_dir / 'single_cell_data' / 'clustering_umap.csv.gz',
        schema={c.CELL_ID_NAME: pl.UInt32, 'leiden_1.00': pl.Utf8, 'UMAP1': pl.Float32, 'UMAP2': pl.Float32},
    )
    clust_umap = clust_umap.rename({'leiden_1.00': 'cluster_id'})  # TODO this is hard coded right now

    metadata = metadata.join(clust_umap, on=c.CELL_ID_NAME, how='left')
    metadata = metadata.with_columns(pl.col('cluster_id').fill_null('unassigned'))

    # 5: load protein data if available
    if smp.includes_protein:
        protein_data = pl.read_csv(smp.data_dir / 'single_cell_data' / 'cell_by_protein.csv.gz')
        metadata = metadata.join(protein_data, on=c.CELL_ID_NAME, how='left')

    return metadata, gex, gene_names


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
