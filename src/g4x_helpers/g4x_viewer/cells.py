import numpy as np
import polars as pl
from numcodecs import Blosc

from .. import io
from ..modules.single_cell import CELL_ID_NAME

DEFAULT_COLOR = '#BFBFBF'

sg_palette = [
    '#72FFAB',
    '#A16CFD',
    '#FF7043',
    '#008FFF',
    '#D32F2F',
    '#7CB342',
    '#7F34BE',
    '#FFCA28',
    '#0C8668',
    '#FB4695',
    '#005EE1',
    '#28EDED',
    '#A17B64',
    '#FFFF58',
    '#BC29AE',
    '#006D8F',
    '#FFBAFF',
    '#FFD091',
    '#5C6BC0',
    '#F490B2',
]


def setup_cell_group(root_group):
    cell_group = root_group.create_group('cells', overwrite=True)

    return cell_group


def write_csr(group, csr, gene_names, compressor=None, chunks=None):
    group.create_array('data', data=csr.data.astype('int16'), compressor=compressor, chunks=chunks)
    group.create_array('indices', data=csr.indices.astype('int32'), compressor=compressor, chunks=chunks)
    group.create_array('indptr', data=csr.indptr.astype('int64'), compressor=compressor, chunks=chunks)
    group.attrs['shape'] = csr.shape
    # TODO find maybe better spot for dtype conversion
    group.create_array('gene_names', data=np.array(gene_names).astype('U12'), compressor=compressor, chunks='auto')


def write_cells(smp, metadata, gex, gene_names, root_group):
    cell_group = setup_cell_group(root_group)
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    # generate polygon vertices
    seg_mask = smp.load_segmentation()
    vertices = extract_vertices(seg_mask)

    metadata = metadata.join(vertices, left_on=CELL_ID_NAME, right_on='label')
    metadata = metadata.sort(CELL_ID_NAME)

    # write arrays for each attribute

    metadata_group = cell_group.create_group('metadata', overwrite=True)

    arr = metadata[CELL_ID_NAME].to_numpy().astype(np.uint32)
    metadata_group.create_array('cell_id', data=arr, compressor=compressor)

    arr = metadata['area_um'].to_numpy().astype(np.uint16)
    metadata_group.create_array('area', data=arr, compressor=compressor)

    # TODO another hard coding that needs removal later
    metadata = metadata.with_columns(cluster_id=pl.col('leiden_0.400').cast(pl.Int8).cast(pl.String).fill_null('-1'))
    metadata_group.attrs['clusterID_colors'] = generate_cluster_palette(metadata['cluster_id'])

    arr = metadata['cluster_id'].to_numpy().astype('U10')
    metadata_group.create_array('cluster_id', data=arr, compressor=compressor)

    arr = metadata['total_counts'].to_numpy().astype(np.uint16)
    metadata_group.create_array('total_counts', data=arr, compressor=compressor)

    arr = metadata['n_genes_by_counts'].to_numpy().astype(np.uint16)
    metadata_group.create_array('total_genes', data=arr, compressor=compressor)

    protein_group = cell_group.create_group('protein', overwrite=True)
    protein_df = metadata.select(*[c for c in metadata.columns if '_intensity_mean' in c])
    arr = protein_df.to_numpy().astype(np.uint16)
    protein_group.create_array('protein_values', data=arr, compressor=compressor)
    protein_names = np.array([s.removesuffix('_intensity_mean') for s in protein_df.columns]).astype('U50')
    protein_group.create_array('protein_names', data=protein_names, compressor=compressor)

    umap = ['X_umap_0.300_0.900_1', 'X_umap_0.300_0.900_2']  # TODO <- remove hard code
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

    polygon_group = cell_group.create_group('polygons', overwrite=True)
    polygon_group.create_array('polygon_offsets', data=offsets, compressor=compressor)
    polygon_group.create_array('polygon_vertices_xy', data=verts_xy, compressor=compressor)

    genes_group = cell_group.create_group('genes', overwrite=True)
    write_csr(genes_group, csr=gex, gene_names=gene_names, compressor=compressor, chunks='auto')


def generate_cluster_palette(clusters: list, max_colors: int = 256) -> dict:
    """
    Generate a color palette mapping for cluster labels.

    This function assigns RGB colors to unique cluster labels using a matplotlib colormap.
    Clusters labeled as "-1" are assigned a default gray color `[191, 191, 191]`.

    The colormap used depends on the number of clusters:
        - `tab10` for ≤10 clusters
        - `tab20` for ≤20 clusters
        - `hsv` for more than 20 clusters, capped by `max_colors`

    Parameters
    ----------
    clusters : list
        A list of cluster identifiers (strings or integers). The special label '-1' is excluded
        from color mapping and handled separately.
    max_colors : int, optional
        Maximum number of colors to use in the HSV colormap. Only used if there are more than
        20 unique clusters. Default is 256.

    Returns
    -------
    dict
        A dictionary mapping each cluster ID (as a string) to a list of three integers
        representing an RGB color in the range [0, 255].

    Examples
    --------
    >>> generate_cluster_palette(['0', '1', '2', '-1'])
    {'0': [31, 119, 180], '1': [255, 127, 14], '2': [44, 160, 44], '-1': [191, 191, 191]}
    """

    import matplotlib.colors as mcolors

    def hex2rgb(hex: str) -> list[int, int, int]:
        return [int(x * 255) for x in mcolors.to_rgb(hex)]

    unique_clusters = [c for c in np.unique(clusters) if c != '-1']
    n_clusters = len(unique_clusters)

    if n_clusters <= 20:
        hex_list = sg_palette

    else:
        from matplotlib.pyplot import get_cmap

        hex_list = get_cmap('hsv', min(max_colors, n_clusters)).colors

    cluster_palette = {}
    for i, cluster in enumerate(unique_clusters):
        cluster_palette[str(cluster)] = hex2rgb(hex_list[i])

    cluster_palette['-1'] = hex2rgb(DEFAULT_COLOR)

    return cluster_palette


def extract_vertices(mask):
    gdf = io.convert.ndarray_to_gdf(mask)

    # Keep only largest polygon per label
    gdf['_area'] = gdf.geometry.area
    gdf = (
        gdf.sort_values('_area', ascending=False)
        .drop_duplicates(subset=CELL_ID_NAME, keep='first')
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

    vertices = pl.from_pandas(gdf[[CELL_ID_NAME, 'vert_x', 'vert_y']]).sort(CELL_ID_NAME)
    return vertices
