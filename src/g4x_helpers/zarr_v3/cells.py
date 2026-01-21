import numpy as np
import polars as pl
from numcodecs import Blosc

from g4x_helpers.modules import segment
from g4x_helpers.zarr_v3.utils import create_array


def setup_cell_group(root_group):
    cell_group = root_group.create_group('cells', overwrite=True)

    return cell_group


def write_csr(group, csr, gene_names, compressor=None, chunks=None):
    create_array(group, 'data', data=csr.data.astype('int16'), compressor=compressor, chunks=chunks)
    create_array(group, 'indices', data=csr.indices.astype('int32'), compressor=compressor, chunks=chunks)
    create_array(group, 'indptr', data=csr.indptr.astype('int64'), compressor=compressor, chunks=chunks)
    group.attrs['shape'] = csr.shape
    # TODO find maybe better spot for dtype conversion
    create_array(group, 'gene_names', data=np.array(gene_names).astype('U12'), compressor=compressor)


def write_cells(smp, metadata, gex, gene_names, root_group):
    cell_group = setup_cell_group(root_group)
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    # generate polygon vertices
    seg_mask = smp.load_segmentation()
    gdf = segment.vectorize_mask(seg_mask)

    # Keep only largest polygon per label
    gdf['_area'] = gdf.geometry.area
    gdf = (
        gdf.sort_values('_area', ascending=False)
        .drop_duplicates(subset='label', keep='first')
        .drop(columns='_area')
        .reset_index(drop=True)
    )

    # Alternative way to keep largest polygon per label
    # idx = gdf.groupby("label")["_area"].idxmax()
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

    vertices = pl.from_pandas(gdf[['label', 'vert_x', 'vert_y']])

    metadata = metadata.join(vertices, left_on='segmentation_label', right_on='label')
    metadata = metadata.sort('segmentation_label')

    # write arrays for each attribute

    metadata_group = cell_group.create_group('metadata', overwrite=True)

    arr = metadata['segmentation_label'].to_numpy().astype(np.uint32)
    create_array(metadata_group, 'cell_id', data=arr, compressor=compressor)

    arr = metadata['area'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'area', data=arr, compressor=compressor)

    arr = metadata['cluster_id'].to_numpy().astype('U10')
    create_array(metadata_group, 'cluster_id', data=arr, compressor=compressor)

    arr = metadata['total_counts'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'total_counts', data=arr, compressor=compressor)

    arr = metadata['n_genes_by_counts'].to_numpy().astype(np.uint16)
    create_array(metadata_group, 'total_genes', data=arr, compressor=compressor)

    protein_group = cell_group.create_group('protein', overwrite=True)
    protein_df = metadata.select(*[c for c in metadata.columns if '_intensity_mean' in c])
    arr = protein_df.to_numpy().astype(np.uint16)
    create_array(protein_group, 'protein_values', data=arr, compressor=compressor)
    protein_names = np.array([s.removesuffix('_intensity_mean') for s in protein_df.columns]).astype('U50')
    create_array(protein_group, 'protein_names', data=protein_names, compressor=compressor)

    umap = ['X_umap_0.300_0.900_1', 'X_umap_0.300_0.900_2']  # TODO <- remove hard code
    arr = metadata.select(umap).rename({umap[0]: 'umap_1', umap[1]: 'umap_2'}).to_numpy().astype(np.float16)
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

    polygon_group = cell_group.create_group('polygons', overwrite=True)
    create_array(polygon_group, 'polygon_offsets', data=offsets, compressor=compressor)
    create_array(polygon_group, 'polygon_vertices_xy', data=verts_xy, compressor=compressor)

    gex_group = cell_group.create_group('gex', overwrite=True)
    write_csr(gex_group, csr=gex, gene_names=gene_names, compressor=compressor, chunks=512)
