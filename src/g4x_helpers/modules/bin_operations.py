import logging
import multiprocessing
from collections import deque
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import polars as pl
from google.protobuf.message import DecodeError
from numba import njit
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance_matrix
from skimage.measure import approximate_polygon
from skimage.morphology import dilation, disk, erosion
from tqdm import tqdm

from .. import utils
from ..g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema
from .resegment import get_cell_ids
from .workflow import workflow

DEFAULT_COLOR = [int(191), int(191), int(191)]


@workflow
def new_bin_core(
    g4x_obj,
    seg_mask: np.ndarray,
    out_dir: str | Path,
    *,
    metadata: str | Path | None = None,
    cluster_key: str | None = None,
    emb_key: str | None = None,
    n_threads: int = 4,
    logger: logging.Logger | None = None,
) -> None:
    # import warnings
    # warnings.filterwarnings(
    #     'ignore',
    #     message='FNV hashing is not implemented in Numba',
    #     category=UserWarning,
    #     module='numba.cpython.old_hashing',
    # )

    logger.info('Creating G4X-viewer bin file.')

    out_tree = utils.OutSchema(out_dir, subdirs=['g4x_viewer'])
    out_file = out_tree.g4x_viewer / f'{g4x_obj.sample_id}_segmentation.bin'

    adata = g4x_obj.load_adata(load_clustering=True, remove_nontargeting=False)

    ### some testing code to simulate missing data
    # adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
    # cells = adata.obs.sample(150000).index
    # adata.obs.loc[cells, cluster_key] = np.nan
    # adata.obs.loc[~adata.obs.index.isin(cells), cluster_key] = 'keep_no_emb'

    # adata.obs.loc[cells, f'{emb_key}_1'] = np.nan
    # adata.obs.loc[cells, f'{emb_key}_2'] = np.nan

    logger.info('Loading clustering information.')
    obs_df = prepare_metadata(
        sample_id=g4x_obj.sample_id, adata=adata, seg_mask=seg_mask, metadata=metadata, cluster_key=cluster_key
    )
    print(f'Number of cells to process: {obs_df.shape[0]}')

    logger.info('Making polygons.')
    cell_ids, pq_args = initialize_segmentation_data(obs_df, seg_mask)

    logger.debug('Adding single-cell metadata.')
    obs_df, gene_names, gex = add_singlecell_info(obs_df, adata, cell_ids, emb_key)

    ## do conversion
    cell_seg_out = CellMasksSchema.CellMasks()
    cell_seg_out.metadata.geneNames.extend(gene_names)

    logger.info('Generating cluster palette.')
    cluster_palette = generate_cluster_palette(obs_df['cluster_id'])

    for cluster_id, color in cluster_palette.items():
        entry = CellMasksSchema.ColormapEntry()
        entry.clusterId = cluster_id
        entry.color.extend(color)
        cell_seg_out.colormap.append(entry)

    protein_values = None
    if g4x_obj.includes_protein:
        protein_names = g4x_obj.proteins
        protein_list = [prot + '_intensity_mean' for prot in protein_names]

        cell_seg_out.metadata.proteinNames.extend(protein_names)

        obs_df = obs_df.with_columns([pl.col(prot).fill_null(0).cast(pl.Int64) for prot in protein_list])
        protein_values = obs_df.select(protein_list).to_numpy()

    logger.info('Refining polygons.')
    with multiprocessing.Pool(processes=n_threads) as pool:
        polygons = pool.starmap(refine_polygon, pq_args)

    logger.info('Adding individual cells.')
    for i, row in enumerate(tqdm(obs_df.iter_rows(named=True), total=obs_df.height, desc='Adding individual cells.')):
        output_mask_data = cell_seg_out.cellMasks.add()
        cell_polygon_pts = [sub_coord for coord in polygons[i] for sub_coord in coord[::-1]]

        output_mask_data.vertices.extend(cell_polygon_pts + cell_polygon_pts[:2])

        output_mask_data.cellId = str(row['cell_id'])
        output_mask_data.area = int(row['area'])
        output_mask_data.totalCounts = int(row['total_counts'])
        output_mask_data.totalGenes = int(row['n_genes_by_counts'])
        output_mask_data.clusterId = str(row['cluster_id'])
        output_mask_data.umapValues.umapX = row['umap_0']
        output_mask_data.umapValues.umapY = row['umap_1']

        start = gex.indptr[i]
        end = gex.indptr[i + 1]
        indices = gex.indices[start:end]
        values = gex.data[start:end]
        output_mask_data.nonzeroGeneIndices.extend(indices.tolist())
        output_mask_data.nonzeroGeneValues.extend(values.astype(int).tolist())

        if protein_values is not None:
            output_mask_data.proteinValues.extend(protein_values[i].astype(int).tolist())

    ## write to file
    with open(out_file, 'wb') as file:
        file.write(cell_seg_out.SerializeToString())

    logger.debug(f'G4X-viewer bin --> {out_file}')


@workflow
def seg_updater(
    bin_file: Path,
    metadata_file: Path,
    out_path: Path,
    *,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    cluster_color_key: str | None = None,
    emb_key: str | None = None,
    logger: logging.Logger,
) -> None:
    ## pre-flight
    if emb_key is None and cluster_key is None:
        logger.warning('neither embedding nor cluster keys were provided, nothing to update.')
        return None

    ## load the bin file
    logger.info(f'Loading {bin_file}.')
    with open(bin_file, 'rb') as f:
        data = f.read()

    cell_masks = CellMasksSchema.CellMasks()
    cell_masks.ParseFromString(data)

    ## load the metadata
    if cellid_key is None:
        logger.info('cellid_key not provided, assuming cell IDs are in first column of metadata.')
        metadata = pd.read_csv(metadata_file, index_col=0, header=0)
    else:
        metadata = pd.read_csv(metadata_file, index_col=None, header=0)
        if cellid_key not in metadata.columns:
            raise KeyError(f'{cellid_key} not a valid column in metadata.')
        metadata.set_index(cellid_key, inplace=True)

    ## check for clustering
    if cluster_key is not None:
        if cluster_key not in metadata.columns:
            raise KeyError(f'{cluster_key} not a valid column in metadata.')
        update_cluster = True
        logger.debug('Updating cluster IDs.')
    else:
        update_cluster = False
        logger.debug('Not updating cluster IDs.')

    ## check for cluster colors
    if cluster_color_key is not None:
        if cluster_key is None:
            raise ValueError('cluster_color_key was provided, but cluster_key was not provided.')
        if cluster_color_key not in metadata.columns:
            raise KeyError(f'{cluster_color_key} not a valid column in metadata.')
        color = metadata[cluster_color_key].iat[0]
        assert color.startswith('#'), 'Cluster colors must be provided as hexcodes.'
        update_cluster_color = True
        logger.debug('Updating cluster colors.')
        cluster_palette = (
            metadata.drop_duplicates(subset=cluster_key)[[cluster_key, cluster_color_key]]
            .set_index(cluster_key)
            .to_dict()[cluster_color_key]
        )
        cluster_palette = {str(k): hex2rgb(v) for k, v in cluster_palette.items()}
    else:
        if cluster_key is not None:
            update_cluster_color = True
            logger.debug('Auto-assigning colors to new clustering.')
            cluster_color_key = 'cluster_color'
            cluster_palette = utils.generate_cluster_palette(metadata[cluster_key].tolist())
            metadata['cluster_color'] = metadata[cluster_key].astype(str).map(cluster_palette).tolist()
        else:
            update_cluster_color = False
            logger.debug('Not updating cluster colors.')

    ## check for embedding
    if emb_key is not None:
        if f'{emb_key}_1' not in metadata.columns or f'{emb_key}_2' not in metadata.columns:
            raise KeyError(f'{emb_key}_1 and {emb_key}_2 are not valid columns in metadata.')
        update_emb = True
        logger.debug('Updating embedding.')
    else:
        update_emb = False
        logger.debug('Not updating embedding.')

    ## Do the actual updating
    logger.info('Updating cells.')
    for cell in tqdm(cell_masks.cellMasks, desc='Updating cell data'):
        current_cellid = cell.cellId
        if current_cellid in metadata.index:
            if update_cluster:
                cell.clusterId = str(metadata.loc[current_cellid, cluster_key])

            # if update_cluster_color and bin_version == 'v2':
            #     # clear out the existing color entries:
            #     cell.ClearField('color')
            #     cell.color.extend(metadata.loc[current_cellid, cluster_color_key])

            if update_emb:
                cell.umapValues.umapX = metadata.loc[current_cellid, f'{emb_key}_1']
                cell.umapValues.umapY = metadata.loc[current_cellid, f'{emb_key}_2']
        else:
            logger.debug(f'{current_cellid} not found in metadata, not updating data for this cell.')

    # clear the entire colormap list:
    if update_cluster_color:
        cell_masks.ClearField('colormap')
        for cluster_id, color in cluster_palette.items():
            entry = CellMasksSchema.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            cell_masks.colormap.append(entry)

    ## Write to file
    logger.debug(f'Writing updated bin file --> {out_path}')
    with open(out_path, 'wb') as file:
        file.write(cell_masks.SerializeToString())


def read_bin_file(bin_file: Path) -> CellMasksSchema.CellMasks:
    with open(bin_file, 'rb') as f:
        data = f.read()

    cell_masks = CellMasksSchema.CellMasks()
    try:
        cell_masks.ParseFromString(data)
        return cell_masks
    except DecodeError:
        print('Failed to parse bin file with current schema. It may need to be updated.')
        return None


def hex2rgb(hex: str) -> list[int, int, int]:
    return [int(x * 255) for x in mcolors.to_rgb(hex)]


def prepare_metadata(sample_id, adata, seg_mask, metadata=None, cluster_key=None):
    obs_df = get_cell_ids(sample_id=sample_id, mask=seg_mask)

    if metadata is None:  # get metadata from adata
        metadata_df = pl.from_pandas(adata.obs)
    else:
        metadata_df = pl.read_csv(metadata)

    obs_df = obs_df.join(metadata_df, left_on='label', right_on='cell_id', how='left').rename({'label': 'cell_id'})

    if not cluster_key:
        print('No cluster_key provided, attempting to infer from metadata.')
        cluster_key = sorted([x for x in obs_df.columns if 'leiden' in x])
        if len(cluster_key) == 0:
            print('Failed to infer cluster_key, no clustering information will be added.')
            clusters_available = False
        else:
            cluster_key = cluster_key[0]
            print(f'Inferred cluster_key: {cluster_key}')
            clusters_available = True
    else:
        clusters_available = cluster_key in obs_df.columns

    if clusters_available:
        obs_df = obs_df.rename({cluster_key: 'cluster_id'}).cast({'cluster_id': pl.Utf8})
        obs_df = obs_df.with_columns(pl.col('cluster_id').fill_null('-1'))
    else:
        obs_df = obs_df.with_columns(pl.lit('-1').alias('cluster_id'))

    obs_df.select('cell_id', 'cluster_id').sort('cluster_id')

    return obs_df


def load_clustered_obs(adata, metadata=None, cluster_key=None):
    if metadata is None:
        if cluster_key in adata.obs.columns:
            clusters_available = True
        else:
            clusters_available = False
        obs_df = adata.obs.copy()
    else:
        clusters_available = True
        clustered_df = pd.read_csv(metadata, index_col=0, header=0)
        if clustered_df.shape[1] > 1:
            assert cluster_key is not None, (
                'ERROR: multiple columns detected in cluster_info, cluster_key must be provided.'
            )
        else:
            cluster_key = clustered_df.columns[0]
        orig_df = adata.obs.copy()

        ## these are cells that were filtered out during clustering
        orig_df = orig_df.loc[list(set(orig_df.index) - set(clustered_df.index)), :].copy()
        for col in list(set(clustered_df.columns) - set(orig_df.columns)):
            orig_df[col] = '-1'

        obs_df = pd.concat([clustered_df, orig_df])
        obs_df.sort_index(inplace=True)

    if clusters_available:
        obs_df['cluster_id'] = obs_df[cluster_key].astype(str)
    else:
        obs_df['cluster_id'] = '-1'

    return obs_df


def initialize_segmentation_data(obs_df, seg_mask):
    ## we create polygons to define the boundaries of each cell mask
    border = get_border(seg_mask)
    seg_mask[border > 0] = 0
    eroded_mask = erosion(seg_mask, disk(1))
    outlines = seg_mask - eroded_mask
    sparse_matrix = csr_matrix(outlines)
    del seg_mask, border, eroded_mask, outlines

    nonzero_values = sparse_matrix.data
    nonzero_row_indices, nonzero_col_indices = sparse_matrix.nonzero()
    sorted_indices = np.argsort(nonzero_values)
    sorted_nonzero_values = nonzero_values[sorted_indices]
    sorted_rows = nonzero_row_indices[sorted_indices]
    sorted_cols = nonzero_col_indices[sorted_indices]

    cell_ids = obs_df['cell_id'].to_list()
    segmentation_labels = obs_df['segmentation_label'].to_list()

    centroid_y = obs_df['cell_x'].to_list()
    centroid_x = obs_df['cell_y'].to_list()

    pq_args = [
        (seg_label, sorted_nonzero_values, sorted_rows, sorted_cols)
        for seg_label, _, _ in zip(segmentation_labels, centroid_x, centroid_y)
    ]

    return cell_ids, pq_args


def add_singlecell_info(obs_df, adata, cell_ids, emb_key):
    # filter adata
    adata = adata[cell_ids]
    gene_names = adata.var_names

    gex = adata.X
    if not issparse(gex):
        print('Converting dense expression matrix to CSR format.')
        gex = csr_matrix(gex)
    else:
        gex = gex.tocsr()

    del adata

    obs_df = obs_df.cast({'total_counts': pl.UInt32, 'n_genes_by_counts': pl.Int32})

    # add_singlecell_area
    area_cols = [c for c in obs_df.columns if 'area' in c]
    if not len(area_cols):
        raise ValueError('No area column found')

    if len(area_cols) == 1:
        area_select = area_cols[0]
    else:
        expanded_area_cols = [f for f in area_cols if 'expanded' in f]
        if len(expanded_area_cols) != 1:
            raise ValueError(
                f'No expanded area column found among multiple area columns: {area_cols}.'
                if not expanded_area_cols
                else f'Multiple expanded area columns found: {expanded_area_cols}'
            )
        print(f'Multiple area columns found. Using {expanded_area_cols[0]}.')
        area_select = expanded_area_cols[0]

    print(f'Using area column: {area_select}')

    obs_df = obs_df.cast({area_select: pl.Int32}).rename({area_select: 'area'})

    def add_umap_coordinates(obs_df, emb_key):
        obs_df = obs_df.with_columns(umap_0=obs_df[f'{emb_key}_1'], umap_1=obs_df[f'{emb_key}_2'])
        obs_df = obs_df.with_columns(pl.col('umap_0').fill_null(float('nan')), pl.col('umap_1').fill_null(float('nan')))
        return obs_df

    if emb_key:
        obs_df = add_umap_coordinates(obs_df, emb_key)
    else:
        emb_keys = [c[:-2] for c in obs_df.columns if c.startswith('X_umap')]
        if len(emb_keys) == 0:
            print('No embedding key available, UMAP coordinates will be set to NaN.')
            obs_df = obs_df.with_columns(umap_0=float('nan'), umap_1=float('nan'))

        else:
            emb_key = emb_keys[0]
            obs_df = add_umap_coordinates(obs_df, emb_key)
            print(f'No embedding key provided. Using: {emb_key}')

    return obs_df, gene_names, gex


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
    from matplotlib.pyplot import get_cmap

    unique_clusters = [c for c in np.unique(clusters) if c != '-1']
    n_clusters = len(unique_clusters)

    if n_clusters <= 10:
        base_cmap = get_cmap('tab10')
    elif n_clusters <= 20:
        base_cmap = get_cmap('tab20')
    else:
        base_cmap = get_cmap('hsv', min(max_colors, n_clusters))

    cluster_palette = {
        str(cluster): [int(255 * c) for c in base_cmap(i / n_clusters)[:3]] for i, cluster in enumerate(unique_clusters)
    }
    cluster_palette['-1'] = DEFAULT_COLOR

    return cluster_palette


def refine_polygon(cell_id, sorted_nonzero_values_ref, sorted_rows_ref, sorted_cols_ref):
    start_idx, end_idx = get_start_stop_idx(sorted_nonzero_values_ref, cell_id)
    points = np.vstack((sorted_rows_ref[start_idx:end_idx], sorted_cols_ref[start_idx:end_idx])).T
    return pointsToSingleSmoothPath(points, tolerance=2.0)


@njit
def get_start_stop_idx(arr, k):
    start_idx = np.searchsorted(arr, k, side='left')
    end_idx = np.searchsorted(arr, k, side='right')
    return start_idx, end_idx


def pointsToSingleSmoothPath(points: np.ndarray, tolerance: float) -> np.ndarray:
    # Calculate the distance matrix
    dist_matrix = distance_matrix(points, points)

    # Create a sparse matrix for the MST calculation
    sparse_matrix = csr_matrix(dist_matrix)

    # Compute the MST
    mst = minimum_spanning_tree(sparse_matrix).toarray()

    adj_list, adj_list_pos = createAdjacencyList_numba(mst)
    adj_list = {row: list(adj_list[row, :pos]) for row, pos in enumerate(adj_list_pos) if pos}

    longest_path = computeLongestPath(adj_list)
    bestPath = indicesToArray(points, longest_path)

    simplified_path = simplify_polygon(bestPath, tolerance=tolerance)

    return simplified_path


def computeLongestPath(adj_list):
    endpoints = returnEndpoints(adj_list)
    longest_path = []
    max_length = 0
    # Use a dictionary to cache paths and avoid recomputation
    path_cache = {}
    # Compute distances between all pairs of endpoints
    for i in range(len(endpoints)):
        for j in range(i + 1, len(endpoints)):
            if (endpoints[i], endpoints[j]) not in path_cache:
                path = bfs_path(endpoints[i], endpoints[j], adj_list)
                path_cache[(endpoints[i], endpoints[j])] = path
            else:
                path = path_cache[(endpoints[i], endpoints[j])]
            if len(path) > max_length:
                max_length = len(path)
                longest_path = path

    return longest_path


def returnEndpoints(adj_list, adjacency=2):
    # Identify endpoints of the MST
    endpoints = [node for node in adj_list if len(adj_list[node]) == adjacency]

    return endpoints


def bfs_path(start, end, adj_list):
    queue = deque([(start, [start])])
    visited = set()
    while queue:
        current, path = queue.popleft()
        if current == end:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in adj_list[current]:
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))
    return []


@njit
def createAdjacencyList_numba(mst):
    """
    Create an adjacency list from a minimum spanning tree (MST) using Numba for performance optimization.

    Parameters:
    mst (numpy.ndarray): The minimum spanning tree represented as a 2D numpy array.

    Returns:
    tuple: A tuple containing:
        - adj_list (numpy.ndarray): An array where each row contains the adjacent nodes for each node.
        - adj_list_pos (numpy.ndarray): An array containing the number of adjacent nodes for each node.
    """
    n = mst.shape[0]
    adj_list = np.zeros((n, n * 2), dtype=np.uint32)
    adj_list_pos = np.zeros(n, dtype=np.uint32)
    for i in range(n):
        for j in range(n):
            if mst[i, j] != 0 or mst[j, i] != 0:
                adj_list[i, adj_list_pos[i]] = j
                adj_list_pos[i] += 1
                adj_list[j, adj_list_pos[j]] = i
                adj_list_pos[j] += 1

    return adj_list, adj_list_pos


def indicesToArray(points, longest_path):
    pth = []

    for j in range(len(longest_path)):
        pth.append([points[longest_path[j], 0], points[longest_path[j], 1]])

    return np.array(pth)


def simplify_polygon(points: np.ndarray, tolerance: float) -> np.ndarray:
    """
    Simplify a series of points representing a polygon using scikit-image's
    approximate_polygon. The tolerance controls how aggressively the polygon
    is simplified (in pixel units).
    """
    if len(points) <= 2:
        return points

    # If the first and last points are not the same, append the first to the end
    # to ensure the polygon is "closed" for approximation (optional).
    if not np.array_equal(points[0], points[-1]):
        points = np.vstack([points, points[0]])

    # Perform polygon simplification
    simplified = approximate_polygon(points, tolerance=tolerance)

    # If approximate_polygon returns only the closed ring, remove the last point
    # to avoid duplication in your pipeline. (approx_polygon returns a closed ring
    # by repeating the first point at the end.)
    if len(simplified) > 2 and np.array_equal(simplified[0], simplified[-1]):
        simplified = simplified[:-1]

    return simplified


def get_border(mask: np.ndarray, s: int = 1) -> np.ndarray:
    d = dilation(mask, disk(s))
    border = (mask != d).astype(np.uint8)
    return border
