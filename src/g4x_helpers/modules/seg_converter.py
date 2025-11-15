import logging
import multiprocessing
import warnings
from collections import deque
from operator import itemgetter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from numba import njit
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance_matrix
from skimage.measure import approximate_polygon
from skimage.morphology import dilation, disk, erosion

from .. import utils
from ..g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema
from .workflow import workflow

DEFAULT_COLOR = [int(191), int(191), int(191)]


@workflow
def seg_converter(
    adata: ad.AnnData,
    seg_mask: np.ndarray,
    out_path: str | Path,
    *,
    metadata: str | Path | None = None,
    cluster_key: str | None = None,
    emb_key: str | None = None,
    protein_list: list[str] | None = None,
    n_threads: int = 4,
    logger: logging.Logger | None = None,
) -> None:
    logger.info('Creating G4X-viewer bin file.')

    warnings.filterwarnings(
        'ignore',
        message='FNV hashing is not implemented in Numba',
        category=UserWarning,
        module='numba.cpython.old_hashing',
    )

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

    ## initialize segmentation data
    ## we create polygons to define the boundaries of each cell mask
    logger.debug('Making polygons.')
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

    ## add single-cell info
    logger.debug('Adding single-cell metadata.')
    cell_ids = obs_df.index.tolist()
    num_cells = len(cell_ids)
    adata = adata[cell_ids]
    gene_names = adata.var_names
    gex = adata.X
    if not issparse(gex):
        logger.debug('Converting dense expression matrix to CSR format.')
        gex = csr_matrix(gex)
    else:
        gex = gex.tocsr()
    centroid_y = obs_df['cell_x'].tolist()
    centroid_x = obs_df['cell_y'].tolist()

    obs_df['total_counts'] = obs_df['total_counts'].astype(int)
    obs_df['n_genes_by_counts'] = obs_df['n_genes_by_counts'].astype(int)

    area_cols = obs_df.filter(like='area').columns
    if not len(area_cols):
        raise ValueError('No area column found')

    if len(area_cols) == 1:
        area_select = area_cols[0]
    else:
        logger.debug(f'Multiple area columns found: {area_cols}. Looking for expanded area column.')
        expanded_area_cols = [f for f in area_cols if 'expanded' in f]
        if len(expanded_area_cols) != 1:
            raise ValueError(
                'No expanded area column found among multiple area columns.'
                if not expanded_area_cols
                else f'Multiple expanded area columns found: {expanded_area_cols}'
            )
        area_select = expanded_area_cols[0]

    logger.debug(f'Using area column: {area_select}')

    obs_df['area'] = obs_df[area_select].astype(int)

    if emb_key:
        obs_df['umap_0'] = obs_df[f'{emb_key}_1']
        obs_df['umap_1'] = obs_df[f'{emb_key}_2']
    else:
        obs_df['umap_0'] = 0
        obs_df['umap_1'] = 0

    if protein_list:
        protein_col_pos = []
        for prot in protein_list:
            protein_col_pos.append(int(np.argwhere([x == prot for x in obs_df.columns])) + 1)
            obs_df[prot] = obs_df[prot].fillna(0).astype(int)
        get_prot_vals = itemgetter(*protein_col_pos)

    ## refine polygons
    logger.debug('Refining polygons.')

    pq_args = [
        (k, cx, cy, sorted_nonzero_values, sorted_rows, sorted_cols)
        for k, cx, cy in zip(np.arange(1, num_cells + 1), centroid_x, centroid_y)
    ]

    with multiprocessing.Pool(processes=n_threads) as pool:
        polygons = pool.starmap(refine_polygon, pq_args)

    ## do conversion
    outputCellSegmentation = CellMasksSchema.CellMasks()
    outputCellSegmentation.metadata.geneNames.extend(gene_names)
    if protein_list:
        protein_names = [x.split('_intensity')[0] for x in protein_list]
        outputCellSegmentation.metadata.proteinNames.extend(protein_names)

    ## add the cluster color map
    if clusters_available:
        obs_df['cluster_id'] = obs_df[cluster_key].astype(str)
        cluster_palette = utils.generate_cluster_palette(obs_df['cluster_id'].tolist())

        for cluster_id, color in cluster_palette.items():
            entry = CellMasksSchema.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            outputCellSegmentation.colormap.append(entry)
    else:
        entry = CellMasksSchema.ColormapEntry()
        entry.clusterId = '-1'
        entry.color.extend(DEFAULT_COLOR)
        outputCellSegmentation.colormap.append(entry)
        obs_df['cluster_id'] = '-1'

    ## add the individual cells
    for i, row in enumerate(obs_df.itertuples(index=True)):
        outputMaskData = outputCellSegmentation.cellMasks.add()
        cellPolygonPoints = [sub_coord for coord in polygons[i] for sub_coord in coord[::-1]]
        _ = outputMaskData.vertices.extend(cellPolygonPoints + cellPolygonPoints[:2])
        area = int(row.area)
        total_counts = int(row.total_counts)
        total_genes = int(row.n_genes_by_counts)
        cluster_id = str(getattr(row, 'cluster_id', '-1'))
        outputMaskData.cellId = str(row[0])
        outputMaskData.area = area
        outputMaskData.totalCounts = total_counts
        outputMaskData.totalGenes = total_genes
        outputMaskData.clusterId = cluster_id
        outputMaskData.umapValues.umapX = row.umap_0
        outputMaskData.umapValues.umapY = row.umap_1

        start = gex.indptr[i]
        end = gex.indptr[i + 1]
        indices = gex.indices[start:end]
        values = gex.data[start:end]
        _ = outputMaskData.nonzeroGeneIndices.extend(indices.tolist())
        _ = outputMaskData.nonzeroGeneValues.extend(values.astype(int).tolist())

        if protein_list:
            _ = outputMaskData.proteinValues.extend(get_prot_vals(row))

    ## write to file
    with open(out_path, 'wb') as file:
        file.write(outputCellSegmentation.SerializeToString())

    logger.debug(f'G4X-viewer bin --> {out_path}')


def bin_file_path(g4x_obj, out_dir: str | Path | None = None) -> Path:
    file_name = f'{g4x_obj.sample_id}.bin'

    if out_dir:
        out_dir = utils.validate_path(out_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)
        outfile = Path(out_dir) / file_name
    else:
        outfile = utils.create_custom_out(out_dir, 'g4x_viewer', file_name)

    return outfile


@njit
def get_start_stop_idx(arr, k):
    start_idx = np.searchsorted(arr, k, side='left')
    end_idx = np.searchsorted(arr, k, side='right')
    return start_idx, end_idx


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


def refine_polygon(k, cx, cy, sorted_nonzero_values_ref, sorted_rows_ref, sorted_cols_ref):
    start_idx, end_idx = get_start_stop_idx(sorted_nonzero_values_ref, k)
    points = np.vstack((sorted_rows_ref[start_idx:end_idx], sorted_cols_ref[start_idx:end_idx])).T
    return pointsToSingleSmoothPath(points, tolerance=2.0)


def get_border(mask: np.ndarray, s: int = 1) -> np.ndarray:
    d = dilation(mask, disk(s))
    border = (mask != d).astype(np.uint8)
    return border
