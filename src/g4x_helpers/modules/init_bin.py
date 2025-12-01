import logging
import multiprocessing
from collections import deque
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from numba import njit
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance_matrix
from skimage.measure import approximate_polygon
from skimage.morphology import dilation, disk, erosion
from tqdm import tqdm

from ..g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema
from .edit_bin import DEFAULT_COLOR, hex2rgb, update_colormap
from .segment import get_cell_ids
from .workflow import OutSchema, workflow

if TYPE_CHECKING:
    from ..models import G4Xoutput


@workflow
def init_bin_file(
    g4x_obj: 'G4Xoutput',
    out_dir: str | Path,
    seg_mask: np.ndarray | None = None,
    *,
    n_threads: int = 4,
    logger: logging.Logger | None = None,
) -> None:
    out_tree = OutSchema(out_dir, subdirs=['g4x_viewer'])
    out_file = out_tree.g4x_viewer / f'{g4x_obj.sample_id}_segmentation.bin'

    adata = g4x_obj.load_adata(load_clustering=True, remove_nontargeting=False)

    if seg_mask is None:
        logger.info('Loading default segmentation mask.')
        seg_mask = g4x_obj.load_segmentation()

    logger.info('Loading clustering information.')
    obs_df = prepare_metadata(sample_id=g4x_obj.sample_id, adata=adata, seg_mask=seg_mask)
    print(f'Number of cells to process: {obs_df.shape[0]}')

    logger.info('Making polygons.')
    cell_ids, pq_args = initialize_segmentation_data(obs_df, seg_mask)

    logger.debug('Adding single-cell metadata.')
    obs_df, gene_names, gex = add_singlecell_info(obs_df, adata, cell_ids)

    ## do conversion
    cell_masks = CellMasksSchema.CellMasks()
    cell_masks.metadata.geneNames.extend(gene_names)

    logger.info('Setting default color.')
    cluster_palette = {}
    cluster_palette['unassigned'] = hex2rgb(DEFAULT_COLOR)

    update_colormap(cell_masks, cluster_palette)

    protein_values = None
    if g4x_obj.includes_protein:
        protein_names = g4x_obj.proteins
        protein_list = [prot + '_intensity_mean' for prot in protein_names]

        cell_masks.metadata.proteinNames.extend(protein_names)

        obs_df = obs_df.with_columns([pl.col(prot).fill_null(0).cast(pl.Int64) for prot in protein_list])
        protein_values = obs_df.select(protein_list).to_numpy()

    logger.info('Refining polygons.')
    with multiprocessing.Pool(processes=n_threads) as pool:
        polygons = pool.starmap(refine_polygon, pq_args)

    logger.info('Adding individual cells.')
    for i, row in enumerate(tqdm(obs_df.iter_rows(named=True), total=obs_df.height, desc='Adding individual cells.')):
        output_mask_data = cell_masks.cellMasks.add()

        cell_polygon_pts = [sub_coord for coord in polygons[i] for sub_coord in coord[::-1]]
        output_mask_data.vertices.extend(cell_polygon_pts + cell_polygon_pts[:2])

        output_mask_data.cellId = str(row['cell_id'])
        output_mask_data.area = int(row['area'])
        output_mask_data.totalCounts = int(row['total_counts'])
        output_mask_data.totalGenes = int(row['n_genes_by_counts'])
        output_mask_data.clusterId = str(row['cluster_id'])
        output_mask_data.umapValues.umapX = float('nan')
        output_mask_data.umapValues.umapY = float('nan')

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
        file.write(cell_masks.SerializeToString())

    logger.debug(f'G4X-viewer bin --> {out_file}')


def prepare_metadata(sample_id: str, adata, seg_mask: np.ndarray) -> pl.DataFrame:
    obs_df = get_cell_ids(sample_id=sample_id, mask=seg_mask)

    metadata_df = pl.from_pandas(adata.obs)
    obs_df = obs_df.join(metadata_df, left_on='label', right_on='cell_id', how='left').rename({'label': 'cell_id'})

    obs_df = obs_df.with_columns(pl.lit('unassigned').alias('cluster_id'))

    return obs_df


def initialize_segmentation_data(obs_df: pl.DataFrame, seg_mask: np.ndarray) -> tuple[list[str], list[tuple]]:
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


def add_singlecell_info(obs_df: pl.DataFrame, adata, cell_ids: list[str]):
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

    return obs_df, gene_names, gex


def refine_polygon(
    cell_id: int, sorted_nonzero_values_ref: np.ndarray, sorted_rows_ref: np.ndarray, sorted_cols_ref: np.ndarray
) -> np.ndarray:
    start_idx, end_idx = get_start_stop_idx(sorted_nonzero_values_ref, cell_id)
    points = np.vstack((sorted_rows_ref[start_idx:end_idx], sorted_cols_ref[start_idx:end_idx])).T
    return pointsToSingleSmoothPath(points, tolerance=2.0)


@njit
def get_start_stop_idx(arr: np.ndarray, k: int) -> tuple[int, int]:
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


def computeLongestPath(adj_list: dict[int, list[int]]) -> list[int]:
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


def returnEndpoints(adj_list: dict[int, list[int]], adjacency: int = 2) -> list[int]:
    # Identify endpoints of the MST
    endpoints = [node for node in adj_list if len(adj_list[node]) == adjacency]

    return endpoints


def bfs_path(start: int, end: int, adj_list: dict[int, list[int]]) -> list[int]:
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
def createAdjacencyList_numba(mst: np.ndarray) -> tuple:
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


def indicesToArray(points: np.ndarray, longest_path: list[int]) -> np.ndarray:
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
