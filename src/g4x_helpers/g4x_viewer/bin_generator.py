from collections import deque
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ray
import random
from numba import njit
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance_matrix
from skimage.measure import approximate_polygon
from skimage.morphology import dilation, disk, erosion
from tqdm import tqdm
from anndata import AnnData
from . import CellMasksSchema_pb2 as CellMasksSchema


def generate_cluster_palette(clusters: list, max_colors: int = 256) -> dict:
    unique_clusters = [c for c in np.unique(clusters) if c != '-1']
    n_clusters = len(unique_clusters)

    if n_clusters <= 10:
        base_cmap = plt.get_cmap('tab10')
    elif n_clusters <= 20:
        base_cmap = plt.get_cmap('tab20')
    else:
        base_cmap = plt.get_cmap('hsv', min(max_colors, n_clusters))

    cluster_palette = {
        str(cluster): [int(255 * c) for c in base_cmap(i / n_clusters)[:3]]
        for i, cluster in enumerate(unique_clusters)
    }
    cluster_palette["-1"] = [int(191), int(191), int(191)]

    return cluster_palette


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


@ray.remote(num_cpus=1)
def refine_polygon(k, cx, cy, sorted_nonzero_values_ref, sorted_rows_ref, sorted_cols_ref):
    start_idx, end_idx = get_start_stop_idx(sorted_nonzero_values_ref, k)
    points = np.vstack((sorted_rows_ref[start_idx:end_idx], sorted_cols_ref[start_idx:end_idx])).T
    return pointsToSingleSmoothPath(points, tolerance=2.0)


def get_border(mask: np.ndarray, s: int = 1) -> np.ndarray:
    d = dilation(mask, disk(s))
    border = (mask != d).astype(np.uint8)
    return border


def seg_converter(
    adata: AnnData,
    seg_mask: np.ndarray,
    outpath: str | Path,
    *,
    cluster_info: str| Path | None = None,
    cluster_key: str | None = None,
    emb_key: str | None = None,
    protein_list: list[str] | None = None,
    n_threads: int = 4
) -> None:

    if cluster_info is None:
        clusters_available = False
        obs_df = adata.obs.copy()
    else:
        clusters_available = True
        clustered_df = pd.read_csv(cluster_info, index_col= 0, header= 0)
        if clustered_df.shape[1] > 1:
            assert cluster_key is not None, "ERROR: multiple columns detected in cluster_info, cluster_key must be provided."
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
    cell_ids = obs_df.index.tolist()
    num_cells = len(cell_ids)
    centroid_y = obs_df['cell_x'].tolist()
    centroid_x = obs_df['cell_y'].tolist()
    areas = obs_df['nuclei_expanded_area'].tolist()
    total_counts = obs_df['total_counts'].tolist()
    total_genes = obs_df['n_genes_by_counts'].tolist()
    if clusters_available:
        clusters = obs_df[cluster_key].tolist()
        cluster_palette = generate_cluster_palette(clusters)
        cluster_colors = obs_df[cluster_key].astype(str).map(cluster_palette).tolist()
    else:
        cmap = plt.get_cmap('hsv', 100)
        clusters = ["-1"] * num_cells
        cluster_colors = [[int(255 * x) for x in cmap(i % 100)[:3]] for i in range(num_cells)]
        random.shuffle(cluster_colors)
    if protein_list:
        prot_vals = list(obs_df[protein_list].to_dict(orient= 'index').values())
    else:
        prot_vals = [{} for _ in range(num_cells)]
    if emb_key:
        umap_x = obs_df[f"{emb_key}_0"].to_numpy()
        umap_y = obs_df[f"{emb_key}_1"].to_numpy()
    else:
        umap_x = np.zeros(num_cells)
        umap_y = np.zeros(num_cells)

    ## refine polygons
    ray.init(
        num_cpus=n_threads,
        logging_level='WARNING',
        include_dashboard=False,
        address='local',
        _node_ip_address='127.0.0.1'
    )

    sorted_nonzero_values_ref = ray.put(sorted_nonzero_values)
    sorted_rows_ref = ray.put(sorted_rows)
    sorted_cols_ref = ray.put(sorted_cols)
    pq_args = [
        (k, cx, cy, sorted_nonzero_values_ref, sorted_rows_ref, sorted_cols_ref)
        for k, cx, cy in zip(np.arange(start=1, stop=len(cell_ids) + 1), centroid_x, centroid_y)
    ]
    polygons = ray.get([refine_polygon.remote(*args) for args in pq_args])

    ray.shutdown()

    ## do conversion
    # logger.debug(f"{sample_id}: Converting data to protobuff format...")
    segmentation_source = {
        'xs': [[y[0] for y in x] for x in polygons],
        'ys': [[y[1] for y in x] for x in polygons],
        'colors': cluster_colors,
        'cell_id': cell_ids,
        'area': areas,
        'centroid_x': centroid_x,
        'centroid_y': centroid_y,
        'total_tx': total_counts,
        'total_genes': total_genes,
        'cluster_id': clusters
    }
    if protein_list:
        segmentation_source['protein'] = prot_vals
    if emb_key:
        segmentation_source['umap_x'] = umap_x
        segmentation_source['umap_y'] = umap_y

    outputCellSegmentation = CellMasksSchema.CellMasks()

    for index in tqdm(range(len(segmentation_source['cell_id']))):
        try:
            cellPolygonPoints = [
                coord
                for pair in zip(segmentation_source['ys'][index], segmentation_source['xs'][index])
                for coord in pair
            ]
            cellPolygonColor = segmentation_source['colors'][index]
            cellId = segmentation_source['cell_id'][index]
            cellTotalCounts = segmentation_source['total_tx'][index]
            cellTotalGenes = segmentation_source['total_genes'][index]
            cellArea = segmentation_source['area'][index]
            clusterId = segmentation_source['cluster_id'][index]
            cellProt = segmentation_source["protein"][index]
            cellUmapX = segmentation_source["umap_x"][index]
            cellUmapY = segmentation_source["umap_y"][index]
        except Exception as e:
            # logger.exception(e)
            pass
        outputMaskData = outputCellSegmentation.cellMasks.add()
        outputMaskData.vertices.extend(cellPolygonPoints + cellPolygonPoints[:2])
        outputMaskData.color.extend(cellPolygonColor)
        outputMaskData.cellId = str(cellId)
        outputMaskData.area = str(cellArea)
        outputMaskData.totalCounts = str(cellTotalCounts)
        outputMaskData.totalGenes = str(cellTotalGenes)
        outputMaskData.clusterId = str(clusterId)
        outputMaskData.umapValues.umapX = cellUmapX
        outputMaskData.umapValues.umapY = cellUmapY
        outputMaskData.proteins.update(cellProt)

    if clusters_available:
        for cluster_id, color in cluster_palette.items():
            entry = CellMasksSchema.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            outputCellSegmentation.colormap.append(entry)
    else:
        entry = CellMasksSchema.ColormapEntry()
        entry.clusterId = '-1'
        entry.color.extend([int(31), int(119), int(180)])
        outputCellSegmentation.colormap.append(entry)   

    ## write to file
    with open(outpath, 'wb') as file:
        file.write(outputCellSegmentation.SerializeToString())
