from collections import deque
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import distance_matrix
from skimage.measure import approximate_polygon
from skimage.morphology import dilation, disk

from .. import utils


def bin_file_path(g4x_obj, out_dir: str | Path | None = None) -> Path:
    file_name = f'{g4x_obj.sample_id}.bin'

    if out_dir:
        out_dir = utils.validate_path(out_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)
        outfile = Path(out_dir) / file_name
    else:
        outfile = utils.create_custom_out(out_dir, 'g4x_viewer', file_name)

    return outfile


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
    unique_clusters = [c for c in np.unique(clusters) if c != '-1']
    n_clusters = len(unique_clusters)

    if n_clusters <= 10:
        base_cmap = plt.get_cmap('tab10')
    elif n_clusters <= 20:
        base_cmap = plt.get_cmap('tab20')
    else:
        base_cmap = plt.get_cmap('hsv', min(max_colors, n_clusters))

    cluster_palette = {
        str(cluster): [int(255 * c) for c in base_cmap(i / n_clusters)[:3]] for i, cluster in enumerate(unique_clusters)
    }
    cluster_palette['-1'] = [int(191), int(191), int(191)]

    return cluster_palette


def hex2rgb(hex: str) -> list[int, int, int]:
    return [int(x * 255) for x in mcolors.to_rgb(hex)]


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
