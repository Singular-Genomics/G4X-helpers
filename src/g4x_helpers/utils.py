from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Optional,
    Tuple,
    Union,
)  # , Any, Callable, Generator, Iterable, Iterator, Literal

if TYPE_CHECKING:
    # This import is only for type checkers (mypy, PyCharm, etc.), not at runtime
    from g4x_helpers.models import G4Xoutput

import logging
import os
import shutil
from functools import lru_cache
from pathlib import Path

import glymur
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


def setup_logger(
    logger_name: str,
    *,
    stream_logger: Optional[bool] = True,
    stream_level: Optional[int] = logging.DEBUG,
    file_logger: Optional[bool] = False,
    file_level: Optional[int] = logging.DEBUG,
    file_mode: Optional[str] = 'a',
    out_dir: Optional[Union[Path, str]] = None,
    format: Optional[str] = None,
    clear_handlers: Optional[bool] = False,
) -> logging.Logger:
    """
    Sets up a logger with configurable stream and file handlers.

    Parameters
    ----------
    logger_name : str
        Name of the logger.
    stream_logger : bool, optional
        Whether to enable logging to the console (default is True).
    stream_level : int, optional
        Logging level for the stream handler (default is logging.DEBUG).
    file_logger : bool, optional
        Whether to enable logging to a file (default is False).
    file_level : int, optional
        Logging level for the file handler (default is logging.DEBUG).
    file_mode : str, optional
        File mode for the log file. Common options: 'a' for append, 'w' for overwrite (default is 'a').
    out_dir : Path or str, optional
        Directory where the log file will be saved. Required if `file_logger` is True.
    format : str, optional
        Custom log message format. If not provided, a default format will be used.
    clear_handlers : bool, optional
        Whether to clear existing handlers from the logger before adding new ones (default is False).

    Returns
    -------
    logging.Logger
        Configured logger instance.

    """

    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    if format is None:
        format = '[%(asctime)s %(levelname)s %(funcName)s %(name)s] %(message)s'
    formatter = logging.Formatter(format)

    ## optionally clear existing handlers
    if clear_handlers:
        logger.handlers.clear()

    if stream_logger:
        h = logging.StreamHandler()
        h.setLevel(stream_level)
        h.setFormatter(formatter)
        logger.addHandler(h)

    if file_logger:
        assert out_dir is not None, 'out_dir must be provided if file_logger is True'
        os.makedirs(out_dir, exist_ok=True)
        fh = logging.FileHandler(f'{out_dir}/{logger_name}.log', mode=file_mode)
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def find_tissue(coords, ubins=100, expand=0.1, order='yx', threshold=0.1, smp_shape=(19200, 15232)):
    roi = []
    idx_order = {'yx': (0, 1), 'xy': (1, 0)}

    # Histogram analysis to determine ROI bounds in each dimension
    for idx in idx_order[order]:
        data = coords[:, idx]
        # Determine number of bins based on a scaling of the maximum coordinate and ubins
        coord_bins = int(data.max() * 0.3125 / ubins)
        counts, bin_edges = np.histogram(data, bins=coord_bins)
        cells = sum(counts)
        deltas = []
        for i in range(len(counts)):
            # Compute percent change between adjacent bins
            d = abs(counts[i] - counts[i - 1])
            mag = (d / cells) * 100
            deltas.append(mag)
        # Find the first bin from the start where change exceeds the threshold
        for i, d in enumerate(deltas):
            if d > threshold:
                low_i = i
                break
        # Find the first bin from the end where change exceeds the threshold
        for i, d in enumerate(reversed(deltas)):
            if d > threshold:
                high_i = len(deltas) - i
                break
        roi.append((int(bin_edges[low_i]), int(bin_edges[high_i])))

    # Use the two ROI boundaries to compute an expanded square ROI
    edge_1 = roi[0][1] - roi[0][0]
    edge_2 = roi[1][1] - roi[1][0]
    edge = max(edge_1, edge_2)
    edge = int(edge * (1 + expand))

    # Compute midpoints for each dimension
    dim_1_mid = int((roi[0][1] + roi[0][0]) / 2)
    dim_2_mid = int((roi[1][1] + roi[1][0]) / 2)

    # Check if the ROI is too large to fit in the allowed boundaries.
    allowed_1 = smp_shape[0]
    allowed_2 = smp_shape[1]

    # Determine scaling factor in each dimension and use the smaller one.
    scale_h = allowed_1 / edge
    scale_w = allowed_2 / edge
    scale_factor = min(scale_w, scale_h, 1)

    if scale_factor < 1:
        edge = int(edge * scale_factor)

    # Recalculate the ROI centered on the computed centroid using the (possibly scaled) edge.
    new_roi_1_min = int(dim_1_mid - edge / 2)
    new_roi_1_max = new_roi_1_min + edge

    new_roi_2_min = int(dim_2_mid - edge / 2)
    new_roi_2_max = new_roi_2_min + edge

    # Shift the ROI if it lies outside the allowed boundaries.
    # For the y dimension:
    if new_roi_1_min < 0:
        new_roi_1_min = 0
        new_roi_1_max = new_roi_1_min + edge
    elif new_roi_1_max > allowed_1:
        new_roi_1_max = allowed_1
        new_roi_1_min = new_roi_1_max - edge

    # For the x dimension:
    if new_roi_2_min < 0:
        new_roi_2_min = 0
        new_roi_2_max = new_roi_2_min + edge
    elif new_roi_2_max > allowed_2:
        new_roi_2_max = allowed_2
        new_roi_2_min = new_roi_2_max - edge

    final_lims = ((new_roi_1_min, new_roi_1_max), (new_roi_2_min, new_roi_2_max))
    final_centroid = ((new_roi_1_min + new_roi_1_max) // 2, (new_roi_2_min + new_roi_2_max) // 2)

    return {'lims': final_lims, 'centroid': final_centroid}


@lru_cache(maxsize=None)
def _load_image_cached(run_base: str, signal: str, thumbnail: bool) -> Tuple[np.ndarray, float, float]:
    """This is cached per (run_base, signal, thumbnail)."""
    base = Path(run_base)

    folder = 'h_and_e' if signal in ('h_and_e', 'nuclear', 'eosin') else 'protein'

    p = base / folder
    suffix = '.jpg' if (thumbnail and signal == 'h_and_e') else ('.png' if thumbnail else '.jp2')

    fname = f'{signal}_thumbnail{suffix}' if thumbnail else f'{signal}.jp2'
    img_file = p / fname

    if img_file.suffix == '.png':
        img = plt.imread(img_file)
    else:
        img = glymur.Jp2k(str(img_file))[:]

    return img


def _get_extent(img):
    shp = img.shape
    return (0, shp[1], shp[0], 0)


def _create_custom_out(sample: 'G4Xoutput', out_dir=None, parent_folder=None, file_name=None):
    if out_dir is None:
        custom_out = sample.run_base / parent_folder / 'custom'
    else:
        custom_out = Path(out_dir) / parent_folder

    if not custom_out.exists():
        os.makedirs(custom_out)
    else:
        shutil.rmtree(custom_out)
        os.makedirs(custom_out)

    outfile = custom_out / file_name
    return outfile


def _normalize(array):
    min_val = array.min()
    max_val = array.max()
    return (array - min_val) / (max_val - min_val) if max_val > min_val else array


def _apply_thresholds(array, tmin, tmax):
    array = np.clip(array, a_min=tmin, a_max=tmax)
    array_normed = _normalize(array)

    return array_normed


def _additive_blend(prot_arrays, marker_colors, clip_mode='hard'):
    rgb_images = []
    for i in range(len(prot_arrays)):
        rgb_images.append(np.stack([prot_arrays[i]] * 3, axis=-1) * np.array(marker_colors[i]))

    blended_image = np.zeros_like(rgb_images[0])

    for img in rgb_images:
        blended_image += img

    if clip_mode == 'hard':
        blended_image = np.clip(blended_image, 0, 255).astype(np.uint8)
    elif clip_mode == 'soft':
        blended_image = _normalize(blended_image) * 255
        blended_image = blended_image.astype(np.uint8)
    elif clip_mode is None:
        pass

    return blended_image


def blend_protein_images(sample, roi=None, query=None, marker_colors=None, cutoffs='global'):
    from .dataclasses.Signal import get_cutoffs_histogram

    # ___> gather images
    signals = []
    for q in query:
        img, vmin, vmax = sample.load_image(q, thumbnail=False, return_th=True)
        pinfo = {'name': q, 'img': img, 'vmin': vmin, 'vmax': vmax}

        signals.append(pinfo)

    # ___> process images
    processed_arrays = {}
    for protein, color in zip(signals, marker_colors):
        if roi is None:
            cropped_array = protein['img']
        else:
            cropped_array = roi._crop_array(protein['img'])

        if cutoffs == 'global':
            vmin, vmax = protein['vmin'], protein['vmax']
        elif cutoffs == 'local':
            vmin, vmax = get_cutoffs_histogram(cropped_array, low_frac=0.5, high_frac=0.995)
        else:
            vmin, vmax = cutoffs

        th_array = _apply_thresholds(cropped_array, vmin, vmax)
        processed_arrays[protein['name']] = {'img': th_array, 'vmin': vmin, 'vmax': vmax, 'color': color}

    arrays = [value['img'] for key, value in processed_arrays.items() if 'img' in value]
    protein_colors = [np.array(mcolors.to_rgb(c)) * 255 for c in marker_colors]
    blended_image = _additive_blend(arrays, protein_colors, clip_mode='hard')

    for key in processed_arrays:
        if 'img' in processed_arrays[key]:
            del processed_arrays[key]['img']

    processed_arrays['blended'] = blended_image

    return processed_arrays
