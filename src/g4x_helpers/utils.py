from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    List,
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
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
from scipy.ndimage import generic_filter
from scipy.stats import median_abs_deviation
from shapely.affinity import scale, translate
from shapely.geometry import Polygon
from skimage import exposure
from skimage.filters import threshold_otsu


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


class Roi:
    def __init__(
        self,
        xlims: Optional[Tuple[int, int]] = None,
        ylims: Optional[Tuple[int, int]] = None,
        extent: Optional[List[int]] = None,
        center: Optional[Tuple[int, int]] = None,
        edge_size: Optional[int] = None,
        polygon: Optional[Polygon] = None,
        name: str = 'ROI',
        sub_rois: Optional[int] = None,
    ):
        self.name = name

        if center and edge_size and not polygon and not (xlims and ylims) and not extent:
            center = (center[0], center[1])
            half_edge = edge_size / 2
            self.xlims = (center[0] - half_edge, center[0] + half_edge)
            self.ylims = (center[1] - half_edge, center[1] + half_edge)

            self.polygon = self._make_polygon()

        if not polygon and (xlims and ylims) or extent:
            xlims = xlims if xlims else extent[0:2]
            ylims = ylims if ylims else extent[2:4]

            self.xlims = xlims
            self.ylims = ylims

            self.polygon = self._make_polygon()

        elif polygon and not (xlims and ylims):
            self.polygon = polygon
            polygon_bounds = self.polygon.bounds
            self.xlims = (polygon_bounds[0], polygon_bounds[2])
            self.ylims = (polygon_bounds[1], polygon_bounds[3])

        # Compute overall width and height
        self.width = self.xlims[1] - self.xlims[0]
        self.height = self.ylims[1] - self.ylims[0]
        self.extent = (self.xlims[0], self.xlims[1], self.ylims[0], self.ylims[1])
        self.extent_array = np.array(self.extent).astype(np.int32)
        self.lims = (self.xlims, self.ylims)
        self.order = 'xy'

        if sub_rois:
            self.sub_rois = self.subtile_roi(n=sub_rois, labels='abc')

    def _make_polygon(self):
        return Polygon(
            [
                (self.xlims[0], self.ylims[0]),
                (self.xlims[1], self.ylims[0]),
                (self.xlims[1], self.ylims[1]),
                (self.xlims[0], self.ylims[1]),
            ]
        )

    def scale(self, factor):
        scaled_roi = scale(self.polygon, xfact=factor, yfact=factor, origin='center')
        sc_roi = Roi(polygon=scaled_roi, name=None)
        return sc_roi

    def translate(self, xoff=0, yoff=0):
        trans_roi = translate(self.polygon, xoff=xoff, yoff=yoff)
        tr_roi = Roi(polygon=trans_roi, name=None)
        return tr_roi

    def subtile_roi(self, n=3, labels='123'):
        scale_factor = 1 / n  # Each sub-ROI should have width and height 1/n of the original

        scaled_roi = self.scale(scale_factor)

        sub_rois = []
        label = 1  # Initialize counter so the top-left sub-ROI becomes 1

        # Reverse the vertical loop: This ensures we start with the top row
        for j in reversed(range(n)):
            for i in range(n):
                # Compute the translation offsets relative to the overall ROI center.
                offset_x = (i - (n - 1) / 2) * (self.width / n)
                offset_y = (j - (n - 1) / 2) * (self.height / n)

                if labels == 'abc':
                    roi_label = chr(65 + (label - 1))
                elif labels == '123':
                    roi_label = str(label)

                sub_roi = scaled_roi.translate(xoff=offset_x, yoff=offset_y)
                sub_roi.name = roi_label
                sub_rois.append(sub_roi)

                label += 1

        return sub_rois

    def _crop_array(self, array, thumbnail=False):
        fct = 1 if not thumbnail else 5

        xlim = (self.extent_array[0:2] / fct).astype(int)
        ylim = (self.extent_array[2:4] / fct).astype(int)

        return array[ylim[0] : ylim[1], xlim[0] : xlim[1]]

    def _limit_ax(self, ax):
        ax.set_xlim(self.xlims)
        ax.set_ylim(self.ylims)

    def add_to_plot(
        self,
        ax,
        color='white',
        lw=2,
        label=True,
        label_position='top_left',
        prefix='ROI ',
        font_size=12,
        pad=0.05,
        fontweight='bold',
    ):
        import matplotlib.patches as patches

        xvals, yvals = self.lims

        x, y = xvals[0], yvals[0]
        w, h = xvals[1] - x, yvals[1] - y

        rect = patches.Rectangle((x, y), w, h, linewidth=lw, edgecolor=color, facecolor='none', zorder=10)
        ax.add_patch(rect)

        # If label is True, add a label at the desired location using relative padding.
        if label:
            if hasattr(self, 'name'):
                label_text = self.name
                if label_text == 'A':
                    label_position = 'top_left'
                elif label_text == 'B':
                    label_position = 'top_right'
                elif label_text == 'C':
                    label_position = 'bottom_left'
                elif label_text == 'D':
                    label_position = 'bottom_right'

            else:
                label_text = f'{prefix}'

            # Calculate relative padding in data units.
            pad_x = pad * w
            pad_y = pad * h

            position_map = {
                'top_left': (lambda x, y, w, h, pad_x, pad_y: (x + pad_x, y + h - pad_y), 'left', 'top'),
                'top_right': (lambda x, y, w, h, pad_x, pad_y: (x + w - pad_x, y + h - pad_y), 'right', 'top'),
                'bottom_left': (lambda x, y, w, h, pad_x, pad_y: (x + pad_x, y + pad_y), 'left', 'bottom'),
                'bottom_right': (lambda x, y, w, h, pad_x, pad_y: (x + w - pad_x, y + pad_y), 'right', 'bottom'),
                'top_center': (lambda x, y, w, h, pad_x, pad_y: (x + w / 2, y + h - pad_y), 'center', 'top'),
                'bottom_center': (lambda x, y, w, h, pad_x, pad_y: (x + w / 2, y + pad_y), 'center', 'bottom'),
                'center': (lambda x, y, w, h, pad_x, pad_y: (x + w / 2, y + h / 2), 'center', 'center'),
            }

            def get_label_position(label_position, x, y, w, h, pad_x, pad_y):
                entry = position_map.get(label_position)
                if entry is None:
                    print('Invalid label_position. Use one of', list(position_map.keys()))
                    return
                position_fn, ha, va = entry
                tx, ty = position_fn(x, y, w, h, pad_x, pad_y)
                return tx, ty, ha, va

            tx, ty, ha, va = get_label_position(label_position, x, y, w, h, pad_x, pad_y)

            # Add the label text to the axis.
            ax.text(
                tx,
                ty,
                label_text,
                color=color,
                fontsize=font_size,
                fontweight=fontweight,
                verticalalignment=va,
                horizontalalignment=ha,
                zorder=11,
            )


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


def get_cutoffs_quantile(image, lower=0.05, upper=0.995):
    vmin, vmax = np.quantile(image, lower), np.quantile(image, upper)
    return vmin, vmax


def get_cutoffs_min_ratio(image, lower=0.05, upper=0.995):
    vmax = np.quantile(image, upper)
    vmin = vmax * lower
    return vmin, vmax


def get_cutoffs_otsu(image, upper_quantile=0.995):
    vmin = threshold_otsu(image)
    vmax = np.quantile(image, upper_quantile)
    return vmin, vmax


def get_cutoffs_robust(image, z_thresh=3.5):
    flat = image.flatten()
    med = np.median(flat)
    mad = median_abs_deviation(flat, scale='normal')
    z_scores = (flat - med) / mad
    filtered = flat[(z_scores > -z_thresh) & (z_scores < z_thresh)]
    return np.min(filtered), np.max(filtered)


def get_cutoffs_histogram(image, low_frac=0.001, high_frac=0.995):
    hist, bin_edges = exposure.histogram(image)
    cdf = np.cumsum(hist) / hist.sum()
    vmin = bin_edges[np.searchsorted(cdf, low_frac)]
    vmax = bin_edges[np.searchsorted(cdf, high_frac)]
    return vmin, vmax


def get_cutoffs_local_percentile(image, min_percentile=1, max_percentile=99, size=51):
    def local_percentile(image, percentile, size=51):
        def func(x):
            return np.percentile(x, percentile)

        return generic_filter(image, func, size=(size, size))

    vmin_map = local_percentile(image, min_percentile, size)  # local low-end, could be noise
    vmax_map = local_percentile(image, max_percentile, size)  # local high-end, reduces outliers

    # clipped = np.clip(image, vmin_map, vmax_map)

    vmin = np.mean(vmin_map)
    vmax = np.mean(vmax_map)
    return vmin, vmax


def apply_clahe(image, kernel_size=(50, 50), clip_limit=0.01):
    return exposure.equalize_adapthist(image, kernel_size=kernel_size, clip_limit=clip_limit)
