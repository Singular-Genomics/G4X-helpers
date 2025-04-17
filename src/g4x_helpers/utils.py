import logging
import os
from pathlib import Path
from typing import Any, Callable, Generator, Iterable, Iterator, List, Literal, Optional, Tuple, Union

import numpy as np
from shapely.affinity import scale, translate
from shapely.geometry import Polygon


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


def downsample_mask(mask, mask_length=None, ax=None, display_size=None, downsample='auto'):
    def get_ax_size(ax):
        fig = ax.get_figure()

        bbox = ax.get_position()
        fig_width, fig_height = fig.get_size_inches()
        width_in = bbox.width * fig_width
        height_in = bbox.height * fig_height

        ax_dim = (width_in, height_in)
        return ax_dim

    if downsample is None:
        return mask

    elif isinstance(downsample, int):
        downsample_fct = int(downsample)

    elif downsample.startswith('auto'):
        if downsample == 'auto':
            target_ppi = 200
        elif downsample.startswith('auto_'):
            target_ppi = int(downsample.split('_')[-1])

        if ax is None and display_size is None:
            raise ValueError("Either ax or display_size must be provided when downsample is 'auto'")
        elif ax is None and display_size:
            panel_size = display_size
        elif ax and display_size is None:
            ax_dim = get_ax_size(ax)
            panel_size = np.array(ax_dim).max()

        if mask_length is None:
            mask_length = np.array(mask.shape).max()

        downsample_fct = int((mask_length / panel_size) / target_ppi)
        # print(f'Downsample factor: {downsample_fct}')

    else:
        raise ValueError("downsample must be None, int or 'auto'")

    if downsample_fct > 0:
        sampled_mask = mask[::downsample_fct, ::downsample_fct].copy()
    else:
        sampled_mask = mask.copy()
    return sampled_mask


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


class Roi:
    def __init__(self, xlims=None, ylims=None, extent=None, center=None, edge_size=None, polygon=None, name='ROI', sub_rois=None):
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
        self.lims = (self.xlims, self.ylims)
        self.order = 'xy'

        if sub_rois:
            self.sub_rois = self.subtile_roi(n=sub_rois, labels='abc')

    def _make_polygon(self):
        return Polygon([(self.xlims[0], self.ylims[0]), (self.xlims[1], self.ylims[0]), (self.xlims[1], self.ylims[1]), (self.xlims[0], self.ylims[1])])

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


def add_rois(ax, rois=list, color='white', lw=2, order='xy', label=True, label_position='top_left', prefix='ROI ', font_size=12, pad=0.05):
    import matplotlib.patches as patches

    if not isinstance(rois, list):
        rois = [rois]

    for idx, roi in enumerate(rois):
        if hasattr(roi, 'lims'):
            xvals, yvals = roi.lims
        else:
            if order == 'xy':
                xvals, yvals = roi[0], roi[1]
            elif order == 'yx':
                yvals, xvals = roi[0], roi[1]
            elif order not in ['xy', 'yx']:
                print("order must be specified as 'xy' or 'yx'")
                return

        x, y = xvals[0], yvals[0]
        w, h = xvals[1] - x, yvals[1] - y

        rect = patches.Rectangle((x, y), w, h, linewidth=lw, edgecolor=color, facecolor='none', zorder=10)
        ax.add_patch(rect)

        # If label is True, add a label at the desired location using relative padding.
        if label:
            if hasattr(roi, 'name'):
                label_text = roi.name
                if label_text == 'A':
                    label_position = 'top_left'
                elif label_text == 'B':
                    label_position = 'top_right'
                elif label_text == 'C':
                    label_position = 'bottom_left'
                elif label_text == 'D':
                    label_position = 'bottom_right'

            else:
                label_text = f'{prefix}{idx + 1}'
            # Calculate relative padding in data units.
            pad_x = pad * w
            pad_y = pad * h

            # Determine text coordinates and alignments based on label_position.
            if label_position == 'top_left':
                tx, ty = x + pad_x, y + h - pad_y
                ha, va = 'left', 'top'
            elif label_position == 'top_right':
                tx, ty = x + w - pad_x, y + h - pad_y
                ha, va = 'right', 'top'
            elif label_position == 'bottom_left':
                tx, ty = x + pad_x, y + pad_y
                ha, va = 'left', 'bottom'
            elif label_position == 'bottom_right':
                tx, ty = x + w - pad_x, y + pad_y
                ha, va = 'right', 'bottom'
            elif label_position == 'top_center':
                tx, ty = x + w / 2, y + h - pad_y
                ha, va = 'center', 'top'
            elif label_position == 'bottom_center':
                tx, ty = x + w / 2, y + pad_y
                ha, va = 'center', 'bottom'
            elif label_position == 'center':
                tx, ty = x + w / 2, y + h / 2
                ha, va = 'center', 'center'
            else:
                print("Invalid label_position. Use one of ['top_left', 'top_right', 'bottom_left', 'bottom_right', 'top_center', 'bottom_center', 'center'].")
                return

            # Add the label text to the axis.
            ax.text(tx, ty, label_text, color=color, fontsize=font_size, fontweight='bold', verticalalignment=va, horizontalalignment=ha, zorder=11)
