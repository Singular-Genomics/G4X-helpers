from __future__ import annotations

from dataclasses import dataclass, field  # , asdict, field
from typing import (
    Any,
    Dict,
    Optional,
)

import g4x_helpers.plotting as g4pl
import g4x_helpers.utils as utils
import glymur
import matplotlib.pyplot as plt
import numpy as np
from g4x_helpers.dataclasses.Roi import Roi
from scipy.ndimage import generic_filter
from scipy.stats import median_abs_deviation
from skimage import exposure
from skimage.filters import threshold_otsu

glymur.set_option('lib.num_threads', 16)


@dataclass
class ImageThumbnail:
    signal_name: Optional[str] = None
    run_base: Optional[str] = None

    def __post_init__(self):
        self.img = utils._load_image_cached(self.run_base, self.signal_name, thumbnail=True)
        self.vmin, self.vmax = 0, 1


@dataclass
class ImageSignal:
    signal_name: Optional[str] = None
    run_base: Optional[str] = None

    vmin: Optional[float] = None
    vmax: Optional[float] = None

    cutoff_method: Optional[str] = 'histogram'
    cutoff_kwargs: Dict[str, Any] = field(default_factory=lambda: {'low_frac': 0.5, 'high_frac': 0.995})

    def __post_init__(self):
        self.thumb = ImageThumbnail(self.signal_name, self.run_base).img

        self.img = utils._load_image_cached(self.run_base, self.signal_name, thumbnail=False)

        if self.vmin is None and self.vmax is None:
            self.vmin, self.vmax = self.get_cutoffs(method=self.cutoff_method, **self.cutoff_kwargs)
        else:
            self.cutoff_method = 'manual'
            self.cutoff_kwargs = {}

        self.extent = utils._get_extent(self.img)
        roi_data = Roi(extent=self.extent, name='data')
        self.rois = {'data': roi_data}

    def set_cutoffs(self, vmin: float, vmax: float):
        self.vmin = vmin
        self.vmax = vmax
        self.cutoff_method = 'manual'

    def get_cutoffs(self, method: str = 'histogram', **kwargs):
        """Get cutoffs for the image using the specified method."""
        if method == 'histogram':
            return get_cutoffs_histogram(self.img, **kwargs)
        else:
            raise ValueError(f'Unknown cutoff method: {method}')

    def show(
        self,
        roi=None,
        scale=1,
        clip=True,
        scale_bar=True,
        colorbar=True,
        thumbnail=False,
        ax=None,
        size=4.5,
        return_ax=False,
    ):
        if ax is None:
            fig, ax = g4pl.montage_plot(1, panel_size=size, layout='compressed', return_fig=True)
        # else:
        #     fig = ax.figure

        ax.set_title(self.signal_name)

        if thumbnail:
            plot_image = self.thumb
            # vmax, vmin = 0, 1
        else:
            plot_image = self.img
            vmin, vmax = self.vmin, self.vmax

        if not clip or thumbnail:
            vmax = vmin = None
        elif isinstance(clip, tuple):
            vmin, vmax = clip

        roi = self.rois['data'] if roi is None else roi

        if scale != 1:
            roi = roi.scale(scale)
            crop_arr = roi._crop_array(plot_image, thumbnail=thumbnail)
        else:
            crop_arr = plot_image

        im = g4pl._show_image(crop_arr, ax, roi, vmin=vmin, vmax=vmax, cmap=g4pl.cm.lipari, scale_bar=scale_bar)

        if colorbar:
            loc = 'right' if thumbnail else 'bottom'
            g4pl._add_colorbar(ax=ax, cax=im, loc=loc)
            # fig.colorbar(im, ax=ax, orientation='vertical', pad=0.025, fraction=0.05, aspect=35)

        ax.set_xticks([])
        ax.set_yticks([])

        if return_ax:
            return ax
        else:
            plt.show()


# region later
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
