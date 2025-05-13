from __future__ import annotations

import warnings
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
from scipy.stats import kurtosis, median_abs_deviation, skew
from skimage import exposure, filters
from skimage.morphology import disk, opening, remove_small_objects

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
        self.shape = self.img.shape

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

    def provide_flat_region_mask(self, flat_region_mask):
        assert flat_region_mask.shape == self.shape

        self.flat_region_mask = flat_region_mask

    def calculate_signal_stats(self, pp_img=None, noise_stats=None, clean_mask=False, return_mask=False):
        if pp_img is None and noise_stats is None:
            # print('No pp_img or noise_values provided, using _pp_snr()')
            pp_img, noise_stats = self._pp_snr()

        noise_values, mu_n, med_n, sigma_n = noise_stats.values()

        report = {}

        signal_mask, thresh = mask_signal_ostu(pp_img)

        if thresh < 3 * sigma_n:
            signal_mask = np.zeros_like(pp_img, dtype=bool)

        if clean_mask:
            signal_mask = _clean_signal_mask(signal_mask, opening_radius=2, min_object_size=64)

        signal_values = pp_img[signal_mask]
        backgrnd_values = pp_img[~signal_mask]

        def _get_basic_stats_safe(values):
            """
            If values is empty, return (nan, nan, nan),
            otherwise delegate to your existing _get_basic_stats(),
            suppressing any RuntimeWarning or SmallSampleWarning.
            """
            if values.size == 0:
                return np.nan, np.nan, np.nan
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # silences RuntimeWarning, SmallSampleWarning, etc.
                return _get_basic_stats(values)

        def safe_div(numer, denom, default=np.nan, eps=1e-12):
            try:
                # if denom is a scalar
                return numer / denom if abs(denom) > eps else default
            except TypeError:
                # denom might be an array
                safe_d = np.where(np.abs(denom) > eps, denom, np.nan)
                return numer / safe_d

        def shannon_entropy(values, bins=256):
            """
            Shannon entropy of intensity distribution, but returns nan if empty.
            """
            if values.size == 0:
                return np.nan
            flat = values.ravel()
            hist, _ = np.histogram(flat, bins=bins, density=True)
            hist = hist[hist > 0]
            return -np.sum(hist * np.log2(hist))

        mu_s, med_s, sigma_s = _get_basic_stats_safe(signal_values)
        mu_b, med_b, sigma_b = _get_basic_stats_safe(backgrnd_values)

        bit_depth = 12
        I_max = 2**bit_depth - 1
        psnr = 10 * np.log10((I_max**2) / (sigma_n**2))

        # standard metrics
        snr = safe_div(mu_s, sigma_n)
        sbr = safe_div(mu_s, mu_b)
        snr_bs = safe_div((mu_s - mu_b), sigma_n)
        cnr = safe_div((mu_s - mu_b), sigma_b)

        rms = np.sqrt(((pp_img - pp_img.mean()) ** 2).mean())

        cvs = safe_div(sigma_s, mu_s)
        cvb = safe_div(sigma_b, mu_b)

        # robust metrics
        # mad_bkg = np.median(np.abs(backgrnd_values - med_b))
        # mad_noise = np.median(np.abs(noise_values - med_n))

        # sbr_robust = safe_div(med_s, med_b)
        # snr_robust = safe_div(med_s, (1.4826 * mad_noise))
        # cnr_robust = safe_div((med_s - med_b), (1.4826 * mad_bkg))

        # rms_robust = 1.4826 * np.median(np.abs(pp_img.flatten() - np.median(pp_img)))

        effect_size = safe_div((mu_s - mu_b), np.sqrt((sigma_s**2 + sigma_b**2) / 2))

        z_factor = 1 - 3 * (sigma_s + sigma_b) / abs(mu_s - mu_b)

        h = shannon_entropy(pp_img)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sk = skew(signal_values, nan_policy='omit')
            ku = kurtosis(signal_values, nan_policy='omit', fisher=True)

        report['snr'] = snr
        report['snr_bs'] = snr_bs
        report['sbr'] = sbr
        report['cnr'] = cnr
        report['rms'] = rms
        report['psnr'] = psnr
        report['effect_size'] = effect_size
        report['z_factor'] = z_factor
        report['h'] = h
        report['sk'] = sk
        report['ku'] = ku
        report['cvs'] = cvs
        report['cvb'] = cvb
        report['mu_s'] = mu_s
        report['mu_b'] = mu_b
        report['mu_n'] = mu_n
        report['sigma_s'] = sigma_s
        report['sigma_b'] = sigma_b
        report['sigma_n'] = sigma_n
        report['med_s'] = med_s
        report['med_b'] = med_b
        report['med_n'] = med_n

        report['thresh'] = thresh

        # report['snr_robust'] = snr_robust
        # report['cnr_robust'] = cnr_robust
        # report['rms_robust'] = rms_robust
        # report['sbr_robust'] = sbr_robust
        # report['mad_bkg'] = mad_bkg
        # report['mad_noise'] = mad_noise

        if return_mask:
            return report, signal_mask
        else:
            return report

    def _pp_snr(self):
        pp_img = np.clip(self.img, 0, self.vmax)

        noise_values = pp_img[self.flat_region_mask]

        pp_bgsub_img = np.maximum(pp_img - np.mean(noise_values), 0)
        mu_n, med_n, sigma_n = _get_basic_stats(noise_values)

        noise_stats = {}
        noise_stats['noise_values'] = noise_values
        noise_stats['mu_n'] = mu_n
        noise_stats['med_n'] = med_n
        noise_stats['sigma_n'] = sigma_n

        return pp_bgsub_img, noise_stats

    def calculate_signal_stats_tiled(self, tile_size=500, clean_mask=False):
        from tqdm import tqdm

        pp_img, noise_stats = self._pp_snr()

        tiles = self._tile_array(pp_img, tile_size=tile_size)

        for tile in tqdm(tiles, desc='Processing'):
            tile_img = tiles[tile]['array']
            tiles[tile]['signal_stats'], signal_mask = self.calculate_signal_stats(
                tile_img, noise_stats, clean_mask=clean_mask, return_mask=True
            )
            tiles[tile]['signal_mask'] = signal_mask

        return tiles

    def _tile_array(self, img, tile_size=1000):
        h, w = img.shape

        cols = np.ceil(h / tile_size)
        rows = np.ceil(w / tile_size)

        tiles = {}
        for c in range(int(cols)):
            for r in range(int(rows)):
                cidx = (c * tile_size), ((c + 1) * tile_size)
                ridx = (r * tile_size), ((r + 1) * tile_size)

                _ = img[cidx[0] : cidx[1], ridx[0] : ridx[1]]

                tiles[(c, r)] = {'array': _, 'coords': (cidx, ridx)}

        return tiles

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

        if ax is not None or return_ax:
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
    vmin = filters.threshold_otsu(image)
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


def mask_signal_ostu(img):
    # 2. Threshold to separate signal from background (Otsu’s method)
    thresh = filters.threshold_otsu(img)
    signal_mask = img > thresh
    return signal_mask, thresh


def mask_signal_quant(img, q=0.5):
    # 2. Threshold to separate signal from background (Otsu’s method)
    thresh = np.quantile(img, q)
    signal_mask = img > thresh
    return signal_mask, thresh


def _clean_signal_mask(signal_mask, opening_radius=2, min_object_size=64):
    # 1. Morphological opening to remove small bumps/bridges
    selem = disk(opening_radius)
    opened = opening(signal_mask, selem)

    # 2. Remove any connected components smaller than min_object_size
    cleaned_mask = remove_small_objects(opened, min_size=min_object_size)

    return cleaned_mask


def _get_basic_stats(values):
    mean = np.mean(values)
    median = np.median(values)
    std = np.std(values)
    return mean, median, std
