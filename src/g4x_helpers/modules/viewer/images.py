from __future__ import annotations

import logging
import warnings

import dask.array as da
import numpy as np
import polars as pl
import zarr
from numcodecs import Blosc
from ome_zarr import scale as oz_scale
from ome_zarr import writer as oz_writer

from ... import io
from ...utils import get_image_shape

log = logging.getLogger(__name__)

saturated_colors = {
    'red': '#FF0000',
    'orange': '#FF8000',
    'yellow': '#FFFF00',
    'chartreuse': '#80FF00',
    'green': '#00FF00',
    'spring_green': '#00FF80',
    'cyan': '#00FFFF',
    'azure': '#007FFF',
    'blue': '#0000FF',
    'violet': '#7F00FF',
    'magenta': '#FF00FF',
    'rose': '#FF0080',
    'white': '#FFFFFF',
}

CHANNELS_SET_A = ['PanCK', 'aSMA', 'CD31', 'CD45']
CHANNELS_SET_B = ['ATPase', 'Isotype', 'cytoplasmic', 'nuclear']
CHANNELS_SET_C = ['HLA-DR', 'CD68', 'CD11c', 'CD20']
CHANNELS_SET_D = ['CD3', 'CD4', 'CD8', 'FOXP3']
CHANNELS_SET_E = ['CD34', 'PanCadherin', 'KI67', 'PD1', 'PDL1']  # < most don't have colors yet

DEFAULT_VISIBLE_CHANNELS = CHANNELS_SET_A

channel_color_map = {
    'aSMA': 'yellow',
    'CD31': 'spring_green',
    'PanCK': 'azure',
    'CD45': 'magenta',
    'HLA-DR': 'orange',
    'CD68': 'violet',
    'CD11c': 'cyan',
    'CD20': 'yellow',
    'CD3': 'spring_green',
    'CD4': 'azure',
    'CD8': 'magenta',
    'FOXP3': 'red',
    'KI67': 'chartreuse',
    'ATPase': 'violet',
    'Isotype': 'orange',
    'cytoplasmicstain': 'rose',
    'nuclearstain': 'white',
}

OMERO_DEFAULT = {
    'label': 'channel',
    'color': 'FFFFFF',
    'active': True,
    'window': {'start': 0, 'end': 1, 'min': 0, 'max': 1},
}


# region main functions
def write_images(
    zarr_path: str,
    images: dict[str, str],
    visible_channels: list[str] | None = None,
    channel_colors: dict[str, str] = {},
    channel_windows: dict[str, str] = {},
    overwrite: bool = True,
    chunk_size: int = 1024,
    use_cache: bool = False,
):

    zarr_path = io.pathval.validate_dir_path(zarr_path)

    mode = 'w' if overwrite else 'r+'
    img_group = zarr.open_group(zarr_path / 'images' / 'multiplex', mode=mode)

    # Prepare dask arrays for each channel
    for name, path in images.items():
        shape = get_image_shape(path)
        arr = io.import_image_dask(path, shape=shape, use_cache=use_cache)

        if arr.ndim == 2:
            arr = arr[None, ...]  # add Z

        images[name] = arr

    if visible_channels is None:
        visible_channels = list(images.keys())[0:4]

    channels = []
    for name, arr in images.items():
        log.debug(f'Processing channel: {name}')

        color = channel_colors.get(name, None)
        if color is None:
            color = saturated_colors[
                list(saturated_colors.keys())[list(images.keys()).index(name) % len(saturated_colors)]
            ]

        active = True if name in visible_channels else False

        window = channel_windows.get(name, None)
        if window is None:
            window = default_window_recipe(arr)
        color = color.removeprefix('#')
        ic = ImageChannel(
            arr, label=name, dtype=np.uint16, omero_attrs={'color': color, 'active': active, 'window': window}
        )
        channels.append(ic)

    write_channel_stack(img_group, channels, chunk_size=chunk_size)


def write_rgb_image(
    zarr_path: str,
    image_name: str,
    image_path: str,
    dtype: np.dtype = np.uint8,
    overwrite: bool = True,
    chunk_size: int = 1024,
    use_cache: bool = False,
):

    zarr_path = io.pathval.validate_dir_path(zarr_path)

    mode = 'w' if overwrite else 'r+'
    img_group = zarr.open_group(zarr_path / 'images' / 'h_and_e', mode=mode)

    shape = get_image_shape(image_path)
    image = io.import_image_dask(image_path, shape=shape, dtype=dtype, use_cache=use_cache)

    if image.ndim == 3 and image.shape[-1] == 3:
        log.debug('Image in YXC format, moving channel axis to front')
        image = da.moveaxis(image, -1, 0)
    elif image.ndim == 3 and image.shape[0] == 3:
        log.debug('Image already in CYX format')
    else:
        raise ValueError(f'Unexpected RGB image shape: {image.shape}')

    c1 = ImageChannel(image[0], label='R', omero_attrs={'color': 'FF0000', 'active': True})
    c2 = ImageChannel(image[1], label='G', omero_attrs={'color': '00FF00', 'active': True})
    c3 = ImageChannel(image[2], label='B', omero_attrs={'color': '0000FF', 'active': True})

    write_channel_stack(img_group, [c1, c2, c3], chunk_size=chunk_size)


def write_channel_stack(
    img_group, channels: list[ImageChannel], chunk_size: int = 256, levels: int = 4, clevel: int = 5
):
    channel_arrays = [ch.img for ch in channels]
    data = da.stack(channel_arrays, axis=0)
    axes = ['c', 'z', 'y', 'x'] if data.ndim == 4 else ['c', 'y', 'x']

    scaler = oz_scale.Scaler(downscale=2, max_layer=levels)
    compressor = Blosc(cname='zstd', clevel=clevel, shuffle=Blosc.SHUFFLE)

    chunks = (1, 1, chunk_size, chunk_size)[-data.ndim :]
    storage_options = {'chunks': chunks, 'compressor': compressor}

    omero = {'channels': [ch.omero for ch in channels]}

    _write_image_withouth_storage_warning(
        data,
        img_group,
        scaler=scaler,
        axes=axes,
        storage_options=storage_options,
        metadata={'omero': omero},
    )


# region private functions
class ImageChannel:
    def __init__(self, img, label: str = 'channel', dtype: np.dtype = None, omero_attrs: dict = {}):
        self.img = da.array(img, dtype=dtype)
        self.label = label
        self.attrs = omero_attrs
        self.attrs['label'] = self.label

        if 'window' not in self.attrs:
            dtype_max = np.iinfo(self.img.dtype).max
            self.attrs['window'] = {'start': 0, 'end': dtype_max, 'min': 0, 'max': dtype_max}

        self.omero = OMERO_DEFAULT.copy()
        self.omero.update(self.attrs)


def default_window_recipe(arr):
    arr_max = int(arr.max().compute())
    clip_vmax = int(da.percentile(arr.ravel(), 99.5).compute())
    clip_vmin = int(clip_vmax * 0.10)
    window = {'min': 0, 'max': arr_max, 'start': clip_vmin, 'end': clip_vmax}
    return window


def _write_image_withouth_storage_warning(*args, **kwargs):
    # ome-zarr currently passes storage params to da.to_zarr via **kwargs;
    # suppress only that known deprecation warning until upstream updates.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            category=FutureWarning,
            message=r'Passing storage-related arguments via \*\*kwargs is deprecated\..*',
        )
        return oz_writer.write_image(*args, **kwargs)


def _determine_visible_channels(available_channels: list[str]) -> list[str]:
    vi_chs = DEFAULT_VISIBLE_CHANNELS
    channels = available_channels

    num_def = len(vi_chs)
    selected = [ch for ch in vi_chs if ch in channels]
    remaining = [ch for ch in channels if ch not in vi_chs]

    for stain in ['nuclearstain', 'cytoplasmicstain']:
        remaining.remove(stain)

    n_missing = num_def - len(selected)

    for i in range(n_missing):
        if remaining:
            selected.append(remaining.pop(0))

    return selected


# region testing
def _add_rgb_astronaut_to_img(data):
    ### Add rgb image to bottom left corner
    from skimage import data as skdata

    astro = skdata.astronaut()
    h, w, _ = astro.shape
    H, W, _ = data.shape

    row0 = H - h
    row1 = H
    col0 = 0
    col1 = w

    # Ensure dtype matches (optional but recommended)
    np_img2 = astro.astype(data.dtype, copy=False)

    # Overwrite that region
    out = data.copy()
    out[row0:row1, col0:col1, :] = da.from_array(np_img2, chunks=(h, w, 3))

    return out


# region metadata handling
def get_channel_metadata(img_group):
    channels = dict(img_group.attrs)['omero']['channels']
    rows = []
    for channel in channels:
        # window = channel.pop('window')
        window = channel['window']
        channel = {k: v for k, v in channel.items() if k != 'window'}

        channel.update(window)
        rows.append(channel)

    return pl.DataFrame(rows).select(['label', 'active', 'color', 'min', 'max', 'start', 'end'])


def apply_channel_metadata(img_group, new_data):
    old_data = get_channel_metadata(img_group)

    if not set(old_data.columns) == set(new_data.columns):
        raise ValueError(
            f'New data columns {new_data.columns} do not match existing channel metadata columns {old_data.columns}'
        )

    existing_labels = old_data['label'].to_list()

    if set(existing_labels) != set(new_data['label']):
        raise ValueError(
            f'New data labels {new_data["label"]} do not match existing channel metadata labels {existing_labels}'
        )

    order_map = {label: i for i, label in enumerate(existing_labels)}

    new_data = (
        new_data.with_columns(pl.col('label').replace_strict(order_map).alias('_order')).sort('_order').drop('_order')
    )

    channels = []
    for row in new_data.iter_rows(named=True):
        window = {k: row[k] for k in ['min', 'max', 'start', 'end']}
        channel_info = {k: row[k] for k in ['label', 'active', 'color']}
        channel_info['window'] = window
        channels.append(channel_info)

    attrs = dict(img_group.attrs)

    omero = attrs['omero']
    omero['channels'] = channels

    multiscales = attrs['multiscales']
    multiscales[0]['metadata']['omero']['channels'] = channels

    img_group.attrs['omero'] = omero
    img_group.attrs['multiscales'] = multiscales
