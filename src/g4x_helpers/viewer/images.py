import warnings
from typing import Literal

import dask
import dask.array as da
import numpy as np
from numcodecs import Blosc
from ome_zarr import scale as oz_scale
from ome_zarr import writer as oz_writer

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
    'cytoplasmic': 'rose',
    'nuclear': 'white',
}

OMERO_DEFAULT = {
    'label': 'channel',
    'color': 'FFFFFF',
    'active': True,
    'window': {'start': 0, 'end': 1, 'min': 0, 'max': 1},
}


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


def write_images(smp, root_group):
    write_muliplex_img(smp, root_group)
    write_he_img(smp, root_group)


def default_window_recipe(arr):
    arr_max = int(arr.max().compute())
    clip_vmax = int(da.percentile(arr.ravel(), 99.5).compute())
    clip_vmin = int(clip_vmax * 0.10)
    window = {'min': 0, 'max': arr_max, 'start': clip_vmin, 'end': clip_vmax}
    return window


def write_muliplex_img(smp, root_group):
    # TODO this should happen in the sample class
    # Sort proteins to have aSMA first and Isotype last
    channel_names = sorted(smp.proteins, key=str.lower)

    dtype = np.uint16
    channel_arrays = []

    # Prepare dask arrays for each channel
    if smp.proteins:
        if 'Isotype' in channel_names:
            channel_names.remove('Isotype')
            channel_names = channel_names + ['Isotype']

        for ch in channel_names:
            arr = _load_image_dask(smp, img_type='protein', protein_name=ch, dtype=dtype)
            channel_arrays.append(arr)

    stain_channels = ['cytoplasmic', 'nuclear']
    channel_names.extend(stain_channels)

    for ch in stain_channels:
        arr = _load_image_dask(smp, img_type=ch, protein_name=None, dtype=dtype)
        channel_arrays.append(arr)

    for i, arr in enumerate(channel_arrays):
        if arr.ndim == 2:
            channel_arrays[i] = arr[None, ...]  # add Z

    # build the channels and their metadata
    channels = []
    for arr, name in zip(channel_arrays, channel_names):
        if name in channel_color_map:
            color = saturated_colors[channel_color_map[name]]
        else:
            color = saturated_colors[list(saturated_colors.keys())[channel_names.index(name) % len(saturated_colors)]]
        # TODO define what to do when t
        active = True if name in DEFAULT_VISIBLE_CHANNELS else False
        window = default_window_recipe(arr)
        color = color.removeprefix('#')
        ic = ImageChannel(
            arr, label=name, dtype=np.uint16, omero_attrs={'color': color, 'active': active, 'window': window}
        )
        channels.append(ic)

    img_group = root_group['images']['multiplex']
    write_channel_stack(img_group, channels)


def write_channel_stack(img_group, channels):
    channel_arrays = [ch.img for ch in channels]
    data = da.stack(channel_arrays, axis=0)
    axes = ['c', 'z', 'y', 'x'] if data.ndim == 4 else ['c', 'y', 'x']

    scaler = _get_default_scaler()
    compressor = _get_default_compressor()

    chunks = (1, 1, 256, 256)[-data.ndim :]
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


def write_he_img(smp, root_group):
    dtype = np.uint8
    image = _load_image_dask(smp, img_type='h_and_e', protein_name=None, dtype=dtype)
    image = _add_rgb_astronaut_to_img(image)

    if image.ndim == 3 and image.shape[-1] == 3:
        image = da.moveaxis(image, -1, 0)
    elif image.ndim == 3 and image.shape[0] == 3:
        pass
    else:
        raise ValueError(f'Unexpected H&E image shape: {image.shape}')

    c1 = ImageChannel(image[0], label='R', omero_attrs={'color': 'FF0000', 'active': True})
    c2 = ImageChannel(image[1], label='G', omero_attrs={'color': '00FF00', 'active': True})
    c3 = ImageChannel(image[2], label='B', omero_attrs={'color': '0000FF', 'active': True})

    img_group = root_group['images']['h_and_e']
    write_channel_stack(img_group, [c1, c2, c3])


# region helpers
def _get_default_compressor():
    return Blosc(cname='zstd', clevel=5, shuffle=Blosc.SHUFFLE)


def _get_default_scaler(levels=4):
    return oz_scale.Scaler(downscale=2, max_layer=levels)


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


def _load_image_dask(
    smp,
    img_type: str = Literal['protein', 'h_and_e', 'nuclear', 'cytoplasmic'],
    protein_name: str | None = None,
    dtype=np.uint16,
    use_cache: bool = False,
):
    shape = smp.shape + (3,) if img_type == 'h_and_e' else smp.shape

    if img_type == 'protein':
        data = da.from_delayed(
            dask.delayed(smp.load_protein_image)(protein=protein_name, use_cache=use_cache),
            shape=shape,
            dtype=dtype,
        )
    elif img_type == 'h_and_e':
        data = da.from_delayed(
            dask.delayed(smp.load_he_image)(use_cache=use_cache),
            shape=shape,
            dtype=dtype,
        )
    elif img_type == 'nuclear':
        data = da.from_delayed(
            dask.delayed(smp.load_nuclear_image)(use_cache=use_cache),
            shape=shape,
            dtype=dtype,
        )
    elif img_type == 'cytoplasmic':
        data = da.from_delayed(
            dask.delayed(smp.load_cytoplasmic_image)(use_cache=use_cache),
            shape=shape,
            dtype=dtype,
        )

    return data


# region testing
def _add_rgb_astronaut_to_img(data):
    ### ### ### ### ### ### ### ### ### ###
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

    ### Add rgb image to bottom left corner
    ### ### ### ### ### ### ### ### ### ###
    return out
