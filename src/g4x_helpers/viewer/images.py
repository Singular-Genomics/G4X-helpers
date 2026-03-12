import warnings
from typing import Literal

import dask
import dask.array as da
import numpy as np
from numcodecs import Blosc
from ome_zarr import scale as oz_scale
from ome_zarr import writer as oz_writer

DEFAULT_VISIBLE_CHANNELS = ['PanCK', 'aSMA', 'CD31', 'CD45']

# TODO apply these properly to channels
sat_cols = [
    'FF0000',
    'FF8000',
    'FFFF00',
    '80FF00',
    '00FF00',
    '00FF80',
    '00FFFF',
    '007FFF',
    '0000FF',
    '7F00FF',
    'FF00FF',
    'FF0080',
]

sat_cols2 = [
    '00FFFF',
    'FF8000',
    '0000FF',
    '80FF00',
    'FF00FF',
    '00FF80',
    'FF0000',
    '007FFF',
    'FFFF00',
    '7F00FF',
    '00FF00',
    'FF0080',
]

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


def write_muliplex_img(smp, root_group):
    # TODO this should happen in the sample class
    # Sort proteins to have aSMA first and Isotype last
    channel_names = sorted(smp.proteins, key=str.lower)

    if 'Isotype' in channel_names:
        channel_names.remove('Isotype')
        channel_names = channel_names + ['Isotype']

    # Prepare dask arrays for each channel
    dtype = np.uint16
    channel_arrays = []
    for ch in channel_names:
        arr = _load_image_dask(smp, img_type='protein', protein_name=ch, dtype=dtype)
        channel_arrays.append(arr)

    for ch in ['nuclear', 'cytoplasmic']:
        arr = _load_image_dask(smp, img_type=ch, protein_name=None, dtype=dtype)
        channel_arrays.append(arr)

    for arr in channel_arrays:
        if arr.ndim == 2:
            arr = arr[None, ...]  # add Z

    # build the channels and their metadata
    channels = []
    for arr, name in zip(channel_arrays, channel_names):
        color = sat_cols2[channel_names.index(name) % len(sat_cols2)]
        active = True if name in DEFAULT_VISIBLE_CHANNELS else False
        arr_max = int(arr.max().compute())
        clip_vmax = float(da.percentile(arr, 99).compute())
        clip_vmin = clip_vmax * 0.05
        window = {'min': 0, 'max': arr_max, 'start': clip_vmin, 'end': clip_vmax}
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
