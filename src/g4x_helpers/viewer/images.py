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


def _get_default_compressor():
    return Blosc(cname='zstd', clevel=5, shuffle=Blosc.SHUFFLE)


def _get_default_scaler(levels=4):
    return oz_scale.Scaler(downscale=2, max_layer=levels)


def _load_image_dask(smp, img_type, channel_name=None, dtype=np.uint16):
    shape = smp.shape + (3,) if img_type == 'h_and_e' else smp.shape

    data = da.from_delayed(
        dask.delayed(smp.load_image_by_type)(image_type=img_type, protein=channel_name, thumbnail=False, cached=False),
        shape=shape,
        dtype=dtype,
    )

    return data


def write_images(smp, root_group):
    # img_group = root_group.create_group('images', overwrite=True)
    # img_group.attrs['axes'] = {'unit': 'micrometer', 'pixel_per_um': c.PIXEL_PER_MICRON}

    write_muliplex_img(smp, root_group)
    write_he_img(smp, root_group)


def write_muliplex_img(smp, root_group):
    # Sort proteins to have aSMA first and Isotype last
    channel_names = sorted(smp.proteins)

    if 'aSMA' in channel_names:
        channel_names.remove('aSMA')
        channel_names = ['aSMA'] + channel_names

    if 'Isotype' in channel_names:
        channel_names.remove('Isotype')
        channel_names = channel_names + ['Isotype']

    # Prepare dask arrays for each channel
    dtype = np.uint16
    channel_arrays = []
    for ch in channel_names:
        arr = _load_image_dask(smp, img_type='protein', channel_name=ch, dtype=dtype)
        channel_arrays.append(arr)

    for ch in ['nuclear', 'eosin']:
        arr = _load_image_dask(smp, img_type=ch, channel_name=None, dtype=dtype)
        channel_arrays.append(arr)

    for arr in channel_arrays:
        if arr.ndim == 2:
            arr = arr[None, ...]  # add Z

    data = da.stack(channel_arrays, axis=0)
    axes = ['c', 'z', 'y', 'x'] if data.ndim == 4 else ['c', 'y', 'x']

    scaler = _get_default_scaler()
    compressor = _get_default_compressor()

    chunks = (1, 1, 256, 256)[-data.ndim :]
    storage_options = {'chunks': chunks, 'compressor': compressor}

    # dtype_max = np.iinfo(dtype).max
    channel_names += ['nuclear', 'eosin']
    omero = {
        'channels': [
            {
                'label': name,
                'color': sat_cols[channel_names.index(name) % len(sat_cols)],
                'active': True if name in DEFAULT_VISIBLE_CHANNELS else False,
                # 'window': {'start': 0, 'end': dtype_max, 'min': 0, 'max': dtype_max},
            }
            for name in channel_names
        ]
    }

    img_group = root_group['images']['multiplex']
    oz_writer.write_image(
        data,
        img_group,
        scaler=scaler,
        axes=axes,
        storage_options=storage_options,
        metadata={'omero': omero},
    )


def write_he_img(smp, root_group):
    channel_names = ['R', 'G', 'B']
    channel_cols = dict(zip(channel_names, ['FF0000', '00FF00', '0000FF']))

    scaler = _get_default_scaler()
    compressor = _get_default_compressor()

    dtype = np.uint8
    data = _load_image_dask(smp, img_type='h_and_e', channel_name=None, dtype=dtype)

    # TODO remove the astronaut at some point =)
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
    data = out
    ### Add rgb image to bottom left corner
    ### ### ### ### ### ### ### ### ### ###

    # Ensure channels are the first axis without reinterpreting memory.
    if data.ndim == 3 and data.shape[-1] == 3:
        data = da.moveaxis(data, -1, 0)
    elif data.ndim == 3 and data.shape[0] == 3:
        pass
    else:
        raise ValueError(f'Unexpected H&E image shape: {data.shape}')

    data = data[:, None, :, :]
    axes = ['c', 'z', 'y', 'x'] if data.ndim == 4 else ['c', 'y', 'x']

    dtype_max = np.iinfo(dtype).max
    omero = {
        'channels': [
            {
                'label': name,
                'window': {'start': 0, 'end': dtype_max, 'min': 0, 'max': dtype_max},
                'color': channel_cols[name],
                'active': True,
            }
            for name in channel_names
        ]
    }

    chunks = (1, 1, 256, 256)[-data.ndim :]
    storage_options = {'chunks': chunks, 'compressor': compressor}

    img_group = root_group['images']['h_and_e']
    oz_writer.write_image(
        data,
        img_group,
        scaler=scaler,
        axes=axes,
        storage_options=storage_options,
        metadata={'omero': omero},
    )
