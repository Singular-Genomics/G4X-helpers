from __future__ import annotations

import colorsys
import logging
import os
import shutil

import zarr

from ... import constants as c
from ... import io

LOGGER = logging.getLogger(__name__)


def setup_viewer_zarr(
    zarr_path: str,
    overwrite: bool = True,
) -> None:

    zarr_path = io.pathval.validate_dir_parent(zarr_path)

    mode = 'w' if overwrite else 'a'
    root_group = zarr.open_group(zarr_path, mode=mode, zarr_version=2)

    img_group = root_group.create_group('images', overwrite=overwrite)
    img_group.attrs['axes'] = {'unit': 'micrometer', 'pixel_per_um': c.PIXEL_PER_MICRON}

    img_group.create_group('multiplex', overwrite=overwrite)
    img_group.create_group('h_and_e', overwrite=overwrite)

    tx_group = root_group.create_group('transcripts', overwrite=overwrite)
    tx_group.attrs['gene_order'] = []
    tx_group.attrs['gene_colors'] = {}
    tx_group.attrs['layer_config'] = {
        'layers': 1,
        'tile_size': 1,
        'layer_height': 1,
        'layer_width': 1,
        'coordinate_order': ['default_x', 'default_y'],
    }

    cell_group = root_group.create_group('cells', overwrite=overwrite)
    cell_group.attrs['segmentation_sources'] = {}
    cell_group.attrs['segmentation_order'] = []

    (zarr_path / 'misc').mkdir(parents=True, exist_ok=True)

    return root_group


def link_viewer_group(smp, group_name: str, overwrite: bool = True):
    target = smp.smp_dir / c.FILE_VIEWER_ZARR / group_name
    link = smp.src.ViewerZarr.p / group_name

    # Compute target relative to the link's parent directory
    relative_target = os.path.relpath(target, start=link.parent)

    if link.exists() and not overwrite:
        raise FileExistsError(f"Zarr group '{link}' already exists. Set overwrite=True to replace it.")

    if link.is_symlink() or link.is_file():
        link.unlink()

    elif link.exists():
        shutil.rmtree(link)

    # Create the symlink
    link.symlink_to(relative_target, target_is_directory=True)


# this function satisfies both zarr 2 and zarr 3 APIs, trying different combinations of parameters until one works
def create_array(group, name, data, compressor=None, chunks=None):
    create = getattr(group, 'create_array', None) or group.create_dataset

    attempts = [
        {'chunks': chunks, 'compressor': compressor},
        {'chunk_shape': chunks, 'compressors': [compressor] if compressor is not None else None},
        {'chunks': chunks, 'compressors': [compressor] if compressor is not None else None},
        {'chunk_shape': chunks, 'compressor': compressor},
        {},
    ]

    last_error = None
    for kwargs in attempts:
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        try:
            return create(name, data=data, **kwargs)
        except TypeError as e:
            last_error = e

    raise last_error


def calculate_chunks(arr, target_mb=4):
    TARGET_BYTES = target_mb * 1024 * 1024  # 4 MiB
    bytes_per_row = arr.dtype.itemsize if arr.ndim == 1 else arr.dtype.itemsize * arr.shape[1]
    row_chunk = max(1, TARGET_BYTES // bytes_per_row)
    return (row_chunk, *arr.shape[1:])


def hex_to_rgb(hex_color, normalized=False):
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    if normalized:
        return tuple(v / 255 for v in rgb)
    return rgb


def rgb_to_hex(rgb, normalized=False):
    if normalized:
        rgb = tuple(int(v * 255) for v in rgb)

    return '#{:02x}{:02x}{:02x}'.format(*rgb)


def hsv_to_hex(h, s, v):
    # wrap hue into [0,1]
    h = h % 1.0
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


# region metdata handling
def get_viewer_group(viewer_dir, group: str, subgroup: str | None = None):
    img_group = zarr.open(viewer_dir / group, mode='r+')

    if subgroup is None:
        return img_group

    if subgroup not in img_group:
        raise ValueError(f'Subgroup "{subgroup}" not found in the data. Available subgroups: {list(img_group.keys())}')

    return img_group[subgroup]
