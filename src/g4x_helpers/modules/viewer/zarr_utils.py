from __future__ import annotations

import colorsys
import logging
import os
import shutil

import zarr

from ... import constants as c
from ... import io
from ...schema.definition import ViewerZarr
from ..workflow import PRESET_SOURCE, reroute_source

LOGGER = logging.getLogger(__name__)
DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


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


def init_viewer_zarr(
    smp,
    *,
    out_dir: str = PRESET_SOURCE,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
) -> None:
    log = logger or LOGGER
    log.info('Running init_viewer_zarr')

    if out_dir == PRESET_SOURCE:
        out_dir = smp.smp_dir if not smp.uses_branch else smp.alt_source
    else:
        out_dir = io.pathval.validate_dir_path(out_dir)

    # out_dir = smp.smp_dir if out_dir == PRESET_SOURCE else io.pathval.validate_dir_path(out_dir)
    reroute_source(smp, out_dir, validator=ViewerZarr, overwrite=overwrite, logger=log)

    mode = 'w' if overwrite else 'a'
    root_group = zarr.open_group(smp.out.ViewerZarr.p, mode=mode, zarr_version=2)

    # write_metadata_defaults
    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

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

    (smp.out.ViewerZarr.p / 'misc').mkdir(parents=True, exist_ok=True)
    if smp.src.QCSummary.path_exists():
        shutil.copy(smp.src.QCSummary.p, smp.out.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')

    return root_group


def link_viewer_group(smp, group_name: str, overwrite: bool = True):
    target = smp.src.ViewerZarr.p / group_name
    link = smp.out.ViewerZarr.p / group_name

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


# def get_gene_metadata(viewer_dir: str):
#     g = zarr.open(viewer_dir, mode='r')
#     cmap = dict(g['transcripts'].attrs)['gene_colors']
#     cmap = {k: rgb_to_hex(v) for k, v in cmap.items()}

#     df = pl.DataFrame(dict(g['transcripts'].attrs)['gene_order'], schema=['gene_id'])
#     df = df.with_columns(pl.col('gene_id').replace(cmap).alias('color'))
#     return df
