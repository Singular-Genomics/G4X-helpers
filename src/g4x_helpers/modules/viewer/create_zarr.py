from __future__ import annotations

import logging
import shutil

import zarr

from ... import c, io
from ...schema.definition import ViewerZarr
from ..workflow import PRESET_SOURCE, reroute_source

LOGGER = logging.getLogger(__name__)
DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


# @g4x_workflow
def init_viewer_zarr(
    smp,
    *,
    out_dir: str = PRESET_SOURCE,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
) -> None:
    log = logger or LOGGER
    log.info('Running init_viewer_zarr')

    out_dir = smp.data_dir if out_dir == PRESET_SOURCE else io.pathval.validate_dir_path(out_dir)
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


def link_viewer_group(smp, out_dir, group_name: str, overwrite: bool = True):
    link = out_dir / smp.src.ViewerZarr.DEFAULT_TARGET_PATH / group_name
    
    if link.exists() or link.is_symlink():
        if not overwrite:
            raise FileExistsError(f'Link {link} already exists and overwrite is set to False.')

        if link.is_symlink() or link.is_file():
            link.unlink()

        else:
            shutil.rmtree(link)

    link.symlink_to(smp.src.ViewerZarr.p, target_is_directory=True)
