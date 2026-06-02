import logging

from .cells import write_cells
from .images import write_he_img, write_multiplex_img
from .transcripts import write_transcripts
from .zarr_utils import init_viewer_zarr, link_viewer_group

LOGGER = logging.getLogger(__name__)


def create_default_viewer(
    smp,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    init_viewer_zarr(smp, overwrite=overwrite, logger=log)

    if smp.uses_branch:
        link_viewer_group(smp, group_name='images', overwrite=overwrite)
    else:
        write_multiplex_img(smp, overwrite=overwrite, logger=log)
        write_he_img(smp, overwrite=overwrite, logger=log)

    write_transcripts(smp, overwrite=overwrite, logger=log)
    write_cells(smp, seg_name='g4x-default', overwrite=overwrite, logger=log)


__all__ = [
    'init_viewer_zarr',
    'write_he_img',
    'write_multiplex_img',
    'write_transcripts',
    'write_cells',
    'link_viewer_group',
]
