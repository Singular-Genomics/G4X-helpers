from .cells import write_cells
from .images import write_he_img, write_multiplex_img
from .manage_zarr import init_viewer_zarr, link_viewer_group
from .transcripts import write_transcripts


def create_default_viewer(smp, logger):
    init_viewer_zarr(smp, logger=logger)
    write_multiplex_img(smp, logger=logger)
    write_he_img(smp, logger=logger)
    write_transcripts(smp, logger=logger)
    write_cells(smp, logger=logger)


__all__ = [
    'init_viewer_zarr',
    'write_he_img',
    'write_multiplex_img',
    'write_transcripts',
    'write_cells',
    'link_viewer_group',
]
