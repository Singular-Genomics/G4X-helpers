from . import cells as cell
from . import images as img
from . import transcripts as tx
from .cells import write_cells
from .create_zarr import init_viewer_zarr, link_viewer_group
from .images import write_he_img, write_muliplex_img
from .transcripts import write_transcripts


def create_default_viewer(smp):
    init_viewer_zarr(smp)
    write_muliplex_img(smp)
    write_he_img(smp)
    write_transcripts(smp)
    write_cells(smp)


__all__ = [
    'cell',
    'img',
    'tx',
    'init_viewer_zarr',
    'write_he_img',
    'write_muliplex_img',
    'write_transcripts',
    'write_cells',
    'link_viewer_group',
]
