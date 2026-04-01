import logging
import shutil

from ... import io
from ...schema.definition import ViewerZarr
from ..workflow import DEFAULT_INPUT, prepare_output
from .cells import prepare_cell_group_input, write_cells
from .images import write_images
from .setup import setup_zarr_tree
from .transcripts import write_transcripts

LOGGER = logging.getLogger(__name__)


def create_viewer_zarr(
    smp,
    out_dir: str = DEFAULT_INPUT,
    overwrite=True,
    protein_list: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> None:
    log = logger or LOGGER
    log.info('Running create_viewer_zarr')

    out_dir = smp.data_dir if out_dir == DEFAULT_INPUT else io.pathval.validate_dir_path(out_dir)
    prepare_output(smp, out_dir, validator=ViewerZarr, overwrite=overwrite, logger=log)

    if protein_list is not None:
        smp.set_proteins(protein_list)

    root_group = setup_zarr_tree(smp.data_dir, overwrite=overwrite)

    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

    write_images(smp, root_group, logger=log)
    write_transcripts(smp, root_group, overwrite=True, logger=log)

    prechew = prepare_cell_group_input(smp)
    write_cells(smp, root_group, prechew=prechew, overwrite=True, logger=log)

    if smp.src.QCSummary.path_exists:
        shutil.copy(smp.src.QCSummary.p, smp.out.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')
