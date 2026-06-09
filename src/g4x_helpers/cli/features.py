import functools
import logging
from typing import Literal

from .. import __version__, io
from .. import constants as c
from .. import logging_utils as logut
from .. import sample_ops as ops
from .. import utils as ut
from ..g4x_output import G4Xoutput

log = logging.getLogger(__name__)


def _create_branch(sample_dir: str, name: str):
    HELPERS_DIR_NAME = 'g4x-helpers'

    sample_dir = io.pathval.validate_dir_path(sample_dir)
    branch_dir = sample_dir / HELPERS_DIR_NAME / name

    if not branch_dir.exists():
        branch_dir.mkdir(parents=True, exist_ok=True)
    return branch_dir


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        smp_dir: str,
        *,
        out_dir: str | None = None,
        verbose: int = 1,
        downstream: bool = True,
        backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
        logger: logging.Logger | None = None,
        **kwargs,
    ):
        smp_dir = io.pathval.validate_dir_path(smp_dir)

        if out_dir is None:
            out_dir = _create_branch(smp_dir, func.__name__)
        else:
            out_dir = io.pathval.validate_dir_path(out_dir)

        if logger is None:
            log_dir = out_dir / 'logs'
            if not log_dir.exists():
                log_dir.mkdir(parents=True, exist_ok=True)

            lvl = logut.verbose_to_level(verbose)
            logger = logut.configure_g4x_logging(
                level=lvl, file_log=True, out_dir=log_dir, append_time=True, file_mode='w'
            )

        compute_backend = io.get_backend(backend)
        compute_eng = f'{compute_backend.kind}'
        compute_eng += ' (auto-detected)' if backend == 'auto' else ''

        d = {
            'sample_dir': f'{smp_dir}',
            'out_dir': f'{out_dir}',
            'downstream': f'{downstream}',
            'verbosity': f'{verbose}',
            'compute_eng': compute_eng,
            'g4x-helpers': f'v{__version__}',
        }

        header = f'Initializing G4X-helpers [{func.__name__}]\n'
        msg = ut.pretty_dict_str(d)
        logut.log_msg_wrapped(header=header, msg=msg, prefix='  ')

        try:
            result = func(
                smp_dir=smp_dir,
                out_dir=out_dir,
                downstream=downstream,
                backend=compute_backend.kind,
                **kwargs,
            )
            logger.info(f'Completed: [{func.__name__}]\n')
            return result

        except Exception as e:
            logger.error(f'{str(e)}')
            raise e

    return wrapper


@_base_command
def redemux(
    smp_dir: str,
    manifest: str,
    *,
    out_dir: str | None = None,
    batch_size: int = c.DEFAULT_BATCH_SIZE,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):

    smp = G4Xoutput(smp_dir, alt_source=out_dir)

    ops.demux(
        smp, manifest=manifest, out_dir=out_dir, overwrite=overwrite, batch_size=batch_size, show_progress=show_progress
    )

    if downstream:
        ops.aggregate(smp, out_dir=out_dir, overwrite=overwrite, backend=backend, show_progress=show_progress)
        ops.sc_process(smp, out_dir=out_dir, backend=backend, overwrite=overwrite)
        ops.viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite)

    return smp


@_base_command
def resegment(
    smp_dir: str,
    segmentation_mask: str,
    *,
    out_dir: str,
    mask_key: str | None = None,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):

    smp = G4Xoutput(smp_dir, alt_source=out_dir)
    ops.aggregate(
        smp,
        cell_mask=segmentation_mask,
        mask_key=mask_key,
        out_dir=out_dir,
        overwrite=overwrite,
        backend=backend,
        show_progress=show_progress,
    )

    if downstream:
        ops.sc_process(smp, out_dir=out_dir, backend=backend, overwrite=overwrite)
        ops.viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite)

    return smp


@_base_command
def migrate(
    smp_dir: str,
    out_dir: str,
    *,
    roi_coords: tuple | None = None,
    downstream: bool = True,
    **kwargs,
) -> None:
    from ..modules import migrate

    migrate.migrate_sample(
        sample_dir=smp_dir,
        out_dir=out_dir,
        roi_coords=roi_coords,
        downstream=downstream,
        **kwargs,
    )


def migrate_check(smp_dir: str) -> None:
    from ..modules import migrate

    migrate.status(sample_dir=smp_dir)


def validate(smp_dir: str, **kwargs):

    from ..schema import FileTree

    ft = FileTree(smp_dir)
    report = ft.validation_report(raise_exception=False)
    print(report)


def cell_metadata(viewer_zarr, import_metadata, export_metadata, segmentation):
    from ..modules.viewer import cells

    if import_metadata is not None and export_metadata is not None:
        raise ValueError('--import-metadata and --export-metadata cannot be used together.')

    seg_group = cells.get_seg_group(viewer_zarr, segmentation)

    if export_metadata is not None:
        meta = cells.get_cell_metadata(seg_group)
        meta.write_csv(export_metadata)
    if import_metadata is not None:
        cells.apply_viewer_metadata(seg_group, import_metadata)
