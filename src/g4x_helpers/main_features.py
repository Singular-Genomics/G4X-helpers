import functools
import logging
from typing import TYPE_CHECKING, Literal

from . import __version__, c, io
from . import logging_utils as logut
from . import utils as ut

if TYPE_CHECKING:
    from .g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        smp: 'G4Xoutput',
        *,
        out_dir: str | None = None,
        n_threads: int = c.DEFAULT_THREADS,
        verbose: int = 1,
        compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
        logger: logging.Logger | None = None,
        **kwargs,
    ):
        func_name = func.__name__

        if out_dir is None:
            out_dir = smp.data_dir
        else:
            out_dir = io.pathval.validate_dir_path(out_dir)
            func_out = out_dir / func_name
            if not func_out.exists():
                func_out.mkdir(parents=True, exist_ok=True)
            out_dir = func_out

        log_dir = out_dir / 'logs'
        if not log_dir.exists():
            log_dir.mkdir(parents=True, exist_ok=True)

        if logger is None:
            # TODO enable append time when testing is complete
            func_log = logut.configure_g4x_logging(
                level=verbose * 10, file_log=True, out_dir=log_dir, append_time=False, file_mode='w'
            )
        else:
            func_log = logger

        backend = io.get_backend(compute_backend)
        compute_eng = f'{backend.kind}'
        compute_eng += ' (detected)' if compute_backend == 'auto' else ''

        d = {
            'sample_dir': f'{smp.data_dir}',
            'out_dir': f'{out_dir}',
            'verbosity': f'{verbose}',
            'n_threads': f'{n_threads}',
            'compute_eng': compute_eng,
            'g4x-helpers': f'v{__version__}',
        }

        header = f'Initializing G4X-helpers [{func_name}]\n'
        msg = ut.pretty_dict_str(d)
        logut.log_msg_wrapped(header=header, msg=msg, prefix='  ', logger=func_log)

        result = func(
            smp=smp,
            out_dir=out_dir,
            n_threads=n_threads,
            compute_backend=backend.kind,
            # logger=func_log,
            **kwargs,
        )

        func_log.info(f'Completed: {func.__name__}')
        return result

    return wrapper


@_base_command
def resegment(
    smp: 'G4Xoutput',
    out_dir: str,
    segmentation_mask: str,
    mask_key: str | None = None,
    overwrite: bool = True,
    show_progress: bool = False,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):
    from .modules import aggregate, single_cell, viewer

    log = kwargs.get('logger', LOGGER)

    aggregate.aggregate_cell_data(
        smp,
        segmentation_mask=segmentation_mask,
        mask_key=mask_key,
        out_dir=out_dir,
        overwrite=overwrite,
        compute_backend=compute_backend,
        show_progress=show_progress,
        logger=log,
    )
    single_cell.process_sc_output(
        smp, out_dir=out_dir, compute_backend=compute_backend, overwrite=overwrite, logger=log
    )

    viewer.create_viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite, logger=log)
    viewer.cells.write_cells(smp, seg_name='g4x-default', overwrite=overwrite, logger=log)


@_base_command
def redemux(
    smp: 'G4Xoutput',
    out_dir: str,
    manifest: str,
    overwrite: bool = True,
    show_progress: bool = False,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):
    from .modules import aggregate, demux, single_cell, viewer

    log = kwargs.get('logger', LOGGER)

    demux.demux_raw_features(
        smp,
        manifest=manifest,
        out_dir=out_dir,
        overwrite=overwrite,
        show_progress=show_progress,
        logger=log,
    )

    aggregate.aggregate_cell_data(
        smp,
        out_dir=out_dir,
        overwrite=overwrite,
        compute_backend=compute_backend,
        show_progress=show_progress,
        logger=log,
    )
    single_cell.process_sc_output(
        smp, out_dir=out_dir, compute_backend=compute_backend, overwrite=overwrite, logger=log
    )

    viewer.create_viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite, logger=log)
    viewer.cells.write_cells(smp, seg_name='g4x-default', overwrite=overwrite, logger=log)


@_base_command
def migrate(
    g4x_obj: 'G4Xoutput',
    restore: bool = False,
    n_threads: int = c.DEFAULT_THREADS,
    *,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .schemas import migration

    if restore:
        migration.restore_backup(data_dir=g4x_obj.data_dir, sample_id=g4x_obj.sample_id, logger=logger)
    else:
        migration.migrate_g4x_data(
            data_dir=g4x_obj.data_dir, sample_id=g4x_obj.sample_id, n_threads=n_threads, logger=logger
        )


@_base_command
def validate(smp: 'G4Xoutput', **kwargs):

    log = kwargs.get('logger', LOGGER)
    report = smp.src.validation_report(raise_exception=False)

    if smp.src.is_valid_all:
        logut.log_msg_wrapped('Sample is valid:\n', report, logger=log, level='info', prefix=' ')
    else:
        logut.log_msg_wrapped('Sample validation failed:\n', report, logger=log, level='error', prefix=' ')
