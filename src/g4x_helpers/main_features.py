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
        out_dir = io.pathval.validate_dir_path(out_dir)

        if out_dir != smp.data_dir:
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
def testo(
    smp: 'G4Xoutput',
    out_dir: str,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    overwrite: bool = True,
    **kwargs,
):
    log = kwargs.get('logger', LOGGER)
    log.info('This is a test command.')
    log.info(f'Output directory: {out_dir}')
    log.info(f'Compute backend: {compute_backend}')
    log.info(f'Overwrite: {overwrite}')


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
def validate(
    g4x_obj: 'G4Xoutput',
    *,
    logger: logging.Logger,
    **kwargs,
):
    from .schemas import validate

    try:
        validate.validate_g4x_data(
            g4x_obj.data_dir, schema_name='base_schema', formats={'sample_id': g4x_obj.sample_id}, report='long'
        )
    except validate.ValidationError as e:
        print('Directory structure validation failed!: ', e)

    try:
        print('\n')
        validate.validate_file_schemas(g4x_obj.data_dir, verbose=True)
    except validate.ValidationError as e:
        print('File schema validation failed!: ', e)


# def _base_command(func):
#     """Decorator to apply standard command initialization logic."""

#     @functools.wraps(func)
#     def wrapper(
#         g4x_obj: 'G4Xoutput',
#         out_dir: str | None = None,
#         n_threads: int = c.DEFAULT_THREADS,
#         verbose: int = 1,
#         file_logger: bool = True,
#         **kwargs,
#     ):
#         func_name = func.__name__

#         if func_name not in ('migrate', 'validate'):
#             g4x_obj.validate()

#         out_dir = g4x_obj.data_dir if out_dir is None else out_dir
#         out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

#         if out_dir != g4x_obj.data_dir:
#             func_out = out_dir / func_name
#             func_out.mkdir(parents=True, exist_ok=True)
#             out_dir = func_out

#         gap = 12
#         click.secho(f'\nStarting: {func_name}\n', bold=True)
#         print_k_v('sample_dir', f'{g4x_obj.data_dir}', gap)
#         print_k_v('out_dir', f'{out_dir}', gap)
#         print_k_v('n_threads', f'{n_threads}', gap)
#         print_k_v('verbosity', f'{verbose}', gap)
#         print_k_v('g4x-helpers', f'v{__version__}', gap)
#         click.echo('')

#         log_dir = g4x_obj.data_dir / 'g4x_helpers' / 'logs'
#         log_dir.mkdir(parents=True, exist_ok=True)
#         if not kwargs.get('logger', None):
#             logger = utils.setup_logger(
#                 logger_name=func_name, out_dir=log_dir, stream_level=verbose, file_logger=file_logger
#             )
#         logger.info(f'Running {func_name} with G4X-helpers v{__version__}')

#         result = func(
#             g4x_obj=g4x_obj,
#             out_dir=out_dir,
#             n_threads=n_threads,
#             verbose=verbose,
#             logger=logger,
#             **kwargs,
#         )

#         click.secho(f'\nCompleted: {func.__name__}', bold=True, fg='green')
#         return result

#     return wrapper
