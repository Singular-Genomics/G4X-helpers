import functools
import logging
from typing import Literal

from . import __version__, io
from . import logging_utils as logut
from . import utils as ut
from .g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        smp_dir: str,
        *,
        out_dir: str | None = None,
        verbose: int = 1,
        downstream: bool = True,
        compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
        logger: logging.Logger | None = None,
        **kwargs,
    ):
        smp_dir = io.pathval.validate_dir_path(smp_dir)

        if out_dir is None:
            out_dir = smp_dir
        else:
            out_dir = io.pathval.validate_dir_path(out_dir)

        if logger is None:
            log_dir = out_dir / 'logs'
            if not log_dir.exists():
                log_dir.mkdir(parents=True, exist_ok=True)

            logger = logut.configure_g4x_logging(
                level='INFO', file_log=True, out_dir=log_dir, append_time=True, file_mode='w'
            )

        backend = io.get_backend(compute_backend)
        compute_eng = f'{backend.kind}'
        compute_eng += ' (auto-detected)' if compute_backend == 'auto' else ''

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
        logut.log_msg_wrapped(header=header, msg=msg, prefix='  ', logger=logger)

        try:
            result = func(
                smp_dir=smp_dir,
                out_dir=out_dir,
                downstream=downstream,
                compute_backend=backend.kind,
                logger=logger,
                **kwargs,
            )
            logger.info(f'Completed: [{func.__name__}]\n')
            return result

        except Exception as e:
            logger.error(f'{str(e)}')
            raise e

    return wrapper


@_base_command
def demux(
    smp_dir: str,
    manifest: str,
    *,
    out_dir: str | None = None,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):
    from .modules import aggregate, demux, single_cell, viewer

    log = kwargs.get('logger', LOGGER)
    smp = G4Xoutput(smp_dir, alt_source=out_dir)

    demux.demux_raw_features(
        smp,
        manifest=manifest,
        out_dir=out_dir,
        overwrite=overwrite,
        show_progress=show_progress,
        logger=log,
    )

    if downstream:
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

        if smp.alt_source is not None:
            viewer.init_viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite, logger=log)
            viewer.link_viewer_group(smp, out_dir, group_name='images', overwrite=overwrite)

        viewer.write_transcripts(smp, overwrite=overwrite, logger=log)
        viewer.write_cells(smp, seg_name='g4x-default', overwrite=overwrite, logger=log)

    return smp


@_base_command
def aggregate(
    smp_dir: str,
    segmentation_mask: str,
    *,
    mask_key: str | None = None,
    out_dir: str,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
):
    from .modules import aggregate, single_cell, viewer

    log = kwargs.get('logger', LOGGER)

    smp = G4Xoutput(smp_dir, alt_source=out_dir)
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

    if downstream:
        single_cell.process_sc_output(
            smp, out_dir=out_dir, compute_backend=compute_backend, overwrite=overwrite, logger=log
        )

        if smp.alt_source is not None:
            viewer.init_viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite, logger=log)
            viewer.link_viewer_group(smp, out_dir, group_name='images', overwrite=overwrite)

        viewer.write_transcripts(smp, overwrite=overwrite, logger=log)
        viewer.write_cells(smp, seg_name='g4x-default', overwrite=overwrite, logger=log)

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
    from .modules import migrate
    
    log = kwargs.get('logger', LOGGER)

    migrate.migrate_sample(
        sample_dir=smp_dir, out_dir=out_dir, roi_coords=roi_coords, downstream=downstream, logger=log
    )


@_base_command
def validate(smp_dir: str, **kwargs):
    log = kwargs.get('logger', LOGGER)
    smp = G4Xoutput(smp_dir)
    report = smp.src.validation_report(raise_exception=False)

    if smp.src.is_valid_all:
        logut.log_msg_wrapped('Sample is valid:\n', report, logger=log, level='info', prefix=' ')
    else:
        logut.log_msg_wrapped('Sample validation failed:\n', report, logger=log, level='error', prefix=' ')
