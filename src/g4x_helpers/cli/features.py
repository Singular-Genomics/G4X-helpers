import functools
import logging
from typing import TYPE_CHECKING, Literal

from .. import __version__, io
from .. import constants as c
from .. import sample_ops as ops
from .. import utils as ut
from ..g4x_sample import G4Xsample

if TYPE_CHECKING:
    from pathlib import Path

log = logging.getLogger(__name__)


def _create_branch(sample_dir: str, name: str) -> 'Path':
    HELPERS_DIR_NAME = 'g4x-helpers'

    sample_dir = io.pathval.validate_dir_path(sample_dir)
    branch_dir = sample_dir / HELPERS_DIR_NAME / name

    if not branch_dir.exists():
        branch_dir.mkdir(parents=True, exist_ok=True)
    return branch_dir


def _base_command(func) -> None:
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

            lvl = ut.verbose_to_level(verbose)
            logger = ut.configure_logging(level=lvl, file_log=True, out_dir=log_dir, append_time=True, file_mode='w')

        compute_backend = io.get_backend(backend)
        compute_eng = f'{compute_backend.kind}'
        compute_eng += ' (auto-detected)' if backend == 'auto' else ''

        d = {
            'sample_dir': f'{smp_dir.resolve()}',
            'out_dir': f'{out_dir.resolve()}',
            'downstream': f'{downstream}',
            'verbosity': f'{verbose}',
            'compute_eng': compute_eng,
            'g4x-helpers': f'v{__version__}',
        }

        header = f'Initializing G4X-helpers [{func.__name__}]\n'
        msg = ut.pretty_dict_str(d)
        ut.log_msg_wrapped(header=header, msg=msg, prefix='  ')

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
    out_dir: str | None = None,
    *,
    batch_size: int = c.DEFAULT_BATCH_SIZE,
    max_ham_dist: int = 2,
    min_delta: int = 2,
    demux_length: int = 15,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
) -> None:

    smp = G4Xsample(smp_dir, alt_source=out_dir)

    ops.demux(
        smp,
        manifest=manifest,
        out_dir=out_dir,
        overwrite=overwrite,
        max_ham_dist=max_ham_dist,
        min_delta=min_delta,
        demux_length=demux_length,
        batch_size=batch_size,
        show_progress=show_progress,
        **kwargs,
    )

    if downstream:
        ops.aggregate(smp, out_dir=out_dir, overwrite=overwrite, backend=backend, show_progress=show_progress)
        ops.sc_process(smp, out_dir=out_dir, backend=backend, overwrite=overwrite)
        ops.viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite)


@_base_command
def resegment(
    smp_dir: str,
    cell_mask: str,
    out_dir: str,
    *,
    mask_key: str | None = None,
    overwrite: bool = True,
    downstream: bool = True,
    show_progress: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
) -> None:

    smp = G4Xsample(smp_dir, alt_source=out_dir)
    ops.aggregate(
        smp,
        cell_mask=cell_mask,
        mask_key=mask_key,
        out_dir=out_dir,
        overwrite=overwrite,
        backend=backend,
        show_progress=show_progress,
    )

    if downstream:
        ops.sc_process(smp, out_dir=out_dir, backend=backend, overwrite=overwrite)
        ops.viewer_zarr(smp, out_dir=out_dir, overwrite=overwrite)


@_base_command
def migrate(
    smp_dir: str,
    out_dir: str,
    *,
    roi_coords: tuple | None = None,
    downstream: bool = True,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
) -> None:

    ops.migrate(
        legacy_dir=smp_dir,
        out_dir=out_dir,
        roi_coords=roi_coords,
        downstream=downstream,
        backend=backend,
        **kwargs,
    )


def migrate_check(smp_dir: str) -> None:
    from ..modules import migrate

    print(migrate.status(legacy_dir=smp_dir))


def validate(smp_dir: str, raw_only: bool = False) -> None:

    from ..schema import FileTree

    ft = FileTree(smp_dir)
    report = ft.validation_report(raw_only=raw_only, raise_exception=False)
    print(report)


# region viewer metadata
def _viewer_metadata_command(func) -> None:

    @functools.wraps(func)
    def wrapper(
        viewer_zarr,
        *,
        import_metadata: str | None = None,
        export_metadata: str | None = None,
        **kwargs,
    ):

        if not import_metadata and not export_metadata:
            raise ValueError('Please provide one of --import-metadata or --export-metadata options')

        if import_metadata is not None and export_metadata is not None:
            raise ValueError('--import-metadata and --export-metadata cannot be used together')

        if import_metadata is not None:
            import_metadata = io.pathval.validate_file_path(import_metadata)
        if export_metadata is not None:
            export_metadata = io.pathval.validate_file_parent(export_metadata)

        result = func(
            viewer_zarr=viewer_zarr,
            import_metadata=import_metadata,
            export_metadata=export_metadata,
            **kwargs,
        )
        return result

    return wrapper


@_viewer_metadata_command
def cell_metadata(
    viewer_zarr,
    *,
    import_metadata: str | None = None,
    export_metadata: str | None = None,
    segmentation: str = 'g4x_default_segmentation',
) -> None:
    from ..modules.viewer import cells
    from ..modules.viewer import zarr_utils as zu

    seg_group = zu.get_viewer_group(viewer_zarr, 'cells', segmentation)

    if export_metadata is not None:
        meta = cells.get_cell_metadata(seg_group)
        meta.write_csv(export_metadata)
    if import_metadata is not None:
        cells.apply_cell_metadata(seg_group, import_metadata)


@_viewer_metadata_command
def image_metadata(
    viewer_zarr,
    *,
    import_metadata: str | None = None,
    export_metadata: str | None = None,
) -> None:
    from ..modules.viewer import images
    from ..modules.viewer import zarr_utils as zu

    img_group = zu.get_viewer_group(viewer_zarr, 'images', 'multiplex')

    if export_metadata is not None:
        meta = images.get_channel_metadata(img_group)
        meta.write_csv(export_metadata)
    if import_metadata is not None:
        images.apply_channel_metadata(img_group, import_metadata)


@_viewer_metadata_command
def transcript_metadata(
    viewer_zarr,
    *,
    import_metadata: str | None = None,
    export_metadata: str | None = None,
) -> None:
    from ..modules.viewer import transcripts
    from ..modules.viewer import zarr_utils as zu

    tx_group = zu.get_viewer_group(viewer_zarr, 'transcripts')

    if export_metadata is not None:
        meta = transcripts.get_tx_metadata(tx_group)
        meta.write_csv(export_metadata)
    if import_metadata is not None:
        transcripts.apply_tx_metadata(tx_group, import_metadata)
