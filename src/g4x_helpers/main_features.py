import functools
import logging
from typing import TYPE_CHECKING

import rich_click as click

from . import __version__, utils

if TYPE_CHECKING:
    from .models import G4Xoutput


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        g4x_obj: 'G4Xoutput',
        out_dir: str | None = None,
        n_threads: int = utils.DEFAULT_THREADS,
        verbose: int = 1,
        file_logger: bool = True,
        **kwargs,
    ):
        func_name = func.__name__

        if func_name != 'migrate':
            g4x_obj.validate()

        out_dir = g4x_obj.data_dir if out_dir is None else out_dir
        out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

        if out_dir != g4x_obj.data_dir:
            func_out = out_dir / func_name
            func_out.mkdir(parents=True, exist_ok=True)
            out_dir = func_out

        gap = 12
        click.secho(f'\nStarting: {func_name}\n', bold=True)
        utils.print_k_v('sample_dir', f'{g4x_obj.data_dir}', gap)
        utils.print_k_v('out_dir', f'{out_dir}', gap)
        utils.print_k_v('n_threads', f'{n_threads}', gap)
        utils.print_k_v('verbosity', f'{verbose}', gap)
        utils.print_k_v('g4x-helpers', f'v{__version__}', gap)
        click.echo('')

        if not kwargs.get('logger', None):
            logger = utils.setup_logger(
                logger_name=func_name, out_dir=out_dir / 'logs', stream_level=verbose, file_logger=file_logger
            )
        logger.info(f'Running {func_name} with G4X-helpers v{__version__}')

        # Pass the initialized logger and parameters to the wrapped function
        result = func(
            g4x_obj=g4x_obj,
            out_dir=out_dir,
            n_threads=n_threads,
            verbose=verbose,
            logger=logger,
            **kwargs,
        )

        click.secho(f'\nCompleted: {func.__name__}', bold=True, fg='green')
        return result

    return wrapper


@_base_command
def resegment(
    g4x_obj: 'G4Xoutput',
    out_dir: str,
    cell_labels: str,
    *,
    labels_key: str | None = None,
    n_threads: int = 4,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .modules import edit_bin, init_bin, segment

    cell_labels = utils.validate_path(cell_labels, must_exist=True, is_dir_ok=False, is_file_ok=True)

    logger.info('Loading segmentation mask.')
    labels = segment.try_load_segmentation(cell_labels=cell_labels, expected_shape=g4x_obj.shape, labels_key=labels_key)

    segment.apply_segmentation(
        g4x_obj=g4x_obj,
        labels=labels,
        out_dir=out_dir,
        skip_protein_extraction=False,  ### TODO <<< RMOVE THIS HACK LATER
        logger=logger,
    )

    init_bin.init_bin_file(g4x_obj=g4x_obj, out_dir=out_dir, logger=logger, n_threads=n_threads)

    bin_file = out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin'
    edit_bin.edit_bin_file(g4x_obj=g4x_obj, bin_file=bin_file, logger=logger)


@_base_command
def redemux(
    g4x_obj: 'G4Xoutput',
    out_dir: str,
    manifest: str,
    *,
    batch_size: int = 1_000_000,
    n_threads: int = 4,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from g4x_helpers.modules import demux, edit_bin, init_bin, segment, transcript_tar

    demux.demux_raw_features(
        g4x_obj=g4x_obj,
        manifest=manifest,
        out_dir=out_dir,
        batch_size=batch_size,
        logger=logger,
    )

    labels = g4x_obj.load_segmentation()

    segment.apply_segmentation(
        g4x_obj=g4x_obj,
        labels=labels,
        out_dir=out_dir,
        skip_protein_extraction=True,
        create_source=True,
        logger=logger,
    )

    init_bin.init_bin_file(g4x_obj=g4x_obj, out_dir=out_dir, n_threads=n_threads, logger=logger)

    bin_file = out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_segmentation.bin'
    edit_bin.edit_bin_file(g4x_obj=g4x_obj, bin_file=bin_file, logger=logger)

    transcript_tar.create_tx_tarfile(
        g4x_obj=g4x_obj,
        out_path=out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_transcripts.tar',
        n_threads=n_threads,
        logger=logger,
    )


@_base_command
def update_bin(
    g4x_obj: 'G4Xoutput',
    out_dir: str,
    metadata: str,
    *,
    bin_file: str | None = None,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    clustercolor_key: str | None = None,
    emb_key: str | None = None,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .modules import edit_bin

    metadata = utils.validate_path(metadata, must_exist=True, is_dir_ok=False, is_file_ok=True)

    edit_bin.edit_bin_file(
        g4x_obj=g4x_obj,
        bin_file=bin_file,
        metadata=metadata,
        bin_out=None,
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        clustercolor_key=clustercolor_key,
        emb_key=emb_key,
        logger=logger,
    )


@_base_command
def new_bin(
    g4x_obj: 'G4Xoutput',  #
    out_dir: str,
    *,
    n_threads: int = 4,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .modules import init_bin

    init_bin.init_bin_file(
        g4x_obj=g4x_obj,
        out_dir=out_dir,
        n_threads=n_threads,
        logger=logger,
    )


@_base_command
def tar_viewer(
    g4x_obj: 'G4Xoutput',  #
    out_dir: str | None = None,
    *,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .modules import viewer_dir

    viewer_dir.package_viewer_dir(
        viewer_dir=g4x_obj.data_dir / 'g4x_viewer',
        sample_id=g4x_obj.sample_id,
        out_dir=out_dir,
        logger=logger,
    )


@_base_command
def migrate(
    g4x_obj: 'G4Xoutput',
    restore: bool = False,
    *,
    logger: logging.Logger,
    **kwargs,
) -> None:
    from .schemas import migration

    if restore:
        migration.restore_backup(sample_base=g4x_obj.data_dir, sample_id=g4x_obj.sample_id, logger=logger)
    else:
        migration.migrate_g4x_data(data_dir=g4x_obj.data_dir, sample_id=g4x_obj.sample_id, logger=logger)
