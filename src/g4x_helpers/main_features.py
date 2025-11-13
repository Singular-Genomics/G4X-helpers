import functools
import logging
from typing import TYPE_CHECKING

import rich_click as click

from . import __version__, utils
from . import workflows as wflw

if TYPE_CHECKING:
    from .models import G4Xoutput


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        g4x_obj: 'G4Xoutput',
        out_dir: str,
        n_threads: int = utils.DEFAULT_THREADS,
        verbose: int = 1,
        file_logger: bool = True,
        **kwargs,
    ):
        func_name = func.__name__

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
    from .modules import segmentation as seg

    cell_labels = utils.validate_path(cell_labels, must_exist=True, is_dir_ok=False, is_file_ok=True)

    labels = seg.try_load_segmentation(
        cell_labels=cell_labels,
        expected_shape=g4x_obj.shape,
        labels_key=labels_key,
    )

    adata, mask = wflw.intersect_segmentation(
        g4x_obj=g4x_obj,
        labels=labels,
        out_dir=out_dir,
        logger=logger,
    )

    # TODO this pattern is repeated in several places, refactor into utility function
    view_dir = out_dir / 'g4x_viewer'
    view_dir.mkdir(parents=True, exist_ok=True)

    wflw.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=view_dir / f'{g4x_obj.sample_id}.bin',
        protein_list=[f'{x}_intensity_mean' for x in g4x_obj.proteins],
        n_threads=n_threads,
        logger=logger,
    )


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
    wflw.redemux(
        g4x_obj=g4x_obj,
        manifest=manifest,
        out_dir=out_dir,
        batch_size=batch_size,
        logger=logger,
    )

    ## now regenerate the secondary files
    logger.info('Regenerating downstream files.')
    g4x_obj.data_dir = out_dir  ## set base directory to the redemux output folder for downstream steps

    # resegment with existing segmentation
    logger.info('Intersecting with existing cell segmentation.')
    adata, mask = wflw.intersect_segmentation(
        g4x_obj=g4x_obj,
        labels=g4x_obj.load_segmentation(),
        out_dir=out_dir,
        logger=logger,
    )

    view_dir = out_dir / 'g4x_viewer'
    view_dir.mkdir(parents=True, exist_ok=True)

    wflw.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=view_dir / f'{g4x_obj.sample_id}.bin',
        protein_list=[f'{x}_intensity_mean' for x in g4x_obj.proteins],
        n_threads=n_threads,
        logger=logger,
    )

    wflw.tx_converter(
        g4x_obj=g4x_obj,
        out_path=view_dir / f'{g4x_obj.sample_id}.tar',
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
    cluster_color_key: str | None = None,
    emb_key: str | None = None,
    logger: logging.Logger,
    **kwargs,
) -> None:
    metadata = utils.validate_path(metadata, must_exist=True, is_dir_ok=False, is_file_ok=True)

    view_dir = out_dir / 'g4x_viewer'
    view_dir.mkdir(parents=True, exist_ok=True)

    wflw.seg_updater(
        bin_file=g4x_obj.data_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}.bin' if bin_file is None else bin_file,
        metadata_file=metadata,
        out_path=view_dir / f'{g4x_obj.sample_id}.bin',
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        cluster_color_key=cluster_color_key,
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
    ## set up the data
    try:
        adata = g4x_obj.load_adata()
        emb_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
        cluster_key = sorted([x for x in adata.obs.columns if 'leiden' in x])[0]
        logger.info('Successfully loaded adata with clustering information.')

    except Exception:
        adata = g4x_obj.load_adata(load_clustering=False)
        emb_key = None
        cluster_key = None
        logger.info('Clustering information was not found, cell coloring will be random.')

    view_dir = out_dir / 'g4x_viewer'
    view_dir.mkdir(parents=True, exist_ok=True)

    wflw.seg_converter(
        adata=adata,
        seg_mask=g4x_obj.load_segmentation(),
        out_path=view_dir / f'{g4x_obj.sample_id}.bin',
        cluster_key=cluster_key,
        emb_key=emb_key,
        protein_list=[f'{x}_intensity_mean' for x in g4x_obj.proteins],
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
    viewer_dir = g4x_obj.data_dir / 'g4x_viewer'
    wflw.tar_viewer(viewer_dir=viewer_dir, out_dir=out_dir, logger=logger)
