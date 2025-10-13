import functools
from pathlib import Path
from typing import TYPE_CHECKING

import rich_click as click

from . import __version__, utils
from . import demux as dmx
from . import segmentation as seg
from . import tar_viewer as tv
from .g4x_viewer import bin_generator as bg

if TYPE_CHECKING:
    from .models import G4Xoutput


def _base_command(func):
    """Decorator to apply standard command initialization logic."""

    @functools.wraps(func)
    def wrapper(
        g4x_out: 'G4Xoutput',
        out_dir: Path | str,
        n_threads: int = 4,
        verbose: int = 1,
        **kwargs,
    ):
        out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)
        func_out = out_dir / func.__name__
        func_out.mkdir(parents=True, exist_ok=True)
        out_dir = func_out

        logger = utils.setup_logger(logger_name=func.__name__, out_dir=out_dir, stream_level=verbose)

        gap = 11
        click.secho(f'\nStarting: {func.__name__}', bold=True)
        utils.print_k_v('sample_dir', f'{g4x_out.run_base}', gap)
        utils.print_k_v('out_dir', f'{out_dir}', gap)
        utils.print_k_v('n_threads', f'{n_threads}', gap)
        utils.print_k_v('verbosity', f'{verbose}', gap)
        utils.print_k_v('version', f'{__version__}', gap)

        # Pass the initialized logger and parameters to the wrapped function
        result = func(
            g4x_out=g4x_out,
            out_dir=out_dir,
            n_threads=n_threads,
            verbose=verbose,
            logger=logger,
            **kwargs,
        )

        click.secho(f'Completed: {func.__name__}\n', bold=True, fg='green')
        return result

    return wrapper


@_base_command
def resegment(
    g4x_out: 'G4Xoutput',
    segmentation_mask: str,
    out_dir: str,
    *,
    segmentation_mask_key: str | None = None,
    n_threads: int = 4,
    verbose: int = 1,
    **kwargs,
) -> None:
    logger = kwargs['logger']
    segmentation_mask = utils.validate_path(segmentation_mask, must_exist=True, is_dir_ok=False, is_file_ok=True)

    labels = seg.try_load_segmentation(
        segmentation_mask=segmentation_mask,
        expected_shape=g4x_out.shape,
        segmentation_mask_key=segmentation_mask_key,
    )

    adata, mask = seg.intersect_segmentation(
        g4x_out=g4x_out,
        labels=labels,
        out_dir=out_dir,
        n_threads=n_threads,
        logger=logger,
    )

    bg.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=bg.bin_file_path(g4x_out, out_dir=out_dir),
        protein_list=[f'{x}_intensity_mean' for x in g4x_out.proteins],
        n_threads=n_threads,
        logger=logger,
    )


@_base_command
def redemux(
    g4x_out: 'G4Xoutput',
    manifest: Path | str,
    out_dir: Path | str,
    batch_size: int = 1_000_000,
    n_threads: int = 4,
    verbose: int = 1,
    **kwargs,
) -> None:
    logger = kwargs['logger']

    ## preflight checks
    manifest = utils.validate_path(manifest, must_exist=True, is_dir_ok=False, is_file_ok=True)
    # out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

    batch_dir = out_dir / 'diagnostics' / 'batches'
    batch_dir.mkdir(parents=True, exist_ok=True)

    ## make output directory with symlinked files from original
    utils.symlink_original_files(g4x_out, out_dir)

    ## update metadata and transcript panel file
    dmx.update_metadata_and_tx_file(g4x_out, manifest, out_dir)

    ## load the new manifest file that we will demux against
    logger.info('Loading manifest file.')
    manifest, probe_dict = dmx.load_manifest(manifest)

    logger.info('Performing re-demuxing.')
    dmx.batched_demuxing(g4x_out, manifest, probe_dict, batch_dir, batch_size)

    ## concatenate results into final csv and parquet
    logger.info('Writing updated transcript table.')
    dmx.concatenate_and_cleanup(batch_dir, out_dir)

    ## now regenerate the secondary files
    logger.info('Regenerating downstream files.')
    ## set run_base to the redemux output folder for downstream steps
    g4x_out.run_base = out_dir

    # resegment with existing segmentation
    logger.info('Intersecting with existing segmentation.')

    labels = g4x_out.load_segmentation()
    adata, mask = seg.intersect_segmentation(
        g4x_out, labels=labels, out_dir=out_dir, n_threads=n_threads, logger=logger
    )

    bg.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=out_dir / 'g4x_viewer' / f'{g4x_out.sample_id}.bin',
        protein_list=[f'{x}_intensity_mean' for x in g4x_out.proteins],
        n_threads=n_threads,
        logger=logger,
    )

    dmx.tx_converter(
        g4x_out,
        out_path=out_dir / 'g4x_viewer' / f'{g4x_out.sample_id}.tar',
        n_threads=n_threads,
        logger=logger,
    )


@_base_command
def update_bin(
    g4x_out: 'G4Xoutput',
    metadata: str,
    out_dir: str,
    *,
    bin_file: str | None = None,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    cluster_color_key: str | None = None,
    emb_key: str | None = None,
    verbose: int = 1,
    **kwargs,
) -> None:
    logger = kwargs['logger']

    metadata = utils.validate_path(metadata, must_exist=True, is_dir_ok=False, is_file_ok=True)
    
    if bin_file:
        bin_file = utils.validate_path(bin_file, must_exist=True, is_dir_ok=False, is_file_ok=True)
    else:
        g4x_out.run_base / 'g4x_viewer' / f'{g4x_out.sample_id}.bin'

    out_path = out_dir / f'{g4x_out.sample_id}.bin'

    bg.seg_updater(
        bin_file=bin_file,
        metadata_file=metadata,
        out_path=out_path,
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        cluster_color_key=cluster_color_key,
        emb_key=emb_key,
        logger=logger,
    )


@_base_command
def new_bin(
    g4x_out: 'G4Xoutput',  #
    out_dir: str,
    n_threads: int = 4,
    verbose: int = 1,
    **kwargs,
) -> None:
    logger = kwargs['logger']

    ## set up the data
    try:
        adata = g4x_out.load_adata()
        emb_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
        cluster_key = sorted([x for x in adata.obs.columns if 'leiden' in x])[0]
        logger.info('Successfully loaded adata with clustering information.')
    except Exception:
        adata = g4x_out.load_adata(load_clustering=False)
        emb_key = None
        cluster_key = None
        logger.info('Clustering information was not found, cell coloring will be random.')

    mask = g4x_out.load_segmentation()

    _ = bg.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=bg.bin_file_path(g4x_out, out_dir=out_dir),
        cluster_key=cluster_key,
        emb_key=emb_key,
        protein_list=[f'{x}_intensity_mean' for x in g4x_out.proteins],
        n_threads=n_threads,
        logger=logger,
    )


@_base_command
def tar_viewer(
    g4x_out: 'G4Xoutput',  #
    viewer_dir: str | None = None,
    out_dir: str | None = None,
    **kwargs,
) -> None:
    logger = kwargs['logger']

    if viewer_dir is None:
        viewer_dir = g4x_out.run_base / 'g4x_viewer'

    viewer_dir = utils.validate_path(viewer_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)

    tv.tar_viewer(out_path=out_dir, viewer_dir=viewer_dir, logger=logger)
