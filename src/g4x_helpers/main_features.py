from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .models import G4Xoutput

from . import demux as dmx
from . import segmentation as seg
from . import tar_viewer as tv
from . import utils
from .g4x_viewer import bin_generator as bg


def redemux(
    g4x_out: 'G4Xoutput',
    manifest: Path | str,
    out_dir: Path | str,
    batch_size: int = 1_000_000,
    n_threads: int = 4,
    verbose: int = 1,
):
    print('Starting redemux')
    logger = utils.setup_logger(logger_name='redemux', out_dir=out_dir, stream_level=verbose)

    ## preflight checks
    manifest = utils.validate_path(manifest, must_exist=True, is_dir_ok=False, is_file_ok=True)
    out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

    batch_dir = out_dir / 'diagnostics' / 'batches'
    batch_dir.mkdir(parents=True, exist_ok=True)

    ## make output directory with symlinked files from original
    dmx.symlink_original_files(g4x_out, out_dir)

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

    ## set run_base to the redemux output folder
    # TODO make sure I understand why this is needed
    g4x_out.run_base = out_dir

    ## now regenerate the secondary files
    logger.info('Regenerating downstream files.')

    # resegment with existing segmentation
    logger.info('Intersecting with existing segmentation.')

    labels = g4x_out.load_segmentation()
    _ = seg.intersect_segmentation(g4x_out, labels=labels, out_dir=out_dir, n_threads=n_threads, logger=logger)

    logger.info('Generating viewer transcript file.')
    _ = dmx.tx_converter(
        g4x_out,
        out_path=out_dir / 'g4x_viewer' / f'{g4x_out.sample_id}.tar',
        n_threads=n_threads,
        logger=logger,
    )

    logger.info('Completed redemux.')


def resegment(
    g4x_out: 'G4Xoutput',
    segmentation_mask: str,
    *,
    out_dir: Path | str,
    segmentation_mask_key: str | None = None,
    n_threads: int = 4,
    verbose: int = 1,
) -> None:
    print('Starting resegmentation')
    logger = utils.setup_logger(logger_name='resegment', out_dir=out_dir, stream_level=verbose)

    # load new segmentation
    labels = seg.try_load_segmentation(segmentation_mask, g4x_out.shape, segmentation_mask_key)

    # run intersection with new segmentation
    seg.intersect_segmentation(g4x_out=g4x_out, labels=labels, out_dir=out_dir, n_threads=n_threads, logger=logger)

    print('Resegmentation complete.')


def update_bin(
    bin_file: str,
    out_dir: str,
    metadata: str,
    *,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    cluster_color_key: str | None = None,
    emb_key: str | None = None,
    verbose: int = 1,
):
    print('Starting update_bin')
    logger = utils.setup_logger(logger_name='update_bin', out_dir=out_dir, stream_level=verbose)

    bin_file = utils.validate_path(bin_file, must_exist=True, is_dir_ok=False, is_file_ok=True)
    out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)
    metadata = utils.validate_path(metadata, must_exist=True, is_dir_ok=False, is_file_ok=True)

    out_dir = out_dir.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    _ = bg.seg_updater(
        bin_file=bin_file,
        metadata_file=metadata,
        out_path=out_dir,
        cellid_key=cellid_key,
        cluster_key=cluster_key,
        cluster_color_key=cluster_color_key,
        emb_key=emb_key,
        logger=logger,
    )


def new_bin(g4x_out: 'G4Xoutput', out_dir: str, n_threads: int = 4, verbose: int = 1) -> None:
    print('Starting new_bin')
    logger = utils.setup_logger(logger_name='new_bin', out_dir=out_dir, stream_level=verbose)

    out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

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

    # TODO need to make decision on whether to allow no-output dir
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        outfile = Path(out_dir) / f'{g4x_out.sample_id}.bin'
    else:
        outfile = g4x_out.run_base / 'g4x_viewer' / f'{g4x_out.sample_id}.bin'

    logger.info('Making G4X-Viewer bin file.')
    _ = bg.seg_converter(
        adata=adata,
        seg_mask=mask,
        out_path=outfile,
        cluster_key=cluster_key,
        emb_key=emb_key,
        protein_list=[f'{x}_intensity_mean' for x in g4x_out.proteins],
        n_threads=n_threads,
        logger=logger,
    )
    logger.debug(f'G4X-Viewer bin --> {outfile}')


def tar_viewer(g4x_out: 'G4Xoutput', out_dir: str, viewer_dir: str | None = None, verbose: int = 1):
    print('Starting tar_viewer')
    logger = utils.setup_logger(logger_name='tar_viewer', out_dir=out_dir, stream_level=verbose)

    out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)
    if viewer_dir is None:
        viewer_dir = g4x_out.run_base / 'g4x_viewer'
    else:
        viewer_dir = utils.validate_path(viewer_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)

    tv.tar_viewer(out_path=out_dir, viewer_dir=viewer_dir, logger=logger)
