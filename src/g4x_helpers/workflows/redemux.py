import logging
from pathlib import Path
from typing import TYPE_CHECKING

from .. import utils
from ..modules import demultiplexing as dmx

if TYPE_CHECKING:
    from ..models import G4Xoutput


def redemux(
    g4x_out: 'G4Xoutput',
    manifest: Path | str,
    out_dir: Path | str,
    *,
    batch_size: int = 1_000_000,
    logger: logging.Logger,
) -> None:
    ## preflight checks
    manifest = utils.validate_path(manifest, must_exist=True, is_dir_ok=False, is_file_ok=True)
    out_dir = utils.validate_path(out_dir, must_exist=False, is_dir_ok=True, is_file_ok=False)

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
