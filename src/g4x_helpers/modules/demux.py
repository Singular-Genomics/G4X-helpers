import logging
import math
import shutil
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from tqdm import tqdm

from .. import c, io
from .. import logging_utils as logut
from ..schema.definition import Manifest, TxTable
from .workflow import PRESET_SOURCE, collect_input, reroute_source

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)


# region main function
# @g4x_workflow
def demux_raw_features(
    smp: 'G4Xoutput',
    manifest: str = PRESET_SOURCE,
    *,
    out_dir: str = PRESET_SOURCE,
    batch_size: int = c.DEFAULT_BATCH_SIZE,
    overwrite: bool = False,
    show_progress: bool | None = None,
    logger: logging.Logger | None = None,
) -> None:

    log = logger or LOGGER

    # 1: Validate and collect input
    log.debug('Validating input and preparing output paths')
    manifest_in = collect_input(smp, manifest, Manifest, logger=log)

    # 2: Validate and prepare output
    out_dir = smp.data_dir if out_dir == PRESET_SOURCE else io.pathval.validate_dir_path(out_dir)

    overwrite_manifest = True if manifest == PRESET_SOURCE else overwrite
    reroute_source(smp, out_dir, validator=Manifest, overwrite=overwrite_manifest, logger=log)
    reroute_source(smp, out_dir, validator=TxTable, overwrite=overwrite, logger=log)

    # if we're using a provided manifest, copy it into the output tree
    if manifest != PRESET_SOURCE:
        log.debug('Copying manifest from provided path into output tree')
        shutil.copy(manifest_in.p, smp.out.Manifest.p)

    # 3: Do the demuxing
    log.info('Starting batched demuxing of raw features')

    manifest = manifest_in.parse()
    show_progress = sys.stderr.isatty() if show_progress is None else show_progress
    try:
        batch_dir = io.pathval.ensure_dir(out_dir / 'demux_batches')
        batched_demuxing(
            feature_table_path=smp.src.RawFeatures.p,
            manifest=manifest,
            batch_dir=batch_dir,
            batch_size=batch_size,
            show_progress=show_progress,
            logger=log,
        )

        # 4: Compile the demuxed transcript table
        logut.log_with_path('Compiling demuxed transcript table from batch-dir:', batch_dir, logger=log)
        tx_table = pl.scan_parquet(list(batch_dir.glob('*.parquet')))
        tx_table = tx_table.filter(pl.col('demuxed')).drop('demuxed')

        # 5: Write the demuxed transcript table
        logut.log_with_path(f'Writing {smp.out.TxTable.name} table:', smp.out.TxTable.p, level='info', logger=log)
        tx_table.sink_csv(smp.out.TxTable.p, compression='gzip')

        # make sure that the gene list is populated with new genes for downstream steps
        smp.set_genes()

    # 6: Remove temporary demux batches
    finally:
        if batch_dir.exists():
            log.debug('Removing temporary demux-batch directory')
            shutil.rmtree(batch_dir)


def batched_demuxing(
    feature_table_path: str,
    manifest: pl.DataFrame,
    batch_dir: str,
    batch_size: int = c.DEFAULT_BATCH_SIZE,
    show_progress: bool | None = None,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    probe_dict = dict(zip(manifest['probe'].to_list(), manifest['gene_name'].to_list()))
    probe_dict['UNDETERMINED'] = 'UNDETERMINED'

    seq_reads = manifest['read'].unique().to_list()
    seq_reads = [int(x.split('_')[-1]) if isinstance(x, str) else x for x in seq_reads]
    manifest_by_read = {seq_read: manifest.filter(pl.col('read') == seq_read) for seq_read in seq_reads}

    num_features = pl.scan_parquet(feature_table_path).select(pl.len()).collect().item()
    num_expected_batches = math.ceil(num_features / batch_size)
    next_progress_pct = 10

    if num_features == 0:
        log.info('Demuxing progress: 100%% (no features to process)')
        return

    LUT = np.zeros((256, 4), dtype=np.float32)
    for base, idx in zip(c.BASE_ORDER, range(4)):
        LUT[ord(base), idx] = 1.0

    for i, feature_batch in tqdm(
        enumerate(stream_features(feature_table_path, batch_size)),
        total=num_expected_batches,
        desc='Demuxing transcripts',
        position=0,
        disable=not show_progress,
    ):
        feature_batch = feature_batch.with_columns(pl.col('TXUID').str.split('_').list.last().cast(int).alias('read'))
        redemuxed_feature_batch = []
        for seq_read in seq_reads:
            feature_batch_read = feature_batch.filter(pl.col('read') == seq_read)
            manifest_read = manifest_by_read[seq_read]

            if len(feature_batch_read) == 0 or len(manifest_read) == 0:
                continue

            seqs = feature_batch_read['sequence'].to_list()
            codes = manifest_read['sequence'].to_list()
            codebook_target_ids = np.array(manifest_read['probe'].to_list())

            hammings = batched_dot_product_hamming_matrix(seqs, codes, lut=LUT, batch_size=batch_size)
            feature_batch_read = demux(hammings, feature_batch_read, codebook_target_ids, probe_dict)
            feature_batch_read = feature_batch_read.drop(['sequence', 'read'])
            redemuxed_feature_batch.append(feature_batch_read)

        batch_dir = Path(batch_dir)
        pl.concat(redemuxed_feature_batch).write_parquet(batch_dir / f'batch_{i}.parquet')

        pct_complete = ((i + 1) * 100) // num_expected_batches
        while pct_complete >= next_progress_pct:
            log.debug('Demuxing progress: %d%% (%d/%d batches)', next_progress_pct, i + 1, num_expected_batches)
            next_progress_pct += 10


def stream_features(
    feature_table_path: str, batch_size: int = c.DEFAULT_BATCH_SIZE, columns: str | list[str] | None = None
) -> Iterator[pl.DataFrame]:
    df = pl.scan_parquet(feature_table_path)
    if columns:
        df = df.select(columns)
    offset = 0
    while True:
        batch = df.slice(offset, batch_size).collect()
        if batch.is_empty():
            break
        yield batch
        offset += batch_size


def demux(
    hammings: np.ndarray,
    reads: pl.DataFrame,
    codebook_target_ids: np.ndarray,
    probe_dict: dict,
    max_ham_dist: int = 2,
    min_delta: int = 2,
) -> pl.DataFrame:
    demuxed = np.zeros(hammings.shape[0], dtype=bool)

    for i in range(max_ham_dist + 1):
        hits = hammings == i
        close_hits = hammings <= (i + min_delta)
        uniquely_hit = hits.sum(axis=1) == 1
        close_hit = close_hits.sum(axis=1) > 1
        pass_filter = uniquely_hit & ~close_hit
        demuxed[pass_filter] = 1

        # logger.info(f"""
        # ... ... {fq.stem}
        # hamming == {i}, min_delta == {min_delta}
        # unique hits = {sum(uniquely_hit)}
        # total cumulative hits within min_delta = {sum(close_hit)}
        # total demuxed (unique hits without another hit within min_delta) = {sum(pass_filter)}
        # """)

    # --- Get best hits ---
    hit_ids = hammings.argmin(axis=1)
    hit_targets = codebook_target_ids[hit_ids]

    transcripts = np.where(demuxed, hit_targets, 'UNDETERMINED')
    transcript_condensed = [probe_dict.get(t, 'UNDETERMINED') for t in transcripts]

    reads = reads.with_columns(
        [
            pl.Series('probe_name', transcripts),
            pl.Series(c.GENE_ID_NAME, transcript_condensed),
            pl.Series('demuxed', demuxed),
        ]
    )

    return reads


def batched_dot_product_hamming_matrix(
    reads: list[str],
    codebook: list[str],
    lut: np.ndarray,
    batch_size: int,
) -> np.ndarray:
    """
    Compute full Hamming distance matrix (N_reads, N_codebook)
    using batched dot-product with one-hot encoding.
    """
    seq_len = len(codebook[0])
    assert all(len(seq) == seq_len for seq in codebook), 'All codebook entries must be same length'

    # One-hot encode the codebook once
    codebook_oh = one_hot_encode_str_array(codebook, seq_len, lut)
    M = len(codebook)

    # Prepare final result
    N = len(reads)
    hamming_matrix = np.empty((N, M), dtype=np.uint8)

    # num_expected_batches = math.ceil(N / batch_size)
    # for i in tqdm(range(0, N, batch_size), total=num_expected_batches, desc='Demuxing batch', position=1, leave=False):
    for i in range(0, N, batch_size):  # , total=num_expected_batches, desc='Demuxing batch', position=1, leave=False):
        batch_reads = reads[i : i + batch_size]
        batch_oh = one_hot_encode_str_array(batch_reads, seq_len, lut)
        matches = batch_oh @ codebook_oh.T
        hamming = seq_len - matches
        hamming_matrix[i : i + len(batch_reads)] = hamming

    return hamming_matrix


def one_hot_encode_str_array(seqs: list[str], seq_len: int, lut: np.ndarray) -> np.ndarray:
    """
    Fast one-hot encoding using LUT.
    Returns: (N, seq_len * 4) float32 array
    """
    N = len(seqs)
    # Flatten all sequences into a byte array and reshape to (N, seq_len)
    arr = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8).reshape(N, seq_len)
    # Apply LUT: arr → (N, seq_len, 4), then flatten
    return lut[arr].reshape(N, seq_len * 4)


# TODO consider re-implementing this function if needed
# def update_metadata_and_tx_file(g4x_obj: 'G4Xoutput', manifest, out_dir):
#     if not manifest == out_dir / 'transcript_panel.csv':
#         shutil.copy(manifest, out_dir / 'transcript_panel.csv')

#     panel_name = manifest.name
#     timestamp = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
#     ## add info to run_meta.json
#     with open(g4x_obj.data_dir / 'run_meta.json', 'r') as f:
#         meta = json.load(f)

#     meta['transcript_panel'] = panel_name
#     meta['redemuxed_timestamp'] = timestamp

#     with open(out_dir / 'run_meta.json', 'w') as f:
#         json.dump(meta, f, indent=2)

#     ## add info to run_meta.json in g4x_viewer
#     with open(g4x_obj.data_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_run_metadata.json', 'r') as f:
#         meta = json.load(f)

#     meta['run_metadata']['transcript_panel'] = panel_name
#     meta['run_metadata']['redemuxed_timestamp'] = timestamp

#     with open(out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_run_metadata.json', 'w') as f:
#         json.dump(meta, f, indent=2)
