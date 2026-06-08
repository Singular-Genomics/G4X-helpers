import logging
import math
import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import polars as pl
from tqdm import tqdm

from .. import constants as c
from .. import io
from .. import utils as ut

log = logging.getLogger(__name__)


# region main function
def demux_raw_features(
    raw_features: pl.LazyFrame,
    manifest: pl.DataFrame,
    *,
    max_ham_dist: int = 2,
    min_delta: int = 2,
    demux_length: int = 15,
    batch_size: int = c.DEFAULT_BATCH_SIZE,
    batch_dir: Path | None = None,
    show_progress: bool = False,
):

    log.info('Starting batched demuxing of raw features')

    if batch_dir is None:
        log.debug('Creating temporary directory for demux batches')
        batch_dir = io.pathval.validate_dir_path(tempfile.mkdtemp(prefix='g4x_demux_tmp_'))

    else:
        batch_dir = io.pathval.validate_dir_path(batch_dir)
        batch_dir = io.pathval.ensure_dir(batch_dir / 'demux_batches')

    ut.log_with_path('Directory for temporary demux batches:', batch_dir)
    try:
        probe_dict = _build_probe_id_to_gene_name(manifest)
        seq_reads, manifest_by_read = _group_manifest_by_read(manifest)

        num_features = raw_features.select(pl.len()).collect().item()
        num_expected_batches = math.ceil(num_features / batch_size)
        next_progress_pct = 10

        lut = _build_base_lut()

        for i, feature_batch in tqdm(
            enumerate(_iter_feature_batches(raw_features, batch_size)),
            total=num_expected_batches,
            desc='Demuxing transcripts',
            position=0,
            disable=not show_progress,
        ):
            feature_batch = feature_batch.with_columns(
                pl.col('TXUID').str.split('_').list.last().cast(int).alias('read_num')
            )
            redemuxed_feature_batch = []
            for seq_read in seq_reads:
                feature_batch_read = _demux_feature_batch(
                    feature_batch=feature_batch,
                    seq_read=seq_read,
                    manifest_by_read=manifest_by_read,
                    probe_dict=probe_dict,
                    lut=lut,
                    demux_length=demux_length,
                    batch_size=batch_size,
                    max_ham_dist=max_ham_dist,
                    min_delta=min_delta,
                )

                redemuxed_feature_batch.append(feature_batch_read)

            demuxed_batch = pl.concat(redemuxed_feature_batch)
            demuxed_batch.write_parquet(batch_dir / f'batch_{i}.parquet')

            if num_expected_batches > 1:
                pct_complete = ((i + 1) * 100) // num_expected_batches
                while pct_complete >= next_progress_pct:
                    log.debug('Demuxing progress: %d%% (%d/%d batches)', next_progress_pct, i + 1, num_expected_batches)
                    next_progress_pct += 10

        return _compile_demuxed_batches(batch_dir)

    finally:
        if batch_dir.exists():
            log.debug('Removing temporary demux-batch directory')
            shutil.rmtree(batch_dir)


# region private functions
def _demux_feature_batch(
    feature_batch: pl.DataFrame,
    *,
    seq_read: int,
    manifest_by_read: dict[int, pl.DataFrame],
    probe_dict: dict[str, str],
    lut: np.ndarray,
    demux_length: int,
    batch_size: int,
    max_ham_dist: int,
    min_delta: int,
) -> pl.DataFrame:
    feature_batch_read = feature_batch.filter(pl.col('read_num') == seq_read)
    manifest_read = manifest_by_read[seq_read]

    if len(feature_batch_read) == 0 or len(manifest_read) == 0:
        return _mark_as_undemuxed(feature_batch_read)

    seqs = feature_batch_read['sequence'].to_list()
    codes = manifest_read['sequence'].to_list()
    seqs = [seq[:demux_length] for seq in seqs]
    codes = [seq[:demux_length] for seq in codes]

    codebook_target_ids = np.array(manifest_read['probe_id'].to_list())

    hammings = _compute_hamming_distance_matrix(seqs, codes, lut=lut, batch_size=batch_size)
    feature_batch_read = _assign_probe_matches(
        hammings=hammings,
        reads=feature_batch_read,
        codebook_target_ids=codebook_target_ids,
        probe_dict=probe_dict,
        max_ham_dist=max_ham_dist,
        min_delta=min_delta,
    )
    feature_batch_read = feature_batch_read.drop(['sequence', 'read_num'])
    return feature_batch_read


def _mark_as_undemuxed(feature_batch_read: pl.DataFrame) -> pl.DataFrame:
    return feature_batch_read.with_columns(
        [
            pl.lit('UNDETERMINED').alias('probe_name'),
            pl.lit('UNDETERMINED').alias(c.GENE_ID_NAME),
            pl.lit(False).alias('demuxed'),
        ]
    ).drop(['sequence', 'read_num'])


def _build_probe_id_to_gene_name(manifest: pl.DataFrame) -> dict[str, str]:
    mapping = dict(zip(manifest['probe_id'].to_list(), manifest['gene_name'].to_list()))
    mapping['UNDETERMINED'] = 'UNDETERMINED'
    return mapping


def _group_manifest_by_read(manifest: pl.DataFrame) -> tuple[list[int], dict[int, pl.DataFrame]]:
    seq_reads = manifest['read_num'].unique().to_list()
    seq_reads = [int(x.split('_')[-1]) if isinstance(x, str) else x for x in seq_reads]
    manifest_by_read = {read: manifest.filter(pl.col('read_num') == read) for read in seq_reads}
    return seq_reads, manifest_by_read


def _build_base_lut() -> np.ndarray:
    lut = np.zeros((256, 4), dtype=np.float32)
    for base, idx in zip(c.BASE_ORDER, range(4)):
        lut[ord(base), idx] = 1.0
    return lut


def _iter_feature_batches(
    raw_features: pl.LazyFrame, batch_size: int = c.DEFAULT_BATCH_SIZE, columns: str | list[str] | None = None
) -> Iterator[pl.DataFrame]:

    if columns:
        raw_features = raw_features.select(columns)
    offset = 0
    while True:
        batch = raw_features.slice(offset, batch_size).collect()
        if batch.is_empty():
            break
        yield batch
        offset += batch_size


def _compile_demuxed_batches(batch_dir: Path) -> pl.DataFrame:
    batch_paths = sorted(batch_dir.glob('batch_*.parquet'), key=lambda path: int(path.stem.split('_')[-1]))
    tx_table = pl.scan_parquet(batch_paths)
    return tx_table.filter(pl.col('demuxed')).drop('demuxed').collect()


def _assign_probe_matches(
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


def _compute_hamming_distance_matrix(
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
    codebook_oh = _one_hot_encode_sequences(codebook, seq_len, lut)
    M = len(codebook)

    # Prepare final result
    N = len(reads)
    hamming_matrix = np.empty((N, M), dtype=np.uint8)

    for i in range(0, N, batch_size):
        batch_reads = reads[i : i + batch_size]
        batch_oh = _one_hot_encode_sequences(batch_reads, seq_len, lut)
        matches = batch_oh @ codebook_oh.T
        hamming = seq_len - matches
        hamming_matrix[i : i + len(batch_reads)] = hamming

    return hamming_matrix


def _one_hot_encode_sequences(seqs: list[str], seq_len: int, lut: np.ndarray) -> np.ndarray:
    """
    Fast one-hot encoding using LUT.
    Returns: (N, seq_len * 4) float32 array
    """
    N = len(seqs)
    # Flatten all sequences into a byte array and reshape to (N, seq_len)
    arr = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8).reshape(N, seq_len)
    # Apply LUT: arr → (N, seq_len, 4), then flatten
    return lut[arr].reshape(N, seq_len * 4)
