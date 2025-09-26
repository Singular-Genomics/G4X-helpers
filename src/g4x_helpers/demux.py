import re
import math
import numpy as np
import polars as pl
from tqdm import tqdm
from pathlib import Path

BASE_ORDER = 'CTGA'
LUT = np.zeros((256, 4), dtype=np.float32)
LUT[ord('C'), 0] = 1.0
LUT[ord('T'), 1] = 1.0
LUT[ord('G'), 2] = 1.0
LUT[ord('A'), 3] = 1.0


def one_hot_encode_str_array(seqs: list[str], seq_len: int) -> np.ndarray:
    """
    Fast one-hot encoding using LUT.
    Returns: (N, seq_len * 4) float32 array
    """
    N = len(seqs)
    # Flatten all sequences into a byte array and reshape to (N, seq_len)
    arr = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8).reshape(N, seq_len)
    # Apply LUT: arr → (N, seq_len, 4), then flatten
    return LUT[arr].reshape(N, seq_len * 4)


def batched_dot_product_hamming_matrix(
    reads: list[str],
    codebook: list[str],
    batch_size: int,
) -> np.ndarray:
    """
    Compute full Hamming distance matrix (N_reads, N_codebook)
    using batched dot-product with one-hot encoding.
    """
    seq_len = len(codebook[0])
    assert all(len(seq) == seq_len for seq in codebook), 'All codebook entries must be same length'

    # One-hot encode the codebook once
    codebook_oh = one_hot_encode_str_array(codebook, seq_len)
    M = len(codebook)

    # Prepare final result
    N = len(reads)
    hamming_matrix = np.empty((N, M), dtype=np.uint8)

    num_expected_batches = math.ceil(N / batch_size)
    for i in tqdm(range(0, N, batch_size), total=num_expected_batches, desc='Demuxing batch', position=1, leave=False):
        batch_reads = reads[i : i + batch_size]
        batch_oh = one_hot_encode_str_array(batch_reads, seq_len)
        matches = batch_oh @ codebook_oh.T
        hamming = seq_len - matches
        hamming_matrix[i : i + len(batch_reads)] = hamming

    return hamming_matrix


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
            pl.Series('transcript_new', transcripts),
            pl.Series('transcript_condensed_new', transcript_condensed),
            pl.Series('demuxed_new', demuxed),
        ]
    )

    return reads


def load_manifest(file_path: str | Path) -> tuple[pl.DataFrame, dict]:
    manifest = pl.read_csv(file_path)
    target_list = manifest['target'].to_list()

    p = re.compile('-[ACGT]{2,30}')
    match = re.findall(pattern=p, string=target_list[0])
    if len(match) == 0:
        probe_dict = {p: '-'.join(p.split('-')[:-1]) for p in target_list}
    else:
        probe_dict = {}
        for probe in target_list:
            match = re.findall(pattern=p, string=probe)
            probe_dict[probe] = probe.split(match[0])[0]
    probe_dict['UNDETERMINED'] = 'UNDETERMINED'

    return manifest, probe_dict
