import json
import math

# import re
import shutil
from datetime import datetime

# from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from tqdm import tqdm

from .. import utils

if TYPE_CHECKING:
    from ..models import G4Xoutput

BASE_ORDER = 'CTGA'
LUT = np.zeros((256, 4), dtype=np.float32)
LUT[ord('C'), 0] = 1.0
LUT[ord('T'), 1] = 1.0
LUT[ord('G'), 2] = 1.0
LUT[ord('A'), 3] = 1.0

primer_read_map = {
    'SP1': 1,
    'm7a': 2,
    'm9a': 3,
    'm6a': 4,
    'm8a': 5,
    'm3a': 6,
}

PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'


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


# def load_manifest(file_path: str | Path) -> tuple[pl.DataFrame, dict]:
#     manifest = pl.read_csv(file_path)
#     target_list = manifest['probe_name'].to_list()

#     p = re.compile('-[ACGT]{2,30}')
#     match = re.findall(pattern=p, string=target_list[0])
#     if len(match) == 0:
#         probe_dict = {p: '-'.join(p.split('-')[:-1]) for p in target_list}
#     else:
#         probe_dict = {}
#         for probe in target_list:
#             match = re.findall(pattern=p, string=probe)
#             probe_dict[probe] = probe.split(match[0])[0]
#     probe_dict['UNDETERMINED'] = 'UNDETERMINED'

#     return manifest, probe_dict


def parse_input_manifest(file_path: str) -> pl.DataFrame:
    """
    Parse strings in `df[col]` of the form '<prefix>-<15mer>-<primer>'
    and add columns: gene_name, sequence, primer, primer_code, read.
    Rows that don't match the pattern get nulls in the new columns.
    """
    df = pl.read_csv(file_path)

    col = 'probe_name'
    if col not in df.columns:
        raise ValueError(f"transcript manifest must contain column '{col}'")

    parsed = df.with_columns(
        [
            pl.col(col).str.extract(PROBE_PATTERN, 1).alias('gene_name'),
            pl.col(col).str.extract(PROBE_PATTERN, 2).alias('sequence'),
            pl.col(col).str.extract(PROBE_PATTERN, 3).alias('primer'),
        ]
    )

    null_count = parsed.null_count()['sequence'][0]
    if null_count > 0:
        print(f'{null_count} probes with invalid sequence format will be ignored:')
        null_seqs = parsed.filter(pl.col('sequence').is_null())['probe_name'].to_list()
        for ns in null_seqs:
            print(f'- {ns}')

    if 'read' in df.columns:
        print("Using 'read' column provided by input manifest.")
        parsed = parsed.drop('read')
    else:
        plist = parsed['primer'].unique().to_list()
        ign_primer = [p for p in plist if p not in primer_read_map]
        if len(ign_primer) > 0:
            print('Warning: the following primer names are not known and will be ignored:')
            for ip in ign_primer:
                print(f'- {ip}')
        parsed = parsed.filter(pl.col('primer').is_in(primer_read_map.keys())).with_columns(
            pl.col('primer').replace(primer_read_map).cast(pl.Int8).alias('read')
        )

    if 'gene_name' in df.columns:
        print("Using 'gene_name' column provided by input manifest.")
        parsed = parsed.drop('gene_name')

    manifest = parsed.join(df, on=col, how='left')

    probe_dict = dict(zip(manifest['probe_name'].to_list(), manifest['gene_name'].to_list()))
    probe_dict['UNDETERMINED'] = 'UNDETERMINED'

    return manifest, probe_dict


def update_metadata_and_tx_file(g4x_obj: 'G4Xoutput', manifest, out_dir):
    ## update metadata and transcript panel file
    shutil.copy(manifest, out_dir / 'transcript_panel.csv')

    with open(out_dir / 'run_meta.json', 'r') as f:
        meta = json.load(f)

    meta['transcript_panel'] = manifest.name
    meta['redemuxed_timestamp'] = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
    (out_dir / 'run_meta.json').unlink()

    with open(out_dir / 'run_meta.json', 'w') as f:
        _ = json.dump(meta, f)
    meta = {'run_metadata': meta}
    (out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_run_metadata.json').unlink()

    with open(out_dir / 'g4x_viewer' / f'{g4x_obj.sample_id}_run_metadata.json', 'w') as f:
        _ = json.dump(meta, f)


def batched_demuxing(g4x_obj: 'G4Xoutput', manifest, probe_dict, batch_dir, batch_size):
    seq_reads = manifest['read'].unique().to_list()
    seq_reads = [int(x.split('_')[-1]) if isinstance(x, str) else x for x in seq_reads]

    num_features = pl.scan_parquet(g4x_obj.feature_table_path).select(pl.len()).collect().item()
    num_expected_batches = math.ceil(num_features / batch_size)
    cols_to_select = [
        'x_coord_shift',
        'y_coord_shift',
        'z',
        'demuxed',
        'transcript_condensed',
        'meanQS',
        'sequence_to_demux',
        'transcript',
        'TXUID',
    ]
    for i, feature_batch in tqdm(
        enumerate(g4x_obj.stream_features(batch_size, cols_to_select)),
        total=num_expected_batches,
        desc='Demuxing transcripts',
        position=0,
    ):
        feature_batch = feature_batch.with_columns(pl.col('TXUID').str.split('_').list.last().cast(int).alias('read'))
        redemuxed_feature_batch = []
        for seq_read in seq_reads:
            feature_batch_read = feature_batch.filter(pl.col('read') == seq_read)
            manifest_read = manifest.filter(pl.col('read') == seq_read)

            if len(feature_batch_read) == 0 or len(manifest_read) == 0:
                continue

            seqs = feature_batch_read['sequence_to_demux'].to_list()
            codes = manifest_read['sequence'].to_list()
            codebook_target_ids = np.array(manifest_read['probe_name'].to_list())

            hammings = batched_dot_product_hamming_matrix(seqs, codes, batch_size=batch_size)
            feature_batch_read = demux(hammings, feature_batch_read, codebook_target_ids, probe_dict)
            redemuxed_feature_batch.append(feature_batch_read)

        pl.concat(redemuxed_feature_batch).write_parquet(batch_dir / f'batch_{i}.parquet')


def concatenate_and_cleanup(batch_dir, out_dir):
    final_tx_table_path = out_dir / 'rna' / 'transcript_table.csv'
    final_tx_pq_path = out_dir / 'diagnostics' / 'transcript_table.parquet'

    utils.delete_existing(final_tx_table_path)
    utils.delete_existing(final_tx_table_path.with_suffix('.csv.gz'))
    utils.delete_existing(final_tx_pq_path)

    tx_table = pl.scan_parquet(list(batch_dir.glob('*.parquet')))
    tx_table.sink_parquet(final_tx_pq_path)

    _ = (
        tx_table.filter(pl.col('demuxed_new'))
        .drop(
            'transcript',
            'transcript_new',
            'transcript_condensed',
            'demuxed',
            'demuxed_new',
            'sequence_to_demux',
            'TXUID',
        )
        .rename(
            {
                'transcript_condensed_new': 'gene_name',
                'x_coord_shift': 'y_pixel_coordinate',
                'y_coord_shift': 'x_pixel_coordinate',
                'z': 'z_level',
                'meanQS': 'confidence_score',
            }
        )
        .sink_csv(final_tx_table_path)
    )

    _ = utils.gzip_file(final_tx_table_path)
    shutil.rmtree(batch_dir)
