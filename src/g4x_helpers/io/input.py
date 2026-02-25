from __future__ import annotations

import json
from functools import lru_cache, wraps
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas
import numpy as np
import polars as pl

from .. import utils
from . import convert

if TYPE_CHECKING:
    pass

PROBE_PATTERN = r'^(.*?)-([ACGT]{2,30})-([^-]+)$'

primer_read_map = {
    'SP1': 1,
    'm7a': 2,
    'm9a': 3,
    'm6a': 4,
    'm8a': 5,
    'm3a': 6,
}


def optionally_cached(func):
    cached_func = lru_cache(maxsize=32)(func)

    @wraps(func)
    def wrapper(*args, use_cache=True, **kwargs):
        if use_cache:
            return cached_func(*args, **kwargs)
        return func(*args, **kwargs)

    wrapper.cache_clear = cached_func.cache_clear
    return wrapper


def build_sample_metadata(sample_sheet, meta_json, sample_id):
    with open(meta_json, 'r') as f:
        run_meta = json.load(f)

    if 'sample_id' not in run_meta:
        run_meta['sample_id'] = sample_id

    sample_info, run_info = _parse_samplesheet(sample_sheet)

    ssheet_info = sample_info.filter(
        pl.col('lane') == int(sample_id[-1]), pl.col('sample_position') == sample_id[0]
    ).to_dicts()[0]

    ri_keys = {'run_name', 'user_name', 'assay'}
    filtered = [d for d in run_info.to_dicts() if d['Key'] in ri_keys]
    run_info = {d['Key']: d['Value'] for d in filtered}

    run_meta.update(run_info)
    run_meta.update(ssheet_info)

    order = [
        'run_name',
        'machine',
        'run_id',
        'platform',
        'assay',
        'user_name',
        'time_of_creation',
        'fc',
        'lane',
        'sample_position',
        'sample_id',
        'tissue_type',
        'block',
        'transcript_panel',
        'transcript_addon',
        'protein_panel',
        'protein_addon',
        'software',
        'software_version',
    ]

    return {k: run_meta[k] for k in order if k in run_meta}


def _parse_samplesheet(path):
    df = pl.read_csv(path, has_header=False, row_index_name='index')

    data_index = df.filter(pl.col('column_1') == '[Data]')['index'][0]

    run_info = (
        df.filter(pl.col('index').is_between(0, data_index, closed='none'))
        .select('column_1', 'column_2')
        .rename({'column_1': 'Key', 'column_2': 'Value'})
    )

    run_info = run_info.with_columns(pl.col('Key').str.to_lowercase().str.replace(' ', '_'))

    sample_info = pl.read_csv(path, skip_rows=data_index + 1)
    sample_info = sample_info.rename({c: c.lower().replace(' ', '_') for c in sample_info.columns})

    # TODO protein_addon is not always the last column ...
    final_col = np.where(np.array(sample_info.columns) == 'protein_addon')[0][0]
    sample_info = sample_info.select(pl.col(sample_info.columns[: final_col + 1])).fill_null('none')

    return sample_info, run_info


# @optionally_cached
def import_segmentation(seg_path: str, expected_shape: tuple[int], labels_key: str | None = None) -> np.ndarray:
    SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}
    ## load new segmentation
    cell_labels = utils.validate_path(seg_path, must_exist=True, is_dir_ok=False, is_file_ok=True)
    suffix = cell_labels.suffix.lower()

    if suffix not in SUPPORTED_MASK_FILETYPES:
        raise ValueError(f'{suffix} is not a supported file type.')

    if suffix == '.npz':
        with np.load(cell_labels) as labels:
            available_keys = list(labels.keys())

            if labels_key:  # if a key is specified
                if labels_key not in labels:
                    raise KeyError(f"Key '{labels_key}' not found in .npz; available keys: {available_keys}")
                seg = labels[labels_key]

            else:
                if len(available_keys) != 1:
                    raise ValueError(
                        f"Found multiple keys in .npz: {available_keys}.\nPlease specify a key using 'labels_key'"
                    )
                seg = labels[available_keys[0]]

    elif suffix == '.npy':
        # .npy: directly returns the array, no context manager available
        if labels_key is not None:
            print('file is .npy, ignoring provided labels_key.')
        seg = np.load(cell_labels, allow_pickle=False)

    elif suffix == '.geojson':
        gdf = geopandas.read_file(cell_labels)

        if labels_key is not None:
            if labels_key not in gdf.columns:
                raise KeyError(f"Column '{labels_key}' not found in GeoJSON; available columns: {gdf.columns.tolist()}")

            # ensure that a column named 'label' exists
            gdf['label'] = gdf[labels_key]

        else:
            if 'label' not in gdf.columns:
                raise ValueError(
                    "No column named 'label' found in GeoJSON. Please specify which column to use for labels via labels_key."
                )

        print('Rasterizing provided GeoDataFrame.')
        seg = convert.rasterize_polygons(gdf=gdf, target_shape=expected_shape)

    # validate shape for final numpy arrays
    if seg.shape != expected_shape:
        raise ValueError(f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}')

    return seg


def parse_input_manifest(file_path: Path, verbose: bool = False) -> pl.DataFrame:
    """
    Parse strings in `df[col]` of the form '<prefix>-<15mer>-<primer>'
    and add columns: gene_name, sequence, primer, primer_code, read.
    Rows that don't match the pattern get nulls in the new columns.
    """
    # print(f'Parsing transcript manifest: {file_path.name}')
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
        if verbose:
            print(f'{null_count} probes with invalid sequence format will be ignored:')
            null_seqs = parsed.filter(pl.col('sequence').is_null())['probe_name'].to_list()
            for ns in null_seqs:
                print(f'- {ns}')

    if 'panel_type' in df.columns:
        parsed = parsed.drop('panel_type')

    if 'read' in df.columns:
        if verbose:
            print("Using 'read' column provided by input manifest.")
        parsed = parsed.drop('read')
    else:
        plist = parsed['primer'].unique().to_list()
        ign_primer = [p for p in plist if p not in primer_read_map]
        if len(ign_primer) > 0:
            if verbose:
                print('Warning: the following primer names are not known and will be ignored:')
                for ip in ign_primer:
                    print(f'- {ip}')
        parsed = parsed.filter(pl.col('primer').is_in(primer_read_map.keys())).with_columns(
            pl.col('primer').replace(primer_read_map).cast(pl.Int8).alias('read')
        )

    if 'gene_name' in df.columns:
        if verbose:
            print("Using 'gene_name' column provided by input manifest.")
        parsed = parsed.drop('gene_name')

    manifest = parsed.join(df, on=col, how='left')

    manifest = manifest.with_columns(
        probe_type=(
            pl.when(pl.col('gene_name').str.to_lowercase() == 'gdna')
            .then(pl.lit('GCP'))
            .when(pl.col('gene_name').str.starts_with('NCS-'))
            .then(pl.lit('NCS'))
            .when(pl.col('gene_name').str.starts_with('NCP-'))
            .then(pl.lit('NCP'))
            .otherwise(pl.lit('targeting'))
        )
    )

    return manifest
