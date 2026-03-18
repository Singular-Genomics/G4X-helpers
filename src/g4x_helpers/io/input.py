from __future__ import annotations

from functools import lru_cache, wraps
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas
import glymur
import numpy as np
import polars as pl
import polars.selectors as cs
import tifffile
from matplotlib.pyplot import imread

from .. import constants, utils
from . import convert

if TYPE_CHECKING:
    from pandas import DataFrame as pdDF
    from polars import DataFrame as plDF
    from polars import LazyFrame as plLF


def optionally_cached(func):
    cached_func = lru_cache(maxsize=32)(func)

    @wraps(func)
    def wrapper(*args, use_cache=True, **kwargs):
        if use_cache:
            return cached_func(*args, **kwargs)
        return func(*args, **kwargs)

    wrapper.cache_clear = cached_func.cache_clear
    return wrapper


# TODO this belongs to migration logic
def _infer_smp_id(sample_dir):
    summary_file = list(sample_dir.glob('summary*.html'))[0].stem
    core_metrics = sample_dir / 'metrics' / 'transcript_core_metrics.csv'

    smp_id_met = pl.read_csv(core_metrics).to_dicts()[0]['sample_id']
    smp_id_sum = summary_file.split('_')[-1]

    if smp_id_met == smp_id_sum:
        smp_id = smp_id_met
    else:
        smp_id = None

    return smp_id


def _parse_samplesheet(ss_path: str):
    ss_df = pl.read_csv(ss_path, has_header=False, row_index_name='index')

    data_index = ss_df.filter(pl.col('column_1') == '[Data]')['index'][0]

    run_info_section = (
        ss_df.filter(pl.col('index').is_between(0, data_index, closed='none'))
        .select('column_1', 'column_2')
        .rename({'column_1': 'Key', 'column_2': 'Value'})
        .cast({'Value': pl.Utf8})
    )

    data_section = pl.read_csv(ss_path, skip_rows=data_index + 1)
    data_section = data_section.rename({k: k.replace('Addon', 'Custom') for k in data_section.columns if 'Addon' in k})
    data_section = data_section.drop(cs.contains('_duplicated'))

    return run_info_section, data_section


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
            pl.col(col).str.extract(constants.PROBE_PATTERN, 1).alias('gene_name'),
            pl.col(col).str.extract(constants.PROBE_PATTERN, 2).alias('sequence'),
            pl.col(col).str.extract(constants.PROBE_PATTERN, 3).alias('primer'),
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
        ign_primer = [p for p in plist if p not in constants.primer_read_map]
        if len(ign_primer) > 0:
            if verbose:
                print('Warning: the following primer names are not known and will be ignored:')
                for ip in ign_primer:
                    print(f'- {ip}')
        parsed = parsed.filter(pl.col('primer').is_in(constants.primer_read_map.keys())).with_columns(
            pl.col('primer').replace(constants.primer_read_map).cast(pl.Int8).alias('read')
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


@optionally_cached
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
        seg = convert.gdf_to_ndarray(gdf=gdf, target_shape=expected_shape)

    # validate shape for final numpy arrays
    if seg.shape != expected_shape:
        raise ValueError(f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}')

    return seg


@optionally_cached
def import_table(file_path: str, lazy: bool = False, columns: list[str] | None = None) -> 'pdDF | plDF | plLF':
    file_path = Path(file_path)
    if lazy:
        if file_path.suffix == '.parquet':
            reads = pl.scan_parquet(file_path)
        else:
            reads = pl.scan_csv(file_path)
    else:
        if file_path.suffix == '.parquet':
            reads = pl.read_parquet(file_path)
        else:
            reads = pl.read_csv(file_path)
    if columns:
        reads = reads.select(columns)

    return reads


@optionally_cached
def import_image(img_path):
    img_path = Path(img_path)
    suffix = img_path.suffix.lower()

    def _read_jp2(path):
        return glymur.Jp2k(str(path))[:]

    def _read_standard_image(path):
        return imread(str(path))

    def _read_tiff(path):
        return tifffile.imread(path)

    readers = {
        '.jp2': _read_jp2,
        '.png': _read_standard_image,
        '.jpg': _read_standard_image,
        '.jpeg': _read_standard_image,
        '.tif': _read_tiff,
        '.tiff': _read_tiff,
    }

    try:
        return readers[suffix](img_path)
    except KeyError:
        raise ValueError(f'Unsupported image format: {suffix}') from None


# def create_sample_g4x(sample_dir, save_file: bool = False):
#     sample_dir = Path(sample_dir)

#     meta_json = sample_dir / 'run_meta.json'
#     sample_sheet = sample_dir / 'samplesheet.csv'

#     with open(meta_json, 'r') as f:
#         run_meta = json.load(f)

#     if 'sample_id' in run_meta:
#         sample_id = run_meta['sample_id']
#     else:
#         sample_id = _infer_smp_id(sample_dir)
#         run_meta['sample_id'] = sample_id

#     run_info, data_section = _parse_samplesheet(sample_sheet)

#     ssheet_info = data_section.filter(
#         pl.col('lane') == int(sample_id[-1]), pl.col('sample_position') == sample_id[0]
#     ).to_dicts()[0]

#     ri_keys = {'run_name', 'user_name', 'assay'}
#     filtered = [d for d in run_info.to_dicts() if d['Key'] in ri_keys]
#     run_info = {d['Key']: d['Value'] for d in filtered}

#     run_meta.update(run_info)
#     run_meta.update(ssheet_info)

#     order = [
#         'run_name',
#         'sample_position',
#         'sample_id',
#         'tissue_type',
#         'block',
#         'assay',
#         'machine',
#         'run_id',
#         'fc',
#         'lane',
#         'platform',
#         'user_name',
#         'time_of_creation',
#         'transcript_panel',
#         'transcript_custom',
#         'protein_panel',
#         'protein_custom',
#         'software',
#         'software_version',
#         'output_version',
#     ]

#     out = {k: str(run_meta[k]) for k in order if k in run_meta}
#     if save_file:
#         with open(sample_dir / 'sample.g4x', 'w') as f:
#             json.dump(out, f, indent=4)

#     return out
