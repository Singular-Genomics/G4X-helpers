from __future__ import annotations

import logging
from functools import lru_cache, wraps
from inspect import signature
from pathlib import Path
from typing import TYPE_CHECKING

import dask
import dask.array as da
import geopandas
import glymur
import numpy as np
import polars as pl
import polars.selectors as cs
import tifffile
from matplotlib.pyplot import imread

from .. import c
from . import convert, pathval

if TYPE_CHECKING:
    from pandas import DataFrame as pdDF
    from polars import DataFrame as plDF
    from polars import LazyFrame as plLF

LOGGER = logging.getLogger(__name__)
log_p_gap = 2 * ' ' + '> '


def optionally_cached(*, maxsize=32, ignore_kwargs=()):
    ignore_kwargs = frozenset(ignore_kwargs)

    def decorator(func):
        sig = signature(func)

        @lru_cache(maxsize=maxsize)
        def cached_call(bound_args_items):
            bound = sig.bind_partial()
            bound.arguments.update(dict(bound_args_items))
            return func(*bound.args, **bound.kwargs)

        @wraps(func)
        def wrapper(*args, use_cache=False, **kwargs):
            if not use_cache:
                return func(*args, **kwargs)

            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            cache_dict = {k: v for k, v in bound.arguments.items() if k not in ignore_kwargs}

            cache_key = tuple(cache_dict.items())
            return cached_call(cache_key)

        wrapper.cache_clear = cached_call.cache_clear
        wrapper.cache_info = cached_call.cache_info
        return wrapper

    return decorator


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


def parse_samplesheet(ss_path: str):
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


def parse_input_manifest(file_path: str, verbose: bool = False):
    manifest = pl.read_csv(file_path)

    # col = 'probe_name'
    # if col not in manifest.columns:
    #     raise ValueError(f"transcript manifest must contain column '{col}'")

    for i, parsed_column in enumerate(['gene_name', 'sequence', 'primer']):
        if parsed_column not in manifest.columns:
            manifest = manifest.with_columns(
                [manifest['probe'].str.extract(c.PROBE_PATTERN, i + 1).alias(parsed_column)]
            )

    null_count = manifest.null_count()['sequence'][0]
    if null_count > 0:
        if verbose:
            print(f'{null_count} probes with invalid sequence format will be ignored:')
            null_seqs = manifest.filter(pl.col('sequence').is_null())['probe'].to_list()
            for ns in null_seqs:
                print(f'- {ns}')

    if 'read' not in manifest.columns:
        plist = manifest['primer'].unique().to_list()
        ign_primer = [p for p in plist if p not in c.primer_read_map]
        if len(ign_primer) > 0:
            if verbose:
                print('Warning: the following primer names are not known and will be ignored:')
                for ip in ign_primer:
                    print(f'- {ip}')
        manifest = manifest.filter(pl.col('primer').is_in(c.primer_read_map.keys())).with_columns(
            pl.col('primer').replace(c.primer_read_map).cast(pl.Int8).alias('read')
        )

    if 'probe_type' not in manifest.columns:
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

    return manifest.select(['probe', 'gene_name', 'sequence', 'primer', 'read', 'probe_type'])


@optionally_cached(maxsize=2)
def import_segmentation(
    seg_path: str, labels_key: str | None = None, *, expected_shape: tuple[int] | None
) -> np.ndarray:

    SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}

    ## load new segmentation
    cell_labels = pathval.validate_file_path(seg_path, must_exist=True)

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
    if expected_shape is not None:
        if seg.shape != expected_shape:
            raise ValueError(f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}')

    return seg


@optionally_cached(maxsize=2)
def import_table(file_path: str, lazy: bool = False, columns: tuple[str] | None = None) -> 'pdDF | plDF | plLF':
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


@optionally_cached(maxsize=8)
def import_image(img_path: str) -> np.ndarray:
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


def import_image_dask(img_path: str, shape: tuple[int], dtype=np.uint16, use_cache: bool = False) -> da.Array:
    return da.from_delayed(
        dask.delayed(import_image)(img_path, use_cache=use_cache),
        shape=shape,
        dtype=dtype,
    )


# def parse_input_manifest(file_path: Path, verbose: bool = False) -> pl.DataFrame:
#     """
#     Parse strings in `df[col]` of the form '<prefix>-<15mer>-<primer>'
#     and add columns: gene_name, sequence, primer, primer_code, read.
#     Rows that don't match the pattern get nulls in the new columns.
#     """
#     # print(f'Parsing transcript manifest: {file_path.name}')
#     df = pl.read_csv(file_path)

#     col = 'probe_name'
#     if col not in df.columns:
#         raise ValueError(f"transcript manifest must contain column '{col}'")

#     parsed = df.with_columns(
#         [
#             pl.col(col).str.extract(c.PROBE_PATTERN, 1).alias('gene_name'),
#             pl.col(col).str.extract(c.PROBE_PATTERN, 2).alias('sequence'),
#             pl.col(col).str.extract(c.PROBE_PATTERN, 3).alias('primer'),
#         ]
#     )

#     null_count = parsed.null_count()['sequence'][0]
#     if null_count > 0:
#         if verbose:
#             print(f'{null_count} probes with invalid sequence format will be ignored:')
#             null_seqs = parsed.filter(pl.col('sequence').is_null())['probe_name'].to_list()
#             for ns in null_seqs:
#                 print(f'- {ns}')

#     if 'panel_type' in df.columns:
#         parsed = parsed.drop('panel_type')

#     if 'read' in df.columns:
#         if verbose:
#             print("Using 'read' column provided by input manifest.")
#         parsed = parsed.drop('read')
#     else:
#         plist = parsed['primer'].unique().to_list()
#         ign_primer = [p for p in plist if p not in c.primer_read_map]
#         if len(ign_primer) > 0:
#             if verbose:
#                 print('Warning: the following primer names are not known and will be ignored:')
#                 for ip in ign_primer:
#                     print(f'- {ip}')
#         parsed = parsed.filter(pl.col('primer').is_in(c.primer_read_map.keys())).with_columns(
#             pl.col('primer').replace(c.primer_read_map).cast(pl.Int8).alias('read')
#         )

#     if 'gene_name' in df.columns:
#         if verbose:
#             print("Using 'gene_name' column provided by input manifest.")
#         parsed = parsed.drop('gene_name')

#     manifest = parsed.join(df, on=col, how='left')

#     manifest = manifest.with_columns(
#         probe_type=(
#             pl.when(pl.col('gene_name').str.to_lowercase() == 'gdna')
#             .then(pl.lit('GCP'))
#             .when(pl.col('gene_name').str.starts_with('NCS-'))
#             .then(pl.lit('NCS'))
#             .when(pl.col('gene_name').str.starts_with('NCP-'))
#             .then(pl.lit('NCP'))
#             .otherwise(pl.lit('targeting'))
#         )
#     )

#     return manifest
