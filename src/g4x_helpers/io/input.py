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
import zarr
from matplotlib.pyplot import imread

from .. import constants as c
from . import convert, pathval

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from pandas import DataFrame as pdDF
    from polars import DataFrame as plDF
    from polars import LazyFrame as plLF


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


def parse_samplesheet(ss_path: str):
    ss_df = pl.read_csv(ss_path, has_header=False, row_index_name='index')

    err_pfx = 'Could not parse samplesheet. '
    if len(ss_df.columns) < 3:
        raise ValueError(f'{err_pfx} Too few columns detected')

    if '[Data]' not in ss_df['column_1']:
        raise ValueError(f"{err_pfx} Missing '[Data]' section-header.")

    data_index = ss_df.filter(pl.col('column_1') == '[Data]')['index'][0]

    run_info_section = (
        ss_df.filter(pl.col('index').is_between(0, data_index, closed='none'))
        .select('column_1', 'column_2')
        .rename({'column_1': 'Key', 'column_2': 'Value'})
        .cast({'Value': pl.Utf8})
    )

    data_section = pl.read_csv(ss_path, skip_rows=data_index + 1)
    data_section = data_section.drop(*([''] if '' in data_section.columns else []), cs.contains('_duplicated'))
    data_section = data_section.rename({c: c.replace('-', ' ') for c in data_section.columns})

    return run_info_section, data_section


def parse_input_manifest(manifest):

    for i, parsed_column in enumerate(['gene_name', 'sequence', 'primer']):
        if parsed_column not in manifest.columns:
            manifest = manifest.with_columns(
                [manifest['probe'].str.extract(c.PROBE_PATTERN, i + 1).alias(parsed_column)]
            )

    null_count = manifest.null_count()['sequence'][0]
    if null_count > 0:
        log.warning(f'{null_count} probes with invalid sequence format will be ignored:')
        null_seqs = manifest.filter(pl.col('sequence').is_null())['probe'].to_list()
        for ns in null_seqs:
            log.warning(f'- {ns}')

    if 'read_num' not in manifest.columns:
        plist = manifest['primer'].unique().to_list()
        ign_primer = [p for p in plist if p not in c.primer_read_map]

        if len(ign_primer) > 0:
            log.warning('Warning: the following primer names are not known and will be ignored:')
            for ip in ign_primer:
                log.warning(f'- {ip}')

        manifest = manifest.filter(pl.col('primer').is_in(c.primer_read_map.keys())).with_columns(
            pl.col('primer').replace(c.primer_read_map).cast(pl.Int8).alias('read_num')
        )

    if 'probe_id' not in manifest.columns:
        manifest = (
            manifest.with_columns(gene_probe_idx=pl.int_range(0, pl.len()).over('gene_name'))
            .with_columns(pl.col('gene_probe_idx').cast(pl.Utf8).str.zfill(4) + '-' + pl.col('gene_name'))
            .rename({'gene_probe_idx': 'probe_id'})
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

    manifest = manifest.with_columns(pl.col('probe_id').fill_null('<not provided>'))
    return manifest.select(['probe', 'probe_id', 'gene_name', 'sequence', 'primer', 'read_num', 'probe_type'])


@optionally_cached(maxsize=2)
def import_segmentation(
    seg_path: str, labels_key: str | None = None, *, expected_shape: tuple[int] | None
) -> np.ndarray:

    SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}

    ## load new segmentation
    cell_labels = pathval.validate_file_path(seg_path)

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
def import_image(
    img_path: str, extent: tuple[int, int, int, int] | None = None, n_threads: int | None = None
) -> np.ndarray:
    img_path = Path(img_path)
    suffix = img_path.suffix.lower()

    if extent is not None:
        extent = np.array(extent, dtype=np.int32)

    def _read_jp2(path, extent=extent, n_threads=n_threads):
        n_threads = c.DEFAULT_THREADS if n_threads is None else n_threads
        glymur.set_option('lib.num_threads', n_threads)

        img_path = str(path)
        if extent is not None:
            x0, x1, y0, y1 = extent
            return glymur.Jp2k(img_path)[y0:y1, x0:x1]
        return glymur.Jp2k(img_path)[:]

    def _read_tiff(path, extent=extent):
        if extent is not None:
            x0, x1, y0, y1 = extent
            store = tifffile.imread(path, aszarr=True)
            z = zarr.open(store, mode='r')
            return z[y0:y1, x0:x1]

        return tifffile.imread(path)

    def _read_standard_image(path):
        return imread(str(path))

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
