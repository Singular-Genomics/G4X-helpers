from __future__ import annotations

import os
from pathlib import Path

import polars as pl


def default_workers(max_workers: int = 16, reserve: int = 1) -> int:
    cpu = os.cpu_count() or 1
    return min(max_workers, max(1, cpu - reserve))


def write_table_to_csv(table: pl.DataFrame | pl.LazyFrame, out_path, source_path: str = None):
    if not isinstance(table, (pl.DataFrame, pl.LazyFrame)):
        raise ValueError(f'Expected a Polars DataFrame or LazyFrame, got {type(table)}')
    out_path = Path(out_path)

    is_lazy = isinstance(table, pl.LazyFrame)
    is_gz = out_path.suffix == '.gz'
    compression = 'gzip' if is_gz else 'uncompressed'

    if source_path is not None:
        if is_lazy and Path(source_path) == out_path:
            table = table.collect()
            is_lazy = False

    if is_lazy:
        table.sink_csv(out_path, compression=compression)
    else:
        table.write_csv(out_path, compression=compression)


def get_image_shape(img_path):
    import glymur
    import tifffile

    if img_path.suffix == '.tiff':
        with tifffile.TiffFile(img_path) as tif:
            series = tif.series[0]
            return series.shape
    elif img_path.suffix == '.jp2':
        return glymur.Jp2k(img_path).shape
    else:
        raise ValueError(f'Unsupported image format: {img_path.suffix}')


def kv_line_gap(key, value, gap=2):
    value = '<undefined>' if not value else value
    line = f'{key:<{gap}}'
    line += ' - '
    line += f'{value}'

    return line


def pretty_dict_str(d):
    max_len = max([len(k) for k in d.keys()])
    msg = ''
    for k, v in d.items():
        msg += kv_line_gap(k, v, gap=max_len) + '\n'
    return msg
