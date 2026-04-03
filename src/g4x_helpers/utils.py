from __future__ import annotations

import os
from pathlib import Path

import polars as pl


def default_workers(max_workers: int = 16, reserve: int = 1) -> int:
    cpu = os.cpu_count() or 1
    return min(max_workers, max(1, cpu - reserve))


def validate_path(path_str, must_exist=True, is_dir_ok=True, is_file_ok=True, resolve_path=False) -> Path:
    path = Path(path_str)  # .expanduser().resolve()
    if resolve_path:
        path = path.resolve()

    if must_exist and not path.exists():
        raise FileNotFoundError(f'Path does not exist: {path}')

    if path.exists():
        if path.is_dir() and not is_dir_ok:
            raise ValueError(f'Expected a file but got a directory: {path}')
        if path.is_file() and not is_file_ok:
            raise ValueError(f'Expected a directory but got a file: {path}')

    if not path.exists() and not must_exist:
        if not is_file_ok:
            path.mkdir(parents=True, exist_ok=True)
        else:
            parent = path.parent
            parent.mkdir(parents=True, exist_ok=True)

    return path


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
