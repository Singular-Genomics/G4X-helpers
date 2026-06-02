from __future__ import annotations

import colorsys

import polars as pl
import zarr

from ... import c

DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


# this function satisfies both zarr 2 and zarr 3 APIs, trying different combinations of parameters until one works
def create_array(group, name, data, compressor=None, chunks=None):
    create = getattr(group, 'create_array', None) or group.create_dataset

    attempts = [
        {'chunks': chunks, 'compressor': compressor},
        {'chunk_shape': chunks, 'compressors': [compressor] if compressor is not None else None},
        {'chunks': chunks, 'compressors': [compressor] if compressor is not None else None},
        {'chunk_shape': chunks, 'compressor': compressor},
        {},
    ]

    last_error = None
    for kwargs in attempts:
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        try:
            return create(name, data=data, **kwargs)
        except TypeError as e:
            last_error = e

    raise last_error


def calculate_chunks(arr, target_mb=4):
    TARGET_BYTES = target_mb * 1024 * 1024  # 4 MiB
    bytes_per_row = arr.dtype.itemsize if arr.ndim == 1 else arr.dtype.itemsize * arr.shape[1]
    row_chunk = max(1, TARGET_BYTES // bytes_per_row)
    return (row_chunk, *arr.shape[1:])


# def get_gene_metadata(viewer_dir: str):
#     g = zarr.open(viewer_dir, mode='r')
#     cmap = dict(g['transcripts'].attrs)['gene_colors']
#     cmap = {k: rgb_to_hex(v) for k, v in cmap.items()}

#     df = pl.DataFrame(dict(g['transcripts'].attrs)['gene_order'], schema=['gene_id'])
#     df = df.with_columns(pl.col('gene_id').replace(cmap).alias('color'))
#     return df


def hex_to_rgb(hex_color, normalized=False):
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    if normalized:
        return tuple(v / 255 for v in rgb)
    return rgb


def rgb_to_hex(rgb, normalized=False):
    if normalized:
        rgb = tuple(int(v * 255) for v in rgb)

    return '#{:02x}{:02x}{:02x}'.format(*rgb)


def hsv_to_hex(h, s, v):
    # wrap hue into [0,1]
    h = h % 1.0
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))
