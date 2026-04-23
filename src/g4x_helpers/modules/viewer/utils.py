from __future__ import annotations

import zarr

from ... import c

DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


def populate_zarr_metadata(
    root_group: zarr.Group,
    sample_metadata: dict | None = None,
    cluster_ids: dict | None = None,
    gene_mtx_shape: tuple | None = None,
    gene_colors: list | None = None,
    tx_layer_config: dict | None = None,
):

    if sample_metadata is not None:
        root_group.attrs['run_metadata'] = {'Sample Information': sample_metadata}

    # if cluster_ids is not None:
    #     root_group['cells']['metadata'].attrs['clusterID_order'] = list(cluster_ids.keys())
    #     root_group['cells']['metadata'].attrs['clusterID_colors'] = cluster_ids

    # if gene_mtx_shape is not None:
    #     root_group['cells']['genes'].attrs['shape'] = gene_mtx_shape

    if gene_colors is not None:
        root_group['transcripts'].attrs['gene_order'] = list(gene_colors.keys())
        root_group['transcripts'].attrs['gene_colors'] = gene_colors

    if tx_layer_config is not None:
        root_group['transcripts'].attrs['layer_config'] = tx_layer_config


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
