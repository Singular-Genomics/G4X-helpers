from __future__ import annotations

import inspect
import shutil
from pathlib import Path

import zarr

from .. import constants as c


def setup_zarr_tree(target_dir: str, store_name: str = c.FILE_VIEWER_ZARR, overwrite: bool = True):
    target_dir = Path(target_dir) / store_name

    if target_dir.exists() and overwrite:
        print(f'Removing existing Zarr store at {target_dir}')
        shutil.rmtree(target_dir)

    root_group = zarr.open_group(target_dir, mode='w', zarr_version=2)

    img_group = root_group.create_group('images', overwrite=overwrite)

    _ = img_group.create_group('multiplex', overwrite=overwrite)
    _ = img_group.create_group('h_and_e', overwrite=overwrite)

    _ = root_group.create_group('transcripts', overwrite=overwrite)

    cell_group = root_group.create_group('cells', overwrite=overwrite)
    _ = cell_group.create_group('metadata', overwrite=overwrite)
    _ = cell_group.create_group('protein', overwrite=overwrite)
    _ = cell_group.create_group('polygons', overwrite=overwrite)
    _ = cell_group.create_group('genes', overwrite=overwrite)

    (target_dir / 'misc').mkdir(parents=True, exist_ok=True)

    write_metadata_defaults(root_group)
    return root_group


def write_metadata_defaults(root_group: zarr.Group):

    # set static image attrs once
    root_group['images'].attrs['axes'] = {'unit': 'micrometer', 'pixel_per_um': c.PIXEL_PER_MICRON}

    tx_layer_info = {
        'layers': 1,
        'tile_size': 1,
        'layer_height': 1,
        'layer_width': 1,
        'coordinate_order': ['default_x', 'default_y'],
    }
    populate_zarr_metadata(
        root_group,
        sample_metadata={'sample_metadata': 'default'},
        cluster_ids={'default': [114, 255, 171]},
        gene_mtx_shape=(0, 0),
        gene_colors={},
        tx_layer_config=tx_layer_info,
    )


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

    if cluster_ids is not None:
        root_group['cells']['metadata'].attrs['clusterID_colors'] = cluster_ids

    if gene_mtx_shape is not None:
        root_group['cells']['genes'].attrs['shape'] = gene_mtx_shape

    if gene_colors is not None:
        root_group['transcripts'].attrs['gene_colors'] = gene_colors

    if tx_layer_config is not None:
        root_group['transcripts'].attrs['layer_config'] = tx_layer_config


def _supports_param(func, name):
    try:
        return name in inspect.signature(func).parameters
    except (TypeError, ValueError):
        return False


def _normalize_chunks(chunks):
    if chunks == 'auto':
        return None
    return chunks


def create_array(group, name, data, compressor=None, chunks=None):
    kwargs = {}
    chunks = _normalize_chunks(chunks)

    create = getattr(group, 'create_array', None)
    if create is None:
        create = group.create_dataset

    if chunks is not None:
        if _supports_param(create, 'chunks'):
            kwargs['chunks'] = chunks
        elif _supports_param(create, 'chunk_shape'):
            kwargs['chunk_shape'] = chunks

    if compressor is not None:
        if _supports_param(create, 'compressors'):
            kwargs['compressors'] = [compressor]
        elif _supports_param(create, 'compressor'):
            kwargs['compressor'] = compressor

    return create(name, data=data, **kwargs)
