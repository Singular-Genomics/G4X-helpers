import shutil
from pathlib import Path

import zarr

from .. import constants as c
from .cells import generate_cluster_palette


def setup_zarr_tree(target_dir: str, store_name: str = c.FILE_VIEWER_ZARR, replace_existing: bool = True):
    target_dir = Path(target_dir) / store_name

    if target_dir.exists() and replace_existing:
        print(f'Removing existing Zarr store at {target_dir}')
        shutil.rmtree(target_dir)

    root_group = zarr.open_group(target_dir, mode='w', zarr_version=2)

    img_group = root_group.create_group('images', overwrite=replace_existing)

    _ = img_group.create_group('multiplex', overwrite=replace_existing)
    _ = img_group.create_group('h_and_e', overwrite=replace_existing)

    _ = root_group.create_group('transcripts', overwrite=replace_existing)

    cell_group = root_group.create_group('cells', overwrite=replace_existing)
    _ = cell_group.create_group('metadata', overwrite=replace_existing)
    _ = cell_group.create_group('protein', overwrite=replace_existing)
    _ = cell_group.create_group('polygons', overwrite=replace_existing)
    _ = cell_group.create_group('genes', overwrite=replace_existing)

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
    setup_zarr_metadata(
        root_group,
        sample_metadata={'sample_metadata': 'default'},
        cluster_ids=['default'],
        gene_mtx_shape=(0, 0),
        gene_colors={},
        tx_layer_config=tx_layer_info,
    )


def setup_zarr_metadata(
    root_group: zarr.Group,
    sample_metadata: dict | None = None,
    cluster_ids: list | None = None,
    gene_mtx_shape: tuple | None = None,
    gene_colors: list | None = None,
    tx_layer_config: dict | None = None,
):

    if sample_metadata is not None:
        root_group.attrs['run_metadata'] = {'Sample Information': sample_metadata}

    if cluster_ids is not None:
        root_group['cells']['metadata'].attrs['clusterID_colors'] = generate_cluster_palette(cluster_ids)

    if gene_mtx_shape is not None:
        root_group['cells']['genes'].attrs['shape'] = gene_mtx_shape

    if gene_colors is not None:
        root_group['transcripts'].attrs['gene_colors'] = gene_colors

    if tx_layer_config is not None:
        root_group['transcripts'].attrs['layer_config'] = tx_layer_config

    # set static image attrs once
    # if 'axes' not in root_group['images'].attrs:
    #     root_group['images'].attrs['axes'] = {'unit': 'micrometer', 'pixel_per_um': c.PIXEL_PER_MICRON}


# meta_path = smp.data_dir / c.FILE_SMP_META
#     with open(meta_path, 'r') as f:
#         run_metadata = json.load(f)

# shutil.copy(smp.data_dir / c.FILE_SUMMARY, target_dir / 'misc' / c.FILE_SUMMARY)
