from __future__ import annotations

import logging
import shutil
from pathlib import Path

import zarr

from ... import c, io
from ... import logging_utils as logut
from ...schema.definition import ViewerZarr
from ..workflow import DEFAULT_INPUT, route_output
from .images import write_images
from .transcripts import write_transcripts
from .utils import populate_zarr_metadata

LOGGER = logging.getLogger(__name__)
DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


def create_viewer_zarr(
    smp,
    tx_table: str = DEFAULT_INPUT,
    manifest: str = DEFAULT_INPUT,
    single_cell_dir: str = DEFAULT_INPUT,
    *,
    out_dir: str = DEFAULT_INPUT,
    overwrite=True,
    store_name: str = DEFAULT_ZARR_NAME,
    protein_list: list[str] | None = None,
    logger: logging.Logger | None = None,
) -> None:
    log = logger or LOGGER
    log.info('Running create_viewer_zarr')

    if protein_list is not None:
        smp.set_proteins(protein_list)

    single_cell_dir = smp.out.SingleCellFolder.p if single_cell_dir == DEFAULT_INPUT else single_cell_dir
    sc_dir = io.pathval.validate_dir_path(single_cell_dir, must_exist=True)

    out_dir = smp.data_dir if out_dir == DEFAULT_INPUT else io.pathval.validate_dir_path(out_dir)
    route_output(smp, out_dir, validator=ViewerZarr, overwrite=overwrite, logger=log)

    root_group = setup_zarr_tree(smp.data_dir, store_name=store_name, overwrite=overwrite)

    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

    if smp.src.QCSummary.path_exists:
        shutil.copy(smp.src.QCSummary.p, smp.out.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')

    write_images(smp, root_group, logger=log)
    write_transcripts(
        smp, root_group, tx_table=tx_table, manifest=manifest, dgex=sc_dir / smp.out.Dgex.n, overwrite=True, logger=log
    )

    return root_group


def setup_zarr_tree(
    target_dir: str, store_name: str = c.FILE_VIEWER_ZARR, overwrite: bool = True, logger: logging.Logger | None = None
) -> zarr.Group:
    log = logger or LOGGER
    target_dir = Path(target_dir) / store_name

    if target_dir.exists() and overwrite:
        logut.log_with_path('Removing existing Zarr store at:', target_dir, logger=log, level='INFO')
        shutil.rmtree(target_dir)

    root_group = zarr.open_group(target_dir, mode='w', zarr_version=2)

    img_group = root_group.create_group('images', overwrite=overwrite)

    img_group.create_group('multiplex', overwrite=overwrite)
    img_group.create_group('h_and_e', overwrite=overwrite)

    root_group.create_group('transcripts', overwrite=overwrite)

    cell_group = root_group.create_group('cells', overwrite=overwrite)
    cell_group.attrs['segmentation_sources'] = {}
    cell_group.attrs['segmentation_order'] = []
    # for seg_name in ['default_segmentation', 'custom_segmentation']:
    #     seg_group = cell_group.create_group(seg_name, overwrite=overwrite)
    #     seg_group.create_group('metadata', overwrite=overwrite)
    #     seg_group.create_group('protein', overwrite=overwrite)
    #     seg_group.create_group('polygons', overwrite=overwrite)
    #     seg_group.create_group('genes', overwrite=overwrite)

    (target_dir / 'misc').mkdir(parents=True, exist_ok=True)

    # write_metadata_defaults
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
    return root_group
