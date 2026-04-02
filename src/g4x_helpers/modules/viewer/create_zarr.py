from __future__ import annotations

import logging
import shutil
from pathlib import Path

import zarr

from ... import c, io
from ...schema.definition import ViewerZarr
from ..workflow import DEFAULT_INPUT, prepare_output
from .cells import process_cell_data, write_cells
from .images import write_images
from .transcripts import write_transcripts
from .utils import write_metadata_defaults

LOGGER = logging.getLogger(__name__)
DEFAULT_ZARR_NAME = c.FILE_VIEWER_ZARR


def create_viewer_zarr(
    smp,
    segmentation_mask: str = DEFAULT_INPUT,
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
    prepare_output(smp, out_dir, validator=ViewerZarr, overwrite=overwrite, logger=log)

    root_group = setup_zarr_tree(smp.data_dir, store_name=store_name, overwrite=overwrite)

    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

    write_images(smp, root_group, logger=log)
    write_transcripts(
        smp, root_group, tx_table=tx_table, manifest=manifest, dgex=sc_dir / smp.out.Dgex.n, overwrite=True, logger=log
    )

    prechew = process_cell_data(
        smp,
        segmentation_mask=segmentation_mask,
        cell_metadata=sc_dir / smp.out.CellMetadata.n,
        cell_x_gene=sc_dir / smp.out.CellxGene.n,
        cell_x_protein=None if not smp.src.pr_detected else sc_dir / smp.out.CellxProt.n,
        clustering_umap=sc_dir / smp.out.ClusteringUmap.n,
        logger=log,
    )
    write_cells(smp, root_group, prechew=prechew, overwrite=True, logger=log)

    if smp.src.QCSummary.path_exists:
        shutil.copy(smp.src.QCSummary.p, smp.out.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')


def setup_zarr_tree(target_dir: str, store_name: str = c.FILE_VIEWER_ZARR, overwrite: bool = True):
    target_dir = Path(target_dir) / store_name

    if target_dir.exists() and overwrite:
        print(f'Removing existing Zarr store at {target_dir}')
        shutil.rmtree(target_dir)

    root_group = zarr.open_group(target_dir, mode='w', zarr_version=2)

    img_group = root_group.create_group('images', overwrite=overwrite)

    img_group.create_group('multiplex', overwrite=overwrite)
    img_group.create_group('h_and_e', overwrite=overwrite)

    root_group.create_group('transcripts', overwrite=overwrite)

    cell_group = root_group.create_group('cells', overwrite=overwrite)
    cell_group.create_group('metadata', overwrite=overwrite)
    cell_group.create_group('protein', overwrite=overwrite)
    cell_group.create_group('polygons', overwrite=overwrite)
    cell_group.create_group('genes', overwrite=overwrite)

    (target_dir / 'misc').mkdir(parents=True, exist_ok=True)

    write_metadata_defaults(root_group)
    return root_group
