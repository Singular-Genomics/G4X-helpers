# Notes
# all functions accept a G4Xoutput object (exception is migrate)
# if no other inputs are provided, all required inputs are loaded from the sample context
# output is written to specific locations relative to the sample directory, or optionally to an alternative directory
# the output paths are registered in the G4Xoutput object for downstream access
# (ie. if demux writes a TxTable to out_dir X, then subsequent functions will use TxTable in out_dir X as default input unless specified otherwise)

import logging
import shutil
from typing import TYPE_CHECKING, Literal

import numpy as np
import polars as pl

from . import constants as c
from . import io
from . import utils as ut
from .g4x_output import G4Xoutput
from .schema import definition as sd

if TYPE_CHECKING:
    from .schema.validator import BaseValidator

log = logging.getLogger(__name__)


# region demux
def demux(
    smp: 'G4Xoutput',
    manifest: str | None = None,
    out_dir: str | None = None,
    *,
    overwrite: bool = True,
    **kwargs,
) -> None:
    from .modules.demux import demux_raw_features

    out_dir = _ingest_out_dir(smp, out_dir)

    manifest_path = smp.src.Manifest.p if manifest is None else manifest

    raw_features_file = _collect_input(smp.src.RawFeatures.p, sd.RawFeatures)
    manifest_file = _collect_input(manifest_path, sd.Manifest)

    tx_table = demux_raw_features(
        raw_features=raw_features_file.load(lazy=True), manifest=manifest_file.load(), **kwargs
    )

    smp.reroute_source(sd.Manifest, out_dir, overwrite=overwrite)
    smp.reroute_source(sd.TxTable, out_dir, overwrite=overwrite)
    tx_table.write_csv(smp.src.TxTable.p, compression='gzip')
    manifest_file.load().write_csv(smp.src.Manifest.p)


# region aggregate
def aggregate(
    smp: 'G4Xoutput',
    cell_mask: str | None = None,
    out_dir: str | None = None,
    *,
    mask_key: str | None = None,
    overwrite: bool = True,
    show_progress: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
) -> None:
    from g4x_helpers.modules import aggregate

    compute_backend = io.get_backend(which=backend)

    out_dir = _ingest_out_dir(smp, out_dir)

    segmentation_path = smp.src.Segmentation.p if cell_mask is None else cell_mask
    segmentation_file = _collect_input(segmentation_path, sd.Segmentation)
    tx_table_file = _collect_input(smp.src.TxTable.p, sd.TxTable)

    save_segmentation = True
    if mask_key is None and segmentation_file.is_default:
        log.info('Loading both nuclei_exp and nuclei masks from default segmentation file')
        segmentation_mask = segmentation_file.load(key='nuclei_exp')
        nuclei_mask = segmentation_file.load(key='nuclei')
        save_segmentation = False
    else:
        mask_key = segmentation_file.available_keys[0] if mask_key is None else mask_key
        log.info('Loading specified segmentation mask with key: %s', mask_key)
        segmentation_mask = segmentation_file.load(key=mask_key)
        nuclei_mask = None

    log.info('Aggregating cell metadata')

    static_columns = {k: v for k, v in smp.smp_meta.items() if k in ['sample_id', 'tissue_type', 'block']}
    static_columns['seg_source'] = 'g4x-default' if segmentation_file.is_default else 'custom'

    cell_meta = aggregate.cell_metadata(
        cell_mask=segmentation_mask,
        nuclei_labels=nuclei_mask,
        static_columns=static_columns,
        show_progress=show_progress,
    )

    # cell_ids for the rest of the aggregation are extracted from the cell metadata to ensure consistency
    cell_ids = cell_meta.select('cell_id').collect()['cell_id'].to_list()

    log.info('Aggregating cell by gene matrix')
    tx_table = tx_table_file.load(lazy=True)
    cell_by_gene, tx_table = aggregate.cell_x_gene(
        tx_table=tx_table,
        cell_mask=segmentation_mask,
        included_cells=cell_ids,
        included_genes=smp.genes,
        return_tx_table=True,
    )

    images = {
        f'{stain}stain': smp.src.HnEDir.mapped_files[stain]
        for stain in smp.stains
        if stain in smp.src.HnEDir.mapped_files
    }

    if smp.src.pr_detected:
        images.update(smp.src.ProteinDir.mapped_files)

    log.info('Aggregating cell by img-signal matrix')
    bead_mask = smp.load_bead_mask()
    cell_by_signal = aggregate.cell_x_signal(
        images=images,
        cell_mask=segmentation_mask,
        bead_mask=bead_mask,
        included_cells=cell_ids,
        backend=compute_backend,
        show_progress=show_progress,
    )

    handles = [f'{k}_intensity_mean' for k in images.keys() if k.endswith('stain')]

    stain_data = cell_by_signal.select([c.CELL_ID_NAME] + handles)
    cell_meta = cell_meta.join(stain_data, on=c.CELL_ID_NAME, how='left').sort(c.CELL_ID_NAME)

    smp.reroute_source(sd.TxTable, out_dir, overwrite=True)
    tx_table.collect().write_csv(smp.src.TxTable.p, compression='gzip')

    smp.reroute_source(sd.CellMetadata, out_dir, overwrite=overwrite)
    cell_meta.sink_csv(smp.src.CellMetadata.p, compression='gzip')

    smp.reroute_source(sd.CellxGene, out_dir, overwrite=overwrite)
    cell_by_gene.sink_csv(smp.src.CellxGene.p, compression='gzip')

    if smp.src.pr_detected:
        log.info('Creating cell by protein matrix')
        cell_by_signal = cell_by_signal.drop(handles)
        smp.reroute_source(sd.CellxProt, out_dir, overwrite=overwrite)
        cell_by_signal.sink_csv(smp.src.CellxProt.p, compression='gzip')

    if save_segmentation:
        smp.reroute_source(sd.Segmentation, out_dir, overwrite=overwrite)
        mask_key = 'custom' if mask_key is None else mask_key
        mask_data = {mask_key: segmentation_mask}
        smp.src.Segmentation.main_key = mask_key
        np.savez(smp.src.Segmentation.p, **mask_data)


# region single cell processing
def sc_process(
    smp: 'G4Xoutput',
    out_dir: str | None = None,
    *,
    overwrite: bool = True,
    omit_correlation: bool = False,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
) -> None:
    from g4x_helpers.modules.single_cell import sc_utils
    from g4x_helpers.modules.single_cell.correlation import run_correlation_analysis
    from g4x_helpers.modules.single_cell.process import (
        default_sc_pipeline,
        dummy_clustering_output,
        dummy_dgex_output,
        init_adata,
        run_dgex,
    )

    compute_backend = io.get_backend(which=backend)

    out_dir = _ingest_out_dir(smp, out_dir)

    adata = init_adata(
        manifest=smp.src.Manifest.parse(),
        cell_metadata=smp.src.CellMetadata.load(),
        cell_x_gene=smp.src.CellxGene.load(),
        cell_x_protein=smp.src.CellxProt.load() if smp.src.pr_detected else None,
    )

    # 1: Befor we modify the adata object, we write cell metadata with QC metrics
    smp.reroute_source(sd.CellMetadata, out_dir, overwrite=True)
    sc_meta = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt32})
    sc_meta.write_csv(smp.src.CellMetadata.p, compression='gzip')

    adata, status = default_sc_pipeline(adata, compute_backend=compute_backend, **kwargs)

    success_clusterings = [k for k in adata.uns.keys() if k.startswith('leiden_')]

    # 4. Generate UMAP and clustering dataframes
    if status != 'success':
        clustering_umap = dummy_clustering_output(adata, failure_code=status)
    else:
        clustering_umap = sc_utils._extract_umap_clustering(adata, cluster_keys=success_clusterings)

    # TODO this try/except is technically duplicated within run_dgex... should consolidate this logic
    try:
        dgex = run_dgex(adata, cluster_keys=success_clusterings, downsample=1000)
        return dgex
    except Exception as e:
        log.warning(f'Failed to run differential gene expression analysis: {e}')
        dgex = dummy_dgex_output('dgex_failed')

    smp.reroute_source(sd.Dgex, out_dir, overwrite=overwrite)
    dgex.write_csv(smp.src.Dgex.p, compression='gzip')

    smp.reroute_source(sd.ClusteringUmap, out_dir, overwrite=overwrite)
    clustering_umap.write_csv(smp.src.ClusteringUmap.p, compression='gzip')

    smp.reroute_source(sd.AdataH5, out_dir, overwrite=overwrite)
    adata.write(smp.src.AdataH5.p)

    ### NOTE unclear if this is ok to run at the very end of the pipeline
    if smp.src.pr_detected and not omit_correlation:
        pr_corr_df, rna_pr_corr_df = run_correlation_analysis(adata)
        pr_corr_df.to_csv(smp.src.AdataH5.p.parent / 'protein_sc_correlation.csv')
        rna_pr_corr_df.to_csv(smp.src.AdataH5.p.parent / 'rna_protein_sc_correlation.csv')


# region migrate
def migrate(
    legacy_dir: str,
    out_dir: str,
    *,
    roi_coords: tuple[float, float, float, float] | None = None,
    downstream: bool = True,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    **kwargs,
) -> None:
    from .modules.migrate import migrate_legacy_raw_data

    legacy_dir = io.pathval.validate_dir_path(legacy_dir)
    out_dir = io.pathval.validate_dir_path(out_dir)

    migrate_legacy_raw_data(legacy_dir, out_dir=out_dir, roi_coords=roi_coords, **kwargs)
    smp = G4Xoutput(smp_dir=out_dir)

    if downstream:
        aggregate(smp, overwrite=True, backend=backend)
        sc_process(smp, overwrite=False, backend=backend)
        viewer_zarr(smp, overwrite=False)

    return smp


# region viewer
def viewer_zarr(
    smp: 'G4Xoutput',
    out_dir: str | None = None,
    *,
    overwrite: bool = True,
    symlink_images: bool = True,
) -> None:
    from g4x_helpers.modules.viewer.zarr_utils import link_viewer_group

    if out_dir is None:
        symlink_images = False

    out_dir = _ingest_out_dir(smp, out_dir)

    viewer_zarr_init(smp, out_dir=out_dir, overwrite=overwrite)

    source_viewer = smp.smp_dir / c.FILE_VIEWER_ZARR
    if source_viewer.exists() and symlink_images:
        link_viewer_group(smp, group_name='images', overwrite=True)
    else:
        viewer_zarr_images(smp, overwrite=overwrite)

    viewer_zarr_transcripts(smp, overwrite=overwrite)
    viewer_zarr_cells(smp, seg_name='g4x-default', overwrite=overwrite)


def viewer_zarr_init(
    smp: 'G4Xoutput',
    out_dir: str | None = None,
    *,
    overwrite: bool = True,
) -> None:
    from g4x_helpers.modules.viewer.zarr_utils import setup_viewer_zarr

    out_dir = _ingest_out_dir(smp, out_dir)

    smp.reroute_source(sd.ViewerZarr, out_dir, overwrite=overwrite)

    root_group = setup_viewer_zarr(smp.src.ViewerZarr.p, overwrite=overwrite)

    # write_metadata_defaults
    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

    if smp.src.QCSummary.path_exists():
        shutil.copy(smp.src.QCSummary.p, smp.src.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')


def viewer_zarr_images(
    smp: 'G4Xoutput',
    *,
    protein_list: list[str] | None = None,
    overwrite: bool = True,
    chunk_size: int = 1024,
) -> None:
    from .modules.viewer import images as viewer_img

    log.debug('Preparing multiplex image')
    images = {
        f'{stain}stain': smp.src.HnEDir.mapped_files[stain]
        for stain in smp.stains
        if stain in smp.src.HnEDir.mapped_files
    }

    if smp.src.pr_detected:
        pr_images = smp.src.ProteinDir.mapped_files
        if protein_list is not None:
            pr_images = {protein: pr_images[protein] for protein in protein_list if protein in pr_images}

        images.update(pr_images)

    if smp.src.pr_detected:
        visible_channels = viewer_img._determine_visible_channels(list(images.keys()))
    else:
        visible_channels = [c.NUCLEAR_STAIN + 'stain']

    channel_colors = {
        name: viewer_img.saturated_colors[viewer_img.channel_color_map[name]]
        for name in images.keys()
        if name in viewer_img.channel_color_map
    }

    viewer_img.write_images(
        smp.src.ViewerZarr.p,
        images=images,
        overwrite=overwrite,
        visible_channels=visible_channels,
        channel_colors=channel_colors,
        chunk_size=chunk_size,
    )

    log.debug('Preparing h_and_e image')
    he_file = smp.src.HnEDir.mapped_files['h_and_e']
    viewer_img.write_rgb_image(
        smp.src.ViewerZarr.p, image_name='h_and_e', image_path=he_file, overwrite=overwrite, chunk_size=chunk_size
    )


def viewer_zarr_transcripts(
    smp: 'G4Xoutput',
    *,
    overwrite: bool = True,
) -> None:
    from .modules.viewer import transcripts as viewer_tx

    viewer_tx.write_transcripts(
        smp.src.ViewerZarr.p,
        tx_table=smp.load_transcript_table(),
        manifest=smp.src.Manifest.parse(),
        data_shape=smp.shape,
        dgex=smp.src.Dgex.load(),
        overwrite=overwrite,
    )


def viewer_zarr_cells(
    smp: 'G4Xoutput',
    *,
    seg_name='g4x-default',
    overwrite: bool = True,
) -> None:
    from .modules.viewer import cells as viewer_cells

    viewer_cells.write_cells(
        smp.src.ViewerZarr.p,
        segmentation_mask=smp.src.Segmentation.load(),
        cell_metadata=smp.src.CellMetadata.load(),
        cell_x_gene=smp.src.CellxGene.load(),
        clustering_umap=smp.src.ClusteringUmap.load(),
        cell_x_protein=smp.src.CellxProt.load() if smp.src.pr_detected else None,
        seg_name=seg_name,
        overwrite=overwrite,
    )


# region private functions
def _ingest_out_dir(smp: 'G4Xoutput', out_dir: str | None) -> str:
    return smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)


def _collect_input(
    path: str,
    validator: 'BaseValidator',
    validate: bool = True,
) -> 'BaseValidator':
    path_valid = io.pathval.validate_file_path(path)
    in_obj = validator(target_path=path_valid)

    ut.log_with_path(f'Using the following file as input {validator.__name__}', in_obj.p, level='debug')

    if validate and not in_obj.is_valid:
        raise ValueError(f'Provided {validator.__name__} is not valid!\n{in_obj.report_validation()}')

    return in_obj
