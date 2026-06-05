import logging
import shutil
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from . import constants as c
from . import io
from .schema import definition as sd

if TYPE_CHECKING:
    from .g4x_output import G4Xoutput
    from .schema.validator import BaseValidator

LOGGER = logging.getLogger(__name__)


def demux(
    smp: 'G4Xoutput',
    manifest: str | None = None,
    out_dir: str | None = None,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
    **kwargs,
) -> None:

    from .modules.demux import demux_raw_features

    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    manifest_path = smp.src.Manifest.p if manifest is None else manifest
    manifest_file = _collect_input(manifest_path, sd.Manifest)
    raw_features_file = smp.src.RawFeatures

    tx_table = demux_raw_features(
        raw_features=raw_features_file.load(lazy=True), manifest=manifest_file.parse(), logger=log, **kwargs
    )

    smp.reroute_source(sd.Manifest, out_dir, overwrite=overwrite)
    smp.reroute_source(sd.TxTable, out_dir, overwrite=overwrite)
    tx_table.write_csv(smp.src.TxTable.p, compression='gzip')
    manifest_file.load().write_csv(smp.src.Manifest.p)


def aggregate(
    smp: 'G4Xoutput',
    segmentation_mask: str | None = None,
    out_dir: str | None = None,
    *,
    mask_key: str | None = None,
    overwrite: bool = True,
    show_progress: bool = False,
    compute_backend: io.ComputeBackend = io.get_backend(which='auto'),
    logger: logging.Logger | None = None,
) -> None:

    from g4x_helpers.modules.aggregate import create_cell_metadata, create_cell_x_gene, create_cell_x_signal

    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    static_columns = {k: v for k, v in smp.smp_meta.items() if k in ['sample_id', 'tissue_type', 'block']}
    static_columns['seg_source'] = 'g4x-default' if segmentation_mask is None else 'custom'
    log.info('Using %s segmentation source', static_columns['seg_source'])

    segmentation_path = smp.src.Segmentation.p if segmentation_mask is None else segmentation_mask
    segmentation_file = _collect_input(segmentation_path, sd.Segmentation)

    if mask_key is None and segmentation_mask is None:
        cell_mask = segmentation_file.load(key='nuclei_exp')
        nuclei_mask = segmentation_file.load(key='nuclei')
    else:
        cell_mask = segmentation_file.load(key=segmentation_file.available_keys[0])
        nuclei_mask = None

    cell_meta = create_cell_metadata(
        segmentation_mask=cell_mask, nuclei_mask=nuclei_mask, static_columns=static_columns, show_progress=show_progress
    )

    # cell_ids for the rest of the aggregation are extracted from the cell metadata to ensure consistency
    cell_ids = cell_meta.select('cell_id').collect()['cell_id'].to_list()

    tx_table = smp.src.TxTable.load(lazy=True)
    cell_by_gene, tx_table = create_cell_x_gene(
        tx_table=tx_table,
        segmentation_mask=cell_mask,
        included_cells=cell_ids,
        return_tx_table=True,
    )

    images = {
        f'{stain}stain': smp.src.HnEDir.mapped_files[stain]
        for stain in smp.stains
        if stain in smp.src.HnEDir.mapped_files
    }

    if smp.src.pr_detected:
        images.update(smp.src.ProteinDir.mapped_files)

    bead_mask = smp.load_bead_mask()
    cell_by_signal = create_cell_x_signal(
        images=images,
        segmentation_mask=cell_mask,
        bead_mask=bead_mask,
        included_cells=cell_ids,
        backend=compute_backend,
        show_progress=show_progress,
    )

    handles = [f'{k}_intensity_mean' for k in images.keys() if k.endswith('stain')]

    stain_data = cell_by_signal.select([c.CELL_ID_NAME] + handles)
    cell_meta = cell_meta.join(stain_data, on=c.CELL_ID_NAME, how='left').sort(c.CELL_ID_NAME)

    smp.reroute_source(sd.TxTable, out_dir, overwrite=overwrite)
    tx_table.collect().write_csv(smp.src.TxTable.p, compression='gzip')

    smp.reroute_source(sd.CellMetadata, out_dir, overwrite=overwrite)
    cell_meta.sink_csv(smp.src.CellMetadata.p, compression='gzip')

    smp.reroute_source(sd.CellxGene, out_dir, overwrite=overwrite)
    cell_by_gene.sink_csv(smp.src.CellxGene.p, compression='gzip')

    if smp.src.pr_detected:
        cell_by_signal = cell_by_signal.drop(handles)
        smp.reroute_source(sd.CellxProt, out_dir, overwrite=overwrite)
        cell_by_signal.sink_csv(smp.src.CellxProt.p, compression='gzip')

    if segmentation_mask is not None:
        smp.reroute_source(sd.Segmentation, out_dir, overwrite=overwrite)
        mask_key = 'custom' if mask_key is None else mask_key
        mask_data = {mask_key: cell_mask}
        smp.src.Segmentation.main_key = mask_key
        np.savez(smp.src.Segmentation.p, **mask_data)


def sc_process(
    smp: 'G4Xoutput',
    out_dir: str | None = None,
    overwrite: bool = True,
    omit_correlation: bool = False,
    logger: logging.Logger | None = None,
    **kwargs,
):

    from g4x_helpers.modules.single_cell import sc_utils
    from g4x_helpers.modules.single_cell.correlation import run_correlation_analysis
    from g4x_helpers.modules.single_cell.process import (
        default_sc_pipeline,
        dummy_clustering_output,
        dummy_dgex_output,
        init_adata,
        run_dgex,
    )

    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    adata = init_adata(
        manifest=smp.src.Manifest.parse(),
        cell_metadata=smp.src.CellMetadata.load(),
        cell_x_gene=smp.src.CellxGene.load(),
        cell_x_protein=smp.src.CellxProt.load(),
    )

    # 1: Befor we modify the adata object, we write cell metadata with QC metrics
    smp.reroute_source(sd.CellMetadata, out_dir, overwrite=overwrite)
    sc_meta = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt32})
    sc_meta.write_csv(smp.src.CellMetadata.p, compression='gzip')

    adata, status = default_sc_pipeline(adata)

    success_clusterings = [k for k in adata.uns.keys() if k.startswith('leiden_')]

    # 4. Generate UMAP and clustering dataframes
    if status != 'success':
        clustering_umap = dummy_clustering_output(adata, failure_code=status)
    else:
        clustering_umap = sc_utils._extract_umap_clustering(adata, cluster_keys=success_clusterings)

    try:
        dgex = run_dgex(adata, cluster_keys=success_clusterings, downsample=1000)
    except Exception as e:
        log.warning(f'Failed to run differential gene expression analysis: {e}')
        dgex = dummy_dgex_output('dgex_failed')

    smp.reroute_source(sd.Dgex, out_dir, overwrite=overwrite)
    dgex.write_csv(smp.src.Dgex.p, compression='gzip')

    smp.reroute_source(sd.ClusteringUmap, out_dir, overwrite=overwrite)
    clustering_umap.write_csv(smp.src.ClusteringUmap.p, compression='gzip')

    smp.reroute_source(sd.AdataH5, out_dir, overwrite=overwrite)
    adata.write(smp.src.AdataH5.p)

    ### unclear if this is ok to run at the very end of the pipeline
    if smp.src.pr_detected and not omit_correlation:
        pr_corr_df, rna_pr_corr_df = run_correlation_analysis(adata)
        pr_corr_df.to_csv(smp.src.AdataH5.p.parent / 'protein_sc_correlation.csv')
        rna_pr_corr_df.to_csv(smp.src.AdataH5.p.parent / 'rna_protein_sc_correlation.csv')


def init_viewer_zarr(
    smp,
    out_dir: str | None = None,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
) -> None:
    from g4x_helpers.modules.viewer.zarr_utils import setup_viewer_zarr

    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    smp.reroute_source(sd.ViewerZarr, out_dir, overwrite=overwrite)

    root_group = setup_viewer_zarr(smp.src.ViewerZarr.p, overwrite=overwrite)

    # write_metadata_defaults
    root_group.attrs['run_metadata'] = {'Sample Information': smp.smp_meta}
    root_group.attrs['smp_info_order'] = list(smp.smp_meta.keys())

    if smp.src.QCSummary.path_exists():
        shutil.copy(smp.src.QCSummary.p, smp.src.ViewerZarr.p / 'misc' / 'summary.html')
    else:
        log.info('QCSummary file does not exist, skipping copy to ViewerZarr.')


def write_viewer_images(
    smp,
    protein_list: list[str] | None = None,
    overwrite: bool = True,
    chunk_size: int = 1024,
    logger: logging.Logger | None = None,
):
    from .modules.viewer import images as viewer_img

    log = logger or LOGGER
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

    viewer_img.write_images_to_zarr(
        smp.src.ViewerZarr.p,
        images=images,
        overwrite=overwrite,
        visible_channels=visible_channels,
        channel_colors=channel_colors,
        chunk_size=chunk_size,
        logger=log,
    )

    he_file = smp.src.HnEDir.mapped_files['h_and_e']
    viewer_img.write_rgb_img(
        smp.src.ViewerZarr.p,
        image_name='h_and_e',
        image_path=he_file,
        overwrite=overwrite,
        chunk_size=chunk_size,
        logger=log,
    )


def write_viewer_transcripts(
    smp,
    overwrite: bool = True,
    logger: logging.Logger | None = None,
):
    from g4x_helpers.modules.viewer.transcripts import write_transcripts_to_zarr

    log = logger or LOGGER
    log.debug('Preparing multiplex image')

    write_transcripts_to_zarr(
        smp.src.ViewerZarr.p,
        tx_table=smp.load_transcript_table(),
        manifest=smp.src.Manifest.parse(),
        data_shape=smp.shape,
        dgex=smp.src.Dgex.load(),
        overwrite=overwrite,
    )


def write_viewer_cells(
    smp,
    seg_name='g4x-default',
    overwrite: bool = True,
    logger: logging.Logger | None = None,
):
    from g4x_helpers.modules.viewer.cells import write_cells_to_zarr

    log = logger or LOGGER
    log.debug('Preparing cells')

    write_cells_to_zarr(
        smp.src.ViewerZarr.p,
        segmentation_mask=smp.src.Segmentation.load(),
        cell_metadata=smp.src.CellMetadata.load(),
        cell_x_gene=smp.src.CellxGene.load(),
        clustering_umap=smp.src.ClusteringUmap.load(),
        cell_x_protein=smp.src.CellxProt.load(),
        seg_name=seg_name,
        overwrite=overwrite,
    )


def _collect_input(
    path: str,
    validator: 'BaseValidator',
    validate: bool = True,
):
    path_valid = io.pathval.validate_file_path(path)
    in_obj = validator(target_path=path_valid)

    if validate and not in_obj.is_valid:
        raise ValueError(f'Provided {validator.__name__} is not valid!\n{in_obj.report_validation()}')

    return in_obj
