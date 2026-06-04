import logging
from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from . import constants as c
from . import io
from .schema import definition as sd

if TYPE_CHECKING:
    from .g4x_output import G4Xoutput
    from .modules.single_cell.filtering import FilterPanel
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

    # cell_ids for the rest of the aggregation must be extracted from the cell metadata to ensure consistency
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
    cell_meta = cell_meta.join(stain_data, on=c.CELL_ID_NAME, how='left')

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


def single_cell(
    smp: 'G4Xoutput',
    out_dir: str | None = None,
    *,
    overwrite: bool = True,
    filter_panel: 'FilterPanel' | None = None,
    n_neighbors: int = 15,
    cluster_attempts: int = 10,
    rnd_st: int = 111,
    compute_backend: io.ComputeBackend = io.get_backend(which='auto'),
    omit_correlation: bool = False,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    out_dir = smp.smp_dir if out_dir is None else io.pathval.validate_dir_path(out_dir)

    backend = io.get_backend(which=compute_backend)

    from g4x_helpers.modules import single_cell as sc_module

    filter_panel = sc_module.filtering._get_default_filter_panel() if filter_panel is None else filter_panel

    manifest = smp.src.Manifest.parse()
    cell_x_gene = smp.src.CellxGene.load()
    cell_metadata = smp.src.CellMetadata.load()
    cell_x_protein = smp.src.CellxProt.load() if smp.src.pr_detected else None

    adata = sc_module.init_adata.init_adata(
        manifest=manifest, cell_metadata=cell_metadata, cell_x_gene=cell_x_gene, cell_x_protein=cell_x_protein
    )

    # 1: Befor we modify the adata object, we write cell metadata with QC metrics
    smp.reroute_source(sd.CellMetadata, out_dir, overwrite=overwrite)
    sc_meta = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt32})
    sc_meta.write_csv(smp.out.CellMetadata.p, compression='gzip')

    # 1. Filter AnnData object
    adata_init = adata.copy()
    adata, cell_summary, gene_summary = sc_module.filtering.filter_adata(
        adata=adata, filter_panel=filter_panel, logger=log
    )

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning('No cells or genes passed the filtering criteria.')
        sc_module.process.write_dummys(adata=adata_init, failure_code='filter_not_passed')
        return
    del adata_init

    if smp.src.pr_detected and not omit_correlation:
        pr_corr_df, rna_pr_corr_df = sc_module.correlation.run_correlation_analysis(adata, logger=log)
        pr_corr_df.to_csv(smp.out.AdataH5.p.parent / 'protein_sc_correlation.csv')
        rna_pr_corr_df.to_csv(smp.out.AdataH5.p.parent / 'rna_protein_sc_correlation.csv')

    # 2. Pre-Processings (CPU/GPU) split path
    try:
        adata = sc_module.process.pre_process_adata(
            adata=adata, n_neighbors=n_neighbors, compute_backend=backend, rnd_st=rnd_st, logger=log
        )
    except Exception as e:
        log.warning(f'Preprocessing failed: {e}')
        sc_module.process.write_dummys(adata=adata, failure_code='preprocessing_failed')
        return

    # 3. Optimize Leiden clusters (CPU/GPU) split path
    success_clusterings = []
    for k, (target_clusters, init_res) in DEFAULT_CLUSTERINGS.items():
        try:
            adata = sc_module.cluster_dgex.optimize_leiden_clusters(
                adata,
                cluster_name=k,
                target_clusters=target_clusters,
                init_res=init_res,
                max_attempts=cluster_attempts,
                compute_backend=backend,
                rnd_st=rnd_st,
                logger=log,
            )
            success_clusterings.append(k)

        except Exception as e:
            log.warning(f'Failed to optimize clusters for {k}: {e}')

    if success_clusterings == []:
        log.warning('No successful clusterings to run differential gene expression analysis.')
        sc_module.process.write_dummys(adata=adata, failure_code='clustering_failed')
        return

    # Move AnnData object to CPU for downstream processing
    if backend.use_gpu:
        backend.rsc.get.anndata_to_CPU(adata)

    # 4. Generate UMAP and clustering dataframes
    umap_df = pl.from_numpy(adata.obsm['X_umap'], schema={'UMAP1': pl.Float32, 'UMAP2': pl.Float32})
    leiden_df = pl.from_pandas(adata.obs[success_clusterings], include_index=True)
    clustering_umap = leiden_df.hstack(umap_df)

    log_with_path(f'Writing {smp.out.ClusteringUmap.name} table:', smp.out.ClusteringUmap.p)
    clustering_umap.write_csv(smp.out.ClusteringUmap.p, compression='gzip')

    # 5. Run differential gene expression analysis
    try:
        dgex = run_dgex(adata, cluster_keys=success_clusterings, downsample=1000, logger=log)
        log_with_path(f'Writing {smp.out.Dgex.name} table:', smp.out.Dgex.p)
        dgex.write_csv(smp.out.Dgex.p, compression='gzip')
    except Exception as e:
        log.warning(f'Failed to run differential gene expression analysis: {e}')
        write_dummys(adata=adata, failure_code='dgex_failed')

    log_with_path(f'Writing {smp.out.AdataH5.name} h5ad:', smp.out.AdataH5.p)
    adata.write(smp.out.AdataH5.p)


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
