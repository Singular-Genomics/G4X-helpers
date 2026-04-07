from __future__ import annotations

import logging
from functools import partial
from typing import TYPE_CHECKING, Literal

import numpy as np
import polars as pl
import scanpy as sc
from anndata import AnnData

from ... import c, io
from ... import logging_utils as logut
from ...schema.definition import AdataH5, CellMetadata, ClusteringUmap, Dgex
from ..workflow import PRESET_SOURCE, reroute_source
from .cluster_dgex import optimize_leiden_clusters, run_dgex
from .filtering import FilterPanel, _get_default_filter_panel, filter_adata
from .init_adata import init_adata

if TYPE_CHECKING:
    from anndata import AnnData

    from ...g4x_output import G4Xoutput


LOGGER = logging.getLogger(__name__)

DEFAULT_CLUSTERINGS = {'leiden_coarse': (6, 0.25), 'leiden_fine': (12, 0.5)}


# region main functions
def process_sc_output(
    smp: 'G4Xoutput',
    adata: 'AnnData' | None = None,
    out_dir: str = PRESET_SOURCE,
    *,
    overwrite: bool = False,
    filter_panel: 'FilterPanel' | None = None,
    n_neighbors: int = 15,
    cluster_attempts: int = 10,
    rnd_st: int = 111,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    backend = io.get_backend(which=compute_backend)

    if adata is None:
        log.info('Initializing AnnData object from raw data')
        adata = init_adata(smp, logger=log)
    else:
        log.info('Using provided AnnData object with %d cells and %d genes', adata.n_obs, adata.n_vars)

    # Validate and prepare output, set up reusable functions
    out_dir = smp.data_dir if out_dir == PRESET_SOURCE else io.pathval.validate_dir_path(out_dir)
    prep_out = partial(reroute_source, smp, out_dir, overwrite=overwrite, logger=log)
    log_with_path = partial(logut.log_with_path, logger=log, level='info')
    write_dummys = partial(write_dummy_clustering_outputs, smp=smp, logger=log)

    prep_out(validator=AdataH5)
    prep_out(validator=Dgex)
    prep_out(validator=ClusteringUmap)
    prep_out(validator=CellMetadata, overwrite=True)

    # 1: Befor we modify the adata object, we write cell metadata with QC metrics
    logut.log_with_path(f'Writing {smp.out.CellMetadata.name} table:', smp.out.CellMetadata.p, level='info', logger=log)
    sc_meta = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt64})
    sc_meta.write_csv(smp.out.CellMetadata.p, compression='gzip')

    # 1. Filter AnnData object
    if filter_panel is None:
        filter_panel = _get_default_filter_panel()

    adata_init = adata.copy()
    adata, cell_summary, gene_summary = filter_adata(adata=adata_init, filter_panel=filter_panel, logger=log)

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning('No cells or genes passed the filtering criteria.')
        write_dummys(adata=adata_init, failure_code='filter_not_passed')
        return
    del adata_init

    # 2. Pre-Processings (CPU/GPU) split path
    try:
        adata = pre_process_adata(adata=adata, n_neighbors=n_neighbors, compute_backend=backend, logger=log)
    except Exception as e:
        log.warning(f'Preprocessing failed: {e}')
        write_dummys(adata=adata, failure_code='preprocessing_failed')
        return

    # 3. Optimize Leiden clusters (CPU/GPU) split path
    success_clusterings = []
    for k, (target_clusters, init_res) in DEFAULT_CLUSTERINGS.items():
        try:
            adata = optimize_leiden_clusters(
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

    if not success_clusterings:
        log.warning('No successful clusterings to run differential gene expression analysis.')
        write_dummys(adata=adata, failure_code='clustering_failed')
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


# region higher-level functions
def pre_process_adata(
    adata: 'AnnData',
    *,
    n_neighbors: int = 15,
    n_pcs: int = 15,
    pca_comps: int = 50,
    umap_key: str = 'X_umap',
    umap_min_dist: float = 0.15,
    rnd_st=777,
    compute_backend: io.ComputeBackend,
    logger: logging.Logger | None = None,
):

    log = logger or LOGGER
    backend_proxy = 'rapids (GPU)' if compute_backend.use_gpu else 'scanpy (CPU)'

    log.info('Pre-processing AnnData object using %s', backend_proxy)

    log.debug('Storing raw counts in adata.layers["counts"]')
    adata.layers['counts'] = adata.X.copy()

    pca_params = {'svd_solver': 'auto', 'n_comps': pca_comps, 'random_state': rnd_st}
    neighbors_params = {'n_neighbors': n_neighbors, 'n_pcs': n_pcs, 'random_state': rnd_st}
    umap_params = {'key_added': umap_key, 'min_dist': umap_min_dist, 'spread': 1, 'random_state': rnd_st}

    if compute_backend.use_gpu:
        rsc = compute_backend.rsc

        log.debug('Moving adata to GPU.')
        adata.X = adata.X.astype(float)
        rsc.get.anndata_to_GPU(adata)

        pp = rsc.pp
        tl = rsc.tl
        neighbors_kwargs = {'algorithm': 'brute', **neighbors_params}
        umap_kwargs = {'init_pos': 'random', **umap_params}
    else:
        pp = sc.pp
        tl = sc.tl
        neighbors_kwargs = neighbors_params
        umap_kwargs = {'init_pos': 'spectral', **umap_params}
        pca_params['svd_solver'] = 'arpack'  # more stable for small datasets on CPU

    steps = [
        ('normalizing total counts', lambda: pp.normalize_total(adata)),
        ('log-transforming data', lambda: pp.log1p(adata)),
        ('computing PCA', lambda: pp.pca(adata, **pca_params)),
        ('computing neighbors', lambda: pp.neighbors(adata, **neighbors_kwargs)),
        ('running UMAP', lambda: tl.umap(adata, **umap_kwargs)),
    ]

    log.info('Normalize -> Log-transform -> PCA -> Neighbors -> UMAP')
    for message, fn in steps:
        log.debug(message)
        fn()

    adata.uns['pre_process_method'] = backend_proxy
    adata.uns['is_gpu'] = compute_backend.use_gpu
    return adata


# region core functions
def write_dummy_clustering_outputs(
    smp,
    *,
    adata,
    success_clusterings: list[str] = [],
    failure_code: str = 'failed',
    logger: logging.Logger | None = None,
) -> None:

    log = logger or LOGGER
    log.info('Filling missing outputs with placeholders due to failure code: %s', failure_code)

    # 1: Clustering / UMAP table
    if 'X_umap' not in adata.obsm:
        log.debug('X_umap was not found, filling with NaNs')
        umap = np.full((adata.n_obs, 2), np.nan)
    else:
        umap = adata.obsm['X_umap']

    # get the adata.obs
    obs_df = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt64})
    dummy_clust = obs_df.select([c.CELL_ID_NAME] + success_clusterings)

    dummy_clust = dummy_clust.with_columns(
        UMAP1=umap[:, 0],
        UMAP2=umap[:, 1],
    )

    # add a default clustering if none were successful
    if not success_clusterings:
        log.debug('No successful clusterings found, adding default label "unassigned" for all cells')
        dummy_clust = dummy_clust.with_columns(leiden=pl.lit('unassigned'))
        success_clusterings.append('leiden')

    dummy_clust = dummy_clust.with_columns(failure_code=pl.lit(failure_code))

    # select the relevant columns for the output
    col_select = [c.CELL_ID_NAME] + success_clusterings + ['UMAP1', 'UMAP2', 'failure_code']
    dummy_clust = dummy_clust.select(col_select)

    # 2: Dgex table <- this is the last step in the pipeline. so when the previous steps fail, we still need to create a placeholder
    dummy_dgex = pl.DataFrame(schema=smp.src.Dgex.SCHEMA)
    null_row = pl.DataFrame({col: pl.Series([None], dtype=dummy_dgex.schema[col]) for col in dummy_dgex.columns})
    dummy_dgex = dummy_dgex.vstack(null_row).with_columns(failure_code=pl.lit(failure_code))

    msg = str(smp.out.ClusteringUmap.p) + '\n'
    msg += str(smp.out.Dgex.p) + '\n'
    msg += str(smp.out.AdataH5.p)

    logut.log_msg_wrapped(header='Writing clustering outputs with placeholders', msg=msg, prefix=logut.PGAP, logger=log)
    dummy_clust.write_csv(smp.out.ClusteringUmap.p, compression='gzip')
    dummy_dgex.write_csv(smp.out.Dgex.p, compression='gzip')
    adata.write(smp.out.AdataH5.p)


# polars implementation
# def reorder_clusters_by_size(clust_umap, key: str, prefix='C') -> pl.DataFrame:
#     clust_sorted = clust_umap.group_by(key).agg(pl.len()).sort('len', descending=True)
#     size_order = clust_sorted[key].to_list()
#     numeric_order = np.arange(len(size_order)).astype(str).tolist()

#     numeric_order = [prefix + str(uid) for uid in numeric_order]
#     order_map = dict(zip(size_order, numeric_order))

#     new_umap = clust_umap.with_columns(pl.col(key).replace(order_map)).sort('seg_cell_id')
#     return new_umap
