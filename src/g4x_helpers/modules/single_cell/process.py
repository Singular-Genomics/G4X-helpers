from __future__ import annotations

import logging
import traceback
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix

from ... import c, io
from ...schema import definition as sd
from . import sc_utils
from .filter_panel import FilterPanel, _get_default_filter_panel

if TYPE_CHECKING:
    from anndata import AnnData


log = logging.getLogger(__name__)

DEFAULT_CLUSTERINGS = {'leiden_coarse': (6, 0.25), 'leiden_fine': (12, 0.5)}


def default_sc_pipeline(
    adata: 'AnnData',
    *,
    n_neighbors: int = 15,
    clusterings: dict[str, tuple[int, float]] = DEFAULT_CLUSTERINGS,
    cluster_attempts: int = 10,
    rnd_st: int = 111,
    compute_backend: io.ComputeBackend = io.get_backend(which='auto'),
) -> tuple['AnnData', str]:

    # 1. Filter AnnData object
    adata_init = adata.copy()

    try:
        adata, _, _ = filter_adata(adata=adata)
    except Exception as e:
        log.warning(f'Filtering failed: {e}')
        return adata_init, 'filter_not_passed'

    del adata_init

    # 2. Pre-Processings (CPU/GPU) split path
    try:
        adata = pre_process_adata(adata=adata, n_neighbors=n_neighbors, compute_backend=compute_backend, rnd_st=rnd_st)
    except Exception as e:
        log.warning(f'Preprocessing failed: {e}')
        return adata, 'preprocessing_failed'

    # 3. Optimize Leiden clusters (CPU/GPU) split path
    success_clusterings = []
    for k, (target_clusters, init_res) in clusterings.items():
        try:
            adata = optimized_clustering(
                adata,
                cluster_name=k,
                target_clusters=target_clusters,
                init_res=init_res,
                max_attempts=cluster_attempts,
                compute_backend=compute_backend,
                rnd_st=rnd_st,
            )
            success_clusterings.append(k)

        except Exception as e:
            log.warning(f'Failed to optimize clusters for {k}: {e}')

    if success_clusterings == []:
        log.warning('No successful clusterings to run differential gene expression analysis.')
        return adata, 'clustering_failed'

    # Move AnnData object to CPU for downstream processing
    if compute_backend.use_gpu:
        compute_backend.rsc.get.anndata_to_CPU(adata)

    return adata, 'success'


def init_adata(
    manifest: pl.DataFrame,
    cell_metadata: pl.DataFrame,
    cell_x_gene: pl.DataFrame,
    cell_x_protein: pl.DataFrame | None = None,
) -> 'AnnData':

    cell_metadata = cell_metadata.sort(c.CELL_ID_NAME)
    cell_metadata = cell_metadata.sort(c.CELL_ID_NAME)

    sanitize_cols = [
        'n_genes_by_counts',
        'log1p_n_genes_by_counts',
        'total_counts',
        'log1p_total_counts',
        'total_counts_ctrl',
        'log1p_total_counts_ctrl',
        'pct_counts_ctrl',
    ]
    for col in sanitize_cols:
        if col in cell_metadata.columns:
            cell_metadata = cell_metadata.drop(col)

    # 3: Ensure that cell IDs match between metadata and expression/protein matrices
    sc_utils._validate_cell_ids(cell_metadata, cell_x_gene, 'CellMetadata', 'CellxGene')

    if cell_x_protein is not None:
        cell_x_protein = cell_x_protein.sort(c.CELL_ID_NAME)
        sc_utils._validate_cell_ids(cell_metadata, cell_x_protein, 'CellMetadata', 'CellxProt')
        cell_x_protein = cell_x_protein.drop(c.CELL_ID_NAME)

    # 3: Set up AnnData object
    log.info('Setting up AnnData components')
    log.debug('Converting cell x gene matrix to sparse format')
    X = cell_x_gene.drop(c.CELL_ID_NAME).to_numpy().astype(np.uint16)
    X = csr_matrix(X)

    log.info('Processing metadata for cells and genes')
    obs_df = cell_metadata.to_pandas().set_index(c.CELL_ID_NAME)
    obs_df.index = obs_df.index.astype(str)

    gene_ids = pl.Series(name=c.GENE_ID_NAME, values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    log.debug('Processing panel type information')
    probe_type = manifest.unique('gene_name').select('gene_name', 'probe_type').rename({'gene_name': c.GENE_ID_NAME})

    var_df = (
        var_df.join(
            probe_type,
            on=c.GENE_ID_NAME,
            how='left',
        )
        .to_pandas()
        .set_index(c.GENE_ID_NAME)
    )
    var_df['probe_type'] = var_df['probe_type'].str.lower()
    var_df.index = var_df.index.astype(str)

    # 4: bring it all together
    log.debug('Initializing AnnData object')
    adata = AnnData(X=X, obs=obs_df, var=var_df)

    if cell_x_protein is not None:
        log.debug('Adding protein data to AnnData object')
        adata.uns['protein_names'] = [col.removesuffix(c.IMG_INTENSITY_HANDLE) for col in cell_x_protein.columns]
        adata.obsm['protein'] = cell_x_protein.to_numpy()

    # 5: Calculate QC metrics
    log.info('Calculating QC metrics')
    adata.var['ctrl'] = adata.var['probe_type'] != 'targeting'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ctrl'], inplace=True, percent_top=None)
    adata.var.drop(columns=['ctrl'], inplace=True)

    adata = sc_utils._sanitize_categorical_columns(adata, threshold=10)
    return adata.copy()


def filter_adata(
    adata: 'AnnData',
    filter_panel: 'FilterPanel' | None = None,
) -> tuple['AnnData', pl.DataFrame, pl.DataFrame]:

    log.info('Filtering adata')

    if filter_panel is None:
        filter_panel = _get_default_filter_panel()

    obs_total, var_total = adata.n_obs, sum(adata.var['probe_type'] == 'targeting')  # adata.n_vars
    cell_summary, gene_summary = filter_panel.filter(adata, apply=True)

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError('No cells or genes remaining after filtering')

    for df, name, n_total in [(cell_summary, 'cells', obs_total), (gene_summary, 'targeting genes', var_total)]:
        n_retained = df.filter(pl.all_horizontal(pl.col('^.*_ok$'))).select('n_total').item()
        retained = n_retained / n_total
        log.info('Retained {:,} ({:.2%}) {} after filtering'.format(n_retained, retained, name))

    return adata, cell_summary, gene_summary


def pre_process_adata(
    adata: 'AnnData',
    *,
    n_neighbors: int = 15,
    n_pcs: int = 15,
    pca_comps: int = 50,
    umap_key: str = 'X_umap',
    umap_min_dist: float = 0.15,
    rnd_st=777,
    compute_backend: io.ComputeBackend = io.get_backend(which='auto'),
):

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


def optimized_clustering(
    adata: 'AnnData',
    *,
    cluster_name='leiden_clusters',
    cluster_prefix='C',
    target_clusters=6,
    init_res=0.5,
    max_attempts=15,
    rnd_st=777,
    compute_backend: io.ComputeBackend = io.get_backend(which='auto'),
):

    backend_proxy = 'rapids (GPU)' if compute_backend.use_gpu else 'scanpy (CPU)'
    log.info('Optimizing Leiden clusters for: %s', cluster_name)

    results = {}

    low_res, high_res = None, None
    res = init_res

    for i in range(max_attempts):
        if compute_backend.use_gpu:
            compute_backend.rsc.tl.leiden(adata, resolution=res, key_added='tmp_leiden', random_state=rnd_st)
        else:
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added='tmp_leiden',
                flavor='igraph',
                n_iterations=-1,
                objective_function='modularity',
                random_state=rnd_st,
            )
        assignment = adata.obs['tmp_leiden'].copy()
        n_clusters = len(assignment.cat.categories)
        delta = n_clusters - target_clusters

        log.debug(f'resolution: {res:.4f}, n-clusters: {n_clusters}, d-ideal: {delta}')

        results[i] = {
            'resolution': res,
            'n_clusters': n_clusters,
            'delta_ideal': delta,
            'assignment': assignment,
        }

        if delta == 0:
            break

        if delta < 0:
            # too few clusters -> need higher resolution
            low_res = res
            if high_res is None:
                res *= 2
            else:
                res = (res + high_res) / 2
        else:
            # too many clusters -> need lower resolution
            high_res = res
            if low_res is None:
                res /= 2
            else:
                res = (low_res + res) / 2

        res = max(res, 1e-6)

    adata.obs.drop(columns=['tmp_leiden'], inplace=True, errors='ignore')

    best_key = min(
        results,
        key=lambda k: (
            abs(results[k]['delta_ideal']),
            abs(results[k]['resolution'] - init_res),
        ),
    )

    column = results[best_key].pop('assignment')
    adata.obs[cluster_name] = column

    results[best_key]['random_state'] = rnd_st
    results[best_key]['n_iterations'] = adata.uns['tmp_leiden']['params']['n_iterations']
    del adata.uns['tmp_leiden']

    adata.uns[cluster_name] = {'params': results[best_key]}

    adata.obs = sc_utils._reorder_clusters_by_size(adata.obs, key=cluster_name, prefix=cluster_prefix)

    adata.uns['clustering_method'] = backend_proxy
    adata.uns['is_gpu'] = compute_backend.use_gpu
    return adata


def run_dgex(adata: 'AnnData', cluster_keys: list[str] = ['leiden'], downsample: int = 1000) -> pl.DataFrame:

    dfList = []
    for leiden in cluster_keys:
        log.info('Running dgex on: %s', leiden)
        adata_downsample = sc_utils._subset_adata_by_key(adata, key=leiden, n_per_group=downsample)
        try:
            sc.tl.rank_genes_groups(
                adata_downsample,
                groupby=leiden,
                use_raw=False,
                method='wilcoxon',
                pts=True,
                key_added=f'{leiden}_rank_genes_groups',
            )
            for g in adata_downsample.obs[leiden].unique():
                tmp = sc.get.rank_genes_groups_df(adata_downsample, group=str(g), key=f'{leiden}_rank_genes_groups')
                tmp['group'] = g
                tmp['leiden_res'] = leiden
                dfList.append(pl.from_pandas(tmp))

        except Exception as e:
            log.warning(f'DGEX failed for {leiden}: {e}')
            log.debug(traceback.format_exc())

    dgex = pl.DataFrame()
    if len(dfList) > 0:
        dgex = pl.concat(dfList).rename(
            {
                'names': c.GENE_ID_NAME,
                'group': 'cluster_id',
                'scores': 'score',
                'logfoldchanges': 'logfoldchange',
                'pvals': 'pval',
                'pvals_adj': 'pval_adj',
            }
        )
        first_cols = ['leiden_res', 'cluster_id']
        last_cols = [c for c in dgex.columns if c not in first_cols]
        dgex = dgex.select(first_cols + last_cols)

    return dgex


def dummy_dgex_output(failure_code):
    dummy_dgex = pl.DataFrame(schema=sd.Dgex.SCHEMA)
    null_row = pl.DataFrame({col: pl.Series([None], dtype=dummy_dgex.schema[col]) for col in dummy_dgex.columns})
    dummy_dgex = dummy_dgex.vstack(null_row).with_columns(failure_code=pl.lit(failure_code))
    return dummy_dgex


def dummy_clustering_output(adata, failure_code: str = 'failed') -> pl.DataFrame:

    log.info('Filling missing outputs with placeholders due to failure code: %s', failure_code)

    # 1: Clustering / UMAP table
    if 'X_umap' not in adata.obsm:
        log.debug('X_umap was not found, filling with NaNs')
        umap = np.full((adata.n_obs, 2), np.nan)
    else:
        umap = adata.obsm['X_umap']

    # get the adata.obs
    obs_df = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt64})
    clusters = [col for col in obs_df.columns if col.startswith('leiden_')]
    dummy_clust = obs_df.select([c.CELL_ID_NAME] + clusters)

    dummy_clust = dummy_clust.with_columns(
        UMAP1=umap[:, 0],
        UMAP2=umap[:, 1],
    )

    # add a default clustering if none were successful
    if clusters == []:
        log.debug('No successful clusterings found, adding default label "unassigned" for all cells')
        dummy_clust = dummy_clust.with_columns(leiden=pl.lit('unassigned'))
        clusters.append('leiden')

    dummy_clust = dummy_clust.with_columns(failure_code=pl.lit(failure_code))

    # select the relevant columns for the output
    col_select = [c.CELL_ID_NAME] + clusters + ['UMAP1', 'UMAP2', 'failure_code']
    dummy_clust = dummy_clust.select(col_select)

    return dummy_clust
