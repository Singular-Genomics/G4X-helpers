from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix

from .. import c, io
from ..schema.definition import CellMetadata, CellxGene, CellxProtein, TranscriptPanel

if TYPE_CHECKING:
    from anndata import AnnData

    from ..g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)

DEFAULT_CLUSTERINGS = {'leiden_coarse': (6, 0.25), 'leiden_fine': (12, 0.5)}

DEFAULT_FILTER_CONFIG = [
    dict(filter_type='cells', alias='area', key='nuclei_area_um', quantile_min=0.01, quantile_max=0.995),
    dict(filter_type='cells', alias='counts', key='total_counts', value_min=10, quantile_max=0.995),
    dict(filter_type='cells', alias='genes', key='n_genes_by_counts', value_min=5),
    dict(filter_type='cells', alias='controls', key='pct_counts_ctrl', value_max=5),
    dict(filter_type='genes', alias='targeting', key='probe_type', subset='targeting'),
    dict(filter_type='genes', alias='coverage', key='n_cells_by_counts', value_min=100),
]


def process_sc_output(
    g4x_obj,
    adata,
    *,
    cluster_attempts: int = 10,
    rnd_st: int = 111,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    # 1. Filter AnnData object
    adata_init = adata.copy()
    adata = filter_adata(adata=adata_init, logger=log)

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning('No cells or genes passed the filtering criteria.')
        write_dummy_clustering_outputs(g4x_obj, adata=adata_init, failure_code='filter_not_passed')
        return
    del adata_init

    # 2. Pre-Processings (CPU/GPU) split path
    backend = io.get_backend(which=compute_backend)

    try:
        adata = pre_process_adata(adata=adata, compute_backend=backend, logger=log)
    except Exception as e:
        log.warning(f'Preprocessing failed: {e}')
        write_dummy_clustering_outputs(g4x_obj, adata=adata, failure_code='preprocessing_failed')
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
        write_dummy_clustering_outputs(g4x_obj, adata=adata, failure_code='clustering_failed')
        return

    # 4. Generate UMAP and clustering dataframes
    umap_df = pl.from_numpy(adata.obsm['X_umap'], schema={'UMAP1': pl.Float32, 'UMAP2': pl.Float32})
    leiden_df = pl.from_pandas(adata.obs[success_clusterings], include_index=True)
    clustering_umap = leiden_df.hstack(umap_df)
    clustering_umap.write_csv(g4x_obj.tree.ClusteringUmap.p, compression='gzip')

    # 5. Run differential gene expression analysis
    try:
        backend.rsc.get.anndata_to_CPU(adata)
        dgex = run_dgex(adata, cluster_keys=success_clusterings, downsample=1000)
        dgex.write_csv(g4x_obj.tree.Dgex.p, compression='gzip')
    except Exception as e:
        log.warning(f'Failed to run differential gene expression analysis: {e}')
        write_dummy_clustering_outputs(g4x_obj, adata=adata, failure_code='dgex_failed')


# region initialize
def init_adata(
    g4x_obj: 'G4Xoutput',
    *,
    tx_panel: pl.DataFrame | str = '__g4x_default__',
    cell_metadata: pl.DataFrame | str = '__g4x_default__',
    cell_x_gene: pl.DataFrame | str = '__g4x_default__',
    cell_x_protein: pl.DataFrame | None = None,
    logger: logging.Logger | None = None,
) -> 'AnnData':

    log = logger or LOGGER
    log.info('Creating AnnData object')

    log.debug('Validating inputs')
    if tx_panel == '__g4x_default__':
        tx_panel = g4x_obj.tree.TranscriptPanel
    else:
        tx_panel = TranscriptPanel(target_path=tx_panel)

    if cell_metadata == '__g4x_default__':
        cell_meta_obj = g4x_obj.tree.CellMetadata
    else:
        cell_meta_obj = CellMetadata(target_path=cell_metadata)

    if cell_x_gene == '__g4x_default__':
        cell_x_gene = g4x_obj.tree.CellxGene
    else:
        cell_x_gene = CellxGene(target_path=cell_x_gene)

    if cell_x_protein is None:
        pass
    elif cell_x_protein == '__g4x_default__':
        cell_x_protein = g4x_obj.tree.CellxProtein
    else:
        cell_x_protein = CellxProtein(target_path=cell_x_protein)

    all_valid = (
        tx_panel.is_valid
        and cell_meta_obj.is_valid
        and cell_x_gene.is_valid
        and (cell_x_protein is None or cell_x_protein.is_valid)
    )
    if not all_valid:
        raise ValueError('One or more input dataframes are invalid')

    tx_panel = tx_panel.parse()
    cell_metadata = cell_meta_obj.load()
    cell_x_gene = cell_x_gene.load()
    if cell_x_protein is not None:
        cell_x_protein = cell_x_protein.load()

    cxg_ok = cell_metadata[c.CELL_ID_NAME].equals(cell_x_gene[c.CELL_ID_NAME])
    if not cxg_ok:
        raise ValueError('Cell IDs in CellMetadata do not match those in CellxGene')

    if cell_x_protein is not None:
        cxp_ok = cell_metadata[c.CELL_ID_NAME].equals(cell_x_protein[c.CELL_ID_NAME])
        if not cxp_ok:
            raise ValueError('Cell IDs in CellMetadata do not match those in CellxProtein')
        cell_x_protein = cell_x_protein.drop(c.CELL_ID_NAME)

    log.debug('Converting cell x gene matrix to sparse format')
    X = cell_x_gene.drop(c.CELL_ID_NAME).to_numpy().astype(np.uint16)
    X = csr_matrix(X)

    log.info('Processing cell metadata')
    obs_df = cell_metadata.to_pandas().set_index(c.CELL_ID_NAME)
    obs_df.index = obs_df.index.astype(str)

    log.info('Processing gene metadata')
    gene_ids = pl.Series(name=c.GENE_ID_NAME, values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    log.debug('Processing panel type information')
    panel_type = tx_panel.unique('gene_name').select('gene_name', 'probe_type').rename({'gene_name': c.GENE_ID_NAME})

    var_df = (
        var_df.join(
            panel_type,
            on=c.GENE_ID_NAME,
            how='left',
        )
        .to_pandas()
        .set_index(c.GENE_ID_NAME)
    )
    var_df.index = var_df.index.astype(str)

    log.info('Initializing AnnData object')
    adata = AnnData(X=X, obs=obs_df, var=var_df)

    if cell_x_protein is not None:
        log.debug('Adding protein data to AnnData object')
        adata.uns['protein_names'] = [c.removesuffix('_intensity_mean') for c in cell_x_protein.columns]
        adata.obsm['protein'] = cell_x_protein.to_numpy()

    log.info('Calculating QC metrics')
    adata.var['ctrl'] = adata.var['probe_type'] != 'targeting'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ctrl'], inplace=True, percent_top=None)
    adata.var.drop(columns=['ctrl'], inplace=True)

    sc_meta = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt64})
    sc_meta.write_csv(cell_meta_obj.p, compression='gzip')

    adata = _sanitize_categorical_columns(adata, threshold=10)
    return adata.copy()


# region processing
def filter_adata(
    adata: 'AnnData',
    filter_panel: 'FilterPanel' | None = None,
    *,
    logger: logging.Logger | None = None,
):

    log = logger or LOGGER

    log.info('Filtering adata')
    n_obs = adata.n_obs
    n_var = adata.n_vars

    if filter_panel is None:
        filter_panel = _get_default_filter_panel()
    cell_summary, gene_summary = filter_panel.filter(adata, apply=True)

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning('No cells or genes remaining after filtering. Returning empty AnnData object.')
        return adata  # , cell_summary, gene_summary

    for df, name, n_total in [(cell_summary, 'cells', n_obs), (gene_summary, 'genes', n_var)]:
        n_retained = df.filter(pl.all_horizontal(pl.col('^.*_ok$'))).select('n_total').item()
        retained = n_retained / n_total
        log.info('Retained {:,} ({:.2%}) {} after filtering'.format(n_retained, retained, name))

    return adata  # , cell_summary, gene_summary


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

    steps = [
        ('normalizing total counts', lambda: pp.normalize_total(adata)),
        ('log-transforming data', lambda: pp.log1p(adata)),
        ('computing PCA', lambda: pp.pca(adata, **pca_params)),
        ('computing neighbors', lambda: pp.neighbors(adata, **neighbors_kwargs)),
        ('running UMAP', lambda: tl.umap(adata, **umap_kwargs)),
    ]

    for message, fn in steps:
        log.info(message)
        fn()

    adata.uns['pre_process_method'] = backend_proxy
    adata.uns['is_gpu'] = compute_backend.use_gpu
    return adata


def optimize_leiden_clusters(
    adata: 'AnnData',
    *,
    cluster_name='leiden_clusters',
    cluster_prefix='C',
    target_clusters=6,
    init_res=0.5,
    max_attempts=15,
    rnd_st=777,
    compute_backend: io.ComputeBackend,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

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

    adata.obs = _reorder_clusters_by_size(adata.obs, key=cluster_name, prefix=cluster_prefix)

    adata.uns['clustering_method'] = backend_proxy
    adata.uns['is_gpu'] = compute_backend.use_gpu
    return adata


def run_dgex(
    adata: 'AnnData',
    cluster_keys: list[str] = ['leiden'],
    downsample: int = 1000,
    logger: logging.Logger | None = None,
) -> pl.DataFrame:
    log = logger or LOGGER

    dfList = []
    for leiden in cluster_keys:
        log.info('Running dgex on %s:', leiden)
        adata_downsample = _subset_adata_by_key(adata, key=leiden, n_per_group=downsample)
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


# region utilities
class FilterMethod:
    def __init__(
        self,
        filter_type: Literal['cells', 'genes'] = 'cells',
        *,
        key: str,
        alias: str | None = None,
        value_min: float | None = None,
        value_max: float | None = None,
        quantile_min: float | None = None,
        quantile_max: float | None = None,
        subset: str | None = None,
    ):
        if filter_type not in {'cells', 'genes'}:
            raise ValueError("filter_type must be 'cells' or 'genes'")

        if value_min is not None and value_max is not None and value_min > value_max:
            raise ValueError('value_min cannot be greater than value_max')

        if quantile_min is not None and quantile_max is not None and quantile_min > quantile_max:
            raise ValueError('quantile_min cannot be greater than quantile_max')

        if value_min is not None and quantile_min is not None:
            raise ValueError('Cannot specify both value_min and quantile_min')

        if value_max is not None and quantile_max is not None:
            raise ValueError('Cannot specify both value_max and quantile_max')

        self.key = key
        self.filter_type = filter_type
        self.value_min = value_min
        self.value_max = value_max
        self.quantile_min = quantile_min
        self.quantile_max = quantile_max
        self.subset = subset
        self.alias = alias if alias is not None else key

    @property
    def axis(self) -> str:
        return 'obs' if self.filter_type == 'cells' else 'var'

    def resolve_thresholds(self, df):
        v_min = self.value_min
        v_max = self.value_max

        if self.quantile_min is not None:
            v_min = df[self.key].quantile(self.quantile_min)

        if self.quantile_max is not None:
            v_max = df[self.key].quantile(self.quantile_max)

        return v_min, v_max

    def filter(self, adata, apply: bool = False):
        df = adata.obs if self.filter_type == 'cells' else adata.var

        if self.subset is not None:
            df = df[df[self.key] == self.subset]

        v_min, v_max = self.resolve_thresholds(df)
        return _filter_axis(adata, axis=self.axis, key=self.key, val_min=v_min, val_max=v_max, apply=apply)


class FilterPanel:
    def __init__(self, filters: list[FilterMethod]):
        self.filters = filters

    def _summarize_filtering(self, n_total: int, results: dict) -> pl.DataFrame:
        series = [pl.Series(name, values) for name, values in results.items()]
        df = pl.DataFrame(series).with_row_index()

        df = df.group_by(results.keys()).agg(pl.len().alias('n_total')).sort('n_total', descending=True)
        df = df.with_columns(((pl.col('n_total') / n_total) * 100).alias('pct'))
        return df

    def filter(self, adata, return_masks: bool = False, apply: bool = False):
        cell_results = {}
        gene_results = {}
        for mth in self.filters:
            if mth.filter_type == 'cells':
                cell_results[f'{mth.alias}_ok'] = mth.filter(adata, apply=False)
                all_passed_cell = np.logical_and.reduce(list(cell_results.values()))
            else:
                gene_results[f'{mth.alias}_ok'] = mth.filter(adata, apply=False)
                all_passed_gene = np.logical_and.reduce(list(gene_results.values()))

        cell_summary = self._summarize_filtering(adata.n_obs, cell_results)
        gene_summary = self._summarize_filtering(adata.n_vars, gene_results)

        if apply:
            adata._inplace_subset_obs(all_passed_cell)
            adata._inplace_subset_var(all_passed_gene)

        if return_masks:
            return (all_passed_cell, cell_summary), (all_passed_gene, gene_summary)

        return cell_summary, gene_summary


def _filter_axis(
    adata,
    axis: str,
    key: str,
    val_min: float | None = None,
    val_max: float | None = None,
    apply: bool = False,
):
    if axis == 'obs':
        df = adata.obs
        size = adata.n_obs
    elif axis == 'var':
        df = adata.var
        size = adata.n_vars
    else:
        raise ValueError("axis must be 'obs' or 'var'")

    if key not in df:
        raise ValueError(f"key '{key}' not found in adata.{axis}")

    values = df[key].to_numpy()
    mask = np.ones(size, dtype=bool)

    if val_min is not None:
        mask &= values >= val_min
    if val_max is not None:
        mask &= values <= val_max

    if apply:
        return adata[mask].copy() if axis == 'obs' else adata[:, mask].copy()
    return mask


def _reorder_clusters_by_size(obs: pd.DataFrame, key: str, prefix: str = 'C') -> pd.DataFrame:
    # Ensure consistent grouping behavior
    clust_sorted = (
        obs.groupby(key, observed=False)  # explicit to silence warning
        .size()
        .sort_values(ascending=False)
    )

    size_order = clust_sorted.index.tolist()

    # Create new labels
    numeric_order = [f'{prefix}{i}' for i in range(len(size_order))]

    # Mapping old -> new labels
    order_map = dict(zip(size_order, numeric_order))

    reordered_obs = obs.copy()

    # Handle categorical safely
    if isinstance(reordered_obs[key].dtype, pd.CategoricalDtype):
        # Convert to string (or object) before replacing
        reordered_obs[key] = reordered_obs[key].astype(str)

    reordered_obs[key] = reordered_obs[key].replace(order_map)

    reordered_obs[key] = pd.Categorical(
        reordered_obs[key],
        categories=numeric_order,  # this defines the order!
        ordered=True,
    )

    return reordered_obs


def write_dummy_clustering_outputs(
    g4x_obj,
    *,
    adata,
    out_dir: str = '__g4x_default__',
    success_clusterings: list[str] = [],
    failure_code: str = 'failed',
) -> None:
    out_dir = g4x_obj.data_dir if out_dir == '__g4x_default__' else Path(out_dir)

    if 'X_umap' in adata.obsm:
        print('X_umap was found')
        umap = adata.obsm['X_umap']
    else:
        print('X_umap was not found, filling with NaNs')
        umap = np.full((adata.n_obs, 2), np.nan)

    obs_df = pl.from_pandas(adata.obs, include_index=True).cast({c.CELL_ID_NAME: pl.UInt64})
    dummy_clust = obs_df.select([c.CELL_ID_NAME] + success_clusterings)

    umap = adata.obsm['X_umap']

    dummy_clust = dummy_clust.with_columns(
        UMAP1=umap[:, 0],
        UMAP2=umap[:, 1],
    )

    if not success_clusterings:
        dummy_clust = dummy_clust.with_columns(leiden=pl.lit('unassigned'))
        success_clusterings.append('leiden')

    dummy_clust = dummy_clust.with_columns(failure_code=pl.lit(failure_code))

    col_select = [c.CELL_ID_NAME] + success_clusterings + ['UMAP1', 'UMAP2', 'failure_code']
    dummy_clust = dummy_clust.select(col_select)
    dummy_clust.write_csv(out_dir / g4x_obj.tree.ClusteringUmap.target_rel, compression='gzip')

    dummy_dgex = pl.DataFrame(schema=g4x_obj.tree.Dgex.SCHEMA)
    null_row = pl.DataFrame({col: pl.Series([None], dtype=dummy_dgex.schema[col]) for col in dummy_dgex.columns})

    dummy_dgex = dummy_dgex.vstack(null_row).with_columns(failure_code=pl.lit('filter_not_passed'))
    dummy_dgex.write_csv(out_dir / g4x_obj.tree.Dgex.target_rel, compression='gzip')


def _sanitize_categorical_columns(adata, threshold: int = 10):
    for df in [adata.obs, adata.var]:
        for col_name in df.columns:
            col = df[col_name]
            if col.nunique() < threshold and col.dtype == 'O':
                df[col_name] = col.astype('category')
    return adata


def _get_default_filter_panel():
    default_filter_panel = FilterPanel(filters=[FilterMethod(**cfg) for cfg in DEFAULT_FILTER_CONFIG])
    return default_filter_panel


def _subset_adata_by_key(adata, key='leiden', n_per_group=1000):

    if key not in adata.obs:
        raise ValueError(f"Key '{key}' not found in adata.obs")

    rng = np.random.default_rng(42)

    sampled_idx = (
        adata.obs.groupby(key, group_keys=False, observed=True)
        .apply(
            lambda x: x.sample(n=min(len(x), n_per_group), random_state=rng.integers(1e9)),
            include_groups=False,  # <-- fixes deprecation warning
        )
        .index
    )

    return adata[sampled_idx].copy()


# polars implementation
# def reorder_clusters_by_size(clust_umap, key: str, prefix='C') -> pl.DataFrame:
#     clust_sorted = clust_umap.group_by(key).agg(pl.len()).sort('len', descending=True)
#     size_order = clust_sorted[key].to_list()
#     numeric_order = np.arange(len(size_order)).astype(str).tolist()

#     numeric_order = [prefix + str(uid) for uid in numeric_order]
#     order_map = dict(zip(size_order, numeric_order))

#     new_umap = clust_umap.with_columns(pl.col(key).replace(order_map)).sort('seg_cell_id')
#     return new_umap


# def create_adata():
# --- a more fine grained breakdown of control metrics ---
# ctrl_types = ['NCP', 'NCS', 'GCP']
# for ctrl_type in ctrl_types:
#     adata.var[ctrl_type.lower()] = adata.var['probe_type'] == ctrl_type

# ctrl_types_lower = [c.lower() for c in ctrl_types]
# sc.pp.calculate_qc_metrics(adata, qc_vars=ctrl_types_lower, inplace=True, percent_top=None)

# adata.var.drop(columns=ctrl_types_lower, inplace=True)

# --- also requires different cell filtering ---
# pct_counts_ncp = filter_by_limits(adata, axis='obs', key='pct_counts_ncp', val_min=0, val_max=5)
# pct_counts_ncs = filter_by_limits(adata, axis='obs', key='pct_counts_ncs', val_min=0, val_max=5)
# pct_counts_gcp = filter_by_limits(adata, axis='obs', key='pct_counts_gcp', val_min=0, val_max=5)
# controls_passed = pct_counts_ncp & pct_counts_ncs & pct_counts_gcp
