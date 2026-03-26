import logging
import traceback
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix

from .. import c, io

if TYPE_CHECKING:
    from anndata import AnnData

    from ..g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)


def create_adata(
    g4x_obj: 'G4Xoutput',
    *,
    cluster_attempts: int = 10,
    rnd_st: int = 111,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    tx_panel = g4x_obj.tree.TranscriptPanel.parse()
    cell_metadata = pl.read_csv(g4x_obj.tree.CellMetadata.p)

    # prot = pl.read_csv(g4x_obj.tree.CellxProtein.p).drop('seg_cell_id')#.sample(fraction=1, shuffle=True)
    # cell_metadata = cell_metadata.hstack(prot)

    cell_x_gene = pl.read_csv(g4x_obj.tree.CellxGene.p)

    adata = init_adata(tx_panel=tx_panel, cell_metadata=cell_metadata, cell_x_gene=cell_x_gene, logger=log)

    backend = io.get_backend(which='auto')
    adata = pre_process_adata(
        adata, rnd_st=rnd_st, pca_comps=50, n_neighbors=15, n_pcs=15, compute_backend=backend, logger=log
    )

    cluster_set_names = ['leiden_coarse', 'leiden_fine']
    adata = optimize_leiden_clusters(
        adata,
        cluster_name=cluster_set_names[0],
        target_clusters=6,
        init_res=0.25,
        max_attempts=cluster_attempts,
        compute_backend=backend,
        rnd_st=rnd_st,
        logger=log,
    )
    adata = optimize_leiden_clusters(
        adata,
        cluster_name=cluster_set_names[1],
        target_clusters=12,
        init_res=0.5,
        max_attempts=cluster_attempts,
        compute_backend=backend,
        rnd_st=rnd_st,
        logger=log,
    )

    if backend.use_gpu:
        log.debug('Fetching adata back to CPU to run DGEX.')
        backend.rsc.get.anndata_to_CPU(adata)

    for df in [adata.obs, adata.var]:
        for col_name in df.columns:
            col = df[col_name]
            if col.nunique() == 1 and col.dtype == 'O':
                df[col_name] = col.astype('category')

    umap_df = pl.from_numpy(adata.obsm['X_umap'], schema={'UMAP1': pl.Float32, 'UMAP2': pl.Float32})
    leiden_df = pl.from_pandas(adata.obs[cluster_set_names], include_index=True)

    clustering_umap = leiden_df.hstack(umap_df)
    clustering_umap.write_csv(g4x_obj.tree.ClusteringUmap.p, compression='gzip')

    adata.write_h5ad(g4x_obj.tree.AdataH5.p)

    dgex = run_dgex(adata, downsample=1000, logger=log)
    dgex.write_csv(g4x_obj.tree.Dgex.p, compression='gzip')

    return adata


# TODO addition of protein values is missing
# region create adata
def init_adata(
    tx_panel: pl.DataFrame,
    cell_x_gene: pl.DataFrame,
    cell_metadata: pl.DataFrame | None = None,
    *,
    filter_obs_area: dict[str, float] = {'q_min': 0.01, 'q_max': 0.995},
    filter_obs_counts: dict[str, float] = {'val_min': 0, 'q_max': 0.995},
    filter_obs_genes: dict[str, float] = {'val_min': 5, 'val_max': None},
    filter_obs_controls: dict[str, float] = {'val_min': 0, 'val_max': 5},
    filter_var_cells: dict[str, float] = {'val_min': 100, 'val_max': None},
    logger: logging.Logger | None = None,
) -> 'AnnData':

    log = logger or LOGGER
    log.info('Creating AnnData object')

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

    log.info('Calculating QC metrics')
    adata.var['ctrl'] = adata.var['probe_type'] != 'targeting'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ctrl'], inplace=True, percent_top=None)
    adata.var.drop(columns=['ctrl'], inplace=True)

    log.info('Filtering adata')
    n_obs = adata.n_obs
    n_var = adata.n_vars

    cell_results = filter_cells(
        adata, area=filter_obs_area, counts=filter_obs_counts, genes=filter_obs_genes, controls=filter_obs_controls
    )
    gene_results = filter_genes(adata, cells=filter_var_cells)

    for df, name, n_total in [(cell_results, 'cells', n_obs), (gene_results, 'genes', n_var)]:
        n_retained = df.filter(pl.all_horizontal(pl.col('^.*_ok$'))).select('n_total').item()
        retained = n_retained / n_total
        log.info('Retained {:,} ({:.2%}) {} after filtering'.format(n_retained, retained, name))

    return adata.copy()


# region processing
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
    downsample: int = 1000,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    # take at most {downsample} cells per cluster, but never more than the smallest cluster
    smallest_cluster = adata.obs['leiden_fine'].value_counts().min()
    min_cells = np.min([downsample, smallest_cluster])

    idx = adata.obs.groupby('leiden_fine', observed=False).sample(min_cells).index
    adata_downsample = adata[idx].copy()

    dfList = []
    for leiden in ['leiden_coarse', 'leiden_fine']:
        log.info('Running dgex on %s:', leiden)
        try:
            # with warnings.catch_warnings():
            #     warnings.simplefilter('ignore', PerformanceWarning)
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


# region filtering
def filter_cells(
    adata: 'AnnData',
    area: dict[str, float],
    counts: dict[str, float],
    genes: dict[str, float],
    controls: dict[str, float],
) -> pl.DataFrame:
    n_obs = adata.n_obs
    size_passed = _filter_by_qtiles(adata, key='nuclei_area_um', axis='obs', q_min=area['q_min'], q_max=area['q_max'])

    max_counts = _filter_by_qtiles(adata, axis='obs', key='total_counts', q_min=0, q_max=counts['q_max'])
    min_counts = _filter_by_limits(adata, axis='obs', key='total_counts', val_min=counts['val_min'], val_max=None)
    counts_passed = max_counts & min_counts

    genes_passed = _filter_by_limits(
        adata, axis='obs', key='n_genes_by_counts', val_min=genes['val_min'], val_max=genes['val_max']
    )

    controls_passed = _filter_by_limits(
        adata, axis='obs', key='pct_counts_ctrl', val_min=controls['val_min'], val_max=controls['val_max']
    )

    all_passed = counts_passed & genes_passed & controls_passed & size_passed
    adata._inplace_subset_obs(all_passed)

    results = {
        'counts_ok': counts_passed,
        'genes_ok': genes_passed,
        'controls_ok': controls_passed,
        'size_ok': size_passed,
    }
    cell_results = _summarize_filtering(n_obs, results)

    return cell_results


def filter_genes(adata: 'AnnData', cells: dict[str, float]) -> pl.DataFrame:
    n_var = adata.n_vars
    targeting_passed = adata.var['probe_type'] == 'targeting'
    coverage_passed = _filter_by_limits(
        adata, axis='var', key='n_cells_by_counts', val_min=cells['val_min'], val_max=cells['val_max']
    )
    genes_all_passed = targeting_passed & coverage_passed
    adata._inplace_subset_var(genes_all_passed)

    results = {'targeting_ok': targeting_passed, 'coverage_ok': coverage_passed}
    gene_results = _summarize_filtering(n_var, results)

    return gene_results


# region utilities
def _filter_by_limits(
    adata: 'AnnData',
    key: str,
    axis: str = 'obs',
    val_min: float | None = None,
    val_max: float | None = None,
    apply: bool = False,
) -> np.ndarray:
    return _filter_axis(adata=adata, axis=axis, key=key, val_min=val_min, val_max=val_max, apply=apply)


def _filter_by_qtiles(
    adata: 'AnnData', key: str, axis: str = 'obs', q_min: float = 0, q_max: float = 1, apply: bool = False
) -> np.ndarray:
    if not (0 <= q_min <= 1 and 0 <= q_max <= 1):
        raise ValueError('q_min and q_max must be between 0 and 1')
    if q_min > q_max:
        raise ValueError('q_min must be <= q_max')

    df = adata.obs if axis == 'obs' else adata.var
    if axis not in {'obs', 'var'}:
        raise ValueError("axis must be 'obs' or 'var'")
    if key not in df:
        raise ValueError(f"key '{key}' not found in adata.{axis}")

    val_min, val_max = df[key].quantile([q_min, q_max]).to_numpy()

    return _filter_axis(adata=adata, axis=axis, key=key, val_min=val_min, val_max=val_max, apply=apply)


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


def _filter_axis(
    adata: 'AnnData',
    axis: str,
    key: str,
    val_min: float | None = None,
    val_max: float | None = None,
    apply: bool = False,
) -> np.ndarray:
    if axis == 'obs':
        df = adata.obs
        values = df[key].to_numpy()
        mask = np.ones(adata.n_obs, dtype=bool)
    elif axis == 'var':
        df = adata.var
        values = df[key].to_numpy()
        mask = np.ones(adata.n_vars, dtype=bool)
    else:
        raise ValueError("axis must be 'obs' or 'var'")

    if key not in df:
        raise ValueError(f"key '{key}' not found in adata.{axis}")

    if val_min is not None:
        mask &= values >= val_min
    if val_max is not None:
        mask &= values <= val_max

    if apply:
        if axis == 'obs':
            return adata[mask].copy()
        else:
            return adata[:, mask].copy()
    return mask


def _summarize_filtering(n_total: int, results: dict) -> pl.DataFrame:
    series = [pl.Series(name, values) for name, values in results.items()]
    df = pl.DataFrame(series).with_row_index()

    df = df.group_by(results.keys()).agg(pl.len().alias('n_total')).sort('n_total', descending=True)
    df = df.with_columns(((pl.col('n_total') / n_total) * 100).alias('pct'))
    return df


# polars implementation
# def reorder_clusters_by_size(clust_umap, key: str, prefix='C') -> pl.DataFrame:
#     clust_sorted = clust_umap.group_by(key).agg(pl.len()).sort('len', descending=True)
#     size_order = clust_sorted[key].to_list()
#     numeric_order = np.arange(len(size_order)).astype(str).tolist()

#     numeric_order = [prefix + str(uid) for uid in numeric_order]
#     order_map = dict(zip(size_order, numeric_order))

#     new_umap = clust_umap.with_columns(pl.col(key).replace(order_map)).sort('seg_cell_id')
#     return new_umap


# def process_adata(adata: 'AnnData'):
#     ## 1. filter on cell size, remove top and bottom 1%
#     ## 2. remove genes not expressed by at least 5% of remaining cells
#     ## 3. filter on total transcripts and unique genes remove top 1% and bottom 5%

#     adata = filter_by_qtiles(adata, axis='obs', key='nuclei_area_um', q_min=0.01, q_max=0.99, apply=True)

#     min_cells = int(0.05 * adata.n_obs)
#     sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)

#     sc.pp.filter_cells(adata, min_counts=10)

#     adata = filter_by_qtiles(adata, axis='obs', key='total_counts', q_min=0.05, q_max=0.99, apply=True)
#     adata = filter_by_qtiles(adata, axis='obs', key='n_genes_by_counts', q_min=0.05, q_max=0.99, apply=True)

#     # normalize data
#     adata.layers['counts'] = adata.X
#     sc.pp.normalize_total(adata)

#     sc.pp.log1p(adata)
#     sc.pp.pca(adata)
#     sc.pp.neighbors(adata)

#     res = 1
#     sc.tl.leiden(
#         adata,
#         resolution=res,
#         objective_function='modularity',
#         n_iterations=-1,
#         random_state=42,
#         flavor='igraph',
#         key_added=f'leiden_{res:0.2f}',
#     )

#     sc.tl.umap(adata)

#     with warnings.catch_warnings():
#         warnings.simplefilter('ignore', PerformanceWarning)
#         sc.tl.rank_genes_groups(
#             adata,
#             groupby=f'leiden_{res:0.2f}',
#             use_raw=False,
#             method='wilcoxon',
#             pts=True,
#             key_added=f'leiden_{res:0.2f}_rank_genes_groups',
#         )
#     return adata


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
