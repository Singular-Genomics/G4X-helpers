import warnings
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from pandas.errors import PerformanceWarning

from .. import constants as c
from .. import io

if TYPE_CHECKING:
    from anndata import AnnData

    from ..g4x_output import G4Xoutput


# TODO addition of protein values is missing
def create_adata(
    g4x_obj: 'G4Xoutput', cell_x_gene: pl.DataFrame, cell_metadata: pl.DataFrame | None = None
) -> 'AnnData':
    from anndata import AnnData
    from scipy.sparse import csr_matrix

    X = cell_x_gene.drop(c.CELL_ID_NAME).to_numpy().astype(np.uint16)
    X = csr_matrix(X)

    # if cell_metadata is None:
    #     cell_metadata = create_cell_metadata(g4x_obj, mask=None)

    obs_df = cell_metadata.to_pandas().set_index(c.CELL_ID_NAME)
    obs_df.index = obs_df.index.astype(str)

    gene_ids = pl.Series(name=c.GENE_ID_NAME, values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    panel_type = (
        io.parse_input_manifest(g4x_obj.tree.TranscriptPanel.path)
        .unique('gene_name')
        .select('gene_name', 'probe_type')
        .rename({'gene_name': c.GENE_ID_NAME})
    )

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

    adata = AnnData(X=X, obs=obs_df, var=var_df)

    adata.var['ctrl'] = adata.var['probe_type'] != 'targeting'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['ctrl'], inplace=True, percent_top=None)
    adata.var.drop(columns=['ctrl'], inplace=True)

    return adata.copy()


def filter_cells(adata: 'AnnData') -> pl.DataFrame:
    n_obs = adata.n_obs
    size_passed = filter_by_qtiles(adata, key='nuclei_area_um', axis='obs', q_min=0.01, q_max=0.995)

    max_counts = filter_by_qtiles(adata, axis='obs', key='total_counts', q_min=0, q_max=0.995)
    min_counts = filter_by_limits(adata, axis='obs', key='total_counts', val_min=1, val_max=None)
    counts_passed = max_counts & min_counts

    genes_passed = filter_by_limits(adata, axis='obs', key='n_genes_by_counts', val_min=5, val_max=None)

    controls_passed = filter_by_limits(adata, axis='obs', key='pct_counts_ctrl', val_min=0, val_max=5)

    all_passed = counts_passed & genes_passed & controls_passed & size_passed
    adata._inplace_subset_obs(all_passed)

    results = {
        'counts_ok': counts_passed,
        'genes_ok': genes_passed,
        'controls_ok': controls_passed,
        'size_ok': size_passed,
    }
    cell_results = _summarize_filtering(n_obs, results)

    n_retained = sum(all_passed)
    retained = n_retained / n_obs
    print(f'Retained {n_retained:,} ({retained:.2%}) cells after filtering')

    return cell_results


def filter_genes(adata: 'AnnData') -> pl.DataFrame:
    n_var = adata.n_vars
    targeting_passed = adata.var['probe_type'] == 'targeting'
    coverage_passed = filter_by_limits(adata, axis='var', key='n_cells_by_counts', val_min=100, val_max=None)
    genes_all_passed = targeting_passed & coverage_passed
    adata._inplace_subset_var(genes_all_passed)

    results = {'targeting_ok': targeting_passed, 'coverage_ok': coverage_passed}
    gene_results = _summarize_filtering(n_var, results)

    n_retained = sum(genes_all_passed)
    retained = n_retained / n_var
    print(f'Retained {n_retained:,} ({retained:.2%}) genes after filtering')

    return gene_results


def filter_by_limits(
    adata: 'AnnData',
    key: str,
    axis: str = 'obs',
    val_min: float | None = None,
    val_max: float | None = None,
    apply: bool = False,
) -> np.ndarray:
    return _filter_axis(adata=adata, axis=axis, key=key, val_min=val_min, val_max=val_max, apply=apply)


def filter_by_qtiles(
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


def optimize_leiden_clusters(
    adata: 'AnnData',
    cluster_name='leiden_clusters',
    cluster_prefix='C',
    ideal_clusters=6,
    init_res=0.5,
    max_attempts=15,
    backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
) -> None:
    step_sizes = {
        50: lambda x: x - (x * 0.5),
        10: lambda x: x - (x * 0.65),
        5: lambda x: x - (x * 0.8),
        2: lambda x: x - (x * 0.9),
    }

    results = {}

    i = 0
    advance = True
    res = init_res
    while advance:
        # TODO proper backend integration
        backend = io.get_backend(preference=backend)
        rsc = backend.rsc

        rsc.tl.leiden(adata, resolution=res, key_added='tmp_leiden', random_state=777)
        leiden_clusters = adata.obs['tmp_leiden']
        cats = leiden_clusters.cat.categories
        n_cluster = len(cats)
        d_ideal = n_cluster - ideal_clusters
        print(f'resolution: {round(res, 3)}, n-clusters: {n_cluster}, d-ideal: {d_ideal}')

        results[i] = {
            'resolution': res,
            'n_clusters': n_cluster,
            'delta_ideal': d_ideal,
            'assignment': leiden_clusters,
        }

        if n_cluster == ideal_clusters:
            advance = False

        x = abs(d_ideal)
        closest = min(step_sizes.keys(), key=lambda k: abs(k - x))
        d_res = step_sizes[closest](res)
        res = res - d_res if d_ideal > 0 else res + d_res
        # print(f'Updated resolution: {res}')
        i += 1
        if i >= max_attempts:
            advance = False

    adata.obs.drop(columns=['tmp_leiden'], inplace=True)

    # collect best leiden result
    best_key = min(
        results,
        key=lambda k: (
            abs(results[k]['delta_ideal']),
            results[k]['n_clusters'],
            results[k]['resolution'],
        ),
    )
    leiden_result = results[best_key]

    adata.obs[cluster_name] = leiden_result['assignment']
    adata.obs = _reorder_clusters_by_size(adata.obs, key=cluster_name, prefix=cluster_prefix)


# region utilities
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


# region processing
def process_adata(adata: 'AnnData'):
    ## 1. filter on cell size, remove top and bottom 1%
    ## 2. remove genes not expressed by at least 5% of remaining cells
    ## 3. filter on total transcripts and unique genes remove top 1% and bottom 5%

    adata = filter_by_qtiles(adata, axis='obs', key='nuclei_area_um', q_min=0.01, q_max=0.99, apply=True)

    min_cells = int(0.05 * adata.n_obs)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)

    sc.pp.filter_cells(adata, min_counts=10)

    adata = filter_by_qtiles(adata, axis='obs', key='total_counts', q_min=0.05, q_max=0.99, apply=True)
    adata = filter_by_qtiles(adata, axis='obs', key='n_genes_by_counts', q_min=0.05, q_max=0.99, apply=True)

    # normalize data
    adata.layers['counts'] = adata.X
    sc.pp.normalize_total(adata)

    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)

    res = 1
    sc.tl.leiden(
        adata,
        resolution=res,
        objective_function='modularity',
        n_iterations=-1,
        random_state=42,
        flavor='igraph',
        key_added=f'leiden_{res:0.2f}',
    )

    sc.tl.umap(adata)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PerformanceWarning)
        sc.tl.rank_genes_groups(
            adata,
            groupby=f'leiden_{res:0.2f}',
            use_raw=False,
            method='wilcoxon',
            pts=True,
            key_added=f'leiden_{res:0.2f}_rank_genes_groups',
        )
    return adata


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
