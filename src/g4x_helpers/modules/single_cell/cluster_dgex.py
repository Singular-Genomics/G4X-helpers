from __future__ import annotations

import logging
import traceback
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from anndata import AnnData

from ... import c, io

if TYPE_CHECKING:
    from anndata import AnnData


LOGGER = logging.getLogger(__name__)


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
        log.info('Running dgex on: %s', leiden)
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
