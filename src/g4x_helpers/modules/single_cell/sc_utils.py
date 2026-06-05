from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
from anndata import AnnData

from ... import constants as c

if TYPE_CHECKING:
    from anndata import AnnData


LOGGER = logging.getLogger(__name__)


def _sanitize_categorical_columns(adata, threshold: int = 10):
    for df in [adata.obs, adata.var]:
        for col_name in df.columns:
            col = df[col_name]
            if col.nunique() < threshold and col.dtype == 'O':
                df[col_name] = col.astype('category')
    return adata


def _validate_cell_ids(df1, df2, name_1, name_2):
    cxg_ok = df1[c.CELL_ID_NAME].equals(df2[c.CELL_ID_NAME])
    if not cxg_ok:
        raise ValueError(f'Cell IDs in {name_1} do not match those in {name_2}')


def _extract_umap_clustering(adata, cluster_keys):
    umap_df = pl.from_numpy(adata.obsm['X_umap'], schema={'UMAP1': pl.Float32, 'UMAP2': pl.Float32})
    leiden_df = pl.from_pandas(adata.obs[cluster_keys], include_index=True)
    clustering_umap = leiden_df.hstack(umap_df).sort(c.CELL_ID_NAME)
    return clustering_umap


def _get_protein_df(adata):
    if 'protein' not in adata.obsm:
        raise ValueError('Protein data not found in adata.obsm["protein"]')
    names_handle = [n + c.IMG_INTENSITY_HANDLE for n in adata.uns['protein_names']]
    prot_df = pd.DataFrame(adata.obsm['protein'], columns=names_handle, index=adata.obs_names)
    return prot_df


def downsample_adata(adata: AnnData, downsample: int, logger: logging.Logger = LOGGER):
    if adata.n_obs > downsample:
        idx = np.random.choice(adata.n_obs, size=downsample, replace=False)
        logger.debug(f'Protein correlation data subsampled to {downsample} cells.')
        return adata[idx, :].copy()
    else:
        logger.debug(f'Protein correlation data has {adata.shape[0]} cells, no subsampling applied.')
        return adata


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


# polars implementation
# def reorder_clusters_by_size(clust_umap, key: str, prefix='C') -> pl.DataFrame:
#     clust_sorted = clust_umap.group_by(key).agg(pl.len()).sort('len', descending=True)
#     size_order = clust_sorted[key].to_list()
#     numeric_order = np.arange(len(size_order)).astype(str).tolist()

#     numeric_order = [prefix + str(uid) for uid in numeric_order]
#     order_map = dict(zip(size_order, numeric_order))

#     new_umap = clust_umap.with_columns(pl.col(key).replace(order_map)).sort('seg_cell_id')
#     return new_umap
