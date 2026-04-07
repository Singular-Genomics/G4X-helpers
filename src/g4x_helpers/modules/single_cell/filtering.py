from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
import polars as pl
from anndata import AnnData

if TYPE_CHECKING:
    from anndata import AnnData


LOGGER = logging.getLogger(__name__)

DEFAULT_FILTER_CONFIG = [
    dict(filter_type='cells', alias='area', key='cell_area_um', quantile_min=0.01, quantile_max=0.995),
    dict(filter_type='cells', alias='intensity', key='nuclearstain_intensity_mean', subset='__notna__'),
    dict(filter_type='cells', alias='counts', key='total_counts', value_min=10, quantile_max=0.995),
    dict(filter_type='cells', alias='genes', key='n_genes_by_counts', value_min=5),
    dict(filter_type='cells', alias='controls', key='pct_counts_ctrl', value_max=5),
    dict(filter_type='genes', alias='targeting', key='probe_type', subset='targeting'),
    dict(filter_type='genes', alias='coverage', key='n_cells_by_counts', value_min=100),
]


def filter_adata(
    adata: 'AnnData',
    filter_panel: 'FilterPanel' | None = None,
    *,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    log.info('Filtering adata')

    if filter_panel is None:
        filter_panel = _get_default_filter_panel()

    obs_total, var_total = adata.n_obs, adata.n_vars
    cell_summary, gene_summary = filter_panel.filter(adata, apply=True)

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning('No cells or genes remaining after filtering. Returning empty AnnData object.')
        return adata, cell_summary, gene_summary

    for df, name, n_total in [(cell_summary, 'cells', obs_total), (gene_summary, 'genes', var_total)]:
        n_retained = df.filter(pl.all_horizontal(pl.col('^.*_ok$'))).select('n_total').item()
        retained = n_retained / n_total
        log.info('Retained {:,} ({:.2%}) {} after filtering'.format(n_retained, retained, name))

    return adata, cell_summary, gene_summary


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

    def _filter_axis(
        self,
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

        # if self.subset is not None:
        #     df = df[df[self.key] == self.subset]

        if self.subset is not None:
            values = df[self.key].to_numpy()

            if self.subset == '__notna__':
                mask = ~pd.isna(values)
            else:
                mask = values == self.subset
            if apply:
                return adata[mask].copy() if self.axis == 'obs' else adata[:, mask].copy()
            return mask

        v_min, v_max = self.resolve_thresholds(df)
        return self._filter_axis(adata, axis=self.axis, key=self.key, val_min=v_min, val_max=v_max, apply=apply)


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
        cell_summary = pl.DataFrame()
        gene_summary = pl.DataFrame()
        cell_results = {}
        gene_results = {}
        for mth in self.filters:
            if mth.filter_type == 'cells':
                cell_results[f'{mth.alias}_ok'] = mth.filter(adata, apply=False)
                all_passed_cell = np.logical_and.reduce(list(cell_results.values()))
            else:
                gene_results[f'{mth.alias}_ok'] = mth.filter(adata, apply=False)
                all_passed_gene = np.logical_and.reduce(list(gene_results.values()))

        if cell_results:
            cell_summary = self._summarize_filtering(adata.n_obs, cell_results)
        if gene_results:
            gene_summary = self._summarize_filtering(adata.n_vars, gene_results)

        if apply:
            adata._inplace_subset_obs(all_passed_cell)
            adata._inplace_subset_var(all_passed_gene)

        if return_masks:
            return (all_passed_cell, cell_summary), (all_passed_gene, gene_summary)

        return cell_summary, gene_summary


def _get_default_filter_panel():
    default_filter_panel = FilterPanel(filters=[FilterMethod(**cfg) for cfg in DEFAULT_FILTER_CONFIG])
    return default_filter_panel


def _format_summary(summary):
    df = summary.with_columns(pl.col('pct').round(3)).rename({'pct': '%'})
    bool_cols = df.select(pl.col(pl.Boolean)).columns

    for col in bool_cols:
        df = df.with_columns(pl.when(pl.col(col)).then(pl.lit('✓')).otherwise(pl.lit('-')).alias(col))

    lines = str(df).splitlines()
    filtered = '\n'.join([line for i, line in enumerate(lines) if i not in (0, 3, 4)])
    return filtered


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
