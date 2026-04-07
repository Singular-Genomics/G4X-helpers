import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import rankdata

from .filtering import FilterMethod

## RNA-protein pairs
PAIRS = {
    'CD3D': 'CD3',
    'CD4': 'CD4',
    'CD8A': 'CD8',
    'HLA-DRA': 'HLA-DR',
    'MS4A1': 'CD20',
    'MKI67': 'KI67',
    'CD274': 'PDL1',
    'PDCD1': 'PD1',
    'CTLA4': 'CTLA4',
    'CD68': 'CD68',
    'ITGAX': 'CD11c',
    'FOXP3': 'FOXP3',
    'PECAM1': 'CD31',
    'NCAM1': 'CD56',
    'VIM': 'VIM',
    'GZMB': 'GZMB',
    'TIGIT': 'TIGIT',
    'CD34': 'CD34',
    'PODXL': 'PODXL',
    'PVALB': 'PVALB',
    'UMOD': 'UMOD',
}

LOGGER = logging.getLogger(__name__)


def run_correlation_analysis(smp, logger: logging.Logger = LOGGER) -> None:
    adata = smp.load_adata(processed=False)
    fm = FilterMethod(filter_type='cells', key='nuclearstain_intensity_mean', subset='__notna__')
    adata = fm.filter(adata, apply=True)

    pr_corr_df = protein_protein(adata, logger)
    rna_corr_df, rna_pr_corr_df = protein_rna(adata, logger)
    return pr_corr_df, rna_pr_corr_df


def protein_protein(adata: AnnData, logger: logging.Logger = LOGGER) -> None:
    mask = _drop_zeros_mask(adata.obsm['protein'])
    prot_df = _get_protein_df(adata[mask])
    logger.info(f'Un-Filtered protein intensity Data Frame: {prot_df.shape=}')

    prot_df = _drop_zero_variance_proteins(prot_df, logger)
    logger.info(f'Filtered protein intensity Data Frame: {prot_df.shape=}')

    return _calculate_correlation(prot_df)


def protein_rna(adata: AnnData, logger: logging.Logger = LOGGER) -> None:
    mask_X = _drop_zeros_mask(adata.X)
    mask_protein = _drop_zeros_mask(adata.obsm['protein'])
    adata = adata[mask_X & mask_protein].copy()

    prot_df = _get_protein_df(adata)
    prot_df = _drop_zero_variance_proteins(prot_df, logger)

    filtered_pairs = {}
    for k, v in PAIRS.items():
        if k in adata.var_names and f'{v}_intensity_mean' in prot_df.columns:
            filtered_pairs[k] = v
    if len(filtered_pairs) == 0:
        logger.warning('No matching protein-RNA pairs in this data')
        return None

    ## get protein data for pairs
    prot_df = prot_df[[f'{x}_intensity_mean' for x in filtered_pairs.values()]].copy()
    logger.info(f'Filtered protein intensity Data Frame: {prot_df.shape=}')

    ## get RNA data
    rna_df = pd.DataFrame(
        adata[:, list(filtered_pairs.keys())].X.toarray(), index=adata.obs_names, columns=list(filtered_pairs.keys())
    )
    logger.info(f'Filtered RNA counts Data Frame: {rna_df.shape=}')

    rna_df = _drop_zero_count_genes(rna_df, logger)
    logger.info(f'Filtered RNA counts Data Frame after removing zero counts: {rna_df.shape=}')

    ## calculate correlation
    rna_corr_df = _calculate_correlation(rna_df)

    ## get ordering for later
    row_order = rna_df.columns.tolist()
    prot_row_order = [f'{filtered_pairs[y]}_intensity_mean' for y in row_order]

    ## merge RNA and Protein and specify order
    final_df = rna_df.merge(prot_df, how='left', left_index=True, right_index=True)
    final_df = final_df[row_order + prot_row_order].copy()

    ## calculate correlation
    ranked_final_df = rankdata(final_df.to_numpy(), axis=0)
    rna_pr_corr_df = _calculate_correlation(ranked_final_df, index=final_df.columns, columns=final_df.columns)

    return rna_corr_df, rna_pr_corr_df


# region helper functions
def _calculate_correlation(df, index=None, columns=None):
    index = index if index is not None else df.index.tolist()
    columns = columns if columns is not None else df.columns.tolist()
    corr = np.corrcoef(df, rowvar=False)
    corr_df = pd.DataFrame(corr, index=columns, columns=columns)
    return corr_df


def _drop_zeros_mask(arr):
    return np.asarray(arr.sum(axis=1)).ravel() > 0


def _get_protein_df(adata):
    if 'protein' not in adata.obsm:
        raise ValueError('Protein data not found in adata.obsm["protein"]')
    names_handle = [n + '_intensity_mean' for n in adata.uns['protein_names']]
    prot_df = pd.DataFrame(adata.obsm['protein'], columns=names_handle, index=adata.obs_names)
    return prot_df


def _drop_zero_variance_proteins(prot_df, logger: logging.Logger = LOGGER) -> None:
    # Constant columns = 0 variance = nan correlation
    col_var = np.var(prot_df, axis=0)
    zero_var_cols = col_var == 0
    zero_var_proteins = list(np.array(prot_df.columns)[zero_var_cols])
    if len(zero_var_proteins) > 0:
        logger.warning(f'Proteins had zero variance in sampled Data Frame: {zero_var_proteins}')

    return prot_df.loc[:, ~zero_var_cols]


def _drop_zero_count_genes(rna_df, logger: logging.Logger = LOGGER) -> None:
    non_zero_counts = rna_df.sum(axis=0) > 0
    zero_count_genes = list(np.array(rna_df.columns)[~non_zero_counts])

    if len(zero_count_genes) > 0:
        logger.warning(f'Genes total 0 counts in sampled cells: {zero_count_genes}')

    return rna_df.loc[:, non_zero_counts]
