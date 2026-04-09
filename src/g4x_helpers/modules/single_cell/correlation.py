import logging

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import rankdata

from ... import constants as c
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


def run_correlation_analysis(
    adata: AnnData, downsample: int = 25_000, logger: logging.Logger = LOGGER
) -> tuple[pd.DataFrame, pd.DataFrame]:

    n_obs_in = adata.n_obs
    fm = FilterMethod(filter_type='cells', key=c.NUC_STAIN_INTENSITY, subset='__notna__')
    adata = fm.filter(adata, apply=True)

    mask_X = _drop_zeros_mask(adata.X)
    mask_protein = _drop_zeros_mask(adata.obsm['protein'])
    adata = adata[mask_X & mask_protein].copy()

    if downsample is not None:
        adata = downsample_adata(adata, downsample=downsample, logger=logger)

    if adata.n_obs < n_obs_in:
        logger.warning(
            f'Filtered out {n_obs_in - adata.n_obs} ({(n_obs_in - adata.n_obs) / n_obs_in:.2%}) cells with no RNA or Protein data'
        )

    pr_corr_df = protein_protein(adata, logger)
    _, rna_pr_corr_df = protein_rna(adata, logger)
    return pr_corr_df, rna_pr_corr_df


def protein_protein(adata: AnnData, logger: logging.Logger = LOGGER) -> None:
    prot_df = _get_protein_df(adata)
    prot_df = _drop_zero_variance_proteins(prot_df, logger)

    if not check_correlation_feasibility(prot_df, 'proteins', logger):
        return create_dummy_output()

    try:
        return _calculate_correlation(prot_df)
    except Exception as e:
        logger.error(f'Error calculating protein-protein correlation: {e}')
        return create_dummy_output()


def protein_rna(adata: AnnData, logger: logging.Logger = LOGGER) -> None:
    prot_df = _get_protein_df(adata)
    prot_df = _drop_zero_variance_proteins(prot_df, logger)

    if not check_correlation_feasibility(prot_df, 'proteins', logger):
        return create_dummy_output(), create_dummy_output()

    filtered_pairs = {}
    for k, v in PAIRS.items():
        if k in adata.var_names and f'{v}{c.IMG_INTENSITY_HANDLE}' in prot_df.columns:
            filtered_pairs[k] = v
    if len(filtered_pairs) == 0:
        logger.warning('No matching protein-RNA pairs in this data. Returning empty correlation matrices.')
        return create_dummy_output(), create_dummy_output()

    ## get protein data for pairs
    prot_df = prot_df[[f'{x}{c.IMG_INTENSITY_HANDLE}' for x in filtered_pairs.values()]].copy()
    if not check_correlation_feasibility(prot_df, 'proteins', logger):
        return create_dummy_output(), create_dummy_output()

    ## get RNA data
    rna_df = pd.DataFrame(
        adata[:, list(filtered_pairs.keys())].X.toarray(), index=adata.obs_names, columns=list(filtered_pairs.keys())
    )
    rna_df = _drop_zero_count_genes(rna_df, logger)

    if not check_correlation_feasibility(rna_df, 'genes', logger):
        return create_dummy_output(), create_dummy_output()

    try:
        ## calculate correlation
        rna_corr_df = _calculate_correlation(rna_df)

        ## get ordering for later
        row_order = rna_df.columns.tolist()
        prot_row_order = [f'{filtered_pairs[y]}{c.IMG_INTENSITY_HANDLE}' for y in row_order]

        ## merge RNA and Protein and specify order
        final_df = rna_df.merge(prot_df, how='left', left_index=True, right_index=True)
        final_df = final_df[row_order + prot_row_order].copy()

        ## calculate correlation
        ranked_final_df = rankdata(final_df.to_numpy(), axis=0)
        rna_pr_corr_df = _calculate_correlation(ranked_final_df, index=final_df.columns, columns=final_df.columns)
    except Exception as e:
        logger.error(f'Error calculating protein-RNA correlation: {e}')
        return create_dummy_output(), create_dummy_output()

    return rna_corr_df, rna_pr_corr_df


# region helper functions
def create_dummy_output():
    df = pd.DataFrame(
        [
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        columns=['not', 'available'],
        index=['not', 'available'],
    )
    return df


def check_correlation_feasibility(df, data_type: str, logger: logging.Logger = LOGGER) -> None:
    if df.shape[0] < 100 or df.shape[1] < 20:
        logger.warning(
            f'Only {df.shape[0]} cells and {df.shape[1]} {data_type} left after filtering. Returning empty correlation matrix.'
        )
        return False
    return True


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
    names_handle = [n + c.IMG_INTENSITY_HANDLE for n in adata.uns['protein_names']]
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


def downsample_adata(adata: AnnData, downsample: int, logger: logging.Logger = LOGGER):
    if adata.n_obs > downsample:
        idx = np.random.choice(adata.n_obs, size=downsample, replace=False)
        logger.debug(f'Protein correlation data subsampled to {downsample} cells.')
        return adata[idx, :].copy()
    else:
        logger.debug(f'Protein correlation data has {adata.shape[0]} cells, no subsampling applied.')
        return adata
