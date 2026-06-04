from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix

from ... import c

if TYPE_CHECKING:
    from anndata import AnnData


LOGGER = logging.getLogger(__name__)


# region main functions
def init_adata(
    manifest: pl.DataFrame,
    cell_metadata: pl.DataFrame,
    cell_x_gene: pl.DataFrame,
    cell_x_protein: pl.DataFrame | None = None,
    logger: logging.Logger | None = None,
) -> 'AnnData':

    log = logger or LOGGER

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
    _validate_cell_ids(cell_metadata, cell_x_gene, 'CellMetadata', 'CellxGene')

    if cell_x_protein is not None:
        cell_x_protein = cell_x_protein.sort(c.CELL_ID_NAME)
        _validate_cell_ids(cell_metadata, cell_x_protein, 'CellMetadata', 'CellxProt')
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

    adata = _sanitize_categorical_columns(adata, threshold=10)
    return adata.copy()


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
