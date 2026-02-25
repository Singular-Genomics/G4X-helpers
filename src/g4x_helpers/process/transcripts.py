from __future__ import annotations

from typing import TYPE_CHECKING

import polars as pl


if TYPE_CHECKING:
    from g4x_helpers.g4x_output import G4Xoutput


def create_cell_x_gene(g4x_obj: G4Xoutput, return_lazy: bool = False) -> tuple[pl.DataFrame, pl.DataFrame]:
    tx_table = g4x_obj.transcript_table_path
    reads = pl.scan_csv(tx_table)

    cell_by_gene = (
        reads.filter(pl.col('cell_id') != 0)
        .group_by('cell_id', 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index='cell_id', on_columns=g4x_obj.genes)
    )

    print('Adding missing cells with zero counts.')
    cell_by_gene = g4x_obj.cell_frame.lazy().join(cell_by_gene, left_on='seg_cell_id', right_on='cell_id', how='left')

    existing = cell_by_gene.collect_schema().names()
    if not set(existing[1:]) == set(g4x_obj.genes):
        raise ValueError('Mismatch between cell_by_gene columns and g4x_obj.genes')

    cell_by_gene = cell_by_gene.select(['seg_cell_id'] + g4x_obj.genes)
    cell_by_gene = cell_by_gene.fill_null(0)

    if return_lazy:
        return cell_by_gene
    return cell_by_gene.collect()
