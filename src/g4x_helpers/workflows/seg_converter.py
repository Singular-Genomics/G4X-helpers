import logging
import multiprocessing
import warnings
from operator import itemgetter
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse
from skimage.morphology import disk, erosion

from ..modules import bin_generation as bg
from ..modules.g4x_viewer.v2 import CellMasksSchema_pb2 as CellMasksSchema_v2
from ..modules.g4x_viewer.v3 import CellMasksSchema_pb2 as CellMasksSchema_v3
from .decorator import workflow

DEFAULT_COLOR = [int(191), int(191), int(191)]


@workflow
def seg_converter(
    adata: ad.AnnData,
    seg_mask: np.ndarray,
    out_path: str | Path,
    *,
    bin_version: Literal['v2', 'v3'] = 'v3',
    metadata: str | Path | None = None,
    cluster_key: str | None = None,
    emb_key: str | None = None,
    protein_list: list[str] | None = None,
    n_threads: int = 4,
    logger: logging.Logger | None = None,
) -> None:
    logger.info('Creating G4X-viewer bin file.')

    warnings.filterwarnings(
        'ignore',
        message='FNV hashing is not implemented in Numba',
        category=UserWarning,
        module='numba.cpython.old_hashing',
    )

    if metadata is None:
        if cluster_key in adata.obs.columns:
            clusters_available = True
        else:
            clusters_available = False
        obs_df = adata.obs.copy()
    else:
        clusters_available = True
        clustered_df = pd.read_csv(metadata, index_col=0, header=0)
        if clustered_df.shape[1] > 1:
            assert cluster_key is not None, (
                'ERROR: multiple columns detected in cluster_info, cluster_key must be provided.'
            )
        else:
            cluster_key = clustered_df.columns[0]
        orig_df = adata.obs.copy()

        ## these are cells that were filtered out during clustering
        orig_df = orig_df.loc[list(set(orig_df.index) - set(clustered_df.index)), :].copy()
        for col in list(set(clustered_df.columns) - set(orig_df.columns)):
            orig_df[col] = '-1'

        obs_df = pd.concat([clustered_df, orig_df])
        obs_df.sort_index(inplace=True)

    ## initialize segmentation data
    ## we create polygons to define the boundaries of each cell mask
    logger.debug('Making polygons.')
    border = bg.get_border(seg_mask)
    seg_mask[border > 0] = 0
    eroded_mask = erosion(seg_mask, disk(1))
    outlines = seg_mask - eroded_mask
    sparse_matrix = csr_matrix(outlines)
    del seg_mask, border, eroded_mask, outlines

    nonzero_values = sparse_matrix.data
    nonzero_row_indices, nonzero_col_indices = sparse_matrix.nonzero()
    sorted_indices = np.argsort(nonzero_values)
    sorted_nonzero_values = nonzero_values[sorted_indices]
    sorted_rows = nonzero_row_indices[sorted_indices]
    sorted_cols = nonzero_col_indices[sorted_indices]

    ## add single-cell info
    logger.debug('Adding single-cell metadata.')
    cell_ids = obs_df.index.tolist()
    num_cells = len(cell_ids)
    adata = adata[cell_ids]
    gene_names = adata.var_names
    gex = adata.X
    if not issparse(gex):
        logger.debug('Converting dense expression matrix to CSR format.')
        gex = csr_matrix(gex)
    else:
        gex = gex.tocsr()
    centroid_y = obs_df['cell_x'].tolist()
    centroid_x = obs_df['cell_y'].tolist()

    obs_df['total_counts'] = obs_df['total_counts'].astype(int)
    obs_df['n_genes_by_counts'] = obs_df['n_genes_by_counts'].astype(int)

    area_cols = obs_df.filter(like='area').columns
    if not len(area_cols):
        raise ValueError('No area column found')

    if len(area_cols) == 1:
        area_select = area_cols[0]
    else:
        logger.debug(f'Multiple area columns found: {area_cols}. Looking for expanded area column.')
        expanded_area_cols = [f for f in area_cols if 'expanded' in f]
        if len(expanded_area_cols) != 1:
            raise ValueError(
                'No expanded area column found among multiple area columns.'
                if not expanded_area_cols
                else f'Multiple expanded area columns found: {expanded_area_cols}'
            )
        area_select = expanded_area_cols[0]

    logger.debug(f'Using area column: {area_select}')

    obs_df['area'] = obs_df[area_select].astype(int)

    if emb_key:
        obs_df['umap_0'] = obs_df[f'{emb_key}_1']
        obs_df['umap_1'] = obs_df[f'{emb_key}_2']
    else:
        obs_df['umap_0'] = 0
        obs_df['umap_1'] = 0

    if protein_list:
        protein_col_pos = []
        for prot in protein_list:
            protein_col_pos.append(int(np.argwhere([x == prot for x in obs_df.columns])) + 1)
            obs_df[prot] = obs_df[prot].fillna(0).astype(int)
        get_prot_vals = itemgetter(*protein_col_pos)

    ## refine polygons
    logger.debug('Refining polygons.')

    pq_args = [
        (k, cx, cy, sorted_nonzero_values, sorted_rows, sorted_cols)
        for k, cx, cy in zip(np.arange(1, num_cells + 1), centroid_x, centroid_y)
    ]

    with multiprocessing.Pool(processes=n_threads) as pool:
        polygons = pool.starmap(bg.refine_polygon, pq_args)

    ## do conversion
    if bin_version == 'v2':
        outputCellSegmentation = CellMasksSchema_v2.CellMasks()
    else:
        outputCellSegmentation = CellMasksSchema_v3.CellMasks()
        outputCellSegmentation.metadata.geneNames.extend(gene_names)
        if protein_list:
            protein_names = [x.split('_intensity')[0] for x in protein_list]
            outputCellSegmentation.metadata.proteinNames.extend(protein_names)

    ## add the cluster color map
    if clusters_available:
        obs_df['cluster_id'] = obs_df[cluster_key].astype(str)
        cluster_palette = bg.generate_cluster_palette(obs_df['cluster_id'].tolist())
        for cluster_id, color in cluster_palette.items():
            if bin_version == 'v2':
                entry = CellMasksSchema_v2.ColormapEntry()
            else:
                entry = CellMasksSchema_v3.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            outputCellSegmentation.colormap.append(entry)
    else:
        if bin_version == 'v2':
            entry = CellMasksSchema_v2.ColormapEntry()
        else:
            entry = CellMasksSchema_v3.ColormapEntry()
        entry.clusterId = '-1'
        entry.color.extend(DEFAULT_COLOR)
        outputCellSegmentation.colormap.append(entry)
        obs_df['cluster_id'] = '-1'

    ## add the individual cells
    for i, row in enumerate(obs_df.itertuples(index=True)):
        outputMaskData = outputCellSegmentation.cellMasks.add()
        cellPolygonPoints = [sub_coord for coord in polygons[i] for sub_coord in coord[::-1]]
        _ = outputMaskData.vertices.extend(cellPolygonPoints + cellPolygonPoints[:2])
        area = int(row.area)
        total_counts = int(row.total_counts)
        total_genes = int(row.n_genes_by_counts)
        cluster_id = str(getattr(row, 'cluster_id', '-1'))
        outputMaskData.cellId = str(row[0])
        if bin_version == 'v2':
            outputMaskData.area = str(area)
            outputMaskData.totalCounts = str(total_counts)
            outputMaskData.totalGenes = str(total_genes)
        else:
            outputMaskData.area = area
            outputMaskData.totalCounts = total_counts
            outputMaskData.totalGenes = total_genes
        outputMaskData.clusterId = cluster_id
        outputMaskData.umapValues.umapX = row.umap_0
        outputMaskData.umapValues.umapY = row.umap_1

        if bin_version == 'v3':
            start = gex.indptr[i]
            end = gex.indptr[i + 1]
            indices = gex.indices[start:end]
            values = gex.data[start:end]
            _ = outputMaskData.nonzeroGeneIndices.extend(indices.tolist())
            _ = outputMaskData.nonzeroGeneValues.extend(values.astype(int).tolist())

        if protein_list:
            if bin_version == 'v2':
                prot_vals = {prot: val for prot, val in zip(protein_list, get_prot_vals(row))}
                outputMaskData.proteins.update(prot_vals)
            else:
                _ = outputMaskData.proteinValues.extend(get_prot_vals(row))

    ## write to file
    with open(out_path, 'wb') as file:
        file.write(outputCellSegmentation.SerializeToString())

    logger.debug(f'G4X-viewer bin --> {out_path}')
