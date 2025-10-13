import logging
import multiprocessing
import random
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from skimage.morphology import disk, erosion
from tqdm import tqdm

from .. import utils
from . import CellMasksSchema_pb2 as CellMasksSchema
from . import bin_utils as butls


def bin_file_path(g4x_out, out_dir: str | Path | None = None) -> Path:
    file_name = f'{g4x_out.sample_id}.bin'
    
    if out_dir:
        out_dir = utils.validate_path(out_dir, must_exist=True, is_dir_ok=True, is_file_ok=False)
        outfile = Path(out_dir) / file_name
    else:
        outfile = utils.create_custom_out(out_dir, 'g4x_viewer', file_name)

    return outfile


def seg_converter(
    adata: ad.AnnData,
    seg_mask: np.ndarray,
    out_path: str | Path,
    *,
    metadata: str | Path | None = None,
    cluster_key: str | None = None,
    emb_key: str | None = None,
    protein_list: list[str] | None = None,
    n_threads: int = 4,
    logger: logging.Logger | None = None,
) -> None:
    
    logger.info('Making G4X-Viewer bin file.')
    
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
    border = butls.get_border(seg_mask)
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
    centroid_y = obs_df['cell_x'].tolist()
    centroid_x = obs_df['cell_y'].tolist()

    if 'area' in obs_df.columns:
        areas = obs_df['area'].tolist()
    else:
        areas = obs_df['nuclei_expanded_area'].tolist()
    total_counts = obs_df['total_counts'].tolist()
    total_genes = obs_df['n_genes_by_counts'].tolist()

    if clusters_available:
        clusters = obs_df[cluster_key].tolist()
        cluster_palette = butls.generate_cluster_palette(clusters)
        cluster_colors = obs_df[cluster_key].astype(str).map(cluster_palette).tolist()
    else:
        cmap = plt.get_cmap('hsv', 100)
        clusters = ['-1'] * num_cells
        cluster_colors = [[int(255 * x) for x in cmap(i % 100)[:3]] for i in range(num_cells)]
        random.shuffle(cluster_colors)

    if protein_list:
        prot_vals = list(obs_df[protein_list].to_dict(orient='index').values())
    else:
        prot_vals = [{} for _ in range(num_cells)]
    if emb_key:
        umap_x = obs_df[f'{emb_key}_1'].to_numpy()
        umap_y = obs_df[f'{emb_key}_2'].to_numpy()
    else:
        umap_x = np.zeros(num_cells)
        umap_y = np.zeros(num_cells)

    ## refine polygons
    logger.debug('Refining polygons.')

    pq_args = [
        (k, cx, cy, sorted_nonzero_values, sorted_rows, sorted_cols)
        for k, cx, cy in zip(np.arange(1, num_cells + 1), centroid_x, centroid_y)
    ]

    with multiprocessing.Pool(processes=n_threads) as pool:
        polygons = pool.starmap(butls.refine_polygon, pq_args)

    ## do conversion
    # logger.debug(f"{sample_id}: Converting data to protobuff format...")
    segmentation_source = {
        'xs': [[y[0] for y in x] for x in polygons],
        'ys': [[y[1] for y in x] for x in polygons],
        'colors': cluster_colors,
        'cell_id': cell_ids,
        'area': areas,
        'centroid_x': centroid_x,
        'centroid_y': centroid_y,
        'total_tx': total_counts,
        'total_genes': total_genes,
        'cluster_id': clusters,
        'umap_x': umap_x,
        'umap_y': umap_y,
        'protein': prot_vals,
    }

    outputCellSegmentation = CellMasksSchema.CellMasks()

    for index in tqdm(range(len(segmentation_source['cell_id'])), desc='Processing cells'):
        try:
            cellPolygonPoints = [
                coord
                for pair in zip(segmentation_source['ys'][index], segmentation_source['xs'][index])
                for coord in pair
            ]
            cellPolygonColor = segmentation_source['colors'][index]
            cellId = segmentation_source['cell_id'][index]
            cellTotalCounts = segmentation_source['total_tx'][index]
            cellTotalGenes = segmentation_source['total_genes'][index]
            cellArea = segmentation_source['area'][index]
            clusterId = segmentation_source['cluster_id'][index]
            cellProt = segmentation_source['protein'][index]
            cellUmapX = segmentation_source['umap_x'][index]
            cellUmapY = segmentation_source['umap_y'][index]
        except Exception as e:
            logger.debug(e)
            pass
        outputMaskData = outputCellSegmentation.cellMasks.add()
        outputMaskData.vertices.extend(cellPolygonPoints + cellPolygonPoints[:2])
        outputMaskData.color.extend(cellPolygonColor)
        outputMaskData.cellId = str(cellId)
        outputMaskData.area = str(cellArea)
        outputMaskData.totalCounts = str(cellTotalCounts)
        outputMaskData.totalGenes = str(cellTotalGenes)
        outputMaskData.clusterId = str(clusterId)
        outputMaskData.umapValues.umapX = cellUmapX
        outputMaskData.umapValues.umapY = cellUmapY
        outputMaskData.proteins.update(cellProt)

    if clusters_available:
        for cluster_id, color in cluster_palette.items():
            entry = CellMasksSchema.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            outputCellSegmentation.colormap.append(entry)
    else:
        entry = CellMasksSchema.ColormapEntry()
        entry.clusterId = '-1'
        entry.color.extend([int(31), int(119), int(180)])
        outputCellSegmentation.colormap.append(entry)

    ## write to file
    with open(out_path, 'wb') as file:
        file.write(outputCellSegmentation.SerializeToString())
    
    logger.debug(f'G4X-Viewer bin --> {out_path}')


def seg_updater(
    bin_file: str | Path,
    metadata_file: str | Path,
    out_path: Path,
    *,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    cluster_color_key: str | None = None,
    emb_key: str | None = None,
    logger: logging.Logger,
) -> None:
    ## pre-flight
    if emb_key is None and cluster_key is None:
        logger.warning('neither embedding nor cluster keys were provided, nothing to update.')
        return None

    ## load the bin file
    logger.info(f'Loading {bin_file}.')
    with open(bin_file, 'rb') as f:
        data = f.read()
    cell_masks = CellMasksSchema.CellMasks()
    cell_masks.ParseFromString(data)

    ## load the metadata
    if cellid_key is None:
        logger.info('cellid_key not provided, assuming cell IDs are in first column of metadata.')
        metadata = pd.read_csv(metadata_file, index_col=0, header=0)
    else:
        metadata = pd.read_csv(metadata_file, index_col=None, header=0)
        if cellid_key not in metadata.columns:
            raise KeyError(f'{cellid_key} not a valid column in metadata.')
        metadata.set_index(cellid_key, inplace=True)

    ## check for clustering
    if cluster_key is not None:
        if cluster_key not in metadata.columns:
            raise KeyError(f'{cluster_key} not a valid column in metadata.')
        update_cluster = True
        logger.debug('Updating cluster IDs.')
    else:
        update_cluster = False
        logger.debug('Not updating cluster IDs.')

    ## check for cluster colors
    if cluster_color_key is not None:
        if cluster_key is None:
            raise ValueError('cluster_color_key was provided, but cluster_key was not provided.')
        if cluster_color_key not in metadata.columns:
            raise KeyError(f'{cluster_color_key} not a valid column in metadata.')
        color = metadata[cluster_color_key].iat[0]
        assert color.startswith('#'), 'Cluster colors must be provided as hexcodes.'
        update_cluster_color = True
        logger.debug('Updating cluster colors.')
        cluster_palette = (
            metadata.drop_duplicates(subset=cluster_key)[[cluster_key, cluster_color_key]]
            .set_index(cluster_key)
            .to_dict()[cluster_color_key]
        )
        cluster_palette = {str(k): butls.hex2rgb(v) for k, v in cluster_palette.items()}
    else:
        if cluster_key is not None:
            update_cluster_color = True
            logger.debug('Auto-assigning colors to new clustering.')
            cluster_color_key = 'cluster_color'
            cluster_palette = butls.generate_cluster_palette(metadata[cluster_key].tolist())
            metadata['cluster_color'] = metadata[cluster_key].astype(str).map(cluster_palette).tolist()
        else:
            update_cluster_color = False
            logger.debug('Not updating cluster colors.')

    ## check for embedding
    if emb_key is not None:
        if f'{emb_key}_1' not in metadata.columns or f'{emb_key}_2' not in metadata.columns:
            raise KeyError(f'{emb_key}_1 and {emb_key}_2 are not valid columns in metadata.')
        update_emb = True
        logger.debug('Updating embedding.')
    else:
        update_emb = False
        logger.debug('Not updating embedding.')

    ## Do the actual updating
    logger.info('Updating cells.')
    for cell in tqdm(cell_masks.cellMasks, desc='Updating cell data'):
        current_cellid = cell.cellId
        if current_cellid in metadata.index:
            if update_cluster:
                cell.clusterId = str(metadata.loc[current_cellid, cluster_key])
            if update_cluster_color:
                # clear out the existing color entries:
                cell.ClearField('color')
                cell.color.extend(metadata.loc[current_cellid, cluster_color_key])
            if update_emb:
                cell.umapValues.umapX = metadata.loc[current_cellid, f'{emb_key}_1']
                cell.umapValues.umapY = metadata.loc[current_cellid, f'{emb_key}_2']
        else:
            logger.debug(f'{current_cellid} not found in metadata, not updating data for this cell.')

    if update_cluster_color:
        # clear the entire colormap list:
        cell_masks.ClearField('colormap')
        for cluster_id, color in cluster_palette.items():
            entry = CellMasksSchema.ColormapEntry()
            entry.clusterId = cluster_id
            entry.color.extend(color)
            cell_masks.colormap.append(entry)

    ## Write to file
    logger.debug(f'Writing updated bin file --> {out_path}')
    with open(out_path, 'wb') as file:
        file.write(cell_masks.SerializeToString())
