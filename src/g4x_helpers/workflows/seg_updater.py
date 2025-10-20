import logging
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from ..modules import bin_generation as bg
from ..modules.g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema
from .decorator import workflow


@workflow
def seg_updater(
    bin_file: Path,
    metadata_file: Path,
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
        cluster_palette = {str(k): bg.hex2rgb(v) for k, v in cluster_palette.items()}
    else:
        if cluster_key is not None:
            update_cluster_color = True
            logger.debug('Auto-assigning colors to new clustering.')
            cluster_color_key = 'cluster_color'
            cluster_palette = bg.generate_cluster_palette(metadata[cluster_key].tolist())
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
