import logging
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.colors as mcolors
import numpy as np
import polars as pl
from tqdm import tqdm

from .. import utils
from ..g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema
from ..schemas.validate import read_bin_file
from .workflow import workflow

if TYPE_CHECKING:
    from ..models import G4Xoutput


DEFAULT_COLOR = '#BFBFBF'

sg_palette = [
    '#72FFAB',
    '#A16CFD',
    '#FF7043',
    '#008FFF',
    '#D32F2F',
    '#7CB342',
    '#7F34BE',
    '#FFCA28',
    '#0C8668',
    '#FB4695',
    '#005EE1',
    '#28EDED',
    '#A17B64',
    '#FFFF58',
    '#BC29AE',
    '#006D8F',
    '#FFBAFF',
    '#FFD091',
    '#5C6BC0',
    '#F490B2',
]


@workflow
def edit_bin_file(
    g4x_obj: 'G4Xoutput',
    bin_file: Path,
    *,
    metadata: Path | None = None,
    cellid_key: str | None = None,
    cluster_key: str | None = None,
    clustercolor_key: str | None = None,
    emb_key: str | None = None,
    bin_out: Path | None = None,
    logger: logging.Logger,
) -> None:
    logger.info(f'Reading bin file: {bin_file}')

    bin_file = utils.validate_path(bin_file, must_exist=True, is_dir_ok=False, is_file_ok=True)

    if not bin_out:
        bin_out = bin_file
    else:
        bin_out = utils.validate_path(bin_out, must_exist=False, is_dir_ok=False, is_file_ok=True)

    cell_masks = read_bin_file(bin_file)

    update_clusters = update_emb = False

    if metadata is None:
        logger.info('No metadata provided. Building from adata source.')

        unused_keys = [k for k, v in locals().items() if k.endswith('_key') and v is not None]
        if len(unused_keys) > 0:
            logger.info(f'Ignoring provided keys: {unused_keys}. Only needed when metadata is provided.')

        logger.info('Loading adata.')
        obs_df = g4x_obj.load_adata(remove_nontargeting=False, load_clustering=True).obs
        obs_df = pl.from_pandas(obs_df)

        if 'cell_id' not in obs_df.columns:
            raise ValueError('Invalid adata. No "cell_id" column found in adata.obs.')

        cellid_key = 'cell_id'
        cluster_key = infer_cluster_key(obs_df)
        emb_key = infer_emb_key(obs_df)
        clustercolor_key = None

        obs_df = obs_df.cast({cluster_key: pl.UInt16})
        update_clusters = cluster_key is not None
        update_emb = emb_key is not None

    if metadata is not None:
        if emb_key is None and cluster_key is None:
            logger.info('neither emb_key nor cluster_key were provided, nothing to update.')
            return None

        obs_df = pl.read_csv(metadata)

        if cellid_key is not None:
            check_key(obs_df, 'cellid_key', cellid_key)
        else:
            logger.info('No cellid_key provided, attempting to infer from "label" column.')
            if 'label' in obs_df.columns:
                cellid_key = 'label'
            else:
                raise KeyError('cellid_key not provided and "label" not found in metadata. please provide cellid_key.')

        if cluster_key is not None:
            logger.info(f'Checking cluster_key: {cluster_key}')
            check_key(obs_df, 'cluster_key', cluster_key)

        if emb_key is not None:
            logger.info(f'Checking emb_key: {emb_key}')
            if f'{emb_key}_1' not in obs_df.columns and f'{emb_key}_2' not in obs_df.columns:
                raise ValueError(
                    f'Matching columns for emb_key: "{emb_key}" missing. \nRequired columns: {emb_key}_1 and {emb_key}_2'
                )

        if clustercolor_key is not None and cluster_key is None:
            raise ValueError('clustercolor_key provided without cluster_key. Please provide cluster_key as well.')
        if clustercolor_key is not None:
            logger.info(f'Checking clustercolor_key: {clustercolor_key}')
            check_key(obs_df, 'clustercolor_key', clustercolor_key)

        update_clusters = cluster_key is not None
        update_emb = emb_key is not None

    logger.info('Preparing metadata for updating.')
    obs_df = obs_df.with_row_index('meta_index')
    bin_cell_ids = [cell.cellId for cell in cell_masks.cellMasks]
    bin_df = pl.Series(name=cellid_key, values=bin_cell_ids).to_frame().with_row_index('bin_index')
    obs_df = obs_df.join(bin_df, on=cellid_key, how='right').sort('meta_index')

    not_present = obs_df.null_count()['meta_index'].item()
    if not_present > 0:
        logger.warning(
            f'{not_present} cells in bin file not found in metadata. These will have missing/default values.'
        )

    logger.info('Selecting final columns and filling missing values.')
    final_cols = {
        'cell_id': (cellid_key, 'unknown_cell'),
        'cluster_id': (cluster_key, '-1'),
        'cluster_color': (clustercolor_key, DEFAULT_COLOR),
        'umap_1': (f'{emb_key}_1', float('nan')),
        'umap_2': (f'{emb_key}_2', float('nan')),
    }

    for new_col, (old_col, fill) in final_cols.items():
        if old_col not in obs_df.columns:
            obs_df = obs_df.with_columns(pl.lit(fill).alias(new_col))
        else:
            if new_col != old_col:
                obs_df = obs_df.rename({old_col: new_col})
            obs_df = obs_df.with_columns(pl.col(new_col).fill_null(fill))

    obs_df = obs_df.select(final_cols.keys())

    if update_clusters:
        logger.info('Building cluster palette.')
        cluster_palette = parse_clustercolor_key(obs_df, clustercolor_key)

    ## Do the actual updating
    logger.info('Updating cells.')
    metadata = obs_df.to_pandas().set_index('cell_id')
    for cell in tqdm(cell_masks.cellMasks, desc='Updating cell data'):
        current_cellid = cell.cellId
        if current_cellid not in metadata.index:
            continue

        if update_clusters:
            cell.clusterId = str(metadata.loc[current_cellid, 'cluster_id'])

        if update_emb:
            cell.umapValues.umapX = metadata.loc[current_cellid, 'umap_1']
            cell.umapValues.umapY = metadata.loc[current_cellid, 'umap_2']

    if update_clusters:
        update_colormap(cell_masks, cluster_palette)

    logger.debug(f'Writing updated bin file --> {bin_out}')
    with open(bin_out, 'wb') as file:
        file.write(cell_masks.SerializeToString())


def check_key(obs_df: pl.DataFrame, key_name: str, key: str = None) -> None:
    if key not in obs_df.columns:
        raise KeyError(f'{key_name} "{key}" not found in metadata columns.')


def parse_clustercolor_key(obs_df: pl.DataFrame, clustercolor_key: str | None) -> dict:
    if clustercolor_key:
        color_mapping = obs_df.select('cluster_id', 'cluster_color').unique('cluster_id').sort('cluster_id')

        if not color_mapping['cluster_color'][0].startswith('#'):
            raise ValueError(f'Colors in column "{clustercolor_key}" must be in hex starting with #.')

        cluster_palette = {}
        for cluster, color in color_mapping.iter_rows():
            cluster_palette[cluster] = hex2rgb(color)

    else:
        cluster_palette = generate_cluster_palette(obs_df['cluster_id'])
    return cluster_palette


def hex2rgb(hex: str) -> list[int, int, int]:
    return [int(x * 255) for x in mcolors.to_rgb(hex)]


def infer_cluster_key(obs_df: pl.DataFrame) -> pl.DataFrame:
    cluster_keys = sorted([x for x in obs_df.columns if 'leiden' in x])
    if len(cluster_keys) == 0:
        print('Failed to infer cluster_key, no clustering information will be added.')
        cluster_key = None
    else:
        cluster_key = cluster_keys[0]
        print(f'Inferred cluster_key: {cluster_key}')

    return cluster_key


def infer_emb_key(obs_df: pl.DataFrame) -> pl.DataFrame:
    emb_keys = [c[:-2] for c in obs_df.columns if c.startswith('X_umap')]

    if len(emb_keys) == 0:
        print('No embedding keys available, UMAP coordinates will be missing.')
        emb_key = None
    else:
        emb_key = emb_keys[0]
        print(f'Inferred embedding_key: {emb_key}')
    return emb_key


def generate_cluster_palette(clusters: list, max_colors: int = 256) -> dict:
    """
    Generate a color palette mapping for cluster labels.

    This function assigns RGB colors to unique cluster labels using a matplotlib colormap.
    Clusters labeled as "-1" are assigned a default gray color `[191, 191, 191]`.

    The colormap used depends on the number of clusters:
        - `tab10` for ≤10 clusters
        - `tab20` for ≤20 clusters
        - `hsv` for more than 20 clusters, capped by `max_colors`

    Parameters
    ----------
    clusters : list
        A list of cluster identifiers (strings or integers). The special label '-1' is excluded
        from color mapping and handled separately.
    max_colors : int, optional
        Maximum number of colors to use in the HSV colormap. Only used if there are more than
        20 unique clusters. Default is 256.

    Returns
    -------
    dict
        A dictionary mapping each cluster ID (as a string) to a list of three integers
        representing an RGB color in the range [0, 255].

    Examples
    --------
    >>> generate_cluster_palette(['0', '1', '2', '-1'])
    {'0': [31, 119, 180], '1': [255, 127, 14], '2': [44, 160, 44], '-1': [191, 191, 191]}
    """

    unique_clusters = [c for c in np.unique(clusters) if c != '-1']
    n_clusters = len(unique_clusters)

    if n_clusters <= 20:
        hex_list = sg_palette

    else:
        from matplotlib.pyplot import get_cmap

        hex_list = get_cmap('hsv', min(max_colors, n_clusters)).colors

    cluster_palette = {}
    for i, cluster in enumerate(unique_clusters):
        cluster_palette[str(cluster)] = hex2rgb(hex_list[i])

    cluster_palette['-1'] = hex2rgb(DEFAULT_COLOR)

    return cluster_palette


def update_colormap(cell_masks: CellMasksSchema.CellMasks, cluster_palette: dict[str, list[int]]) -> None:
    cell_masks.ClearField('colormap')
    for cluster_id, color in cluster_palette.items():
        entry = CellMasksSchema.ColormapEntry()
        entry.clusterId = cluster_id
        entry.color.extend(color)
        cell_masks.colormap.append(entry)
