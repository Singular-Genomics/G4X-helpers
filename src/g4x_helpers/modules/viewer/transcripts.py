import colorsys
import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import polars as pl
from numcodecs import Blosc

from ... import c
from ... import logging_utils as logut
from ...schema.definition import Dgex, Manifest, TxTable
from ..workflow import PRESET_SOURCE, collect_input
from .utils import create_array, populate_zarr_metadata

if TYPE_CHECKING:
    from zarr.hierarchy import Group as zGroup

    from ...g4x_output import G4Xoutput

LOGGER = logging.getLogger(__name__)


def write_transcripts(
    smp: 'G4Xoutput',
    root_group: 'zGroup',
    tx_table: str = PRESET_SOURCE,
    manifest: str = PRESET_SOURCE,
    dgex: str = PRESET_SOURCE,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
) -> None:

    log = LOGGER or logger
    log.info('Preparing transcript data')

    # 1: load inputs
    txtable_in = collect_input(smp, tx_table, validator=TxTable, logger=log)
    gene_metadata = get_gene_metadata(smp, manifest=manifest, dgex=dgex, logger=log)

    # 2: load tx table and filter to relevant columns
    aggregation_level = c.GENE_ID_NAME
    keep_cols = ['x_pixel_coordinate', 'y_pixel_coordinate', c.CELL_ID_NAME, aggregation_level]
    df = txtable_in.load(lazy=False).select(keep_cols)

    # 3: check that all genes in tx table have metadata
    gene_list = df[c.GENE_ID_NAME].unique().sort().to_list()

    unavailable_genes = [g for g in gene_list if g not in gene_metadata]
    if unavailable_genes:
        raise ValueError(f'The following genes are missing from the gene metadata: {unavailable_genes}')

    # 4: construct pyramid
    pyramid, tile_specs = build_tx_pyramid(
        image_resolution=smp.shape,
        total_points=df.height,
        target_points_per_tile=5000,
        min_tile_size=256,
    )

    msg = f'\n{tile_specs}\n'
    for level, specs in pyramid.items():
        msg += f'Level {level}: tile_size: {specs["tile_size"]} - scale: {specs["scale_fct"]}\n'

    logut.log_msg_wrapped('Tile specs:', msg, logger=log, level='debug')

    # 5: assign tiles to tx and construct tile dataframes
    pyramid = construct_tile_dfs(df, pyramid)

    # 6: populate attrs
    layer_config = {
        'layers': len(pyramid) - 1,
        'tile_size': pyramid[len(pyramid) - 1]['tile_size'],
        'layer_height': smp.shape[0],
        'layer_width': smp.shape[1],
        'coordinate_order': ['x_pixel_coordinate', 'y_pixel_coordinate'],
    }

    gene_colors = {k: v['color'] for k, v in gene_metadata.items()}
    populate_zarr_metadata(root_group, gene_colors=gene_colors, tx_layer_config=layer_config)

    log.info('Writing transcript data')
    tx_group = root_group['transcripts']

    write_tx_zarr(tx_group, pyramid, overwrite=overwrite)


def build_tx_pyramid(
    image_resolution: tuple[int, int], total_points: int, target_points_per_tile: int = 5000, min_tile_size: int = 256
) -> dict[int, dict[str, Any]]:

    tile_specs = choose_square_tiling(
        image_resolution_hw=image_resolution,
        total_points=total_points,
        target_points_per_tile=target_points_per_tile,
        min_tile_size=min_tile_size,
    )

    pyramid = {}
    level = 0
    tile_size = tile_specs['tile_size']

    # build pyramid until last tile size exceeds image resolution
    while tile_size <= max(image_resolution):
        scale = 2**level
        tile_size = tile_specs['tile_size'] * scale
        sampling_fct = 1 / (scale**2)
        pyramid[level] = {'tile_size': tile_size, 'scale_fct': scale, 'sampling_fct': sampling_fct}
        level += 1

    return pyramid, tile_specs


def choose_square_tiling(
    image_resolution_hw: tuple[int, int], total_points: int, target_points_per_tile: int, min_tile_size=64
):
    tiles_needed = np.ceil(total_points / target_points_per_tile).astype(int)

    H, W = image_resolution_hw  # (H, W)
    s = max(1, int(np.sqrt((W * H) / tiles_needed)))  # start guess

    def _num_tiles(s):
        nx = np.ceil(W / s)
        ny = np.ceil(H / s)
        return nx, ny

    while True:
        nx, ny = _num_tiles(s)
        if nx * ny >= tiles_needed:
            break
        s -= 1

    if s < min_tile_size:
        s = min_tile_size
        nx, ny = _num_tiles(s)

    res = {'tile_size': int(s), 'nx': int(nx), 'ny': int(ny), 'grid_total': int(nx * ny)}
    return res


def construct_tile_dfs(df: pl.DataFrame, pyramid: dict[int, dict[str, Any]]) -> dict[int, dict[str, Any]]:
    for level, specs in pyramid.items():
        tile_size, _, sampling_fct = specs.values()
        df_smp = (
            df.sample(fraction=sampling_fct, with_replacement=False)
            .with_columns(
                (pl.col('x_pixel_coordinate') / tile_size).cast(pl.Int32).alias('tile_x'),
                (pl.col('y_pixel_coordinate') / tile_size).cast(pl.Int32).alias('tile_y'),
            )
            .sort('tile_y', 'tile_x')
        )
        pyramid[level].update({'tile_df': df_smp})
    return pyramid


def write_tx_zarr(
    tx_group: 'zGroup',
    pyramid: dict[int, dict[str, Any]],
    overwrite: bool = False,
    logger: logging.Logger | None = None,
):
    log = LOGGER or logger
    compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)

    for level in pyramid:
        log.debug(f'Writing transcript data for level {level}')
        tile_df = pyramid[level]['tile_df']
        all_coords = tile_df.select(['x_pixel_coordinate', 'y_pixel_coordinate']).to_numpy()
        all_cell_ids = tile_df.select(c.CELL_ID_NAME).to_numpy()
        all_gene_names = tile_df.select(c.GENE_ID_NAME).to_numpy()
        all_tile_ids = tile_df.select(['tile_y', 'tile_x']).to_numpy()

        tile_ids = tile_df.unique(['tile_y', 'tile_x'], maintain_order=True).select(['tile_y', 'tile_x']).to_numpy()

        for tile in tile_ids:
            y = str(tile[0]).zfill(2)
            x = str(tile[1]).zfill(2)
            tile_path = f'p{level}/y{y}/x{x}'
            tile_group = tx_group.require_group(tile_path)

            idx = np.where((all_tile_ids == tile).all(axis=1))[0]

            coords = all_coords[idx].astype(np.int32)
            gene_names = all_gene_names[idx].astype('U10')
            cell_ids = all_cell_ids[idx].astype(np.int32)

            for key, arr in [('position', coords), ('gene_name', gene_names), ('cell_id', cell_ids)]:
                if overwrite and key in tile_group:
                    del tile_group[key]
                create_array(tile_group, key, data=arr, compressor=compressor)


def get_gene_metadata(smp, manifest, dgex, logger: logging.Logger | None = None):
    log = logger or LOGGER

    manifest_in = collect_input(smp, manifest, validator=Manifest, logger=log)
    dgex_in = collect_input(smp, dgex, validator=Dgex, validate=False, logger=log)

    tx_panel = manifest_in.load()

    if dgex_in.is_valid:
        dgex = dgex_in.load()

        leiden_result = (
            dgex.unique(['leiden_res', 'cluster_id'])
            .group_by(['leiden_res'])
            .agg(pl.len())
            .sort('len')
            .head(1)['leiden_res']
            .item()
        )

        dgex = dgex.filter(pl.col('leiden_res') == leiden_result)
        dgex_fil = dgex.filter(pl.col('score') > 0)

        assignments = assign_colors_to_clusters(dgex_fil)
        colors = complete_panel_colors(tx_panel=tx_panel, assignments=assignments)

        gene_metadata = {}
        for g in colors.iter_rows(named=True):
            gene_metadata[g['gene_id']] = {'color': hex_to_rgb(g['hex'])}
    else:
        gene_list = tx_panel['gene_name'].unique().sort().to_list()

        N = len(gene_list)  # number of colors
        colors = np.random.randint(0, 256, size=(N, 3), dtype=np.uint8)
        gene_metadata = {g: {'color': tuple(c.tolist())} for g, c in zip(gene_list, colors)}
        gene_metadata

    return gene_metadata


# region colors
def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')

    return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))


def hsv_to_hex(h, s, v):
    # wrap hue into [0,1]
    h = h % 1.0

    r, g, b = colorsys.hsv_to_rgb(h, s, v)

    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


def _normalize_range(df, column, out_range=(0, 1)):
    result = df.with_columns((out_range[0] + (out_range[1] - out_range[0]) * pl.col(column)).alias(column))
    return result


def rank_normalized(df, column, new_column, out_range=(0, 1)):
    result = df.with_columns(((pl.col(column).rank(method='average') - 1) / (pl.len() - 1)).alias(new_column))
    result = _normalize_range(result, new_column, out_range)
    return result


def column_normalized(df, column, new_column, out_range=(0, 1)):
    result = df.with_columns(
        pl.when(pl.col(column).max() == pl.col(column).min())
        .then(0.0)
        .otherwise((pl.col(column) - pl.col(column).min()) / (pl.col(column).max() - pl.col(column).min()))
        .alias(new_column)
    )
    result = _normalize_range(result, new_column, out_range)
    return result


def assign_colors(df, hue=180, hue_spread=0.1, sat_range=(0, 1), val_range=(0, 1)):

    hue_n = hue / 360
    hue_range = (hue_n - hue_spread / 2, hue_n + hue_spread / 2)

    subset = df.sort('gene_id').with_row_index(name='rank')

    subset = rank_normalized(subset, 'rank', 'hue', out_range=hue_range)
    subset = column_normalized(subset, 'score', 'sat', out_range=sat_range)
    subset = rank_normalized(subset, 'pct_nz_group', 'val', out_range=val_range)

    # subset = subset.with_columns(
    #     pl.struct(['hue', 'sat', 'val']).map_elements(lambda x: hsv_to_hex(x['hue'], x['sat'], x['val'])).alias('hex')
    # )

    return subset


def assign_colors_to_clusters(df_fil):
    clusters = df_fil['cluster_id'].unique().sort().to_list()
    duplications = df_fil.group_by('gene_id').agg(pl.len()).sort('len', descending=True)

    singles = duplications.filter(pl.col('len') == 1)
    multis = duplications.filter(pl.col('len') > 1)

    assignments = df_fil.filter(pl.col('gene_id').is_in(singles['gene_id'].implode()))

    for gene in multis['gene_id'].to_list():
        match = df_fil.filter(pl.col('gene_id') == gene).sort('score', descending=True).head(1)
        assignments = assignments.vstack(match)

    step = 360 / len(clusters)
    spread = 1 / len(clusters) / 1.75

    final_subset = pl.DataFrame()
    for i, cluster in enumerate(clusters):
        hue = i * step
        subset = assignments.filter(pl.col('cluster_id') == cluster)
        subset = assign_colors(subset, hue=hue, hue_spread=spread, sat_range=(0.25, 1.0), val_range=(0.25, 1.0))
        final_subset = final_subset.vstack(subset)
    return final_subset


def complete_panel_colors(tx_panel, assignments):
    complete_panel = (
        tx_panel.unique(['probe_type', 'gene_name'])
        .sort('probe_type', descending=True)
        .select(['gene_name', 'probe_type'])
    )

    missing = complete_panel.filter(~pl.col('gene_name').is_in(assignments['gene_id'].implode()))

    missing = missing.with_columns(pl.lit(0.0).alias('hue'), pl.lit(0.0).alias('sat'))
    missing = missing.with_row_index(name='rank')

    missing = rank_normalized(missing, 'rank', 'val', out_range=(0.0, 1.0))

    colors = assignments.select('gene_id', 'hue', 'sat', 'val')
    missing_colors = missing.select('gene_name', 'hue', 'sat', 'val').rename({'gene_name': 'gene_id'})

    colors = colors.vstack(missing_colors)
    colors = colors.with_columns(
        pl.struct(['hue', 'sat', 'val']).map_elements(lambda x: hsv_to_hex(x['hue'], x['sat'], x['val'])).alias('hex')
    )
    return colors
