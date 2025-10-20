from __future__ import annotations

import logging
import math
import multiprocessing as mp
import os
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import polars as pl

from ..modules import tx_generation as tg
from .decorator import workflow

if TYPE_CHECKING:
    from ..models import G4Xoutput

mp.set_start_method('spawn', force=True)


@workflow
def tx_converter(
    g4x_obj: 'G4Xoutput',
    out_path: Path,
    *,
    aggregation_level: str = 'gene',
    n_threads: int = 4,
    sampling_fraction: float = 0.2,
    logger: logging.Logger,
) -> None:
    logger.info('Generating viewer transcript file.')

    ## prelims
    IMAGE_RESOLUTION = g4x_obj.shape
    if g4x_obj.platform == 'g4x-2lane':
        MIN_TILE_SIZE = 1028
    else:
        MIN_TILE_SIZE = 512
    out_dir = g4x_obj.data_dir / 'g4x_viewer_temp'
    os.makedirs(out_dir, exist_ok=True)

    ## get transcript table
    if aggregation_level == 'probe':
        tx_column = 'transcript'
    else:
        tx_column = 'gene_name'
    keep_cols = ['x_pixel_coordinate', 'y_pixel_coordinate', 'cell_id', tx_column]
    df = g4x_obj.load_transcript_table(lazy=True, columns=keep_cols)

    ## make colormap
    unique_genes = df.select(tx_column).collect().unique().to_series().to_list()
    num_genes = len(unique_genes)
    base_cmap = plt.get_cmap('hsv', num_genes)
    palette = {gene: [int(255 * c) for c in base_cmap(i / num_genes)[:3]] for i, gene in enumerate(unique_genes)}

    df = df.with_columns(
        (pl.col('x_pixel_coordinate') / MIN_TILE_SIZE).cast(pl.Int32).alias('tile_y_coord'),
        (pl.col('y_pixel_coordinate') / MIN_TILE_SIZE).cast(pl.Int32).alias('tile_x_coord'),
        pl.col(tx_column).replace_strict(palette, default=[127, 127, 127]).alias('color'),
        pl.concat_list(['x_pixel_coordinate', 'y_pixel_coordinate']).alias('position'),
    )

    num_tiles_x, num_tiles_y = tg.GetPyramidTilesConfigData(IMAGE_RESOLUTION, MIN_TILE_SIZE)
    NUMBER_OF_LEVELS = tg.GetPyramidLevelsConfigData(IMAGE_RESOLUTION, MIN_TILE_SIZE)

    min_zoom_tiles_x = math.ceil(num_tiles_x / (pow(2, NUMBER_OF_LEVELS - 1)))
    min_zoom_tiles_y = math.ceil(num_tiles_y / (pow(2, NUMBER_OF_LEVELS - 1)))

    logger.info(f"""
        Final configurations for parsing:
            Image resolution: {IMAGE_RESOLUTION[0]} x {IMAGE_RESOLUTION[1]}
            Number of max zoom tiles: X = {num_tiles_x} | Y = {num_tiles_y}
            Number of min zoom tiles: X = {min_zoom_tiles_x} | Y = {min_zoom_tiles_y}
            Number of levels: {NUMBER_OF_LEVELS}
    """)

    tg.save_configuration_file(out_dir, IMAGE_RESOLUTION, MIN_TILE_SIZE, NUMBER_OF_LEVELS, palette)
    logger.info('Parsing and classifying tiles...')

    pq_args = []
    for level_index in reversed(range(NUMBER_OF_LEVELS + 1)):
        ## subsampling factor
        sampling_factor = sampling_fraction ** (NUMBER_OF_LEVELS - level_index)

        # factor for computing tile coordinates at this level
        scaling_factor = 2 ** (NUMBER_OF_LEVELS - level_index)
        current_tile_size = MIN_TILE_SIZE * scaling_factor
        x_num_of_tiles = math.ceil(IMAGE_RESOLUTION[0] / current_tile_size)
        y_num_of_tiles = math.ceil(IMAGE_RESOLUTION[1] / current_tile_size)

        # Ensure even numbers of tiles
        if x_num_of_tiles % 2 != 0:
            x_num_of_tiles += 1
        if y_num_of_tiles % 2 != 0:
            y_num_of_tiles += 1

        for tile_x_index in range(x_num_of_tiles):
            tileOutputDirPath = out_dir / f'{level_index}' / f'{tile_x_index}'
            os.makedirs(tileOutputDirPath, exist_ok=True)
            pq_args.append([tileOutputDirPath, y_num_of_tiles, scaling_factor])

            _ = (
                df.filter(
                    ((pl.col('tile_x_coord') // scaling_factor) == tile_x_index),
                )
                .select(['position', 'color', tx_column, 'cell_id', 'tile_y_coord'])
                .collect()
                .sample(fraction=sampling_factor)
                .write_parquet(tileOutputDirPath / 'tmp.parquet')
            )

    ctx = mp.get_context('spawn')
    with ctx.Pool(processes=n_threads) as pool:
        pool.starmap(tg._process_x, pq_args)

    logger.info('Tarring up.')
    if out_path.exists() or out_path.is_symlink():
        out_path.unlink()
    _ = tg.create_tar_from_directory(out_dir, out_path)
    shutil.rmtree(out_dir, ignore_errors=True)
