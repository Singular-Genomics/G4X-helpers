from __future__ import annotations
from typing import TYPE_CHECKING
import json
import logging
import math
import os
import shutil
import tarfile
import multiprocessing as mp
from pathlib import Path
import matplotlib.pyplot as plt
import polars as pl

import g4x_helpers.utils as utils

mp.set_start_method('spawn', force=True)

if TYPE_CHECKING:
    # This import is only for type checkers (mypy, PyCharm, etc.), not at runtime
    from g4x_helpers.models import G4Xoutput


def create_tar_from_directory(directory_path: str | Path, archive_name: str | Path) -> None:
    if not os.path.isdir(directory_path):
        raise ValueError(f'The directory {directory_path} does not exist.')
    with tarfile.open(archive_name, 'w') as tar:
        tar.add(directory_path, arcname=os.path.basename(directory_path))


def GetPyramidTilesConfigData(image_resolution: tuple[int, int], tile_size: int) -> tuple[int, int]:
    num_tiles_x = math.ceil(image_resolution[0] / tile_size)
    num_tiles_y = math.ceil(image_resolution[1] / tile_size)

    if not num_tiles_x % 2 == 0:
        num_tiles_x += 1

    if not num_tiles_y % 2 == 0:
        num_tiles_y += 1

    return num_tiles_x, num_tiles_y


def GetPyramidLevelsConfigData(image_resolution: tuple[int, int], tile_size: int, min_tiles_number: int = 16) -> int:
    min_num_levels = 0

    current_level = 1
    while min_num_levels == 0:
        level_tile_size = tile_size * pow(2, current_level)
        level_tiles_x = math.ceil(image_resolution[0] / level_tile_size)
        level_tiles_y = math.ceil(image_resolution[1] / level_tile_size)
        if not level_tiles_x % 2 == 0:
            level_tiles_x += 1
        if not level_tiles_y % 2 == 0:
            level_tiles_y += 1
        if level_tiles_x * level_tiles_y <= min_tiles_number:
            min_num_levels = current_level
            break
        current_level += 1

    return min_num_levels


def save_configuration_file(
    outputDirPath: str, image_resolution: tuple[int, int], min_tile_size: int, number_of_levels: int, palette: dict
) -> None:
    start_tile_size = min_tile_size * pow(2, number_of_levels)

    config_data = {
        'layer_height': image_resolution[0],
        'layer_width': image_resolution[1],
        'layers': number_of_levels,
        'tile_size': start_tile_size,
        'color_map': [{'gene_name': key, 'color': value} for key, value in palette.items()],
    }

    with open(f'{outputDirPath}/config.json', 'w') as json_file:
        json.dump(config_data, json_file, indent=2)


def _process_x(tileOutputDirPath: Path, y_num_of_tiles: int, scaling_factor: int) -> None:
    from . import MetadataSchema_pb2 as MetadataSchema

    for tile_y_index in range(y_num_of_tiles):
        outputTileData = MetadataSchema.TileData()

        df_current = (
            pl.scan_parquet(tileOutputDirPath / 'tmp.parquet')
            .filter(((pl.col('tile_y_coord') // scaling_factor) == tile_y_index))
            .drop('tile_y_coord')
            .collect()
        )
        # Iterate over rows directly
        ## this can potentially be done lazily with this PR: https://github.com/pola-rs/polars/pull/23980
        for position, color, gene, cell_id in df_current.iter_rows():
            outputPointData = outputTileData.pointsData.add()
            _ = outputPointData.position.extend(position)
            _ = outputPointData.color.extend(color)
            outputPointData.geneName = gene
            outputPointData.cellId = str(cell_id)

        with open(f'{tileOutputDirPath}/{tile_y_index}.bin', 'wb') as file:
            _ = file.write(outputTileData.SerializeToString())


def tx_converter(
    sample: 'G4Xoutput',
    out_path: str | Path,
    *,
    aggregation_level: str = 'gene',
    n_threads: int = 4,
    sampling_fraction: float = 0.2,
    logger: logging.Logger | None = None,
    log_level: int = logging.DEBUG,
) -> None:
    if logger is None:
        logger = utils.setup_logger(
            'tx_converter',
            stream_logger=True,
            stream_level=log_level,
            file_logger=True,
            file_level=logging.DEBUG,
            out_dir=out_path.parent,
            clear_handlers=True,
        )

    ## prelims
    IMAGE_RESOLUTION = sample.shape
    if sample.platform == 'g4x-2lane':
        MIN_TILE_SIZE = 1028
    else:
        MIN_TILE_SIZE = 512
    out_dir = sample.run_base / 'g4x_viewer_temp'
    os.makedirs(out_dir, exist_ok=True)

    ## get transcript table
    if aggregation_level == 'probe':
        tx_column = 'transcript'
    else:
        tx_column = 'gene_name'
    keep_cols = ['x_pixel_coordinate', 'y_pixel_coordinate', 'cell_id', tx_column]
    df = sample.load_transcript_table(lazy=True, columns=keep_cols)

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

    num_tiles_x, num_tiles_y = GetPyramidTilesConfigData(IMAGE_RESOLUTION, MIN_TILE_SIZE)
    NUMBER_OF_LEVELS = GetPyramidLevelsConfigData(IMAGE_RESOLUTION, MIN_TILE_SIZE)

    min_zoom_tiles_x = math.ceil(num_tiles_x / (pow(2, NUMBER_OF_LEVELS - 1)))
    min_zoom_tiles_y = math.ceil(num_tiles_y / (pow(2, NUMBER_OF_LEVELS - 1)))

    logger.info(f"""
        Final configurations for parsing:
            Image resolution: {IMAGE_RESOLUTION[0]} x {IMAGE_RESOLUTION[1]}
            Number of max zoom tiles: X = {num_tiles_x} | Y = {num_tiles_y}
            Number of min zoom tiles: X = {min_zoom_tiles_x} | Y = {min_zoom_tiles_y}
            Number of levels: {NUMBER_OF_LEVELS}
    """)

    save_configuration_file(out_dir, IMAGE_RESOLUTION, MIN_TILE_SIZE, NUMBER_OF_LEVELS, palette)
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
        pool.starmap(_process_x, pq_args)

    logger.info('Tarring up.')
    if out_path.exists() or out_path.is_symlink():
        out_path.unlink()
    _ = create_tar_from_directory(out_dir, out_path)
    shutil.rmtree(out_dir, ignore_errors=True)
