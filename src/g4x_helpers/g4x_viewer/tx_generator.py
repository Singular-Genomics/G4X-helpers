from __future__ import annotations
from typing import TYPE_CHECKING
import json
import logging
import multiprocessing
import math
import os
import shutil
import tarfile
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from tqdm import tqdm

import g4x_helpers.utils as utils
from . import MetadataSchema_pb2 as MetadataSchema

if TYPE_CHECKING:
    # This import is only for type checkers (mypy, PyCharm, etc.), not at runtime
    from g4x_helpers.models import G4Xoutput


class NumpyEncoder(json.JSONEncoder):
    """Custom encoder for numpy data types"""

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):
            return int(obj)

        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)

        elif isinstance(obj, (np.complex_, np.complex64, np.complex128)):
            return {'real': obj.real, 'imag': obj.imag}

        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()

        elif isinstance(obj, (np.bool_)):
            return bool(obj)

        elif isinstance(obj, (np.void)):
            return None

        return json.JSONEncoder.default(self, obj)


def create_tar_from_directory(directory_path, archive_name):
    # Ensure the directory path exists
    if not os.path.isdir(directory_path):
        raise ValueError(f'The directory {directory_path} does not exist.')

    with tarfile.open(archive_name, 'w') as tar:
        # Add the directory to the tar file
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


def GenerateInitialTiles(
    number_of_tiles_x: int, number_of_tiles_y: int, levels: int, df: pl.DataFrame, logger: logging.Logger
) -> dict:
    # Initialize empty data structure
    tile_x_data = {}
    for tile_x_index in range(number_of_tiles_x):
        logger.info(f'Processing tile_x = {tile_x_index}')
        tile_y_data = {}
        for tile_y_index in range(number_of_tiles_y):
            logger.info(f'Processing tile_y = {tile_y_index}')
            tile_y_data.update(
                {
                    tile_y_index: {
                        'points': (
                            df.filter(pl.col('tile_x_coord') == tile_x_index, pl.col('tile_y_coord') == tile_y_index)
                            .select(['position', 'color', 'transcript_condensed', 'cell_id'])
                            .to_dicts()
                        )
                    }
                }
            )

        tile_x_data.update({tile_x_index: tile_y_data})

    return {levels: tile_x_data}


def GeneratePyramidData(
    parsed_csv_data: dict, levels: int, image_resolution: tuple[int, int], tile_size: int, logger: logging.Logger
) -> dict:
    for level_index in reversed(range(levels)):
        current_tile_size = tile_size * pow(2, (levels - level_index))
        x_num_of_tiles = math.ceil(image_resolution[0] / current_tile_size)
        y_num_of_tiles = math.ceil(image_resolution[1] / current_tile_size)

        logger.info(f"""
            --- Generating data for pyramid level: {level_index}
            --- Current tile size: {current_tile_size}
            --- Number of tiles: X = {x_num_of_tiles} | Y = {y_num_of_tiles}
        """)

        if not x_num_of_tiles % 2 == 0:
            x_num_of_tiles += 1
        if not y_num_of_tiles % 2 == 0:
            y_num_of_tiles += 1

        tile_x_data = {}
        for tile_x_index in tqdm(range(x_num_of_tiles)):
            tile_y_data = {}
            for tile_y_index in range(y_num_of_tiles):
                try:
                    sub_tile_0_points = np.array(
                        parsed_csv_data[level_index + 1][tile_x_index * 2][tile_y_index * 2]['points']
                    )
                    sub_tile_1_points = np.array(
                        parsed_csv_data[level_index + 1][(tile_x_index * 2) + 1][tile_y_index * 2]['points']
                    )
                    sub_tile_2_points = np.array(
                        parsed_csv_data[level_index + 1][(tile_x_index * 2) + 1][(tile_y_index * 2) + 1]['points']
                    )
                    sub_tile_3_points = np.array(
                        parsed_csv_data[level_index + 1][tile_x_index * 2][(tile_y_index * 2) + 1]['points']
                    )

                    sparse_points = (
                        sub_tile_0_points[
                            np.random.choice(len(sub_tile_0_points), int(len(sub_tile_0_points) * 0.2), replace=False)
                        ].tolist()
                        + sub_tile_1_points[
                            np.random.choice(len(sub_tile_1_points), int(len(sub_tile_1_points) * 0.2), replace=False)
                        ].tolist()
                        + sub_tile_2_points[
                            np.random.choice(len(sub_tile_2_points), int(len(sub_tile_2_points) * 0.2), replace=False)
                        ].tolist()
                        + sub_tile_3_points[
                            np.random.choice(len(sub_tile_3_points), int(len(sub_tile_3_points) * 0.2), replace=False)
                        ].tolist()
                    )

                    tile_y_data.update(
                        {
                            tile_y_index: {
                                'points': sparse_points,
                            }
                        }
                    )
                except Exception as e:
                    logger.error(e)
                    logger.error(traceback.format_exc())
                    tile_y_data.update({tile_y_index: {'points': []}})
                    continue

            tile_x_data.update({tile_x_index: tile_y_data})
        parsed_csv_data.update({level_index: tile_x_data})

    return parsed_csv_data


def write_tile_file_from_df(
    df_tile: pl.DataFrame, level: int, tile_x: int, tile_y: int, outputDirPath: str, outputDirName: str
) -> None:
    """
    Write .bin file directly from a filtered Polars DataFrame.
    No intermediate dictionary creation is done here.

    df_tile must contain columns: position, color, transcript_condensed, cell_id.
    """
    fullOutputDirPath = f'{outputDirPath}/{outputDirName}'
    tileOutputDirPath = f'{fullOutputDirPath}/{level}/{tile_x}'
    os.makedirs(tileOutputDirPath, exist_ok=True)

    outputTileData = MetadataSchema.TileData()

    # Iterate over rows directly
    for position, color, gene, cell_id in df_tile.select(
        ['position', 'color', 'transcript_condensed', 'cell_id']
    ).iter_rows():
        outputPointData = outputTileData.pointsData.add()
        try:
            outputPointData.position.extend(position)
            outputPointData.color.extend(color)
            outputPointData.geneName = gene
            outputPointData.cellId = str(cell_id)
        except Exception as e:
            print(e)
            continue

    with open(f'{tileOutputDirPath}/{tile_y}.bin', 'wb') as file:
        file.write(outputTileData.SerializeToString())


def process_x_tile(df_current, tile_x_index, y_num_of_tiles, factor, level_index, outputDirPath, outputDirName) -> None:
    for tile_y_index in range(y_num_of_tiles):
        try:
            # Filter DataFrame for points belonging to this tile at the current level.
            df_tile = df_current.filter(((pl.col('tile_y_coord') // factor) == tile_y_index))

            # Write the tile file directly from the df
            write_tile_file_from_df(
                df_tile=df_tile,
                level=level_index,
                tile_x=tile_x_index,
                tile_y=tile_y_index,
                outputDirPath=outputDirPath,
                outputDirName=outputDirName,
            )

        except Exception:
            # Write empty tile in case of error
            empty_df = df_current.head(0)  # empty dataframe with same schema
            write_tile_file_from_df(
                df_tile=empty_df,
                level=level_index,
                tile_x=tile_x_index,
                tile_y=tile_y_index,
                outputDirPath=outputDirPath,
                outputDirName=outputDirName,
            )


def save_multi_files(
    outputMetadataValues: dict, outputDirPath: str, outputDirName: str, logger: logging.Logger
) -> None:
    fullOutputDirPath = f'{outputDirPath}/{outputDirName}'

    for level, level_data in outputMetadataValues.items():
        for y_coord, y_coord_data in level_data.items():
            for x_coord, x_coord_data in y_coord_data.items():
                try:
                    points = x_coord_data['points']
                    outputTileData = MetadataSchema.TileData()

                    for point in points:
                        outputPointData = outputTileData.pointsData.add()

                        try:
                            outputPointData.position.extend(point['position'])
                            outputPointData.color.extend(point['color'])
                            outputPointData.geneName = point['transcript_condensed']
                            outputPointData.cellId = f'{point["cell_id"]}'
                        except Exception as e:
                            logger.error(e)
                            logger.error(traceback.format_exc())
                            continue

                    tileOutputDirPath = f'{fullOutputDirPath}/{level}/{y_coord}/'
                    os.makedirs(tileOutputDirPath, exist_ok=True)
                    with open(f'{tileOutputDirPath}/{x_coord}.bin', 'wb') as file:
                        file.write(outputTileData.SerializeToString())

                except Exception as e:
                    logger.error(e)
                    logger.error(traceback.format_exc())
                    continue


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


def tx_converter(
    sample: 'G4Xoutput',
    out_path: str | Path,
    *,
    aggregation_level: str = 'gene',
    n_threads: int = 4,
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
        tx_column = 'transcript_condensed'
    keep_cols = ['x_coord_shift', 'y_coord_shift', 'cell_id', tx_column]
    df = sample.load_transcript_table(columns=keep_cols)

    ## make colormap
    unique_genes = df[tx_column].unique().to_list()
    num_genes = len(unique_genes)
    base_cmap = plt.get_cmap('hsv', num_genes)
    palette = {gene: [int(255 * c) for c in base_cmap(i / num_genes)[:3]] for i, gene in enumerate(unique_genes)}

    ## actual processing
    df = (
        df.with_columns(
            (pl.col('x_coord_shift') / MIN_TILE_SIZE).cast(pl.Int32).alias('tile_x_coord'),
            (pl.col('y_coord_shift') / MIN_TILE_SIZE).cast(pl.Int32).alias('tile_y_coord'),
            pl.col(tx_column).replace(palette, default=[127, 127, 127]).alias('color'),
            pl.concat_list(['y_coord_shift', 'x_coord_shift']).alias('position'),
        )
        .rename({tx_column: 'transcript_condensed'})
        .collect()
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

    df_current = df
    sampling_fraction = 0.2

    for level_index in reversed(range(NUMBER_OF_LEVELS + 1)):
        if level_index < NUMBER_OF_LEVELS and sampling_fraction < 1.0:
            df_current = df_current.sample(fraction=sampling_fraction)

        # factor for computing tile coordinates at this level
        factor = 2 ** (NUMBER_OF_LEVELS - level_index)
        current_tile_size = MIN_TILE_SIZE * factor
        x_num_of_tiles = math.ceil(IMAGE_RESOLUTION[0] / current_tile_size)
        y_num_of_tiles = math.ceil(IMAGE_RESOLUTION[1] / current_tile_size)

        # Ensure even numbers of tiles
        if x_num_of_tiles % 2 != 0:
            x_num_of_tiles += 1
        if y_num_of_tiles % 2 != 0:
            y_num_of_tiles += 1

        logger.info(f"""
            --- Generating data for pyramid level: {level_index}
            --- Current tile size: {current_tile_size}
            --- Number of tiles: X = {x_num_of_tiles} | Y = {y_num_of_tiles}
            --- factor: {factor}
        """)

        # Create output directory for this level
        os.makedirs(out_dir / f'{level_index}', exist_ok=True)

        pq_args = []
        for tile_x_index in range(x_num_of_tiles):
            df_tile = df_current.filter(((pl.col('tile_x_coord') // factor) == tile_x_index))
            pq_args.append((df_tile, tile_x_index, y_num_of_tiles, factor, level_index, out_dir, sample.sample_id))
        with multiprocessing.Pool(processes=n_threads) as pool:
            _ = pool.starmap(process_x_tile, pq_args)

    logger.info('Tarring up.')
    if out_path.exists() or out_path.is_symlink():
        out_path.unlink()
    _ = create_tar_from_directory(out_dir, out_path)
    shutil.rmtree(out_dir, ignore_errors=True)
