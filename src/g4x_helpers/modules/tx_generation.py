from __future__ import annotations

import json
import math
import os
import tarfile
from pathlib import Path

import polars as pl


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
    from ..modules.g4x_viewer import MetadataSchema_pb2 as MetadataSchema

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
