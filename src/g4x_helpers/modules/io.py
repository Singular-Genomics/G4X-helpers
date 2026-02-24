from __future__ import annotations

from typing import TYPE_CHECKING

import geopandas
import numpy as np
import polars as pl
from tqdm import tqdm

from .. import utils

if TYPE_CHECKING:
    from geopandas.geodataframe import GeoDataFrame


def build_sample_metadata(p, meta):
    meta = meta.copy()
    smp_id = meta['sample_position'] + meta['lane'][-2:]

    sample_info, run_info = _parse_samplesheet(p)

    ssheet_info = sample_info.filter(
        pl.col('lane') == int(smp_id[-1]), pl.col('sample_position') == smp_id[0]
    ).to_dicts()[0]
    ri_keys = {'run_name', 'user_name', 'assay'}
    filtered = [d for d in run_info.to_dicts() if d['Key'] in ri_keys]
    run_info = {d['Key']: d['Value'] for d in filtered}

    meta.update(run_info)
    meta.update(ssheet_info)

    order = [
        'run_name',
        'machine',
        'run_id',
        'platform',
        'assay',
        'user_name',
        'time_of_creation',
        'fc',
        'lane',
        'sample_position',
        'tissue_type',
        'block',
        'transcript_panel',
        'transcript_addon',
        'protein_panel',
        'protein_addon',
        'software',
        'software_version',
    ]

    return {k: meta[k] for k in order if k in meta}


def _parse_samplesheet(path):
    df = pl.read_csv(path, has_header=False, row_index_name='index')

    data_index = df.filter(pl.col('column_1') == '[Data]')['index'][0]

    run_info = (
        df.filter(pl.col('index').is_between(0, data_index, closed='none'))
        .select('column_1', 'column_2')
        .rename({'column_1': 'Key', 'column_2': 'Value'})
    )

    run_info = run_info.with_columns(pl.col('Key').str.to_lowercase().str.replace(' ', '_'))

    sample_info = pl.read_csv(path, skip_rows=data_index + 1)
    sample_info = sample_info.rename({c: c.lower().replace(' ', '_') for c in sample_info.columns})

    # TODO protein_addon is not always the last column ...
    final_col = np.where(np.array(sample_info.columns) == 'protein_addon')[0][0]
    sample_info = sample_info.select(pl.col(sample_info.columns[: final_col + 1])).fill_null('none')

    return sample_info, run_info


def import_segmentation(seg_path: str, expected_shape: tuple[int], labels_key: str | None = None) -> np.ndarray:
    SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}
    ## load new segmentation
    cell_labels = utils.validate_path(seg_path, must_exist=True, is_dir_ok=False, is_file_ok=True)
    suffix = cell_labels.suffix.lower()

    if suffix not in SUPPORTED_MASK_FILETYPES:
        raise ValueError(f'{suffix} is not a supported file type.')

    if suffix == '.npz':
        with np.load(cell_labels) as labels:
            available_keys = list(labels.keys())

            if labels_key:  # if a key is specified
                if labels_key not in labels:
                    raise KeyError(f"Key '{labels_key}' not found in .npz; available keys: {available_keys}")
                seg = labels[labels_key]

            else:
                if len(available_keys) != 1:
                    raise ValueError(
                        f"Found multiple keys in .npz: {available_keys}.\nPlease specify a key using 'labels_key'"
                    )
                seg = labels[available_keys[0]]

    elif suffix == '.npy':
        # .npy: directly returns the array, no context manager available
        if labels_key is not None:
            print('file is .npy, ignoring provided labels_key.')
        seg = np.load(cell_labels, allow_pickle=False)

    elif suffix == '.geojson':
        gdf = geopandas.read_file(cell_labels)

        if labels_key is not None:
            if labels_key not in gdf.columns:
                raise KeyError(f"Column '{labels_key}' not found in GeoJSON; available columns: {gdf.columns.tolist()}")

            # ensure that a column named 'label' exists
            gdf['label'] = gdf[labels_key]

        else:
            if 'label' not in gdf.columns:
                raise ValueError(
                    "No column named 'label' found in GeoJSON. Please specify which column to use for labels via labels_key."
                )

        print('Rasterizing provided GeoDataFrame.')
        seg = rasterize_polygons(gdf=gdf, target_shape=expected_shape)

    # validate shape for final numpy arrays
    if seg.shape != expected_shape:
        raise ValueError(f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}')

    return seg


def rasterize_polygons(gdf: 'GeoDataFrame', target_shape: tuple) -> np.ndarray:
    from rasterio.features import rasterize
    from rasterio.transform import Affine

    height, width = target_shape
    transform = Affine.identity()

    # wrap the zip in tqdm; total=len(gdf) gives a proper progress bar
    wrapped = tqdm(zip(gdf.geometry, gdf['label']), total=len(gdf), desc='Rasterizing polygons')
    # feed that wrapped iterator into rasterize
    shapes = ((geom, int(lbl)) for geom, lbl in wrapped)

    label_array = rasterize(
        shapes=shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,  # background value
        dtype='int32',
    )

    return label_array
