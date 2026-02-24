from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from skimage import measure
from tqdm import tqdm

# from .segment import get_cell_ids

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput

PX_SIZE_MICRONS = 0.3125


def get_cell_ids(sample_id: str, mask) -> pl.DataFrame:
    print('Getting cell-IDs from segmentation mask.')
    seg_ids = np.unique(mask[mask != 0])
    cell_ids = (
        pl.Series(name='segmentation_label', values=seg_ids)
        .to_frame()
        .with_columns((f'{sample_id}-' + pl.col('segmentation_label').cast(pl.String)).alias('label'))
        .select(['label', 'segmentation_label'])
    )
    return cell_ids


def get_mask_properties(cell_ids: pl.DataFrame | None, mask: np.ndarray) -> pl.DataFrame:
    print('Getting regionprops.')
    props = measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar
    for prop in tqdm(props, desc='Extracting mask properties'):
        label = prop.label  # The label (mask id)

        cell_y, cell_x = prop.centroid  # coordinate order: 'yx' (row, col)

        px_to_um_area = PX_SIZE_MICRONS**2
        area_um = prop.area * px_to_um_area  # prop.area is in pixels

        prop_dict.append(
            {
                'seg_cell_id': label,
                'area_um': area_um,
                'cell_x': cell_x,
                'cell_y': cell_y,
            }
        )
    schema = {
        'seg_cell_id': pl.Int32,
        'area_um': pl.Float32,
        'cell_x': pl.Float32,
        'cell_y': pl.Float32,
    }
    prop_dict_df = pl.DataFrame(prop_dict, schema=schema).sort('seg_cell_id')

    # prop_dict_df = cell_ids.join(prop_dict_df, on='segmentation_label', how='left')
    return prop_dict_df


def build_g4x_cell_properties(g4x_obj: 'G4Xoutput') -> pl.DataFrame:
    # first load expanded mask
    mask = g4x_obj.load_segmentation(expanded=True)
    cell_ids = get_cell_ids(g4x_obj.sample_id, mask)
    cell_props_expanded = get_mask_properties(cell_ids, mask)

    # now load non-expanded mask
    mask = g4x_obj.load_segmentation(expanded=False)
    cell_props_nuc = get_mask_properties(cell_ids, mask)
    del mask

    df = cell_props_nuc.rename({'area_um': 'nuclei_area_um'})
    df_exp = cell_props_expanded.rename({'area_um': 'wholecell_area_um'}).drop(
        ['segmentation_label', 'cell_x', 'cell_y']
    )

    segmentation_props = df.join(df_exp, on='label').select(
        ['label', 'segmentation_label', 'cell_x', 'cell_y', 'nuclei_area_um', 'wholecell_area_um']
    )

    return segmentation_props


def build_custom_cell_properties(g4x_obj: 'G4Xoutput', mask: np.ndarray) -> pl.DataFrame:
    cell_ids = get_cell_ids(g4x_obj.sample_id, mask)
    cell_props = get_mask_properties(cell_ids, mask)
    del mask

    segmentation_props = cell_props.select(['label', 'segmentation_label', 'cell_x', 'cell_y', 'area_um'])

    return segmentation_props
