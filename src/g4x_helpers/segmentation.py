# from g4x_helpers.models import G4Xoutput
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

from pathlib import Path
import anndata as ad
import numpy as np
import polars as pl
import scanpy as sc
import shapely.affinity
import skimage.measure
from geopandas import GeoDataFrame
from rasterio.features import rasterize
from rasterio.transform import Affine
from shapely.geometry import Polygon
from skimage.measure._regionprops import RegionProperties

if TYPE_CHECKING:
    # This import is only for type checkers (mypy, PyCharm, etc.), not at runtime
    from g4x_helpers.models import G4Xoutput


def _make_cell_metadata(segmentation_props: pl.DataFrame, cell_by_protein: pl.DataFrame) -> pl.DataFrame:
    cell_metadata = segmentation_props.join(cell_by_protein, on='segmentation_cell_id', how='left')
    return cell_metadata


def _make_cell_by_gene(segmentation_props: pl.DataFrame, reads_new_labels: pl.DataFrame) -> pl.DataFrame:
    cell_by_gene = (
        reads_new_labels.filter(pl.col('segmentation_label') != 0)
        .group_by('segmentation_cell_id', 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index='segmentation_cell_id')
    )

    cell_by_gene = (
        segmentation_props.select('segmentation_cell_id')
        .join(cell_by_gene, on='segmentation_cell_id', how='left')
        .fill_null(0)
    )
    return cell_by_gene


def _make_adata(cell_by_gene: pl.DataFrame, cell_metadata: pl.DataFrame) -> ad.AnnData:
    X = cell_by_gene.drop('segmentation_cell_id').to_numpy()

    obs_df = cell_metadata.drop('segmentation_label').to_pandas()
    obs_df.index = obs_df.index.astype(str)

    adata = ad.AnnData(X=X, obs=obs_df)

    adata.var['gene_id'] = cell_by_gene.columns[1:]
    adata.var['modality'] = 'tx'

    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)

    return adata


def get_mask_properties(sample: 'G4Xoutput', mask: np.ndarray) -> pl.DataFrame:
    props = skimage.measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid
    for prop in props:
        label = prop.label  # The label (mask id)
        area = prop.area  # Area: count of pixels
        centroid = prop.centroid  # Centroid: (row, col) coordinates

        if sample._get_coord_order() == 'yx':
            cell_y = centroid[0]
            cell_x = centroid[1]
        else:
            cell_x = centroid[0]
            cell_y = centroid[1]

        prop_dict.append(
            {
                'segmentation_label': label,
                'segmentation_cell_id': f'{sample.sample_id}-{label}',
                'area': area,
                'cell_x': cell_x,
                'cell_y': cell_y,
            }
        )
        # prop_dict.append({'segmentation_label': label, 'area': area, 'cell_x': cell_x, 'cell_y': cell_y})

    prop_dict_df = pl.DataFrame(prop_dict)
    # prop_dict_df = prop_dict_df.with_columns(pl.col('segmentation_label')).cast(pl.Utf8).cast(pl.Categorical)
    return prop_dict_df


def assign_tx_to_mask_labels(sample: 'G4Xoutput', mask: np.ndarray) -> pl.DataFrame:
    reads = sample.load_transcript_table(return_polars=True)

    if sample._get_coord_order() == 'yx':
        coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    else:
        coord_order = ['x_pixel_coordinate', 'y_pixel_coordinate']

    tx_coords = reads[coord_order].to_numpy()
    cell_ids = mask[tx_coords[:, 0].astype(int), tx_coords[:, 1].astype(int)]

    reads = reads.with_columns(pl.lit(cell_ids).alias('segmentation_label'))
    reads = reads.with_columns(
        (f'{sample.sample_id}-' + pl.col('segmentation_label').cast(pl.String)).alias('segmentation_cell_id')
    )

    return reads


def image_mask_intensity_extraction(
    mask: np.ndarray,
    img: np.ndarray,
    *,
    bead_mask: Optional[np.ndarray] = None,
    label_prefix: Optional[str] = '',
    channel_label: Optional[str] = 'intensity_mean',
    lazy: Optional[bool] = True,
) -> Union[pl.DataFrame, pl.LazyFrame]:
    mask_flat = mask.ravel()
    img_flat = img.ravel()

    ## optional masking with beads
    if bead_mask is not None:
        bead_mask_flat = ~bead_mask.ravel()
        mask_flat = mask_flat[bead_mask_flat]
        img_flat = img_flat[bead_mask_flat]

    intensity_means = np.bincount(mask_flat, weights=img_flat)[1:] / np.bincount(mask_flat)[1:]

    if lazy:
        prop_tab = pl.LazyFrame(
            {'label': np.arange(start=1, stop=mask_flat.max() + 1), channel_label: intensity_means}
        ).with_columns((pl.lit(label_prefix) + pl.col('label').cast(pl.Utf8)).cast(pl.Categorical).alias('label'))
    else:
        prop_tab = pl.DataFrame(
            {'label': np.arange(start=1, stop=mask_flat.max() + 1), channel_label: intensity_means}
        ).with_columns((pl.lit(label_prefix) + pl.col('label').cast(pl.Utf8)).cast(pl.Categorical).alias('label'))

    return prop_tab


def join_frames(lf_left, lf_right):
    return lf_left.join(lf_right, on='label', how='inner')


def _region_props_to_polygons(region_props: RegionProperties) -> list[Polygon]:
    mask = np.pad(region_props.image, 1)
    contours = skimage.measure.find_contours(mask, 0.5)

    # shapes with <= 3 vertices, i.e. lines, can't be converted into a polygon
    polygons = [Polygon(contour[:, [1, 0]]) for contour in contours if contour.shape[0] >= 4]

    yoff, xoff, *_ = region_props.bbox
    return [shapely.affinity.translate(poly, xoff, yoff) for poly in polygons]


def _vectorize_mask(mask: np.ndarray, nudge: bool = True) -> GeoDataFrame:
    if mask.max() == 0:
        return GeoDataFrame(geometry=[])

    regions = skimage.measure.regionprops(mask)

    polygons_list = [_region_props_to_polygons(region) for region in regions]
    geoms = [poly for polygons in polygons_list for poly in polygons]
    labels = [region.label for i, region in enumerate(regions) for _ in range(len(polygons_list[i]))]
    gdf = GeoDataFrame({'label': labels}, geometry=geoms)

    if nudge:
        gdf['geometry'] = gdf['geometry'].translate(xoff=-0.5, yoff=-0.5)

    return gdf


def rasterize_polygons(gdf: GeoDataFrame, target_shape: tuple) -> np.ndarray:
    height, width = target_shape
    transform = Affine.identity()

    shapes = ((geom, int(lbl)) for geom, lbl in zip(gdf.geometry, gdf['label']))

    label_array = rasterize(
        shapes=shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,  # background value
        dtype='int32',
    )

    return label_array


def extract_image_signals(
    sample: 'G4Xoutput',
    mask: np.ndarray,
    lazy: bool = False,
    include_channels: list = None,
    exclude_channels: list = None,
) -> Union[pl.DataFrame, pl.LazyFrame]:
    if include_channels:
        signal_list = include_channels
    else:
        signal_list = ['nuclear', 'eosin'] + sample.proteins

    if isinstance(exclude_channels, (str, type(None))):
        exclude_channels = [exclude_channels]

    signal_list = [item for item in signal_list if item not in exclude_channels]

    for i, signal_name in enumerate(signal_list):
        print(f'Extracting {signal_name} signal...')
        if signal_name not in exclude_channels:
            if signal_name in ('nuclear', 'eosin'):
                # parent_directory = 'h_and_e'
                channel_name_map = {'nuclear': 'nuclearstain', 'eosin': 'cytoplasmicstain'}
            else:
                # parent_directory = 'protein'
                channel_name_map = {protein: protein for protein in sample.proteins}

            signal_img = sample.load_image('PanCK', thumbnail=False)
            # signal_img = sample._get_image(parent_directory=parent_directory, pattern=f'{signal_name}.jp2')
            ch_label = f'{channel_name_map[signal_name]}_intensity_mean'

            prop_tab = image_mask_intensity_extraction(
                mask,
                signal_img,
                bead_mask=None,
                label_prefix=f'{sample.sample_id}-',
                channel_label=ch_label,
                lazy=lazy,
            )

            if i == 0:
                signal_df = prop_tab
            else:
                signal_df = join_frames(signal_df, prop_tab)

    signal_df = signal_df.cast({'label': pl.String}).rename({'label': 'segmentation_cell_id'})

    return signal_df


def _create_custom_out(sample: 'G4Xoutput', out_dir=None, parent_folder=None, file_name=None):
    if out_dir is None:
        custom_out = sample.run_base / parent_folder / 'custom_segmentation'
    else:
        custom_out = Path(out_dir) / parent_folder

    # Ensure the directory exists
    custom_out.mkdir(parents=True, exist_ok=True)

    # Prepare output file path
    outfile = custom_out / file_name

    # If the output file already exists, delete it
    if outfile.exists():
        outfile.unlink()

    return outfile