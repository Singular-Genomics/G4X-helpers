from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import anndata as ad
import geopandas
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
from tqdm import tqdm

if TYPE_CHECKING:
    from geopandas.geodataframe import GeoDataFrame

    from g4x_helpers.models import G4Xoutput

SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}


def try_load_segmentation(
    cell_labels: Path, expected_shape: tuple[int], labels_key: str | None = None
) -> np.ndarray | geopandas.GeoDataFrame:
    ## load new segmentation
    cell_labels = Path(cell_labels)
    suffix = cell_labels.suffix.lower()

    if cell_labels.suffix not in SUPPORTED_MASK_FILETYPES:
        raise ValueError(f'{cell_labels.suffix} is not a supported file type.')

    if suffix == '.npz':
        with np.load(cell_labels) as labels:
            available_keys = list(labels.keys())

            if labels_key:  # if a key is specified
                if labels_key not in labels:
                    raise KeyError(f"Key '{labels_key}' not found in .npz; available keys: {available_keys}")
                seg = labels[labels_key]

            else:  # if no key specified,
                if len(labels) == 1:  # and only one key is available, use that key
                    seg = labels[available_keys[0]]
                else:  # and multiple keys are available, raise an error
                    if len(labels) > 1:
                        raise ValueError(
                            f"Multiple keys found in .npz: {available_keys}.\nPlease specify a key using 'labels_key'"
                        )
                    seg = labels

    elif suffix == '.npy':
        # .npy: directly returns the array, no context manager available
        if labels_key is not None:
            print('file is .npy, ignoring provided labels_key.')
        seg = np.load(cell_labels)

    elif suffix == '.geojson':
        seg = geopandas.read_file(cell_labels)

        if labels_key is not None:
            if labels_key not in seg.columns:
                raise KeyError(f"Column '{labels_key}' not found in GeoJSON; available columns: {seg.columns.tolist()}")

            # ensure that a coliumn named 'label' exists
            seg['label'] = seg[labels_key]

        else:
            if 'label' not in seg.columns:
                raise ValueError(
                    "No column named 'label' found in GeoJSON. Please specify which column to use for labels via labels_key."
                )

    # only validate shape for numpy arrays
    if isinstance(seg, np.ndarray):
        assert seg.shape == expected_shape, (
            f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}'
        )

    return seg


def _make_cell_metadata(segmentation_props: pl.DataFrame, cell_by_protein: pl.DataFrame) -> pl.DataFrame:
    cell_metadata = segmentation_props.join(cell_by_protein, on='cell_id', how='left')
    return cell_metadata


def _make_cell_by_gene(segmentation_props: pl.DataFrame, reads_new_labels: pl.DataFrame) -> pl.DataFrame:
    cell_by_gene = (
        reads_new_labels.filter(pl.col('segmentation_label') != 0)
        .group_by('cell_id', 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index='cell_id')
    )

    cell_by_gene = segmentation_props.select('cell_id').join(cell_by_gene, on='cell_id', how='left').fill_null(0)
    return cell_by_gene


def _make_adata(cell_by_gene: pl.DataFrame, cell_metadata: pl.DataFrame) -> ad.AnnData:
    X = cell_by_gene.drop('cell_id').to_numpy()

    obs_df = cell_metadata.drop('segmentation_label').to_pandas()
    obs_df.index = obs_df.index.astype(str)

    adata = ad.AnnData(X=X, obs=obs_df)

    adata.var['gene_id'] = cell_by_gene.columns[1:]
    adata.var['modality'] = 'tx'
    adata.obs_names = adata.obs['cell_id']
    adata.var_names = adata.var['gene_id']

    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)

    return adata


def get_mask_properties(sample: 'G4Xoutput', mask: np.ndarray) -> pl.DataFrame:
    props = skimage.measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar
    for prop in tqdm(props, desc='Extracting mask properties'):
        label = prop.label  # The label (mask id)
        area = prop.area  # Area: count of pixels
        centroid = prop.centroid  # Centroid: (row, col)

        # assuming coordinate order: 'yx':
        cell_y, cell_x = centroid

        prop_dict.append(
            {
                'segmentation_label': label,
                'cell_id': f'{sample.sample_id}-{label}',
                'area': area,
                'cell_x': cell_x,
                'cell_y': cell_y,
            }
        )

    prop_dict_df = pl.DataFrame(prop_dict)
    return prop_dict_df


def assign_tx_to_mask_labels(sample: 'G4Xoutput', mask: np.ndarray) -> pl.DataFrame:
    reads = sample.load_transcript_table()

    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']

    tx_coords = reads[coord_order].to_numpy()
    cell_ids = mask[tx_coords[:, 0].astype(int), tx_coords[:, 1].astype(int)]

    reads = reads.with_columns(pl.lit(cell_ids).alias('segmentation_label'))
    reads = reads.with_columns((f'{sample.sample_id}-' + pl.col('segmentation_label').cast(pl.String)).alias('cell_id'))

    return reads


def image_mask_intensity_extraction(
    mask: np.ndarray,
    img: np.ndarray,
    *,
    bead_mask: np.ndarray | None = None,
    label_prefix: str = '',
    channel_label: str = 'intensity_mean',
    lazy: bool = True,
) -> pl.DataFrame | pl.LazyFrame:
    mask_flat = mask.ravel()
    img_flat = img.ravel()

    ## optional masking with beads
    if bead_mask is not None:
        bead_mask_flat = ~bead_mask.ravel()
        mask_flat = mask_flat[bead_mask_flat]
        img_flat = img_flat[bead_mask_flat]

    # intensity_means = np.bincount(mask_flat, weights=img_flat)[1:] / np.bincount(mask_flat)[1:]
    counts = np.bincount(mask_flat)[1:]
    sums = np.bincount(mask_flat, weights=img_flat)[1:]

    # Avoid division by zero:
    intensity_means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts != 0)

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

    geoms = []
    labels = []
    # Wrap the iteration in tqdm to show progress
    for region in tqdm(regions, desc='Vectorizing regions'):
        polys = _region_props_to_polygons(region)
        geoms.extend(polys)
        # add the region label once per polygon
        labels.extend([region.label] * len(polys))

    gdf = GeoDataFrame({'label': labels}, geometry=geoms)

    if nudge:
        # GeoSeries.translate works elementwise
        gdf['geometry'] = gdf['geometry'].translate(xoff=-0.5, yoff=-0.5)

    return gdf


def rasterize_polygons(gdf: GeoDataFrame, target_shape: tuple) -> np.ndarray:
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


def extract_image_signals(
    sample: G4Xoutput, mask: np.ndarray, lazy: bool = False, signal_list: list[str] | None = None, cached: bool = False
) -> pl.DataFrame | pl.LazyFrame:
    if signal_list is None:
        signal_list = ['nuclear', 'eosin'] + sample.proteins

    channel_name_map = {protein: protein for protein in sample.proteins}
    channel_name_map['nuclear'] = 'nuclearstain'
    channel_name_map['eosin'] = 'cytoplasmicstain'

    for signal_name in tqdm(signal_list, desc='Extracting protein signal'):
        if signal_name not in ['nuclear', 'eosin']:
            image_type = 'protein'
            protein = signal_name
        else:
            image_type = signal_name
            protein = None

        signal_img = sample.load_image_by_type(image_type, thumbnail=False, protein=protein, cached=cached)

        ch_label = f'{channel_name_map[signal_name]}_intensity_mean'

        prop_tab = image_mask_intensity_extraction(
            mask,
            signal_img,
            bead_mask=None,
            label_prefix=f'{sample.sample_id}-',
            channel_label=ch_label,
            lazy=lazy,
        )

        if signal_name == signal_list[0]:
            signal_df = prop_tab
        else:
            signal_df = join_frames(signal_df, prop_tab)

    signal_df = signal_df.cast({'label': pl.String}).rename({'label': 'cell_id'})

    return signal_df
