from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import shapely.affinity
import skimage.measure
from geopandas import GeoDataFrame
from shapely.geometry import Polygon
from skimage.measure._regionprops import RegionProperties
from tqdm import tqdm

if TYPE_CHECKING:
    from geopandas.geodataframe import GeoDataFrame


def rasterize_polygons(gdf: GeoDataFrame, target_shape: tuple) -> np.ndarray:
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


def vectorize_mask(mask: np.ndarray, nudge: bool = True) -> GeoDataFrame:
    if mask.max() == 0:
        return GeoDataFrame(geometry=[])

    regions = skimage.measure.regionprops(mask)

    geoms = []
    labels = []
    # Wrap the iteration in tqdm to show progress
    for region in tqdm(regions, desc='Vectorizing regions'):
        polys = region_props_to_polygons(region)
        geoms.extend(polys)
        # add the region label once per polygon
        labels.extend([region.label] * len(polys))

    gdf = GeoDataFrame({'label': labels}, geometry=geoms)

    if nudge:
        # GeoSeries.translate works elementwise
        gdf['geometry'] = gdf['geometry'].translate(xoff=-0.5, yoff=-0.5)

    return gdf


def region_props_to_polygons(region_props: RegionProperties) -> list[Polygon]:
    mask = np.pad(region_props.image, 1)
    contours = skimage.measure.find_contours(mask, 0.5)

    # shapes with <= 3 vertices, i.e. lines, can't be converted into a polygon
    polygons = [Polygon(contour[:, [1, 0]]) for contour in contours if contour.shape[0] >= 4]

    yoff, xoff, *_ = region_props.bbox
    return [shapely.affinity.translate(poly, xoff, yoff) for poly in polygons]
