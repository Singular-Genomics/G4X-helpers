from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

import glymur
import imageio.v3 as iio
import numpy as np
import shapely.affinity
import skimage.measure
import tifffile as tiff
from shapely.geometry import Polygon
from skimage import exposure, transform
from skimage.measure._regionprops import RegionProperties
from tqdm import tqdm

from .. import constants as c

if TYPE_CHECKING:
    from geopandas.geodataframe import GeoDataFrame


img_settings = {
    'rgb': {'axes': 'YXC', 'tn_clip': None},
    'grey': {'axes': 'YX', 'tn_clip': 0.99},
}


def gdf_to_ndarray(gdf: 'GeoDataFrame', target_shape: tuple) -> np.ndarray:
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


def ndarray_to_gdf(mask: np.ndarray, nudge: bool = True) -> 'GeoDataFrame':
    from geopandas import GeoDataFrame

    if mask.max() == 0:
        return GeoDataFrame(geometry=[])

    def _region_props_to_polygons(region_props: RegionProperties) -> list[Polygon]:
        mask = np.pad(region_props.image, 1)
        contours = skimage.measure.find_contours(mask, 0.5)

        # shapes with <= 3 vertices, i.e. lines, can't be converted into a polygon
        polygons = [Polygon(contour[:, [1, 0]]) for contour in contours if contour.shape[0] >= 4]

        yoff, xoff, *_ = region_props.bbox
        return [shapely.affinity.translate(poly, xoff, yoff) for poly in polygons]

    regions = skimage.measure.regionprops(mask)

    geoms = []
    labels = []
    # Wrap the iteration in tqdm to show progress
    for region in tqdm(regions, desc='Vectorizing regions'):
        polys = _region_props_to_polygons(region)
        geoms.extend(polys)
        # add the region label once per polygon
        labels.extend([region.label] * len(polys))

    gdf = GeoDataFrame({c.CELL_ID_NAME: labels}, geometry=geoms)

    if nudge:
        gdf['geometry'] = gdf['geometry'].translate(xoff=-0.5, yoff=-0.5)

    return gdf


def jp2_to_ometiff(
    in_file: str | Path,
    out_file: str | Path,
    img_type: Literal['rgb', 'grey', 'auto'] = 'auto',
    create_thumb: bool = True,
    report_size: bool = True,
):
    in_file = Path(in_file)
    out_file = Path(out_file)

    if create_thumb:
        thumb_dir = in_file.parent / 'thumbs'
        thumb_suffix = '_thumbnail'
        thumb_dir.mkdir(exist_ok=True)

    if not in_file.exists():
        raise FileNotFoundError(f"Input file '{in_file}' does not exist.")
    if in_file.suffix != '.jp2':
        raise ValueError(f"Input file '{in_file}' is not a JP2 file.")

    img = glymur.Jp2k(in_file)[:]

    # --- auto detect ---
    if img_type == 'auto':
        if img.ndim == 2:
            img_type = 'grey'
        elif img.ndim == 3 and img.shape[2] == 1:
            img_type = 'grey'
            img = img[..., 0]  # squeeze channel
        elif img.ndim == 3 and img.shape[2] == 3:
            img_type = 'rgb'
        else:
            raise ValueError(f'Cannot infer image type from shape {img.shape}')

    if img_type not in img_settings:
        raise ValueError(f"Invalid img_type '{img_type}'. Expected one of {list(img_settings.keys())}.")

    settings = img_settings[img_type]

    if img_type == 'rgb' and (img.ndim != 3 or img.shape[2] != 3):
        raise ValueError(f'Expected RGB image with shape (H, W, 3), got {img.shape}.')

    # Write OME-TIFF
    tiff.imwrite(
        out_file,
        img,
        ome=True,
        bigtiff=True,
        compression='zstd',
        metadata={
            'axes': settings['axes'],
            'PhysicalSizeX': c.PIXEL_SIZE_MICRONS,
            'PhysicalSizeY': c.PIXEL_SIZE_MICRONS,
            'PhysicalSizeXUnit': 'µm',
            'PhysicalSizeYUnit': 'µm',
        },
    )

    if create_thumb:
        q = settings['tn_clip']
        img_thumb = create_img_thumbnail(img, downsample_ratio=0.2, clip_q=q, anti_aliasing=False)

        fname = out_file.name.split('.')[0]
        out_file_thumb = out_file.with_stem(f'{fname}{thumb_suffix}').with_suffix('.png')
        out_file_thumb = thumb_dir / out_file_thumb.name
        iio.imwrite(out_file_thumb, img_thumb)

    if report_size:
        file_size_i = (in_file).stat().st_size
        file_size_o = (out_file).stat().st_size
        delta = file_size_o - file_size_i

        i = _human_readable_size(file_size_i)
        o = _human_readable_size(file_size_o)

        size_pct = (file_size_o / file_size_i * 100) - 100

        sign = '+' if delta > 0 else '-'
        print(f'Size difference: {sign}{abs(size_pct):.2f}% | jp2={i}, tiff={o}')

        if create_thumb:
            print(f'thumbnail size: {_human_readable_size((out_file_thumb).stat().st_size)}')


def create_img_thumbnail(
    img: np.ndarray,
    downsample_ratio: float = 0.25,
    clip_q: float | None = None,
    anti_aliasing: bool = False,
    dtype: type = np.uint8,
):
    """Create a smaller version of the input image with optional intensity clipping."""

    new_shape = (np.array(img.shape) * downsample_ratio).astype(int)
    new_shape = (new_shape[0], new_shape[1])
    img_resize = transform.resize(img, new_shape, anti_aliasing=anti_aliasing, preserve_range=False)

    if clip_q is not None:
        p_clip = np.quantile(img_resize, clip_q)
        in_range = (img_resize.min(), p_clip)
    else:
        in_range = 'image'

    img_rescale = exposure.rescale_intensity(img_resize, in_range=in_range, out_range=dtype)
    return img_rescale


def _human_readable_size(num_bytes: int, decimals: int = 2) -> str:
    """
    Convert a file size in bytes to a human-readable string.
    """
    units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

    size = float(num_bytes)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f'{size:.{decimals}f} {unit}'
        size /= 1024
