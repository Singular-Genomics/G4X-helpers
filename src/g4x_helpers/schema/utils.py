import numpy as np
import polars as pl

from .. import io
from ..roi import Roi
from . import definition as sd


def detect_assay_type(smp_meta):
    tx_detected = smp_meta.get('transcript_panel', None) is not None
    pr_detected = smp_meta.get('protein_panel', None) is not None

    assay_type = (
        'combined'
        if tx_detected and pr_detected
        else 'tx_only'
        if tx_detected
        else 'pr_only'
        if pr_detected
        else 'undefined'
    )
    return assay_type


def migrate_image(
    migrator: sd.ImgDirectoryValidator,
    img_name: str,
    out_dir: str,
    img_type='auto',
    roi: Roi | None = None,
) -> None:
    io.convert.jp2_to_ometiff(
        in_file=migrator.mapped_files[img_name],
        out_file=f'{out_dir}/{img_name}.ome.tiff',
        img_type=img_type,
        create_thumb=True,
        report_size=False,
        extent=roi.extent_array if roi else None,
    )


def crop_tx_features(df: pl.LazyFrame, roi: Roi) -> pl.LazyFrame:
    df = df.filter(
        pl.col('x_pixel_coordinate').is_between(*roi.xlims, closed='left'),
        pl.col('y_pixel_coordinate').is_between(*roi.ylims, closed='left'),
    ).with_columns(
        pl.col('x_pixel_coordinate') - roi.xlims[0],
        pl.col('y_pixel_coordinate') - roi.ylims[0],
    )
    return df


def crop_bead_mask(bead_mask: sd.BeadMask, roi: Roi) -> dict[str, np.ndarray]:

    cropped = {}
    arr = bead_mask.load()
    cropped[bead_mask.DEFAULT_KEY] = roi.crop_array(arr)
    return cropped


def crop_segmentations(segmentation: sd.Segmentation, roi: Roi) -> dict[str, np.ndarray]:

    cleaned_masks = {}

    if segmentation.main_key == 'nuclei_exp':
        segmentation.main_key = 'nuclei'

    main_seg = segmentation.load()
    main_crop = roi.crop_array(main_seg)
    main_crop_cleaned = _remove_boundary_labels(main_crop)

    keep_labels = np.unique(main_crop_cleaned)
    if len(keep_labels) == 1 and keep_labels[0] == 0:
        raise ValueError('No valid segmentation labels remain after cropping.')

    cleaned_masks[segmentation.main_key] = main_crop_cleaned

    for key in segmentation.available_keys:
        if key == segmentation.main_key:
            continue

        sub_seg = segmentation.load(key=key)
        sub_crop = roi.crop_array(sub_seg)
        mask = ~np.isin(sub_crop, keep_labels)
        sub_crop_cleaned = sub_crop.copy()
        sub_crop_cleaned[mask] = 0

        sub_labels = np.unique(sub_crop_cleaned)
        if not np.isin(keep_labels, sub_labels).all():
            raise ValueError(
                f'Sub-segmentation "{key}" contains labels not present in main segmentation after cropping.'
            )
        cleaned_masks[key] = sub_crop_cleaned

    return cleaned_masks


def _remove_boundary_labels(labels, inplace=False):
    arr = labels if inplace else labels.copy()

    boundary_labels = np.unique(
        np.concatenate(
            [  # Collect all labels touching the boundary
                arr[0, :],  # top row
                arr[-1, :],  # bottom row
                arr[:, 0],  # left column
                arr[:, -1],  # right column
            ]
        )
    )
    # Create a mask for boundary labels and set them to 0
    mask = np.isin(arr, boundary_labels)
    arr[mask] = 0

    return arr
