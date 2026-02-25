from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from tqdm import tqdm


if TYPE_CHECKING:
    from g4x_helpers.g4x_output import G4Xoutput


def create_cell_x_protein(
    g4x_obj: 'G4Xoutput',
    mask: np.ndarray,
    signal_list: list[str] | None = None,
    cached: bool = False,
) -> pl.DataFrame | pl.LazyFrame:
    # if signal_list is None:
    #     signal_list = ['nuclear', 'eosin'] + g4x_obj.proteins

    print(f'Creating cell x protein matrix for {len(signal_list)} signals.')

    channel_name_map = {protein: protein for protein in signal_list}
    channel_name_map['nuclear'] = 'nuclearstain'
    channel_name_map['eosin'] = 'cytoplasmicstain'

    # TODO return here when bead masking is implemented
    # bead_mask = g4x_obj.load_bead_mask()
    # bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
    mask_flat = mask.ravel()

    cell_frame = g4x_obj.cell_frame

    for signal_name in tqdm(signal_list, desc='Extracting protein signal'):
        if signal_name not in ['nuclear', 'eosin']:
            image_type = 'protein'
            protein = signal_name
        else:
            image_type = signal_name
            protein = None

        signal_img = g4x_obj.load_image_by_type(image_type, thumbnail=False, protein=protein, cached=cached)

        ch_label = f'{channel_name_map[signal_name]}_intensity_mean'

        intensities = image_intensity_extraction(
            signal_img,
            mask_flat=mask_flat,
            bead_mask_flat=None,
        )

        cell_frame = cell_frame.with_columns(pl.Series(name=ch_label, values=intensities))

    return cell_frame


def image_intensity_extraction(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
) -> np.ndarray:
    img_flat = img.ravel()

    # Optional bead masking
    if bead_mask_flat is not None:
        bead_mask_flat = ~bead_mask_flat
        mask_flat = mask_flat[bead_mask_flat]
        img_flat = img_flat[bead_mask_flat]

    # Remove zeros and find unique labels
    mask_nonzero = mask_flat > 0
    labels = mask_flat[mask_nonzero]
    pixels = img_flat[mask_nonzero]

    _, inv = np.unique(labels, return_inverse=True)
    # `inv` now contains remapped labels 0..n_unique-1

    # Compute sums and counts using remapped labels
    sums = np.bincount(inv, weights=pixels)
    counts = np.bincount(inv)

    # Safe divide
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts != 0)

    return means
