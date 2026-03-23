import logging
import sys
from functools import lru_cache
from typing import TYPE_CHECKING, Literal

import numpy as np
import polars as pl
from skimage.measure import regionprops
from tqdm import tqdm

from .. import constants as c
from .. import io

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput


LOGGER = logging.getLogger(__name__)


# region metadata
def create_cell_metadata(
    g4x_obj: 'G4Xoutput',
    mask: np.ndarray | None,
    logger: logging.Logger | None = None,
    show_progress: bool | None = None,
):
    log = logger or LOGGER

    seg_source = 'g4x-default' if mask is None else 'custom'

    log.info('Creating cell metadata from %s segmentation source', seg_source)

    cell_meta = _cell_frame(g4x_obj)

    cell_meta = cell_meta.with_columns(
        pl.lit(g4x_obj.sample_id).alias('sample_id'),
        pl.lit(g4x_obj.tissue_type).alias('tissue_type'),
        pl.lit(g4x_obj.block).alias('block'),
    ).collect()

    # TODO this might break if someone changes the .npz file
    if mask is None:
        mask_props = extract_cell_props_g4x(g4x_obj, show_progress=show_progress)
    else:
        mask_props = extract_cell_props(mask, show_progress=show_progress)

    mask_props = mask_props.with_columns(pl.lit(seg_source).alias('seg_source'))

    if mask_props[c.CELL_ID_NAME].equals(cell_meta[c.CELL_ID_NAME]):
        cell_meta = cell_meta.join(mask_props, on=c.CELL_ID_NAME, how='left')
    else:
        raise ValueError('The CELL_ID columns in cell_meta and mask_props do not match.')

    log.info('Cell metadata created successfully')
    return cell_meta


def extract_cell_props(
    mask: np.ndarray,
    mask_name: str | None = None,
    show_progress: bool | None = None,
) -> pl.DataFrame:
    props = regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar

    if show_progress is None:
        show_progress = sys.stderr.isatty()

    prefix = f' {mask_name} ' if mask_name else ' '
    for prop in tqdm(props, desc=f'Extracting{prefix}mask properties', disable=not show_progress):
        label = prop.label  # The label (mask id)

        cell_y, cell_x = prop.centroid  # coordinate order: 'yx' (row, col)

        px_to_um_area = c.PIXEL_SIZE_MICRONS**2
        area_um = prop.area * px_to_um_area  # prop.area is in pixels

        prop_dict.append(
            {
                c.CELL_ID_NAME: label,
                'area_um': area_um,
                c.CELL_COORD_X: cell_x,
                c.CELL_COORD_Y: cell_y,
            }
        )
    schema = {
        c.CELL_ID_NAME: pl.Int32,
        'area_um': pl.Float32,
        c.CELL_COORD_X: pl.Float32,
        c.CELL_COORD_Y: pl.Float32,
    }
    mask_props = pl.DataFrame(prop_dict, schema=schema).sort(c.CELL_ID_NAME)
    return mask_props


def extract_cell_props_g4x(
    g4x_obj: 'G4Xoutput',
    show_progress: bool | None = None,
) -> pl.DataFrame:
    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.tree.Segmentation.path, expected_shape=g4x_obj.shape, labels_key='nuclei'
    )
    mask_props_nuc = extract_cell_props(mask=seg_mask, mask_name='nuclei', show_progress=show_progress)

    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.tree.Segmentation.path, expected_shape=g4x_obj.shape, labels_key='nuclei_exp'
    )
    mask_props_exp = extract_cell_props(mask=seg_mask, mask_name='nuclei_exp', show_progress=show_progress)
    del seg_mask

    mask_props_nuc = mask_props_nuc.rename({'area_um': 'nuclei_area_um'})
    mask_props_exp = mask_props_exp.rename({'area_um': 'wholecell_area_um'}).drop([c.CELL_COORD_X, c.CELL_COORD_Y])
    mask_props = mask_props_nuc.join(mask_props_exp, on=c.CELL_ID_NAME)

    mask_props = mask_props.select(
        [c.CELL_ID_NAME, c.CELL_COORD_X, c.CELL_COORD_Y, 'nuclei_area_um', 'wholecell_area_um']
    )

    return mask_props


# region transcripts
def create_cell_x_gene(
    g4x_obj: 'G4Xoutput',
    return_lazy: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    log = logger or LOGGER
    log.info('Creating cell x gene matrix')

    reads = pl.scan_csv(g4x_obj.tree.TranscriptTable.path)

    cell_by_gene = (
        reads.filter(pl.col(c.CELL_ID_NAME) != 0)
        .group_by(c.CELL_ID_NAME, 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index=c.CELL_ID_NAME, on_columns=g4x_obj.genes)
    )

    # Adding missing cells with zero counts
    cell_by_gene = _cell_frame(g4x_obj).join(cell_by_gene, left_on='seg_cell_id', right_on=c.CELL_ID_NAME, how='left')

    existing = cell_by_gene.collect_schema().names()
    if not set(existing[1:]) == set(g4x_obj.genes):
        raise ValueError('Mismatch between cell_by_gene columns and g4x_obj.genes')

    cell_by_gene = cell_by_gene.select(['seg_cell_id'] + g4x_obj.genes)
    cell_by_gene = cell_by_gene.fill_null(0)

    if return_lazy:
        return cell_by_gene

    log.info('Cell x gene matrix created successfully')
    return cell_by_gene.collect()


def intersect_tx_with_cells(
    tx_table: pl.DataFrame | pl.LazyFrame, mask: np.ndarray, column_name: str = c.CELL_ID_NAME
) -> pl.DataFrame:
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = tx_table.select(coord_order).collect().to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]
    tx_table = tx_table.with_columns(pl.lit(cell_ids).alias(column_name))
    return tx_table


def intersect_tx_with_cells_g4x(g4x_obj: 'G4Xoutput', tx_table: pl.DataFrame | None = None) -> pl.DataFrame:
    if tx_table is None:
        tx_table = g4x_obj.load_transcript_table(lazy=True)

    # Assign transcripts to segmentation labels
    mask = g4x_obj.load_segmentation(expanded=True)
    tx_table = intersect_tx_with_cells(tx_table, mask, column_name=c.CELL_ID_NAME)

    mask = g4x_obj.load_segmentation(expanded=False)
    tx_table = intersect_tx_with_cells(tx_table, mask, column_name='in_nucleus')
    del mask

    return tx_table


# region protein
def create_cell_x_protein(
    g4x_obj: 'G4Xoutput',
    mask: np.ndarray,
    signal_list: list[str] | None = None,
    suffix: str = '_intensity_mean',
    logger: logging.Logger | None = None,
    **kwargs,
) -> pl.LazyFrame:

    log = logger or LOGGER

    if signal_list is None:
        signal_list = g4x_obj.stains + g4x_obj.proteins

    log.info('Creating cell x protein matrix for %d signals.', len(signal_list))

    bead_mask = g4x_obj.load_bead_mask()
    bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
    mask_flat = mask.ravel()

    channel_df = _cell_frame(g4x_obj, lazy=True)

    for signal_name in tqdm(signal_list, desc='Extracting protein signal'):
        log.debug('Processing signal: %s', signal_name)
        if signal_name == c.NUCLEAR_STAIN:
            signal_img = g4x_obj.load_nuclear_image()
            signal_name += 'stain'
        elif signal_name == c.CYTOPLASMIC_STAIN:
            signal_img = g4x_obj.load_cytoplasmic_image()
            signal_name += 'stain'
        else:
            signal_img = g4x_obj.load_protein_image(protein=signal_name)

        ch_label = f'{signal_name}{suffix}'

        intensity_df = image_intensity_extraction(
            signal_img, mask_flat=mask_flat, bead_mask_flat=bead_mask_flat, ch_label=ch_label, **kwargs
        )
        channel_df = channel_df.join(intensity_df, on=c.CELL_ID_NAME, how='left')

    log.info('Cell x protein matrix created successfully')
    return channel_df


def image_intensity_extraction(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
    ch_label: str = 'img_intensity_mean',
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
) -> tuple[np.ndarray, np.ndarray]:

    backend = io.get_backend(which=compute_backend)

    if backend.use_gpu:
        unique_labels, means = _image_intensity_extraction_gpu(img, mask_flat, bead_mask_flat, backend=backend)
    else:
        unique_labels, means = _image_intensity_extraction_cpu(img, mask_flat, bead_mask_flat)

    result = pl.LazyFrame(
        {
            c.CELL_ID_NAME: unique_labels,
            ch_label: means,
        }
    )
    return result


def _image_intensity_extraction_cpu(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    img_flat = img.ravel()

    all_labels = np.unique(mask_flat[mask_flat > 0])
    means = np.full(all_labels.shape, np.nan, dtype=float)

    if bead_mask_flat is not None:
        keep = ~bead_mask_flat
        work_mask = mask_flat[keep]
        work_img = img_flat[keep]
    else:
        work_mask = mask_flat
        work_img = img_flat

    nz = work_mask > 0
    labels = work_mask[nz]
    pixels = work_img[nz]

    if labels.size == 0:
        return all_labels, means

    present_labels, inv = np.unique(labels, return_inverse=True)
    sums = np.bincount(inv, weights=pixels)
    counts = np.bincount(inv)
    present_means = sums / counts

    idx = np.searchsorted(all_labels, present_labels)
    means[idx] = present_means

    return all_labels, means


def _image_intensity_extraction_gpu(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
    gpu_backend='ComputeBackend',
) -> tuple[np.ndarray, np.ndarray]:

    cp = gpu_backend.cp

    img_flat = cp.asarray(img).ravel()
    mask_flat_gpu = cp.asarray(mask_flat)

    all_labels = cp.unique(mask_flat_gpu[mask_flat_gpu > 0])
    means = cp.full(all_labels.shape, cp.nan, dtype=cp.float32)

    if bead_mask_flat is not None:
        keep = ~cp.asarray(bead_mask_flat)
        work_mask = mask_flat_gpu[keep]
        work_img = img_flat[keep]
    else:
        work_mask = mask_flat_gpu
        work_img = img_flat

    nz = work_mask > 0
    labels = work_mask[nz]
    pixels = work_img[nz]

    if labels.size == 0:
        return cp.asnumpy(all_labels), cp.asnumpy(means)

    present_labels, inv = cp.unique(labels, return_inverse=True)
    sums = cp.bincount(inv, weights=pixels)
    counts = cp.bincount(inv)

    counts_safe = cp.where(counts != 0, counts, 1)
    present_means = sums / counts_safe

    idx = cp.searchsorted(all_labels, present_labels)
    means[idx] = present_means

    return cp.asnumpy(all_labels), cp.asnumpy(means)


# region utilities
@lru_cache(maxsize=1)
def _cell_frame(g4x_obj: 'G4Xoutput', lazy: bool = True) -> int:
    lf = pl.LazyFrame(g4x_obj.cell_labels, schema={c.CELL_ID_NAME: pl.Int32}).sort(c.CELL_ID_NAME)
    return lf if lazy else lf.collect()
