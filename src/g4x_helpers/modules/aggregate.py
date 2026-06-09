import logging
import sys

import numpy as np
import polars as pl
from skimage.measure import regionprops
from tqdm import tqdm

from .. import constants as c
from .. import io

log = logging.getLogger(__name__)


# region main functions
def cell_metadata(
    cell_mask: np.ndarray,
    nuclei_labels: np.ndarray | None = None,
    static_columns: dict = {'sample_id': 'not-specified'},
    *,
    show_progress: bool = False,
) -> pl.LazyFrame:

    cell_meta = extract_cell_props(cell_mask, show_progress=show_progress)

    if nuclei_labels is not None:
        cell_meta = cell_meta.drop([c.CELL_COORD_X, c.CELL_COORD_Y])
        nuc_meta = extract_cell_props(nuclei_labels, show_progress=show_progress)
        nuc_meta = nuc_meta.rename({c.CELL_AREA_NAME: c.NUC_AREA_NAME})
        cell_meta = cell_meta.join(nuc_meta, on=c.CELL_ID_NAME, how='left')

    cell_meta = cell_meta.with_columns([pl.lit(value).alias(key) for key, value in static_columns.items()])

    col_order = (
        [c.CELL_ID_NAME]
        + list(static_columns.keys())
        + [c.CELL_COORD_X, c.CELL_COORD_Y]
        + [col for col in cell_meta.collect_schema().names() if col.endswith('_area_um')]
    )

    cell_meta = cell_meta.select(col_order).sort(c.CELL_ID_NAME)

    return cell_meta


def cell_x_gene(
    tx_table: pl.LazyFrame,
    cell_mask: np.ndarray,
    included_cells: list | None = None,
    included_genes: list | None = None,
    return_tx_table: bool = False,
) -> pl.LazyFrame:

    tx_table = intersect_cells_with_tx(tx_table, cell_mask)

    existing_gene_ids = tx_table.select(c.GENE_ID_NAME).unique().sort(c.GENE_ID_NAME).collect().to_series().to_list()

    if included_genes is not None:
        _report_comparison(requested=included_genes, existing=existing_gene_ids, data_type='genes')
    else:
        included_genes = existing_gene_ids

    cell_by_gene = (
        tx_table.filter(pl.col(c.CELL_ID_NAME) != 0)
        .group_by(c.CELL_ID_NAME, c.GENE_ID_NAME)
        .agg(pl.len().alias('counts'))
        .sort(c.GENE_ID_NAME)
        .pivot(on=c.GENE_ID_NAME, values='counts', index=c.CELL_ID_NAME, on_columns=included_genes)
    )

    # Adding missing cells with zero counts
    if included_cells is not None:
        existing_cell_ids = (
            cell_by_gene.select(c.CELL_ID_NAME).unique().sort(c.CELL_ID_NAME).collect().to_series().to_list()
        )

        _report_comparison(requested=included_cells, existing=existing_cell_ids, data_type='cells')

        all_cells = pl.LazyFrame(included_cells, schema={c.CELL_ID_NAME: pl.UInt32})
        cell_by_gene = all_cells.join(cell_by_gene, on=c.CELL_ID_NAME, how='left')

    # ensure final table order mathes the gene order in the input table
    cell_by_gene = cell_by_gene.select([c.CELL_ID_NAME] + included_genes)

    # fill missing values with zeros (i.e. genes not detected in a cell)
    cell_by_gene = cell_by_gene.fill_null(0).sort(c.CELL_ID_NAME)

    if return_tx_table:
        return cell_by_gene, tx_table

    return cell_by_gene


def cell_x_signal(
    images: dict,
    cell_mask: np.ndarray,
    *,
    included_cells: list[int] | None = None,
    bead_mask: np.ndarray | None = None,
    suffix: str = c.IMG_INTENSITY_HANDLE,
    show_progress: bool = False,
    backend: io.ComputeBackend = io.get_backend(which='auto'),
) -> pl.LazyFrame:

    signal_list = list(images.keys())
    log.debug('Intersecting cells with signals: %s', signal_list)

    if bead_mask is not None:
        bead_mask_flat = bead_mask.ravel()
    else:
        log.warning('Bead mask not provided. Proceeding without excluding beads from signal extraction.')
        bead_mask_flat = None

    mask_flat = cell_mask.ravel()

    for signal_name in tqdm(signal_list, desc='Extracting image signals', disable=not show_progress):
        log.debug('Processing signal: %s', signal_name)
        signal_img = io.import_image(images[signal_name])

        ch_label = f'{signal_name}{suffix}'

        intensity_df = intersect_cells_with_img(
            signal_img, mask_flat=mask_flat, bead_mask_flat=bead_mask_flat, ch_label=ch_label, backend=backend
        )
        if signal_name == signal_list[0]:
            channel_df = intensity_df
        else:
            channel_df = channel_df.join(intensity_df, on=c.CELL_ID_NAME, how='left')

    if included_cells is not None:
        existing_cell_ids = (
            channel_df.select(c.CELL_ID_NAME).unique().sort(c.CELL_ID_NAME).collect().to_series().to_list()
        )
        _report_comparison(requested=included_cells, existing=existing_cell_ids, data_type='cells')

        all_cells = pl.LazyFrame(included_cells, schema={c.CELL_ID_NAME: pl.UInt32})
        channel_df = all_cells.join(channel_df, on=c.CELL_ID_NAME, how='left')

    channel_df = channel_df.sort(c.CELL_ID_NAME)
    return channel_df


# region core functions
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
                c.CELL_COORD_X: cell_x,
                c.CELL_COORD_Y: cell_y,
                c.CELL_AREA_NAME: area_um,
            }
        )
    schema = {
        c.CELL_ID_NAME: pl.Int32,
        c.CELL_COORD_X: pl.Float32,
        c.CELL_COORD_Y: pl.Float32,
        c.CELL_AREA_NAME: pl.Float32,
    }
    mask_props = pl.LazyFrame(prop_dict, schema=schema).sort(c.CELL_ID_NAME)
    if mask_name:
        mask_props = mask_props.rename(
            {col: col.replace('cell', mask_name) for col in schema.keys() if col.startswith('cell_')}
        )
    return mask_props.sort(c.CELL_ID_NAME)


def intersect_cells_with_tx(
    tx_table: pl.LazyFrame, mask: np.ndarray, column_name: str = c.CELL_ID_NAME
) -> pl.LazyFrame:
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = tx_table.select(coord_order).collect().to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]
    tx_table = tx_table.with_columns(pl.lit(cell_ids).alias(column_name))
    return tx_table


def intersect_cells_with_img(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
    ch_label: str = 'img_mean',
    backend: io.ComputeBackend = io.get_backend(which='auto'),
) -> pl.LazyFrame:

    if backend.use_gpu:
        unique_labels, means = _image_intensity_extraction_gpu(
            img=img, mask_flat=mask_flat, bead_mask_flat=bead_mask_flat, gpu_backend=backend
        )
    else:
        unique_labels, means = _image_intensity_extraction_cpu(
            img=img, mask_flat=mask_flat, bead_mask_flat=bead_mask_flat
        )

    result = pl.LazyFrame(
        {
            c.CELL_ID_NAME: unique_labels,
            ch_label: means,
        }
    )
    return result


# region private functions
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
    batch_size_pixels: int = 8_000_000,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-label mean image intensity using bounded GPU batches.
    """
    img_flat = np.asarray(img).ravel()
    mask_flat = np.asarray(mask_flat).ravel()

    if img_flat.size != mask_flat.size:
        raise ValueError(
            f'img and mask_flat must contain the same number of pixels: {img_flat.size=} != {mask_flat.size=}'
        )

    if bead_mask_flat is not None:
        bead_mask_flat = np.asarray(bead_mask_flat).ravel()

        if bead_mask_flat.size != mask_flat.size:
            raise ValueError(
                f'bead_mask_flat and mask_flat must contain the same number of pixels: '
                f'{bead_mask_flat.size=} != {mask_flat.size=}'
            )

    if batch_size_pixels <= 0:
        raise ValueError(f'batch_size_pixels must be positive: {batch_size_pixels=}')

    all_labels = np.unique(mask_flat[mask_flat > 0])
    means = np.full(all_labels.shape, np.nan, dtype=np.float32)

    if all_labels.size == 0:
        return all_labels, means

    sums_total = np.zeros(all_labels.shape, dtype=np.float64)
    counts_total = np.zeros(all_labels.shape, dtype=np.int64)

    n_pixels = mask_flat.size

    for start in range(0, n_pixels, batch_size_pixels):
        stop = min(start + batch_size_pixels, n_pixels)

        img_flat_batch = img_flat[start:stop]
        mask_flat_batch = mask_flat[start:stop]

        if bead_mask_flat is not None:
            bead_mask_flat_batch = bead_mask_flat[start:stop]
        else:
            bead_mask_flat_batch = None

        present_labels, sums, counts = _image_intensity_extraction_gpu_batch(
            img_flat_batch=img_flat_batch,
            mask_flat_batch=mask_flat_batch,
            bead_mask_flat_batch=bead_mask_flat_batch,
            gpu_backend=gpu_backend,
        )

        if present_labels.size == 0:
            continue

        idx = np.searchsorted(all_labels, present_labels)

        sums_total[idx] += sums
        counts_total[idx] += counts

    valid = counts_total > 0
    means[valid] = (sums_total[valid] / counts_total[valid]).astype(np.float32)

    return all_labels, means


def _image_intensity_extraction_gpu_batch(
    img_flat_batch: np.ndarray,
    mask_flat_batch: np.ndarray,
    bead_mask_flat_batch: np.ndarray | None = None,
    gpu_backend='ComputeBackend',
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute partial per-label sums and counts for one pixel batch on GPU.

    Returns
    -------
    present_labels : np.ndarray
        Labels present in this batch after bead exclusion.
    sums : np.ndarray
        Sum of image intensities for each present label.
    counts : np.ndarray
        Pixel count for each present label.
    """

    cp = gpu_backend.cp

    img_flat_gpu = cp.asarray(img_flat_batch)
    mask_flat_gpu = cp.asarray(mask_flat_batch)

    if bead_mask_flat_batch is not None:
        bead_mask_flat_gpu = cp.asarray(bead_mask_flat_batch)
        keep = ~bead_mask_flat_gpu
        work_mask = mask_flat_gpu[keep]
        work_img = img_flat_gpu[keep]
    else:
        work_mask = mask_flat_gpu
        work_img = img_flat_gpu

    nz = work_mask > 0

    if not bool(cp.any(nz)):
        return (
            np.empty(0, dtype=mask_flat_batch.dtype),
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.int64),
        )

    labels = work_mask[nz]
    pixels = work_img[nz]

    present_labels, inv = cp.unique(labels, return_inverse=True)
    sums = cp.bincount(inv, weights=pixels.astype(cp.float64))
    counts = cp.bincount(inv)

    return (
        cp.asnumpy(present_labels),
        cp.asnumpy(sums).astype(np.float64, copy=False),
        cp.asnumpy(counts).astype(np.int64, copy=False),
    )


def _report_comparison(requested: list, existing: list, data_type: str = 'genes'):
    missing_in_reference = set(existing) - set(requested)
    missing_in_data = set(requested) - set(existing)

    if missing_in_data:
        log.warning(f'{len(missing_in_data)} {data_type} are missing in the data. They will be filled with zeros.')
    if missing_in_reference:
        log.warning(
            f'{len(missing_in_reference)} {data_type} are in the data but were not in the reference. They will be missing from the final output.'
        )


def _add_artifically_large_beads(beads, sq_size=500):
    """ "
    For testing purposes, creates a large "bead" that covers cells entirely.
    In those cases we expect cells with NaN in the output
    """

    half = sq_size // 2
    rows, cols = beads.shape
    center_row, center_col = rows // 2, cols // 2

    # set center block to True
    beads[center_row - half : center_row + half, center_col - half : center_col + half] = True
    return beads
