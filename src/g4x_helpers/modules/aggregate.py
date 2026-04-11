import logging
import sys
from functools import partial
from typing import TYPE_CHECKING, Literal

import numpy as np
import polars as pl
from skimage.measure import regionprops
from tqdm import tqdm

from .. import c, io
from .. import logging_utils as logut
from ..schema.definition import CellMetadata, CellxGene, CellxProt, Segmentation, TxTable
from .workflow import PRESET_SOURCE, collect_input, reroute_source

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput


DEFAULT_MASK_KEY = 'nuclei_exp'
LOGGER = logging.getLogger(__name__)


# region main function
# @g4x_workflow
def aggregate_cell_data(
    smp: 'G4Xoutput',
    segmentation_mask: str = PRESET_SOURCE,
    mask_key: str | None = DEFAULT_MASK_KEY,
    *,
    out_dir: str = PRESET_SOURCE,
    tx_table: str = PRESET_SOURCE,
    overwrite: bool = False,
    gene_list: list[str] = PRESET_SOURCE,
    protein_list: list[str] = PRESET_SOURCE,
    show_progress: bool | None = None,
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
    logger: logging.Logger | None = None,
) -> None:

    log = logger or LOGGER
    backend = io.get_backend(which=compute_backend)

    seg_source = 'g4x-default' if segmentation_mask == PRESET_SOURCE else 'custom'
    log.info('Using %s segmentation source', seg_source)

    # 1: Validate and collect input
    log.debug('Validating input and preparing output paths')
    txtable_in = collect_input(smp, tx_table, TxTable, logger=log)
    segment_in = collect_input(smp, segmentation_mask, Segmentation, validate=False, logger=log)

    # 2: Import segmentation mask (keys will be validated here)
    mask = io.import_segmentation(segment_in.p, labels_key=mask_key, expected_shape=smp.shape)

    # 3: Validate and prepare output
    out_dir = smp.data_dir if out_dir == PRESET_SOURCE else io.pathval.validate_dir_path(out_dir)

    log_with_path = partial(logut.log_with_path, logger=log, level='info')
    route_out = partial(reroute_source, smp, out_dir, overwrite=overwrite, logger=log)

    route_out(validator=CellMetadata)
    route_out(validator=CellxGene)

    overwrite_segmentation = True if segmentation_mask == PRESET_SOURCE else overwrite
    route_out(validator=Segmentation, overwrite=overwrite_segmentation)

    # 4: Create cell metadata
    log.info('Creating cell metadata')
    cell_metadata = create_cell_metadata(
        smp,
        segmentation_mask=mask,
        seg_source=seg_source,
        show_progress=show_progress,
        backend=backend,
        logger=log,
    )

    # 5: Create cell x gene matrix
    log.info('Creating cell x gene matrix')
    cell_x_gene, tx_table_intersected = create_cell_x_gene(
        smp=smp,
        tx_table=txtable_in.load(),
        segmentation_mask=mask,
        gene_labels=gene_list,
        logger=log,
    )

    if segmentation_mask == PRESET_SOURCE:
        log.debug('Adding nuclei properties to cell metadata and tx-table')
        cell_metadata = add_nuclei_properties(smp, cell_metadata, show_progress=show_progress)
        tx_table_intersected = intersect_tx_with_cells(
            tx_table_intersected, smp.load_segmentation(expanded=False), column_name='in_nucleus'
        )

    # 6: Create cell x protein matrix (optional)
    if smp.src.pr_detected:
        log.info('Creating cell x protein matrix')

        route_out(validator=CellxProt)

        if protein_list != PRESET_SOURCE:
            smp.set_proteins(protein_list)

        cell_x_protein = create_cell_x_signal(
            smp=smp,
            mask=mask,
            signal_list=smp.proteins,
            show_progress=show_progress,
            backend=backend,
            logger=log,
        )

        log_with_path(f'Writing {smp.out.CellxProt.name} table:', smp.out.CellxProt.p)
        cell_x_protein.sink_csv(smp.out.CellxProt.p, compression='gzip')

    # 7: Write the demuxed transcript table
    overwrite_txtable = True if tx_table == PRESET_SOURCE else overwrite
    route_out(validator=TxTable, overwrite=overwrite_txtable)
    log_with_path(f'Writing {smp.out.TxTable.name} table:', smp.out.TxTable.p)
    tx_table_intersected.write_csv(smp.out.TxTable.p, compression='gzip')

    # 8: Write the cell x gene matrix
    log_with_path(f'Writing {smp.out.CellxGene.name} table:', smp.out.CellxGene.p)
    cell_x_gene.sink_csv(smp.out.CellxGene.p, compression='gzip')

    # 9: Write the cell metadata table
    log_with_path(f'Writing {smp.out.CellMetadata.name} table:', smp.out.CellMetadata.p)
    cell_metadata.sink_csv(smp.out.CellMetadata.p, compression='gzip')

    # 10: Write the segmentation mask (if using custom)
    if segmentation_mask != PRESET_SOURCE:
        log_with_path(f'Writing {smp.out.Segmentation.name} mask:', smp.out.Segmentation.p)
        mask_data = {mask_key: mask}
        smp.out.Segmentation.main_key = mask_key
        np.savez(smp.out.Segmentation.p, **mask_data)


# region high-level functions
def create_cell_metadata(
    smp: 'G4Xoutput',
    segmentation_mask: np.ndarray,
    *,
    cell_frame: pl.DataFrame | None = None,
    show_progress: bool | None = None,
    return_lazy: bool = True,
    seg_source: str = 'g4x-default',
    backend: io.ComputeBackend = io.get_backend(which='auto'),
    logger: logging.Logger | None = None,
):
    # log = logger or LOGGER

    cell_frame = _cell_frame(segmentation_mask)

    cell_meta = cell_frame.with_columns(
        pl.lit(smp.sample_id).alias('sample_id'),
        pl.lit(smp.tissue_type).alias('tissue_type'),
        pl.lit(smp.block).alias('block'),
        pl.lit(seg_source).alias('seg_source'),
    )  # .collect()

    mask_props = extract_cell_props(segmentation_mask, show_progress=show_progress)

    stain_intensities = create_cell_x_signal(
        smp=smp,
        mask=segmentation_mask,
        signal_list=smp.stains,
        show_progress=show_progress,
        backend=backend,
        logger=logger,
    )

    del segmentation_mask

    mask_props = mask_props.join(stain_intensities, on=c.CELL_ID_NAME, how='left')

    mask_cells = mask_props.select(c.CELL_ID_NAME).collect()
    cell_meta_cells = cell_meta.select(c.CELL_ID_NAME).collect()
    if mask_cells.equals(cell_meta_cells):
        cell_meta = cell_meta.join(mask_props, on=c.CELL_ID_NAME, how='left')
    else:
        raise ValueError('The CELL_ID columns in cell_meta and mask_props do not match.')

    cell_meta = cell_meta.sort(c.CELL_ID_NAME)

    if return_lazy:
        return cell_meta
    return cell_meta.collect()


def create_cell_x_gene(
    smp: 'G4Xoutput',
    tx_table: pl.LazyFrame,
    segmentation_mask: str | np.ndarray = PRESET_SOURCE,
    gene_labels: str | np.ndarray = PRESET_SOURCE,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:

    # log = logger or LOGGER

    tx_table = intersect_tx_with_cells(tx_table, segmentation_mask)
    cell_frame = _cell_frame(segmentation_mask)

    if gene_labels == PRESET_SOURCE:
        gene_labels = smp.genes

    cell_by_gene = (
        tx_table.filter(pl.col(c.CELL_ID_NAME) != 0)
        .group_by(c.CELL_ID_NAME, c.GENE_ID_NAME)
        .agg(pl.len().alias('counts'))
        .sort(c.GENE_ID_NAME)
        .pivot(on=c.GENE_ID_NAME, values='counts', index=c.CELL_ID_NAME, on_columns=gene_labels)
    )

    # Adding missing cells with zero counts
    cell_by_gene = cell_frame.join(cell_by_gene.lazy(), on=c.CELL_ID_NAME, how='left')

    existing = cell_by_gene.collect_schema().names()
    if not set(existing[1:]) == set(gene_labels):
        raise ValueError('Mismatch between cell_by_gene columns and gene_labels')

    cell_by_gene = cell_by_gene.select([c.CELL_ID_NAME] + gene_labels)
    cell_by_gene = cell_by_gene.fill_null(0).sort(c.CELL_ID_NAME)

    if return_lazy:
        return cell_by_gene, tx_table
    return cell_by_gene.collect(), tx_table.collect()


def create_cell_x_signal(
    smp: 'G4Xoutput',
    mask: np.ndarray,
    *,
    signal_list: list[str],
    suffix: str = c.IMG_INTENSITY_HANDLE,
    show_progress: bool | None = None,
    return_lazy: bool = True,
    backend: io.ComputeBackend = io.get_backend(which='auto'),
    logger: logging.Logger | None = None,
    **kwargs,
) -> pl.LazyFrame:

    log = logger or LOGGER
    log.debug('Intersecting cells with signals: %s', signal_list)

    bead_mask = smp.load_bead_mask()
    # bead_mask = _add_artifically_large_beads(bead_mask)
    if bead_mask is not None:
        bead_mask_flat = bead_mask.ravel()
    else:
        log.warning('Bead mask not found. Proceeding without excluding beads from signal extraction.')
        bead_mask_flat = None

    mask_flat = mask.ravel()

    channel_df = _cell_frame(mask)

    if show_progress is None:
        show_progress = sys.stderr.isatty()

    for signal_name in tqdm(signal_list, desc='Extracting image signals', disable=not show_progress):
        log.debug('Processing signal: %s', signal_name)
        if signal_name == c.NUCLEAR_STAIN:
            signal_img = smp.load_nuclear_image()
            signal_name += 'stain'
        elif signal_name == c.CYTOPLASMIC_STAIN:
            signal_img = smp.load_cytoplasmic_image()
            signal_name += 'stain'
        else:
            signal_img = smp.load_protein_image(protein=signal_name)

        ch_label = f'{signal_name}{suffix}'

        intensity_df = image_intensity_extraction(
            signal_img,
            mask_flat=mask_flat,
            bead_mask_flat=bead_mask_flat,
            ch_label=ch_label,
            backend=backend,
            **kwargs,
        )
        channel_df = channel_df.join(intensity_df, on=c.CELL_ID_NAME, how='left')

    channel_df = channel_df.sort(c.CELL_ID_NAME)
    if return_lazy:
        return channel_df
    return channel_df.collect()


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
    return mask_props


def add_nuclei_properties(smp, cell_metadata, show_progress=True):
    segmentation_mask = smp.load_segmentation(expanded=False)
    mask_props_nuc = extract_cell_props(mask=segmentation_mask, mask_name='nuclei', show_progress=show_progress)
    mask_props_nuc = mask_props_nuc.rename({c.CELL_AREA_NAME: c.NUC_AREA_NAME})

    cell_metadata = cell_metadata.drop([c.CELL_COORD_X, c.CELL_COORD_Y]).join(mask_props_nuc, on=c.CELL_ID_NAME)

    col_order = [
        c.CELL_ID_NAME,
        'sample_id',
        'tissue_type',
        'block',
        'seg_source',
        c.CELL_COORD_X,
        c.CELL_COORD_Y,
        c.NUC_AREA_NAME,
        c.CELL_AREA_NAME,
        c.NUC_STAIN_INTENSITY,
        c.CYT_STAIN_INTENSITY,
    ]
    return cell_metadata.select(col_order)


def intersect_tx_with_cells(
    tx_table: pl.DataFrame | pl.LazyFrame, mask: np.ndarray, column_name: str = c.CELL_ID_NAME
) -> pl.DataFrame:
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = tx_table.select(coord_order).to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]
    tx_table = tx_table.with_columns(pl.lit(cell_ids).alias(column_name))
    return tx_table


def image_intensity_extraction(
    img: np.ndarray,
    mask_flat: np.ndarray,
    bead_mask_flat: np.ndarray | None = None,
    ch_label: str = 'img_mean',
    backend: io.ComputeBackend = io.get_backend(which='auto'),
) -> tuple[np.ndarray, np.ndarray]:

    if backend.use_gpu:
        unique_labels, means = _image_intensity_extraction_gpu(img, mask_flat, bead_mask_flat, gpu_backend=backend)
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


def _cell_frame(segmentation_mask: np.ndarray, lazy: bool = True) -> int:
    cell_labels = np.unique(segmentation_mask)
    cell_labels = cell_labels[cell_labels != 0]
    lf = pl.LazyFrame(cell_labels, schema={c.CELL_ID_NAME: pl.Int32}).sort(c.CELL_ID_NAME)
    return lf if lazy else lf.collect()


def _add_artifically_large_beads(beads, sq_size=500):
    half = sq_size // 2
    rows, cols = beads.shape
    center_row, center_col = rows // 2, cols // 2

    # set center block to True
    beads[center_row - half : center_row + half, center_col - half : center_col + half] = True
    return beads
