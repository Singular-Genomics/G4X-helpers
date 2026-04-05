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
from .workflow import DEFAULT_INPUT, collect_input, route_output

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput


DEFAULT_MASK_KEY = 'nuclei_exp'
LOGGER = logging.getLogger(__name__)


# @g4x_workflow
# region main function
def aggregate_cell_data(
    smp: 'G4Xoutput',
    segmentation_mask: str = DEFAULT_INPUT,
    mask_key: str | None = DEFAULT_MASK_KEY,
    *,
    out_dir: str = DEFAULT_INPUT,
    tx_table: str = DEFAULT_INPUT,
    overwrite: bool = False,
    gene_list: list[str] = DEFAULT_INPUT,
    protein_list: list[str] = DEFAULT_INPUT,
    show_progress: bool | None = None,
    logger: logging.Logger | None = None,
) -> None:

    log = logger or LOGGER
    seg_source = 'g4x-default' if segmentation_mask == DEFAULT_INPUT else 'custom'
    log.info('Running aggregate_cell_data with %s segmentation source', seg_source)

    # 1: Validate and collect input
    txtable_in = collect_input(smp, tx_table, TxTable, logger=log)
    segment_in = collect_input(smp, segmentation_mask, Segmentation, validate=False, logger=log)

    # 2: Import segmentation mask (keys will be validated here)
    mask = io.import_segmentation(segment_in.p, labels_key=mask_key, expected_shape=smp.shape)

    # 3: Validate and prepare output
    out_dir = smp.data_dir if out_dir == DEFAULT_INPUT else io.pathval.validate_dir_path(out_dir)

    log_with_path = partial(logut.log_with_path, logger=log, level='info')
    route_out = partial(route_output, smp, out_dir, overwrite=overwrite, logger=log)

    route_out(validator=CellMetadata)
    route_out(validator=CellxGene)
    route_out(validator=TxTable, overwrite=True)  # Always overwrite TxTable since we're updating the input table

    overwrite_segmentation = True if segmentation_mask == DEFAULT_INPUT else overwrite
    route_out(validator=Segmentation, overwrite=overwrite_segmentation)

    # 4: Create cell metadata
    cell_metadata = create_cell_metadata(
        smp,
        segmentation_mask=mask,
        seg_source=seg_source,
        show_progress=show_progress,
        logger=log,
    )

    # 5: Create cell x gene matrix
    cell_x_gene, tx_table_intersected = create_cell_x_gene(
        smp=smp,
        tx_table=txtable_in.load(lazy=False),
        segmentation_mask=mask,
        gene_labels=gene_list,
        logger=log,
    )

    if segmentation_mask == DEFAULT_INPUT:
        cell_metadata = add_nuclei_properties(smp, cell_metadata, show_progress=show_progress)
        tx_table_intersected = intersect_tx_with_cells(
            tx_table_intersected, smp.load_segmentation(expanded=False), column_name='in_nucleus'
        )

    # 6: Create cell x protein matrix (optional)
    if smp.src.pr_detected:
        route_out(validator=CellxProt)
        if protein_list != DEFAULT_INPUT:
            smp.set_proteins(protein_list)

        cell_x_protein = create_cell_x_signal(
            smp=smp,
            mask=mask,
            signal_list=smp.proteins,
            show_progress=show_progress,
            logger=log,
        )

        log_with_path(f'Writing {smp.out.CellxProt.name} table:', smp.out.CellxProt.p)
        cell_x_protein.sink_csv(smp.out.CellxProt.p, compression='gzip')

    # 7: Write the demuxed transcript table
    log_with_path(f'Writing {smp.out.TxTable.name} table:', smp.out.TxTable.p)
    tx_table_intersected.write_csv(smp.out.TxTable.p, compression='gzip')

    # 8: Write the cell x gene matrix
    log_with_path(f'Writing {smp.out.CellxGene.name} table:', smp.out.CellxGene.p)
    cell_x_gene.sink_csv(smp.out.CellxGene.p, compression='gzip')

    # 9: Write the cell metadata table
    log_with_path(f'Writing {smp.out.CellMetadata.name} table:', smp.out.CellMetadata.p)
    cell_metadata.sink_csv(smp.out.CellMetadata.p, compression='gzip')

    # 10: Write the segmentation mask (if using custom)
    if segmentation_mask != DEFAULT_INPUT:
        log_with_path(f'Writing {smp.out.Segmentation.name} mask:', smp.out.Segmentation.p)
        mask_data = {mask_key: mask}
        smp.out.Segmentation.main_key = mask_key
        np.savez(smp.out.Segmentation.p, **mask_data)

    log.info('Completed aggregating cell data')


# region high-level functions
def create_cell_metadata(
    smp: 'G4Xoutput',
    segmentation_mask: np.ndarray,
    cell_frame: pl.DataFrame | None = None,
    show_progress: bool | None = None,
    return_lazy: bool = True,
    seg_source: str = 'g4x-default',
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER
    log.debug('Creating cell metadata')

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
    log.info('Cell metadata created successfully')

    if return_lazy:
        return cell_meta
    return cell_meta.collect()


def create_cell_x_gene(
    smp: 'G4Xoutput',
    tx_table: pl.LazyFrame,
    segmentation_mask: str | np.ndarray = DEFAULT_INPUT,
    gene_labels: str | np.ndarray = DEFAULT_INPUT,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:

    log = logger or LOGGER
    log.debug('Creating cell x gene matrix')

    tx_table = intersect_tx_with_cells(tx_table, segmentation_mask)
    cell_frame = _cell_frame(segmentation_mask)

    if gene_labels == DEFAULT_INPUT:
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

    log.debug('Cell x gene matrix created successfully')

    if return_lazy:
        return cell_by_gene, tx_table
    return cell_by_gene.collect(), tx_table.collect()


def create_cell_x_signal(
    smp: 'G4Xoutput',
    mask: np.ndarray,
    signal_list: list[str],
    suffix: str = '_intensity_mean',
    show_progress: bool | None = None,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
    **kwargs,
) -> pl.LazyFrame:

    log = logger or LOGGER
    log.debug('Creating cell x signal matrix for %d signals.', len(signal_list))

    bead_mask = smp.load_bead_mask()
    bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
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
            signal_img, mask_flat=mask_flat, bead_mask_flat=bead_mask_flat, ch_label=ch_label, **kwargs
        )
        channel_df = channel_df.join(intensity_df, on=c.CELL_ID_NAME, how='left')

    channel_df = channel_df.sort(c.CELL_ID_NAME)
    log.debug('Cell x signal matrix created successfully')
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
    mask_props_nuc = mask_props_nuc.rename({c.CELL_AREA_NAME: 'nuclei_area_um'})

    cell_metadata = cell_metadata.drop([c.CELL_COORD_X, c.CELL_COORD_Y]).join(mask_props_nuc, on=c.CELL_ID_NAME)

    col_order = [
        c.CELL_ID_NAME,
        'sample_id',
        'tissue_type',
        'block',
        'seg_source',
        c.CELL_COORD_X,
        c.CELL_COORD_Y,
        'nuclei_area_um',
        c.CELL_AREA_NAME,
        'nuclearstain_intensity_mean',
        'cytoplasmicstain_intensity_mean',
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
    ch_label: str = 'img_intensity_mean',
    compute_backend: Literal['cpu', 'gpu', 'auto'] = 'auto',
) -> tuple[np.ndarray, np.ndarray]:

    backend = io.get_backend(which=compute_backend)

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
