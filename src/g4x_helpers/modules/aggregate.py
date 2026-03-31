import logging
import sys
from typing import TYPE_CHECKING, Literal

import numpy as np
import polars as pl
from skimage.measure import regionprops
from tqdm import tqdm

from .. import c, io
from ..logging_utils import PGAP

# from .workflow import g4x_workflow

if TYPE_CHECKING:
    from ..g4x_output import G4Xoutput

DEFAULT_INPUT = '__g4x_default__'
LOGGER = logging.getLogger(__name__)


# @g4x_workflow
def aggregate_cell_data(
    g4x_obj: 'G4Xoutput',
    segmentation_mask: str = DEFAULT_INPUT,
    out_dir: str = DEFAULT_INPUT,
    *,
    tx_table: str = DEFAULT_INPUT,
    gene_list: list[str] = DEFAULT_INPUT,
    protein_list: list[str] = DEFAULT_INPUT,
    mask_key: str | None = None,
    override: bool = False,
    show_progress: bool | None = None,
    logger: logging.Logger | None = None,
) -> None:

    log = logger or LOGGER
    log.info('Aggregating cell data')

    if out_dir == DEFAULT_INPUT:
        out_dir = g4x_obj.data_dir
    else:
        out_dir = io.pathval.validate_dir_path(out_dir)

    log.info('Output directory data:\n%s%s', PGAP, out_dir)

    tx_table_out = out_dir / g4x_obj.tree.TranscriptTable.p
    metadata_out = out_dir / g4x_obj.tree.CellMetadata.p
    cxg_out = out_dir / g4x_obj.tree.CellxGene.p
    cxp_out = out_dir / g4x_obj.tree.CellxProtein.p

    for path in [tx_table_out, metadata_out, cxg_out, cxp_out]:
        if path.exists() and not override:
            raise FileExistsError(f'File {path} already exists and override is False.')
        io.pathval.ensure_parent_dir(path)

    protein_list = g4x_obj.proteins if protein_list == DEFAULT_INPUT else protein_list
    if protein_list:
        unavailable_proteins = [p for p in protein_list if p not in g4x_obj.proteins]
        if unavailable_proteins:
            raise ValueError(f'The following requested proteins are not available in this data: {unavailable_proteins}')

    if tx_table == DEFAULT_INPUT:
        log.debug('Using default transcript table from G4X-output')
        tx_table_path = g4x_obj.tree.TranscriptTable.p
    else:
        tx_table_path = io.pathval.validate_file_path(tx_table)
        log.debug('Using transcript table from:\n%s%s', PGAP, tx_table_path)

    if segmentation_mask == DEFAULT_INPUT:
        log.debug('Using default segmentation mask from G4X-output')
        mask = None
        cell_frame = _cell_frame(segmentation_mask=g4x_obj.load_segmentation(expanded=False))
    else:
        segmentation_mask_path = io.pathval.validate_file_path(segmentation_mask)
        log.info('Importing segmentation mask from:\n%s%s', PGAP, segmentation_mask_path)
        mask = io.import_segmentation(
            segmentation_mask_path,
            expected_shape=g4x_obj.shape,
            labels_key=mask_key,
            logger=log,
        )
        cell_frame = _cell_frame(segmentation_mask=mask)

    # log.info('Creating cell metadata')
    cell_metadata = create_cell_metadata(
        g4x_obj,
        cell_frame=cell_frame,
        segmentation_mask=mask,
        show_progress=show_progress,
        logger=log,
    )

    # log.info('Creating cell x gene matrix')
    tx_table_df = pl.scan_csv(tx_table_path)
    cell_x_gene, tx_table_intersected = create_cell_x_gene(
        g4x_obj=g4x_obj,
        tx_table=tx_table_df,
        segmentation_mask=segmentation_mask,
        gene_labels=gene_list,
        logger=log,
    )

    # log.info('Creating cell x protein matrix')

    if protein_list:
        cell_mask = g4x_obj.load_segmentation(expanded=True) if segmentation_mask == DEFAULT_INPUT else mask
        cell_x_protein = create_cell_x_signal(
            g4x_obj=g4x_obj,
            mask=cell_mask,
            signal_list=protein_list,
            logger=log,
        )

        log.debug('Writing cell x protein matrix to CSV:\n%s%s', PGAP, cxp_out)
        cell_x_protein.sink_csv(cxp_out, compression='gzip')

    log.debug('Writing transcript table to CSV:\n%s%s', PGAP, tx_table_out)
    if out_dir == g4x_obj.data_dir:
        tx_table_intersected.collect().write_csv(tx_table_out, compression='gzip')
    else:
        tx_table_intersected.sink_csv(tx_table_out, compression='gzip')

    log.debug('Writing cell metadata to CSV:\n%s%s', PGAP, metadata_out)
    cell_metadata.sink_csv(metadata_out, compression='gzip')
    log.debug('Writing cell x gene matrix to CSV:\n%s%s', PGAP, cxg_out)
    cell_x_gene.sink_csv(cxg_out, compression='gzip')

    log.info('Completed aggregating cell data')


# region metadata
def create_cell_metadata(
    g4x_obj: 'G4Xoutput',
    segmentation_mask: np.ndarray | None,
    cell_frame: pl.DataFrame | None = None,
    show_progress: bool | None = None,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
):
    log = logger or LOGGER

    seg_source = 'g4x-default' if segmentation_mask is None else 'custom'

    log.info('Creating cell metadata from %s segmentation source', seg_source)

    if cell_frame is None:
        cell_frame = _cell_frame(g4x_obj.cell_labels)

    cell_meta = cell_frame.with_columns(
        pl.lit(g4x_obj.sample_id).alias('sample_id'),
        pl.lit(g4x_obj.tissue_type).alias('tissue_type'),
        pl.lit(g4x_obj.block).alias('block'),
        pl.lit(seg_source).alias('seg_source'),
    )  # .collect()

    # TODO this might break if someone changes the .npz file
    if segmentation_mask is None:
        segmentation_mask = g4x_obj.load_segmentation(expanded=False)
        mask_props_nuc = extract_cell_props(mask=segmentation_mask, mask_name='nuclei', show_progress=show_progress)

        segmentation_mask = g4x_obj.load_segmentation(expanded=True)
        mask_props_exp = extract_cell_props(mask=segmentation_mask, mask_name='nuclei_exp', show_progress=show_progress)

        mask_props_nuc = mask_props_nuc.rename({'area_um': 'nuclei_area_um'})
        mask_props_exp = mask_props_exp.rename({'area_um': 'wholecell_area_um'}).drop([c.CELL_COORD_X, c.CELL_COORD_Y])
        mask_props = mask_props_nuc.join(mask_props_exp, on=c.CELL_ID_NAME)

        mask_props = mask_props.select(
            [c.CELL_ID_NAME, c.CELL_COORD_X, c.CELL_COORD_Y, 'nuclei_area_um', 'wholecell_area_um']
        )

        stain_intensities = create_cell_x_signal(
            g4x_obj=g4x_obj,
            mask=segmentation_mask,
            signal_list=g4x_obj.stains,
            logger=logger,
        )

    else:
        mask_props = extract_cell_props(segmentation_mask, show_progress=show_progress)

        stain_intensities = create_cell_x_signal(
            g4x_obj=g4x_obj,
            mask=segmentation_mask,
            signal_list=g4x_obj.stains,
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
    mask_props = pl.LazyFrame(prop_dict, schema=schema).sort(c.CELL_ID_NAME)
    return mask_props


# region transcripts
def create_cell_x_gene(
    g4x_obj: 'G4Xoutput',
    tx_table: pl.LazyFrame,
    segmentation_mask: str | np.ndarray = DEFAULT_INPUT,
    gene_labels: str | np.ndarray = DEFAULT_INPUT,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:

    log = logger or LOGGER
    log.info('Creating cell x gene matrix')

    if segmentation_mask == DEFAULT_INPUT:
        mask = g4x_obj.load_segmentation(expanded=True)
        tx_table = intersect_tx_with_cells(tx_table, mask, column_name=c.CELL_ID_NAME)

        mask = g4x_obj.load_segmentation(expanded=False)
        tx_table = intersect_tx_with_cells(tx_table, mask, column_name='in_nucleus')
        cell_frame = _cell_frame(mask)
        del mask
    else:
        tx_table = intersect_tx_with_cells(tx_table, segmentation_mask)
        cell_frame = _cell_frame(segmentation_mask)

    if gene_labels == DEFAULT_INPUT:
        gene_labels = g4x_obj.genes

    cell_by_gene = (
        tx_table.filter(pl.col(c.CELL_ID_NAME) != 0)
        .group_by(c.CELL_ID_NAME, c.GENE_ID_NAME)
        .agg(pl.len().alias('counts'))
        .sort(c.GENE_ID_NAME)
        .pivot(on=c.GENE_ID_NAME, values='counts', index=c.CELL_ID_NAME, on_columns=gene_labels)
    )

    # Adding missing cells with zero counts
    cell_by_gene = cell_frame.join(cell_by_gene, on=c.CELL_ID_NAME, how='left')

    existing = cell_by_gene.collect_schema().names()
    if not set(existing[1:]) == set(gene_labels):
        raise ValueError('Mismatch between cell_by_gene columns and gene_labels')

    cell_by_gene = cell_by_gene.select([c.CELL_ID_NAME] + gene_labels)
    cell_by_gene = cell_by_gene.fill_null(0).sort(c.CELL_ID_NAME)

    log.info('Cell x gene matrix created successfully')

    if return_lazy:
        return cell_by_gene, tx_table
    return cell_by_gene.collect(), tx_table.collect()


def intersect_tx_with_cells(
    tx_table: pl.DataFrame | pl.LazyFrame, mask: np.ndarray, column_name: str = c.CELL_ID_NAME
) -> pl.DataFrame:
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = tx_table.select(coord_order).collect().to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]
    tx_table = tx_table.with_columns(pl.lit(cell_ids).alias(column_name))
    return tx_table


# region protein
def create_cell_x_signal(
    g4x_obj: 'G4Xoutput',
    mask: np.ndarray,
    signal_list: list[str],
    suffix: str = '_intensity_mean',
    show_progress: bool | None = None,
    return_lazy: bool = True,
    logger: logging.Logger | None = None,
    **kwargs,
) -> pl.LazyFrame:

    log = logger or LOGGER
    log.info('Creating cell x protein matrix for %d signals.', len(signal_list))

    bead_mask = g4x_obj.load_bead_mask()
    bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
    mask_flat = mask.ravel()

    channel_df = _cell_frame(mask)

    if show_progress is None:
        show_progress = sys.stderr.isatty()

    for signal_name in tqdm(signal_list, desc='Extracting image signals', disable=not show_progress):
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

    channel_df = channel_df.sort(c.CELL_ID_NAME)
    log.info('Cell x protein matrix created successfully')
    if return_lazy:
        return channel_df
    return channel_df.collect()


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


# region utilities
def _cell_frame(segmentation_mask: np.ndarray, lazy: bool = True) -> int:
    cell_labels = np.unique(segmentation_mask)
    cell_labels = cell_labels[cell_labels != 0]
    lf = pl.LazyFrame(cell_labels, schema={c.CELL_ID_NAME: pl.Int32}).sort(c.CELL_ID_NAME)
    return lf if lazy else lf.collect()
