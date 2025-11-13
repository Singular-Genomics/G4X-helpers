import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from geopandas import GeoDataFrame

from .. import utils
from ..modules import segmentation as seg
from .decorator import workflow

if TYPE_CHECKING:
    from ..models import G4Xoutput


@workflow
def intersect_segmentation(
    g4x_obj: 'G4Xoutput',
    labels: np.ndarray | GeoDataFrame,
    out_dir: Path,
    *,
    exclude_channels: list[str] | None = None,
    logger: logging.Logger,
) -> None:
    logger.info(f'Using provided output directory: {out_dir}')

    signal_list = ['nuclear', 'eosin'] + g4x_obj.proteins

    if exclude_channels is not None:
        logger.info(f'Not processing channels: {", ".join(exclude_channels)}')
        if isinstance(exclude_channels, str):
            exclude_channels = [exclude_channels]

        signal_list = [item for item in signal_list if item not in exclude_channels]
    else:
        logger.info('Processing all channels.')

    if isinstance(labels, GeoDataFrame):
        logger.info('Rasterizing provided GeoDataFrame.')
        mask = seg.rasterize_polygons(gdf=labels, target_shape=g4x_obj.shape)
    else:
        mask = labels

    logger.info('Extracting mask properties.')
    segmentation_props = seg.get_mask_properties(g4x_obj, mask)

    logger.info('Assigning transcripts to mask labels.')
    reads_new_labels = seg.assign_tx_to_mask_labels(g4x_obj, mask)

    logger.info('Extracting protein signal.')
    cell_by_protein = seg.extract_image_signals(g4x_obj, mask, signal_list=signal_list)

    logger.info('Building output data structures.')
    cell_metadata = seg._make_cell_metadata(segmentation_props, cell_by_protein)
    cell_by_gene = seg._make_cell_by_gene(segmentation_props, reads_new_labels)
    adata = seg._make_adata(cell_by_gene, cell_metadata)

    probe_types = g4x_obj.load_adata(remove_nontargeting=False, load_clustering=False).var['probe_type']
    adata.var = adata.var.merge(probe_types, left_index=True, right_index=True, how='left')

    logger.info(f'Saving output files to {out_dir}')

    outfile = utils.create_custom_out(out_dir, 'segmentation', 'segmentation_mask.npz')
    logger.debug(f'segmentation mask --> {outfile}')
    np.savez(outfile, cell_labels=mask)

    outfile = utils.create_custom_out(out_dir, 'rna', 'transcript_table.csv')
    logger.debug(f'transcript table --> {outfile}.gz')
    reads_new_labels.write_csv(outfile)
    _ = utils.gzip_file(outfile, remove_original=True)

    outfile = utils.create_custom_out(out_dir, 'single_cell_data', 'cell_by_transcript.csv')
    logger.debug(f'cell x transcript --> {outfile}.gz')
    cell_by_gene.write_csv(outfile)
    _ = utils.gzip_file(outfile, remove_original=True)

    outfile = utils.create_custom_out(out_dir, 'single_cell_data', 'cell_by_protein.csv')
    logger.debug(f'cell x protein --> {outfile}.gz')
    cell_by_protein.write_csv(outfile)
    _ = utils.gzip_file(outfile, remove_original=True)

    outfile = utils.create_custom_out(out_dir, 'single_cell_data', 'feature_matrix.h5')
    logger.debug(f'single-cell h5 --> {outfile}')
    adata.write_h5ad(outfile)

    outfile = utils.create_custom_out(out_dir, 'single_cell_data', 'cell_metadata.csv.gz')
    logger.debug(f'cell metadata --> {outfile}')
    adata.obs.to_csv(outfile, compression='gzip')

    return adata, mask
