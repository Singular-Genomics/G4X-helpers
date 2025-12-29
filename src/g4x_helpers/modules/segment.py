from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas
import numpy as np
import polars as pl
import scanpy as sc
import shapely.affinity
import skimage.measure
from anndata import AnnData
from geopandas import GeoDataFrame
from rasterio.features import rasterize
from rasterio.transform import Affine
from shapely.geometry import Polygon
from skimage.measure._regionprops import RegionProperties
from tqdm import tqdm

from .. import utils
from .workflow import OutSchema, workflow

if TYPE_CHECKING:
    from geopandas.geodataframe import GeoDataFrame

    from g4x_helpers.models import G4Xoutput


@workflow
def apply_segmentation(
    g4x_obj: 'G4Xoutput',
    labels: np.ndarray,
    
    out_dir: Path,
    *,
    tx_table: Path | None = None,
    skip_protein_extraction: bool = False,
    signal_list: list[str] | None = None,
    logger: logging.Logger,
    **kwargs,
):
    out_tree = OutSchema(out_dir, subdirs=['single_cell_data', 'rna', 'segmentation'])

    cell_x_gene_out = out_tree.single_cell_data / 'cell_by_transcript.csv'
    cell_x_protein_out = out_tree.single_cell_data / 'cell_by_protein.csv'
    cell_metadata_out = out_tree.single_cell_data / 'cell_metadata.csv'
    feature_matrix_out = out_tree.single_cell_data / 'feature_matrix.h5'
    tx_table_out = out_tree.rna / 'transcript_table.csv'
    seg_mask_out = out_tree.segmentation / 'segmentation_mask.npz'

    create_source = kwargs.get('create_source', False)
    nuc_mask = g4x_obj.load_segmentation(expanded=False) if create_source else None

    logger.info('Creating cell_by_transcript table.')
    if tx_table is None:
        tx_table = g4x_obj.transcript_table_path
    cell_x_gene, tx_table = create_cell_x_gene(g4x_obj=g4x_obj, tx_table=tx_table, mask=labels, nucleus_mask=nuc_mask)

    logger.info(f'Writing cell_by_transcript matrix to: {cell_x_gene_out.relative_to(out_dir)}.')
    utils.write_csv_gz(df=cell_x_gene, path=cell_x_gene_out)

    logger.info(f'Writing transcript table after adding cell-ids to: {tx_table_out.relative_to(out_dir)}.')
    utils.write_csv_gz(df=tx_table, path=tx_table_out)

    if signal_list is None:
        if g4x_obj.includes_protein:
            signal_list = ['nuclear', 'eosin'] + g4x_obj.proteins
        else:
            signal_list = ['nuclear', 'eosin']
        logger.info(f'Processing image channels: {signal_list}.')
    else:
        logger.info('Processing all image channels.')

    if skip_protein_extraction and g4x_obj.includes_protein:
        logger.warning('Skipping protein extraction. Will try existing cell_by_protein table.')
        
        cxp_path = g4x_obj.data_dir / 'single_cell_data' / 'cell_by_protein.csv.gz'
        cxp_short = Path(cxp_path.parent.name) / cxp_path.name
        print(f'Cell_x_protein table not provided, loading from: {cxp_short}')
        cell_x_protein = pl.read_csv(cxp_path)

        if signal_list is not None:
            logger.info(f'Selecting channels: {signal_list}')
            signal_list = ['nuclear', 'eosin'] + g4x_obj.proteins
            cell_x_protein = cell_x_protein.select(['label'] + [f'{signal_name}_intensity_mean' for signal_name in signal_list])
    else:
        logger.info('Gathering image intensities.')
    
        cell_x_protein = create_cell_x_protein(
            g4x_obj=g4x_obj,
            mask=labels,
            signal_list=signal_list,
            cached=False,
        )

        if g4x_obj.includes_protein:
            logger.info(f'Writing cell_by_protein matrix to: {cell_x_protein_out.relative_to(out_dir)}.')
            utils.write_csv_gz(df=cell_x_protein, path=cell_x_protein_out)

    logger.info('Creating adata and metadata.')
    cell_metadata_pre = create_cell_metadata(
        g4x_obj, cell_x_protein=cell_x_protein, mask=labels, create_source=create_source
    )

    adata = create_adata(
        g4x_obj=g4x_obj,
        cell_x_gene=cell_x_gene,
        cell_metadata=cell_metadata_pre,
    )

    logger.info(f'Writing adata to: {feature_matrix_out.relative_to(out_dir)}.')
    adata.write_h5ad(feature_matrix_out)

    logger.info(f'Writing cell metadata table to: {cell_metadata_out.relative_to(out_dir)}.')
    cell_metadata = adata.obs.rename(columns={'cell_id': 'label'}).set_index('label')
    utils.write_csv_gz(df=cell_metadata, path=cell_metadata_out)

    logger.info(f'Writing segmentation mask to: {seg_mask_out.relative_to(out_dir)}.')
    np.savez(seg_mask_out, cell_labels=labels)


def try_load_segmentation(cell_labels: str, expected_shape: tuple[int], labels_key: str | None = None) -> np.ndarray:
    SUPPORTED_MASK_FILETYPES = {'.npy', '.npz', '.geojson'}
    ## load new segmentation
    cell_labels = utils.validate_path(cell_labels, must_exist=True, is_dir_ok=False, is_file_ok=True)
    suffix = cell_labels.suffix.lower()

    if suffix not in SUPPORTED_MASK_FILETYPES:
        raise ValueError(f'{suffix} is not a supported file type.')

    if suffix == '.npz':
        with np.load(cell_labels) as labels:
            available_keys = list(labels.keys())

            if labels_key:  # if a key is specified
                if labels_key not in labels:
                    raise KeyError(f"Key '{labels_key}' not found in .npz; available keys: {available_keys}")
                seg = labels[labels_key]

            else:
                if len(available_keys) != 1:
                    raise ValueError(
                        f'Expected exactly one key in .npz but found {len(available_keys)}: {available_keys}. '
                        "Please specify a key using 'labels_key'."
                    )
                seg = labels[available_keys[0]]

    elif suffix == '.npy':
        # .npy: directly returns the array, no context manager available
        if labels_key is not None:
            print('file is .npy, ignoring provided labels_key.')
        seg = np.load(cell_labels, allow_pickle=False)

    elif suffix == '.geojson':
        gdf = geopandas.read_file(cell_labels)

        if labels_key is not None:
            if labels_key not in gdf.columns:
                raise KeyError(f"Column '{labels_key}' not found in GeoJSON; available columns: {gdf.columns.tolist()}")

            # ensure that a column named 'label' exists
            gdf['label'] = gdf[labels_key]

        else:
            if 'label' not in gdf.columns:
                raise ValueError(
                    "No column named 'label' found in GeoJSON. Please specify which column to use for labels via labels_key."
                )

        print('Rasterizing provided GeoDataFrame.')
        seg = rasterize_polygons(gdf=gdf, target_shape=expected_shape)

    # validate shape for final numpy arrays
    if seg.shape != expected_shape:
        raise ValueError(f'provided mask shape {seg.shape} does not match G4X sample shape {expected_shape}')

    return seg


def get_cell_ids(sample_id: str, mask) -> pl.DataFrame:
    print('Getting cell-IDs from segmentation mask.')
    seg_ids = np.unique(mask[mask != 0])
    cell_ids = (
        pl.Series(name='segmentation_label', values=seg_ids[seg_ids != 0])
        .to_frame()
        .with_columns((f'{sample_id}-' + pl.col('segmentation_label').cast(pl.String)).alias('label'))
        .select(['label', 'segmentation_label'])
    )
    return cell_ids


def get_mask_properties(cell_ids: pl.DataFrame, mask: np.ndarray) -> pl.DataFrame:
    print('Getting regionprops.')
    props = skimage.measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar
    for prop in tqdm(props, desc='Extracting mask properties'):
        label = prop.label  # The label (mask id)
        area = prop.area  # Area: count of pixels
        centroid = prop.centroid  # Centroid: (row, col)

        # assuming coordinate order: 'yx':
        cell_y, cell_x = centroid

        px_to_um_area = 0.3125**2

        prop_dict.append(
            {
                'segmentation_label': label,
                'area_um': area * px_to_um_area,
                'cell_x': cell_x,
                'cell_y': cell_y,
            }
        )
    schema = {
        'segmentation_label': pl.Int32,
        'area_um': pl.Float32,
        'cell_x': pl.Float32,
        'cell_y': pl.Float32,
    }
    prop_dict_df = pl.DataFrame(prop_dict, schema=schema)

    prop_dict_df = cell_ids.join(prop_dict_df, on='segmentation_label', how='left')
    return prop_dict_df


def create_cell_metadata(
    g4x_obj: G4Xoutput,
    cell_x_protein: pl.DataFrame | None = None,
    mask: np.ndarray | None = None,
    create_source: bool = False,
) -> pl.DataFrame:
    if create_source:
        cell_props = build_g4x_cell_properties(g4x_obj)
    else:
        cell_props = build_custom_cell_properties(g4x_obj, mask)

    if cell_x_protein is None:
        cxp_path = g4x_obj.data_dir / 'single_cell_data' / 'cell_by_protein.csv.gz'
        cxp_short = Path(cxp_path.parent.name) / cxp_path.name
        print(f'Cell_x_protein table not provided, loading from: {cxp_short}')
        cell_x_protein = pl.read_csv(cxp_path)
        mergeable = cell_x_protein['label'].sort().equals(cell_props['label'].sort())
        if not mergeable:
            raise ValueError(
                f'Cell labels in {cxp_path} do not match those in cell properties. Please provide cell_x_protein DataFrame directly.'
            )

    cell_metadata = cell_props.join(cell_x_protein, on='label')
    return cell_metadata


def build_g4x_cell_properties(g4x_obj: G4Xoutput) -> pl.DataFrame:
    mask = g4x_obj.load_segmentation(expanded=True)
    cell_ids = get_cell_ids(g4x_obj.sample_id, mask)
    cell_props_expanded = get_mask_properties(cell_ids, mask)

    mask = g4x_obj.load_segmentation(expanded=False)
    cell_props_nuc = get_mask_properties(cell_ids, mask)
    del mask

    df = cell_props_nuc.rename({'area_um': 'nuclei_area_um'})
    df_exp = cell_props_expanded.rename({'area_um': 'nuclei_expanded_area_um'}).drop(
        ['segmentation_label', 'cell_x', 'cell_y']
    )

    segmentation_props = df.join(df_exp, on='label').select(
        ['label', 'segmentation_label', 'cell_x', 'cell_y', 'nuclei_area_um', 'nuclei_expanded_area_um']
    )

    return segmentation_props


def build_custom_cell_properties(g4x_obj: G4Xoutput, mask: np.ndarray) -> pl.DataFrame:
    cell_ids = get_cell_ids(g4x_obj.sample_id, mask)
    cell_props = get_mask_properties(cell_ids, mask)
    del mask

    segmentation_props = cell_props.select(['label', 'segmentation_label', 'cell_x', 'cell_y', 'area_um'])

    return segmentation_props


def create_cell_x_gene(
    g4x_obj: G4Xoutput,
    tx_table: Path,
    mask: np.ndarray,
    nucleus_mask: np.ndarray | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    reads = pl.scan_csv(tx_table)

    print('Assigning transcripts to segmentation labels.')
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = reads.select(coord_order).collect().to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]

    if nucleus_mask is not None:
        print('Assigning transcripts to nucleus labels.')
        nuc_ids = nucleus_mask[tx_coords[:, 0], tx_coords[:, 1]]
        reads = reads.with_columns(pl.lit(nuc_ids).alias('in_nucleus'))

    if nucleus_mask is None:
        if 'in_nucleus' in reads.collect_schema().names():
            reads = reads.drop('in_nucleus')

    reads = reads.with_columns(pl.lit(cell_ids).alias('segmentation_label')).with_columns(
        cell_id=(g4x_obj.sample_id + '-' + pl.col('segmentation_label').cast(pl.Utf8))
    )

    print('Creating cell x gene matrix.')
    cell_by_gene = (
        reads.filter(pl.col('segmentation_label') != 0)
        .group_by('segmentation_label', 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .collect()
        .pivot(on='gene_name', values='counts', index='segmentation_label')
    )

    print('Adding missing cells with zero counts.')
    cell_ids = get_cell_ids(g4x_obj.sample_id, mask)
    cell_by_gene = (
        cell_ids.join(cell_by_gene, on='segmentation_label', how='left')  # .select('segmentation_label')
        .fill_null(0)
        .drop('segmentation_label')
    )

    print('Adding missing genes with zero counts.')
    zero_count_probes = (
        utils.parse_input_manifest(file_path=g4x_obj.data_dir / 'transcript_panel.csv')
        .unique('gene_name')
        .filter(~pl.col('gene_name').is_in(cell_by_gene.columns))['gene_name']
        .to_list()
    )

    for probe in zero_count_probes:
        cell_by_gene = cell_by_gene.with_columns(pl.lit(None, dtype=pl.UInt32).alias(probe))

    ordered = [cell_by_gene.columns[0]] + sorted(cell_by_gene.columns[1:])
    cell_by_gene = cell_by_gene.select(ordered).fill_null(0)

    return cell_by_gene, reads.drop('segmentation_label')


def create_cell_x_protein(
    g4x_obj: G4Xoutput,
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

    cell_x_protein = get_cell_ids(g4x_obj.sample_id, mask).drop('segmentation_label')

    for signal_name in tqdm(signal_list, desc='Extracting protein signal'):
        if signal_name not in ['nuclear', 'eosin']:
            image_type = 'protein'
            protein = signal_name
        else:
            image_type = signal_name
            protein = None

        signal_img = g4x_obj.load_image_by_type(image_type, thumbnail=False, protein=protein, cached=cached)

        ch_label = f'{channel_name_map[signal_name]}_intensity_mean'

        intensities = image_mask_intensity_extraction_v2(
            signal_img,
            mask_flat=mask_flat,
            bead_mask_flat=None,
        )

        cell_x_protein = cell_x_protein.with_columns(pl.Series(name=ch_label, values=intensities))

    return cell_x_protein


def create_adata(
    g4x_obj: G4Xoutput,
    cell_x_gene: pl.DataFrame,
    cell_metadata: pl.DataFrame,
) -> AnnData:
    X = cell_x_gene.drop('label').to_numpy().astype(np.uint16)

    obs_df = cell_metadata.to_pandas()  # .set_index('label')
    obs_df.index = obs_df.index.astype(str)

    gene_ids = pl.Series(name='gene_id', values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    panel_type = (
        utils.parse_input_manifest(g4x_obj.data_dir / 'transcript_panel.csv')
        .unique('gene_name')
        .select('gene_name', 'probe_type')
        .rename({'gene_name': 'gene_id'})
    )

    var_df = var_df.join(
        panel_type,
        on='gene_id',
        how='left',
    ).to_pandas()
    var_df.index = var_df.index.astype(str)

    adata = AnnData(X=X, obs=obs_df, var=var_df)
    adata.obs = adata.obs.rename(columns={'label': 'cell_id'})  # if you still want this

    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)

    return adata


# def image_mask_intensity_extraction(
#     img: np.ndarray,
#     mask_flat: np.ndarray,
#     bead_mask_flat: np.ndarray | None = None,
# ) -> pl.DataFrame | pl.LazyFrame:
#     img_flat = img.ravel()

#     ## optional masking with beads
#     if bead_mask_flat is not None:
#         bead_mask_flat = ~bead_mask_flat
#         mask_flat = mask_flat[bead_mask_flat]
#         img_flat = img_flat[bead_mask_flat]

#     # intensity_means = np.bincount(mask_flat, weights=img_flat)[1:] / np.bincount(mask_flat)[1:]
#     counts = np.bincount(mask_flat)[1:]
#     sums = np.bincount(mask_flat, weights=img_flat)[1:]

#     # Avoid division by zero:
#     intensity_means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts != 0)
#     return intensity_means


def image_mask_intensity_extraction_v2(
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

    unique_labels, inv = np.unique(labels, return_inverse=True)
    # `inv` now contains remapped labels 0..n_unique-1

    # Compute sums and counts using remapped labels
    sums = np.bincount(inv, weights=pixels)
    counts = np.bincount(inv)

    # Safe divide
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts != 0)

    return means


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


def rasterize_polygons(gdf: GeoDataFrame, target_shape: tuple) -> np.ndarray:
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
