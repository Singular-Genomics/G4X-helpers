from typing import TYPE_CHECKING

import numpy as np
import polars as pl
import scanpy as sc
from skimage import measure
from tqdm import tqdm

from .. import constants, io

if TYPE_CHECKING:
    from anndata import AnnData

    from ..g4x_output import G4Xoutput


CELL_ID_NAME = 'seg_cell_id'
GENE_ID_NAME = 'gene_id'
CELL_COORD_X = 'cell_x'
CELL_COORD_Y = 'cell_y'


def cell_frame(g4x_obj: 'G4Xoutput', lazy: bool = True) -> int:
    lf = pl.LazyFrame(g4x_obj.cell_labels, schema={CELL_ID_NAME: pl.Int32}).sort(CELL_ID_NAME)
    return lf if lazy else lf.collect()


# TODO addition of protein values is missing
def create_adata(g4x_obj: 'G4Xoutput', cell_metadata: pl.DataFrame, cell_x_gene: pl.DataFrame):
    from anndata import AnnData
    from scipy.sparse import csr_matrix

    X = cell_x_gene.drop(CELL_ID_NAME).to_numpy().astype(np.uint16)
    X = csr_matrix(X)

    obs_df = cell_metadata.to_pandas().set_index(CELL_ID_NAME)
    obs_df.index = obs_df.index.astype(str)

    gene_ids = pl.Series(name=GENE_ID_NAME, values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    panel_type = (
        io.parse_input_manifest(g4x_obj.data_dir / 'transcript_panel.csv')
        .unique('gene_name')
        .select('gene_name', 'probe_type')
        .rename({'gene_name': GENE_ID_NAME})
    )

    var_df = (
        var_df.join(
            panel_type,
            on=GENE_ID_NAME,
            how='left',
        )
        .to_pandas()
        .set_index(GENE_ID_NAME)
    )

    # # TODO I think the correct thing to do here is to set the index to the gene_id
    # var_df.index = var_df.index.astype(str)
    # .set_index('label')

    adata = AnnData(X=X, obs=obs_df, var=var_df)
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
    return adata.copy()


def process_adata(adata: 'AnnData'):
    ## 1. filter on cell size, remove top and bottom 1%
    ## 2. remove genes not expressed by at least 5% of remaining cells
    ## 3. filter on total transcripts and unique genes remove top 1% and bottom 5%

    adata = filter_by_quantiles(adata, obs_key='nuclei_area_um', quantiles=(0.01, 0.99))

    min_cells = int(0.05 * adata.n_obs)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)

    adata = filter_by_quantiles(adata, obs_key='total_counts', quantiles=(0.05, 0.99))
    adata = filter_by_quantiles(adata, obs_key='n_genes_by_counts', quantiles=(0.05, 0.99))

    # normalize data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)

    res = 1
    sc.tl.leiden(
        adata,
        resolution=res,
        objective_function='modularity',
        n_iterations=-1,
        random_state=42,
        flavor='igraph',
        key_added=f'leiden_{res:0.2f}',
    )

    sc.tl.umap(adata)

    sc.tl.rank_genes_groups(
        adata,
        groupby=f'leiden_{res:0.2f}',
        use_raw=False,
        method='wilcoxon',
        pts=True,
        key_added=f'leiden_{res:0.2f}_rank_genes_groups',
    )
    return adata


def create_cell_metadata(g4x_obj: 'G4Xoutput', mask: np.ndarray | None):
    cell_meta = init_cell_metadata(g4x_obj).collect()

    source = 'custom' if mask is not None else 'g4x'

    # TODO this might break if someone changes the .npz file
    if mask is None:
        mask_props = extract_cell_props_g4x(g4x_obj)
    else:
        mask_props = extract_cell_props(mask)

    mask_props = mask_props.with_columns(pl.lit(source).alias('source'))

    if mask_props[CELL_ID_NAME].equals(cell_meta[CELL_ID_NAME]):
        cell_meta = cell_meta.join(mask_props, on=CELL_ID_NAME, how='left')
    else:
        raise ValueError('The CELL_ID columns in cell_meta and mask_props do not match.')

    return cell_meta


def init_cell_metadata(g4x_obj: 'G4Xoutput'):
    df = cell_frame(g4x_obj)

    df = df.with_columns(
        pl.lit(g4x_obj.sample_id).alias('sample_id'),
        pl.lit(g4x_obj.tissue_type).alias('tissue_type'),
        pl.lit(g4x_obj.block).alias('block'),
    )
    return df


def extract_cell_props(mask: np.ndarray, mask_name: str | None = None) -> pl.DataFrame:
    props = measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar

    prefix = f' {mask_name} ' if mask_name else ' '
    for prop in tqdm(props, desc=f'Extracting{prefix}mask properties'):
        label = prop.label  # The label (mask id)

        cell_y, cell_x = prop.centroid  # coordinate order: 'yx' (row, col)

        px_to_um_area = constants.PIXEL_SIZE_MICRONS**2
        area_um = prop.area * px_to_um_area  # prop.area is in pixels

        prop_dict.append(
            {
                CELL_ID_NAME: label,
                'area_um': area_um,
                CELL_COORD_X: cell_x,
                CELL_COORD_Y: cell_y,
            }
        )
    schema = {
        CELL_ID_NAME: pl.Int32,
        'area_um': pl.Float32,
        CELL_COORD_X: pl.Float32,
        CELL_COORD_Y: pl.Float32,
    }
    prop_dict_df = pl.DataFrame(prop_dict, schema=schema).sort(CELL_ID_NAME)
    return prop_dict_df


def extract_cell_props_g4x(g4x_obj: 'G4Xoutput') -> pl.DataFrame:
    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.segmentation_path, expected_shape=g4x_obj.shape, labels_key='nuclei'
    )
    mask_props_nuc = extract_cell_props(mask=seg_mask, mask_name='nuclei')

    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.segmentation_path, expected_shape=g4x_obj.shape, labels_key='nuclei_exp'
    )
    mask_props_exp = extract_cell_props(mask=seg_mask, mask_name='nuclei_exp')
    del seg_mask

    mask_props_nuc = mask_props_nuc.rename({'area_um': 'nuclei_area_um'})
    mask_props_exp = mask_props_exp.rename({'area_um': 'wholecell_area_um'}).drop([CELL_COORD_X, CELL_COORD_Y])
    mask_props = mask_props_nuc.join(mask_props_exp, on=CELL_ID_NAME)

    mask_props = mask_props.select([CELL_ID_NAME, CELL_COORD_X, CELL_COORD_Y, 'nuclei_area_um', 'wholecell_area_um'])

    return mask_props


def create_cell_x_gene(g4x_obj: 'G4Xoutput', return_lazy: bool = False) -> tuple[pl.DataFrame, pl.DataFrame]:
    reads = pl.scan_csv(g4x_obj.transcript_table_path)

    cell_by_gene = (
        reads.filter(pl.col(CELL_ID_NAME) != 0)
        .group_by(CELL_ID_NAME, 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index=CELL_ID_NAME, on_columns=g4x_obj.genes)
    )

    # Adding missing cells with zero counts
    cell_by_gene = cell_frame(g4x_obj).join(cell_by_gene, left_on='seg_cell_id', right_on=CELL_ID_NAME, how='left')

    existing = cell_by_gene.collect_schema().names()
    if not set(existing[1:]) == set(g4x_obj.genes):
        raise ValueError('Mismatch between cell_by_gene columns and g4x_obj.genes')

    cell_by_gene = cell_by_gene.select(['seg_cell_id'] + g4x_obj.genes)
    cell_by_gene = cell_by_gene.fill_null(0)

    if return_lazy:
        return cell_by_gene
    return cell_by_gene.collect()


def create_cell_x_protein(
    g4x_obj: 'G4Xoutput',
    mask: np.ndarray,
    signal_list: list[str] | None = None,
    cached: bool = False,
) -> pl.DataFrame | pl.LazyFrame:
    # if signal_list is None:
    #     signal_list = ['nuclear', 'eosin'] + g4x_obj.proteins

    print(f'Creating cell x protein matrix for {len(signal_list)} signals.')

    # TODO eosin channel has been renamed
    channel_name_map = {protein: protein for protein in signal_list}
    channel_name_map['nuclear'] = 'nuclearstain'
    channel_name_map['eosin'] = 'cytoplasmicstain'

    # TODO return here when bead masking is implemented
    # bead_mask = g4x_obj.load_bead_mask()
    # bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
    mask_flat = mask.ravel()

    # TODO make lazy implementation
    cf = cell_frame(g4x_obj, lazy=False)

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

        cf = cf.with_columns(pl.Series(name=ch_label, values=intensities))

    return cf


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


def filter_by_quantiles(adata: 'AnnData', obs_key: str, quantiles: tuple[float, float] = (0.0, 1.0)):
    if obs_key not in adata.obs:
        raise ValueError(f"obs_key '{obs_key}' not found in adata.obs")

    qs = adata.obs[obs_key].quantile([quantiles[0], quantiles[1]]).to_numpy(dtype=int)

    adata = adata[(adata.obs[obs_key] >= qs[0]) & (adata.obs[obs_key] <= qs[1])]
    return adata.copy()


def intersect_tx_with_cells(tx_table, mask):
    coord_order = ['y_pixel_coordinate', 'x_pixel_coordinate']
    tx_coords = tx_table.select(coord_order).collect().to_numpy().astype(int)
    cell_ids = mask[tx_coords[:, 0], tx_coords[:, 1]]
    return cell_ids
