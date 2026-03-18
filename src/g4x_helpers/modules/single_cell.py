import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from pandas.errors import PerformanceWarning
from skimage import measure
from tqdm import tqdm

from .. import constants as c
from .. import io

if TYPE_CHECKING:
    from anndata import AnnData

    from ..g4x_output import G4Xoutput

rapids = io.check_rapids()


def cell_frame(g4x_obj: 'G4Xoutput', lazy: bool = True) -> int:
    lf = pl.LazyFrame(g4x_obj.cell_labels, schema={c.CELL_ID_NAME: pl.Int32}).sort(c.CELL_ID_NAME)
    return lf if lazy else lf.collect()


# region create files
def create_cell_metadata(g4x_obj: 'G4Xoutput', mask: np.ndarray | None):
    cell_meta = init_cell_metadata(g4x_obj).collect()

    source = 'custom' if mask is not None else 'g4x-default'

    # TODO this might break if someone changes the .npz file
    if mask is None:
        mask_props = extract_cell_props_g4x(g4x_obj)
    else:
        mask_props = extract_cell_props(mask)

    mask_props = mask_props.with_columns(pl.lit(source).alias('seg_source'))

    if mask_props[c.CELL_ID_NAME].equals(cell_meta[c.CELL_ID_NAME]):
        cell_meta = cell_meta.join(mask_props, on=c.CELL_ID_NAME, how='left')
    else:
        raise ValueError('The CELL_ID columns in cell_meta and mask_props do not match.')

    return cell_meta


def create_cell_x_gene(g4x_obj: 'G4Xoutput', return_lazy: bool = False) -> tuple[pl.DataFrame, pl.DataFrame]:
    reads = pl.scan_csv(g4x_obj.tree.TranscriptTable.path)

    cell_by_gene = (
        reads.filter(pl.col(c.CELL_ID_NAME) != 0)
        .group_by(c.CELL_ID_NAME, 'gene_name')
        .agg(pl.len().alias('counts'))
        .sort('gene_name')
        .pivot(on='gene_name', values='counts', index=c.CELL_ID_NAME, on_columns=g4x_obj.genes)
    )

    # Adding missing cells with zero counts
    cell_by_gene = cell_frame(g4x_obj).join(cell_by_gene, left_on='seg_cell_id', right_on=c.CELL_ID_NAME, how='left')

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
    nuc_name = Path(c.NUC_IMG).name.split('.')[0]
    cyt_name = Path(c.CYT_IMG).name.split('.')[0]
    channel_name_map[nuc_name] = nuc_name + 'stain'
    channel_name_map[cyt_name] = cyt_name + 'stain'

    # TODO return here when bead masking is implemented
    # bead_mask = g4x_obj.load_bead_mask()
    # bead_mask_flat = bead_mask.ravel() if bead_mask is not None else None
    mask_flat = mask.ravel()

    # TODO make lazy implementation
    cf = cell_frame(g4x_obj, lazy=False)

    for signal_name in tqdm(signal_list, desc='Extracting protein signal'):
        if signal_name == nuc_name:
            signal_img = g4x_obj.load_nuclear_image(use_cache=cached)
        elif signal_name == cyt_name:
            signal_img = g4x_obj.load_cytoplasmic_image(use_cache=cached)
        else:
            signal_img = g4x_obj.load_protein_image(protein=signal_name, use_cache=cached)

        ch_label = f'{channel_name_map[signal_name]}_intensity_mean'

        intensities = image_intensity_extraction(
            signal_img,
            mask_flat=mask_flat,
            bead_mask_flat=None,
        )

        cf = cf.with_columns(pl.Series(name=ch_label, values=intensities))

    return cf


# TODO addition of protein values is missing
def create_adata(g4x_obj: 'G4Xoutput', cell_x_gene: pl.DataFrame, cell_metadata: pl.DataFrame | None = None):
    from anndata import AnnData
    from scipy.sparse import csr_matrix

    X = cell_x_gene.drop(c.CELL_ID_NAME).to_numpy().astype(np.uint16)
    X = csr_matrix(X)

    # if cell_metadata is None:
    #     cell_metadata = create_cell_metadata(g4x_obj, mask=None)

    obs_df = cell_metadata.to_pandas().set_index(c.CELL_ID_NAME)
    obs_df.index = obs_df.index.astype(str)

    gene_ids = pl.Series(name=c.GENE_ID_NAME, values=cell_x_gene.columns[1:])
    var_df = pl.DataFrame(gene_ids).with_columns(pl.lit('tx').alias('modality'))

    panel_type = (
        io.parse_input_manifest(g4x_obj.tree.TranscriptPanel.path)
        .unique('gene_name')
        .select('gene_name', 'probe_type')
        .rename({'gene_name': c.GENE_ID_NAME})
    )

    var_df = (
        var_df.join(
            panel_type,
            on=c.GENE_ID_NAME,
            how='left',
        )
        .to_pandas()
        .set_index(c.GENE_ID_NAME)
    )
    var_df.index = var_df.index.astype(str)

    adata = AnnData(X=X, obs=obs_df, var=var_df)

    ctrl_types = ['NCP', 'NCS', 'GCP']
    for ctrl_type in ctrl_types:
        adata.var[ctrl_type.lower()] = adata.var['probe_type'] == ctrl_type

    ctrl_types_lower = [c.lower() for c in ctrl_types]
    sc.pp.calculate_qc_metrics(adata, qc_vars=ctrl_types_lower, inplace=True, percent_top=None)

    adata.var.drop(columns=ctrl_types_lower, inplace=True)

    return adata.copy()


# region processing
def process_adata(adata: 'AnnData'):
    ## 1. filter on cell size, remove top and bottom 1%
    ## 2. remove genes not expressed by at least 5% of remaining cells
    ## 3. filter on total transcripts and unique genes remove top 1% and bottom 5%

    adata = filter_by_qtiles(adata, obs_key='nuclei_area_um', quantiles=(0.01, 0.99))

    min_cells = int(0.05 * adata.n_obs)
    sc.pp.filter_genes(adata, min_cells=min_cells, inplace=True)

    sc.pp.filter_cells(adata, min_counts=10)

    adata = filter_by_qtiles(adata, obs_key='total_counts', quantiles=(0.05, 0.99))
    adata = filter_by_qtiles(adata, obs_key='n_genes_by_counts', quantiles=(0.05, 0.99))

    # normalize data
    adata.layers['counts'] = adata.X
    print(adata.obs['total_counts'].min())

    sc.pp.normalize_total(adata)
    return adata

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

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PerformanceWarning)
        sc.tl.rank_genes_groups(
            adata,
            groupby=f'leiden_{res:0.2f}',
            use_raw=False,
            method='wilcoxon',
            pts=True,
            key_added=f'leiden_{res:0.2f}_rank_genes_groups',
        )
    return adata


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


def extract_cell_props(mask: np.ndarray, mask_name: str | None = None) -> pl.DataFrame:
    props = measure.regionprops(mask)

    prop_dict = []
    # Loop through each region to get the area and centroid, with a progress bar

    prefix = f' {mask_name} ' if mask_name else ' '
    for prop in tqdm(props, desc=f'Extracting{prefix}mask properties'):
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
    prop_dict_df = pl.DataFrame(prop_dict, schema=schema).sort(c.CELL_ID_NAME)
    return prop_dict_df


def extract_cell_props_g4x(g4x_obj: 'G4Xoutput') -> pl.DataFrame:
    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.tree.Segmentation.path, expected_shape=g4x_obj.shape, labels_key='nuclei'
    )
    mask_props_nuc = extract_cell_props(mask=seg_mask, mask_name='nuclei')

    seg_mask = io.import_segmentation(
        seg_path=g4x_obj.tree.Segmentation.path, expected_shape=g4x_obj.shape, labels_key='nuclei_exp'
    )
    mask_props_exp = extract_cell_props(mask=seg_mask, mask_name='nuclei_exp')
    del seg_mask

    mask_props_nuc = mask_props_nuc.rename({'area_um': 'nuclei_area_um'})
    mask_props_exp = mask_props_exp.rename({'area_um': 'wholecell_area_um'}).drop([c.CELL_COORD_X, c.CELL_COORD_Y])
    mask_props = mask_props_nuc.join(mask_props_exp, on=c.CELL_ID_NAME)

    mask_props = mask_props.select(
        [c.CELL_ID_NAME, c.CELL_COORD_X, c.CELL_COORD_Y, 'nuclei_area_um', 'wholecell_area_um']
    )

    return mask_props


# region methods
def init_cell_metadata(g4x_obj: 'G4Xoutput'):
    df = cell_frame(g4x_obj)

    df = df.with_columns(
        pl.lit(g4x_obj.sample_id).alias('sample_id'),
        pl.lit(g4x_obj.tissue_type).alias('tissue_type'),
        pl.lit(g4x_obj.block).alias('block'),
    )
    return df


# def filter_by_quantiles(adata: 'AnnData', obs_key: str, quantiles: tuple[float, float] = (0.0, 1.0)):
#     if obs_key not in adata.obs:
#         raise ValueError(f"obs_key '{obs_key}' not found in adata.obs")

#     qs = adata.obs[obs_key].quantile([quantiles[0], quantiles[1]]).to_numpy(dtype=int)

#     adata = adata[(adata.obs[obs_key] >= qs[0]) & (adata.obs[obs_key] <= qs[1])]
#     return adata.copy()


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


# def reorder_clusters_by_size(clust_umap, key: str, prefix='C') -> pl.DataFrame:
#     clust_sorted = clust_umap.group_by(key).agg(pl.len()).sort('len', descending=True)
#     size_order = clust_sorted[key].to_list()
#     numeric_order = np.arange(len(size_order)).astype(str).tolist()

#     numeric_order = [prefix + str(uid) for uid in numeric_order]
#     order_map = dict(zip(size_order, numeric_order))

#     new_umap = clust_umap.with_columns(pl.col(key).replace(order_map)).sort('seg_cell_id')
#     return new_umap


def reorder_clusters_by_size(obs: pd.DataFrame, key: str, prefix: str = 'C') -> pd.DataFrame:
    # Ensure consistent grouping behavior
    clust_sorted = (
        obs.groupby(key, observed=False)  # explicit to silence warning
        .size()
        .sort_values(ascending=False)
    )

    size_order = clust_sorted.index.tolist()

    # Create new labels
    numeric_order = [f'{prefix}{i}' for i in range(len(size_order))]

    # Mapping old -> new labels
    order_map = dict(zip(size_order, numeric_order))

    reordered_obs = obs.copy()

    # Handle categorical safely
    if isinstance(reordered_obs[key].dtype, pd.CategoricalDtype):
        # Convert to string (or object) before replacing
        reordered_obs[key] = reordered_obs[key].astype(str)

    reordered_obs[key] = reordered_obs[key].replace(order_map)

    reordered_obs[key] = pd.Categorical(
        reordered_obs[key],
        categories=numeric_order,  # this defines the order!
        ordered=True,
    )

    return reordered_obs


def _filter_axis(
    adata,
    axis: str,
    key: str,
    val_min: float | None = None,
    val_max: float | None = None,
):
    if axis == 'obs':
        df = adata.obs
        values = df[key].to_numpy()
        mask = np.ones(adata.n_obs, dtype=bool)
    elif axis == 'var':
        df = adata.var
        values = df[key].to_numpy()
        mask = np.ones(adata.n_vars, dtype=bool)
    else:
        raise ValueError("axis must be 'obs' or 'var'")

    if key not in df:
        raise ValueError(f"key '{key}' not found in adata.{axis}")

    if val_min is not None:
        mask &= values >= val_min
    if val_max is not None:
        mask &= values <= val_max

    if axis == 'obs':
        return adata[mask].copy()
    else:
        return adata[:, mask].copy()


def filter_by_limits(
    adata,
    key: str,
    axis: str = 'obs',
    val_min: float | None = None,
    val_max: float | None = None,
):
    return _filter_axis(
        adata=adata,
        axis=axis,
        key=key,
        val_min=val_min,
        val_max=val_max,
    )


def filter_by_qtiles(
    adata,
    key: str,
    axis: str = 'obs',
    q_min: float = 0,
    q_max: float = 1,
):
    if not (0 <= q_min <= 1 and 0 <= q_max <= 1):
        raise ValueError('q_min and q_max must be between 0 and 1')
    if q_min > q_max:
        raise ValueError('q_min must be <= q_max')

    df = adata.obs if axis == 'obs' else adata.var
    if axis not in {'obs', 'var'}:
        raise ValueError("axis must be 'obs' or 'var'")
    if key not in df:
        raise ValueError(f"key '{key}' not found in adata.{axis}")

    val_min, val_max = df[key].quantile([q_min, q_max]).to_numpy()

    return _filter_axis(
        adata=adata,
        axis=axis,
        key=key,
        val_min=val_min,
        val_max=val_max,
    )
