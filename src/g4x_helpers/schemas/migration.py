import shutil
from pathlib import Path

import anndata as ad
import polars as pl

from .. import utils

files_to_migrate = {
    'bin': ('g4x_viewer/{sample_id}.bin', 'g4x_viewer/{sample_id}_segmentation.bin', False),
    'tar': ('g4x_viewer/{sample_id}.tar', 'g4x_viewer/{sample_id}_transcripts.tar', False),
    'ome.tiff': ('g4x_viewer/{sample_id}.ome.tiff', 'g4x_viewer/{sample_id}_multiplex.ome.tiff', False),
    'raw_tx_v2': ('diagnostics/transcript_table.parquet', 'rna/raw_features.parquet', True),
    'raw_tx_v3': ('rna/transcript_table.parquet', 'rna/raw_features.parquet', True),
    'tx_table': ('rna/transcript_table.csv.gz', 'rna/transcript_table.csv.gz', True),
    'tx_panel': ('transcript_panel.csv', 'transcript_panel.csv', True),
    'feat_mtx': ('single_cell_data/feature_matrix.h5', 'single_cell_data/feature_matrix.h5', True),
    'cell_meta': ('single_cell_data/cell_metadata.csv.gz', 'single_cell_data/cell_metadata.csv.gz', True),
}

parquet_rename = {
    'TXUID': 'TXUID',
    'sequence_to_demux': 'sequence',
    'meanQS': 'confidence_score',
    'x_coord_shift': 'y_pixel_coordinate',
    'y_coord_shift': 'x_pixel_coordinate',
    'z': 'z_level',
    'transcript': 'probe_name',
    'transcript_condensed': 'gene_name',
}

parquet_col_order = ['TXUID', 'sequence', 'confidence_score', 'x_pixel_coordinate', 'y_pixel_coordinate', 'z_level']


def migrate_file(sample_base, sample_id, old_path, new_path, backup):
    backup_path = sample_base / 'migration_backup'
    backup_path.mkdir(exist_ok=True, parents=True)

    old_path = old_path.format(sample_id=sample_id)
    new_path = new_path.format(sample_id=sample_id)

    old_path_full = sample_base / old_path.format(sample_id=sample_id)
    new_path_full = sample_base / new_path.format(sample_id=sample_id)

    moved = False
    if old_path_full.exists():
        print('Moving:', old_path, '->', new_path)
        if backup:
            shutil.copy2(old_path_full, backup_path / old_path_full.name)

        shutil.move(old_path_full, new_path_full)
        moved = True
    return moved


def migrate_sample_files(sample_base, sample_id):
    moved_files = {}
    for f, paths in files_to_migrate.items():
        old_path, new_path, backup = paths
        moved_files[f] = migrate_file(sample_base, sample_id, old_path, new_path, backup)
    return moved_files


def put_back(sample_base, sample_id):
    for f in files_to_migrate:
        old_path, new_path, has_backup = files_to_migrate[f]
        if has_backup:
            file_name = Path(old_path).name.format(sample_id=sample_id)
            new_path = sample_base / 'migration_backup' / file_name
        migrate_file(sample_base, sample_id, str(new_path), old_path, False)
        if f == 'raw_tx_v2':
            (sample_base / files_to_migrate[f][1]).unlink()


def infer_parquet_schema(sample_base):
    raw_features_path = sample_base / files_to_migrate['raw_tx_v3'][1]
    parquet_lf = pl.scan_parquet(raw_features_path)
    parquet_cols = parquet_lf.collect_schema().names()

    if {'cell_id', 'gene_name', 'probe_name', 'demuxed'} <= set(parquet_cols):
        parquet_shema = 'v3'
    elif {'x_coord_shift', 'y_coord_shift', 'z', 'transcript', 'transcript_condensed'} <= set(parquet_cols):
        parquet_shema = 'v2'
    else:
        parquet_shema = 'unknown'

    if parquet_shema == 'v2':
        parquet_cols = parquet_lf.collect_schema().names()

        for c in parquet_cols:
            if c in parquet_rename:
                print(f'Renaming {c} to {parquet_rename[c]}')
                parquet_lf = parquet_lf.rename({c: parquet_rename[c]})

    return parquet_lf, parquet_shema


def migrate_tx_panel(sample_base, parquet_lf):
    tx_panel_path = sample_base / files_to_migrate['tx_panel'][1]
    tx_panel_old = pl.read_csv(tx_panel_path)
    tx_panel_cols = tx_panel_old.columns

    if {'target_condensed'} <= set(tx_panel_cols):
        # tx_panel_schema = 'v2'
        print('Building new transcript panel csv from migrated transcript table')
        tx_panel = tx_panel_from_parquet(parquet_lf, tx_panel_old)
        tx_panel.write_csv(sample_base / 'transcript_panel.csv')
    elif {'probe_name'} <= set(tx_panel_cols):
        # tx_panel_schema = 'v3'
        print('Transcript panel already in v3 schema. No migration needed.')
    else:
        # tx_panel_schema = 'unknown'
        print('Could not infer transcript panel schema. No migration performed.')


def tx_panel_from_parquet(parquet_lf, tx_panel_old):
    tx_panel = (
        parquet_lf.filter(pl.col('probe_name') != 'UNDETERMINED')
        .unique('probe_name')
        .sort('probe_name')
        .select('probe_name', 'gene_name')
        .collect()
    )

    tx_panel = (
        tx_panel.unique('probe_name')
        .join(tx_panel_old, left_on='gene_name', right_on='target_condensed', how='right')
        .rename({'target_condensed': 'gene_name'})
        .sort('panel_type', 'gene_name')
        .fill_null('NOT_CONVERTED')
    )

    unconverted = tx_panel.filter(pl.col('probe_name') == 'NOT_CONVERTED')['gene_name']
    if len(unconverted) > 0:
        print(f'Warning: {len(unconverted)} gene_names failed to map to probe_names:')
        print(unconverted.to_list())
        print('They were listed in original transcript_panel.csv but not found in transcript data.')

    return tx_panel


def write_parquet_in_order(parquet_lf, sample_base):
    raw_features_path = sample_base / files_to_migrate['raw_tx_v3'][1]
    tmp_file = raw_features_path.with_name('tmp.parquet')
    parquet_lf.select(parquet_col_order).sink_parquet(tmp_file)
    shutil.move(tmp_file, raw_features_path)


def migrate_adata(sample_base):
    adata_path = sample_base / files_to_migrate['feat_mtx'][1]
    adata = ad.read_h5ad(adata_path)

    if {'expanded_cell_x', 'expanded_cell_y'} <= set(adata.obs.columns):
        adata_schema = 'v2'

        adata.obs['cell_x'], adata.obs['cell_y'] = (
            adata.obs['cell_y'],
            adata.obs['cell_x'],
        )

        adata.obs = adata.obs.drop(columns=['expanded_cell_x', 'expanded_cell_y'])
        px_to_um_area = 0.3125**2

        adata.obs['nuclei_area'] = adata.obs['nuclei_area'] * px_to_um_area
        adata.obs['nuclei_expanded_area'] = adata.obs['nuclei_expanded_area'] * px_to_um_area

        adata.obs = adata.obs.rename(
            columns={'nuclei_area': 'nuclei_area_um', 'nuclei_expanded_area': 'nuclei_expanded_area_um'}
        )
    elif {'centroid-0', 'centroid-1'} <= set(adata.obs.columns):
        adata_schema = 'v3'

        adata.obs = adata.obs.drop(columns=['centroid-0', 'centroid-1'])
    else:
        adata_schema = 'unknown'
        print('could not detect an older adata version')

    if adata_schema != 'unknown':
        print(f'detected {adata_schema} adata schema')
        print('writing new adata and cell metadata')
        adata.write_h5ad(adata_path)

        cols = adata.obs.columns
        first_cols = ['nuclearstain_intensity_mean', 'cytoplasmicstain_intensity_mean']
        prot_cols = first_cols + [c for c in cols if '_intensity_mean' in c and c not in first_cols]
        cell_meta_cols = ['cell_id', 'cell_x', 'cell_y', 'nuclei_area_um', 'nuclei_expanded_area_um']
        last_cols = [c for c in cols if c not in cell_meta_cols and c not in prot_cols]

        new_cols = cell_meta_cols + prot_cols + last_cols

        cell_metadata = adata.obs[new_cols].rename(columns={'cell_id': 'label'})

        utils.write_csv_gz(cell_metadata, (sample_base / files_to_migrate['cell_meta'][1]).with_suffix(''))
