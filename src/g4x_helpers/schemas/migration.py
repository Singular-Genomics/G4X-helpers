import importlib.resources as resources
import logging
import os
import shutil
from pathlib import Path

import polars as pl

from .. import utils
from ..models import G4Xoutput
from ..modules import edit_bin_file, init_bin_file
from ..modules.workflow import workflow
from . import validate
from .constants import files_to_migrate

parquet_col_order = ['TXUID', 'sequence', 'confidence_score', 'x_pixel_coordinate', 'y_pixel_coordinate', 'z_level']


# register schemas
path = resources.files('g4x_helpers.schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont


@workflow
def migrate_g4x_data(
    sample_base: Path,
    sample_id: str,
    logger: logging.Logger,
):
    logger.info('Migrating to latest G4X-data schema')

    backup_size = total_size_gb(files_to_migrate, sample_id=sample_id, base_path=sample_base)
    logger.info(f'Creating backup of {backup_size:.2f} GB before proceeding')

    logger.info('Updating file locations')
    migrate_sample_files(sample_base, sample_id=sample_id)

    logger.info('Validating results...')
    valid_schema = validate.validate_g4x_data(
        path=sample_base, schema_name='base_schema', formats={'sample_id': sample_id}, report=False
    )

    if not valid_schema:
        logger.error('Migration failed to produce correct G4X-data schema.')
    else:
        logger.info('Successfully updated file locations! Checking file structures...')
        if not validate.validate_file_schemas(sample_base):
            logger.info('Some file structures need to be updated.')
            update_files = True
        else:
            update_files = False
            logger.info('File structures are up to date! Migration complete.')

    if update_files:
        parquet_lf, parquet_shema = validate.infer_parquet_schema(sample_base)
        logger.info(f'parquet schema is: {parquet_shema}')

        tx_panel_df, tx_panel_schema = validate.infer_tx_panel_schema(sample_base)
        logger.info(f'tx_panel schema is: {tx_panel_schema}')

        bin_file_schema = validate.infer_bin_schema(sample_base)
        logger.info(f'bin file schema is: {bin_file_schema}')

        adata, adata_schema = validate.infer_adata_schema(sample_base)
        logger.info(f'anndata schema is: {adata_schema}')

        if tx_panel_schema != 'valid' and tx_panel_schema != 'unknown':
            logger.info('Building new transcript_panel.csv from migrated transcript table')
            tx_panel = tx_panel_from_parquet(parquet_lf, tx_panel_df)
            tx_panel.write_csv(sample_base / 'transcript_panel.csv')

        if parquet_shema != 'valid' and parquet_shema != 'unknown':
            logger.info('Bringing parquet file in correct order')
            write_parquet_in_order(parquet_lf, sample_base)

        if adata_schema != 'valid' and adata_schema != 'unknown':
            logger.info('Bringing anndata file in correct order')
            write_adata_in_order(adata, sample_base)

        if bin_file_schema != 'valid' and bin_file_schema != 'unknown':
            smp = G4Xoutput(sample_base)
            init_bin_file(g4x_obj=smp, out_dir=sample_base, logger=logger)
            edit_bin_file(
                g4x_obj=smp, bin_file=sample_base / 'g4x_viewer' / f'{sample_id}_segmentation.bin', logger=logger
            )

        logger.info('Re-validating file structures...')
        if not validate.validate_file_schemas(sample_base):
            raise RuntimeError('Migration failed to produce correct file structures.')
        else:
            logger.info('File structures are up to date! Migration complete.')


@workflow
def restore_backup(
    sample_base: Path,
    sample_id: str,
    logger: logging.Logger,
) -> None:
    backup_dir = sample_base / 'migration_backup'
    if backup_dir.exists():
        if any(backup_dir.iterdir()):
            put_back(sample_base, sample_id)
            shutil.rmtree(backup_dir)
        else:
            logger.info('No backed up files found. Nothing to restore.')

    else:
        logger.info('No backup found. Nothing to restore.')


def total_size_gb(files_dict, sample_id=None, base_path='.'):
    total_bytes = 0

    for key, (src, _, include) in files_dict.items():
        if not include:
            continue

        # Format path if sample_id is used
        if sample_id is not None:
            src = src.format(sample_id=sample_id)

        full_path = os.path.join(base_path, src)

        try:
            total_bytes += os.path.getsize(full_path)
        except FileNotFoundError:
            pass  # skip missing files

    return round(total_bytes / (1024**3), 2)


def migrate_file(sample_base, sample_id, old_path, new_path, backup):
    backup_path = sample_base / 'migration_backup'
    backup_path.mkdir(exist_ok=True, parents=True)

    old_path = old_path.format(sample_id=sample_id)
    new_path = new_path.format(sample_id=sample_id)

    old_path_full = sample_base / old_path.format(sample_id=sample_id)
    new_path_full = sample_base / new_path.format(sample_id=sample_id)

    moved = False
    if old_path_full.exists():
        print('Moving:', old_path_full.relative_to(sample_base), '->', new_path)
        if backup:
            shutil.copy2(old_path_full, backup_path / old_path_full.name)

        shutil.move(old_path_full, new_path_full)
        moved = True
    return moved


def migrate_sample_files(sample_base, sample_id):
    backup_dir = sample_base / 'migration_backup'
    if backup_dir.exists() and any(backup_dir.iterdir()):
        raise RuntimeError(
            'Non-empty backup directory detected. Will not proceed with migration to avoid overwriting existing backups. \n'
            'Please restore from backup via "g4x-helpers migrate --restore" or remove the backup directory before migrating again.'
        )
    else:
        moved_files = {}
        for f, paths in files_to_migrate.items():
            old_path, new_path, backup = paths
            moved_files[f] = migrate_file(sample_base, sample_id, old_path, new_path, backup)


def put_back(sample_base, sample_id):
    for f in files_to_migrate:
        old_path, new_path, has_backup = files_to_migrate[f]

        if has_backup:
            file_name = Path(old_path).name.format(sample_id=sample_id)
            new_path = sample_base / 'migration_backup' / file_name
            if not new_path.exists():
                continue

        migrate_file(sample_base, sample_id, str(new_path), old_path, False)
        if f == 'raw_tx_v2':
            (sample_base / files_to_migrate[f][1]).unlink()


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


def write_adata_in_order(adata, sample_base: Path) -> None:
    adata_path = sample_base / files_to_migrate['feat_mtx'][1]
    cell_meta_path = sample_base / files_to_migrate['cell_meta'][1]

    cols = adata.obs.columns
    first_cols = ['nuclearstain_intensity_mean', 'cytoplasmicstain_intensity_mean']
    prot_cols = first_cols + [c for c in cols if '_intensity_mean' in c and c not in first_cols]
    cell_meta_cols = ['cell_id', 'cell_x', 'cell_y', 'nuclei_area_um', 'nuclei_expanded_area_um']
    last_cols = [c for c in cols if c not in cell_meta_cols and c not in prot_cols]

    new_cols = cell_meta_cols + prot_cols + last_cols

    adata.write_h5ad(adata_path)
    cell_metadata = adata.obs[new_cols].rename(columns={'cell_id': 'label'})
    utils.write_csv_gz(cell_metadata, cell_meta_path.with_suffix(''))
