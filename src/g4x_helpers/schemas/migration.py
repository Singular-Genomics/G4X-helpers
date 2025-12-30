import importlib.resources as resources
import json
import logging
import os
import re
import shutil
from pathlib import Path

import polars as pl

from .. import utils
from ..models import G4Xoutput
from ..modules import create_tx_tarfile, edit_bin_file, init_bin_file
from ..modules.workflow import workflow
from . import validate

parquet_col_order = ['TXUID', 'sequence', 'confidence_score', 'x_pixel_coordinate', 'y_pixel_coordinate', 'z_level']


# register schemas
path = resources.files('g4x_helpers.schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont


class MigrationError(Exception):
    pass


migration_targets = [
    ('g4x_viewer/{sample_id}.bin', 'g4x_viewer/{sample_id}_segmentation.bin', True),
    ('g4x_viewer/{sample_id}.tar', 'g4x_viewer/{sample_id}_transcripts.tar', True),
    ('g4x_viewer/{sample_id}.ome.tiff', 'g4x_viewer/{sample_id}_multiplex.ome.tiff', False),
    ('diagnostics/transcript_table.parquet', 'rna/raw_features.parquet', True),
    ('rna/transcript_table.parquet', 'rna/raw_features.parquet', True),
    ('rna/transcript_table.csv.gz', None, True),
    ('transcript_panel.csv', None, True),
    ('single_cell_data/feature_matrix.h5', None, True),
    ('single_cell_data/cell_metadata.csv.gz', None, True),
    ('g4x_viewer/{sample_id}_run_metadata.json', None, True),
]

backup_loc = 'g4x_helpers/migration_backup'


class MigrationTarget:
    def __init__(self, root: str, smp_id: str, source: str, target: str, backup: bool = True):
        self.root = Path(root)

        target = source if target is None else target
        source = source.format(sample_id=smp_id)
        target = target.format(sample_id=smp_id)

        self.is_rename = True if source != target else False
        self.rename_only = self.is_rename and backup is False

        self.source_file = self.root / source
        self.target_file = self.root / target

        self.backup = backup
        self.backup_dir = self.root / backup_loc
        self.backup_file = self.backup_dir / self.source_file.relative_to(self.root)

        self.source_file_short = self.source_file.relative_to(self.root)
        self.target_file_short = self.target_file.relative_to(self.root)
        self.backup_file_short = self.backup_file.relative_to(self.root)

    def check_source_exists(self) -> bool:
        return self.source_file.exists()

    def check_target_exists(self) -> bool:
        return self.target_file.exists()

    def check_backup_exists(self) -> bool:
        return self.backup_file.exists()

    def backup_source(self):
        if not self.backup:
            return
        if self.backup_file.exists():
            print(f'Backup file "{self.backup_file_short}" already exists. Will not overwrite.')
            return
        if self.source_file.exists():
            self.backup_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(self.source_file, self.backup_file)

    def restore_backup(self):
        if not self.backup:
            if self.check_target_exists() and self.is_rename:
                self.target_file.rename(self.source_file)
            return
        if not self.backup_file.exists():
            print(f'Can not restore file. Backup missing: "{self.backup_file_short}"')
            return
        self.target_file.unlink(missing_ok=False)
        self.backup_file.rename(self.source_file)

    def migrate(self):
        if not self.source_file.exists():
            print(f'Source file "{self.source_file_short}" does not exist. Cannot migrate.')
            return
        self.source_file.rename(self.target_file)

    def __repr__(self):
        return f'MigrationTarget: "{self.target_file_short}"'


@workflow
def migrate_g4x_data(
    data_dir: Path,
    sample_id: str,
    n_threads: int = utils.DEFAULT_THREADS,
    *,
    logger: logging.Logger,
):
    try:
        valid_schema = validate.validate_g4x_data(
            path=data_dir, schema_name='base_schema', formats={'sample_id': sample_id}, report=False
        )
        valid_files = validate.validate_file_schemas(data_dir)
        if not valid_schema or not valid_files:
            raise validate.ValidationError

        logger.info('All validations passed successfully. No migration needed.')
        return
    except Exception:
        logger.info('Migrating to latest G4X-data schema')

    backup_size = total_size_gb(migration_targets, sample_id=sample_id, base_path=data_dir)
    logger.info(f'Creating backup of {backup_size:.2f} GB before proceeding')

    ### STEP 1: Back up and migrate files to new locations
    logger.info('Updating file locations')
    mts = collect_targets(data_dir=data_dir, sample_id=sample_id)
    if any([mt.check_backup_exists() for mt in mts]):
        raise RuntimeError(
            'Non-empty backup directory detected. Will not proceed with migration to avoid overwriting existing backups. \n'
            'Please restore from backup via "g4x-helpers migrate --restore" or remove the backup directory if not needed.'
        )

    mt_sources = [mt for mt in mts if mt.check_source_exists()]

    for mt in mt_sources:
        mt.backup_source()
        logger.info(f'Migrating: {mt.source_file_short} -> {mt.target_file_short}')
        mt.migrate()

    logger.info('Validating results...')
    try:
        valid_schema = validate.validate_g4x_data(
            path=data_dir, schema_name='base_schema', formats={'sample_id': sample_id}, report=False
        )
    except validate.ValidationError:
        msg = (
            'Migration failed to produce correct G4X-data schema\n'
            'Your data was restored to its original state.\n'
            'Please run "g4x-helpers validate" to see which files are blocking validation.\n'
            'Contact care@singulargenomics.com for support'
        )
        restore(data_dir=data_dir, sample_id=sample_id, logger=logger)
        raise MigrationError(msg)

    logger.info('Successfully updated file locations! Checking file structures...')

    ### STEP 2: Validate file structures after migration
    logger.info('Checking for invalid NaN values in run_meta.json')
    fix_json_nan(data_dir / 'g4x_viewer' / f'{sample_id}_run_metadata.json')

    if not validate.validate_file_schemas(data_dir):
        logger.info('Some file structures need to be updated.')
        update_files = True
    else:
        logger.info('File structures are up to date! Migration complete.')
        return

    ### STEP 3: Update file structures as needed
    if update_files:
        parquet_lf, parquet_shema = validate.infer_parquet_schema(data_dir)
        logger.info(f'parquet schema is: {parquet_shema}')

        tx_panel_df, tx_panel_schema = validate.infer_tx_panel_schema(data_dir)
        logger.info(f'tx_panel schema is: {tx_panel_schema}')

        bin_file_schema = validate.infer_bin_schema(data_dir)
        logger.info(f'bin file schema is: {bin_file_schema}')

        adata, adata_schema = validate.infer_adata_schema(data_dir)
        logger.info(f'anndata schema is: {adata_schema}')

        if tx_panel_schema != 'valid' and tx_panel_schema != 'unknown':
            logger.info('Building new transcript_panel.csv from migrated transcript table')
            tx_panel = tx_panel_from_parquet(parquet_lf, tx_panel_df)
            tx_panel.write_csv(data_dir / 'transcript_panel.csv')

        if parquet_shema != 'valid' and parquet_shema != 'unknown':
            logger.info('Bringing parquet file in correct order')
            write_parquet_in_order(parquet_lf, data_dir)

        if adata_schema != 'valid' and adata_schema != 'unknown':
            logger.info('Bringing anndata file in correct order')
            write_adata_in_order(adata, data_dir)

        if bin_file_schema != 'valid' and bin_file_schema != 'unknown':
            smp = G4Xoutput(data_dir)
            init_bin_file(g4x_obj=smp, out_dir=data_dir, n_threads=n_threads, logger=logger)
            edit_bin_file(
                g4x_obj=smp, bin_file=data_dir / 'g4x_viewer' / f'{sample_id}_segmentation.bin', logger=logger
            )
            logger.info('Re-created bin file with updated schema.')
            create_tx_tarfile(
                g4x_obj=smp,
                out_path=data_dir / 'g4x_viewer' / f'{sample_id}_transcripts.tar',
                n_threads=n_threads,
                logger=logger,
            )

        logger.info('Re-validating file structures...')
        if not validate.validate_file_schemas(data_dir, verbose=True):
            msg = (
                'Migration failed to produce correct file structures.\n'
                'Your data was restored to its original state.\n'
                'Contact care@singulargenomics.com for support'
            )

            restore(data_dir=data_dir, sample_id=sample_id, logger=logger)
            raise MigrationError(msg)
        else:
            logger.info('File structures are up to date! Migration complete.')
            return True


@workflow
def restore_backup(
    data_dir: Path,
    sample_id: str,
    logger: logging.Logger,
):
    backup_dir = Path(data_dir) / backup_loc
    if backup_dir.exists():
        mts = collect_targets(data_dir=data_dir, sample_id=sample_id)
        mt_backups = [mt for mt in mts if mt.check_backup_exists() or mt.rename_only]

        if len(mt_backups) > 0:
            for mt in mt_backups:
                logger.info(f'Restoring: {mt.source_file_short}')
                mt.restore_backup()
        else:
            logger.info('No backed up files found. Nothing to restore.')
        shutil.rmtree(Path(data_dir) / backup_loc)
    else:
        logger.info('No backup found. Nothing to restore.')


def restore(
    data_dir: Path,
    sample_id: str,
    logger: logging.Logger,
):
    mts = collect_targets(data_dir=data_dir, sample_id=sample_id)
    mt_backups = [mt for mt in mts if mt.check_backup_exists() or mt.rename_only]

    if len(mt_backups) > 0:
        for mt in mt_backups:
            logger.info(f'Restoring: {mt.source_file_short}')
            mt.restore_backup()
    else:
        logger.info('No backed up files found. Nothing to restore.')
    shutil.rmtree(Path(data_dir) / backup_loc)


def collect_targets(data_dir: Path, sample_id: str):
    return [
        MigrationTarget(root=data_dir, smp_id=sample_id, source=mt[0], target=mt[1], backup=mt[2])
        for mt in migration_targets
    ]


def total_size_gb(files, sample_id=None, base_path='.'):
    total_bytes = 0

    for src, _, include in files:
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


def fix_json_nan(path: Path) -> bool:
    with open(path) as f:
        txt = f.read()

    has_nan = bool(re.search(r'(?<!")\bNaN\b(?!")', txt))
    if has_nan:
        print('Found bare NaN. fixing...')

        # Replace bare NaN with "NaN"
        txt_fixed = txt.replace(' NaN', ' "NaN"').replace(':NaN', ':"NaN"')

        # Parse JSON safely now
        data = json.loads(txt_fixed)

        with open(path, 'w') as f:
            json.dump(data, f, indent=2)


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
    raw_features_path = sample_base / 'rna' / 'raw_features.parquet'
    tmp_file = raw_features_path.with_name('tmp.parquet')
    parquet_lf.select(parquet_col_order).sink_parquet(tmp_file)
    shutil.move(tmp_file, raw_features_path)


def write_adata_in_order(adata, sample_base: Path) -> None:
    adata_path = sample_base / 'single_cell_data' / 'feature_matrix.h5'
    cell_meta_path = sample_base / 'single_cell_data' / 'cell_metadata.csv.gz'

    cols = adata.obs.columns
    first_cols = ['nuclearstain_intensity_mean', 'cytoplasmicstain_intensity_mean']
    prot_cols = first_cols + [c for c in cols if '_intensity_mean' in c and c not in first_cols]
    cell_meta_cols = ['cell_id', 'cell_x', 'cell_y', 'nuclei_area_um', 'nuclei_expanded_area_um']
    last_cols = [c for c in cols if c not in cell_meta_cols and c not in prot_cols]

    new_cols = cell_meta_cols + prot_cols + last_cols

    adata.write_h5ad(adata_path)
    cell_metadata = adata.obs[new_cols].rename(columns={'cell_id': 'label'})
    utils.write_csv_gz(cell_metadata, cell_meta_path.with_suffix(''))
