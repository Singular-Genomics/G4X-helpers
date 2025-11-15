import importlib.resources as resources
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from anndata import read_h5ad
from pathschema import validate

if TYPE_CHECKING:
    from ..models import G4Xoutput

# register schemas
path = resources.files('g4x_helpers.schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont

col_rename = {
    'x_coord_shift': 'y_pixel_coordinate',
    'y_coord_shift': 'x_pixel_coordinate',
    'z': 'z_level',
    'demuxed': 'demuxed',
    'transcript': 'probe_name',
    'transcript_condensed': 'gene_name',
    'meanQS': 'confidence_score',
    'cell_id': 'cell_id',
    'sequence_to_demux': 'sequence',
    'TXUID': 'TXUID',
}

new_feature_table = 'rna/transcript_table_new.parquet'
new_feature_matrix = 'single_cell_data/feature_matrix_new.h5'
new_cell_metadata = 'single_cell_data/cell_metadata_new.csv.gz'
new_transcript_panel = 'transcript_panel_new.csv'


def check_and_convert_sample(smp: 'G4Xoutput'):
    print('Checking schema version for sample:')
    schema = infer_schema_version(smp)
    if schema == 'invalid':
        print('Schema is invalid. Cannot proceed with this sample.')
    elif schema == 'v2':
        print('Detected schema version v2. Nothing to do.')
    elif schema == 'v1':
        if confirm_v1_schema(smp):
            print('Confirmed schema version v1. Proceeding with conversion to v2...')
            convert_v1_to_v2(smp, delete_old=False)


def infer_schema_version(smp: 'G4Xoutput'):
    is_v1 = validate_g4x_data(path=smp.data_dir, schema_name='v1', report=False, formats={'sample_id': smp.sample_id})
    is_v2 = validate_g4x_data(path=smp.data_dir, schema_name='v2', report=False, formats={'sample_id': smp.sample_id})

    if is_v1 and is_v2:
        result = 'ambiguous'
    elif is_v1 and not is_v2:
        result = 'v1'
    elif is_v2 and not is_v1:
        result = 'v2'
    elif not is_v1 and not is_v2:
        result = 'invalid'

    return result


def confirm_v1_schema(smp: 'G4Xoutput') -> bool:
    old_feature_table = smp.data_dir / 'diagnostics' / 'transcript_table.parquet'
    old_feature_matrix = smp.data_dir / 'single_cell_data' / 'feature_matrix.h5'

    cols = pl.scan_parquet(old_feature_table).collect_schema().names()

    parquet_v1 = adata_v1 = False
    if all(k in cols for k in ['x_coord_shift', 'y_coord_shift', 'transcript_condensed']):
        parquet_v1 = True

    cols = read_h5ad(old_feature_matrix).obs.columns
    if all(k in cols for k in ['expanded_cell_x', 'expanded_cell_y', 'nuclei_expanded_area']):
        adata_v1 = True

    if adata_v1 and parquet_v1:
        return True

    return False


def convert_v1_to_v2(smp: 'G4Xoutput', delete_old: bool = False):
    new_feature_table = smp.data_dir / 'rna/transcript_table_new.parquet'
    new_feature_matrix = smp.data_dir / 'single_cell_data/feature_matrix_new.h5'
    new_cell_metadata = smp.data_dir / 'single_cell_data/cell_metadata_new.csv.gz'

    print('Migrating transcript_table.parquet schema: v1 -> v2')

    (
        pl.scan_parquet(smp.data_dir / 'diagnostics' / 'transcript_table.parquet')
        .rename(col_rename)
        .sink_parquet(new_feature_table)
    )

    print('Building new transcript panel csv from migrated transcript table')
    df = (
        pl.read_parquet(new_feature_table)
        .filter(pl.col('probe_name') != 'UNDETERMINED')
        .unique('probe_name')
        .sort('probe_name')
    )
    df.select('probe_name', 'gene_name').write_csv(smp.data_dir / new_transcript_panel)

    print('Migrating feature_matrix.h5 and cell_metadata.csv.gz schema: v1 -> v2')
    adata = smp.load_adata(remove_nontargeting=False, load_clustering=False)

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

    adata.write_h5ad(new_feature_matrix)
    adata.obs.to_csv(new_cell_metadata)

    if delete_old:
        (smp.data_dir / 'diagnostics' / 'transcript_table.parquet').unlink()


def _print_details(path, result, errors_only=True):
    gap = 21
    for err_path, errors in result.errors_by_path.items():
        err_path = Path(err_path)

        relative = err_path.relative_to(path)

        if not errors_only:
            if len(errors) == 0:
                relative = 'root' if str(relative) == '.' else relative
                print(f'{"path valid":<{gap}} - {relative}')

        for err in errors:
            if err.startswith('Missing required file'):
                err_full = err
                err = 'Missing required file'
                relative = relative / err_full.removeprefix('Missing required file ')

            print(f'{err:<{gap}} - {relative}')


def validate_g4x_data(
    path,
    schema_name: str,
    formats: dict | None = None,
    report: str | None = 'short',
):
    path = Path(path)

    if formats is None:
        formats = {'sample_id': 'A01'}

    with open(schemas[schema_name], 'r') as f:
        schema = f.read()

    schema = schema.format(**formats)
    result = validate(path, schema)

    ok = not result.has_error()

    # store callables, not results
    reports = {
        True: {
            'short': lambda: _print_details(path, result, errors_only=True),
            'long': lambda: _print_details(path, result, errors_only=False),
        },
        False: {
            'short': lambda: _print_details(path, result, errors_only=True),
            'long': lambda: _print_details(path, result, errors_only=False),
        },
    }

    if report:
        try:
            reports[ok][report]()  # <-- call the function
        except KeyError:
            raise ValueError(f'Unknown report type: {report!r}')

    return ok
