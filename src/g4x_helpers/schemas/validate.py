import importlib.resources as resources
import re
from pathlib import Path

import anndata as ad
import polars as pl
from google.protobuf.message import DecodeError
from pathschema import validate

from .. import utils
from ..g4x_viewer import CellMasksSchema_pb2 as CellMasksSchema

probe_re = re.compile(utils.PROBE_PATTERN)

# register schemas
path = resources.files('g4x_helpers.schemas')
schemas = {}
for cont in path.iterdir():
    schemas[cont.name.removesuffix('.txt')] = cont


class ValidationError(Exception):
    pass


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


def validate_g4x_data(
    path,
    schema_name: str,
    formats: dict | None = None,
    report: str | None = 'short',
):
    path = Path(path)

    if formats is None:
        formats = {'sample_id': 'unknown_sample'}

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

    if not ok:
        raise ValidationError(
            'Output does not conform to latest G4X-data schema. \n'
            'If your data was generated with a previous version of G4X-software, '
            'you can migrate it to the latest schema with: "g4x-helpers migrate" \n'
            'Please ensure that all files and folder are in their original configurations.'
        )
    return True


def infer_sample_id(data_dir) -> str:
    data_dir = Path(data_dir)
    summs = list(data_dir.glob('summary_*.html'))

    failure = 'Could not determine sample_id:'
    if len(summs) == 0:
        raise ValidationError(f'{failure} "summary_{{sample_id}}.html" missing')
    elif len(summs) > 1:
        raise ValidationError(f'{failure} Multiple "summary_{{sample_id}}.html" files found')
    else:
        summary_file = summs[0]
        sample_id = summary_file.stem.split('_')[-1]
    return sample_id


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


def validate_file_schemas(sample_base, verbose: bool = False) -> bool:
    if verbose:
        print('Validating file schemas...')

    _, parquet_shema = infer_parquet_schema(sample_base)
    if verbose:
        print(f'raw_features.parquet: {parquet_shema}')

    _, tx_panel_schema = infer_tx_panel_schema(sample_base)
    if verbose:
        print(f'transcript_panel.csv: {tx_panel_schema}')

    _, adata_schema = infer_adata_schema(sample_base)
    if verbose:
        print(f'feature_matrix.h5: {adata_schema}')

    bin_file_schema = infer_bin_schema(sample_base)
    if verbose:
        print(f'segmentation .bin file: {bin_file_schema}')

    if parquet_shema == tx_panel_schema == adata_schema == bin_file_schema == 'valid':
        if verbose:
            print('All files conform to latest G4X-data schema!')
        return True
    else:
        # print('Some files do not conform to latest G4X-data schema:')
        return False


def infer_parquet_schema(sample_base):
    raw_features_path = sample_base / 'rna' / 'raw_features.parquet'
    if not raw_features_path.exists():
        raise ValidationError('raw_features.parquet file not found.')
    parquet_lf = pl.scan_parquet(raw_features_path)
    parquet_cols = parquet_lf.collect_schema().names()

    if set(parquet_cols) == {
        'TXUID',
        'sequence',
        'confidence_score',
        'x_pixel_coordinate',
        'y_pixel_coordinate',
        'z_level',
    }:
        parquet_shema = 'valid'
    elif {'cell_id', 'gene_name', 'probe_name', 'demuxed'} <= set(parquet_cols):
        parquet_shema = 'invalid'
    elif {'x_coord_shift', 'y_coord_shift', 'z', 'transcript', 'transcript_condensed'} <= set(parquet_cols):
        parquet_shema = 'need_rename'
    else:
        parquet_shema = 'unknown'

    if parquet_shema == 'need_rename':
        parquet_cols = parquet_lf.collect_schema().names()
        for c in parquet_cols:
            if c in parquet_rename:
                # print(f'Renaming {c} to {parquet_rename[c]}')
                parquet_lf = parquet_lf.rename({c: parquet_rename[c]})
        parquet_shema = 'invalid'

    return parquet_lf, parquet_shema


def infer_bin_schema(sample_base):
    p = sample_base / 'g4x_viewer'
    bin_files = list(p.glob('*.bin'))

    if len(bin_files) > 1:
        print('Multiple bin files found. Cannot infer schema.')
        return 'unknown'
    elif len(bin_files) == 0:
        print('No bin file found. Nothing to update.')
        return 'unknown'

    bin_file = bin_files[0]

    bin_file_schema = 'invalid'
    if bin_file.exists():
        cell_masks = read_bin_file(bin_file)
        if not cell_masks:
            bin_file_schema = 'invalid'
        else:
            bin_file_schema = 'valid'
    else:
        print('No bin file found. Nothing to update.')
    return bin_file_schema


def is_valid_probe(s: str) -> bool:
    return bool(probe_re.match(s))


def infer_tx_panel_schema(sample_base: Path) -> str:
    tx_panel_path = sample_base / 'transcript_panel.csv'
    if not tx_panel_path.exists():
        raise ValidationError('transcript_panel.csv file not found.')
    tx_panel_old = pl.read_csv(tx_panel_path)
    tx_panel_cols = tx_panel_old.columns

    tx_panel_schema = 'unknown'
    if {'target_condensed'} <= set(tx_panel_cols):
        tx_panel_schema = 'invalid'
    elif {'probe_name'} <= set(tx_panel_cols):
        if all(
            [
                is_valid_probe(row['probe_name'])
                for row in tx_panel_old.iter_rows(named=True)
                if row['probe_name'] != 'NOT_CONVERTED'
            ]
        ):
            tx_panel_schema = 'valid'
    return tx_panel_old, tx_panel_schema


def infer_adata_schema(sample_base: Path) -> str:
    adata_path = sample_base / 'single_cell_data' / 'feature_matrix.h5'

    if not adata_path.exists():
        raise ValidationError('feature_matrix.h5 file not found.')

    adata = ad.read_h5ad(adata_path)

    cols = set(adata.obs.columns)

    has_expanded = {'expanded_cell_x', 'expanded_cell_y'} <= cols
    has_centroid = {'centroid-0', 'centroid-1'} <= cols
    has_um_areas = {'nuclei_area_um', 'nuclei_expanded_area_um'} <= cols

    adata_schema = 'invalid'
    if not has_expanded and not has_centroid and has_um_areas:
        adata_schema = 'valid'

    elif has_expanded and not has_um_areas:
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
    elif has_centroid:
        adata.obs = adata.obs.drop(columns=['centroid-0', 'centroid-1'])
    else:
        adata_schema = 'unknown'
        print('could not detect adata version')

    if adata_schema == 'invalid':
        cols = adata.obs.columns
        first_cols = ['nuclearstain_intensity_mean', 'cytoplasmicstain_intensity_mean']
        prot_cols = first_cols + [c for c in cols if '_intensity_mean' in c and c not in first_cols]
        cell_meta_cols = ['cell_id', 'cell_x', 'cell_y', 'nuclei_area_um', 'nuclei_expanded_area_um']
        last_cols = [c for c in cols if c not in cell_meta_cols and c not in prot_cols]

        new_cols = cell_meta_cols + prot_cols + last_cols
        adata.obs = adata.obs[new_cols]

    return adata, adata_schema


def read_bin_file(bin_file: Path) -> CellMasksSchema.CellMasks:
    with open(bin_file, 'rb') as f:
        data = f.read()

    # print('Parsing bin file with current schema.')
    cell_masks = CellMasksSchema.CellMasks()
    try:
        # print('Attempting to parse bin file.')
        cell_masks.ParseFromString(data)
        return cell_masks
    except DecodeError:
        print('Failed to parse bin file with current schema. It may need to be updated.')
        return None
