import json

import polars as pl

from .. import __version__
from ..schema.definition import SampleG4X, SampleSheet
from .input import parse_samplesheet

EXPECTED_KEYS_RUN_META = [
    'machine',
    'run_id',
    'platform',
    'fc',
    'lane',
    'time_of_creation',
    'software',
    'software_version',
    'transcript_panel',
    'protein_panel',
]

EXPECTED_KEYS_SS_DATA_SECTION = SampleSheet.EXPECTED_KEYS_DATA_SECTION
EXPECTED_KEYS_SS_RUN_SECTION = SampleSheet.EXPECTED_KEYS_RUN_SECTION


def create_sample_g4x(sample_id: str, run_meta: dict, ssheet: str, out_path: str | None = None) -> dict:

    # create sample_g4x dictionary with default values
    DEFAULT_VALUE = '<not-provided>'
    sample_g4x = {k: DEFAULT_VALUE for k in SampleG4X.KEYS}
    sample_g4x['sample_id'] = sample_id
    sample_g4x['output_version'] = __version__

    # sanitize run_meta input
    run_meta = _format_keys(run_meta)

    missing_keys = [k for k in EXPECTED_KEYS_RUN_META if k not in run_meta]
    for key in missing_keys:
        print(f"Warning: Expected key '{key}' is missing from run_meta.")

    run_meta_sanitized = {k: v for k, v in run_meta.items() if k in EXPECTED_KEYS_RUN_META and v is not None}

    drop_panels_ssheet = False
    if 'transcript_panel' in run_meta_sanitized or 'protein_panel' in run_meta_sanitized:
        drop_panels_ssheet = True

    ### extract samplesheet information into a dictionary
    sample_sheet_info = _extract_sample_sheet_info(sample_id, ssheet)
    sample_sheet_info = _format_keys(sample_sheet_info)
    sample_sheet_info = {k: v for k, v in sample_sheet_info.items() if v is not None}

    ### create a single value for custom panels
    sample_sheet_info = _handle_ssheet_panel_names(sample_sheet_info)

    if drop_panels_ssheet:
        # if the panels were already provided in run_meta, drop them from the sample_sheet_info to avoid conflicts
        sample_sheet_info.pop('transcript_panel', None)
        sample_sheet_info.pop('protein_panel', None)

    ### merge the sanitized run_meta and sample_sheet_info into the sample_g4x dictionary
    sample_g4x = _merge_existing_keys(run_meta_sanitized, sample_g4x)
    sample_g4x = _merge_existing_keys(sample_sheet_info, sample_g4x)

    sample_g4x = _format_custom_panels(sample_g4x)

    if out_path is not None:
        with open(out_path, 'w') as f:
            json.dump(sample_g4x, f, indent=4)

    return sample_g4x


def _extract_sample_sheet_info(sample_id: str, ssheet: str) -> dict:
    run_info_section, data_section = parse_samplesheet(ssheet)

    # sanitize columns and filter the data_section for the specific sample_id
    def format_keys(key):
        return key.replace(' ', '_').replace('-', '_').lower()

    exp_lower = [format_keys(c) for c in EXPECTED_KEYS_SS_DATA_SECTION]
    data_section = data_section.rename({c: format_keys(c) for c in data_section.columns})

    missing_keys = [c for c in exp_lower if c not in data_section.columns]
    for key in missing_keys:
        data_section = data_section.with_columns(pl.lit(None).alias(key))

    data_section = data_section.select(exp_lower)

    data_section_dict = data_section.filter(
        pl.col('lane') == int(sample_id[-1]), pl.col('sample_position') == sample_id[0]
    ).to_dicts()[0]

    # sanitize run_info_section
    ri_keys = EXPECTED_KEYS_SS_RUN_SECTION
    filtered = [d for d in run_info_section.to_dicts() if d['Key'] in ri_keys]
    run_info = {d['Key']: d['Value'] for d in filtered}

    ### merge the fields that were extracted from the samplesheet
    return {**run_info, **data_section_dict}


def _handle_ssheet_panel_names(sample_sheet_info: dict) -> dict:
    tx_panel = sample_sheet_info.pop('transcript_panel', None)
    tx_panel_cust = sample_sheet_info.pop('transcript_custom', None)
    pr_panel = sample_sheet_info.pop('protein_panel', None)
    pr_panel_cust = sample_sheet_info.pop('protein_custom', None)

    def _parse_ssheet_columns(panel, cust_panel):
        result = {}
        if cust_panel is not None:
            result['name'] = cust_panel
            result['is_custom'] = True
        else:
            result['name'] = panel
            result['is_custom'] = False
        return result

    sample_sheet_info['transcript_panel'] = _parse_ssheet_columns(tx_panel, tx_panel_cust)
    sample_sheet_info['protein_panel'] = _parse_ssheet_columns(pr_panel, pr_panel_cust)
    return sample_sheet_info


def _format_custom_panels(sample_g4x: dict) -> dict:

    def _format(panel):
        prefix = '[custom] ' if panel.get('is_custom', False) else ''
        name = f'{prefix}{panel["name"]}'
        return name

    for key in ['transcript_panel', 'protein_panel']:
        panel = sample_g4x.get(key, None)
        if panel is not None and isinstance(panel, dict) and 'name' in panel:
            sample_g4x[key] = _format(sample_g4x[key])

    return sample_g4x


def _format_keys(dict: dict) -> dict:
    return {k.replace(' ', '_').lower(): v for k, v in dict.items()}


def _merge_existing_keys(from_dict: dict, to_dict: dict) -> dict:
    updated = {k: from_dict.get(k, v) for k, v in to_dict.items()}
    return updated
