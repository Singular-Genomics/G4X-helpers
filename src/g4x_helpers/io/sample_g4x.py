import json

import polars as pl

from g4x_helpers.io.input import _parse_samplesheet
from g4x_helpers.io.validate import SampleMetadata, SampleSheet

from .. import __version__

EXPECTED_KEYS_RUN_META = [
    'machine',
    'run_id',
    'platform',
    'fc',
    'lane',
    'time_of_creation',
    'software',
    'software_version',
]

EXPECTED_KEYS_SS_DATA_SECTION = SampleSheet.EXPECTED_KEYS_DATA_SECTION
EXPECTED_KEYS_SS_RUN_SECTION = SampleSheet.EXPECTED_KEYS_RUN_SECTION


def create_sample_g4x(sample_id: str, run_meta: dict, ssheet: str, out_path: str | None = None) -> dict:

    # create sample_g4x dictionary with default values
    DEFAULT_VALUE = '<unavailable>'
    sample_g4x = {k: DEFAULT_VALUE for k in SampleMetadata.KEYS}
    sample_g4x['sample_id'] = sample_id
    sample_g4x['output_version'] = __version__

    # sanitize run_meta input
    run_meta = _format_keys(run_meta)

    missing_keys = [k for k in EXPECTED_KEYS_RUN_META if k not in run_meta]
    for key in missing_keys:
        print(f"Warning: Expected key '{key}' is missing from run_meta.")

    run_meta_sanitized = {k: v for k, v in run_meta.items() if k in EXPECTED_KEYS_RUN_META and v is not None}

    ### extract samplesheet information into a dictionary
    sample_sheet_info = _extract_sample_sheet_info(sample_id, ssheet)
    sample_sheet_info = _format_keys(sample_sheet_info)
    sample_sheet_info = {k: v for k, v in sample_sheet_info.items() if v is not None}

    ### create a single value for custom panels
    sample_g4x = _handle_custom_panels(sample_g4x, sample_sheet_info)

    ### merge the sanitized run_meta and sample_sheet_info into the sample_g4x dictionary
    sample_g4x = _merge_existing_keys(run_meta_sanitized, sample_g4x)
    sample_g4x = _merge_existing_keys(sample_sheet_info, sample_g4x)

    if out_path is not None:
        with open(out_path, 'w') as f:
            json.dump(sample_g4x, f, indent=4)

    return sample_g4x


def _extract_sample_sheet_info(sample_id: str, ssheet: str) -> dict:
    run_info_section, data_section = _parse_samplesheet(ssheet)

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


def _handle_custom_panels(sample_g4x: dict, sample_sheet_info: dict) -> dict:
    tx_panel = sample_sheet_info.pop('transcript_panel', None)
    tx_panel_cust = sample_sheet_info.pop('transcript_custom', None)
    pr_panel = sample_sheet_info.pop('protein_panel', None)
    pr_panel_cust = sample_sheet_info.pop('protein_custom', None)

    if tx_panel_cust is not None:
        sample_g4x['transcript_panel'] = f'[custom] {tx_panel_cust}'
    else:
        sample_g4x['transcript_panel'] = tx_panel
    if pr_panel_cust is not None:
        sample_g4x['protein_panel'] = f'[custom] {pr_panel_cust}'
    else:
        sample_g4x['protein_panel'] = pr_panel

    return sample_g4x


def _format_keys(dict: dict) -> dict:
    return {k.replace(' ', '_').lower(): v for k, v in dict.items()}


def _merge_existing_keys(from_dict: dict, to_dict: dict) -> dict:
    updated = {k: from_dict.get(k, v) for k, v in to_dict.items()}
    return updated
