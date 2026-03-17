import json
from pathlib import Path
from typing import Literal

import numpy as np
import polars as pl

from .. import constants as c
from ..io.input import _parse_samplesheet


def kv_line_gap(key, value, gap=2):
    value = '<undefined>' if not value else value
    line = f'{key:<{gap}}'
    line += ' - '
    line += f'{value}'

    return line


def pretty_dict_str(d):
    max_len = max([len(k) for k in d.keys()])
    msg = ''
    for k, v in d.items():
        msg += kv_line_gap(k, v, gap=max_len) + '\n'
    return msg


class FileTree:
    def __init__(self, sample_dir):
        self.smp_dir = Path(sample_dir)

        meta_validator = SampleMetadata(self.smp_dir)

        self.raw_validators = [meta_validator]
        self.secondary_validators = []

        self.tx_detected = False
        self.pr_detected = False

        if meta_validator.is_valid:
            self.detect_assay_type(meta_validator.target_path_resolved)

            self.raw_validators.extend(
                [
                    SampleSheet(self.smp_dir),
                    Segmentation(self.smp_dir),
                    BeadMask(self.smp_dir),
                    HnEDir(self.smp_dir),
                ]
            )

            self.secondary_validators.extend(
                [
                    QCSummary(self.smp_dir),
                    ViewerZarr(self.smp_dir),
                    SingleCellFolder(self.smp_dir),
                    CellMetadata(self.smp_dir),
                    AdataH5(self.smp_dir),
                ]
            )

            if self.tx_detected:
                self.raw_validators.extend(
                    [
                        TranscriptPanel(self.smp_dir),
                        RawFeatures(self.smp_dir),
                    ]
                )
                self.secondary_validators.extend(
                    [
                        TranscriptTable(self.smp_dir),
                        CellxGene(self.smp_dir),
                    ]
                )

            if self.pr_detected:
                self.raw_validators.extend(
                    [
                        ProteinPanel(self.smp_dir),
                        ProteinDir(self.smp_dir),
                    ]
                )
                self.secondary_validators.extend(
                    [
                        CellxProtein(self.smp_dir),
                    ]
                )

        self.validators = self.raw_validators + self.secondary_validators

        for v in self.validators:
            setattr(self, v.name, v)

    def detect_assay_type(self, meta_path):
        with open(meta_path, 'r') as f:
            smp_meta = json.load(f)

        # TODO ensure these hooks are a good choice
        # inspect tx and pr presence from valid metadata
        self.tx_detected = 'transcript_panel' in smp_meta
        self.pr_detected = 'protein_panel' in smp_meta

        self.assay_type = (
            'combined'
            if self.tx_detected and self.pr_detected
            else 'tx_only'
            if self.tx_detected
            else 'pr_only'
            if self.pr_detected
            else 'undefined'
        )

    @property
    def is_valid_raw(self):
        return all([v.is_valid for v in self.raw_validators])

    @property
    def is_valid_all(self):
        return all([v.is_valid for v in self.validators])

    @property
    def errors(self):
        return {v.name: v.validation() for v in self.validators if not v.is_valid}

    @property
    def reports(self):
        return {v.name: v.validation() for v in self.validators}

    def _val_report_minimal(self, raw_only: bool = True, report_pass: bool = True, raise_exception: bool = True):
        gate = self.is_valid_raw if raw_only else self.is_valid_all
        what = 'raw' if raw_only else 'all'
        validation_title = f'G4X-[{what} data] validation'
        if not gate:
            msg = f'{validation_title} failed for:\n'
            msg += f'{self.smp_dir}\n'

            errs = pretty_dict_str(self.errors)
            msg += f'\nErrors:\n{errs}'
            if raise_exception:
                raise ValidationError(msg)
            else:
                print(msg)
        else:
            if report_pass:
                print(f'{validation_title} passed for:\n{self.smp_dir}')

    def _val_report_verbose(self, raw_only: bool = False, raise_exception: bool = True):

        if self.SampleMetadata.path_exists:
            print(f'Detected G4X-metadata file:\n{self.SampleMetadata.target_path_resolved}')
            print(f'assay type: {self.assay_type}')
            print('\n> Validating required raw data ...')
        else:
            self.SampleMetadata.report_validation()

        for v in self.raw_validators:
            v.report_validation()

        if not raw_only:
            print('\n> Validating secondary data ...')
            for v in self.secondary_validators:
                v.report_validation()
        print('')
        self._val_report_minimal(raw_only=raw_only, raise_exception=raise_exception)

    def validation_report(
        self,
        format: Literal['verbose', 'minimal'] = 'verbose',
        raw_only: bool = False,
        report_pass: bool = True,
        raise_exception: bool = True,
    ):
        if format == 'verbose':
            self._val_report_verbose(raw_only=raw_only, raise_exception=raise_exception)
        elif format == 'minimal':
            self._val_report_minimal(raw_only=raw_only, report_pass=report_pass, raise_exception=raise_exception)
        else:
            raise ValueError(f"Invalid format: {format}. Expected 'full' or 'minimal'.")


class ValidationError(Exception):
    pass


class DataValidator:
    DEFAULT_TARGET_PATH = '.'
    TYPE = 'file'
    VALIDATION_TESTS = ()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        tests = []
        for base in reversed(cls.__mro__[1:]):
            tests.extend(getattr(base, 'VALIDATION_TESTS', ()))

        for name, value in cls.__dict__.items():
            if getattr(value, '_is_validation_test', False):
                tests.append(name)

        cls.VALIDATION_TESTS = tuple(dict.fromkeys(tests))

    def __init__(self, smp_dir, target_path: str | None = None, validate_absence: bool = False):
        self.smp_dir = Path(smp_dir)
        self._target_path = Path(target_path) if target_path is not None else Path(type(self).DEFAULT_TARGET_PATH)
        self.name = type(self).__name__  # + 'Validator'
        self.validate_absence = validate_absence
        self.backup_dir = self.smp_dir / 'g4x_helpers' / 'migration_backup'

    @property
    def target_path(self):
        return self._target_path

    @target_path.setter
    def target_path(self, value):
        self._target_path = Path(value)

    @property
    def target_path_resolved(self):
        rel_path = self.target_path if not self.has_wildcard else self.resolve_wildcard()
        return (self.smp_dir / rel_path).resolve()

    @property
    def path(self):
        return self.target_path_resolved

    @property
    def has_wildcard(self):
        parts = self.target_path.parts
        return any('*' in part or '?' in part or '[' in part for part in parts)

    @property
    def path_exists(self):
        return self.target_path_resolved.exists()

    def resolve_wildcard(self):
        if self.has_wildcard:
            matches = sorted(self.smp_dir.glob(str(self.target_path)))

            if len(matches) != 1:
                # print(
                #     f'{self.name}: Could not resolve wildcard, expected exactly 1 match for {self.target_path!s}, found {len(matches)}: {matches}'
                # )
                return self.target_path
            return matches[0].relative_to(self.smp_dir)
        else:
            return self.target_path

    @property
    def is_valid(self):
        return all(self.validation().values())

    def validation(self):
        validation_results = {
            # 'resolved_path': self.target_path,  #
            'path_exists': self.path_exists,
        }

        if self.VALIDATION_TESTS and self.path_exists:
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)()

        if self.validate_absence:
            validation_results['absent'] = not self.path_exists

        return validation_results

    def report_validation(self):
        val_name = self.name.removesuffix('Validator')
        resolved_target_name = self.target_path_resolved.relative_to(self.smp_dir)
        if not self.is_valid:
            code = '[!invalid]'
            msg = f'{code} {val_name} {self.TYPE}: {resolved_target_name}\n'
            msg += f'{len(code) * " "} reason: {self.validation()}'
            print(msg)
        else:
            print(f'[valid] {val_name} {self.TYPE}: {resolved_target_name}')


def validation_test(func):
    func._is_validation_test = True
    return func


class SampleMetadata(DataValidator):
    DEFAULT_TARGET_PATH = c.SMP_META
    TYPE = 'file'

    KEYS = [
        'run_name',
        'sample_position',
        'sample_id',
        'tissue_type',
        'block',
        'assay',
        'machine',
        'run_id',
        'fc_layout',
        'fc',
        'lane',
        'platform',
        'user_name',
        'user_email',
        'run_notes',
        'time_of_creation',
        'transcript_panel',
        'protein_panel',
        'software',
        'software_version',
        'output_version',
    ]

    def _try_load_json(self):
        try:
            with open(self.target_path_resolved, 'r') as f:
                meta = json.load(f)
        except Exception as _:
            return {}
        return meta

    @validation_test
    def is_readable(self):
        return bool(self._try_load_json())

    @validation_test
    def correct_schema(self):
        if self.is_readable():
            smp_meta = self._try_load_json()
            return set(self.KEYS).issubset(set(smp_meta.keys()))
        return False


class Segmentation(DataValidator):
    DEFAULT_TARGET_PATH = c.SEG_MASK
    TYPE = 'file'

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['nuclei', 'nuclei_exp']


class BeadMask(DataValidator):
    DEFAULT_TARGET_PATH = c.BEAD_MASK
    TYPE = 'file'

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['bead_mask']


class QCSummary(DataValidator):
    DEFAULT_TARGET_PATH = c.SUMMARY
    TYPE = 'file'


class SampleSheet(DataValidator):
    DEFAULT_TARGET_PATH = c.SSHEET
    TYPE = 'file'

    EXPECTED_KEYS_RUN_SECTION = [
        'Date',
        'Run Name',
        'User Name',
        'User Email',
        'Workflow',
        'Assay',
        'Run Notes',
        'Stage1',
        'Stage2',
        'FC Layout',
    ]

    EXPECTED_KEYS_DATA_SECTION = [
        'Lane',
        'Sample Position',
        'Tissue Type',
        'Block',
        'Transcript Panel',
        'Protein Panel',
        'Transcript Custom',
        'Protein Custom',
    ]

    def _try_parse_samplesheet(self):
        try:
            res = _parse_samplesheet(self.target_path_resolved)
        except Exception as _:
            return None
        return res

    @validation_test
    def is_parsable(self):
        return bool(self._try_parse_samplesheet())

    @validation_test
    def correct_keys(self):
        run_section, data_section = self._try_parse_samplesheet()
        if run_section is None or data_section is None:
            return False

        run_section_keys = run_section['Key']
        data_section_keys = data_section.columns

        run_section_valid = all(key in run_section_keys for key in self.EXPECTED_KEYS_RUN_SECTION)
        data_section_valid = all(key in data_section_keys for key in self.EXPECTED_KEYS_DATA_SECTION)

        return run_section_valid and data_section_valid


class RawFeatures(DataValidator):
    DEFAULT_TARGET_PATH = c.RAW_FEATURES
    TYPE = 'file'

    SCHEMA = [
        'TXUID',
        'sequence',
        'confidence_score',
        'y_pixel_coordinate',
        'x_pixel_coordinate',
        'z_level',
    ]

    @validation_test
    def correct_schema(self):
        lf = pl.scan_parquet(self.target_path_resolved)
        lf_schema = lf.collect_schema().names()
        return lf_schema == self.SCHEMA


class TranscriptTable(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_TX_TABLE
    TYPE = 'file'


class ViewerZarr(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_VIEWER_ZARR
    TYPE = 'file'


class TranscriptPanel(DataValidator):
    DEFAULT_TARGET_PATH = c.TX_PANEL
    TYPE = 'file'

    SCHEMA = ['probe_name', 'gene_name', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)


class ProteinPanel(DataValidator):
    DEFAULT_TARGET_PATH = c.PR_PANEL
    TYPE = 'file'

    SCHEMA = ['target', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    @validation_test
    def folder_present(self):
        folder = ProteinPanel(self.smp_dir)
        return folder.path_exists


class ProteinDir(DataValidator):
    DEFAULT_TARGET_PATH = c.PR_DIR
    TYPE = 'directory'

    IMG_SUFFIXES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

    def __init__(self, smp_dir):
        super().__init__(smp_dir)
        self.panel = ProteinPanel(self.smp_dir)
        proteins = pl.read_csv(self.panel.target_path_resolved)['target'].to_list()

        self.existing_files = {}
        for pr in proteins:
            img = None
            for suffix in self.IMG_SUFFIXES:
                candidate = (self.target_path_resolved / pr).with_suffix(suffix)
                if candidate.exists():
                    img = candidate
                    break
            self.existing_files[pr] = img  # is not None

    @validation_test
    def has_panel(self):
        return self.panel.is_valid

    @validation_test
    def images_match_panel(self):
        return all(self.existing_files.values())


class HnEDir(DataValidator):
    DEFAULT_TARGET_PATH = c.HE_DIR
    TYPE = 'directory'

    IMG_SUFFIXES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

    IMAGES = [
        c.NUC_IMG,
        c.CYT_IMG,
        c.HNE_IMG,
    ]

    def __init__(self, smp_dir):
        super().__init__(smp_dir)

        self.existing_files = {}
        for img in self.IMAGES:
            img_path = Path(f'{img}__missing__')
            for suffix in self.IMG_SUFFIXES:
                candidate = (self.smp_dir / img).with_suffix(suffix)
                if candidate.exists():
                    img_path = candidate

            self.existing_files[img] = img_path  # is not None

    @validation_test
    def images_present(self):
        existing_files = {}
        for img_name, img_path in self.existing_files.items():
            existing_files[img_name] = img_path.exists()
        return all(existing_files.values())


class CellMetadata(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_METADATA
    TYPE = 'file'


class CellxGene(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_GENE
    TYPE = 'file'


class CellxProtein(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_PROTEIN
    TYPE = 'file'


class AdataH5(DataValidator):
    DEFAULT_TARGET_PATH = c.FILE_FEAT_MTX
    TYPE = 'file'


class SingleCellFolder(DataValidator):
    DEFAULT_TARGET_PATH = c.SINGLE_CELL_DIR
    TYPE = 'directory'

    SUB_VALIDATORS = [
        CellMetadata('.'),
        CellxGene('.'),
        CellxProtein('.'),
        AdataH5('.'),
    ]

    def __init__(self, smp_dir):
        super().__init__(smp_dir, self.DEFAULT_TARGET_PATH)
        for val in self.SUB_VALIDATORS:
            val.smp_dir = self.smp_dir
            setattr(self, val.name, val)

    @validation_test
    def files_present(self):
        existing_files = {}
        for val in self.SUB_VALIDATORS:
            existing_files[val.name] = val.is_valid
        return all(existing_files.values())
