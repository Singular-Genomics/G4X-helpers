import json
from pathlib import Path

import numpy as np
import polars as pl

from .. import constants as c


class FileTree:
    def __init__(self, sample_dir):
        self.smp_dir = sample_dir
        meta_validator = SampleMetadata(sample_dir)

        self.raw_validators = [
            meta_validator,
            SampleSheet(sample_dir),
            Segmentation(sample_dir),
            BeadMask(sample_dir),
            QCSummary(sample_dir),
            HnEFolder(sample_dir),
        ]

        self.tx_detected = False
        self.pr_detected = False
        self.tx_validators = []
        self.pr_validators = []

        # TODO ensure these hooks are a good choice
        if meta_validator.is_valid:
            # inspect tx and pr presence from valid metadata
            with open(meta_validator.target_path_resolved, 'r') as f:
                smp_meta = json.load(f)

            if 'transcript_panel' in smp_meta:
                tx_val = TranscriptPanel(sample_dir)
                rw_val = RawFeatures(sample_dir)
                self.tx_validators = [tx_val, rw_val]
                self.tx_detected = True

            if 'protein_panel' in smp_meta:
                px_val = ProteinPanel(sample_dir)
                pf_val = ProteinFolder(sample_dir)
                self.pr_validators = [px_val, pf_val]
                self.pr_detected = True

        self.raw_validators.extend(self.tx_validators)
        self.raw_validators.extend(self.pr_validators)

        self.secondary_validators = [
            TranscriptTable(sample_dir),
            ViewerZarr(sample_dir),
            SingleCellFolder(sample_dir),
        ]

        self.validators = self.raw_validators + self.secondary_validators

        for v in self.validators:
            setattr(self, v.name, v)

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

    def validation_report_minimal(self, validate_all: bool = False, report_pass: bool = True):
        gate = self.is_valid_raw if not validate_all else self.is_valid_all
        if not gate:
            msg = 'G4X-data validation failed for:\n'
            msg += f'{self.smp_dir}\n'
            msg += f'errors: {self.errors}'
            raise ValidationError(msg)
        else:
            if report_pass:
                print('G4X-data validation passed for:\n', self.smp_dir)

    def validation_report_full(self):

        if self.SampleMetadata.path_exists:
            print('Detected G4X-metadata file: ', self.SampleMetadata.target_path_resolved)
            print('\n> Validating directory ...')
        else:
            self.SampleMetadata.report_validation()

        for v in self.raw_validators:
            v.report_validation()

        if self.tx_detected:
            print('\n> Detected transcript data, validating ...')
            for v in self.tx_validators:
                v.report_validation()

        if self.pr_detected:
            print('\n> Detected protein data, validating ...')
            for v in self.pr_validators:
                v.report_validation()

        print('\n> Validating secondary data ...')
        for v in self.secondary_validators:
            v.report_validation()


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
                raise ValueError(
                    f'{self.name}: Could not resolve wildcard, expected exactly 1 match for {self.target_path!s}, found {len(matches)}: {matches}'
                )
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
            validation_results['required_absent'] = not self.path_exists

        return validation_results

    def report_validation(self):
        val_name = self.name.removesuffix('Validator')
        resolved_target_name = self.target_path_resolved.relative_to(self.smp_dir)
        if not self.is_valid:
            print(f'[!invalid] {val_name} {self.TYPE}: {resolved_target_name} Reason: {self.validation()}')
        else:
            print(f'[valid] {val_name} {self.TYPE}: {resolved_target_name}')


def validation_test(func):
    func._is_validation_test = True
    return func


class SampleMetadata(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SMP_META
    TYPE = 'file'

    SCHEMA = [
        'sample_id',
        'run_name',
        'machine',
        'run_id',
        'platform',
        'fc',
        'lane',
        'sample_position',
        'time_of_creation',
        'software',
        'software_version',
    ]

    @validation_test
    def correct_schema(self):
        with open(self.target_path_resolved, 'r') as f:
            smp_meta = json.load(f)

        return set(self.SCHEMA).issubset(set(smp_meta.keys()))


class Segmentation(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SEG_MASK
    TYPE = 'file'

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['nuclei', 'nuclei_exp']


class BeadMask(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_BEAD_MASK
    TYPE = 'file'

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['bead_mask']


class QCSummary(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SUMMARY
    TYPE = 'file'


class SampleSheet(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SSHEET
    TYPE = 'file'


class RawFeatures(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_RAW_FEATURES
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
    DEFAULT_TARGET_PATH = c.REQUIRED_TX_PANEL
    TYPE = 'file'

    SCHEMA = ['probe_name', 'gene_name', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(lf_names).issubset(self.SCHEMA)


class ProteinPanel(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_PR_PANEL
    TYPE = 'file'

    SCHEMA = ['target', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(lf_names).issubset(self.SCHEMA)

    @validation_test
    def folder_present(self):
        folder = ProteinPanel(self.smp_dir)
        return folder.path_exists


class ProteinFolder(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_PR_DIR
    TYPE = 'directory'

    IMG_SUFFIX = c.REQUIRED_PR_SUFFIX

    def __init__(self, smp_dir):
        super().__init__(smp_dir)
        self.panel = ProteinPanel(self.smp_dir)

    @validation_test
    def has_panel(self):
        return self.panel.is_valid

    @validation_test
    def images_match_panel(self):
        required_proteins = pl.read_csv(self.panel.target_path_resolved)['target'].to_list()

        existing_files = {}
        for pr in required_proteins:
            img = (self.target_path_resolved / pr).with_suffix(self.IMG_SUFFIX)
            existing_files[pr] = img.exists()
        return all(existing_files.values())


class HnEFolder(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_HE_DIR
    TYPE = 'directory'

    REQUIRED_IMAGES = [
        c.REQUIRED_NUC_IMG,
        c.REQUIRED_CYT_IMG,
        c.REQUIRED_HnE_IMG,
    ]

    @validation_test
    def images_present(self):
        existing_files = {}
        for img_path in self.REQUIRED_IMAGES:
            img = self.smp_dir / img_path
            existing_files[img_path] = img.exists()
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


class SingleCellFolder(DataValidator):
    DEFAULT_TARGET_PATH = c.DIRECTORY_SINGLE_CELL
    TYPE = 'directory'

    SUB_VALIDATORS = [
        CellMetadata('.'),
        CellxGene('.'),
        CellxProtein('.'),
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
