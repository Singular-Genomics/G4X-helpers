import json
from pathlib import Path

import numpy as np
import polars as pl

from .. import constants as c


class ValidationError(Exception):
    pass


class DataValidator:
    DEFAULT_TARGET_PATH = '.'
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
        self.name = type(self).__name__ + 'Validator'
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
        validation_results = {'path_exists': self.path_exists}

        if self.VALIDATION_TESTS and self.path_exists:
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)()

        if self.validate_absence:
            validation_results['required_absent'] = not self.path_exists

        return validation_results

    def report_validation(self):
        val_name = self.name.removesuffix('Validator')
        is_file = self.target_path_resolved.is_file()
        proxy = 'file' if is_file else 'directory'
        if not self.is_valid:
            raise ValidationError(
                f'G4X-helpers requires a {val_name} {proxy}: {self.target_path.name}.\nReason: {self.validation()}'
            )
        else:
            print(f'Found valid {val_name} {proxy}: {self.target_path.name}')


def validation_test(func):
    func._is_validation_test = True
    return func


class SampleMetadata(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SMP_META

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

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['nuclei', 'nuclei_exp']


class QCSummary(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SUMMARY


class SampleSheet(DataValidator):
    DEFAULT_TARGET_PATH = c.REQUIRED_SSHEET


class RawFeatures(DataValidator):
    DEFAULT_TARGET_PATH = 'rna/raw_features.parquet'

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


class TranscriptPanel(DataValidator):
    DEFAULT_TARGET_PATH = 'transcript_panel.csv'

    SCHEMA = ['probe_name', 'gene_name', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(lf_names).issubset(self.SCHEMA)


class ProteinPanel(DataValidator):
    DEFAULT_TARGET_PATH = 'protein_panel.csv'

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
    DEFAULT_TARGET_PATH = 'protein'

    IMG_SUFFIX = '.jp2'

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
    DEFAULT_TARGET_PATH = 'h_and_e'

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
