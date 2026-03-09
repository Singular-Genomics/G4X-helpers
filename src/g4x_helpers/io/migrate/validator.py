import json
from pathlib import Path

import numpy as np
import polars as pl


class DataValidator:
    target_path = '.'
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

    def __init__(self, smp_dir, target_path: str | None = None):
        self.smp_dir = Path(smp_dir)
        self.target_path = Path(target_path) if target_path is not None else Path(self.target_path)
        self.name = type(self).__name__ + 'Validator'

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

        return validation_results


def validation_test(func):
    func._is_validation_test = True
    return func


class SampleMeta(DataValidator):
    target_path = 'run_meta.json'

    SCHEMA = ['run_id']

    @validation_test
    def correct_schema(self):
        with open(self.target_path_resolved, 'r') as f:
            smp_meta = json.load(f)

        return set(self.SCHEMA).issubset(set(smp_meta.keys()))


class Segmentation(DataValidator):
    target_path = 'segmentation/segmentation_mask.npz'

    @validation_test
    def correct_keys(self):
        path = self.target_path_resolved
        return list(np.load(path).keys()) == ['nuclei', 'nuclei_exp']


class QCSummary(DataValidator):
    target_path = 'summary_*.html'


class SampleSheet(DataValidator):
    target_path = 'samplesheet.csv'


class RawFeatures(DataValidator):
    target_path = 'rna/raw_features.parquet'

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
    target_path = 'transcript_panel.csv'

    SCHEMA = ['probe_name', 'gene_name', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path_resolved)
        lf_names = lf.collect_schema().names()

        return set(lf_names).issubset(self.SCHEMA)


class ProteinPanel(DataValidator):
    target_path = 'protein_panel.csv'

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
    target_path = 'protein'

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
