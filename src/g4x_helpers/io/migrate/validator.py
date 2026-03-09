from pathlib import Path

import numpy as np
import polars as pl


class DataValidator:
    VALIDATION_TESTS = []

    def __init__(self, smp_dir, target_path: str, name: str | None = None):
        self.name = name or '_unnamed-validator_'
        self.smp_dir = Path(smp_dir)
        self.target_path = target_path

    @property
    def target_path_full(self):
        return self.smp_dir / self.target_path

    @property
    def path_exists(self):
        return self.target_path_full.exists()

    @property
    def is_valid(self):
        return all(self.validation().values())

    def validation(self):
        validation_results = {'path_exists': self.path_exists}

        if self.VALIDATION_TESTS and self.path_exists:
            for test_name in self.VALIDATION_TESTS:
                validation_results[test_name] = getattr(self, test_name)(self.target_path_full)

        return validation_results


class SegmentationValidator(DataValidator):
    target_path = 'segmentation/segmentation_mask.npz'

    def correct_keys(self, path):
        return list(np.load(path).keys()) == ['nuclei', 'nuclei_exp']

    VALIDATION_TESTS = ['correct_keys']

    def __init__(self, smp_dir):
        super().__init__(smp_dir, self.target_path, type(self).__name__)


class QCSummaryValidator(DataValidator):
    target_path = 'summary_*.html'

    def __init__(self, smp_dir):
        matched = list(smp_dir.glob(self.target_path))
        if len(matched) == 1:
            self.target_path = matched[0].relative_to(smp_dir)
        else:
            raise ValueError('Expected exactly one summary HTML file, found {}'.format(len(matched)))

        super().__init__(smp_dir, self.target_path, type(self).__name__)


class SampleSheetValidator(DataValidator):
    target_path = 'samplesheet.csv'

    def __init__(self, smp_dir):
        super().__init__(smp_dir, self.target_path, type(self).__name__)


class RawFeaturesValidator(DataValidator):
    target_path = 'rna/raw_features.parquet'

    SCHEMA = [
        'TXUID',
        'sequence',
        'confidence_score',
        'y_pixel_coordinate',
        'x_pixel_coordinate',
        'z_level',
    ]

    VALIDATION_TESTS = ['correct_schema']

    def correct_schema(self, path):
        lf = pl.scan_parquet(path)
        lf_schema = lf.collect_schema().names()
        return lf_schema == self.SCHEMA

    def __init__(self, smp_dir):
        super().__init__(smp_dir, self.target_path, type(self).__name__)
