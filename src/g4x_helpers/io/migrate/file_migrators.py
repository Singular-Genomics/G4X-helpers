import shutil

import polars as pl

from ..input import build_sample_metadata
from .migrator import DataMigrator


class SampleMetaMigrator(DataMigrator):
    name = 'SampleMetaMigrator'
    probes = {
        'current': 'sample.g4x',
        'legacy_version_1': 'run_meta.json',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes, name=self.name)

    def _migrate_legacy_version_1(self):
        _ = build_sample_metadata(self.smp_dir, save_file=True)


class CytoplasmicImageMigrator(DataMigrator):
    name = 'CytoplasmicImageMigrator'
    probes = {
        'current': 'h_and_e/cytoplasmic.jp2',
        'legacy_version_1': 'h_and_e/eosin.jp2',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes, name=self.name)

    def _migrate_legacy_version_1(self):

        self.relocate(self.detected_path_full, self.target_path_full, how='move')


class RawFeaturesMigrator(DataMigrator):
    name = 'RawFeaturesMigrator'
    probes = {
        'current': 'rna/raw_features.parquet',
        'legacy_version_2': 'rna/transcript_table.parquet',
        'legacy_version_1': 'diagnostics/transcript_table.parquet',
    }

    desired_schema = [
        'TXUID',
        'sequence',
        'confidence_score',
        'y_pixel_coordinate',
        'x_pixel_coordinate',
        'z_level',
    ]

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes, name=self.name)

    @property
    def is_valid(self):
        # start with base behavior if useful
        base_valid = super().is_valid
        # add further checks
        lf = pl.scan_parquet(self.detected_path_full)
        lf_schema = lf.collect_schema().names()
        correct_lf_schema = lf_schema == self.desired_schema
        return base_valid and correct_lf_schema

    def _migrate_legacy_version_2(self):

        lf = pl.scan_parquet(self.detected_path_full)

        # bring the columns in the desired order
        lf = lf.select(self.desired_schema)
        lf.sink_parquet(self.target_path_full)

    def _migrate_legacy_version_1(self):
        schema_remap = {
            'sequence_to_demux': 'sequence',
            'meanQS': 'confidence_score',
            'x_coord_shift': 'y_pixel_coordinate',
            'y_coord_shift': 'x_pixel_coordinate',
            'z': 'z_level',
            'transcript': 'probe_name',
            'transcript_condensed': 'gene_name',
        }

        lf = pl.scan_parquet(self.detected_path_full)
        lf_schema = lf.collect_schema().names()
        lf = lf.rename({col: schema_remap[col] for col in schema_remap if col in lf_schema})

        # bring the columns in the desired order
        lf = lf.select(self.desired_schema)
        lf.sink_parquet(self.target_path_full)


class DiagnosticsMigrator(DataMigrator):
    name = 'DiagnosticsMigrator'
    probes = {
        'current': DataMigrator.ABSENT_SENTINEL,
        'legacy_version_1': 'diagnostics',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes, name=self.name)

    def _migrate_legacy_version_1(self):

        legacy_path = self.smp_dir / self.probes['legacy_version_1']
        shutil.rmtree(legacy_path)
