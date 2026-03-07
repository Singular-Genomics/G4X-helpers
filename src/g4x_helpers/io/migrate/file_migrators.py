import json
import shutil

import polars as pl

from ..input import build_sample_metadata
from .migrator import LegacyFileMigrator


class SampleMetaMigrator(LegacyFileMigrator):
    probes = {
        'current': 'sample.g4x',
        'legacy_version_1': 'run_meta.json',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes)

    def _migrate_legacy_version_1(self, keep_detected: bool = False):
        data = build_sample_metadata(self.smp_dir)

        with open(self.smp_dir / 'sample.g4x', 'w') as f:
            json.dump(data, f, indent=4)

        legacy_file = self.smp_dir / self.probes['legacy_version_1']
        if not keep_detected:
            legacy_file.unlink()


class CytoplasmicImageMigrator(LegacyFileMigrator):
    probes = {
        'current': 'h_and_e/cytoplasmic.jp2',
        'legacy_version_1': 'h_and_e/eosin.jp2',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes)

    def _migrate_legacy_version_1(self, keep_detected: bool = False):
        how = 'copy' if keep_detected else 'move'
        self.copy_or_move(self.detected_file_full, self.target_file_full, how)


class RawFeaturesMigrator(LegacyFileMigrator):
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
        super().__init__(sample_dir, self.probes)

    def _migrate_legacy_version_2(self, keep_detected: bool = False):

        lf = pl.scan_parquet(self.detected_file_full)

        # bring the columns in the desired order
        lf = lf.select(self.desired_schema)
        lf.sink_parquet(self.target_file_full)

        legacy_file = self.smp_dir / self.probes['legacy_version_2']
        if not keep_detected:
            legacy_file.unlink()

    def _migrate_legacy_version_1(self, keep_detected: bool = False):
        schema_remap = {
            'sequence_to_demux': 'sequence',
            'meanQS': 'confidence_score',
            'x_coord_shift': 'y_pixel_coordinate',
            'y_coord_shift': 'x_pixel_coordinate',
            'z': 'z_level',
            'transcript': 'probe_name',
            'transcript_condensed': 'gene_name',
        }

        lf = pl.scan_parquet(self.detected_file_full)
        lf_schema = lf.collect_schema().names()
        lf = lf.rename({col: schema_remap[col] for col in schema_remap if col in lf_schema})

        # bring the columns in the desired order
        lf = lf.select(self.desired_schema)
        lf.sink_parquet(self.target_file_full)

        legacy_file = self.smp_dir / self.probes['legacy_version_1']
        if not keep_detected:
            legacy_file.unlink()


class DiagnosticsMigrator(LegacyFileMigrator):
    probes = {
        'current': '__ABSENT__',
        'legacy_version_1': 'diagnostics',
    }

    def __init__(self, sample_dir):
        super().__init__(sample_dir, self.probes)

    def _migrate_legacy_version_1(self, **kwargs):
        legacy_path = self.smp_dir / self.probes['legacy_version_1']
        shutil.rmtree(legacy_path)
