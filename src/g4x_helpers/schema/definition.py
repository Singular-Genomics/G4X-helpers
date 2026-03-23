import json
from pathlib import Path

import numpy as np
import polars as pl

from .. import c
from ..io.input import _parse_samplesheet
from .validator import DataValidator, validation_test


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
