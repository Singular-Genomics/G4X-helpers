import json
from pathlib import Path

import numpy as np
import polars as pl

from .. import c, io
from .validator import BaseValidator, validation_test


# region root
class SampleG4X(BaseValidator):
    DEFAULT_TARGET_PATH = c.SMP_META

    KEYS = [
        'run_name',
        'sample_id',
        'tissue_type',
        'block',
        'assay',
        'machine',
        'run_id',
        'fc_layout',
        'fc',
        'lane',
        'sample_position',
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
            with open(self.target_path, 'r') as f:
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

    def load(self):
        if self.is_valid:
            with open(self.target_path, 'r') as f:
                smp_meta = json.load(f)
            return smp_meta
        return self.validation()


class QCSummary(BaseValidator):
    DEFAULT_TARGET_PATH = c.SUMMARY


class SampleSheet(BaseValidator):
    DEFAULT_TARGET_PATH = c.SSHEET

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
            res = io.parse_samplesheet(self.target_path)
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


class Manifest(BaseValidator):
    DEFAULT_TARGET_PATH = c.TX_PANEL

    SCHEMA = ['probe', 'gene_name', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    def parse(self):
        if self.is_valid:
            return io.parse_input_manifest(self.target_path)
        return self.validation()


class ProteinPanel(BaseValidator):
    DEFAULT_TARGET_PATH = c.PR_PANEL

    SCHEMA = ['target', 'panel_type']

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    @validation_test
    def folder_present(self):
        folder = ProteinPanel(self.root)
        return folder.path_exists


# region masks
class Segmentation(BaseValidator):
    DEFAULT_TARGET_PATH = c.SEG_MASK
    DEFAULT_KEYS = ['nuclei', 'nuclei_exp']
    _main_key = 'nuclei_exp'

    @property
    def main_key(self):
        return self._main_key

    @main_key.setter
    def main_key(self, value):
        self._main_key = value
        self.DEFAULT_KEYS = [value]

    @validation_test
    def correct_keys(self):
        path = self.target_path
        return set(np.load(path).keys()) == set(self.DEFAULT_KEYS)

    def load(self, key: str | None = None):
        key = self.main_key if key is None else key
        return io.import_segmentation(seg_path=self.target_path, labels_key=key, expected_shape=None, use_cache=True)


class BeadMask(BaseValidator):
    DEFAULT_TARGET_PATH = c.BEAD_MASK

    @validation_test
    def correct_keys(self):
        path = self.target_path
        return list(np.load(path).keys()) == ['bead_mask']


# region rna
class RawFeatures(BaseValidator):
    DEFAULT_TARGET_PATH = c.RAW_FEATURES

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
        lf = pl.scan_parquet(self.target_path)
        lf_schema = lf.collect_schema().names()
        return lf_schema == self.SCHEMA


class TxTable(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_TX_TABLE

    def load(self, lazy: bool = False):
        if self.is_valid:
            return io.import_table(self.target_path, lazy=lazy)
        else:
            return self.validation()


# region protein
class ProteinDir(BaseValidator):
    DEFAULT_TARGET_PATH = c.PR_DIR

    IMG_SUFFIXES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

    def __init__(self, root):
        super().__init__(root=root)
        self.panel = ProteinPanel(root=self.root)
        proteins = pl.read_csv(self.panel.target_path)['target'].to_list()

        self.existing_files = {}
        for pr in proteins:
            img = None
            for suffix in self.IMG_SUFFIXES:
                candidate = (self.target_path / pr).with_suffix(suffix)
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


# region hne
class HnEDir(BaseValidator):
    DEFAULT_TARGET_PATH = c.HE_DIR

    IMG_SUFFIXES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

    IMAGES = [
        c.NUC_IMG,
        c.CYT_IMG,
        c.HNE_IMG,
    ]

    def __init__(self, root):
        super().__init__(root=root)

        self.existing_files = {}
        for img in self.IMAGES:
            img_path = Path(f'{img}__missing__')
            for suffix in self.IMG_SUFFIXES:
                candidate = (self.root / img).with_suffix(suffix)
                if candidate.exists():
                    img_path = candidate

            self.existing_files[img] = img_path  # is not None

    @validation_test
    def images_present(self):
        existing_files = {}
        for img_name, img_path in self.existing_files.items():
            existing_files[img_name] = img_path.exists()
        return all(existing_files.values())


# region single cell
class CellMetadata(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_METADATA
    SCHEMA = [c.CELL_ID_NAME]

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    def load(self, lazy: bool = False):
        if self.is_valid:
            return io.import_table(self.target_path, lazy=lazy)
        else:
            return self.validation()


class CellxGene(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_GENE
    SCHEMA = [c.CELL_ID_NAME]

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    def load(self, lazy: bool = False):
        if self.is_valid:
            return io.import_table(self.target_path, lazy=lazy)
        else:
            return self.validation()


class CellxProt(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_PROTEIN
    SCHEMA = [c.CELL_ID_NAME]

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA).issubset(lf_names)

    def load(self, lazy: bool = False):
        if self.is_valid:
            return io.import_table(self.target_path, lazy=lazy)
        else:
            return self.validation()


class AdataH5(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_FEAT_MTX


class ClusteringUmap(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_CLUSTERING_UMAP

    def load(self, lazy: bool = False):
        if self.is_valid:
            return io.import_table(self.target_path, lazy=lazy)
        else:
            return self.validation()


class Dgex(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_DGEX

    SCHEMA = {
        'leiden_res': pl.String,
        'cluster_id': pl.String,
        'gene_id': pl.String,
        'score': pl.Float64,
        'logfoldchange': pl.Float64,
        'pval': pl.Float64,
        'pval_adj': pl.Float64,
        'pct_nz_group': pl.Float64,
        'pct_nz_reference': pl.Float64,
    }

    # import_method
    def _try_load(self, lazy: bool = False):
        try:
            res = io.import_table(self.target_path, lazy=lazy)
        except Exception as _:
            return None
        return res

    @validation_test
    def correct_schema(self):
        lf = pl.scan_csv(self.target_path)
        lf_names = lf.collect_schema().names()

        return set(self.SCHEMA.keys()).issubset(lf_names)

    @validation_test
    def is_readable(self):
        return self._try_load() is not None

    @validation_test
    def has_clusters(self):
        df = self._try_load()
        if df is not None:
            max_clusters = df.unique(['leiden_res', 'cluster_id']).group_by(['leiden_res']).agg(pl.len())['len'].max()
        return max_clusters > 2

    def load(self, lazy: bool = False):
        if self.is_valid:
            return self._try_load(lazy=lazy)
        else:
            return self.validation()


class SingleCellFolder(BaseValidator):
    DEFAULT_TARGET_PATH = c.SINGLE_CELL_DIR

    SUB_VALIDATORS = [
        CellMetadata(root='.'),
        CellxGene(root='.'),
        CellxProt(root='.'),
        AdataH5(root='.'),
        Dgex(root='.'),
        ClusteringUmap(root='.'),
    ]

    def __init__(self, root=None, target_path=None):
        super().__init__(root=root, target_path=target_path or self.DEFAULT_TARGET_PATH)
        for val in self.SUB_VALIDATORS:
            val.root = self.root
            setattr(self, val.name, val)

    @property
    def existing_files(self):
        existing_files = {}
        for val in self.SUB_VALIDATORS:
            existing_files[val.name] = val.is_valid
        return existing_files

    @validation_test
    def files_present(self):
        return all(self.existing_files.values())


# region viewer
class ViewerZarr(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_VIEWER_ZARR
