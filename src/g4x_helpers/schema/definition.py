import json

import anndata as ad
import numpy as np
import polars as pl

from .. import constants as c
from .. import io
from .validator import (
    BaseValidator,
    FileValidator,
    ImgDirectoryValidator,
    TableValidator,
    validation_test,
)


# region root
class SampleG4X(FileValidator):
    PRIMARY = True
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

    @validation_test
    def correct_schema(self):
        smp_meta = self.load()
        return set(self.KEYS).issubset(set(smp_meta.keys()))

    @validation_test
    def is_loadable(self):
        try:
            self.load()
            return True
        except Exception:
            return False

    def _load_method(self):
        with open(self.target_path, 'r') as f:
            smp_meta = json.load(f)
        return smp_meta


class QCSummary(FileValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.SUMMARY


class SampleSheet(TableValidator):
    PRIMARY = True
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
    ]

    def parse(self):
        return io.parse_samplesheet(self.target_path)

    @validation_test
    def correct_keys(self):
        run_section, data_section = self.parse()
        if run_section is None or data_section is None:
            return False

        run_section_keys = run_section['Key']
        data_section_keys = data_section.columns

        run_section_valid = all(key in run_section_keys for key in self.EXPECTED_KEYS_RUN_SECTION)
        data_section_valid = all(key in data_section_keys for key in self.EXPECTED_KEYS_DATA_SECTION)

        return run_section_valid and data_section_valid


class Manifest(TableValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.TX_PANEL

    SCHEMA = {'probe': pl.String}

    def parse(self):
        return io.parse_input_manifest(self.target_path)


# region masks
class Segmentation(FileValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.SEG_MASK
    DEFAULT_KEYS = ['nuclei', 'nuclei_exp']
    _main_key = 'nuclei_exp'

    @property
    def available_keys(self):
        return list(np.load(self.target_path).keys())

    @property
    def main_key(self):
        return self._main_key

    @main_key.setter
    def main_key(self, value):
        self._main_key = value
        self.DEFAULT_KEYS = [value]

    # @validation_test
    # def correct_keys(self):
    #     return set(self.available_keys) == set(self.DEFAULT_KEYS)

    def _load_method(self, key: str | None = None):
        key = self.main_key if key is None else key
        return io.import_segmentation(seg_path=self.target_path, labels_key=key, expected_shape=None, use_cache=True)


class BeadMask(FileValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.BEAD_MASK
    DEFAULT_KEY = 'bead_mask'

    @validation_test
    def correct_keys(self):
        path = self.target_path
        return list(np.load(path).keys()) == [self.DEFAULT_KEY]

    def _load_method(self):
        return np.load(self.target_path)[self.DEFAULT_KEY]


# region rna
class RawFeatures(TableValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.RAW_FEATURES

    SCHEMA = {
        'TXUID': pl.String,
        'sequence': pl.String,
        'confidence_score': pl.Float64,
        'y_pixel_coordinate': pl.Float64,
        'x_pixel_coordinate': pl.Float64,
        'z_level': pl.Float64,
    }


class TxTable(TableValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_TX_TABLE

    SCHEMA = {
        'TXUID': pl.String,
        'confidence_score': pl.Float64,
        'y_pixel_coordinate': pl.Float64,
        'x_pixel_coordinate': pl.Float64,
        'z_level': pl.Float64,
        'gene_id': pl.String,
        # 'probe_name': pl.String,
    }


# region single cell
class AdataH5(BaseValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_FEAT_MTX

    @property
    def has_qc(self) -> bool:
        qc_cols = [
            'n_genes_by_counts',
            'log1p_n_genes_by_counts',
            'total_counts',
            'log1p_total_counts',
            'total_counts_ctrl',
            'log1p_total_counts_ctrl',
            'pct_counts_ctrl',
        ]

        ad_cols = ad.read_h5ad(self.target_path, backed='r').obs.columns

        return set(qc_cols).issubset(set(ad_cols))

    def load(self):
        return ad.read_h5ad(self.target_path)


class CellMetadata(TableValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_CELL_METADATA

    SCHEMA = {
        c.CELL_ID_NAME: pl.String,
        'sample_id': pl.String,
        'tissue_type': pl.String,
        'block': pl.String,
        'seg_source': pl.String,
        c.CELL_COORD_X: pl.String,
        c.CELL_COORD_Y: pl.String,
        c.CELL_AREA_NAME: pl.String,
        c.NUC_STAIN_INTENSITY: pl.String,
        c.CYT_STAIN_INTENSITY: pl.String,
        # c.NUC_AREA_NAME: pl.String, # don't want to require this since it will be missing for custom segmentations
    }


class CellxGene(TableValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_GENE

    SCHEMA = {c.CELL_ID_NAME: pl.String}


class CellxProt(TableValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_PROTEIN

    SCHEMA = {c.CELL_ID_NAME: pl.String}


class ClusteringUmap(TableValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_CLUSTERING_UMAP

    SCHEMA = {
        c.CELL_ID_NAME: pl.String,
        'UMAP1': pl.Float32,
        'UMAP2': pl.Float32,
    }


class Dgex(TableValidator):
    PRIMARY = False
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

    @validation_test
    def has_clusters(self):
        df = self.load(lazy=False)
        max_clusters = df.unique(['leiden_res', 'cluster_id']).group_by(['leiden_res']).agg(pl.len())['len'].max()
        return max_clusters > 2


class SingleCellFolder(BaseValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.SINGLE_CELL_DIR

    # TODO create a factory method to generate these validators
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


# region protein
class ProteinPanel(TableValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.PR_PANEL
    SCHEMA = {'target': pl.String, 'panel_type': pl.String}

    @validation_test
    def folder_present(self):
        folder = ProteinDir(root=self.root)
        return folder.path_exists()

    @property
    def proteins(self):
        return self.load()['target'].to_list()


class ProteinDir(ImgDirectoryValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.PR_DIR
    VALID_IMG_TYPES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]
    EXPECTED_DIRS = {'thumbs'}

    @property
    def FILE_MAP(self):
        suffix = self.img_type
        return {f.name.removesuffix(suffix): [f.name] for f in self.existing_files() if f.name.endswith(suffix)}

    @property
    def panel(self):
        return ProteinPanel(root=self.root)

    @property
    def proteins(self):
        if self.panel.is_valid:
            return self.panel.proteins
        else:
            return []

    @validation_test
    def has_panel(self):
        return self.panel.is_valid

    @validation_test
    def images_match_panel(self):
        return set(self.mapped_files.keys()).issuperset(set(self.proteins))


# region hne
class HnEDir(ImgDirectoryValidator):
    PRIMARY = True
    DEFAULT_TARGET_PATH = c.HE_DIR

    FILE_MAP = {
        c.CYTOPLASMIC_STAIN: [f'{c.CYTOPLASMIC_STAIN}.ome.tiff', f'{c.CYTOPLASMIC_STAIN}.jp2'],
        c.H_AND_E: [f'{c.H_AND_E}.ome.tiff', f'{c.H_AND_E}.jp2'],
        c.NUCLEAR_STAIN: [f'{c.NUCLEAR_STAIN}.ome.tiff', f'{c.NUCLEAR_STAIN}.jp2'],
    }

    EXPECTED_DIRS = {'thumbs'}


# region viewer
class ViewerZarr(BaseValidator):
    PRIMARY = False
    DEFAULT_TARGET_PATH = c.FILE_VIEWER_ZARR
