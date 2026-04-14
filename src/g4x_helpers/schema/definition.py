import json

import anndata as ad
import numpy as np
import polars as pl

from .. import c, io
from .validator import BaseValidator, FileValidator, FolderValidator, TableValidator, validation_test


# region root
class SampleG4X(FileValidator):
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

    def _load_method(self):
        with open(self.target_path, 'r') as f:
            smp_meta = json.load(f)
        return smp_meta


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


class Manifest(TableValidator):
    DEFAULT_TARGET_PATH = c.TX_PANEL

    SCHEMA = {'probe': pl.String}

    def parse(self):
        return io.parse_input_manifest(self.target_path)


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
        return folder.path_exists()


# region masks
class Segmentation(FileValidator):
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

    def _load_method(self, key: str | None = None):
        key = self.main_key if key is None else key
        return io.import_segmentation(seg_path=self.target_path, labels_key=key, expected_shape=None, use_cache=True)


class BeadMask(FileValidator):
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
    DEFAULT_TARGET_PATH = c.FILE_TX_TABLE

    SCHEMA = {
        'TXUID': pl.String,
        'confidence_score': pl.Float64,
        'y_pixel_coordinate': pl.Float64,
        'x_pixel_coordinate': pl.Float64,
        'z_level': pl.Float64,
        'probe_name': pl.String,
        'gene_id': pl.String,
    }


# region single cell
class AdataH5(BaseValidator):
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
    DEFAULT_TARGET_PATH = c.FILE_CELL_METADATA

    SCHEMA = {
        c.CELL_ID_NAME: pl.String,
        'sample_id': pl.String,
        'tissue_type': pl.String,
        'block': pl.String,
        'seg_source': pl.String,
        c.CELL_COORD_X: pl.String,
        c.CELL_COORD_Y: pl.String,
        # c.NUC_AREA_NAME: pl.String,
        c.CELL_AREA_NAME: pl.String,
        c.NUC_STAIN_INTENSITY: pl.String,
        c.CYT_STAIN_INTENSITY: pl.String,
    }


class CellxGene(TableValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_GENE

    SCHEMA = {c.CELL_ID_NAME: pl.String}


class CellxProt(TableValidator):
    DEFAULT_TARGET_PATH = c.FILE_CELL_X_PROTEIN

    SCHEMA = {c.CELL_ID_NAME: pl.String}


class ClusteringUmap(TableValidator):
    DEFAULT_TARGET_PATH = c.FILE_CLUSTERING_UMAP

    SCHEMA = {
        c.CELL_ID_NAME: pl.String,
        'UMAP1': pl.Float32,
        'UMAP2': pl.Float32,
    }


class Dgex(TableValidator):
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
class HnEDir(FolderValidator):
    DEFAULT_TARGET_PATH = c.HE_DIR
    EXPECTED_FILES = {
        f'{c.CYTOPLASMIC_STAIN}.ome.tiff',
        f'{c.H_AND_E}.ome.tiff',
        f'{c.NUCLEAR_STAIN}.ome.tiff',
    }

    ALT_FILES = {
        f'{c.CYTOPLASMIC_STAIN}.jp2',
        f'{c.H_AND_E}.jp2',
        f'{c.NUCLEAR_STAIN}.jp2',
    }

    EXPECTED_DIRS = {'thumbs'}

    def get_img(self, query: str):
        search_pool = self.ALT_FILES if self.alt_present else self.EXPECTED_FILES
        fname = next(x for x in search_pool if query in x)
        return self.p / fname


# class HnEDir(BaseValidator):
#     DEFAULT_TARGET_PATH = c.HE_DIR

#     IMG_SUFFIXES = [c.PREFERRED_IMG_SUFFIX, c.ALT_IMG_SUFFIX]

#     IMAGES = [
#         c.NUC_IMG,
#         c.CYT_IMG,
#         c.HNE_IMG,
#     ]

#     def __init__(self, root):
#         super().__init__(root=root)

#         self.existing_files = {}
#         for img in self.IMAGES:
#             img_path = Path(f'{img}__missing__')
#             for suffix in self.IMG_SUFFIXES:
#                 candidate = (self.root / img).with_suffix(suffix)
#                 if candidate.exists():
#                     img_path = candidate

#             self.existing_files[img] = img_path  # is not None

#     @validation_test
#     def images_present(self):
#         existing_files = {}
#         for img_name, img_path in self.existing_files.items():
#             existing_files[img_name] = img_path.exists()
#         return all(existing_files.values())


# region viewer
class ViewerZarr(BaseValidator):
    DEFAULT_TARGET_PATH = c.FILE_VIEWER_ZARR
