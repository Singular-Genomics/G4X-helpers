import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Literal

import glymur
import numpy as np
import pandas as pd
import polars as pl
import tifffile
from anndata import AnnData, read_h5ad
from matplotlib.pyplot import imread

from . import utils
from .schemas import validate


@dataclass()
class G4Xoutput:
    """
    Container for managing and processing data from a G4X run.

    This class initializes and loads metadata, image dimensions, transcript and protein panels for downstream analysis of G4X output data.
    It provides methods to load images, segmentations, transcript data, and interact with single-cell and spatial analysis pipelines.

    Parameters
    ----------
    data_dir : str
        Path pointing to a G4X sample output directory. This should contain all run-related files including metadata,
        segmentation masks, panels, and feature matrices.
    sample_id : str, optional
        The sample ID to associate with the run. If not provided, it will be inferred from the `data_dir` path.

    Attributes
    ----------
    run_meta : dict
        Metadata dictionary loaded from `run_meta.json`.
    shape : tuple
        Image shape (height, width) as inferred from the segmentation mask.
    transcript_panel_dict : dict
        Mapping of transcript genes to panel types, if the transcript panel is present.
    protein_panel_dict : dict
        Mapping of proteins to panel types, if the protein panel is present.
    genes : list of str
        List of transcript gene names (only if transcript panel exists).
    proteins : list of str
        List of protein names (only if protein panel exists).

    Notes
    -----
    On instantiation, this class performs the following:
      - Loads metadata from `run_meta.json`.
      - Loads the shape of the segmentation mask.
      - Parses transcript and protein panel files (if present).
    """

    data_dir: str
    sample_id: str | None = None

    def __post_init__(self):
        self.data_dir = Path(self.data_dir)

        if self.sample_id is None:
            self.sample_id = self.infer_sample_id()

        with open(self.data_dir / 'run_meta.json', 'r') as f:
            self.run_meta = json.load(f)

        self.set_meta_attrs()

        if self.transcript_panel:
            try:
                transcript_panel = (
                    utils.parse_input_manifest(self.data_dir / 'transcript_panel.csv')
                    .unique('gene_name')
                    .sort('gene_name')
                )
                self.transcript_panel_dict = dict(zip(transcript_panel['gene_name'], transcript_panel['panel_type']))
                self.genes = transcript_panel['gene_name'].to_list()
            except Exception:
                self.transcript_panel_dict = {}
                self.genes = []

        if self.protein_panel:
            protein_panel = pl.read_csv(self.data_dir / 'protein_panel.csv').sort('target')
            self.protein_panel_dict = dict(zip(protein_panel['target'], protein_panel['panel_type']))
            self.proteins = protein_panel['target'].to_list()

    # region dunder
    def __repr__(self):
        machine_num = self.machine.removeprefix('g4-').lstrip('0')
        mac_run_id = f'G{machine_num.zfill(2)}-{self.run_id}'
        gap = 16
        repr_string = f'G4X-data @ {self.data_dir}\n'

        shp = (np.array(self.shape) * 0.3125) / 1000

        repr_string += f'{"Sample ID":<{gap}} - {self.sample_id} of {mac_run_id}, {self.fc}\n'
        repr_string += f'{"imaged area":<{gap}} - ({shp[1]:.2f} x {shp[0]:.2f}) mm\n'
        repr_string += f'{"software version":<{gap}} - {self.software_version}\n\n'

        panels = [
            ('Transcript panel', len(self.genes), 'genes', self.genes)
            if self.includes_transcript
            else (None, 0, '', []),
            ('Protein panel', len(self.proteins), 'proteins', self.proteins)
            if self.includes_protein
            else (None, 0, '', []),
        ]

        # Step 1: compute lengths of "<count> <label>"
        pre_bracket_lengths = [
            len(str(count)) + 2 + len(label)  # e.g., "128 genes"
            for (_, count, label, _) in panels
        ]

        # Step 2: max width to align the `[`
        max_pre = max(pre_bracket_lengths)

        def format_panel(title, count, label, items):
            return f'{title:<{gap}} - {count} {label:<{max_pre - len(str(count)) - 1}}[{", ".join(items[0:5])} ... ]\n'

        if self.includes_transcript:
            repr_string += format_panel(*panels[0])

        if self.includes_protein:
            repr_string += format_panel(*panels[1])

        return repr_string

    # region properties
    @property
    def includes_protein(self) -> bool:
        return self.protein_panel != []

    @property
    def includes_transcript(self) -> bool:
        return self.transcript_panel != []

    @property
    def transcript_table_path(self) -> Path:
        return self.data_dir / 'rna' / 'transcript_table.csv.gz'

    @property
    def feature_table_path(self) -> Path:
        return self.data_dir / 'rna' / 'raw_features.parquet'

    @property
    def segmentation_path(self) -> Path:
        return self.data_dir / 'segmentation' / 'segmentation_mask.npz'

    @property
    def feature_mtx_path(self) -> Path:
        return self.data_dir / 'single_cell_data' / 'feature_matrix.h5'

    @property
    def bead_mask_path(self) -> Path:
        return self.data_dir / 'protein' / 'bead_mask.npz'

    def set_meta_attrs(self):
        static_attrs = [
            'platform',
            'machine',
            'run_id',
            'fc',
            'lane',
            'software',
            'software_version',
            'transcript_panel',
            'protein_panel',
        ]

        for k in static_attrs:
            setattr(self, k, self.run_meta.get(k, None))

        self.shape = glymur.Jp2k(self.data_dir / 'h_and_e' / 'nuclear.jp2').shape

    def infer_sample_id(self) -> str:
        summs = list(self.data_dir.glob('summary_*.html'))

        failure = 'Could not determine sample_id:'
        if len(summs) == 0:
            raise validate.ValidationError(f'{failure} "summary_{{sample_id}}.html" missing')
        elif len(summs) > 1:
            raise validate.ValidationError(f'{failure} Multiple "summary_{{sample_id}}.html" files found')
        else:
            summary_file = summs[0]
            sample_id = summary_file.stem.removeprefix('summary_')
        return sample_id

    def validate(self, details: bool = False) -> None:
        report_style = 'short' if details else False
        _ = validate.validate_g4x_data(
            self.data_dir, schema_name='base_schema', formats={'sample_id': self.sample_id}, report=report_style
        )

    # region methods
    def load_adata(self, *, remove_nontargeting: bool = True, load_clustering: bool = True) -> AnnData:
        adata = read_h5ad(self.feature_mtx_path)

        adata.obs_names = adata.obs['cell_id']
        adata.var_names = adata.var['gene_id']

        adata.obs['sample_id'] = adata.uns['sample_id'] = self.sample_id
        adata.uns['software_version'] = self.software_version

        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()

        if load_clustering:
            df = pd.read_csv(self.data_dir / 'single_cell_data' / 'clustering_umap.csv.gz', index_col=0, header=0)
            adata.obs = adata.obs.merge(df, how='left', left_index=True, right_index=True)
            umap_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
            adata.obsm['X_umap'] = adata.obs[[f'{umap_key}_1', f'{umap_key}_2']].to_numpy(dtype=float)

            # convert clustering columns to categorical
            for col in adata.obs.columns:
                if 'leiden' in col:
                    adata.obs[col] = adata.obs[col].astype('category')

        adata.obs_names = f'{self.sample_id}-' + adata.obs['cell_id'].str.split('-').str[1]
        return adata

    def load_image_by_type(
        self,
        image_type: Literal['protein', 'h_and_e', 'nuclear', 'eosin'],
        *,
        thumbnail: bool = False,
        protein: str | None = None,
        cached: bool = False,
    ) -> np.ndarray:
        if image_type == 'protein':
            if not self.protein_panel:
                print('No protein results.')
                return None
            if protein is None or protein not in self.proteins:
                print(f'{protein} not in protein panel.')
                return None
            pattern = f'{protein}_thumbnail.*' if thumbnail else f'{protein}.*'
            directory = 'protein'
        else:
            pattern_base = {'h_and_e': 'h_and_e', 'nuclear': 'nuclear', 'eosin': 'eosin'}.get(image_type)

            if not pattern_base:
                print(f'Unknown image type: {image_type}')
                return None

            pattern = f'{pattern_base}_thumbnail.*' if thumbnail else f'{pattern_base}.*'
            directory = 'h_and_e'

        if cached:
            return self._load_image_cached(self.data_dir, directory, pattern)
        else:
            return self._load_image(self.data_dir, directory, pattern)

    def load_protein_image(self, protein: str, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('protein', thumbnail=thumbnail, protein=protein, cached=cached)

    def load_he_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('h_and_e', thumbnail=thumbnail, cached=cached)

    def load_nuclear_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('nuclear', thumbnail=thumbnail, cached=cached)

    def load_eosin_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('eosin', thumbnail=thumbnail, cached=cached)

    def load_segmentation(self, expanded: bool = True, key: str | None = None) -> np.ndarray:
        from .modules.segment import try_load_segmentation

        arr = np.load(self.segmentation_path)
        available_keys = list(arr.keys())

        if 'nuclei' and 'nuclei_exp' in available_keys:
            key = 'nuclei_exp' if expanded else 'nuclei'

        return try_load_segmentation(cell_labels=self.segmentation_path, expected_shape=self.shape, labels_key=key)

    def load_bead_mask(self) -> np.ndarray:
        return np.load(self.bead_mask_path)['bead_mask']

    # TODO: this is not used anywhere, consider removing. it's also broken
    def load_feature_table(
        self,
        *,
        return_polars: bool = True,
        lazy: bool = False,
        columns: list[str] | None = None,
        alt_file_path: str | None = None,
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = self.feature_table_path if alt_file_path is None else alt_file_path
        return self._load_table(file_path, return_polars, lazy, columns)

    def load_transcript_table(
        self,
        *,
        return_polars: bool = True,
        lazy: bool = False,
        columns: list[str] | None = None,
        alt_file_path: str | None = None,
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = self.transcript_table_path if alt_file_path is None else alt_file_path
        return self._load_table(file_path, return_polars, lazy, columns)

    def list_content(self, subdir=None):
        if subdir is None:
            subdir = ''

        list_path = self.data_dir / subdir
        output = os.listdir(list_path)

        contents = {'dirs': [], 'files': []}
        for item in output:
            if os.path.isdir(list_path / item):
                contents['dirs'].append(item)
            if os.path.isfile(list_path / item):
                contents['files'].append(item)

        return contents

    # region internal
    @staticmethod
    def _load_image_base(data_dir: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        img_path = next((data_dir / parent_directory).glob(pattern), None)
        if img_path is None:
            raise FileNotFoundError(f'No file matching {pattern} found.')
        if img_path.suffix == '.jp2' or img_path.suffix == '.jpg':
            img = glymur.Jp2k(img_path)[:]
        elif img_path.suffix == '.png':
            img = imread(img_path)
        else:
            img = tifffile.imread(img_path)
        return img

    @staticmethod
    @lru_cache(maxsize=None)
    def _load_image_cached(data_dir: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        return G4Xoutput._load_image_base(data_dir, parent_directory, pattern)

    @staticmethod
    def _load_image(data_dir: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        return G4Xoutput._load_image_base(data_dir, parent_directory, pattern)

    def _load_table(
        self, file_path: str, return_polars: bool = True, lazy: bool = False, columns: list[str] | None = None
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = Path(file_path)
        if lazy:
            if file_path.suffix == '.parquet':
                reads = pl.scan_parquet(file_path)
            else:
                reads = pl.scan_csv(file_path)
        else:
            if file_path.suffix == '.parquet':
                reads = pl.read_parquet(file_path)
            else:
                reads = pl.read_csv(file_path)
        if columns:
            reads = reads.select(columns)
        if not return_polars:
            reads = reads.collect().to_pandas()
        return reads

    def _clear_image_cache(self):
        """Evict all cached images so subsequent calls re-read from disk."""
        self._load_image.cache_clear()
