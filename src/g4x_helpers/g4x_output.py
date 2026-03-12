import json
import os
from pathlib import Path

import glymur
import numpy as np
import pandas as pd
import polars as pl
import tifffile
from anndata import AnnData, read_h5ad

from . import io


class G4Xoutput:
    """
    Container for managing and processing data from a G4X run.

    This class initializes and loads metadata, image dimensions, transcript and protein panels for downstream analysis of G4X output data.
    It provides methods to load images, segmentations, transcript data, and interact with single-cell and spatial analysis pipelines.

    """

    def __init__(self, data_dir: str):
        self.data_dir = Path(data_dir)
        self.tree = io.FileTree(self.data_dir)

        self.tree.validation_report(format='minimal', raw_only=True, report_pass=False, raise_exception=False)

        with open(self.tree.SampleMetadata.path, 'r') as f:
            self.smp_meta = json.load(f)

        nuc_img = list((self.data_dir / 'h_and_e').glob('nuclear.*'))[0]
        if nuc_img.suffix == '.tiff':
            with tifffile.TiffFile(nuc_img) as tif:
                series = tif.series[0]
                self.shape = series.shape
        elif nuc_img.suffix == '.jp2':
            self.shape = glymur.Jp2k(nuc_img).shape
        else:
            raise ValueError(f'Unsupported image format: {nuc_img.suffix}')

        self.set_meta_attrs()

        self.genes = []
        if self.tree.tx_detected:
            tx_panel = io.parse_input_manifest(self.tree.TranscriptPanel.path)
            tx_panel = tx_panel.sort(by=['probe_type', 'probe_name'], descending=[True, False])
            self.genes = tx_panel['gene_name'].unique(maintain_order=True).to_list()

        self.proteins = []
        if self.tree.pr_detected:
            protein_panel = pl.read_csv(self.tree.ProteinPanel.path)
            protein_panel.sort(by=['panel_type', pl.col('target').str.to_lowercase()], descending=[True, False])
            self.proteins = protein_panel['target'].to_list()

    ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
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
            ('Transcript panel', len(self.genes), 'genes', self.genes) if self.tree.tx_detected else (None, 0, '', []),
            ('Protein panel', len(self.proteins), 'proteins', self.proteins)
            if self.tree.pr_detected
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

        if self.tree.tx_detected:
            repr_string += format_panel(*panels[0])

        if self.tree.pr_detected:
            repr_string += format_panel(*panels[1])

        return repr_string

    ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
    def set_meta_attrs(self):
        static_attrs = [
            'platform',
            'machine',
            'run_id',
            'fc',
            'lane',
            'software_version',
            'sample_id',
            'block',
            'tissue_type',
            # 'assay',
            # 'run_name',
        ]

        for k in static_attrs:
            setattr(self, k, self.smp_meta.get(k, 'unknown'))

    @property
    def cell_labels(self) -> np.ndarray:
        nuc_mask = io.import_segmentation(
            seg_path=self.tree.Segmentation.path,
            expected_shape=self.shape,
            labels_key='nuclei',  # TODO figure out what to do with custom segmentations
        )
        nuc_labels = np.unique(nuc_mask)

        return nuc_labels[nuc_labels != 0]

    @property
    def is_demuxed(self):
        result = False if not self.tree.tx_detected else self.tree.TranscriptTable.is_valid
        return result

    @property
    def is_aggregated(self):
        cxg = True if not self.tree.tx_detected else self.tree.CellxGene.is_valid
        cxp = True if not self.tree.pr_detected else self.tree.CellxProtein.is_valid
        return cxg and cxp

    @property
    def is_scprocessed(self):
        # umap = self.tree.clustering_umap_path.exists()
        # meta = self.tree.cell_metadata_path.exists()
        fm = self.tree.AdataH5.is_valid
        return fm
        # return umap and meta and fm

    @property
    def is_viewer(self):
        return self.tree.ViewerZarr.is_valid

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

    def load_protein_image(self, protein: str, use_cache: bool = False) -> np.ndarray:
        img_path = self.tree.ProteinFolder.existing_files.get(protein)
        if img_path is None:
            print(f'Protein image for {protein} not found.')
            return None

        return io.load_image(img_path=img_path, use_cache=use_cache)

    def load_he_image(self, use_cache: bool = False) -> np.ndarray:
        img_path = self.tree.HnEFolder.path / 'h_and_e.ome.tiff'
        return io.load_image(img_path=img_path, use_cache=use_cache)

    def load_nuclear_image(self, use_cache: bool = False) -> np.ndarray:
        img_path = self.tree.HnEFolder.path / 'nuclear.ome.tiff'
        return io.load_image(img_path=img_path, use_cache=use_cache)

    def load_cytoplasmic_image(self, use_cache: bool = False) -> np.ndarray:
        img_path = self.tree.HnEFolder.path / 'cytoplasmic.ome.tiff'
        return io.load_image(img_path=img_path, use_cache=use_cache)

    def load_segmentation(self, expanded: bool = True, key: str | None = None) -> np.ndarray:
        from .modules.segment import try_load_segmentation

        arr = np.load(self.tree.Segmentation.path)
        available_keys = list(arr.keys())

        if 'nuclei' and 'nuclei_exp' in available_keys:
            key = 'nuclei_exp' if expanded else 'nuclei'

        return try_load_segmentation(cell_labels=self.tree.Segmentation.path, expected_shape=self.shape, labels_key=key)

    def load_bead_mask(self) -> np.ndarray:
        return np.load(self.tree.BeadMask.path)['bead_mask']

    # TODO: this is not used anywhere, consider removing. it's also broken
    def load_feature_table(
        self,
        *,
        return_polars: bool = True,
        lazy: bool = False,
        columns: list[str] | None = None,
        alt_file_path: str | None = None,
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = self.tree.RawFeatures.path if alt_file_path is None else alt_file_path
        return self._load_table(file_path, return_polars, lazy, columns)

    def load_transcript_table(
        self,
        *,
        return_polars: bool = True,
        lazy: bool = False,
        columns: list[str] | None = None,
        alt_file_path: str | None = None,
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = self.tree.TranscriptTable.path if alt_file_path is None else alt_file_path
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
