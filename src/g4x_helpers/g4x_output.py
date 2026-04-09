import json
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from . import c, io, schema, ut

if TYPE_CHECKING:
    from polars import DataFrame as plDF
    from polars import LazyFrame as plLF

LOGGER = logging.getLogger(__name__)


class G4Xoutput:
    """
    Container for managing and processing data from a G4X run.

    This class initializes and loads metadata, image dimensions, transcript and protein panels for downstream analysis of G4X output data.
    It provides methods to load images, segmentations, transcript data, and interact with single-cell and spatial analysis pipelines.

    """

    def __init__(self, data_dir: str, use_cache: bool = False):
        self.data_dir = Path(data_dir)
        self.src = schema.FileTree(self.data_dir)
        self.out = self.src.copy()
        self.use_cache = use_cache

        self.src.validation_report(format='minimal', raw_only=True, report_pass=False, raise_exception=False)

        with open(self.src.SampleG4X.p, 'r') as f:
            self.smp_meta = json.load(f)

        nuc_img = self.src.HnEDir.existing_files['h_and_e/nuclear']
        self.shape = ut.get_image_shape(nuc_img)

        self.set_meta_attrs()
        self.cache = {}

        self.stains = [c.NUCLEAR_STAIN, c.CYTOPLASMIC_STAIN]
        self.set_genes()
        self.set_proteins()

    ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
    # region dunder
    def __repr__(self):
        machine_num = self.machine.removeprefix('g4-').lstrip('0')
        mac_run_id = f'G{machine_num.zfill(2)}-{self.run_id}'
        gap = 16
        repr_string = f'G4X-data @ {self.data_dir}\n'

        shp = (np.array(self.shape) * 0.3125) / 1000

        repr_string += f'{"Sample ID":<{gap}} - {self.sample_id} of {mac_run_id}, {self.fc}\n'
        repr_string += f'{"tissue, block":<{gap}} - {self.tissue_type}, {self.block}\n'
        repr_string += f'{"imaged area":<{gap}} - ({shp[1]:.2f} x {shp[0]:.2f}) mm\n'
        repr_string += f'{"software version":<{gap}} - {self.software_version}\n\n'

        panels = [
            ('Transcript panel', len(self.genes), 'genes', self.genes) if self.src.tx_detected else (None, 0, '', []),
            ('Protein panel', len(self.proteins), 'proteins', self.proteins)
            if self.src.pr_detected
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

        if self.src.tx_detected:
            repr_string += format_panel(*panels[0])

        if self.src.pr_detected:
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
        ]

        for k in static_attrs:
            setattr(self, k, self.smp_meta.get(k, 'unknown'))

    def set_genes(self, genes: list[str] | None = None):
        self.genes = []
        if genes is not None:
            self.genes = genes
            return

        if self.src.tx_detected:
            tx_panel = self.out.Manifest.load()
            tx_panel = tx_panel.sort(by=['probe_type', 'probe'], descending=[True, False])
            self.genes = tx_panel['gene_name'].unique(maintain_order=True).to_list()

    def set_proteins(self, proteins: list[str] | None = None):
        self.available_proteins = []
        if self.src.pr_detected:
            protein_panel = pl.read_csv(self.src.ProteinPanel.p)
            protein_panel.sort(by=['panel_type', pl.col('target').str.to_lowercase()], descending=[True, False])
            self.available_proteins = self.sort_proteins(protein_panel['target'].to_list())

        self.proteins = []
        if proteins is not None and len(proteins) > 0:
            unavailable_proteins = [p for p in proteins if p not in self.available_proteins]
            if unavailable_proteins:
                raise ValueError(
                    f'The following requested proteins are not available in this data: {unavailable_proteins}'
                )

            self.proteins = proteins
            return

        self.proteins = self.available_proteins.copy()

    @property
    def is_demuxed(self):
        return False if not self.src.tx_detected else self.src.TxTable.is_valid

    @property
    def is_aggregated(self):
        met = self.src.CellMetadata.is_valid
        cxg = True if not self.src.tx_detected else self.src.CellxGene.is_valid
        cxp = True if not self.src.pr_detected else self.src.CellxProt.is_valid
        return met and cxg and cxp

    @property
    def is_scprocessed(self):
        return self.src.SingleCellFolder.is_valid

    @property
    def is_viewer(self):
        return self.src.ViewerZarr.is_valid

    # region methods
    def load_adata(self, *, processed: bool = True, load_clustering: bool = False, remove_nontargeting: bool = False):
        if processed:
            adata = self.src.AdataH5.load()
            leiden_cols = [c for c in adata.obs.columns if c.startswith('leiden')]
            adata.obs = adata.obs.drop(columns=leiden_cols)
        else:
            from .modules.single_cell import init_adata

            adata = init_adata(self)

        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()

        if load_clustering and self.src.ClusteringUmap.is_valid:
            df = self.src.ClusteringUmap.load().cast({'cell_id': pl.Utf8}).to_pandas().set_index('cell_id')
            adata.obs = adata.obs.merge(df, how='left', left_index=True, right_index=True)

            if not processed:
                adata.obsm['X_umap'] = adata.obs[['UMAP1', 'UMAP2']].to_numpy(dtype=np.float32)

        return adata

    def _return_image(
        self,
        img_path: str,
        dask: bool = False,
        shape: tuple[int] | None = None,
        use_cache: bool = False,
        dtype: np.dtype = np.uint16,
    ) -> np.ndarray:
        if dask:
            return io.import_image_dask(img_path=img_path, shape=shape or self.shape, dtype=dtype, use_cache=use_cache)
        else:
            return io.import_image(img_path=img_path, use_cache=use_cache)

    def load_nuclear_image(self, dask: bool = False, use_cache: bool = False) -> np.ndarray:
        img_path = self.src.HnEDir.existing_files['h_and_e/nuclear']
        return self._return_image(img_path=img_path, dask=dask, use_cache=use_cache)

    def load_cytoplasmic_image(self, dask: bool = False, use_cache: bool = False) -> np.ndarray:
        img_path = self.src.HnEDir.existing_files['h_and_e/cytoplasmic']
        return self._return_image(img_path=img_path, dask=dask, use_cache=use_cache)

    def load_he_image(self, dask: bool = False, use_cache: bool = False) -> np.ndarray:
        img_path = self.src.HnEDir.existing_files['h_and_e/h_and_e']
        return self._return_image(
            img_path=img_path, shape=self.shape + (3,), dask=dask, dtype=np.uint8, use_cache=use_cache
        )

    def load_protein_image(self, protein: str, dask: bool = False, use_cache: bool = False) -> np.ndarray:
        img_path = self.src.ProteinDir.existing_files.get(protein)
        if img_path is None:
            print(f'Protein image for {protein} not found.')
            return None

        return self._return_image(img_path=img_path, dask=dask, use_cache=use_cache)

    def load_segmentation(self, expanded: bool = True, key: str = False) -> np.ndarray:
        key = 'nuclei_exp' if expanded else 'nuclei'
        return io.import_segmentation(
            seg_path=self.src.Segmentation.p, expected_shape=self.shape, labels_key=key, use_cache=self.use_cache
        )

    def load_bead_mask(self) -> np.ndarray:
        if self.src.BeadMask.is_valid:
            return np.load(self.src.BeadMask.p)['bead_mask']
        return None

    def load_feature_table(
        self,
        *,
        lazy: bool = False,
        columns: list[str] | None = None,
        use_cache: bool = False,
    ) -> 'plDF | plLF':
        file_path = self.src.RawFeatures.p
        return io.import_table(file_path, lazy=lazy, columns=columns, use_cache=use_cache)

    def load_transcript_table(
        self,
        *,
        lazy: bool = False,
        columns: list[str] | None = None,
        use_cache: bool = False,
    ) -> 'plDF | plLF':
        file_path = self.src.TxTable.p
        return io.import_table(file_path, lazy=lazy, columns=columns, use_cache=use_cache)

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

    def sort_proteins(self, proteins):
        proteins = sorted(proteins, key=str.lower)

        if 'Isotype' in proteins:
            proteins.remove('Isotype')
            proteins = proteins + ['Isotype']
        return proteins
