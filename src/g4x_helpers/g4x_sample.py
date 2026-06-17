import json
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from . import constants as c
from . import io, schema
from . import utils as ut

if TYPE_CHECKING:
    from anndata import AnnData

    from .schema.validator import BaseValidator

log = logging.getLogger(__name__)


class G4Xsample:
    """
    Container for managing and processing data from a G4X run.

    This class initializes and loads metadata, image dimensions, transcript and protein panels for downstream analysis of G4X output data.
    It provides methods to load images, segmentations, transcript data, and interact with single-cell and spatial analysis pipelines.

    """

    def __init__(
        self, smp_dir: str, alt_source: str | None = None, use_cache: bool = False, validate: bool = True
    ) -> None:
        self.smp_dir = Path(smp_dir)
        self.alt_source = Path(alt_source) if alt_source is not None else None
        self.src = schema.FileTree(self.smp_dir, alt_source=alt_source)
        self.use_cache = use_cache

        if self.alt_source:
            self.alt = schema.FlatTree(self.alt_source)
        else:
            self.alt = self.src

        if validate:
            self.src.validation_report(format='minimal', raw_only=True, report_pass=False, raise_exception=True)

        self.set_meta_attrs()
        self.cache = {}

        self.stains = [c.NUCLEAR_STAIN, c.CYTOPLASMIC_STAIN]

    ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ### --- ###
    # region dunder
    def __repr__(self):
        machine_num = self.machine.removeprefix('g4-').lstrip('0')
        mac_run_id = f'G{machine_num.zfill(2)}-{self.run_id}'
        gap = 16
        repr_string = f'G4X-data @ {self.smp_dir}\n'

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
    def set_meta_attrs(self) -> None:
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

    @property
    def smp_meta(self) -> dict:
        with open(self.src.SampleG4X.p, 'r') as f:
            smp_meta = json.load(f)
        return smp_meta

    @property
    def shape(self) -> tuple[int, int]:
        # TODO add load() method for images
        nuc_img = self.src.HnEDir.get_img(c.NUCLEAR_STAIN)
        return ut.get_image_shape(nuc_img)

    @property
    def genes(self) -> list[str]:
        if self.src.tx_detected:
            tx_panel = self.src.Manifest.parse()
            tx_panel = tx_panel.sort(by=['probe_type', 'gene_name'], descending=[True, False])
            genes = tx_panel['gene_name'].unique(maintain_order=True).to_list()
            return genes
        return []

    @property
    def proteins(self) -> list[str]:

        def _sort_proteins(proteins: list[str]) -> list[str]:
            proteins = sorted(proteins, key=str.lower)

            if 'Isotype' in proteins:
                proteins.remove('Isotype')
                proteins = proteins + ['Isotype']
            return proteins

        if self.src.pr_detected:
            protein_panel = pl.read_csv(self.src.ProteinPanel.p)
            protein_panel.sort(by=['panel_type', pl.col('target').str.to_lowercase()], descending=[True, False])
            return _sort_proteins(protein_panel['target'].to_list())
        return []

    @property
    def uses_branch(self) -> bool:
        return self.alt_source is not None and self.alt_source != self.smp_dir

    @property
    def is_demuxed(self) -> bool:
        return False if not self.src.tx_detected else self.alt.TxTable.is_valid

    @property
    def is_aggregated(self) -> bool:
        met = self.alt.CellMetadata.is_valid
        cxg = True if not self.src.tx_detected else self.alt.CellxGene.is_valid
        cxp = True if not self.src.pr_detected else self.alt.CellxProt.is_valid
        return met and cxg and cxp

    @property
    def is_scprocessed(self) -> bool:
        return self.src.SingleCellFolder.is_valid

    @property
    def is_viewer(self) -> bool:
        return self.src.ViewerZarr.is_valid

    # region methods
    def load_adata(
        self, *, processed: bool = True, load_clustering: bool = False, remove_nontargeting: bool = False
    ) -> 'AnnData':
        if processed:
            adata = self.src.AdataH5.load()
            leiden_cols = [c for c in adata.obs.columns if c.startswith('leiden')]
            adata.obs = adata.obs.drop(columns=leiden_cols)
        else:
            from .modules.single_cell.process import init_adata

            adata = init_adata(
                manifest=self.src.Manifest.load(),
                cell_metadata=self.src.CellMetadata.load(),
                cell_x_gene=self.src.CellxGene.load(),
                cell_x_protein=self.src.CellxProt.load() if self.src.pr_detected else None,
            )

        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()

        if load_clustering and self.src.ClusteringUmap.is_valid:
            df = self.src.ClusteringUmap.load().cast({'cell_id': pl.Utf8}).to_pandas().set_index('cell_id')
            adata.obs = adata.obs.merge(df, how='left', left_index=True, right_index=True)

            if not processed:
                adata.obsm['X_umap'] = adata.obs[['UMAP1', 'UMAP2']].to_numpy(dtype=np.float32)

                leiden_cols = adata.obs.columns[adata.obs.columns.str.startswith('leiden_')]
                adata.obs[leiden_cols] = adata.obs[leiden_cols].fillna(c.UNASSIGNED_CELL)

        return adata

    def _return_image(
        self,
        img_path: str,
        dask: bool = False,
        shape: tuple[int] | None = None,
        use_cache: bool | None = None,
        dtype: np.dtype = np.uint16,
    ) -> np.ndarray:
        use_cache = self.use_cache if use_cache is None else use_cache
        if dask:
            return io.import_image_dask(img_path=img_path, shape=shape or self.shape, dtype=dtype, use_cache=use_cache)
        else:
            return io.import_image(img_path=img_path, use_cache=use_cache)

    def load_nuclear_image(self, dask: bool = False, **kwargs) -> np.ndarray:
        img_path = self.src.HnEDir.get_img(c.NUCLEAR_STAIN)
        return self._return_image(img_path=img_path, dask=dask, **kwargs)

    def load_cytoplasmic_image(self, dask: bool = False, **kwargs) -> np.ndarray:
        img_path = self.src.HnEDir.get_img(c.CYTOPLASMIC_STAIN)
        return self._return_image(img_path=img_path, dask=dask, **kwargs)

    def load_he_image(self, dask: bool = False, **kwargs) -> np.ndarray:
        img_path = self.src.HnEDir.get_img(c.H_AND_E)
        return self._return_image(img_path=img_path, shape=self.shape + (3,), dask=dask, dtype=np.uint8, **kwargs)

    def load_protein_image(self, protein: str, dask: bool = False, **kwargs) -> np.ndarray:
        if not self.src.pr_detected:
            print('No protein data available')
            return None

        if protein not in self.src.ProteinDir.mapped_files:
            print(f'Protein image for {protein} not found.')
            return None

        img_path = self.src.ProteinDir.get_img(protein)

        return self._return_image(img_path=img_path, dask=dask, **kwargs)

    def load_segmentation(self, expanded: bool = True, key: str = False) -> np.ndarray:
        key = 'nuclei_exp' if expanded else 'nuclei'
        return io.import_segmentation(
            seg_path=self.src.Segmentation.p, expected_shape=self.shape, labels_key=key, use_cache=self.use_cache
        )

    def load_bead_mask(self) -> np.ndarray:
        if self.src.BeadMask.is_valid:
            return np.load(self.src.BeadMask.p)['bead_mask']
        return None

    def list_content(self, subdir=None) -> dict:
        if subdir is None:
            subdir = ''

        list_path = self.smp_dir / subdir
        output = os.listdir(list_path)

        contents = {'dirs': [], 'files': []}
        for item in output:
            if os.path.isdir(list_path / item):
                contents['dirs'].append(item)
            if os.path.isfile(list_path / item):
                contents['files'].append(item)

        return contents

    def reroute_source(self, validator: 'BaseValidator', out_dir: str, overwrite: bool = False) -> None:

        out_obj = getattr(self.src, validator.__name__)
        out_obj.root = Path(out_dir)
        io.pathval.ensure_parent_dir(out_obj.p)

        if out_obj.path_exists() and not overwrite:
            raise RuntimeError(
                f'Operation aborted! {validator.__name__} already exists at:\n{ut.PGAP}{out_obj.p}\nUse overwrite=True to ignore this.',
            )

        suffix = 'overriding existing file' if out_obj.path_exists() else 'creating new file'
        ut.log_with_path(f'Using the following path for {validator.__name__} output ({suffix}):', out_obj.p)
