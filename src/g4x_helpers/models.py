import glob
import json
import logging
import os
import re
from dataclasses import InitVar, asdict, dataclass, field
from pathlib import Path
from typing import Any, Generator, Iterable, Iterator, List, Literal, Optional, Tuple, Union

import anndata
import glymur
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import tifffile
from packaging import version

import g4x_helpers.utils as utils
from g4x_helpers.utils import setup_logger

glymur.set_option('lib.num_threads', 16)


@dataclass()
class G4XOutput:
    run_base: Union[Path, str]
    sample_id: str = None
    log_level: InitVar[int] = logging.INFO

    def __post_init__(self, log_level):
        self.run_base = Path(self.run_base)

        if self.sample_id is None:
            self.sample_id = self.run_base.name

        ## this only sets up a stream handler, no file handlers
        _ = self.setup_logger(stream_level=log_level, clear_handlers=True)

        with open(self.run_base / 'run_meta.json', 'r') as f:
            self.run_meta = json.load(f)

        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / 'transcript_panel.csv', index_col=0, header=0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())

        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / 'protein_panel.csv', index_col=0, header=0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())

        self._get_shape()
        self._get_cell_coords()
        self._get_tissue_crop()

        roi_data = utils.Roi(extent=self.extent, name='data')
        roi_tissue = utils.Roi(extent=self.tissue_extent, name='tissue', sub_rois=2)
        self.rois = {
            'data': roi_data,
            'tissue': roi_tissue,
            'A': roi_tissue.sub_rois[0],
            'B': roi_tissue.sub_rois[1],
            'C': roi_tissue.sub_rois[2],
            'D': roi_tissue.sub_rois[3],
        }

    def __repr__(self):
        return f"""
        G4XOutput for {self.sample_id} located at {self.run_base}.
            contains transcript: {self.includes_transcript}
            contains protein: {self.includes_protein}
        """
    # region properties
    @property
    def machine(self) -> str:
        return self.run_meta.get('machine', None)

    @property
    def run_id(self) -> str:
        return self.run_meta.get('run_id', None)

    @property
    def fc(self) -> str:
        return self.run_meta.get('fc', None)

    @property
    def lane(self) -> str:
        return self.run_meta.get('lane', None)

    @property
    def platform(self) -> str:
        return self.run_meta.get('platform', None)

    @property
    def transcript_panel(self) -> dict:
        return self.run_meta.get('transcript_panel', None)

    @property
    def protein_panel(self) -> dict:
        return self.run_meta.get('protein_panel', None)

    @property
    def software(self) -> str:
        return self.run_meta.get('software', None)

    @property
    def software_version(self) -> str:
        return self.run_meta.get('software_version', None)

    @property
    def logger(self) -> logging.Logger:
        return logging.getLogger(f'{self.sample_id}_G4XOutput')

    @property
    def includes_transcript(self) -> bool:
        return self.transcript_panel is not None

    @property
    def includes_protein(self) -> bool:
        return self.protein_panel is not None

    # region methods
    def setup_logger(
        self,
        *,
        stream_logger: Optional[bool] = True,
        stream_level: Optional[int] = logging.INFO,
        file_logger: Optional[bool] = False,
        file_level: Optional[int] = logging.INFO,
        clear_handlers: Optional[bool] = True,
    ) -> None:
        _ = setup_logger(
            f'{self.sample_id}_G4XOutput',
            stream_logger=stream_logger,
            stream_level=stream_level,
            file_logger=file_logger,
            file_level=file_level,
            clear_handlers=clear_handlers,
        )

    def thumb(self, size=3, ax=None, crop=True, scale=1, show_tissue_bounds=True, background=True):
        if ax is None:
            fig, ax = plt.subplots(figsize=(size, size), dpi=300, layout='compressed')

        ax.set_title(self.sample_id)

        if not hasattr(self, 'he_thumb'):
            self.he_thumb = self._get_image(parent_directory='h_and_e', pattern='h_and_e_thumbnail.jpg')
        plot_image = self.he_thumb.copy()

        roi = self.rois['data'] if crop == 'data' else self.rois[crop]

        # if background:
        #     if isinstance(background, bool):
        #         background = '#F7F3F9'  # HE Image background color
        #     ax.set_facecolor(background)

        if scale != 1:
            roi = roi.scale(scale)

        display_size = np.max([roi.width / 5, roi.height / 5])
        plot_image = utils.downsample_mask(plot_image, ax=ax, mask_length=display_size, downsample='auto')

        ax.imshow(plot_image, extent=self.rois['data'].extent, origin='lower')

        if show_tissue_bounds:
            if crop == 'data':
                crop_roi = self.rois['tissue']
                utils.add_rois(ax, [crop_roi], color='orange', label=False)

        ax.set_xlim(roi.xlims)
        ax.set_ylim(roi.ylims)
        ax.set_xticks([])
        ax.set_yticks([])

    def load_adata(self, *, remove_nontargeting: Optional[bool] = True, load_clustering: Optional[bool] = True) -> anndata.AnnData:
        adata = sc.read_h5ad(self.run_base / 'single_cell_data' / 'feature_matrix.h5')

        # TODO this may conflict with the last line?!
        # adata.obs_names = adata.obs['cell_id']
        adata.var_names = adata.var['gene_id']

        adata.obs['sample_id'] = adata.uns['sample_id'] = self.sample_id
        adata.uns['software_version'] = self.software_version

        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()

        if load_clustering:
            try:
                df = pd.read_csv(self.run_base / 'single_cell_data' / 'clustering_umap.csv.gz', index_col=0, header=0)
                adata.obs = adata.obs.merge(df, how='left', left_index=True, right_index=True)
                umap_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
                adata.obsm['X_umap'] = adata.obs[[f'{umap_key}_1', f'{umap_key}_2']].to_numpy(dtype=float)
            except Exception as e:
                self.logger.exception('No single-cell clustering results found.')
                self.logger.exception(e)

        adata.obs_names = f'{self.sample_id}-' + adata.obs['cell_id'].str.split('-').str[1]
        return adata

    def load_protein_image(self, protein: str) -> np.ndarray:
        if not self.protein_panel:
            self.logger.error('No protein results.')
            return None
        assert protein in self.protein_panel, print(f'{protein} not in protein panel.')
        return self._get_image('protein', f'{protein}.*')

    def load_he_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'h_and_e.*')

    def load_nuclear_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'nuclear.*')

    def load_eosin_image(self) -> np.ndarray:
        return self._get_image('h_and_e', 'eosin.*')

    def load_segmentation(self, expanded: Optional[bool] = True) -> np.ndarray:
        arr = np.load(self.run_base / 'segmentation' / 'segmentation_mask.npz')
        if expanded:
            return arr['nuclei_exp']
        else:
            return arr['nuclei']

    def load_transcript_table(self, return_polars: Optional[bool] = False) -> Union[pd.DataFrame, pl.DataFrame]:
        reads = pl.read_csv(self.run_base / 'rna' / 'transcript_table.csv.gz')
        if not return_polars:
            reads = reads.to_pandas()
        return reads

    def list_contents(self, subdir=None):
        if subdir is None:
            subdir = ''

        list_path = self.run_base / subdir
        output = os.listdir(list_path)

        contents = {'dirs': [], 'files': []}
        for item in output:
            if os.path.isdir(list_path / item):
                contents['dirs'].append(item)
            if os.path.isfile(list_path / item):
                contents['files'].append(item)

        return contents
    
    # region internal
    def _get_image(self, parent_directory: str, pattern: str) -> np.ndarray:
        img_path = next((self.run_base / parent_directory).glob(pattern), None)
        if img_path is None:
            self.logger.error('Image file does not exist.')
            return None
        if img_path.suffix in ('.jp2', '.jpg'):
            img = glymur.Jp2k(img_path)[:]
        else:
            img = tifffile.imread(img_path)
        return img

    def _get_coord_order(self, verbose: bool = False) -> str:
        critical_version = version.parse('2.11.1')

        detected_caretta_version = version.parse(self.software_version)
        if verbose:
            print(f'Detected Caretta version: {detected_caretta_version}')

        if detected_caretta_version >= critical_version:
            return 'yx'
        else:
            return 'xy'

    def _get_shape(self):
        image_path = self.run_base / 'h_and_e' / 'nuclear.jp2'

        image = glymur.Jp2k(image_path)
        self.shape = image.shape
        self.extent = (0, self.shape[1], 0, self.shape[0])

    def _get_cell_coords(self):
        cell_meta_path = self.run_base / 'single_cell_data' / 'cell_metadata.csv.gz'

        if self._get_coord_order() == 'yx':
            coord_select = ['cell_y', 'cell_x']
        else:
            coord_select = ['cell_x', 'cell_y']

        self.cell_coords = pl.read_csv(cell_meta_path).select(coord_select).to_numpy()

    def _get_tissue_crop(self, ubins=100, threshold=0.1, expand=0.1, order='yx', **kwargs):
        def lims_to_extent(lims, lims_order='yx'):
            if lims_order == 'yx':
                y0, y1 = lims[0]
                x0, x1 = lims[1]
            elif lims_order == 'xy':
                x0, x1 = lims[0]
                y0, y1 = lims[1]

            return (x0, x1, y0, y1)

        lims = utils.find_tissue(self.cell_coords, ubins=ubins, threshold=threshold, expand=expand, order=order, smp_shape=self.shape, **kwargs)
        self.tissue_lims = lims['lims']
        self.tissue_center = lims['centroid']
        self.tissue_extent = lims_to_extent(lims['lims'])
