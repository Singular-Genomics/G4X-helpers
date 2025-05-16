import json
import logging
import os
from dataclasses import InitVar, dataclass  # , asdict, field
from pathlib import Path
from typing import Optional, Tuple, Union  # , Any, Generator, Iterable, Iterator, List, Literal

import anndata as ad
import glymur
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from geopandas.geodataframe import GeoDataFrame
from packaging import version
from functools import lru_cache
import matplotlib.pyplot as plt

import g4x_helpers.segmentation as reseg
import g4x_helpers.utils as utils

glymur.set_option('lib.num_threads', 16)


@dataclass()
class G4Xoutput:
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

        self._populate_attrs()
        self._get_shape()
        
        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / 'transcript_panel.csv', index_col=0, header=0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())

        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / 'protein_panel.csv', index_col=0, header=0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())

    

    @property
    def logger(self) -> logging.Logger:
        return logging.getLogger(f'{self.sample_id}_G4XOutput')

    def _populate_attrs(self):
        for key in (
            'machine',
            'run_id',
            'fc',
            'lane',
            'platform',
            'software',
            'software_version',
            'transcript_panel',
            'protein_panel',
        ):
            setattr(self, key, self.run_meta.get(key))

        self.includes_transcript = False if self.transcript_panel == [] else True
        self.includes_protein = False if self.protein_panel == [] else True

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
        _ = utils.setup_logger(
            f'{self.sample_id}_G4XOutput',
            stream_logger=stream_logger,
            stream_level=stream_level,
            file_logger=file_logger,
            file_level=file_level,
            clear_handlers=clear_handlers,
        )



    def load_adata(
        self, *, remove_nontargeting: Optional[bool] = True, load_clustering: Optional[bool] = True
    ) -> ad.AnnData:
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
        assert protein in self.proteins, print(f'{protein} not in protein panel.')
        return self.load_image(signal=protein, thumbnail=False, return_th=False)

    def load_he_image(self) -> np.ndarray:
        return self.load_image(signal='h_and_e', thumbnail=False, return_th=False)

    def load_nuclear_image(self) -> np.ndarray:
        return self.load_image(signal='nuclear', thumbnail=False, return_th=False)

    def load_eosin_image(self) -> np.ndarray:
        return self.load_image(signal='eosin', thumbnail=False, return_th=False)

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

    def list_content(self, subdir=None):
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


    def run_g4x_segmentation(
        self, labels: GeoDataFrame | np.ndarray, out_dir=None, include_channels=None, exclude_channels=None
    ):
        mask = labels

        if isinstance(mask, GeoDataFrame):
            print('Rasterizing provided GeoDataFrame')
            mask = reseg.rasterize_polygons(gdf=mask, target_shape=self.shape)

        print('Extracting mask properties')
        segmentation_props = reseg.get_mask_properties(self, mask)

        print('Assigning transcripts to mask labels')
        reads_new_labels = reseg.assign_tx_to_mask_labels(self, mask)

        print('Extracting image signals')
        cell_by_protein = reseg.extract_image_signals(
            self, mask, include_channels=include_channels, exclude_channels=exclude_channels
        )

        print('Building output data structures')
        cell_metadata = reseg._make_cell_metadata(segmentation_props, cell_by_protein)
        cell_by_gene = reseg._make_cell_by_gene(segmentation_props, reads_new_labels)
        adata = reseg._make_adata(cell_by_gene, cell_metadata)

        if out_dir:
            print(f'Saving output files to {out_dir}')
        else:
            print(f'No output directory specified, saving to ["custom"] directories in {self.run_base}.')

        outfile = reseg._create_custom_out(self, out_dir, 'segmentation', 'segmentation_mask.npz')
        np.savez(outfile, cell_labels=mask)

        outfile = reseg._create_custom_out(self, out_dir, 'rna', 'transcript_table.csv.gz')
        reads_new_labels.write_csv(outfile)
        
        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_transcript.csv.gz')
        cell_by_gene.write_csv(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_protein.csv.gz')
        cell_by_protein.write_csv(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'feature_matrix.h5')
        adata.write_h5ad(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_metadata.csv.gz')
        adata.obs.to_csv(outfile)


    def load_image(
        self, signal: str, thumbnail: bool = False) -> Tuple[np.ndarray, float, float]:
        img = _load_image_cached(self.run_base, signal, thumbnail=thumbnail)
        return img

    # region internal
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


    def _clear_image_cache(self):
        """Evict all cached images so subsequent calls re-read from disk."""
        _load_image_cached.cache_clear()

    # region dunder
    def __repr__(self):
        machine_num = self.machine.removeprefix('g4-').lstrip('0')
        mac_run_id = f'G{machine_num.zfill(2)}-{self.run_id}'

        repr_string = f'G4X_output @ {self.run_base}\n'
        repr_string += f'Sample: \033[1m{self.sample_id}\033[0m of {mac_run_id}, {self.fc}\n'

        shp = (np.array(self.shape) * 0.3125) / 1000
        repr_string += f'imaged area: ({shp[1]:.2f} x {shp[0]:.2f}) mm\n\n'

        if self.includes_transcript:
            repr_string += f'Transcript panel with {len(self.genes)} genes\t[{", ".join(self.genes[0:5])} ... ]\n'

        if self.includes_protein:
            repr_string += f'Protein panel with {len(self.proteins)} proteins\t[{", ".join(self.proteins[0:5])} ... ]\n'

        return repr_string
    
@lru_cache(maxsize=None)
def _load_image_cached(run_base: str, signal: str, thumbnail: bool) -> Tuple[np.ndarray, float, float]:
    """This is cached per (run_base, signal, thumbnail)."""
    base = Path(run_base)

    folder = 'h_and_e' if signal in ('h_and_e', 'nuclear', 'eosin') else 'protein'

    p = base / folder
    suffix = '.jpg' if (thumbnail and signal == 'h_and_e') else ('.png' if thumbnail else '.jp2')

    fname = f'{signal}_thumbnail{suffix}' if thumbnail else f'{signal}.jp2'
    img_file = p / fname

    if img_file.suffix == '.png':
        img = plt.imread(img_file)
    else:
        img = glymur.Jp2k(str(img_file))[:]

    return img