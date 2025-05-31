import json
import logging
import os
from dataclasses import InitVar, dataclass  # , asdict, field
from pathlib import Path
from typing import List, Optional, Union

import anndata as ad
import glymur
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from geopandas.geodataframe import GeoDataFrame
from packaging import version

import g4x_helpers.g4x_viewer.bin_generator as bin_gen
import g4x_helpers.segmentation as reseg
import g4x_helpers.utils as utils

glymur.set_option('lib.num_threads', 24)


@dataclass()
class G4Xoutput:
    run_base: Path | str
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

        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / 'transcript_panel.csv', index_col=0, header=0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())

        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / 'protein_panel.csv', index_col=0, header=0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())

    # region properties
    @property
    def logger(self) -> logging.Logger:
        return logging.getLogger(f'{self.sample_id}_G4XOutput')

    # region methods
    def setup_logger(
        self,
        *,
        stream_logger: bool = True,
        stream_level: int = logging.INFO,
        file_logger: bool = False,
        file_level: int = logging.INFO,
        clear_handlers: bool = True,
    ) -> None:
        _ = utils.setup_logger(
            f'{self.sample_id}_G4XOutput',
            stream_logger=stream_logger,
            stream_level=stream_level,
            file_logger=file_logger,
            file_level=file_level,
            clear_handlers=clear_handlers,
        )

    def load_adata(self, *, remove_nontargeting: bool = True, load_clustering: bool = True) -> ad.AnnData:
        adata = sc.read_h5ad(self.run_base / 'single_cell_data' / 'feature_matrix.h5')

        adata.obs_names = adata.obs['cell_id']
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
        assert protein in self.proteins, self.logger.error(f'{protein} not in protein panel.')
        return self.load_image(signal=protein)

    def load_he_image(self) -> np.ndarray:
        return self.load_image(signal='h_and_e')

    def load_nuclear_image(self) -> np.ndarray:
        return self.load_image(signal='nuclear')

    def load_eosin_image(self) -> np.ndarray:
        return self.load_image(signal='eosin')

    def load_segmentation(self, expanded: bool = True) -> np.ndarray:
        arr = np.load(self.run_base / 'segmentation' / 'segmentation_mask.npz')
        if expanded:
            return arr['nuclei_exp']
        else:
            return arr['nuclei']

    def load_transcript_table(self, return_polars: bool = False) -> pd.DataFrame | pl.DataFrame:
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

    def intersect_segmentation(
        self,
        labels: GeoDataFrame | np.ndarray,
        *,
        out_dir: str | Path | None = None,
        exclude_channels: Optional[List[str]] = None,
        gen_bin_file: bool = True,
        n_threads: int = 4,
    ) -> None:
        mask = labels

        signal_list = ['nuclear', 'eosin'] + self.proteins

        if exclude_channels is not None:
            self.logger.info(f'Not processing channels: {", ".join(exclude_channels)}')
            if isinstance(exclude_channels, str):
                exclude_channels = [exclude_channels]

            signal_list = [item for item in signal_list if item not in exclude_channels]
        else:
            self.logger.info('Processing all channels.')
            

        if out_dir is None:
            self.logger.warning('out_dir was not specified, so files will be updated in-place.')
            out_dir = self.run_base
        else:
            self.logger.info(f'Using provided output directory: {out_dir}')
            out_dir = Path(out_dir)

        if isinstance(mask, GeoDataFrame):
            self.logger.info('Rasterizing provided GeoDataFrame')
            mask = reseg.rasterize_polygons(gdf=mask, target_shape=self.shape)

        self.logger.info('Extracting mask properties')
        segmentation_props = reseg.get_mask_properties(self, mask)

        self.logger.info('Assigning transcripts to mask labels')
        reads_new_labels = reseg.assign_tx_to_mask_labels(self, mask)

        self.logger.info('Extracting image signals')
        cell_by_protein = reseg.extract_image_signals(self, mask, signal_list=signal_list)

        self.logger.info('Building output data structures')
        cell_metadata = reseg._make_cell_metadata(segmentation_props, cell_by_protein)
        cell_by_gene = reseg._make_cell_by_gene(segmentation_props, reads_new_labels)
        adata = reseg._make_adata(cell_by_gene, cell_metadata)

        if out_dir:
            self.logger.info(f'Saving output files to {out_dir}')
        else:
            self.logger.info(f'No output directory specified, saving to ["custom"] directories in {self.run_base}.')

        outfile = reseg._create_custom_out(self, out_dir, 'segmentation', 'segmentation_mask_updated.npz')
        self.logger.debug(f'segmentation mask --> {outfile}')
        np.savez(outfile, cell_labels=mask)

        outfile = reseg._create_custom_out(self, out_dir, 'rna', 'transcript_table.csv.gz')
        self.logger.debug(f'transcript table --> {outfile}')
        reads_new_labels.write_csv(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_transcript.csv.gz')
        self.logger.debug(f'cell x transcript --> {outfile}')
        cell_by_gene.write_csv(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_by_protein.csv.gz')
        self.logger.debug(f'cell x protein --> {outfile}')
        cell_by_protein.write_csv(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'feature_matrix.h5')
        self.logger.debug(f'single-cell h5 --> {outfile}')
        adata.write_h5ad(outfile)

        outfile = reseg._create_custom_out(self, out_dir, 'single_cell_data', 'cell_metadata.csv.gz')
        self.logger.debug(f'cell metadata --> {outfile}')
        adata.obs.to_csv(outfile)

        protein_only_list = [p for p in signal_list if p not in ['nuclear', 'eosin']]
        if gen_bin_file:
            self.logger.info('Making G4X-Viewer bin file.')
            outfile = reseg._create_custom_out(self, out_dir, 'g4x_viewer', f'{self.sample_id}.bin')
            _ = bin_gen.seg_converter(
                adata=adata,
                seg_mask=mask,
                outpath=outfile,
                protein_list=[f'{x}_intensity_mean' for x in protein_only_list],
                n_threads=n_threads,
            )
            self.logger.debug(f'G4X-Viewer bin --> {outfile}')

    def load_image(self, signal: str, thumbnail: bool = False) -> tuple[np.ndarray, float, float]:
        img = utils.load_image_cached(self.run_base, signal, thumbnail=thumbnail)
        return img

    # region internal
    def _get_shape(self):
        img_path = self.run_base / 'h_and_e/nuclear.jp2'
        img = glymur.Jp2k(img_path)
        return img.shape

    def _get_coord_order(self, verbose: bool = False) -> str:
        critical_version = version.parse('2.11.1')

        detected_caretta_version = version.parse(self.software_version)
        if verbose:
            self.logger.debug(f'Detected Caretta version: {detected_caretta_version}')

        if detected_caretta_version >= critical_version:
            return 'yx'
        else:
            return 'xy'

    def _clear_image_cache(self):
        """Evict all cached images so subsequent calls re-read from disk."""
        utils.load_image_cached.cache_clear()

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
            setattr(self, key, self.run_meta.get(key, None))

        self.includes_transcript = False if self.transcript_panel == [] else True
        self.includes_protein = False if self.protein_panel == [] else True
        self.shape = self._get_shape()

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
