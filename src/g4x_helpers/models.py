import json
import logging
import os
from collections.abc import Iterator
from dataclasses import InitVar, dataclass  # , asdict, field
from functools import lru_cache
from pathlib import Path
from typing import Literal

import anndata as ad
import glymur
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import tifffile
from geopandas.geodataframe import GeoDataFrame

import g4x_helpers.g4x_viewer.bin_generator as bin_gen
import g4x_helpers.segmentation as reseg
import g4x_helpers.utils as utils

glymur.set_option('lib.num_threads', 8)


@dataclass()
class G4Xoutput:
    """
    Container for managing and processing data from a G4X run.

    This class initializes and loads metadata, image dimensions, transcript and protein panels,
    and sets up logging for downstream analysis of G4X output data. It provides methods to load
    images, segmentations, transcript data, and interact with single-cell and spatial analysis pipelines.

    Parameters
    ----------
    run_base : str or Path
        The base directory path of the G4X run. This should contain all run-related files including metadata,
        segmentation masks, panels, and feature matrices.
    sample_id : str, optional
        The sample ID to associate with the run. If not provided, it will be inferred from the `run_base` path.
    out_dir : str or Path, optional
        The output directory to use for writing logs and intermediate files. Defaults to `run_base` if not provided.
    log_level : int, optional
        The logging level for the stream logger (default is `logging.INFO`).

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
      - Sets up logging with optional stream and file loggers.
      - Loads metadata from `run_meta.json`.
      - Loads the shape of the segmentation mask.
      - Parses transcript and protein panel files (if present).
    """

    run_base: Path | str
    sample_id: str | None = None
    out_dir: InitVar[Path | str | None] = None
    log_level: InitVar[int] = logging.INFO

    def __post_init__(self, out_dir: Path | str | None, log_level: int):
        self.run_base = Path(self.run_base)

        if self.sample_id is None:
            self.sample_id = self.run_base.name

        _ = self.setup_logger(
            stream_logger=True,
            stream_level=log_level,
            file_logger=True,
            file_level=logging.DEBUG,
            out_dir=out_dir,
            clear_handlers=True,
        )

        _ = self._validate_run_base()

        with open(self.run_base / 'run_meta.json', 'r') as f:
            self.run_meta = json.load(f)

        self.shape = utils.npzGetShape(self.run_base / 'segmentation' / 'segmentation_mask.npz', 'nuclei')

        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / 'transcript_panel.csv', index_col=0, header=0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())

        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / 'protein_panel.csv', index_col=0, header=0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())

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

    # region properties
    @property
    def fc(self) -> str:
        return self.run_meta.get('fc', None)

    @property
    def feature_table_path(self) -> Path:
        return self.run_base / 'diagnostics' / 'transcript_table.parquet'

    @property
    def includes_protein(self) -> bool:
        return self.protein_panel != []

    @property
    def includes_transcript(self) -> bool:
        return self.transcript_panel != []

    @property
    def lane(self) -> str:
        return self.run_meta.get('lane', None)

    @property
    def logger(self) -> logging.Logger:
        if not hasattr(self, '_logger') or self._logger is None:
            self._logger = logging.getLogger(self.logger_name)
        return self._logger

    @property
    def logger_name(self) -> str:
        return f'{self.sample_id}_G4XOutput'

    @property
    def machine(self) -> str:
        return self.run_meta.get('machine', None)

    @property
    def platform(self) -> str:
        return self.run_meta.get('platform', None)

    @property
    def protein_panel(self) -> dict:
        return self.run_meta.get('protein_panel', None)

    @property
    def run_id(self) -> str:
        return self.run_meta.get('run_id', None)

    @property
    def software(self) -> str:
        return self.run_meta.get('software', None)

    @property
    def software_version(self) -> str:
        return self.run_meta.get('software_version', None)

    @property
    def transcript_panel(self) -> dict:
        return self.run_meta.get('transcript_panel', None)

    @property
    def transcript_table_path(self) -> Path:
        return self.run_base / 'rna' / 'transcript_table.csv.gz'

    # region methods
    def setup_logger(
        self,
        *,
        stream_logger: bool = True,
        stream_level: int = logging.INFO,
        file_logger: bool = False,
        file_level: int = logging.INFO,
        out_dir: Path | str | None = None,
        clear_handlers: bool = True,
    ) -> None:
        if out_dir is None:
            out_dir = os.getcwd()
        _ = utils.setup_logger(
            f'{self.sample_id}_G4XOutput',
            stream_logger=stream_logger,
            stream_level=stream_level,
            file_logger=file_logger,
            file_level=file_level,
            out_dir=out_dir,
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
            df = pd.read_csv(self.run_base / 'single_cell_data' / 'clustering_umap.csv.gz', index_col=0, header=0)
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
                self.logger.error('No protein results.')
                return None
            if protein is None or protein not in self.proteins:
                self.logger.error(f'{protein} not in protein panel.')
                return None
            pattern = f'{protein}_thumbnail.*' if thumbnail else f'{protein}.*'
            directory = 'protein'
        else:
            pattern_base = {'h_and_e': 'h_and_e', 'nuclear': 'nuclear', 'eosin': 'eosin'}.get(image_type)

            if not pattern_base:
                self.logger.error(f'Unknown image type: {image_type}')
                return None

            pattern = f'{pattern_base}_thumbnail.*' if thumbnail else f'{pattern_base}.*'
            directory = 'h_and_e'

        if cached:
            return self._load_image_cached(self.run_base, directory, pattern)
        else:
            return self._load_image(self.run_base, directory, pattern)

    def load_protein_image(self, protein: str, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('protein', thumbnail=thumbnail, protein=protein, cached=cached)

    def load_he_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('h_and_e', thumbnail=thumbnail, cached=cached)

    def load_nuclear_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('nuclear', thumbnail=thumbnail, cached=cached)

    def load_eosin_image(self, thumbnail: bool = False, cached: bool = False) -> np.ndarray:
        return self.load_image_by_type('eosin', thumbnail=thumbnail, cached=cached)

    def load_segmentation(self, expanded: bool = True) -> np.ndarray:
        arr = np.load(self.run_base / 'segmentation' / 'segmentation_mask.npz')
        if expanded:
            return arr['nuclei_exp']
        else:
            return arr['nuclei']

    def _load_table(
        self, file_path: str | Path, return_polars: bool = True, lazy: bool = False, columns: list[str] | None = None
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        file_path = Path(file_path)
        if lazy:
            if file_path.suffix == 'parquet':
                reads = pl.scan_parquet(file_path)
            else:
                reads = pl.scan_csv(file_path)
        else:
            if file_path.suffix == 'parquet':
                reads = pl.read_parquet(file_path)
            else:
                reads = pl.read_csv(file_path)
        if columns:
            reads = reads.select(columns)
        if not return_polars:
            reads = reads.collect().to_pandas()
        return reads

    def load_feature_table(
        self, *, return_polars: bool = True, lazy: bool = False, columns: list[str] | None = None
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        return self._load_table(self.feature_table_path, return_polars, lazy, columns)

    def load_transcript_table(
        self, *, return_polars: bool = True, lazy: bool = False, columns: list[str] | None = None
    ) -> pd.DataFrame | pl.DataFrame | pl.LazyFrame:
        return self._load_table(self.transcript_table_path, return_polars, lazy, columns)

    def stream_features(self, batch_size: int, columns: str | list[str] | None = None) -> Iterator[pl.DataFrame]:
        df = pl.scan_parquet(self.feature_table_path)
        if columns:
            df = df.select(columns)
        offset = 0
        while True:
            batch = df.slice(offset, batch_size).collect()
            if batch.is_empty():
                break
            yield batch
            offset += batch_size

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
        exclude_channels: list[str] | None = None,
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
            self.logger.info('Rasterizing provided GeoDataFrame.')
            mask = reseg.rasterize_polygons(gdf=mask, target_shape=self.shape)

        self.logger.info('Extracting mask properties.')
        segmentation_props = reseg.get_mask_properties(self, mask)

        self.logger.info('Assigning transcripts to mask labels.')
        reads_new_labels = reseg.assign_tx_to_mask_labels(self, mask)

        self.logger.info('Extracting protein signal.')
        cell_by_protein = reseg.extract_image_signals(self, mask, signal_list=signal_list)

        self.logger.info('Building output data structures.')
        cell_metadata = reseg._make_cell_metadata(segmentation_props, cell_by_protein)
        cell_by_gene = reseg._make_cell_by_gene(segmentation_props, reads_new_labels)
        adata = reseg._make_adata(cell_by_gene, cell_metadata)

        if out_dir:
            self.logger.info(f'Saving output files to {out_dir}')
        else:
            self.logger.info(f'No output directory specified, saving to ["custom"] directories in {self.run_base}.')

        outfile = self._create_custom_out(out_dir, 'segmentation', 'segmentation_mask_updated.npz')
        self.logger.debug(f'segmentation mask --> {outfile}')
        np.savez(outfile, cell_labels=mask)

        outfile = self._create_custom_out(out_dir, 'rna', 'transcript_table.csv')
        self.logger.debug(f'transcript table --> {outfile}.gz')
        reads_new_labels.write_csv(outfile)
        _ = utils.gzip_file(outfile, remove_original=True)

        outfile = self._create_custom_out(out_dir, 'single_cell_data', 'cell_by_transcript.csv')
        self.logger.debug(f'cell x transcript --> {outfile}.gz')
        cell_by_gene.write_csv(outfile)
        _ = utils.gzip_file(outfile, remove_original=True)

        outfile = self._create_custom_out(out_dir, 'single_cell_data', 'cell_by_protein.csv')
        self.logger.debug(f'cell x protein --> {outfile}.gz')
        cell_by_protein.write_csv(outfile)
        _ = utils.gzip_file(outfile, remove_original=True)

        outfile = self._create_custom_out(out_dir, 'single_cell_data', 'feature_matrix.h5')
        self.logger.debug(f'single-cell h5 --> {outfile}')
        adata.write_h5ad(outfile)

        outfile = self._create_custom_out(out_dir, 'single_cell_data', 'cell_metadata.csv.gz')
        self.logger.debug(f'cell metadata --> {outfile}')
        adata.obs.to_csv(outfile, compression='gzip')

        protein_only_list = [p for p in signal_list if p not in ['nuclear', 'eosin']]
        if gen_bin_file:
            self.logger.info('Making G4X-Viewer bin file.')
            outfile = self._create_custom_out(out_dir, 'g4x_viewer', f'{self.sample_id}.bin')
            _ = bin_gen.seg_converter(
                adata=adata,
                seg_mask=mask,
                out_path=outfile,
                protein_list=[f'{x}_intensity_mean' for x in protein_only_list],
                n_threads=n_threads,
                logger=self.logger,
            )
            self.logger.debug(f'G4X-Viewer bin --> {outfile}')

    # region internal
    @staticmethod
    def _load_image_base(run_base: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        img_path = next((run_base / parent_directory).glob(pattern), None)
        if img_path is None:
            raise FileNotFoundError(f'No file matching {pattern} found.')
        if img_path.suffix == '.jp2' or img_path.suffix == '.jpg':
            img = glymur.Jp2k(img_path)[:]
        elif img_path.suffix == '.png':
            img = plt.imread(img_path)
        else:
            img = tifffile.imread(img_path)
        return img

    @staticmethod
    @lru_cache(maxsize=None)
    def _load_image_cached(run_base: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        return G4Xoutput._load_image_base(run_base, parent_directory, pattern)

    @staticmethod
    def _load_image(run_base: str, parent_directory: str, pattern: str) -> tuple[np.ndarray, float, float]:
        return G4Xoutput._load_image_base(run_base, parent_directory, pattern)

    def _get_shape(self):
        img = self.load_he_image()
        return img.shape

    def _clear_image_cache(self):
        """Evict all cached images so subsequent calls re-read from disk."""
        self._load_image.cache_clear()

    def _validate_run_base(self):
        """check that all expected outputs are present."""
        self.logger.debug('Validating run_base.')

        required_paths = [self.run_base / 'run_meta.json', self.run_base / 'single_cell_data' / 'feature_matrix.h5']

        for p in required_paths:
            if not p.is_file():
                self.logger.error(f'{p} does not exist.')
                raise FileNotFoundError(f'{p} does not exist.')

    def _create_custom_out(
        self,
        out_dir: Path | str,
        parent_folder: Path | str | None = None,
        file_name: str | None = None,
    ) -> Path:
        custom_out = Path(out_dir) / parent_folder

        # Ensure the directory exists
        custom_out.mkdir(parents=True, exist_ok=True)

        # Prepare output file path
        outfile = custom_out / file_name

        utils.delete_existing(outfile)

        return outfile
