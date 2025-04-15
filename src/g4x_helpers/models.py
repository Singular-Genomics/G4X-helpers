import os
import re
import json
import glob
import glymur
import logging
import anndata
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from pathlib import Path
from dataclasses import dataclass, field, asdict, InitVar
from typing import List, Tuple, Union, Iterator, Iterable, Any, Generator, Optional, Literal

from g4x_helpers.utils import setup_logger

glymur.set_option("lib.num_threads", 16)

@dataclass()
class G4XOutput():

    run_base: Union[Path, str]
    sample_id: str
    log_level: InitVar[int] = logging.INFO

    def __post_init__(self, log_level):
        self.run_base = Path(self.run_base)

        ## this only sets up a stream handler, no file handlers
        _= self.setup_logger(stream_level= log_level, clear_handlers= True)

        with open(self.run_base / "run_meta.json", "r") as f:
            self.run_meta = json.load(f)
        
        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / "transcript_panel.csv", index_col= 0, header= 0)
            self.transcript_panel_dict = transcript_panel.to_dict()['panel_type']
            self.genes = list(self.transcript_panel_dict.keys())
        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / "protein_panel.csv", index_col= 0, header= 0)
            self.protein_panel_dict = protein_panel.to_dict()['panel_type']
            self.proteins = list(self.protein_panel_dict.keys())


    def __repr__(self):
        return f"""
        G4XOutput for {self.sample_id} located at {self.run_base}.
            contains transcript: {self.includes_transcript}
            contains protein: {self.includes_protein}
        """


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
      return logging.getLogger(f"{self.sample_id}_G4XOutput")
    @property
    def includes_transcript(self) -> bool:
        return self.transcript_panel is not None
    @property
    def includes_protein(self) -> bool:
        return self.protein_panel is not None
    
    
    def setup_logger(
        self,
        *,
        stream_logger: Optional[bool] = True,
        stream_level: Optional[int] = logging.INFO,
        file_logger: Optional[bool] = False,
        file_level: Optional[int] = logging.INFO,
        clear_handlers: Optional[bool] = True
    ) -> None:
        _= setup_logger(
            f"{self.sample_id}_G4XOutput",
            stream_logger= stream_logger,
            stream_level= stream_level,
            file_logger= file_logger,
            file_level= file_level,
            clear_handlers= clear_handlers
        )

    def load_adata(
        self,
        *,
        remove_nontargeting: Optional[bool] = True,
        load_clustering: Optional[bool] = True
    ) -> anndata.AnnData:
        adata = sc.read_h5ad(self.run_base / "single_cell_data" / "feature_matrix.h5")
        adata.obs_names = adata.obs['cell_id']
        adata.var_names = adata.var['gene_id']
        if remove_nontargeting:
            adata = adata[:, adata.var.query(" probe_type == 'targeting' ").index].copy()
        if load_clustering:
            try:
                df = pd.read_csv(self.run_base / "single_cell_data" / "clustering_umap.csv.gz", index_col= 0, header= 0)
                adata.obs = adata.obs.merge(df, how= 'left', left_index= True, right_index= True)
                umap_key = '_'.join(sorted([x for x in adata.obs.columns if 'X_umap' in x])[0].split('_')[:-1])
                adata.obsm['X_umap'] = adata.obs[[f"{umap_key}_1", f"{umap_key}_2"]].to_numpy(dtype= float)
            except Exception as e:
                self.logger.exception(f"No single-cell clustering results found.")
                self.logger.exception(e)
        adata.obs_names = f"{self.sample_id}-" + adata.obs['cell_id'].str.split('-').str[0]
        return adata

    def _get_image(self, parent_directory: str, pattern: str) -> np.ndarray:
        img_path = next((self.run_base / parent_directory).glob(pattern), None)
        if img_path is None:
            self.logger.error("Image file does not exist.")
            return None
        if img_path.suffix == ".jp2":
            img = glymur.Jp2k(img_path)[:]
        else:
            img = tifffile.imread(img_path)
        return img

    def load_protein_image(self, protein: str) -> np.ndarray:
        if not self.protein_panel:
            self.logger.error("No protein results.")
            return None
        assert protein in self.protein_panel, print(f"{protein} not in protein panel.")
        return _get_image("protein", f"{protein}.*")
    def load_he_image(self) -> np.ndarray:
        return _get_image("h_and_e", "h_and_e.*")
    def load_nuclear_image(self) -> np.ndarray:
        return _get_image("h_and_e", "nuclear.*")
    def load_eosin_image(self) -> np.ndarray:
        return _get_image("h_and_e", "eosin.*")
    
    def load_segmentation(self, expanded: Optional[bool] = True) -> np.ndarray:
        arr = np.load(self.run_base / "segmentation" / "segmentation_mask.npz")
        if expanded:
            return arr["nuclei_exp"]
        else:
            return arr["nuclei"]
        
    def load_transcript_table(self, return_polars: Optional[bool] = False) -> Union[pd.DataFrame, pl.DataFrame]:
        df = pl.read_csv(self.run_base / "rna" / "transcript_table.csv.gz")
        if not return_polars:
            df = df.to_pandas()
        return d