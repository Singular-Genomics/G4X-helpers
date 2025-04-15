import os
import re
import json
import glob
import glymur
import anndata
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
from pathlib import Path
from dataclasses import dataclass, field, asdict, InitVar
from typing import List, Tuple, Union, Iterator, Iterable, Any, Generator, Optional, Literal

glymur.set_option("lib.num_threads", 16)

@dataclass(kw_only= True)
class G4XOutput():

    run_base: Union[Path, str]

    def __post_init__(self):
        self.run_base = Path(self.run_base)

        with open(self.run_base / "run_meta.json", "r") as f:
            run_meta = json.load(f)
            for key, value in run_meta.items():
                setattr(self, key, value)
        
        if self.transcript_panel:
            transcript_panel = pd.read_csv(self.run_base / "transcript_panel.csv", index_col= 0, header= 0)
            self.transcript_panel = transcript_panel.to_dict()['panel_type']
            self.target_genes = list(self.transcript_panel.keys())
        if self.protein_panel:
            protein_panel = pd.read_csv(self.run_base / "protein_panel.csv", index_col= 0, header= 0)
            self.protein_panel = protein_panel.to_dict()['panel_type']
            self.target_proteins = list(self.protein_panel.keys())

    def load_adata(
        self,
        *,
        remove_nontargeting: Optional[bool] = True,
        load_clustering: Optional[bool] = True
    ) -> anndata.AnnData
        adata = sc.read_h5ad(self.run_base / "single_cell_data" / "feature_matrix.h5")
        adata.obs_names = adata.obs["cell_id"]
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
                print(f"No clustering results found. {e}")
        return adata

    def _get_image(self, parent_directory: str, pattern: str) -> np.ndarray:
        img_path = next((self.run_base / parent_directory).glob(pattern), None)
        if img_path is None:
            print("Image file does not exist.")
            return None
        if img_path.suffix == ".jp2":
            img = glymur.Jp2k(img_path)[:]
        else:
            img = tifffile.imread(img_path)
        return img

    def load_protein_image(self, protein: str) -> np.ndarray
        if not self.protein_panel:
            print("No protein results.")
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
        return df