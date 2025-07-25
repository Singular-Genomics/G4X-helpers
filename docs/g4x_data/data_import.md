<br>

# G4X data import

The multi-modal output of the G4X spatial sequencer comprises images, tables and annotated data matrices, which allow deep exploration of your sample. The `single_cell_data` folder in the G4X output contains the final processed form of the data, after transcript and image signals have been aggregated for each each segmented cell. There are several excellent open-source tools available that enable the full stack of analytical needs to gain biological insight from this data.  

#### Below we illustrate data import strategies for [Python](#if-you-are-working-in-python) and [R](#if-you-are-working-in-r) users:

<br>

## if you are working in Python
---
Some excellent tools are:

+ `scanpy`: full-featured single-cell analysis suite
+ `rapids-singlecell`: expands scanpy with GPU-support
+ `squidpy`: functionality for the analysis of spatial data
+ `spatial-data`: analysis and management of spatial data (requires data conversion)

<br>

All of these tools accept, or incorporate [`anndata`](https://anndata.readthedocs.io/en/stable/) objects as the basic representation of your single-cell data. There are several ways to load your data as an `anndata` object into your Python session. Below we will illustrate a few methods that will produce equivalent outputs.


### 1. G4X-helpers

If you have G4X-helpers installed, you can use the `G4Xoutput()` class to access the anndata of your sample via the `load_adata()` method.


```py
import g4x_helpers as g4x

run_base = '/path/to/g4x_output/sample_x1'

sample = g4x.G4Xoutput(run_base=run_base)
adata = sample.load_adata(remove_nontargeting=False, load_clustering=False) 
```

!!!note "Note: `load_adata()` options"
    
    Two options are provided that will impact what will be loaded. In the above example we are overriding the default so that the output matches the other methods.  
    `remove_nontargeting`: bool (default=True)  
    `load_clustering`: bool (default=True)


### 2. scanpy

You can achieve the same result by pointing scanpy's `read_h5ad()` function the `feature_matrix.h5` in your G4X output directory

```py
from pathlib import Path
import scanpy as sc

run_base = Path('/path/to/g4x_output/sample_x1')
ad_file = run_base / 'single_cell_data' / 'feature_matrix.h5'

adata = sc.read_h5ad(ad_file)
```

!!!note
    If you have G4X-helpers installed, then `scanpy` will be available as one of its dependencies. If not, please refer to the [scanpy](https://scanpy.readthedocs.io/en/stable/) documentation for installation guides.


### 3. building from raw data

If you do not wish to use the pre-generated `feature_matrix.h5` you can replicate a similar object by reading the counts and cell metadata into an `anndata` object.

```py
from pathlib import Path
import anndata as ad
import pandas as pd
from scipy import sparse

run_base = Path('/path/to/g4x_output/sample_x1')

txcounts_path = run_base / 'single_cell_data' / 'cell_by_transcript.csv.gz'
metadata_path = run_base / 'single_cell_data' / 'cell_metadata.csv.gz'

counts = pd.read_csv(txcounts_path, index_col='label')
X = sparse.csr_matrix(counts.values)

metadata = pd.read_csv(metadata_path, index_col='label')

adata = ad.AnnData(X=X, obs=metadata)
adata.var_names = counts.columns
```

!!!note
    If you have G4X-helpers installed, then `anndata`, `pandas` and `scipy` will be available as dependencies.

<br>

## if you are working in R
---

The most feature rich package for single-cell analysis in R is:

+ `Seurat`: full-featured single-cell analysis suite with spatial analysis capabilities

To work with your data in Seurat, it needs to be loaded into a `SeuratObject`, which is an annotated data structure.

```R
library('Seurat')
 
run_base = c('/path/to/g4x_output/sample_x1')

txcounts_path = file.path(run_base, "single_cell_data/cell_by_transcript.csv.gz")
metadata_path = file.path(run_base, "single_cell_data/cell_metadata.csv.gz")

counts <- read.csv(txcounts_path, row.names = 1) 
counts <- t(counts) # transpose to match Seurat input
 
metadata <- read.csv(metadata_path, row.names = 1) 
 
sobj <- CreateSeuratObject(counts = counts, meta.data = metadata)
```

--8<-- "_partials/end_cap.md"
