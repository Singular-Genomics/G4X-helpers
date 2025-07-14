# G4X data import


### G4X-helpers (python)
```py
import g4x_helpers as g4x

run_base = '/path/to/g4x_output'

sample = g4x.G4Xoutput(run_base=run_base)
```

### Scanpy (python)
```py
from pathlib import Path
import scanpy as sc

run_base = Path('/path/to/g4x_output')
ad_file = run_base / 'single_cell_data' / 'feature_matrix.h5'

adata = sc.read_h5ad(ad_file)
```


### Seurat (R)
```R
library('Seurat')
 
run_base = c('/path/to/g4x_output')

txcounts_path = file.path(run_base, "single_cell_data/cell_by_transcript.csv.gz")
metadata_path = file.path(run_base, "single_cell_data/cell_metadata.csv.gz")

counts <- read.csv(txcounts_path, row.names = 1) 
counts <- t(counts) # transpose to match Seurat input
 
metadata <- read.csv(metadata_path, row.names = 1) 
 
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
```