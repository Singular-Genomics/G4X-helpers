<br>

# G4X output structure

### Sample directory tree
---
Directory structure depends on run type.

=== "Transcript & Protein"

    ```
    <sample_root>
    │
    ├── diagnostics
    │   └── transcript_table.parquet
    │
    ├── g4x_viewer 
    │   ├── <sample_1>.bin
    │   ├── <sample_1>.ome.tiff              
    │   ├── <sample_1>.tar
    │   ├── <sample_1>_HE.ome.tiff
    │   ├── <sample_1>_nuclear.ome.tiff
    │   └── <sample_1>_run_metadata.json
    │
    ├── h_and_e 
    │   ├── eosin.jp2
    │   ├── eosin_thumbnail.png
    │   ├── h_and_e.jp2
    │   ├── h_and_e_thumbnail.jpg
    │   ├── nuclear.jp2
    │   └── nuclear_thumbnail.png
    │
    ├── metrics 
    │   ├── core_metrics.csv
    │   └── per_area_metrics.csv
    │
    ├── protein                             
    │   ├── <protein_1>.jp2
    │   ├── <protein_1>_thumbnail.png
    │   ├── <protein_2>.jp2
    │   ├── <protein_2>_thumbnail.png
    │   └── …
    │
    ├── protein_panel.csv                   
    │
    ├── rna
    │   └── transcript_table.csv.gz
    │
    ├── run_meta.json
    │
    ├── samplesheet.csv
    │
    ├── segmentation
    │   └── segmentation_mask.npz
    │
    ├── single_cell_data
    │   ├── cell_by_protein.csv.gz          
    │   ├── cell_by_transcript.csv.gz
    │   ├── cell_metadata.csv.gz
    │   ├── clustering_umap.csv.gz
    │   ├── dgex.csv.gz
    │   └── feature_matrix.h5
    │
    ├── summary_<sample_1>.html
    └── transcript_panel.csv
    ```

=== "Transcript"

    ```
    <sample_root>
    │
    ├── diagnostics
    │   └── transcript_table.parquet
    │
    ├── g4x_viewer
    │   ├── <sample_1>.bin
    │   ├── <sample_1>.tar
    │   ├── <sample_1>_HE.ome.tiff
    │   ├── <sample_1>_nuclear.ome.tiff
    │   └── <sample_1>_run_metadata.json
    │
    ├── h_and_e
    │   ├── eosin.jp2
    │   ├── eosin_thumbnail.png
    │   ├── h_and_e.jp2
    │   ├── h_and_e_thumbnail.jpg
    │   ├── nuclear.jp2
    │   └── nuclear_thumbnail.png
    │
    ├── metrics
    │   ├── core_metrics.csv
    │   └── per_area_metrics.csv
    │
    ├── rna
    │   └── transcript_table.csv.gz
    │
    ├── run_meta.json
    │
    ├── samplesheet.csv
    │
    ├── segmentation
    │   └── segmentation_mask.npz
    │
    ├── single_cell_data
    │   ├── cell_by_transcript.csv.gz
    │   ├── cell_metadata.csv.gz
    │   ├── clustering_umap.csv.gz
    │   ├── dgex.csv.gz
    │   └── feature_matrix.h5
    │
    ├── summary_<sample_1>.html
    └── transcript_panel.csv
    ```

<br>

## Sample sub-directory reference

---

### root of sample_folder 
> `run_meta.json`  
> `samplesheet.csv`  
> `summary_<sample_1>.html`  
> `transcript_panel.csv`  
> `protein_panel.csv`  *(only when protein is included)*



### diagnostics
> `transcript_table.parquet`



### g4x_viewer
> `<sample_1>.bin`  
> `<sample_1>.ome.tiff` *(only when protein is included)*  
> `<sample_1>.tar`  
> `<sample_1>_HE.ome.tiff`  
> `<sample_1>_nuclear.ome.tiff`  
> `<sample_1>_run_metadata.json`  



### h_and_e
> `eosin.jp2`  
> `eosin_thumbnail.png`  
> `h_and_e.jp2`  
> `h_and_e_thumbnail.jpg`  
> `nuclear.jp2`  
> `nuclear_thumbnail.png`  



### metrics
> `core_metrics.csv`  
> `per_area_metrics.csv`  



### protein *(only when protein is included)*
> `<protein_1>.jp2`  
> `<protein_1>_thumbnail.png`  
> `<protein_2>.jp2`  
> `<protein_2>_thumbnail.png`  
> `(additional proteins)`



### rna
> `transcript_table.csv.gz`  



### segmentation
> `segmentation_mask.npz`  



### single_cell_data
> `cell_by_protein.csv.gz`  *(only when protein is included)*  
> `cell_by_transcript.csv.gz`  
> `cell_metadata.csv.gz`  
> `clustering_umap.csv.gz`  
> `dgex.csv.gz`  
> `feature_matrix.h5`  


--8<-- "_core/_partials/end_cap.md"
