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
    │   ├── <sample_id>.bin
    │   ├── <sample_id>.ome.tiff              
    │   ├── <sample_id>.tar
    │   ├── <sample_id>_HE.ome.tiff
    │   ├── <sample_id>_nuclear.ome.tiff
    │   └── <sample_id>_run_metadata.json
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
    ├── summary_<sample_id>.html
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
    │   ├── <sample_id>.bin
    │   ├── <sample_id>.tar
    │   ├── <sample_id>_HE.ome.tiff
    │   ├── <sample_id>_nuclear.ome.tiff
    │   └── <sample_id>_run_metadata.json
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
    ├── summary_<sample_id>.html
    └── transcript_panel.csv
    ```

<br>

## Sample sub-directory reference

---

### root of sample_folder 
> `run_meta.json:` JSON file containing versioning information for the panels used, analysis pipelines, and specific sequencer on which the experiment was run.

> `samplesheet.csv:` CSV file containing detailed run information. Details the experimental design, flow cell layout, tissue type, panel utilized, etc. This file is useful when performing analyses to identify which tissue blocks are which. 

> `summary_<sample_id>.html:` HTML file which gives a high level overview of the experiment outputs, quality, and performance for the individual tissue block. For more information, see Summary HTML (link). 

> `transcript_panel.csv:` CSV file containing a full list of all targeted genes in this experiment and the panel(s) which they originated from.

> `protein_panel.csv:` CSV file containing a full list of all targeted proteins in this experiment and the panel(s) which they originated from. *Multiomics runs only*



### /diagnostics/
> `transcript_table.parquet:` Parquet file containing all decoded and non-decoded transcripts and associated metadata (e.g. spatial coordinate, gene identity, cell identity (if assigned to a cell), quality score, sequence).



### /g4x_viewer/
> `<sample_id>.bin:` Binary file containing the segmentation mask for the stitched image. Can be easily read in Python with numpy. For visualizing the stored array, see Visualizing a Segmentation Mask (link).

> `<sample_id>.ome.tiff:` Multidimensional OME-TIFF image file. On windows, this may appear as `<sample_id>.ome`. This image contains aggregated images for all protein targets as well as nuclear stain. Can be loaded into any standard ome.tiff readers, including our G4X Viewer, and napri. For more information on our Viewer, see G4X Viewer (link)  *Multiomics runs only.*. 

> `<sample_id>_HE.ome.tiff:` OME-TIFF image file containing the fH&E stain images. Can be loaded into any standard OME-TIFF readers, including our G4X Viewer and napri. For more information on our Viewer, see G4X Viewer (link). 

> `<sample_id>_nuclear.ome.tiff:` OME-TIFF image file containing the nuclear stain images. Can be loaded into any standard OME-TIFF readers, including our G4X Viewer and napri. For more information on our Viewer, see G4X Viewer (link). 

> `<sample_id>_run_metadata.json:` JSON file containing much of the same information as the run_meta.json along with extra core metrics information (such as tissue area, total tx, etc).

> `<sample_id>.tar:` Tarball containing all other files from this directory bundled into one file. This can be loaded into the G4X Viewer directly with the “single file upload” option to avoid dragging each file individually. May take longer to load than the individual files due to needing to untar the components before displaying on the viewer. 



### /h_and_e/
> `eosin.jp2:` Full-sized eosin stained image used for analysis purposes for selected tissue block.

> `eosin_thumbnail.png:` Downsampled image from the .jp2 file for easier viewing of the eosin stain for selected tissue block.

> `h_and_e.jp2:` Full-sized fH&E image used for analysis purposes for selected tissue block.

> `h_and_e_thumbnail.jpg:` Downsampled image from the .jp2 file for easier viewing of the fH&E stain for selected tissue block.

> `nuclear.jp2:` Full-sized nuclear stain image used for analysis purposes for selected tissue block.

> `nuclear_thumbnail.png:` Downsampled image from the .jp2 file for easier viewing of the nuclear stain for selected tissue block.
  


### /metrics/
> `core_metrics.csv:` CSV file containing a set of core metrics for the tissue block including total transcripts, total area, number of cells and more.

> `per_area_metrics.csv:` CSV file containing a set of per-area metrics for the tissue block (coordinate location, nunmber of transcripts, and number of cells), separated out into images from before the images were stitched together into one whole block. 



### /protein/ *(only when protein is included)*
> `<protein_name>.jp2:` Full-sized image used for analysis purposes. Shows the `<protein_name>` stain for selected tissue block.

> `<protein_name>_thumbnail.png:` Downsampled image of the .jp2 file for easier viewing. Shows the `<protein_name>` stain for selected tissue block.




### /rna/
> `transcript_table.csv.gz:` CSV file containing a transcript table showing all transcripts identified on the whole tissue block. Contains coordinate information, z-layer, gene identity, and cell_id fields. All transcript here are high confidence transcripts post-filtering and processing.



### /segmentation/
> `segmentation_mask.npz:` Compressed numpy file containing the segmentation mask array. This can be easily read with the numpy.load() function. An example of how to visualize and plot this can be found here (link).



### /single_cell_data/
> `cell_by_protein.csv.gz:` Gzipped CSV file in a cell x protein intensity format. Each entry in the table is the average protein intensity for a given protein in a given cell. *Multiomics runs only.*

> `cell_by_transcript.csv.gz:` Gzipped CSV file in a cell x transcript format. Each entry in the table is the counts for a given transcript in a given cell.

> `cell_metadata.csv.gz:` Gzipped CSV file containing the metadata associated with each cell, including cell_id, ... This is needed to laund a Seurat object and perform downstream analyses. For more information, see `(link)`.

??? note "Expand to see column descriptions"
    | name   | description   |
    |------|----------|
    | `label`  | Cell ID  |
    | `<protein>_intensity_mean` | Fluorescence intensity mean for a given protein |
    | `cell_id` | Cell ID for a given cell |
    | `cell_x/y` | Spatial X/Y coordinate for the nuclear segmentation centroid |
    | `expanded_cell_x/y` | Spatial X/Y coordinate for the expanded nuclear segmentation centroid |
    | `log1p_n_genes_by_counts` | Log number of unique genes detected |
    | `log1p_total_counts` | Log number of total transcripts |
    | `n_genes_by_counts` | Number of unique genes detected |
    | `nuclei_area` | Area of the nuclear segmentation |
    | `nuclei_expanded_area` | Area of the expanded nuclear segmentation |
    | `total_counts` | Total transcript counts |

> `clustering_umap.csv.gz:` Gzipped CSV file containing the matrix of cell cluster annotations and UMAP coordinates for each cell. This is used to visualize the clustering of cells in 2D for a given leiden resolution or UMAP embedding setting.

??? note "Expand to see column descriptions"
    | column   | description   |
    |-------|----------|
    | `label`  | Cell ID  |
    | `leiden_<resolution>`  | Leiden cluster identity for the cell at the specified resolution (0.2-1.0)  |
    | `X_umap_<min_dist>_<spread>_<axis>`  | UMAP coordinate for the given cell with the given min_dist and spread for a given axis (1 is typically x/UMAP1, 2 is typically y/UMAP2) |

> `dgex.csv.gz:` Gzipped CSV file containing the differential gene expression (DGEx) results for the selected tissue block. Columns are detailed below.
??? note "Expand to see column descriptions"
    | column   | description   |
    |-------|----------|
    | `names`  | Gene symbol  |
    | `scores`  | Z-score from Wilcoxon rank-sum test  |
    | `logfoldchanges`  | LogFoldChange for the given cluster compared to *all* other clusters combined |
    | `pvals`  | P-value from Wilcoxon rank-sum test |
    | `pvals_adj`  | Adjusted P-value  |
    | `pct_nz_group`  | Percentage of non-zero values in the given cluster.  |
    | `pct_nz_reference`  | Percentage of non-zero values in all cells outside the given cluster.  |
    | `group`  | Leiden cluster identity  |
    | `leiden_res`  | Leiden clustering resolution that this entry is derived from (1 per gene per cluster) |

> `feature_matrix.h5:` H5ad file containing the full cell by gene matrix as well as a wide array of metadata and annotations you might want to use for downstream analysis. This file can be loaded into Python and run through scanpy or a number of other pipelines. See `data import`. The annotations in this file are detailed below.

??? note "Expand to see annotation layer descriptions"
    | name   | layer  | dimensions   | description   |
    |------|----|----------|----------|
    | `<protein>_intensity_mean` | obs | ncell x  1 | Fluorescence intensity mean for a given protein per cell |
    | `cell_id` | obs | ncell x 1 | Cell ID for a given cell |
    | `cell_x/y` | obs | ncell x 1 | Spatial X/Y coordinate for the nuclear segmentation centroid |
    | `expanded_cell_x/y` | obs | ncell x 1 | Spatial X/Y coordinate for the expanded nuclear segmentation centroid per cell |
    | `log1p_n_genes_by_counts` | obs | ncell x 1 | Log number of unique genes detected per cell |
    | `log1p_total_counts` | obs  | ncell x 1 | Log number of total transcripts per cell |
    | `n_genes_by_counts` | obs  | ncell x 1 | Number of unique genes detected per cell |
    | `nuclei_area` | obs | ncell x 1 | Area of the nuclear segmentation per cell |
    | `nuclei_expanded_area` | obs | ncell x 1  | Area of the expanded nuclear segmentation per cell |
    | `total_counts` | obs | ncell x 1 | Total transcript counts per cell |
    | `gene_id` | var | ngenes x 1 | Gene symbol |
    | `log1p_mean_counts` | var | ngenes x 1 | Log mean transcript counts across all cells |
    | `log1p_total_counts` | var | ngenes x 1 | Log total transcript counts across all cells |
    | `mean_counts` | var | ngenes x 1 | Mean counts of each transcript across all cells |
    | `modality` | var | ngenes x 1 | G4X modality (transcript or protein) |
    | `n_cells_by_counts` | var | ngenes x 1 | Number of cells with counts of each transcript |
    | `pct_dropout_by_counts` | var | ngenes x 1 |  |
    | `probe_type` | var | ngenes x 1 | Type of probe: Negative control probe/sequence (NCP/NCS) or tanscript targeting (targeting) |
    | `total_counts` | var | ngenes x 1 | Total transcript counts per gene |




--8<-- "_partials/end_cap.md"
