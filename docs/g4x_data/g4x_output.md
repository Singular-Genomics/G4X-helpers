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
    |   ├── protein_core_metrics.csv
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

> `run_meta.json:` JSON file containing versioning information for the panels used and analysis pipelines as well as the ID for the sequencer on which the experiment was run.

> `samplesheet.csv:` CSV file containing detailed run information. Details the experimental design, flow cell layout, tissue type, panel utilized, etc. This file is useful for analysis, as it designates where each tissue section is positioned on the flow cell.

> `summary_<sample_id>.html:` HTML file which gives a high level overview of the experiment outputs, data quality, and performance for the selected tissue block.

> `transcript_panel.csv:` CSV file containing a full list of all targeted genes in this experiment and the panel(s) which they originated from.

> `protein_panel.csv:` CSV file containing a full list of all targeted proteins in this experiment and the panel(s) which they originated from. *Multiomics runs only*

<br>

### /diagnostics/

> `transcript_table.parquet:` Parquet file containing all decoded and non-decoded transcripts and associated metadata (e.g. spatial coordinate, gene identity, cell identity (if assigned to a cell), quality score, sequence). Parquet files can be loaded in Python using the polars, fastparquet, pandas, and pyarrow packages.

??? note "Expand to see column descriptions"
    | column | type | description |
    |-------|-----|----------|
    | `x_coord_shift` | *float* | The x coordinate for the transcript (shifted to global coordinates) |
    | `y_coord_shift` | *float* | The y coordinate for the transcript (shifted to global coordinates) |
    | `z` | *int* | z layer of identified transcript |
    | `demuxed` | *bool* | Whether or not the transcript was demultiplexed |
    | `transcript_condensed` | *str* | Shortened name of transcript |
    | `meanQS` | *float* | Mean quality score for the transcript |
    | `cell_id` | *uint* | Cell ID  |
    | `sequence_to_demux` | *str* | Sequence identified that will be demultiplexed |
    | `transcript` | *str* | Long form transcript name (specific to single probe) |
    | `TXUID` | *str* | Unique identifier for the transcript |

<br> 

### /g4x_viewer/

!!! tip
    The G4X Viewer is a web-based tool for visualizing and exploring G4X data. All files in this directory are designed to be loaded into and explored with the G4X Viewer. For more information on how to use the G4X Viewer, see [G4X Viewer](../g4x_viewer/index.md).

> `<sample_id>.bin:` Binary file containing the segmentation mask for the stitched image. Can be easily read in Python with numpy.

> `<sample_id>.ome.tiff:` Multidimensional OME-TIFF image file. On windows, this may appear as `<sample_id>.ome`. This image contains aggregated images for all protein targets as well as nuclear stain. Can be loaded into any standard ome.tiff readers, including our G4X Viewer, and napari. *Multiomics runs only.*.

> `<sample_id>_HE.ome.tiff:` OME-TIFF image file containing the fH&E stain images. Can be loaded into any standard OME-TIFF readers, including our G4X Viewer and napari.

> `<sample_id>_nuclear.ome.tiff:` OME-TIFF image file containing the nuclear stain images. Can be loaded into any standard OME-TIFF readers, including our G4X Viewer and napari.

> `<sample_id>_run_metadata.json:` JSON file containing much of the same information as the run_meta.json along with extra core metrics information (such as tissue area, total tx, etc).

> `<sample_id>.tar:` Tarball containing all other files from this directory bundled into one file. This can be loaded into the G4X Viewer directly with the “single file upload” option to avoid dragging each file individually. May take longer to load than the individual files due to needing to untar the components before displaying on the Viewer.

<br>

### /metrics/

> `core_metrics.csv:` CSV file containing a set of core metrics for the tissue block including total transcripts, total area, number of cells and more.

> `protein_core_metrics.csv:` CSV file containing a set of core protein metrics for the tissue block including SNR, background intensity, and Fisher's exact scores for the coocurrence of the protein signal with its associated transcript signal (`<protein>_fisher_score`) and a random background (`<protein>_fisher_score_background`). These scores indicate the likelihood of the signal being true signal compared to the measured background. *Multiomics runs only.*

> `per_area_metrics.csv:` CSV file containing a set of per-area metrics for the tissue block (coordinate location, number of transcripts, and number of cells), separated out into images from before the images were stitched together into one whole block.

<br>

### /h_and_e/

!!! tip
    The .jp2 images in this folder and the `/protein/` folder are suitable to use for both nuclear and cytoplasmic segmentation. For more information on how you might do this, see [segment data](../g4x_tutorials/segment_data.md).

> `eosin.jp2:` Full-sized eosin stained JPEG image used for analysis purposes for selected tissue block.

> `eosin_thumbnail.png:` Downsampled PNG image from the .jp2 file for easier viewing of the eosin stain for selected tissue block.

> `h_and_e.jp2:` Full-sized fH&E JPEG image used for analysis purposes for selected tissue block.

> `h_and_e_thumbnail.jpg:` Downsampled PNG image from the .jp2 file for easier viewing of the fH&E stain for selected tissue block.

> `nuclear.jp2:` Full-sized nuclear stained JPEG image used for analysis purposes for selected tissue block.

> `nuclear_thumbnail.png:` Downsampled PNG image from the .jp2 file for easier viewing of the nuclear stain for selected tissue block.

<br>

### /protein/ *(Only in Multiomics runs)*

> `<protein_name>.jp2:` Full-sized JPEG image used for analysis purposes. Shows the `<protein_name>` stain for selected tissue block.

> `<protein_name>_thumbnail.png:` Downsampled PNG image of the .jp2 file for easier viewing. Shows the `<protein_name>` stain for selected tissue block.

<br>

### /rna/

> `transcript_table.csv.gz:` CSV file containing a transcript table showing all demuxed transcripts on the whole tissue block. Contains coordinate information, z-layer, gene identity, and cell_id fields. All transcripts here are high confidence transcripts post-filtering and processing.

<br>

### /segmentation/

> `segmentation_mask.npz:` Compressed numpy array file containing the segmentation mask. This can be easily read with the `numpy.load()` function.

<br>

### /single_cell_data/

> `cell_by_protein.csv.gz:` Gzipped CSV file in a cell x protein intensity format. Each entry in the table is the average protein intensity for a given protein in a given cell. *Multiomics runs only.*

> `cell_by_transcript.csv.gz:` Gzipped CSV file in a cell x transcript format. Each entry in the table is the counts for a given transcript in a given cell.

> `cell_metadata.csv.gz:` Gzipped CSV file containing the metadata associated with each cell, including cell_id, protein mean intensity, and transcript counts per cell, per transcript species. This is needed to launch a Seurat object and perform downstream analyses. For more information, see [data import](./data_import.md).

??? note "Expand to see column descriptions"
    | name | type | description |
    |------|-----|----------|
    | `label` | *str* | Cell ID |
    | `<protein>_intensity_mean` | *float* | Fluorescence intensity mean for a given protein |
    | `cell_id` | *str* | Cell ID for a given cell |
    | `cell_x/y` | *float* | Spatial X/Y coordinate for the nuclear segmentation centroid |
    | `expanded_cell_x/y` | *float* | Spatial X/Y coordinate for the expanded nuclear segmentation centroid |
    | `log1p_n_genes_by_counts` | *float* | Log number of unique genes detected |
    | `log1p_total_counts` | *float* | Log number of total transcripts |
    | `n_genes_by_counts` | *int* | Number of unique genes detected |
    | `nuclei_area` | *int* | Area of the nuclear segmentation |
    | `nuclei_expanded_area` | *int* | Area of the expanded nuclear segmentation |
    | `total_counts` | *int* | Total transcript counts |

> `clustering_umap.csv.gz:` Gzipped CSV file containing the matrix of cell cluster annotations and UMAP coordinates for each cell. This is used to visualize the clustering of cells in 2D for a given leiden resolution or UMAP embedding setting.

??? note "Expand to see column descriptions"
    | column | type | description |
    |-------|-----|----------|
    | `label` | *str* | Cell ID  |
    | `leiden_<resolution>` | *int* | Leiden cluster identity for the cell at the specified resolution (0.2-1.0)  |
    | `X_umap_<min_dist>_<spread>_<axis>` | *float* | UMAP coordinate for the given cell with the given min_dist and spread for a given axis (1 is typically x/UMAP1, 2 is typically y/UMAP2) |

> `dgex.csv.gz:` Gzipped CSV file containing the differential gene expression (DGEx) results for the selected tissue block. Columns are detailed below.

??? note "Expand to see column descriptions"
    | column  | type | description  |
    |-------|-----|----------|
    | `names`  | *str* | Gene symbol  |
    | `scores`  | *float* | Z-score from Wilcoxon rank-sum test  |
    | `logfoldchanges`  | *float* | LogFoldChange for the given cluster compared to *all* other clusters combined |
    | `pvals`  | *float* | P-value from Wilcoxon rank-sum test |
    | `pvals_adj`  | *float* | Adjusted P-value  |
    | `pct_nz_group`  | *float* | Percentage of non-zero values in the given cluster.  |
    | `pct_nz_reference`  | *float* | Percentage of non-zero values in all cells outside the given cluster.  |
    | `group`  | *int* | Leiden cluster identity  |
    | `leiden_res`  | *str* | Leiden clustering resolution that this entry is derived from (1 per gene per cluster) |

> `feature_matrix.h5:` H5ad file containing the full cell by gene matrix as well as a wide array of metadata and annotations you might want to use for downstream analysis. This file can be loaded into Python and run through scanpy or a number of other pipelines. See [data import](./data_import.md). The annotations in this file are detailed below.

??? note "Expand to see annotation layer descriptions"
    | name   | layer  | *str* | dimensions   | description   |
    |------|----|-----|----------|----------|
    | `<protein>_intensity_mean` | *float* | obs | ncell x  1 | Fluorescence intensity mean for a given protein per cell |
    | `cell_id` | *str* | obs | ncell x 1 | Cell ID for a given cell |
    | `cell_x/y` | *float* | obs | ncell x 1 | Spatial X/Y coordinate for the nuclear segmentation centroid |
    | `expanded_cell_x/y` | *float* | obs | ncell x 1 | Spatial X/Y coordinate for the expanded nuclear segmentation centroid per cell |
    | `log1p_n_genes_by_counts` | *float* | obs | ncell x 1 | Log number of unique genes detected per cell |
    | `log1p_total_counts` | *float* | obs  | ncell x 1 | Log number of total transcripts per cell |
    | `n_genes_by_counts` | *int* | obs  | ncell x 1 | Number of unique genes detected per cell |
    | `nuclei_area` | *float* | obs | ncell x 1 | Area of the nuclear segmentation per cell |
    | `nuclei_expanded_area` | *float* | obs | ncell x 1  | Area of the expanded nuclear segmentation per cell |
    | `total_counts` | *int* | obs | ncell x 1 | Total transcript counts per cell |
    | `gene_id` | *str* | var | ngenes x 1 | Gene symbol |
    | `log1p_mean_counts` | *float* | var | ngenes x 1 | Log mean transcript counts across all cells |
    | `log1p_total_counts` | *float* | var | ngenes x 1 | Log total transcript counts across all cells |
    | `mean_counts` | *float* | var | ngenes x 1 | Mean counts of each transcript across all cells |
    | `modality` | *str* | var | ngenes x 1 | G4X modality (transcript or protein) |
    | `n_cells_by_counts` | *int* | var | ngenes x 1 | Number of cells with counts of each transcript |
    | `pct_dropout_by_counts` | *float* | var | ngenes x 1 |  |
    | `probe_type` | *str* | var | ngenes x 1 | Type of probe: Negative control probe/sequence (NCP/NCS) or transcript targeting (targeting) |
    | `total_counts` | *int* | var | ngenes x 1 | Total transcript counts per gene |
 
 <br>

 --8<-- "_partials/end_cap.md"
