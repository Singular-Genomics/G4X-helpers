# G4X-helpers
Helper models and functions for working with G4X data.

# 📑 Table of Contents

- [G4X-helpers](#g4x-helpers)
- [📑 Table of Contents](#-table-of-contents)
- [🚀 Installation](#-installation)
  - [📦 Source Installation / CLI Usage](#-source-installation--cli-usage)
    - [Step 1: Prepare a Python environment](#step-1-prepare-a-python-environment)
    - [Step 2: Clone and install `g4x-helpers`](#step-2-clone-and-install-g4x-helpers)
    - [Step 3: Verify installation](#step-3-verify-installation)
- [🧬 `resegment`: Reprocess G4X output with a new segmentation](#-resegment-reprocess-g4x-output-with-a-new-segmentation)
  - [`resegment` Usage](#resegment-usage)
  - [`resegment` Argument Descriptions](#resegment-argument-descriptions)
    - [`--run_base` (required)](#--run_base-required)
    - [`--segmentation_mask` (required)](#--segmentation_mask-required)
    - [`--sample_id` (optional)](#--sample_id-optional)
    - [`--out_dir` (optional)](#--out_dir-optional)
    - [`--segmentation_mask_key` (optional)](#--segmentation_mask_key-optional)
    - [`--threads` (optional)](#--threads-optional)
    - [`--verbose` (optional)](#--verbose-optional)
- [🎨 `update_bin`: Update G4X-Viewer `.bin` file with new metadata](#-update_bin-update-g4x-viewer-bin-file-with-new-metadata)
  - [`update_bin` Usage](#update_bin-usage)
  - [`update_bin` Argument Descriptions](#update_bin-argument-descriptions)
    - [`--bin_file` (required)](#--bin_file-required)
    - [`--out_path` (required)](#--out_path-required)
    - [`--metadata` (required)](#--metadata-required)
    - [`--cellid_key` (optional)](#--cellid_key-optional)
    - [`--cluster_key` (optional)](#--cluster_key-optional)
    - [`--cluster_color_key` (optional)](#--cluster_color_key-optional)
    - [`--emb_key` (optional)](#--emb_key-optional)
    - [`--verbose` (optional)](#--verbose-optional-1)
- [🧪 `new_bin`: Generate G4X-Viewer `.bin` Files from G4X Sample Output](#-new_bin-generate-g4x-viewer-bin-files-from-g4x-sample-output)
  - [`new_bin` Usage](#new_bin-usage)
  - [`new_bin` Argument Descriptions](#new_bin-argument-descriptions)
    - [`--run_base` (required)](#--run_base-required-1)
    - [`--out_dir` (optional)](#--out_dir-optional-1)
    - [`--threads` (optional)](#--threads-optional-1)
    - [`--verbose` (optional)](#--verbose-optional-2)
- [📦 `tar_viewer`: Package G4X-Viewer Folder for Distribution](#-tar_viewer-package-g4x-viewer-folder-for-distribution)
  - [`tar_viewer` Usage](#tar_viewer-usage)
  - [`tar_viewer` Argument Description](#tar_viewer-argument-description)
    - [`--viewer_dir` (required)](#--viewer_dir-required)


# 🚀 Installation

## 📦 Source Installation / CLI Usage

`g4x-helpers` can be installed and run directly as a Python package.

### Step 1: Prepare a Python environment

- Install **conda**, **miniconda**, or **mamba**
- Create the environment:

```bash
conda create -n g4x-helpers_env python=3.10
```

- Activate the environment:

```bash
conda activate g4x-helpers_env
```

### Step 2: Clone and install `g4x-helpers`

- Clone the repository:

```bash
git clone git@github.com:Singular-Genomics/G4X-helpers.git
```

- Navigate to the repo directory:

```bash
cd g4x-helpers
```

- Install the package:

```bash
pip install .
```

### Step 3: Verify installation

After installation, you can call the following commands from any terminal and help statements should be printed:

```bash
resegment --help
update_bin --help
new_bin --help
tar_viewer --help
```
---

📘 **Note:**  
Examples of API (Python) usage can be found in the `api_demo.ipynb` notebook included in the repository.


# 🧬 `resegment`: Reprocess G4X output with a new segmentation

Replaces or updates the segmentation mask in a G4X run and regenerates all downstream single-cell data and `.bin` files.

## `resegment` Usage

```bash
resegment \
  --run_base /path/to/G4X/output \
  --segmentation_mask /path/to/new_mask.npz \
  [--sample_id SAMPLE123] \
  [--out_dir /path/to/output_dir] \
  [--segmentation_mask_key nuclei_exp] \
  [--threads 8] \
  [--verbose 2]
```
---

## `resegment` Argument Descriptions

### `--run_base` (required)
> **Type:** `str`  
> **Description:**  
> Path to the G4X sample output folder (the base directory for the run). This directory must contain required files such as `run_meta.json`, segmentation masks, and panel files.

### `--segmentation_mask` (required)
> **Type:** `str`  
> **Description:**  
> Path to the new segmentation mask file. Supported formats include `.npy`, `.npz`, and `.geojson`. This file will be used to replace the existing mask for transcript and protein signal assignment.

### `--sample_id` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Optional sample identifier. If not provided, the sample ID will be inferred from the name of the `run_base` directory.

### `--out_dir` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Directory to write the updated segmentation and downstream output files. If not provided, existing files in the `run_base` directory will be overwritten in-place.

### `--segmentation_mask_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Specifies the identifier for segmentation labels when loading mask data:
> - **If using a `.npz` file**:  
> Provide the name of the array within the archive that corresponds to the segmentation mask  
> (required if multiple arrays are stored).
> 
> - **If using a `.geojson` file**:  
> By default, cell labels are expected in a column named `label`.  
> Use this argument to override and select a different column as the label source.

### `--threads` (optional)
> **Type:** `int`  
> **Default:** `4`  
> **Description:**  
> Number of threads to use for segmentation, signal extraction, and downstream computation.

### `--verbose` (optional)
> **Type:** `int`  
> **Default:** `1`  
> **Description:**  
> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG
> This affects how much information is printed to the console during execution.

---

# 🎨 `update_bin`: Update G4X-Viewer `.bin` file with new metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

## `update_bin` Usage

```bash
update_bin \
  --bin_file /path/to/sample.bin \
  --out_path /path/to/updated_sample.bin \
  --metadata /path/to/metadata.csv \
  [--cellid_key CellID] \
  [--cluster_key Cluster] \
  [--cluster_color_key ClusterColor] \
  [--emb_key UMAP] \
  [--verbose 2]
```
---

## `update_bin` Argument Descriptions

### `--bin_file` (required)
> **Type:** `str`  
> **Description:**  
> Path to the existing `.bin` file used by the G4X-Viewer. This file will be updated with metadata from the provided CSV file.

### `--out_path` (required)
> **Type:** `str`  
> **Description:**  
> Output path where the updated `.bin` file will be saved. The directory will be created if it doesn't already exist.

### `--metadata` (required)
> **Type:** `str`  
> **Description:**  
> Path to a CSV file containing metadata to update the `.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

### `--cellid_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Name of the column in the metadata file that contains cell IDs matching those in the `.bin` file. If not provided, the first column in the metadata will be used.

### `--cluster_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster_color_key` is used.

### `--cluster_color_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster_key`.

### `--emb_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Prefix for embedding coordinates. The script will look for two columns: `{emb_key}_0` and `{emb_key}_1`. Used to embed cells in UMAP/tSNE/etc. space.

### `--verbose` (optional)
> **Type:** `int`  
> **Default:** `1`  
> **Description:**  
> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG
> This affects how much information is printed to the console during execution.

---

# 🧪 `new_bin`: Generate G4X-Viewer `.bin` Files from G4X Sample Output

This tool will create a new `.bin` segmentation file compatible with the G4X-Viewer using the processed output from a G4X run. This is typically only needed to update older outputs to newer versions of the `.bin` format.

---

## `new_bin` Usage

```bash
new_bin \
  --run_base /path/to/G4X_sample_output \
  [--out_dir /path/to/output_folder] \
  [--threads 8] \
  [--verbose 2]
```
---

## `new_bin` Argument Descriptions

### `--run_base` (required)
> **Type:** `str`  
> **Description:**  
> Path to the base directory containing G4X sample output. This folder should include `adata` and segmentation files generated by a previous G4X analysis run.

### `--out_dir` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Output directory where the new `.bin` file will be saved. If not provided, the file will be written in-place to a default path within the `run_base` directory (`g4x_viewer/{sample_id}.bin`).

### `--threads` (optional)
> **Type:** `int`  
> **Default:** `4`  
> **Description:**  
> Number of threads to use for processing. Increase this value to speed up the `.bin` generation process on multicore machines.

### `--verbose` (optional)
> **Type:** `int`  
> **Default:** `1`  
> **Description:**  
> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG  
> This affects how much information is printed to the console during execution.

---

# 📦 `tar_viewer`: Package G4X-Viewer Folder for Distribution

Tars a G4X-Viewer folder by validating and organizing key viewer assets (e.g., `.bin`, `.ome.tiff`, `.tar`, etc.), generating the required `dataset.config.json` metadata file, and creating a `.tar` archive ready for use with the Single-File upload option in the G4X-Viewer.

---

## `tar_viewer` Usage

```bash
tar_viewer \
  --viewer_dir /path/to/g4x_viewer_folder
```
---

## `tar_viewer` Argument Description

### `--viewer_dir` (required)
> **Type:** `str`  
> **Description:**  
> Path to the G4X-Viewer folder that contains the data and metadata to be packaged. This folder must contain:
> - A single `.bin` file
> - A `{sample_id}.ome.tiff` image file
> - A `{sample_id}_run_metadata.json` file
> - A `{sample_id}.tar` transcript file
> - A `{sample_id}_HE.ome.tiff` file (H&E), which will be moved into a subdirectory `h_and_e/`

---