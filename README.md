# G4X-helpers
Helper models and functions for working with G4X data.

# 🧬 `resegment`: Reprocess G4X Output with a New Segmentation

Replaces or updates the segmentation mask in a G4X run and regenerates all downstream single-cell data and `.bin` files.

## ✅ Usage

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

# 🎨 `update_bin`: Update G4X-Viewer `.bin` File Metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

## ✅ Usage

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


# 🧬 `resegment` – Argument Descriptions

### `--run_base` (required)
> **Type:** `str`  
> **Description:**  
> Path to the G4X sample output folder (the base directory for the run). This directory must contain required files such as `run_meta.json`, segmentation masks, and panel files.

---

### `--segmentation_mask` (required)
> **Type:** `str`  
> **Description:**  
> Path to the new segmentation mask file. Supported formats include `.npy`, `.npz`, and `.geojson`. This file will be used to replace the existing mask for transcript and protein signal assignment.

---

### `--sample_id` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Optional sample identifier. If not provided, the sample ID will be inferred from the name of the `run_base` directory.

---

### `--out_dir` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Directory to write the updated segmentation and downstream output files. If not provided, existing files in the `run_base` directory will be overwritten in-place.

---

### `--segmentation_mask_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> If the segmentation mask is provided as a `.npz` file, this key specifies which array inside the archive should be used as the mask.

---

### `--threads` (optional)
> **Type:** `int`  
> **Default:** `4`  
> **Description:**  
> Number of threads to use for segmentation, signal extraction, and downstream computation.

---

### `--verbose` (optional)
> **Type:** `int`  
> **Default:** `1`  
> **Description:**  
> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG

---

# 🎨 `update_bin` – Argument Descriptions

### `--bin_file` (required)
> **Type:** `str`  
> **Description:**  
> Path to the existing `.bin` file used by the G4X-Viewer. This file will be updated with metadata from the provided CSV file.

---

### `--out_path` (required)
> **Type:** `str`  
> **Description:**  
> Output path where the updated `.bin` file will be saved. The directory will be created if it doesn't already exist.

---

### `--metadata` (required)
> **Type:** `str`  
> **Description:**  
> Path to a CSV file containing metadata to update the `.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

---

### `--cellid_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Name of the column in the metadata file that contains cell IDs matching those in the `.bin` file. If not provided, the first column in the metadata will be used.

---

### `--cluster_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster_color_key` is used.

---

### `--cluster_color_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster_key`.

---

### `--emb_key` (optional)
> **Type:** `str`  
> **Default:** `None`  
> **Description:**  
> Prefix for embedding coordinates. The script will look for two columns: `{emb_key}_0` and `{emb_key}_1`. Used to embed cells in UMAP/tSNE/etc. space.

---

### `--verbose` (optional)
> **Type:** `int`  
> **Default:** `1`  
> **Description:**  
> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG



# 🚀 Installation

## ⚙️ Installation via `pyproject.toml`

If installed as a Python package with CLI entry points:

```toml
[project.scripts]
resegment = "g4x_helpers.entrypoint:launch_resegment"
update_bin = "g4x_helpers.entrypoint:launch_update_bin"
```

After installation, you can call:

```bash
resegment --help
update_bin --help
```

from any terminal to use these tools.

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

- Change into the repo directory:

```bash
cd g4x-helpers
```

- Install the package:

```bash
pip install .
```

- Or for an editable install (recommended for development):

```bash
pip install -e .
```

---

📘 **Note:**  
Examples of API (Python) usage can be found in the `api_demo.ipynb` notebook included in the repository.
