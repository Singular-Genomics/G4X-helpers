#  `--tar_viewer`
# Package G4X-Viewer Folder for Distribution

Tars a G4X-Viewer folder by validating and organizing key viewer assets (e.g., `.bin`, `.ome.tiff`, `.tar`, etc.), generating the required `dataset.config.json` metadata file, and creating a `.tar` archive ready for use with the Single-File upload option in the G4X-Viewer.

## Usage

```bash
$ tar_viewer --viewer_dir /path/to/g4x_viewer_folder
```

<br>

## Argument Descriptions
---

### `--viewer_dir` (required)
> **Type:** `str`  
> **Description:**  
> Path to the G4X-Viewer folder that contains the data and metadata to be packaged. This folder must contain:
> - A single `.bin` file
> - A `{sample_id}.ome.tiff` image file
> - A `{sample_id}_run_metadata.json` file
> - A `{sample_id}.tar` transcript file
> - A `{sample_id}_HE.ome.tiff` file (H&E), which will be moved into a subdirectory `h_and_e/`

<br>
