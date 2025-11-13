<br>

# `tar_viewer`
#### Package G4X-viewer Folder for Distribution

Tars a G4X-viewer folder by validating and organizing key viewer assets (e.g., `.bin`, `.ome.tiff`, `.tar`, etc.), generating the required `dataset.config.json` metadata file, and creating a `.tar` archive ready for use with the Single-File upload option in the G4X-viewer.


## Usage
![`g4x-helpers tar_viewer --help`](../img/tar_viewer-help.svg)


## argument description
---

#### `--g4x-data`: (*type:* `str`)

> Path to the G4X sample output directory. By default the command packages the `g4x_viewer` folder inside this directory.

#### `--viewer-dir`: (*type:* `str`  *default:* `None`)

> Optional override for the viewer folder that should be packaged. When omitted, the command uses the `g4x_viewer` directory found within `--g4x-data`. The folder must contain:
> - A single `.bin` file
> - A `{sample_id}.ome.tiff` image file
> - A `{sample_id}_run_metadata.json` file
> - A `{sample_id}.tar` transcript file
> - A `{sample_id}_HE.ome.tiff` file (H&E), which will be moved into a subdirectory `h_and_e/`

#### `--output`: (*type:* `str`  *default:* `None`)

> Directory where the packaged archive should be written. Defaults to the `g4x_viewer` folder within `--g4x-data`.

#### `--sample-id`: (*type:* `str`  *default:* `None`)

> Optional sample identifier override used when naming the generated archive.

<br>
--8<-- "_core/_partials/global_options_note.md"
<br>
--8<-- "_core/_partials/end_cap.md"
