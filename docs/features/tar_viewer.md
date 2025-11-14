<br>

# `tar_viewer`
#### Package G4X-viewer Folder for Distribution

Tars a G4X-viewer folder by validating and organizing key viewer assets (e.g., `.bin`, `.ome.tiff`, `.tar`, etc.), generating the required `dataset.config.json` metadata file, and creating a `.tar` archive ready for use with the Single-File upload option in the G4X-viewer.


## Usage
![`g4x-helpers tar_viewer --help`](../img/tar_viewer-help.svg)

--8<-- "_core/_partials/global_options_note.md"

## argument description
---
--8<-- "_core/_partials/arg_g4x_data.md"


#### `--viewer-dir`: (*type:* `str`  *default:* `None`)

> Optional override for the viewer folder that should be packaged. When omitted, the command uses the `g4x_viewer` directory found within `--g4x-data`.  

> The folder must contain:  
> - `{sample_id}.bin`  
> - `{sample_id}.ome.tiff`  
> - `{sample_id}_run_metadata.json`  
> - `{sample_id}.tar`  
> - `{sample_id}_HE.ome.tiff`  

--8<-- "_core/_partials/arg_output.md"

--8<-- "_core/_partials/arg_smp_id.md"

<br>
--8<-- "_core/_partials/end_cap.md"
