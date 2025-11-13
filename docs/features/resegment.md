<br>

# `resegment`
#### Reprocess G4X-output with a new segmentation

Replaces or updates the segmentation mask in a G4X run and regenerates all downstream single-cell data and `.bin` files.

## Usage
![`g4x-helpers resegment --help`](../img/resegment-help.svg)

## argument descriptions
---
### required
#### `--g4x-data`: (*type:* `str`)

> Path to the G4X sample output directory. This folder must contain required files such as `run_meta.json`, segmentation masks, panel files, and feature tables.

#### `--cell-labels`: (*type:* `str`)

> Path to the new segmentation mask file. Supported formats include `.npy`, `.npz`, and `.geojson`. This file will be used to replace the existing mask for transcript and protein signal assignment.

### optional 
#### `--sample-id`: (*type:* `str`  *default:* `None`)

> Optional sample identifier. If not provided, it is inferred from the name of the `--g4x-data` directory.

#### `--output`: (*type:* `str`  *default:* `None`)

> Directory that receives the updated segmentation and downstream output. When omitted, the command modifies artifacts inside the `--g4x-data` folder in place.

#### `--labels-key`: (*type:* `str`  *default:* `None`)

> Specifies the identifier for segmentation labels when loading mask data:

> - **If using a `.npz` file**:  
> Provide the name of the array within the archive that corresponds to the segmentation mask (required if multiple arrays are stored).

> - **If using a `.geojson` file**:  
> By default, cell labels are expected in a column named `label`.  
> Use this argument to override and select a different column as the label source.

<br>
--8<-- "_core/_partials/global_options_note.md"
<br>
--8<-- "_core/_partials/end_cap.md"
