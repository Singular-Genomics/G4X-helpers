<br>

# `resegment`
#### Reprocess G4X-output with a new segmentation

Replaces or updates the segmentation mask in a G4X run and regenerates all downstream single-cell data and `.bin` files.

---

## Usage
![`g4x-helpers resegment --help`](../img/resegment-help.svg)

--8<-- "_partials/global_options_note.md"

## Arguments
---
### required
--8<-- "_partials/arg_g4x_data.md"

#### `--cell-labels`: (*type:* `str`)

> Path to the new segmentation mask file. Supported formats include `.npy`, `.npz`, and `.geojson`. This file will be used to replace the existing mask for transcript and protein signal assignment.

### optional 

#### `--labels-key`: (*type:* `str`  *default:* `None`)

> Specifies the identifier for segmentation labels when loading mask data:

> - **If using a `.npz` file**:  
> Provide the name of the array within the archive that corresponds to the segmentation mask (required if multiple arrays are stored).

> - **If using a `.geojson` file**:  
> By default, cell labels are expected in a column named `label`.  
> Use this argument to override and select a different column as the label source.

--8<-- "_partials/arg_in_place.md"

<br>
--8<-- "_core/_partials/end_cap.md"
