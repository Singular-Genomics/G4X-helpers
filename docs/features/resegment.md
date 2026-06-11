<br>

# `resegment`
#### Reprocess G4X-output with a new segmentation

Takes new cell-labels from a custom segmentation output and re-assigns transcripts  
 and protein signals to those cells. The operation recreates single-cell outputs and 
 initializes a new G4X-viewer zarr store. Does not regenerate metrics.

---

## Usage
![`g4x-helpers resegment --help`](../img/resegment-help.svg)

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--cell-labels` [required]
_type_ : <span class="acc-2-code">`file path`</span>  
_example_  : `path/to/segmentation.npz`
> Path to a new segmentation mask file containing cell labels. It will be used to create new single-cell outputs by aggregating transcript and protein data on those labels. The extent of the labels must match the shape of the original data.  
> Supported formats are: `(.npy, .npz, .geojson)`

---

### `--labels-key`
_type_ : <span class="acc-2-code">`string`</span>  
_example_  : `cell_id`
> Specifies the identifier for segmentation labels when loading mask data:

> **If using a `.npz` file**:  
> Provide the name of the array within the archive that corresponds to the segmentation mask (required if multiple arrays are stored).

> **If using a `.geojson` file**:  
> By default, cell labels are expected in a column named `label`.  
> Use this argument to override and select a different column as the label source.

---

--8<-- "_partials/arg_no_downstream.md"

--8<-- "_partials/arg_branch.md"

<br>
--8<-- "_core/_partials/end_cap.md"
