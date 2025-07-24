<br>

# `--resegment`
#### Reprocess G4X output with a new segmentation

Replaces or updates the segmentation mask in a G4X run and regenerates all downstream single-cell data and `.bin` files.

--8<-- "_partials/section_break.md"

## usage

``` bash
$ resegment
  --run_base /path/to/G4X/output 
  --segmentation_mask /path/to/new_mask.npz 
  
  # ─── optional ───
  --sample_id <sample_id> 
  --out_dir <output_dir> 
  --segmentation_mask_key <mask_array_name> 
  --threads <n_threads> 
  --verbose <level>
```

## Argument Descriptions

### required
---
#### `--run_base`: (*type:* `str`)
 
> Path to the G4X sample output folder (the base directory for the run). This directory must contain required files such as `run_meta.json`, segmentation masks, and panel files.

#### `--segmentation_mask`: (*type:* `str`)

> Path to the new segmentation mask file. Supported formats include `.npy`, `.npz`, and `.geojson`. This file will be used to replace the existing mask for transcript and protein signal assignment.

### optional
---
#### `--sample_id`: (*type:* `str`  *default:* `None`)

> Optional sample identifier. If not provided, the sample ID will be inferred from the name of the `run_base` directory.

#### `--out_dir`: (*type:* `str`  *default:* `None`)

> Directory to write the updated segmentation and downstream output files. If not provided, existing files in the `run_base` directory will be overwritten in-place.

#### `--segmentation_mask_key`: (*type:* `str`  *default:* `None`)

> Specifies the identifier for segmentation labels when loading mask data:

> - **If using a `.npz` file**:  
> Provide the name of the array within the archive that corresponds to the segmentation mask (required if multiple arrays are stored).
 
> - **If using a `.geojson` file**:  
> By default, cell labels are expected in a column named `label`.  
> Use this argument to override and select a different column as the label source.

#### `--threads`: (*type:* `int`  *default:* `4`)

> Number of threads to use for segmentation, signal extraction, and downstream computation.

#### `--verbose`: (*type:* `int`  *default:* `1`)

> Logging verbosity level.  
> This affects how much information is printed to the console during execution.  

> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG

<br>
