<br>

# `redemux`
#### Reprocess G4X-output with a new transcript manifest

Replaces or updates the transcript manifest in a G4X run, reassigns transcripts via demultiplexing, and regenerates all downstream single-cell data and `.tar` viewer files.

## Usage
![`g4x-helpers redemux --help`](../img/redemux-help.svg)

--8<-- "_core/_partials/global_options_note.md"

## argument descriptions
---
### required
--8<-- "_core/_partials/arg_g4x_data.md"

#### `--manifest`: (*type:* `str`)

> Path to the new transcript manifest for demuxing.  
> The manifest must be a **3-column CSV** with the following header:
> ```
> target,sequence,read
> ```

> - **`target`**: The gene or feature identifier to assign.  
> - **`sequence`**: The nucleotide sequence associated with the target.  
> - **`read`**: The read index (e.g. `read_1` or `2`) for which the sequence is valid.

### optional
#### `--batch-size`: (*type:* `int`  *default:* `1,000,000`)

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

--8<-- "_core/_partials/arg_output.md"

--8<-- "_core/_partials/arg_smp_id.md"

<br>
--8<-- "_core/_partials/end_cap.md"
