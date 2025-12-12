<br>

# `redemux`
#### Reprocess G4X-output with a new transcript manifest

Replaces or updates the transcript manifest in a G4X run, reassigns transcripts via demultiplexing, and regenerates all downstream single-cell data and `.tar` viewer files.

---

## Usage
![`g4x-helpers redemux --help`](../img/redemux-help.svg)

--8<-- "_partials/global_options_note.md"

## Arguments
---
### required
--8<-- "_partials/arg_g4x_data.md"

#### `--manifest`: (*type:* `str`)

> Path to the new transcript manifest for demuxing.  
> Must contain a `probe_name` column with entries formatted as `<gene>-<sequence>-<primer>`. Optional `gene_name` or `read` columns are respected if present; otherwise they are derived from `probe_name`. Invalid probe names are ignored.

### optional
#### `--batch-size`: (*type:* `int`  *default:* `1,000,000`)

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

--8<-- "_partials/arg_in_place.md"

<br>
--8<-- "_core/_partials/end_cap.md"
