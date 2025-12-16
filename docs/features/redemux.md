<br>

# `redemux`
#### Reprocess G4X-output with a new transcript manifest

Replaces or updates the transcript manifest in a G4X run, reassigns transcripts via demultiplexing, and regenerates all downstream single-cell data and `.tar` viewer files.

---

## Usage
![`g4x-helpers redemux --help`](../img/redemux-help.svg)

--8<-- "_partials/global_options_note.md"

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--manifest`
_type_ : <span class="acc-2-code">`file path`</span>  
_example_  : `path/to/transcript_panel.csv`

> Path to the new transcript manifest for demuxing.  
> Must contain a `probe_name` column with entries formatted as `<gene>-<sequence>-<primer>`. Optional `gene_name` or `read` columns are respected if present; otherwise they are derived from `probe_name`. Invalid probe names are ignored.

---

### `--batch-size`
_type_ : <span class="acc-2-code">`integer`</span>  
_default_  : `1.000.000`

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

---

--8<-- "_partials/arg_in_place.md"

<br>
--8<-- "_core/_partials/end_cap.md"
