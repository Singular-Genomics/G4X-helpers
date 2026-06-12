<br>

# `redemux`
#### Reprocess G4X-data with a new transcript manifest

Generates a new [`transcript_table.csv.gz`](https://docs.singulargenomics.com/g4x_data/output_files/rna/#transcript_tablecsvgz) by demultiplexing the raw feature data against a provided list of probe sequences and mapping each feature to its corresponding target gene. It then proceeds to regenerate single-cell outputs and initializes a new G4X-viewer zarr store. Does not regenerate metrics.  

---

## Usage
![`g4x-helpers redemux --help`](../img/redemux-help.svg)

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--manifest`
_type_ : <span class="acc-2-code">`file path`</span>  
_example_  : `path/to/transcript_panel.csv`

> Path to the new transcript manifest for demuxing.  
> Must contain a `probe` column with entries formatted as `<gene>-<sequence>-<primer>`. Optional `gene_name` or `read_num` columns are respected if present; otherwise they are derived from `probe`. Invalid probe names are ignored.

---

### `--batch-size`
_type_ : <span class="acc-2-code">`integer`</span>  
_default_  : `1.000.000`

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

---

--8<-- "_partials/arg_no_downstream.md"

--8<-- "_partials/arg_branch.md"

<br>
--8<-- "_core/_partials/end_cap.md"
