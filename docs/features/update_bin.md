<br>

# `update_bin`
#### Update G4X-viewer `.bin` file with new metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

---

## Usage
![`g4x-helpers update_bin --help`](../img/update_bin-help.svg)

--8<-- "_partials/global_options_note.md"

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--metadata`
_type_ : <span class="acc-2-code">`file path`</span>  
_example_  : `path/to/cell_metadata.csv`

> Path to a CSV file containing metadata to update the `segmentation.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

---

### `--cellid-key`
_type_ : <span class="acc-2-code">`string`</span>  
_default_  : `cell_id`

> Name of the column in the metadata file that contains cell IDs matching those in the `segmentation.bin` file.

---

### `--cluster-key`
_type_ : <span class="acc-2-code">`string`</span>  
_default_  : `None`

> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster-color-key` is used.

---

### `--cluster-color-key`
_type_ : <span class="acc-2-code">`string`</span>  
_default_  : `None`

> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster-key`.

---

### `--emb-key`
_type_ : <span class="acc-2-code">`string`</span>  
_default_  : `None`

> Prefix for embedding coordinates. The command looks for two columns: `{emb_key}_1` and `{emb_key}_2`. Used to embed cells in UMAP/tSNE/etc. space.

---

--8<-- "_partials/arg_in_place.md"

<br>
--8<-- "_core/_partials/end_cap.md"
