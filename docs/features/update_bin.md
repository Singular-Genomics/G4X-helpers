<br>

# `update_bin`
#### Update G4X-viewer `.bin` file with new metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

## Usage
![`g4x-helpers update_bin --help`](../img/update_bin-help.svg)

--8<-- "_partials/global_options_note.md"

## argument descriptions
---
### required 
--8<-- "_partials/arg_g4x_data.md"


#### `--metadata`: (*type:* `str`)

> Path to a CSV file containing metadata to update the `.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

### optional 

--8<-- "_partials/arg_in_place.md"

#### `--cellid-key`: (*type:* `str`  *default:* `cell_id`)

> Name of the column in the metadata file that contains cell IDs matching those in the `.bin` file. Defaults to `cell_id`.

#### `--cluster-key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster-color-key` is used.

#### `--cluster-color-key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster-key`.

#### `--emb-key`: (*type:* `str`  *default:* `None`)

> Prefix for embedding coordinates. The command looks for two columns: `{emb_key}_1` and `{emb_key}_2`. Used to embed cells in UMAP/tSNE/etc. space.

<br>
--8<-- "_core/_partials/end_cap.md"
