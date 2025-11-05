<br>

# `update_bin`
#### Update G4X-viewer `.bin` file with new metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

## usage 
---

```bash
$ g4x-helpers update_bin 
  --g4x-data /path/to/G4X/output
  --metadata /path/to/metadata.csv

  # ─── optional ───
  --output <output_dir>
  --sample-id <sample_id>
  --cellid-key <cellid_key>
  --cluster-key <cluster_key>
  --cluster-color-key <cluster_color_key>
  --emb-key <emb_key>
```

## argument descriptions
---
### required 

#### `--g4x-data`: (*type:* `str`)

> Path to the base directory containing G4X sample output (including the `g4x_viewer` folder with the current `.bin` file).

#### `--metadata`: (*type:* `str`)

> Path to a CSV file containing metadata to update the `.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

### optional 

#### `--output`: (*type:* `str`  *default:* `None`)

> Directory where the updated viewer assets will be written. If omitted, the existing viewer folder inside `--g4x-data` is updated in place.

#### `--sample-id`: (*type:* `str`  *default:* `None`)

> Optional override for the sample identifier used when naming the output `.bin`.

#### `--cellid-key`: (*type:* `str`  *default:* `None`)

> Name of the column in the metadata file that contains cell IDs matching those in the `.bin` file. If not provided, the first column in the metadata will be used.

#### `--cluster-key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster-color-key` is used.

#### `--cluster-color-key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster-key`.

#### `--emb-key`: (*type:* `str`  *default:* `None`)

> Prefix for embedding coordinates. The command looks for two columns: `{emb_key}_1` and `{emb_key}_2`. Used to embed cells in UMAP/tSNE/etc. space.

!!! note
    Control threads and verbosity on the root command:  
    `g4x-helpers --threads 8 update_bin --g4x-data ...`

<br>
