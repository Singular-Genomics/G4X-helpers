
# `--update_bin`
#### Update G4X-Viewer `.bin` file with new metadata

Updates an existing `.bin` segmentation file with new clustering, embedding, or color metadata from a CSV file.

--8<-- "_partials/section_break.md"

## usage 

```bash
$ update_bin 
  --bin_file /path/to/sample.bin
  --out_path /path/to/updated_sample.bin
  --metadata /path/to/metadata.csv

  # ─── optional ───
  --cellid_key <cellid_key>
  --cluster_key <cluster_key>
  --cluster_color_key <cluster_color_key>
  --emb_key <emb_key>
  --verbose <level>
```

## Argument Descriptions

### required
---
#### `--bin_file`: (*type:* `str`)
 
> Path to the existing `.bin` file used by the G4X-Viewer. This file will be updated with metadata from the provided CSV file.

#### `--out_path`: (*type:* `str`)

> Output path where the updated `.bin` file will be saved. The directory will be created if it doesn't already exist.

#### `--metadata`: (*type:* `str`)

> Path to a CSV file containing metadata to update the `.bin` file. This file must include a header row and should contain cell IDs, optional cluster assignments, colors, and embeddings.

### optional
---
#### `--cellid_key`: (*type:* `str`  *default:* `None`)

> Name of the column in the metadata file that contains cell IDs matching those in the `.bin` file. If not provided, the first column in the metadata will be used.

#### `--cluster_key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides cluster assignments for each cell. Required if `--cluster_color_key` is used.

#### `--cluster_color_key`: (*type:* `str`  *default:* `None`)

> Column name in the metadata that provides RGB or hex colors for each cluster. Must be used in conjunction with `--cluster_key`.

#### `--emb_key`: (*type:* `str`  *default:* `None`)

> Prefix for embedding coordinates. The script will look for two columns: `{emb_key}_0` and `{emb_key}_1`. Used to embed cells in UMAP/tSNE/etc. space.

#### `--verbose`: (*type:* `int`  *default:* `1`)

> Logging verbosity level:
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG
> This affects how much information is printed to the console during execution.

<br>
