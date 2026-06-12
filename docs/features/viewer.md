<br>

# `viewer`
#### Modify metadata in a `g4x-viewer.zarr` store

Exports or imports metadata tables for an existing `g4x-viewer.zarr`.  
Each viewer data layer has its own subcommand:

- `images` for multiplex image channel metadata
- `cells` for cell labels, label colors, and UMAP coordinates
- `transcripts` for gene colors and gene display order

The recommended workflow is to export the current metadata, edit the exported CSV, and import the modified file back into the same `g4x-viewer.zarr`.

---

## Usage
![`g4x-helpers viewer --help`](../img/viewer-help.svg)

--8<-- "_partials/args_optns.md"

---

### `VIEWER-ZARR`
_type_ : <span class="acc-2-code">`zarr store`</span>  
_example_ : `g4x-data/g4x-viewer.zarr`

> Positional path to an existing `g4x-viewer.zarr` directory.
> The path is provided before the `images`, `cells`, or `transcripts` subcommand.

---

### `--export-metadata` / `-e`
_type_ : <span class="acc-2-code">`file path`</span>  
_example_ : `image_metadata.csv`

> Export the selected layer metadata to a CSV file.
> Use this file as the safest starting point for edits, because import checks expect the same layer-specific rows and columns.

---

### `--import-metadata` / `-i`
_type_ : <span class="acc-2-code">`file path`</span>  
_example_ : `image_metadata_edited.csv`

> Import a CSV file and update the selected layer metadata in `g4x-viewer.zarr`.

!!! warning
    `--import-metadata` and `--export-metadata` cannot be used together in the same command.

---

## `images`
#### Modify multiplex image channel metadata

Use `images` to export or update channel labels, visibility, colors, and intensity windows for the multiplex image layer.

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr images \
    --export-metadata image_metadata.csv
```

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr images \
    --import-metadata image_metadata.csv
```

The metadata CSV must contain these columns:

| column | description |
| --- | --- |
| `label` | Channel label shown in the viewer. |
| `active` | Whether the channel is enabled by default. |
| `color` | Channel color, using the exported color format. |
| `min` | Minimum value for the channel intensity range. |
| `max` | Maximum value for the channel intensity range. |
| `start` | Lower value for the initially selected channel window. |
| `end` | Upper value for the initially selected channel window. |

Imported image metadata must contain the same labels as the existing image layer.

---

## `cells`
#### Modify cell labels, colors, and UMAP coordinates

Use `cells` to export or update cell-level metadata for a segmentation group.

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr cells \
    --export-metadata cell_metadata.csv
```

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr cells \
    --import-metadata cell_metadata.csv
```

The metadata CSV must include:

| column | description |
| --- | --- |
| `cell_id` | Cell identifier. Values and order must match the existing viewer data. |
| `UMAP1` | First UMAP coordinate. |
| `UMAP2` | Second UMAP coordinate. |
| `<label>` | One or more categorical cell-label columns. |
| `<label>_color` | Optional color column paired with a label column. |

If multiple segmentation groups are present, select one with `--segmentation`:

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr cells \
    --segmentation g4x_default_segmentation \
    --export-metadata cell_metadata.csv
```

---

## `transcripts`
#### Modify transcript gene colors

Use `transcripts` to export or update gene colors for the transcript layer.

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr transcripts \
    --export-metadata transcript_metadata.csv
```

```bash
g4x-helpers viewer path/to/g4x-viewer.zarr transcripts \
    --import-metadata transcript_metadata.csv
```

The metadata CSV must contain these columns:

| column | description |
| --- | --- |
| `gene_id` | Gene identifier. The imported file must contain the same genes as the existing transcript layer. |
| `color` | Gene color as a hex color string. |

The row order in the imported file becomes the viewer gene order.

<br>
--8<-- "_core/_partials/end_cap.md"
