<br>

# `migrate`
#### Upgrade legacy G4X-data layouts

Migrates legacy G4X-data into the current G4X-data schema used by G4X-helpers and G4X-viewer. The command copies and converts supported raw data and metadata into a new output directory, then optionally runs downstream processing to create single-cell output and a `g4x-viewer.zarr` compatible with G4X-viewer v4.

---

## Usage
![`g4x-helpers migrate --help`](../img/migrate-help.svg)

--8<-- "_partials/args_optns.md"

---

--8<-- "_partials/arg_g4x_data.md"

---

### `--out-dir` / `-o` [required]
_type_ : <span class="acc-2-code">`directory`</span>  
_example_ : `path/to/migrated_sample`

> Output directory for migrated data.
> The directory must already exist and must not already contain a `sample.g4x` file.

---

### `--check` / `-c`
_type_ : <span class="acc-2-code">`flag`</span>  
_default_ : `not set`

> Check whether the legacy input can be migrated.
> This reports migration status without copying or converting data.

---

### `--roi`
_type_ : <span class="acc-2-code">`four integers`</span>  
_example_ : `--roi 0 0 4096 4096`

> Restrict migration to a rectangular region of interest.
> Coordinates are provided as `x0 y0 x1 y1` in pixel coordinates.

---

--8<-- "_partials/arg_no_downstream.md"

---

## Examples

Check whether a legacy output can be migrated:

```bash
g4x-helpers migrate /path/to/legacy_sample \
    --out-dir /path/to/migrated_sample \
    --check
```

Migrate a full legacy output and generate downstream outputs:

```bash
g4x-helpers migrate /path/to/legacy_sample \
    --out-dir /path/to/migrated_sample
```

Migrate only raw data and metadata:

```bash
g4x-helpers migrate /path/to/legacy_sample \
    --out-dir /path/to/migrated_sample \
    --no-downstream
```

Migrate a region of interest:

```bash
g4x-helpers migrate /path/to/legacy_sample \
    --out-dir /path/to/migrated_sample \
    --roi 0 0 4096 4096
```

<br>
--8<-- "_core/_partials/end_cap.md"
