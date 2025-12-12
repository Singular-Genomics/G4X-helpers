<br>

# `migrate`
#### Upgrade legacy G4X-output layouts

Brings older G4X outputs up to the latest schema used by G4X-viewer and G4X-helpers. The command renames or relocates legacy files and regenerates viewer artifacts when needed.

Creates a backup in `<G4X-DATA>/g4x_helpers/migration_backup` before making any modifications.

---

## Usage
![`g4x-helpers migrate --help`](../img/migrate-help.svg)

--8<-- "_partials/global_options_note.md"

## Arguments
---
### required
--8<-- "_partials/arg_g4x_data.md"

### optional
#### `--restore`: (*flag*  *default:* `False`)

> Restore files from an existing migration backup (if present) instead of performing a migration. The backup directory is removed after a successful restore.

<br>
--8<-- "_core/_partials/end_cap.md"
