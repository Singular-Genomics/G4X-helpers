<br>

# `redemux`
#### Reprocess G4X-output with a new transcript manifest

Replaces or updates the transcript manifest in a G4X run, reassigns transcripts via demultiplexing, and regenerates all downstream single-cell data and `.tar` viewer files.

## usage
---
```bash
$ g4x-helpers redemux
  --g4x-data /path/to/G4X/output
  --manifest /path/to/manifest.csv

  # ─── optional ───
  --batch-size <n_transcripts>
  --output <output_dir>
  --sample-id <sample_id>
```

## argument descriptions
---
### required
#### `--g4x-data`: (*type:* `str`)

> Path to the G4X sample output folder (the base directory for the run).  
> This directory must contain the required files such as `run_meta.json`, transcript panel, and feature tables.

#### `--manifest`: (*type:* `str`)

> Path to the new transcript manifest for demuxing.  
> The manifest must be a **3-column CSV** with the following header:
> ```
> target,sequence,read
> ```

> - **`target`**: The gene or feature identifier to assign.  
> - **`sequence`**: The nucleotide sequence associated with the target.  
> - **`read`**: The read index (e.g. `read_1` or `2`) for which the sequence is valid.

### optional
#### `--batch-size`: (*type:* `int`  *default:* `1,000,000`)

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

#### `--output`: (*type:* `str`  *default:* `None`)

> Output directory where the re-demuxed files will be written.  
> - If not provided, files in the `--g4x-data` directory will be updated **in place**.  
> - If provided, the directory will be created (if it does not exist) and symlinked to the original run files (excluding specific diagnostic files).

#### `--sample-id`: (*type:* `str`  *default:* `None`)

> Optional sample identifier override used for naming downstream outputs.

!!! note
    Supply runtime controls such as `--threads` and `--verbose` on the root command:  
    `g4x-helpers --threads 8 redemux --g4x-data ...`

<br>
