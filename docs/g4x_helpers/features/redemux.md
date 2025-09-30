# `--redemux`
#### Reprocess G4X output with a new transcript manifest

Replaces or updates the transcript manifest in a G4X run, reassigns transcripts via demultiplexing, and regenerates all downstream single-cell data and `.tar` viewer files.

## usage
---
```bash
$ redemux
  --run_base /path/to/G4X/output
  --manifest /path/to/manifest.csv

  # ─── optional ───
  --batch_size <n_transcripts>
  --out_dir <output_dir>
  --threads <n_threads>
  --verbose <level>
```

## argument descriptions
---
### required
#### `--run_base`: (*type:* `str`)

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
#### `--batch_size`: (*type:* `int`  *default:* `1,000,000`)

> Number of transcripts to process per batch during demultiplexing.  
> Larger batch sizes may improve performance but increase memory usage.

#### `--out_dir`: (*type:* `str`  *default:* `None`)

> Output directory where the redemuxed files will be written.  
> - If not provided, files in the `run_base` directory will be updated **in place**.  
> - If provided, the directory will be created (if it does not exist) and symlinked to the original run files (excluding specific diagnostic files).

#### `--threads`: (*type:* `int`  *default:* `4`)

> Number of threads to use for processing, including transcript reassignment, segmentation intersection, and viewer file generation.

#### `--verbose`: (*type:* `int`  *default:* `1`)

> Logging verbosity level.  
> - `0` = WARNING  
> - `1` = INFO  
> - `2` = DEBUG  

<br>
