# <span class="index-cat-header">usage</span>


<!-- begin section -->
## verify installation

To start using G4X-helpers, ensure that you have activated the environment in which the package was installed. 
If your installation of G4X-helpers was successful, you can call the following commands from any terminal and help statements should be printed:


+ `resegment --help`  
+ `update_bin --help`  
+ `new_bin --help`  
+ `tar_viewer --help`


!!!tip
    [further down](#example-using-resegment-from-the-cli) you can see the output of such a `--help` statement

---
<br>
<!-- end section -->

<!-- begin section -->
## CLI Usage

G4X-helpers offers several tools for common post-processing needs, each of which is described in detail in the [CLI features](./features/index.md) section.

### Basic pattern
every command follows standard CLI syntax:
```
$ command_name --option_1 VALUE_OPTION_1 --option_2 VALUE_OPTION_2
```

To use these functions, you will need the contents of a [G4X output folder](../g4x_data/g4x_output.md), which is referred to as `run_base`, and in some cases some externally generated files.


### Example: using `--resegment` from the CLI
In the following example we want to apply a custom cell segmentation to our data, which is done via the `--resegment` tool.

In the [feature reference](./features/resegment.md) for this command, we can see that there are required and optional arguments.
We can also call `resegment --help` to see which arguments the tool accepts:

```
usage: resegment [-h] --run_base RUN_BASE --segmentation_mask SEGMENTATION_MASK [--sample_id SAMPLE_ID] [--out_dir OUT_DIR]
                [--segmentation_mask_key SEGMENTATION_MASK_KEY] [--threads THREADS] [--verbose VERBOSE]

options:
-h, --help            show this help message and exit
--run_base RUN_BASE   Path to G4X sample output folder
    ...

```

In our scenario, we will provide the following values:

| argument | value | type |
| --- | --- | --- |
| `--run_base` | /path/to/g4x_output/directory/sample_id | directory |
| `--segmentation_mask` | /path/to/custom/segmentation/seg_mask.npz | .npz file |
| `--sample_id` | sample_A01 | string |
| `--out_dir` | /path/to/resegmentation/output | directory |

The full command to run this operation is:

```
resegment \
    --run_base /path/to/g4x_output/directory/sample_id \
    --segmentation_mask /path/to/custom/segmentation/seg_mask.npz \
    --sample_id sample_A01 \
    --out_dir /path/to/resegmentation/output 
```
<!-- end section -->
<br>

---

<!-- begin section -->
## Docker Usage

If you don’t want to install `g4x-helpers` locally, you can run all tools from the Docker image.

!!! info "Why use Docker?"
    - **No local installs**: skip creating a Python environment.  
    - **Reproducibility**: everyone uses the same environment.  
    - **Isolated**: nothing leaks into (or depends on) your system Python.  
    - **Great for HPC/servers**: just bind your data directories.

### Prerequisites

- Docker installed and running.
- Access to the G4X-helpers image (e.g. `ghcr.io/singular-genomics/g4x-helpers:<tag>`).  
  In your command, replace `<tag>` with a real version or `latest`.

---

### Basic pattern

When running via Docker, the pattern is similar to invoking a CLI command, with the exception that all paths and the Docker image need to be specified.
Inside the container, always reference the mounted paths (e.g. /data/...), not the host paths ([see full example below](#example-using-resegment-in-docker)).


```bash
docker run --rm \
  -v /host/path/to/g4x_output:/data/g4x_output \
  -v /host/path/to/other_inputs:/data/inputs \
  -v /host/path/to/outputs:/data/outputs \
  ghcr.io/singular-genomics/g4x-helpers:latest \
  command_name --option_1 VALUE_OPTION_1 --option_2 VALUE_OPTION_2
```

+ `-v host:container` mounts your folder to let the container see your files. 
+ `--rm` cleans up the container after it exits.+
+ `ghcr.io/singular-genomics/g4x-helpers:latest` uses the latest version of G4X-helpers


### Example: using `--resegment` in Docker

Here we run the same --resegment command as in the [cli-example](#example-using-resegment-from-the-cli) above

```bash
docker run --rm \
  -v /path/to/g4x_output/directory:/data/g4x_output \
  -v /path/to/custom/segmentation:/data/inputs \
  -v /path/to/resegmentation/output:/data/outputs \
  ghcr.io/your-org/g4x-helpers:latest \
  resegment \
    --run_base /data/g4x_output/sample_id \
    --segmentation_mask /data/inputs/seg_mask.npz \
    --sample_id sample_A01 \
    --out_dir /data/outputs/resegmentation
```
<!-- end section -->
<br>

---

<!-- begin section -->
