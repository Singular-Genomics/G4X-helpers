<br>

# Usage

G4X-helpers offers several tools for common post-processing needs, each of which is described in detail in the [CLI features](../features/index.md) section. On this page, we will demonstrate how to execute these tools either from the CLI or through a Docker image (depending on how you [installed G4X-helpers](../installation/index.md)).

## Quickstart
If you have installed G4X-helpers and your environment is functional, you can get an overview of the basic usage by calling the main entrypoint `g4x-helpers --help`.  
Please refer to the sections below for details on how to run a subcommand.

![`g4x-helpers --help`](../img/main-help.svg)

<br>

## CLI usage
---
### Activating your environment
Before getting started, you need to ensure that the environment in which you installed G4X-helpers is activated. How this is done depends on your installation method. We recommend [installing it into a Conda environment](../installation/source.md#step-2-install-the-package). 

=== "<span style="font-size:1rem">conda</span>"

    If you created an environment with the name `g4x-helpers_env`, you will activate it with the following command:

    ```bash
    conda activate g4x-helpers_env
    ```

    Once activated, your shell prompt should change to include the environment name:

    ```bash
    (g4x-helpers_env) $
    ```
    You can further verify that you’re using the correct interpreter:

    ```bash
    (g4x-helpers_env) $ which python
    /Users/you/miniconda3/envs/g4x-helpers_env/bin/python
    ```

=== "<span style="font-size:1rem">uv</span>"

    ```bash
    source .venv/bin/activate
    ```

    Once activated, your shell prompt should change to include the environment name (`g4x_helpers` is the default):
    
    ```bash
    (g4x-helpers) $
    ```

    You can further verify that you’re using the correct interpreter:

    ```bash
    (g4x-helpers) $ which python
    /Users/You/G4X-helpers/.venv/bin/python
    ```
<br>

### Basic command pattern
---
```
g4x-helpers [GLOBAL OPTIONS] <command> --option-1 VALUE --option-2 VALUE
```

The global options (`--threads`, `--verbose`, `--version`) live on the `g4x-helpers` command itself and should be provided **before** the sub-command. Each command then defines its own required and optional arguments. All workflows assume you have a single-sample [G4X output directory](https://docs.singulargenomics.com/g4x_data/g4x_output/) and, in some cases, additional input files (e.g. manifests or metadata tables).

#### Global options

| option | default | description |
| --- | --- | --- |
| `--threads` / `-t` | system-dependent | Number of worker threads used by CPU-heavy steps (segmentation, demultiplexing, conversions). Increase to speed up processing or lower it when sharing a server. |
| `--verbose` / `-v` | `1` | Controls logging detail in the terminal. `0` prints warnings only, `1` shows progress-level logs, `2` enables detailed debugging output. |
| `--version` | — | Prints the installed `g4x-helpers` version and exits, which is helpful for troubleshooting and documentation. |
| `--help` / `-h` | — | Standard Click help flag that prints contextual usage information for the command you attach it to. |

!!! tip
    Always place global options **before** the sub-command: `g4x-helpers --threads 8 resegment ...`

#### Discovering CLI commands

- List every available sub-command and a short description:

    ```bash
    g4x-helpers --help
    ```

- Show detailed usage and arguments for a specific command (e.g. `redemux`):

    ```bash
    g4x-helpers redemux --help
    ```

- Confirm which package version your environment is using:

    ```bash
    g4x-helpers --version
    ```


### Example: running `resegment` from the CLI

Below we apply a custom cell segmentation to existing data using the `resegment` sub-command.

First, inspect the available options:

```
$ g4x-helpers resegment --help

Usage: g4x-helpers [OPTIONS] resegment [OPTIONS]

  Reprocess G4X-output with a new segmentation

Options:
  -i, --g4x-data PATH      Directory containing G4X-data for a single sample  [required]
  -s, --sample-id TEXT     Sample ID (used for naming outputs)
  -o, --output PATH        Output directory used. Will edit G4X-data in-place if not provided.
  --cell-labels PATH       File containing cell segmentation labels. supported file types: [.npy, .npz, .geojson]
                           [required]
  --labels-key TEXT        Key/column in npz/geojson where labels should be taken from (optional)
  -h, --help               Show this message and exit.
```

For this run we provide:

| argument | value | type |
| --- | --- | --- |
| `--g4x-data` | /path/to/g4x_output/directory/sample_id | directory |
| `--cell-labels` | /path/to/custom/segmentation/seg_mask.npz | .npz file |
| `--sample-id` | sample_A01 | string |
| `--output` | /path/to/resegmentation/output | directory |

The full command is:

```
g4x-helpers resegment \
    --g4x-data /path/to/g4x_output/directory/sample_id \
    --cell-labels /path/to/custom/segmentation/seg_mask.npz \
    --sample-id sample_A01 \
    --output /path/to/resegmentation/output
```

If the command succeeds you'll see the `resegment` progress log in your terminal.

<br>

## Docker usage
---
### Pull the G4X-helpers image

If you haven't done so already (maybe you skipped the [Docker setup](../installation/docker.md) section) you will first need to pull the G4X-helper Docker image. This is done via `docker pull`, which downloads the image (and all of its layers) from GitHub’s Container Registry ( ghcr.io ) to your local Docker cache.

```bash
docker pull ghcr.io/singular-genomics/g4x-helpers:latest
```

If you want to run a specific version, please refer to the [Docker setup](../installation/docker.md) section for details.

### Basic command pattern

When running via Docker, the pattern is similar to [invoking a CLI command](#basic-command-pattern), with the exception that all paths and the Docker image need to be specified.
Inside the container, always reference the mounted paths (e.g. `/data/...`), not the host paths ([see full example below](#example-using-resegment-in-docker)).


```bash
docker run --rm \
  -v /host/path/to/g4x_output:/data/g4x_output \
  -v /host/path/to/other_inputs:/data/inputs \
  -v /host/path/to/outputs:/data/outputs \
  ghcr.io/singular-genomics/g4x-helpers:latest \
  g4x-helpers [GLOBAL OPTIONS] <command> --option-1 VALUE --option-2 VALUE
```

+ `-v host:container` mounts your folder to let the container see your files. 
+ `--rm` cleans up the container after it exits.
+ `ghcr.io/singular-genomics/g4x-helpers:latest` uses the latest version of G4X-helpers


### Example: using `resegment` in Docker

Here we run the same `resegment` command as in the [CLI example](#example-running-resegment-from-the-cli) above.

```bash
docker run --rm \
  -v /path/to/g4x_output/directory:/data/g4x_output \
  -v /path/to/custom/segmentation:/data/inputs \
  -v /path/to/resegmentation/output:/data/outputs \
  ghcr.io/your-org/g4x-helpers:latest \
  g4x-helpers resegment \
    --g4x-data /data/g4x_output/sample_id \
    --cell-labels /data/inputs/seg_mask.npz \
    --sample-id sample_A01 \
    --output /data/outputs/resegmentation
```

## Python module
---
!!!note
    G4X-helpers can also be imported as a Python module, in which case it can provide some convenience functions to access and interact with your data.  
    Documentation for this will follow in a separate section.

<br>


--8<-- "_core/_partials/end_cap.md"
