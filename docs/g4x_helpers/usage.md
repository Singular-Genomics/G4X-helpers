<br>

# usage

G4X-helpers offers several tools for common post-processing needs, each of which is described in detail in the [CLI features](./features/index.md) section. On this page we will demonstrate how to execute these tools either from the CLI or through a Docker image (depending on how you [installed G4X-helpers](./installation/index.md)).
!!!note
    G4X-helpers can also be imported as a Python module, in which case it can provide some convenience functions to access and interact with your data. Documentation for this will follow in a separate section.

<br>

## CLI Usage
---
### activating your environment
Before getting started, you need to ensure that the environment in which you installed G4X-helpers is activated. How this is done depends on your installation method. We recommend [installing it into a Conda environment](./installation/source.md#step-2-install-the-package). 

If you created an environment with the name `g4x-helpers_env`, you will activate it with the following command:

```bash
$ conda activate g4x-helpers_env
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


### basic command pattern
every command follows standard CLI syntax:
```
$ command_name --option_1 VALUE_OPTION_1 --option_2 VALUE_OPTION_2
```

To use these functions, you will need the contents of a [G4X output folder](../g4x_data/g4x_output.md) (which is referred to as `run_base`) and in some cases some externally generated files.


### example: using `--resegment` from the CLI
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

If your command was successful, you will start seeing the output log of the `resegment` process.

<br>

## Docker Usage
---
### pull the G4X-helpers image

If you haven't done so already (maybe you skipped the [Docker setup](./installation/docker.md) section) you will first need to pull the G4X-helper Docker image. This is done via `docker pull`, which downloads the image (and all of its layers) from GitHub’s Container Registry ( ghcr.io ) to your local Docker cache.

```bash
$ docker pull ghcr.io/singular-genomics/g4x-helpers:latest
```

If you want to run a specific version, please refer to the [Docker setup](./installation/docker.md) section for details.


### basic command pattern

When running via Docker, the pattern is similar to [invoking a CLI command](#basic-command-pattern), with the exception that all paths and the Docker image need to be specified.
Inside the container, always reference the mounted paths (e.g. `/data/...`), not the host paths ([see full example below](#example-using-resegment-in-docker)).


```bash
docker run --rm \
  -v /host/path/to/g4x_output:/data/g4x_output \
  -v /host/path/to/other_inputs:/data/inputs \
  -v /host/path/to/outputs:/data/outputs \
  ghcr.io/singular-genomics/g4x-helpers:latest \
  command_name --option_1 VALUE_OPTION_1 --option_2 VALUE_OPTION_2
```

+ `-v host:container` mounts your folder to let the container see your files  
+ `--rm` cleans up the container after it exits  
+ `ghcr.io/singular-genomics/g4x-helpers:latest` uses the latest version of G4X-helpers  


### example: using `--resegment` in Docker

Here we run the same `--resegment` command as in the [cli-example](#example-using-resegment-from-the-cli) above

```bash
docker run --rm \
  -v /path/to/g4x_output/directory:/data/g4x_output \
  -v /path/to/custom/segmentation:/data/inputs \
  -v /path/to/resegmentation/output:/data/outputs \
  ghcr.io/singular-genomics/g4x-helpers:latest \
  resegment \
    --run_base /data/g4x_output/sample_id \
    --segmentation_mask /data/inputs/seg_mask.npz \
    --sample_id sample_A01 \
    --out_dir /data/outputs/resegmentation
```

--8<-- "_core/_partials/end_cap.md"
