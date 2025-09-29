<br>

# Docker setup

If you prefer working with Docker or you don’t want to install G4X-helpers locally, you can run its CLI tools from the [Docker image](https://github.com/Singular-Genomics/G4X-helpers/pkgs/container/g4x-helpers) that is published along with the repository.

If you have Docker already installed (or Apptainer, Podman ... ) and are familiar with using it, you can continue to the G4X-helpers [Docker usage](../usage/index.md#docker-usage) section.  

Otherwise, please read on to get started.

!!! info "Why use Docker?"
    - **no local installs**: skip creating a Python environment.  
    - **reproducibility**: everyone uses the same environment.  
    - **isolated**: nothing leaks into (or depends on) your system Python.  
    - **great for HPC/servers**: just bind your data directories.

<br>

## Step 1: Install Docker
---
Docker is available for most platforms. Please refer to the Docker installation guide for [Docker Engine](https://docs.docker.com/engine/) (Linux), or [Docker Desktop](https://docs.docker.com/desktop/) (MacOS/Windows). If you are not sure if Docker is already installed, you can simply call `docker --version` in your terminal.

<br>

## Step 2: Pull the G4X-helpers image
---

The docker pull command downloads the G4X-helpers container image (and all of its layers) from GitHub’s Container Registry ( ghcr.io ) to your local Docker cache.

```bash
docker pull ghcr.io/singular-genomics/g4x-helpers:latest
```

Once the image is pulled, you can start new containers from it instantly — even when you’re offline — without re-fetching data from the registry.

!!!tip
    `latest` will always retrieve the most recently published build.  
    You can replace it with an explicit tag, which locks you to a reproducible version.  
    e.g. `ghcr.io/singular-genomics/g4x-helpers:v0.4.11`  
    
if the pull was successful, you can confirm that the image is available by calling

```bash
docker image ls ghcr.io/singular-genomics/g4x-helpers
```

You are now ready to use the G4X-helpers Docker image. Please refer to the [Docker usage](../usage/index.md#docker-usage) section to learn how to run the tools.

<br>

## Updating & cleaning up docker images
---

You can update your image by running `docker pull` again

```bash
docker pull ghcr.io/singular-genomics/g4x-helpers:latest
```
It will only download the layers that have changed since your last pull, giving you the newest build.

Over time, old images, stopped containers and dangling volumes can consume disk space.  

A useful one-shot cleanup command is:

```bash
docker system prune -a
```

For more general information on how to utilize the features of CLI Docker, please refer to the [Docker documentation](https://docs.docker.com/engine/reference/commandline/cli/).

--8<-- "_partials/end_cap.md"
