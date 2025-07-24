<br>

# Docker setup

If you don’t want to install `G4X-helpers` locally, you can run its CLI tools from the Docker image.

!!! info "Why use Docker?"
    - **No local installs**: skip creating a Python environment.  
    - **Reproducibility**: everyone uses the same environment.  
    - **Isolated**: nothing leaks into (or depends on) your system Python.  
    - **Great for HPC/servers**: just bind your data directories.

### prerequisites

- Docker installed and running.
- Access to the G4X-helpers image (e.g. `ghcr.io/singular-genomics/g4x-helpers:<tag>`).  
  In your command, replace `<tag>` with a real version or `latest`.




## 1. Install Docker

Docker is availabe for most plaforms. Please refer to the Docker installation guide for [Docker Engine](https://docs.docker.com/engine/) (Linux), or [Docker Desktop](https://docs.docker.com/desktop/) (MacOS/Windows). If you are not sure if Docker is already installed, you can simply call `docker --version` in your terminal.

⸻

## 2. Pull the G4X-helpers image

docker pull ghcr.io/singular-genomics/g4x-helpers:<tag>

Replace <tag> with a version or use latest.

⸻

## 3. Run a one-off command

Execute any tool or open a shell inside the container:

```bash
# Replace /path/on/host with your local data path
docker run --rm \
  -v /path/on/host:/data \
  -w /data \
  -it ghcr.io/singular-genomics/g4x-helpers:<tag> bash
```

This mounts your data into /data in the container and drops you into a shell.

⸻

## 4. Common Docker options
	•	--rm: Automatically remove the container when it exits.
	•	-v HOST:CONTAINER: Bind mount a host directory.
	•	-w PATH: Set the working directory inside the container.
	•	-e VAR=VALUE: Pass environment variables into the container.
	•	--entrypoint: Override the default entrypoint.

Example running a g4x-helpers tool directly:

docker run --rm \
  -v $(pwd):/data \
  -w /data \
  ghcr.io/singular-genomics/g4x-helpers:<tag> \
  g4x-process --input sample.fastq --output results/




⸻

## 6. Updating & cleanup
	•	Update image:

docker pull ghcr.io/singular-genomics/g4x-helpers:<tag>

	•	Remove unused images/containers:

docker system prune -a



⸻

## 7. Troubleshooting
	•	Permission denied (socket):
	•	On Linux, add your user to the docker group or use sudo.

sudo usermod -aG docker $USER


	•	Image pull fails:
	•	Check login to GitHub Container Registry:


gh auth login –with-token 
```
	•	Container can’t see files:
	•	Ensure correct host path in -v and that files have appropriate permissions.

--8<-- "_partials/end_cap.md"
