<br>

# Docker setup

If you don’t want to install `g4x-helpers` locally, you can run all tools from the Docker image.

!!! info "Why use Docker?"
    - **No local installs**: skip creating a Python environment.  
    - **Reproducibility**: everyone uses the same environment.  
    - **Isolated**: nothing leaks into (or depends on) your system Python.  
    - **Great for HPC/servers**: just bind your data directories.

### prerequisites

- Docker installed and running.
- Access to the G4X-helpers image (e.g. `ghcr.io/singular-genomics/g4x-helpers:<tag>`).  
  In your command, replace `<tag>` with a real version or `latest`.
