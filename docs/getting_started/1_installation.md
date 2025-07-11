# 🚀 Installation

!!! tip

    The installation script may be inspected before use:

    === "macOS and Linux"

        ```bash
        $ curl -LsSf https://astral.sh/uv/install.sh | less
        ```

    === "Windows"

        ```bash
        PS> powershell -c "irm https://astral.sh/uv/install.ps1 | more"
        ```

    Alternatively, the installer or binaries can be downloaded directly from [GitHub](#github-releases).


## 📦 Source Installation / CLI Usage

`g4x-helpers` can be installed and run directly as a Python package.

### Step 1: Prepare a Python environment

- Install **conda**, **miniconda**, or **mamba**
- Create the environment:

```bash
$ conda create -n g4x-helpers_env python=3.10
```

- Activate the environment:

```bash
$ conda activate g4x-helpers_env
```

### Step 2: Install openJPEG

- Install openJPEG.
```bash
apt-get install libopenjp2-7-dev
```
OR
```bash
conda install -c conda-forge openjpeg
```
- Verify openJPEG version is ≥ 2.2.0.
```bash
python -c "import glymur; print(glymur.version.openjpeg_version)"
```

### Step 3: Clone and install `g4x-helpers`

- Clone the repository:

```bash
git clone git@github.com:Singular-Genomics/G4X-helpers.git
```

- Navigate to the repo directory:

```bash
cd G4x-helpers
```

- Install the package:

```bash
pip install .
```

### Step 4: Verify installation

After installation, you can call the following commands from any terminal and help statements should be printed:

```bash
resegment --help
update_bin --help
new_bin --help
tar_viewer --help
```