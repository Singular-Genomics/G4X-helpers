<br>

# source installation

This page explains how you can install `G4X-helpers` on your local machine. If you’re happy running the tool via `Docker`, you can skip these installation steps and head to the [Docker setup](./docker.md) section.  

**You only need a local install if you want to:**

- use the `G4X-helpers` CLI directly on your machine
- use the package beyond the CLI features exposed in the Docker image
- develop or debug the codebase
- integrate pieces of the library into your own Python environment


If these use cases apply to you, please read on

<br>

## Step 1: clone the [`G4X-helpers`](https://github.com/Singular-Genomics/G4X-helpers) repository
---

```bash
git clone https://github.com/singular-genomics/g4x-helpers
```

navigate to the repo directory:
```bash
cd G4X-helpers
```

!!!note
    All subsequent steps assume that you are in the G4X-helpers directory. You can confirm this in your terminal via `pwd`.
    You should see a path ending in: `/G4X-helpers`.
    ```bash
    > pwd
    /path/to/current/directory/.../G4X-helpers
    ```

<br>

## Step 2: install the package
---

!!! note
    
    G4X-helpers depends on [`Glymur`](https://glymur.readthedocs.io/en/v0.14.2/) and [`OpenJPEG >= 2.2.0`](https://www.openjpeg.org/) for multi-threaded image loading. On some systems, Glymur may require additional configuration to [detect OpenJPEG](https://glymur.readthedocs.io/en/v0.14.2/detailed_installation.html).
    
    
    Due to those limitations **we strongly suggest installing OpenJPEG, and thus G4X-helpers, via `conda`**, as it reliably provides compatible and properly linked dependencies across platforms.
    Other installation methods (e.g., Homebrew or manual builds) may lead to issues such as Glymur reporting version `0.0.0` or failing to load JPEG 2000 (.jp2) images.  

<br>

=== "<span style="font-size:1rem">Conda (recommended)</span>"

    #### Create a `conda` environment

    Install **miniconda**, **conda**, or **mamba** ([instructions](https://www.anaconda.com/docs/getting-started/miniconda/install)).
    <br>

    If this is your first time using conda:

    ```bash
    conda init
    ```

    Create the environment:

    ```bash
    conda create -n g4x-helpers_env python=3.10
    ```

    Activate the environment:

    ```bash
    conda activate g4x-helpers_env
    ```

    Install the package:
    
    ```bash
    pip install .
    ```

=== "<span style="font-size:1rem">pip</span>"
    
    #### Install into your current python environment via `pip`

    !!!warning
        `pip` does not create or manage virtual environments, so installing through `pip install .` will require that your local Python version is compatible with the package dependencites (`Python >= 3.10`).
    
    ```bash
    pip install .
    ```

=== "<span style="font-size:1rem">uv</span>"
    #### Create a `venv` using [uv](https://docs.astral.sh/uv/)

    ```bash
    uv sync
    ```  
    activate the environment  

    ```bash
    source .venv/bin/activate
    ```

<br>

## Step 3: Verify OpenJPEG installation
---


After installation of `G4X-helpers`, confirm that Glymur recognizes OpenJPEG via:

```bash
python -c "import glymur; print(glymur.version.openjpeg_version)"
```


!!! success  "Success: Version >=2.2.0"
    The output shows a version string >2.2.0
    ```
    2.4.1
    ```

!!! Failure "Failure: Version <2.2.0"
    ```
    0.0.0
    ```

If your OpenJPEG is installed and was detected successfully, you can now proceed to [use](../usage/usage.md) G4X-helpers.

If Glymur does not detect a compatible OpenJPEG version, please follow the steps below. If this happens, we **strongly suggest** performing the G4X-helpers installation and the following steps in a `conda` environment.

If the above installation was done using `uv` and Glymur failed to detect OpenJPEG, you will likely need to do the next steps with Homebrew, since conda installations are not visible to uv virtual environments.

Hints on other systems are provided, but not supported! You can find further details in the [Glymur documentation](https://glymur.readthedocs.io/en/v0.14.2/detailed_installation.html) on advanced installation methods.

<br>

## Step 4: Install OpenJPEG

=== "<span style="font-size:1rem">conda (recommended)</span>"

    Install openjpeg through `conda`

    ```bash
    conda install -c conda-forge openjpeg
    ```

    Return to Step 3 to verify that Glymur detects OpenJPEG with a version >=2.2.0.

=== "<span style="font-size:1rem">Homebrew</span>"
        
    Install OpenJPEG and pkg-config
    
    ```bash
    brew install openjpeg pkg-config
    ```

    Create a Glymur config directory

    ```bash
    mkdir -p ~/.config/glymur
    ```
    
    Add the Homebrew installed OpenJPEG path to a glymurrc file in the config
    
    ```bash
    printf "[library]\nopenjp2 = /opt/homebrew/lib/libopenjp2.dylib\n" > ~/.config/glymur/glymurrc
    ```

    Return to Step 3 to verify that Glymur detects OpenJPEG with a version >=2.2.0.

<br>

## Step 5: Verify installation
---

To start using G4X-helpers, ensure that you have activated the environment in which the package was installed. 
If your installation of G4X-helpers was successful, you can call the following commands from any terminal and help statements should be printed:


+ `resegment --help`  
+ `update_bin --help`  
+ `new_bin --help`  
+ `tar_viewer --help`


!!!tip
    You can see the expected output of the `--help` statements in the [CLI usage](../usage/usage.md/#cli-usage) section.


--8<-- "_partials/end_cap.md"
