# 🚀 Installation


## Step 1: clone the [`g4x-helpers`](https://github.com/Singular-Genomics/G4X-helpers) repository


```bash
$ git clone git@github.com:Singular-Genomics/G4X-helpers.git
```

navigate to the repo directory:
```bash
$ cd G4X-helpers
```

!!!note
    all following steps assume that you are in the G4X-helpers directory. You can confirm this via:
    ```bash
    $ pwd
    ```
    you should see a path ending in: G4X-helpers
    ```bash
    $ /path/to/current/directory/.../G4X-helpers
    ```
    
---
<br>

## Step 2: install the package

!!! note
    
    G4X-helpers depends on [`Glymur`](https://glymur.readthedocs.io/en/v0.14.2/) and [`OpenJPEG >= 2.2.0`](https://www.openjpeg.org/) for multi-threaded image loading. On some systems, Glymur may require advanced configuration to [properly detect OpenJPEG](https://glymur.readthedocs.io/en/v0.14.2/detailed_installation.html).
    
    
    Due to those limitations **we strongly suggest installing OpenJPEG, and thus G4X-helpers, via `Conda`**, as it reliably provides compatible and properly linked dependencies across platforms.
    Other installation methods (e.g., Homebrew or manual builds) may lead to issues such as Glymur reporting version `0.0.0` or failing to load JPEG 2000 images.  

    learn how to [install and verify OpenJPEG](#step-3-verify-openjpeg-installation) for use with Glymur

<br>

=== "<span style="font-size:1rem">Conda (recommended)</span>"

    ### Create a `Conda` environment

    Install **miniconda**, **conda**, or **mamba**. ([intstructions](https://www.anaconda.com/docs/getting-started/miniconda/install))  
    <br>
    create the environment:

    ```bash
    $ conda create -n g4x-helpers_env python=3.10
    
    # if this is your first time using conda ...
    # $ conda init
    ```

    activate the environment:

    ```bash
    $ conda activate g4x-helpers_env
    ```

    ### Install the package:
    
    ```bash
    $ pip install .
    ```

=== "<span style="font-size:1rem">pip</span>"
    
    ### Install into your current python environment via `pip`
    
    ```bash
    $ pip install .
    ```

=== "<span style="font-size:1rem">uv</span>"
    ### Create a `venv` using [uv](https://docs.astral.sh/uv/)

    ```bash
    $ uv sync
    ```  
    activate the environment  

    ```bash
    $ source .venv/bin/activate
    ```

---


## Step 3: verify OpenJPEG installation

After installation of `G4X-helpers`, you can confirm that Glymur recognizes OpenJPEG via:

```bash
$ python -c "import glymur; print(glymur.version.openjpeg_version)"
```

!!! success 
    The output shows a correct version string
    ```
    2.4.1
    ```

An OpenJPEG version above `2.2.0` is detected. You can now proceed to [use](./usage.md) G4X-helpers

!!! warning 
    ```
    # output
    0.0.0
    ```

Glymur does not detect OpenJPEG and reports `0.0.0` or other error.  
In this case we **strongly suggest** performing the G4X-helpers installation and the following steps to enable OpenJPEG in a Conda enviroment.  
Hints on other systems are provided, but not supported! You can find further details in the [Glymur documentation](https://glymur.readthedocs.io/en/v0.14.2/detailed_installation.html) on advanced installation methods.


## Step 4: install OpenJPEG

=== "Conda (recommended)"

    Inside your `Conda` environment:

    ```bash
    conda install -c conda-forge openjpeg
    ```

=== "Other systems"
    
    This has been tested on MacOS with Homebrew:
    
    install OpenJPEG and pkg-config
    
    ```bash
    $ brew install openjpeg pkg-config
    ```

    create a Glymur config directory

    ```bash
    mkdir -p ~/.config/glymur
    ```
    
    add the Homebrew installed OpenJPEG path to a glymurrc file in the config
    
    ```bash
    printf "[library]\nopenjp2 = /opt/homebrew/lib/libopenjp2.dylib\n" > ~/.config/glymur/glymurrc
    ```

---
<br>
