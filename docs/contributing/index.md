[G4X-helpers]: https://github.com/Singular-Genomics/G4X-helpers
[mkdocs]: https://squidfunk.github.io/mkdocs-material/
[pre-commit]: https://pre-commit.com
[GitHub's CLI]: https://cli.github.com
[scanpy docs]: https://scanpy.readthedocs.io/en/stable/dev/code.html#contributing-code
[pandas]: https://pandas.pydata.org/pandas-docs/stable/development/index.html
[MDAnalysis]: https://userguide.mdanalysis.org/stable/contributing.html

# <span class="g4x-header">Contributing</span>

Contributions to `G4X-helpers` are welcome!  
This section provides some guidelines and tips to follow when you wish to enrich the project with your own code.


## Index
+ [Development workflow](#development-workflow)  
+ [Working with git](#working-with-git)
+ [Building and managing your environment](#building-and-managing-your-environment)
+ [Documentation](#documentation)
+ [Releasing a new version](#documentation)


!!!info
    Parts of these guidelines have been adapted from the [scanpy docs][], which in turn built on the work done by [pandas][] and [MDAnalysis][].
    We highly recommend checking out these excellent guides to learn more.

---
<br>

# Development workflow

1. [Fork](#forking-and-cloning) the `G4X-helpers` repository to your own GitHub account
2. Create an [environment](#building-and-managing-your-environment) with all dev-dependencies
3. Create a [new branch](#creating-a-branch-for-your-feature) for your feature or bugfix
4. Commit your contribution to the codebase
5. Update and check the [documentation](#documentation)
6. [Open a PR](#opening-a-pull-request) back to the main repository

---
<br>



# Working with `git`

This section of the docs covers our practices for working with `git` on our codebase.  
For more in-depth guides, we can recommend a few sources:

[Atlassian's git tutorial](https://www.atlassian.com/git/tutorials)
: Beginner friendly introductions to the git command line interface

[Setting up git for GitHub](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/set-up-git)
: Configuring git to work with your GitHub user account


## Forking and cloning

To get the code, and be able to push changes back to the main project, you'll need to (1) fork the repository on github and (2) clone the repository to your local machine.

This is very straight forward if you're using [GitHub's CLI][]:

```console
$ gh repo fork Singular-Genomics/G4X-helpers --clone --remote
```

This will fork the repo to your github account, create a clone of the repo on your current machine, add our repository as a remote, and set the `main` development branch to track our repository.

To do this manually, first make a fork of the repository by clicking the "fork" button on our main github package. Then, on your machine, run:

```bash
$ # Clone your fork of the repository (substitute in your username)
$ git clone https://github.com/{your-username}/G4X-helpers.git

$ # Enter the cloned repository
$ cd G4X-helpers

$ # Add our repository as a remote
$ git remote add upstream https://github.com/Singular-Genomics/G4X-helpers.git
$ # git branch --set-upstream-to "upstream/main"
```

## Creating a branch for your feature

All development should occur in branches dedicated to the particular work being done.
Additionally, unless you are a maintainer, all changes should be directed at the `main` branch.
You can create a branch with:

```bash
$ git checkout main                 # Starting from the main branch
$ git pull                          # Syncing with the repo
$ git switch -c {your-branch-name}  # Making and changing to the new branch
```


## Opening a pull request

When you're ready to have your code reviewed, push your changes up to your fork:

```bash
$ # The first time you push the branch, you'll need to tell git where
$ git push --set-upstream origin {your-branch-name}

$ # After that, just use
$ git push
```

And open a pull request by going to the [main repo](https://github.com/Singular-Genomics/G4X-helpers) and clicking *`New pull request`*.  
GitHub may also prompt you to open PRs for recently pushed branches.

We'll try and get back to you soon!

---
<br>



# Building and managing your environment

## Installing project dependencies

It is recommended to develop your feature in an isolated virtual environment.
There are many environment managers available for Python (conda, pyenv, Virtualenv ...)

We recommend using [uv](https://docs.astral.sh/uv/), which can manage your virtual environment and use the project's `uv.lock` file to replicate all dependencies from exact sources.

After installing uv, you can build the environment by calling:
```bash
$ uv sync
```

A folder named `.venv` will be created. It holds the correct python version and all project dependencies. It will also install necessary development tools like `ruff`, `mkdocs`, `pre-commit`, `bump-my-version`.

You can now activate this environment with:
```bash
$ source .venv/bin/activate
```

<br>

## Using pre-commit hooks

We use [pre-commit][] to run various checks on new code.  


In order for it to attach to your commits automatically, you need to install it once after building your environment. 

```bash
# from the root of the repo.
$ pre-commit install
```

While most rules will be applied automatically, some checks may prevent your code from being commited. The pre-commit output will help you indentify which sections need to be addressed.

If you choose not to run the hooks on each commit, you can run them manually with  
`pre-commit run --files={your files}`.

!!!note
    If your environment manager did not install pre-commit as a dependency, you can do so via

    ```bash
    $ pip install pre-commit
    ```

<br>

## Code formatting and linting
We use [Ruff](https://docs.astral.sh/ruff) to format and lint the `G4X-helpers` codebase. Ruff is a project dependency and its rules are configured in `ruff.toml`. It  will be invoked on all code contributions via pre-commit hooks (see above) but you can also run it manually via `ruff check`.


<br>


# Documentation

## docstrings

We prefer the numpydoc style for writing docstrings.
We'd primarily suggest looking at existing docstrings for examples, but the [napolean guide to numpy style docstrings][] is also a great source.
If you're unfamiliar with the reStructuredText (rST) markup format, check out the [Sphinx rST primer][].

Look at [`sc.tl.leiden`](https://github.com/scverse/scanpy/blob/350c3424d2f96c4a3a7bb3b7d0428d38d842ebe8/src/scanpy/tools/_leiden.py#L49-L120) as an example of a complete doctring.

[napolean guide to numpy style docstrings]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy
[sphinx rst primer]: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html


### `Params` section

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
These will automatically populate the docstrings on import, and when the documentation is built.

Use the python standard library types (defined in {mod}`collections.abc` and {mod}`typing` modules) for containers, e.g.
{class}`~collections.abc.Sequence`s (like `list`),
{class}`~collections.abc.Iterable`s (like `set`), and
{class}`~collections.abc.Mapping`s (like `dict`).
Always specify what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so: `Literal['elem-1', 'elem-2']`.

### `Returns` section

There are three types of return sections – prose, tuple, and a mix of both.

1. Prose is for simple cases.
2. Tuple return sections are formatted like parameters. Other than in numpydoc, each tuple is first characterized by the identifier and *not* by its type. Provide type annotation in the function header.
3. Mix of prose and tuple is relevant in complicated cases, e.g. when you want to describe that you *added something as annotation to an \`AnnData\` object*.

#### Examples

For simple cases, use prose as in [`sc.pp.normalize_total`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html#scanpy.pp.normalize_total):

```rst
Returns
-------
Returns dictionary with normalized copies of `adata.X` and `adata.layers`
or updates `adata` with normalized versions of the original
`adata.X` and `adata.layers`, depending on `inplace`.
```

For tuple return values, you can use the standard numpydoc way of populating it,
e.g. as in [`sc.pp.calculate_qc_metrics`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics).
Do not add types in the docstring, but specify them in the function signature:

```python
def myfunc(...) -> tuple[int, str]:
    """
    ...
    Returns
    -------
    one_identifier
        Description.
    second_identifier
        Description 2.
    """
    ...
```

Many functions also just modify parts of the passed AnnData object, like e.g. [`scanpy.tl.dpt`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.dpt.html#scanpy.tl.dpt).
You can then combine prose and lists to best describe what happens:

```rst
Returns
-------
Depending on `copy`, returns or updates `adata` with the following fields.

If `n_branchings==0`, no field `dpt_groups` will be written.

dpt_pseudotime : :class:`~pandas.Series` (`adata.obs`, dtype `float`)
    Array of dim (number of samples) that stores the pseudotime of each
    cell, that is, the DPT distance with respect to the root cell.
dpt_groups : :class:`pandas.Series` (`adata.obs`, dtype `category`)
    Array of dim (number of samples) that stores the subgroup id ('0',
    '1', ...) for each cell. The groups  typically correspond to
    'progenitor cells', 'undecided cells' or 'branches' of a process.
```
