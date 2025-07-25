[G4X-helpers]: https://github.com/Singular-Genomics/G4X-helpers
[mkdocs]: https://squidfunk.github.io/mkdocs-material/
[pre-commit]: https://pre-commit.com
[GitHub's CLI]: https://cli.github.com
[scanpy docs]: https://scanpy.readthedocs.io/en/stable/dev/code.html#contributing-code
[pandas]: https://pandas.pydata.org/pandas-docs/stable/development/index.html
[MDAnalysis]: https://userguide.mdanalysis.org/stable/contributing.html


<!-- excluding this header from toc -->
<h3>Contributions to G4X-helpers are welcome!</h3>  

This section provides some guidelines and tips to follow when you wish to enrich the project with your own code or documentation.

## Index
---
+ [Development workflow](#development-workflow)  
+ [Working with git](#working-with-git)
+ [Building and managing your environment](#building-and-managing-your-environment)
+ [Documentation](#documentation)
+ [Releasing a new version](#releasing-and-versioning)

<br>

!!!info
    Parts of these guidelines have been adapted from the [scanpy docs][], which in turn built on the work done by [pandas][] and [MDAnalysis][].  
    We highly recommend checking out these excellent guides to learn more.

<br>

## Development workflow
---

The life-cycle of a new feature or other contribution should follow this pattern:

1. [Fork](#forking-and-cloning) the `G4X-helpers` repository to your own GitHub account
2. Create an [environment](#building-and-managing-your-environment) with all dev-dependencies
3. Create a [new branch](#creating-a-branch-for-your-feature) for your feature or bugfix
4. [Commit](#committing-your-code) your contribution to the codebase
5. Update and check the [documentation](#documentation)
6. [Open a PR](#opening-a-pull-request) back to the main repository

<br>

## Working with `git`
---

This section of the docs covers our practices for working with `git` on our codebase.  
For more in-depth guides, we can recommend a few sources:

[Atlassian's git tutorial](https://www.atlassian.com/git/tutorials)
: Beginner friendly introductions to the git command line interface

[Setting up git for GitHub](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/set-up-git)
: Configuring git to work with your GitHub user account


### Forking and cloning

To get the code, and be able to push changes back to the main project, you'll need to (1) fork the repository on github and (2) clone the repository to your local machine.

This is very straight forward if you're using [GitHub's CLI][]:

```console
$ gh repo fork Singular-Genomics/G4X-helpers --clone --remote
```

This will fork the repo to your github account, create a clone of the repo on your current machine, add our repository as a remote, and set the `main` development branch to track our repository.

To do this manually, first make a fork of the repository by clicking the "fork" button on our main github package. Then, on your machine, run:

```bash
# Clone your fork of the repository (substitute in your username)
$ git clone https://github.com/{your-username}/G4X-helpers.git

# Enter the cloned repository
$ cd G4X-helpers

# Add our repository as a remote
$ git remote add upstream https://github.com/Singular-Genomics/G4X-helpers.git
```

### Creating a branch for your feature

All development should occur in branches dedicated to the particular work being done.
Additionally, unless you are a maintainer, all changes should be directed at the `main` branch.
You can create a branch with:

```bash
$ git checkout main                 # Starting from the main branch
$ git pull                          # Syncing with the repo
$ git switch -c {your-branch-name}  # Making and changing to the new branch
```

### Committing your code

Keep commits small, focused, and well-described. This makes code review easier and history clearer.
When you are ready, add the files that belong to your commit:

```bash
$ git status              # See what changed
$ git add -p              # Interactively stage only the hunks you want

# or everything
$ git add .
```

Write a clear commit-message. Use an imperative, one-line summary (≤ 72 chars).  

```bash
$ git commit -m "Fix NPE in sample parser when header is missing"
```

!!!tip
    Need to fix the last commit before pushing?
    `git commit --amend` lets you change the message or add more files.

### Opening a pull request

When you're ready to have your code reviewed, push your changes up to your fork:

```bash
# The first time you push the branch, you'll need to tell git where
$ git push --set-upstream origin {your-branch-name}

# After that, just use
$ git push
```

Then open a pull request by going to the [main repo](https://github.com/Singular-Genomics/G4X-helpers) and clicking *`New pull request`*.  
GitHub may also prompt you to open PRs for recently pushed branches.

!!!info
    It is important to summarize your changes in the description of the PR so that they get included in the next change-log

*We'll try and get back to you soon!*

<br>

## Building and managing your environment
---

### Installing project dependencies

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

### Using pre-commit hooks

We use [pre-commit][] to run various checks on new code.  


In order for it to attach to your commits automatically, you need to install it once after building your environment. 

```bash
# from the root of the repo.
$ pre-commit install
```

While most rules will be applied automatically, some checks may prevent your code from being committed. The pre-commit output will help you identify which sections need to be addressed.

If you choose not to run the hooks on each commit, you can run them manually with  
`pre-commit run --files={your files}`.

!!!note
    If your environment manager did not install pre-commit as a dependency, you can do so via:

    ```bash
    $ pip install pre-commit
    ```

<br>

### Code formatting and linting
We use [Ruff](https://docs.astral.sh/ruff) to format and lint the `G4X-helpers` codebase. Ruff is a project dependency and its rules are configured in `ruff.toml`. It  will be invoked on all code contributions via pre-commit hooks (see above) but you can also run it manually via `ruff check`.

<br>

## Documentation
---

### docstrings

We prefer the numpydoc style for writing docstrings.
We'd primarily suggest looking at existing docstrings for examples, but the [napolean guide to numpy style docstrings][] is also a great source.
If you're unfamiliar with the reStructuredText (rST) markup format, check out the [Sphinx rST primer][].

Look at [`sc.tl.leiden`](https://github.com/scverse/scanpy/blob/350c3424d2f96c4a3a7bb3b7d0428d38d842ebe8/src/scanpy/tools/_leiden.py#L49-L120) as an example of a complete doctring.

[napolean guide to numpy style docstrings]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy
[sphinx rst primer]: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html


#### `Params` section

The `Params` abbreviation is a legit replacement for `Parameters`.

To document parameter types use type annotations on function parameters.
These will automatically populate the docstrings on import, and when the documentation is built.

Use the python standard library types (defined in `collections.abc` and `typing` modules) for containers, e.g.  

+ `collections.abc.Sequence`s (like `list`),
+ `collections.abc.Iterable`s (like `set`), and
+ `collections.abc.Mapping`s (like `dict`).

Always specify what these contain, e.g. `{'a': (1, 2)}` → `Mapping[str, Tuple[int, int]]`.
If you can’t use one of those, use a concrete class like `AnnData`.
If your parameter only accepts an enumeration of strings, specify them like so: `Literal['elem-1', 'elem-2']`.

#### `Returns` section

+ Function returns nothing? Use None.
+ Single object: pd.DataFrame — description on the next line.
+ Multiple values:
+ Prefer a named tuple/dataclass.
+ Otherwise list each element on its own line:

```
Returns
-------
norm : AnnData
    Normalized copy of the input.
stats : dict
    Summary statistics (mean, var, n_cells).
```
<br>


## Releasing and versioning
---

Versioning and release tagging in G4X-helpers is handled through `bump-my-version` and maintainers handle version bumps and publishing releases.  

Please do not change the project version or the changelog. 

Just submit your code/docs and they will be incorporated into our release workflow.

--8<-- "_partials/end_cap.md"
