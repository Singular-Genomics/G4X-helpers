# Changelog

## `latest` — [2025-09-29]

- G4Xoutput out_dir now defaults to cwd instead of run_base
- removed bump-my-version and replacing with uv-ship
- change all segmentation_cell_id to cell_id (#34)
- here it is
- try to fix the progress bar
- remove temporary folder
- Merge branch 'dev' of github.com:Singular-Genomics/G4X-helpers into dev
- break up re-demux by read. there may still be some issues with symlinks getting overwritten.
- updated _core
- conversion of leiden clusters when load_clustering=true
- updated ruff
- Redemux tool (#30)
- enforce cluster color as hex-codes
- Docs restructuring (#33)
- Merge pull request #28 from Singular-Genomics/cellid_key_bugfix
- actually pass cellid_key to seg_updater
- relaxing glymur dependency
- relaxed protobuf requirement for internal alignment

## `v0.4.14` — [2025-08-12]

#### Fixes:
- `tar_viewer` will now exit leaving the source folder unchanged

#### Improvements:
- `tar_viewer` now requires an out_path


## `v0.4.13` – (2025-08-11)

#### Docs:
- changelog now included in docs
- url updates
- typos

#### Package:
- deployment of multi-arch packages
- replaced dependency `ray` with `multiprocessing`
- set required `uv` version


## `v0.4.12` — [2025-07-25]

#### Docs:
- Added Documentation for G4X-helpers and G4X output

#### Package:
- Implemented bump-my-version for handling version updates and tagging
- Added pre-commit hooks for code quality checks
- Ruff for linting and formatting (invoked via pre-commit)
- Github actions for automated package builds and docs deployment

#### Improvements:
- Added npz util to speed up `G4Xoutput()` initialization

#### Fixes:
- Incorrect gzip compression on some output csvs

#### Housekeeping:
- Updated and trimmed dependencies
- Cleaned up .gitignore file
- Cleaned up README.md


## `v0.2.4` - `v0.4.11` (unreleased)

- Release preparation
- Add Dockerfile
- Update dependency information


## `v0.2.3` — [2025-06-25]

- Add tar_viewer tool to tar up a G4X-viewer folder for the single-file upload option.


## `v0.2.2` — [2025-06-24]

- Add new_bin tool to more quickly generate a new bin file


## `v0.2.1` — [2025-06-02]

- Bug fixes for MVP functionality


## `v0.2.0` — [2025-05-16]

- Add CLI tools for re-segmentation and updating bin files with clustering/embedding information
