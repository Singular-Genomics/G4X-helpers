# Changelog

## [2025-12-09] — `v2.0.2`

- docs: update to match global layout
- docs: set up section indices correctly
- docs: correct readme reference
- pkg: clean up changelog
- pkg: pre-commit

## [2025-12-04] — `v2.0.1`

- fix: docker release auth issue

## [2025-12-04] — `v2.0.0`

#### Overview:
- New "migrate" function to update older datasets to be compatible with the latest G4X-viewer and G4X-helpers versions.
- Added schema validation to ensure data integrity and compatibility.
- Implementation of an initial test framework for the package.
- Re-worked CLI-api to be more modular and easier to extend.

#### Changes:
- new "in_place" output flag replaces "out_dir" and ensures predictable output behavior.
- standardized log location
- new logging formatter
- refactor "workflows" into modules
- adapt G4Xoutput to latest schema
- refactored stream_features method to redemux module
- new write_csv_gz method
- add bead_mask loading to G4Xoutput
- removed deprecated G4X-viewer schema
- new dependency: pathschema
- updates dependencies
- update docs content and match structure to main site
- update changelog

## [2025-12-03] — `v1.0.2`

- fix: update uv.lock

## [2025-12-01] — `v1.0.1`

- fix: get_shape fast via glymur

## [2025-11-13] — `v1.0.0`

- docs: add partials for shared options
- pkg: add rich-codex as docs dependency
- fix: load_segmentation method now allows custom keys
- feat: redemux parses transcript panel for input flexibility
- docs: update for new CLI
- feat: numpy version agnostic npzGetshape
- pkg: repr update
- pkg: implement new cli and re-factor
- pkg: set python 3.12 as default version

## [2025-10-01] — `v0.5.2`

- fix: missing tests in publish workflow

## [2025-10-01] — `v0.5.1`

- pkg: added PyPI publish workflow
- fix: README logo link

## [2025-10-01] — `v0.5.0`

- docs: formatting
- docs: spell check
- docs: integrate redemux docs
- docs: Docs restructuring (#33)
- docs: updated _core docs

- feat: Redemux tool (#30)
- feat: progress bar

- fix: G4Xoutput out_dir now defaults to cwd instead of run_base
- fix: enforce cluster color as hex-codes
- fix: auto-conversion of clusters to 'categorical' when load_clustering=true
- fix: change all segmentation_cell_id to cell_id (#34)
- fix: cellid_key - bugfix

- pkg: updated ruff
- pkg: relaxing glymur dependency
- pkg: relaxed protobuf requirement for internal alignment
- pkg: opening python requirements to include 3.11
- pkg: removed bump-my-version and replacing with uv-ship

##  [2025-08-12] — `v0.4.14`

#### Fixes:
- `tar_viewer` will now exit leaving the source folder unchanged

#### Improvements:
- `tar_viewer` now requires an out_path


## [2025-08-11] — `v0.4.13`

#### Docs:
- changelog now included in docs
- url updates
- typos

#### Package:
- deployment of multi-arch packages
- replaced dependency `ray` with `multiprocessing`
- set required `uv` version


##  [2025-07-25] — `v0.4.12`

#### Docs:
- Added Documentation for G4X-helpers and G4X-output

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


## [unreleased changes]

- Release preparation
- Add Dockerfile
- Update dependency information
- Add tar_viewer tool to tar up a G4X-viewer folder for the single-file upload option.
- Add new_bin tool to more quickly generate a new bin file
- Bug fixes for MVP functionality
- Add CLI tools for re-segmentation and updating bin files with clustering/embedding information
