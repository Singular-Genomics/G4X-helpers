# Changelog

## [Unreleased] — `v4.0.0`

#### Overview

G4X-helpers v4 is a major architectural release designed to support the latest G4X-data v4 format and the new zarr based G4X-viewer input.
This release reorganizes the Python API, processing modules, CLI, schema validation, and G4X-viewer generation. It introduces a more consistent workflow for processing, validating, migrating, and modifying G4X datasets.

#### Highlights
- continued support to migrate legacy data to latest v4 format
- begin phase-out of jp2 files (will be replaced with ome.tiff during migration)
- new functions to inspect and edit G4X-viewer input
- now handles single-cell processing steps after redemux and resegment operations
- added optional GPU acceleration through the `gpu` dependency extra.
- added region-of-interest support for processing cropped datasets.
- expanded support for transcript-only, protein-only, and combined assay data.

#### CLI

- Reorganized the CLI around five supported commands:
  - `redemux`
  - `resegment`
  - `migrate`
  - `validate`
  - `viewer`
- Added `migrate --check` for checking migration compatibility without modifying data.
- Added `migrate --roi` for migrating a selected region.
- Added branch-based output handling for `redemux` and `resegment`.
- Added `--no-downstream` for skipping downstream processing.
- Added `validate --raw-only` for validating only the inputs required by G4X-helpers.
- Added viewer metadata import and export commands for images, cells, and transcripts.
- Improved CLI help messages, progress reporting, error handling, and version reporting.

#### Data processing

- Refactored demultiplexing, aggregation, and single-cell processing into independent modules.
- Added configurable demultiplexing parameters and batch handling.
- Added RNA and protein correlation analysis.
- Added assay-aware image, transcript, and cell output generation.
- Improved handling of missing bead masks, empty transcript data, failed clustering, and transcript-only runs.

#### Validation and migration

- Replaced the legacy schema implementation with dedicated file, directory, and dataset validators.
- Added more detailed validation reporting for raw and processed G4X-data.
- Reworked legacy-data migration around dedicated migrator classes.
- Improved samplesheet, manifest, transcript-table, protein-directory, and `sample.g4x` validation.
- Added migration support for cropped regions and updated G4X-viewer layouts.

#### G4X-viewer

- New G4X-viewer Zarr creation for images, transcripts, cell masks, and metadata.
- Added metadata import and export for existing viewer stores.
- Added support for multiple segmentations and assay-aware channel ordering.
- Improved image conversion, tiling, chunking, and Zarr compatibility.
- Improved transcript visualization, polygon generation, and color controls.

#### Reliability and release infrastructure

- Expanded CLI, processing, smoke, transcript-only, and protein-only tests.
- Added strict documentation, package, wheel, and multi-architecture Docker preflight checks.
- Added AMD64 and ARM64 container validation.
- Refreshed dependencies to address reported security vulnerabilities.
- Removed development and documentation dependencies from production container images.
- Updated GitHub Actions to use Node.js 24-compatible action versions.


## [2026-04-13] — `v3.0.2`

- fix: load adata from correct location

## [2026-04-03] — `v3.0.1`

- update changelog
- remove "staging" as release-branch in uv-ship

## [2026-04-03] — `v3.0.0`

- feat: compatibility with "cytoplasmic" image label
- fix: ensure new_bin feature handles init_bin logic
- fix: resegment not updating bin file
- fix: ensure new gene names propagate to sc-out during redemux
- fix: cannot save adata object

## [2025-12-30] — `v2.1.2`

- Fix/auto revert migrate (#68)

## [2025-12-29] — `v2.1.1`

- fix: migration tripping over file updates despite fail

## [2025-12-29] — `v2.1.0`

- feat: "validate" function for detailed schema output validation (#67)
- fix: re-demux on tx-only runs (#66)

## [2025-12-16] — `v2.0.3`

- fix: prevent validation fail for tx-only runs (#64)
- fix: code scanning alert no. 3: Workflow does not contain permissions (#62)
- pkg: incorporate test data from S3 for full feature tests (#61)
- docs: layout updates and cleaner changelog format (#63)

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
