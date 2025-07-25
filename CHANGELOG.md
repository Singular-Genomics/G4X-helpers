
# Changelog

## v0.4.12 – (2025-07-25)

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

## v0.2.4 - v0.4.11 (unreleased)

- Release preparation
- Add Dockerfile
- Update dependency information

## v0.2.3 – (2025-06-25)

- Add tar_viewer tool to tar up a G4X-Viewer folder for the single-file upload option.

## v0.2.2 – (2025-06-24)

- Add new_bin tool to more quickly generate a new bin file

## v0.2.1 – (2025-06-02)

- Bug fixes for MVP functionality

## v0.2.0 – (2025-05-16)

- Add CLI tools for resegmentation and updating bin files with clustering/embedding information
