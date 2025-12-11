#!/usr/bin/env bash
set -euo pipefail

# Resolve to the directory where this script lives (./tests)
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "$SCRIPT_DIR"
cd ./datasets/test_data

uv run g4x-helpers -v 2 migrate . 
uv run g4x-helpers -v 2 tar_viewer . 
uv run g4x-helpers -v 2 new_bin . 
uv run g4x-helpers -v 2 update_bin . --metadata single_cell_data/clustering_umap.csv.gz --cluster-key leiden_0.400 --cellid-key label 
uv run g4x-helpers -v 2 resegment . --cell-labels segmentation/segmentation_mask.npz --labels-key nuclei_exp 
uv run g4x-helpers -v 2 redemux . --manifest transcript_panel.csv 
