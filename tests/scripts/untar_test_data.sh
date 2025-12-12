#!/usr/bin/env bash
set -euo pipefail

# Resolve to the directory where this script lives and then move to tests/
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
TESTS_DIR="${SCRIPT_DIR%/*}"

cd "$TESTS_DIR/datasets"

rm -r test_data || true
mkdir -p test_data

echo "Extracting test data..."
tar -xf "test_data.tar" -C test_data
echo "Contents of ./tests/test_data:"
ls -lh test_data
