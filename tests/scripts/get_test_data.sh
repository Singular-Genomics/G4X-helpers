#!/usr/bin/env bash

# Usage note:
# ./tests/get_test_data.sh                       # downloads the default demo_data_1x1 dataset
# ./tests/get_test_data.sh another_dataset_id    # downloads a different dataset by ID

set -euo pipefail

# Resolve to the directory where this script lives and then move to tests/
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
TESTS_DIR="${SCRIPT_DIR%/*}"

cd "$TESTS_DIR"

# Configurable dataset ID (defaults to demo_data.tar)
DATASET_ID="${1:-demo_data.tar}"
S3_BUCKET_ROOT="s3://sg-data-portal-public"

# Clean up from previous runs (should be optional)
rm -rf datasets
mkdir -p datasets

echo "Downloading ${DATASET_ID} from S3..."
aws s3 cp "${S3_BUCKET_ROOT}/datasets/${DATASET_ID}" "datasets/${DATASET_ID}" #--no-progress

# Clean up config and checksum file
mv "datasets/${DATASET_ID}" datasets/test_data.tar
