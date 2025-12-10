#!/usr/bin/env bash

# Usage note:
# ./tests/get_demo_data.sh                       # downloads the default demo_data_1x1 dataset
# ./tests/get_demo_data.sh another_dataset_id    # downloads a different dataset by ID

set -euo pipefail

# Resolve to the directory where this script lives (./tests)
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "$SCRIPT_DIR"

# Configurable dataset ID (defaults to demo_data_1x1)
DATASET_ID="${1:-demo_data_1x1}"

S3_BUCKET_ROOT="s3://sg-data-portal-public"
CONFIGS_PREFIX="${S3_BUCKET_ROOT}/configs"

# Clean up from previous runs (optional; drop if you want to keep old data)
rm -rf test_data datasets

mkdir -p datasets

echo "Downloading datasets.json from S3..."
aws s3 cp "${CONFIGS_PREFIX}/datasets.json" "datasets/config.json"

echo "Collecting MD5 checksum for dataset ${DATASET_ID}..."
jq -r \
  ".datasets[] | select(.id==\"${DATASET_ID}\") | \"\(.dataChecksum)  \(.objectKey)\"" \
  datasets/config.json > datasets/test_data.md5

file_key=$(awk '{print $2}' datasets/test_data.md5)
basename=${file_key##*/}

echo "Downloading ${basename} from S3..."
aws s3 cp "${S3_BUCKET_ROOT}/${file_key}" ${file_key} --no-progress

echo "Verifying MD5 checksum..."
md5sum -c datasets/test_data.md5

# Clean up config and checksum file
rm datasets/config.json datasets/test_data.md5

# Match original layout: rename datasets -> test_data
mv datasets test_data

echo "Extracting ${basename}..."
cd test_data
tar -xf "${basename}"
rm "${basename}"
echo "Contents of ./tests/test_data:"
ls -lh .
