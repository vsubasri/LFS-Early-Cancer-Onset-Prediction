#!/bin/bash
#
# Main script to run the modeling pipeline
#

# Set the base directory to the project root
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC_DIR="${BASE_DIR}/src"
DATA_DIR="${BASE_DIR}/data"
RESULTS_DIR="${DATA_DIR}/results"
CHECKPOINT_DIR="${BASE_DIR}/checkpoint"

# Ensure directories exist
mkdir -p ${RESULTS_DIR}
mkdir -p ${CHECKPOINT_DIR}

echo "Starting modeling pipeline..."

# Run preprocessing
echo "Running preprocessing..."
${SRC_DIR}/scripts/run_preprocessing.sh

# Run model training
echo "Running model training..."
Rscript ${SRC_DIR}/model/runSingleModel.R \
  --input-data="${DATA_DIR}/processed/preprocessed_data.RData" \
  --output-dir="${RESULTS_DIR}" \
  --checkpoint-dir="${CHECKPOINT_DIR}" \
  --seed=42

echo "Model execution completed"
echo "Results saved to ${RESULTS_DIR}" 