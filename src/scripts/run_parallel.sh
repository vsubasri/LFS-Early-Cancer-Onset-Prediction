#!/bin/bash
#
# Script to run the model in parallel with multiple seeds
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

# Parameters
NUM_SEEDS=10
MAX_PARALLEL_JOBS=4

echo "Starting parallel model execution with ${NUM_SEEDS} seeds..."

# Create an array to store background PIDs
pids=()

# Run models with different seeds in parallel
for seed in $(seq 1 ${NUM_SEEDS}); do
    # Create seed-specific output directory
    SEED_DIR="${RESULTS_DIR}/seed_${seed}"
    mkdir -p ${SEED_DIR}
    
    echo "Starting model with seed ${seed}..."
    
    # Run the model with this seed
    Rscript ${SRC_DIR}/model/runSingleModel.R \
        --input-data="${DATA_DIR}/processed/preprocessed_data.RData" \
        --output-dir="${SEED_DIR}" \
        --checkpoint-dir="${CHECKPOINT_DIR}/seed_${seed}" \
        --seed=${seed} &
    
    # Store the PID
    pids+=($!)
    
    # If we've reached the maximum number of parallel jobs, wait for one to finish
    if [ ${#pids[@]} -ge ${MAX_PARALLEL_JOBS} ]; then
        wait -n
        # Remove completed processes from the array
        for i in "${!pids[@]}"; do
            if ! kill -0 ${pids[$i]} 2>/dev/null; then
                unset pids[$i]
            fi
        done
        # Recreate array to remove gaps
        pids=("${pids[@]}")
    fi
done

# Wait for remaining jobs to finish
wait

echo "All model executions completed"

# Combine results
echo "Combining results..."
Rscript ${SRC_DIR}/model/combine_results.R \
    --input-dir="${RESULTS_DIR}" \
    --output-file="${RESULTS_DIR}/combined_results.RData"

echo "Results saved to ${RESULTS_DIR}/combined_results.RData" 