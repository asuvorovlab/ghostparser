#!/bin/bash

# Ensure $1 and $2 are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <path_to_simulations> <number_of_simulations>"
    exit 1
fi

# Convert SIM_DIR to an absolute path
SIM_DIR="$(cd "$1" && pwd)"
export SIM_DIR  # Export so it's available to parallel processes
NUM_SIMULATIONS="$2"

# Get the absolute path of the current script
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "SCRIPT_DIR: $SCRIPT_DIR"
export SCRIPT_DIR
echo "SIM_DIR: $SIM_DIR"

# Function to process a single simulation
process_simulation() {
    i="$1"
    SIM_PATH="${SIM_DIR}/simulation_${i}"
    COMBINED_TREES="${SIM_PATH}/combined.trees"
    SUMMARY_STATS="${SIM_PATH}/summary_stats.csv"

    echo "Processing simulation $i"
    echo "SIM_PATH: $SIM_PATH"
    echo "COMBINED_TREES: $COMBINED_TREES"
    echo "SUMMARY_STATS: $SUMMARY_STATS"

    # Skip if summary_stats.csv already exists
    if [ -f "$SUMMARY_STATS" ]; then
        echo "Skipping simulation $i: summary_stats.csv already exists."
        return 0
    fi

    # Combine tree files and delete other files
    if [ -d "$SIM_PATH" ]; then
        cat "${SIM_PATH}"/*treefile > "$COMBINED_TREES"
        find "$SIM_PATH" -type f ! -name "combined.trees" -delete
    else
        echo "Error: SIM_PATH does not exist: $SIM_PATH"
        return 1
    fi

    # Run Python script to get absolute branch lengths
    if [ -f "$SCRIPT_DIR/relative_lengths.py" ]; then
        python3 "$SCRIPT_DIR/relative_lengths.py" "$COMBINED_TREES" "$SUMMARY_STATS"
    else
        echo "Error: relative_lengths.py not found in $SCRIPT_DIR"
        return 1
    fi
}

# Export the function so it works in the parallel subshell
export -f process_simulation

# Run the simulations in parallel
seq 1 "$NUM_SIMULATIONS" | parallel -j "$(nproc)" process_simulation

