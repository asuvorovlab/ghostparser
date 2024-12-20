#!/bin/bash

# Export input and output directories as arguments
export INPUT_DIR="$1"
export OUTPUT_DIR="$2"

# Detect the number of available CPUs for parallel processing
NUM_CPUS=${SLURM_CPUS_PER_TASK:-$(nproc)}

# Create the output directory and a dedicated treefile directory if they do not exist
mkdir -p "$OUTPUT_DIR/treefiles"

# Function to process each file
process_file() {
    local file="$1"
    local new_name=$(basename "$file" | cut -d"." -f1)
    local unique_name="${new_name}_$(date +%s%N)"  # Unique filename to avoid conflicts
    
    echo "Processing file: $file"
    
    # Run linsi and iqtree, checking for errors at each step
    if linsi "$file" > "$OUTPUT_DIR/treefiles/${unique_name}_aligned.fa"; then
        echo "Alignment successful for: $file"
    else
        echo "Error in linsi for: $file" >&2
        return 1
    fi
    
    if iqtree -s "$OUTPUT_DIR/treefiles/${unique_name}_aligned.fa" -m TEST > "$OUTPUT_DIR/treefiles/${unique_name}_iqtree.log" 2>&1; then
        echo "IQ-TREE successful for: $file"
    else
        echo "Error in IQ-TREE for: $file" >&2
        return 1
    fi
}

# Export the function for use in GNU Parallel
export -f process_file  

# Use GNU Parallel to process all files in the input directory
find "$INPUT_DIR" -type f | parallel --jobs "$NUM_CPUS" --bar process_file {}

# Wait for all processes to complete
wait

# Concatenate all treefiles into a single file in the output directory
cat "$OUTPUT_DIR/treefiles"/*treefile > "$OUTPUT_DIR/calculated_trees.txt"

echo "All tasks completed successfully."
