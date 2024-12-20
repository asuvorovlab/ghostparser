#!/bin/bash

# Input and output directories are passed as positional arguments (not using 'export')
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Detect the number of available CPUs for parallel processing, defaulting to nproc if not in SLURM
NUM_CPUS=${SLURM_CPUS_PER_TASK:-$(nproc)}

# Ensure the output directories exist
mkdir -p "$OUTPUT_DIR/treefiles"

# Function to process each file
process_file() {
    local file="$1"
    local new_name=$(basename "$file" | cut -d"." -f1)
    local unique_name="${new_name}_$(uuidgen)"  # Use UUID to ensure uniqueness
    local output_file_path="$OUTPUT_DIR/treefiles/${unique_name}_aligned.fa"
    
    echo "Processing file: $file"
    
    # Run linsi and iqtree, checking for errors at each step
    if linsi "$file" > "$output_file_path"; then
        echo "Alignment successful for: $file"
    else
        echo "Error in linsi for: $file" >&2
        return 1
    fi
    
    if iqtree -s "$output_file_path" -m TEST > "$OUTPUT_DIR/treefiles/${unique_name}_iqtree.log" 2>&1; then
        echo "IQ-TREE successful for: $file"
    else
        echo "Error in IQ-TREE for: $file" >&2
        return 1
    fi
}

# Export the function for use in GNU Parallel
export -f process_file  

# Ensure the PATH is set correctly for commands like linsi, iqtree, parallel, etc.
export PATH="$PATH:/usr/local/bin:/usr/bin"

# Use GNU Parallel to process all files in the input directory
find "$INPUT_DIR" -type f | parallel --jobs "$NUM_CPUS" --bar --joblog "$OUTPUT_DIR/parallel.log" process_file {}

# Wait for all processes to complete (NOTE: parallel already waits, so this is optional)
wait

# Concatenate all treefiles into a single file in the output directory
cat "$OUTPUT_DIR/treefiles"/*treefile > "$OUTPUT_DIR/calculated_trees.txt"

# Set file permissions to ensure no permission issues
chmod -R 775 "$OUTPUT_DIR/treefiles"

echo "All tasks completed successfully."

