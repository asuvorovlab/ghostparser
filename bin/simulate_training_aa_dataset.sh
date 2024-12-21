#!/bin/bash
#set -e  # Exit immediately on error

# Check if the input directory and number of simulations are provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <directory_with_python_files> <number_of_simulations>"
    exit 1
fi

# Set input arguments
export SIM_DIR="$1"
export NUM_SIM="$2"
export SCRIPT_DIR=$(dirname "$0")

# Move each Python file to its own directory if not already done
for file in "$SIM_DIR"/*.py; do
    i=$(basename "$file" .py)
    if [[ ! -d "$SIM_DIR/$i" ]]; then
        mkdir -p "$SIM_DIR/$i"
        mv "$file" "$SIM_DIR/$i"
    fi
done

# Define the function to be parallelized
run_simulation() {
    local i="$1"
    local x="$2"

    # Paths for the script and output files
    local script_path="$SIM_DIR/$i/$i.py"
    local output_file="$SIM_DIR/$i/tree${x}.txt"
    local rescaled_file="$SIM_DIR/$i/rescaled.tre"
    local nuc_alignment="$SIM_DIR/$i/nuc_alignment_$x"
    local nuc_trimmed="$SIM_DIR/$i/nuc_alignment_trimmed_$x"
    local nuc_realigned="$SIM_DIR/$i/nuc_realigned_$x"
    local nuc_tree="$SIM_DIR/$i/nuc_tree_$x"
    local lock_file="$output_file.lock"

    if [[ -f "$script_path" ]]; then

        # Step 1: Run the Python simulation if output file does not exist
        (
            flock -n 200
            if [[ ! -s "$output_file" ]]; then
                echo "Running simulation for $script_path (Tree ${x})"
                python "$script_path" | grep "&R" | cut -d" " -f7- | \
                    sed 's/n0/A/g' | sed 's/n1/B/g' | sed 's/n2/C/g' | sed 's/n3/Ghost/g' | sed 's/n4/Out/g' > "$output_file"

                # Check if the output file is empty
                if [[ ! -s "$output_file" ]]; then
                    echo "Warning: No valid tree output for $script_path"
                    return
                fi
            else
                echo "Checkpoint: Tree output already exists for $output_file"
            fi
        ) 200>"$lock_file"

        # Step 2: Rescale the tree if not already done
        (
            flock -n 201
            if [[ ! -s "$rescaled_file" ]]; then
                echo "Rescaling tree for $output_file"
                python "$SCRIPT_DIR/rescale_tree.py" "$output_file" "$rescaled_file"
            else
                echo "Checkpoint: Rescaled tree already exists for $rescaled_file"
            fi
        ) 201>"$rescaled_file.lock"

        # Step 3: Run IQ-TREE simulations if alignment file does not exist
        (
            flock -n 202
            if [[ ! -s "$nuc_tree.treefile" ]]; then
                echo "Running nucleotide alignment simulation for $rescaled_file"
                iqtree --alisim "$nuc_alignment" -t "$rescaled_file" -m LG --indel .03,.1 --indel-size NB{5/20},POW{1.5/5} --length 350 --out-format fasta 
                seqkit grep -v -p "Ghost" "$nuc_alignment.unaligned.fa" > "$nuc_trimmed"
                if [[ ! -s "$nuc_realigned" ]]; then
                    linsi "$nuc_trimmed" > "$nuc_realigned"
                fi
                iqtree -s "$nuc_realigned" -m TEST --prefix "$nuc_tree"
            else
                echo "Checkpoint: Nucleotide alignment already exists for $nuc_tree"
            fi
        ) 202>"$nuc_tree.lock"

    else
        echo "Error: Script not found: $script_path"
    fi
}

# Export the function for use with GNU Parallel
export -f run_simulation

# Run the simulations in parallel for each subdirectory
parallel --env SIM_DIR --env SCRIPT_DIR -j "${SLURM_CPUS_PER_TASK:-1}" run_simulation ::: $(find "$SIM_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;) ::: $(seq 1 "$NUM_SIM")

