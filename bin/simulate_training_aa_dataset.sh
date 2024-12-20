#!/bin/bash

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
for i in $(ls "$SIM_DIR" | sed 's/.py//g'); do
    if [[ ! -d "$SIM_DIR/$i" ]]; then
        mkdir -p "$SIM_DIR/$i"
        mv "$SIM_DIR/$i.py" "$SIM_DIR/$i"
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
    local aa_alignment="$SIM_DIR/$i/aa_alignment_$x"
    local aa_trimmed="$SIM_DIR/$i/aa_alignment_trimmed_$x"
    local aa_realigned="$SIM_DIR/$i/aa_realigned_$x"
    local aa_tree="$SIM_DIR/$i/aa_tree_$x"
    # Uncomment if using amino acid alignment
    # local aa_alignment="$SIM_DIR/$i/aa_alignment_$x"
    # local nuc_tree="$SIM_DIR/$i/aa_tree_$x"

    # Check if the script exists
    if [[ -f "$script_path" ]]; then
        # Step 1: Run the Python simulation if output file does not exist
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

        # Step 2: Rescale the tree if not already done
        if [[ ! -s "$rescaled_file" ]]; then
            echo "Rescaling tree for $output_file"
            python $SCRIPT_DIR/rescale_tree.py "$output_file" "$rescaled_file"
        else
            echo "Checkpoint: Rescaled tree already exists for $rescaled_file"
        fi

        # Step 3: Run IQ-TREE simulations if alignment file does not exist
#        if [[ ! -s "$nuc_tree.treefile" ]]; then
 #           echo "Running nucleotide alignment simulation for $rescaled_file"
#############Right now the indel model follows the Zipfian distribution as defined in INDELible, with an empirical exponent of 1.5 and a max indel size of 5. There are on average###################
  #          iqtree --alisim "$nuc_alignment" -t "$rescaled_file" -m GTR+G4 --indel .03,.1 --indel-size NB{5/20},POW{1.5/5} --length 1000 --out-format fasta 
#	    seqkit grep -v -p "Ghost" "$nuc_alignment.unaligned.fa" > "$nuc_trimmed"
#	    linsi "$nuc_trimmed" > "$nuc_realigned"
 #           iqtree -s "$nuc_realigned" -m TEST --prefix "$nuc_tree" 
            

#        else
 #           echo "Checkpoint: Nucleotide alignment already exists for $nuc_alignment"
  #      fi

        # Uncomment if running amino acid alignments
         if [[ ! -s "$aa_tree.treefile" ]]; then
             echo "Running amino acid alignment simulation for $rescaled_file"

             iqtree --alisim "$aa_alignment" -t "$rescaled_file" -m LG --length 350 --indel .03,.1 --indel-size NB{5/20},POW{1.5/5} --out-format fasta
             seqkit grep -v -p "Ghost" "$aa_alignment.unaligned.fa" > "$aa_trimmed"
             linsi "$aa_trimmed" > "$aa_realigned"
             iqtree -s "$aa_realigned" -m TEST --prefix "$aa_tree" 


         else
             echo "Checkpoint: Amino acid alignment already exists for $aa_alignment"
         fi

    else
        echo "Error: Script not found: $script_path"
    fi
}

# Export the function for use with GNU Parallel
export -f run_simulation

# Run the simulations in parallel for each subdirectory
for i in $(ls "$SIM_DIR"); do
    parallel -j "$(nproc)" run_simulation ::: "$i" ::: $(seq 1 "$NUM_SIM")
done

