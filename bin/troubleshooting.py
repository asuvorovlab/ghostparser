#!/usr/bin/env python3

import argparse
import os
import random
import sys
from pathlib import Path
import subprocess

# Argument parser
parser = argparse.ArgumentParser(description="Run ghostbuster alignment and tree generation")
parser.add_argument("--working_directory", type=str, required=True, help="Path for working directory for ghostbuster to create.")
parser.add_argument("--input_directory", type=str, required=True, help="Path to the input directory containing files for alignment.")
parser.add_argument("--A_taxon", type=str, required=True, help="Name of the A taxon to be replaced with 'A'.")
parser.add_argument("--B_taxon", type=str, required=True, help="Name of the B taxon to be replaced with 'B'.")
parser.add_argument("--C_taxon", type=str, required=True, help="Name of the C taxon to be replaced with 'C'.")
parser.add_argument("--Out_taxon", type=str, required=True, help="Name of the outgroup taxon to be replaced with 'Out'.")
args = parser.parse_args()
#define functions
def replace_taxa_names(input_file, output_file, a_taxon, b_taxon, c_taxon, out_taxon):
    """
    Replace taxon names in Newick trees with A, B, C, and Out and save the result.

    Parameters:
    input_file (str): Path to the input file containing Newick trees.
    output_file (str): Path to save the modified Newick trees.
    a_taxon (str): Taxon name to replace with 'A'.
    b_taxon (str): Taxon name to replace with 'B'.
    c_taxon (str): Taxon name to replace with 'C'.
    out_taxon (str): Taxon name to replace with 'Out'.
    """
    try:
        # Read the input file
        with open(input_file, 'r') as infile:
            trees = infile.readlines()

        # Replace taxon names
        modified_trees = []
        for tree in trees:
            tree = tree.strip()
            tree = tree.replace(a_taxon, 'A')
            tree = tree.replace(b_taxon, 'B')
            tree = tree.replace(c_taxon, 'C')
            tree = tree.replace(out_taxon, 'Out')
            modified_trees.append(tree)

        # Write the modified trees to the output file
        with open(output_file, 'w') as outfile:
            outfile.write('\n'.join(modified_trees))

        print(f"Modified trees saved to {output_file}")

    except Exception as e:
        print(f"Error: {e}")


# Resolve paths
working_dir = Path(args.working_directory).resolve()
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
generate_alignment_trees = os.path.join(script_dir, "generate_alignment_and_trees_only.sh")

# Run the bash script
subprocess.run(["bash", generate_alignment_trees, args.input_directory, str(working_dir)], check=True)
recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
input_trees_path = os.path.join(working_dir, "calculated_trees.txt")
replace_taxa_names(
input_trees_path,
recoded_tree_file,
args.A_taxon,
args.B_taxon,
args.C_taxon,
args.Out_taxon,
)
output_file=f"{args.working_directory}/input_tree_stats.csv"
process_trees_python = os.path.join(script_dir, "relative_lengths.py")
subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
tuned_model_file = f"{working_dir}/tuned_model.RData"
ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)
