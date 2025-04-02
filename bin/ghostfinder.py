#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run GhostBuster workflow by generating intermediate taxa files, "
                    "invoking ghostbuster.py, and filtering the output."
    )
    parser.add_argument("--out_taxa", required=True,
                        help="File with outgroup taxa, one per line.")
    parser.add_argument("--ghostbuster_output", required=True,
                        help="GhostBuster output file.")
    parser.add_argument("--input_trees", required=True,
                        help="File containing the gene trees, one per line.")
    parser.add_argument("--output_file", required=True,
                        help="Path to the final output file.")
    return parser.parse_args()

def generate_taxa_files(ghostbuster_output, work_dir):
    """
    Process ghostbuster_output to create three files:
      - GhostFinder_A_taxa.txt: from token2 (second token) of field1
      - GhostFinder_B_taxa.txt: from token3 (third token) of field1
      - additional_taxa_to_exclude.txt: from token1 (first token) of field1
    (Only lines with field2 equal to "Evidence of unsampled introgression" are used.)
    """
    a_set = set()
    b_set = set()
    additional_set = set()
    
    with open(ghostbuster_output, 'r') as fin:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            # Split by tab; expect at least two fields.
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            # Use strip to remove any extra spaces.
            if fields[1].strip() == "Evidence of unsampled introgression":
                tokens = fields[0].split()
                if len(tokens) >= 2:
                    a_set.add(tokens[1])
                if len(tokens) >= 3:
                    b_set.add(tokens[2])
                if len(tokens) >= 1:
                    additional_set.add(tokens[0])
    
    a_file = os.path.join(work_dir, "GhostFinder_A_taxa.txt")
    b_file = os.path.join(work_dir, "GhostFinder_B_taxa.txt")
    additional_file = os.path.join(work_dir, "additional_taxa_to_exclude.txt")
    
    with open(a_file, 'w') as fout:
        for taxon in sorted(a_set):
            fout.write(f"{taxon}\n")
    with open(b_file, 'w') as fout:
        for taxon in sorted(b_set):
            fout.write(f"{taxon}\n")
    with open(additional_file, 'w') as fout:
        for taxon in sorted(additional_set):
            fout.write(f"{taxon}\n")
    
    return a_file, b_file, additional_file

def generate_ghostfinder_c_taxa(a_file, b_file, additional_file, out_taxa, input_trees, work_dir):
    """
    Build GhostFinder_C_taxa.txt by:
      1. Building an exclusion set from:
         - All taxa from GhostFinder_B_taxa.txt,
         - additional_taxa_to_exclude.txt,
         - GhostFinder_A_taxa.txt,
         - And (if out_taxa is literally named "Out.txt") taxa from that file.
      2. Then scanning the out_taxa and input_trees files:
         - Replace any of the characters ( ) , ; with spaces.
         - Split the line into tokens and remove any branch lengths (remove colon and following text).
         - If a token is nonempty, not in the exclusion set, and not solely numeric,
           it is included in the output.
    """
    exclusion = set()
    
    def add_file_to_exclusion(filename):
        try:
            with open(filename, 'r') as fin:
                for line in fin:
                    token = line.strip()
                    if token:
                        exclusion.add(token)
        except Exception as e:
            print(f"Error reading {filename}: {e}", file=sys.stderr)
    
    # Process the first three files.
    add_file_to_exclusion(b_file)
    add_file_to_exclusion(additional_file)
    add_file_to_exclusion(a_file)
    # Also add tokens from out_taxa if its basename is "Out.txt"
    if os.path.basename(out_taxa) == "Out.txt":
        add_file_to_exclusion(out_taxa)
    
    result_tokens = set()
    
    def process_file(filename):
        with open(filename, 'r') as fin:
            for line in fin:
                # Replace the characters ( ) , ; with a space.
                line_clean = re.sub(r"[(),;]", " ", line)
                for token in line_clean.split():
                    # Remove branch lengths (anything after a colon).
                    token = token.split(":", 1)[0].strip()
                    # Skip empty tokens, tokens in exclusion, or tokens that are purely numbers.
                    if token and token not in exclusion and not re.fullmatch(r"\d+", token):
                        result_tokens.add(token)
    
    process_file(out_taxa)
    process_file(input_trees)
    
    c_file = os.path.join(work_dir, "GhostFinder_C_taxa.txt")
    with open(c_file, 'w') as fout:
        for token in sorted(result_tokens):
            fout.write(f"{token}\n")
    
    return c_file

def run_ghostbuster(script_dir, out_taxa, a_file, b_file, c_file, input_trees, full_output_file):
    """
    Run ghostbuster.py from the same directory as this script.
    Pass along the necessary file arguments.
    """
    ghostbuster_script = os.path.join(script_dir, "ghostbuster.py")
    cmd = [
        sys.executable,  # Use the same Python interpreter.
        ghostbuster_script,
        "--out_taxa", out_taxa,
        "--A_taxa", a_file,
        "--B_taxa", b_file,
        "--C_taxa", c_file,
        "--input_trees", input_trees,
        "--output_file", full_output_file
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running ghostbuster.py:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

def filter_output(full_output_file, final_output_file):
    """
    Read the ghostbuster output and filter only lines containing " sampled".
    Write these lines to the final output file.
    """
    with open(full_output_file, 'r') as fin, open(final_output_file, 'w') as fout:
        for line in fin:
            if " sampled" in line: # Also add syntax that proportion of AC is greater than BC!!!
                fout.write(line)

def remove_intermediate_files(files):
    for file in files:
        try:
            os.remove(file)
        except Exception as e:
            print(f"Warning: could not remove file {file}: {e}", file=sys.stderr)

def main():
    args = parse_args()
    # Use the current working directory for intermediate files.
    work_dir = os.getcwd()
    # Get the directory where this script resides (for locating ghostbuster.py)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Step 1: Generate intermediate taxa files from ghostbuster_output.
    a_file, b_file, additional_file = generate_taxa_files(args.ghostbuster_output, work_dir)
    
    # Step 2: Build GhostFinder_C_taxa.txt using the out_taxa and input_trees files.
    c_file = generate_ghostfinder_c_taxa(a_file, b_file, additional_file,
                                         args.out_taxa, args.input_trees, work_dir)
    
    # Step 3: Run ghostbuster.py, outputting to Full_GhostFinder_output.txt.
    full_output_file = os.path.join(work_dir, "Full_GhostFinder_output.txt")
    run_ghostbuster(script_dir, args.out_taxa, a_file, b_file, c_file,
                    args.input_trees, full_output_file)
    
    # Step 4: Filter the ghostbuster output to only include lines with " sampled"
    filter_output(full_output_file, args.output_file)
    
    # Step 5: Remove all intermediate files.
   # intermediate_files = [a_file, b_file, additional_file, c_file, full_output_file]
   # remove_intermediate_files(intermediate_files)

if __name__ == "__main__":
    main()
