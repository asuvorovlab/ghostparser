#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys
import tempfile
import concurrent.futures
import threading
from functools import partial

# Global lock to serialize access to the shared tree_stats file.
lock = threading.Lock()

def usage_error():
    sys.stderr.write("Usage: {} --out_taxa <value> --A_taxa <value> --B_taxa <value> --C_taxa <value> --input_trees <value> --output_file <value> --threads <value> (optional)\n".format(sys.argv[0]))
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--out_taxa", help="File containing each of the Outgroup taxa, one per line")
    parser.add_argument("--A_taxa", help="File containing each of the A taxa, one per line")
    parser.add_argument("--B_taxa", help="File containing each of the B taxa, one per line")
    parser.add_argument("--C_taxa", help="File containing each of the C taxa, one per line")
    parser.add_argument("--input_trees", help="File containing the gene trees, one per line")
    parser.add_argument("--output_file", help="User specified path to an output file")
    # New --threads option; default to using all available threads.
    parser.add_argument("--threads", type=int, default=os.cpu_count(),
                        help="Number of threads to use (default: all available threads)")
    args, unknown = parser.parse_known_args()
    if unknown:
        sys.stderr.write("Error: Unknown option {}\n".format(" ".join(unknown)))
        usage_error()
    if (not args.out_taxa or not args.A_taxa or not args.B_taxa or 
        not args.C_taxa or not args.input_trees or not args.output_file):
        sys.stderr.write("Error: Missing required arguments.\n")
        usage_error()
    return args

def read_file_lines(filepath):
    with open(filepath, "r") as f:
        # Skip empty lines and remove whitespace
        return [line.strip() for line in f if line.strip()]

def check_taxa_in_trees(taxa_files, input_trees_file):
    # Read the entire input trees file into a single string
    with open(input_trees_file, "r") as f:
        trees_content = f.read()
    missing_taxa = []
    for taxa_file in taxa_files:
        for taxon in read_file_lines(taxa_file):
            # Using regex for whole-word matching
            if not re.search(r'\b{}\b'.format(re.escape(taxon)), trees_content):
                missing_taxa.append(taxon)
    if missing_taxa:
        sys.stderr.write("Error: The following taxa were not found in any input tree:\n")
        for taxon in missing_taxa:
            sys.stderr.write(" - {}\n".format(taxon))
        sys.exit(1)

def generate_combinations(C_taxa_file, B_taxa_file, A_taxa_file, out_taxa_file):
    # Read the taxa from each file
    taxa_C = read_file_lines(C_taxa_file)
    taxa_B = read_file_lines(B_taxa_file)
    taxa_A = read_file_lines(A_taxa_file)
    taxa_out = read_file_lines(out_taxa_file)

    combinations = []
    # Nested loops in the same order as the bash script:
    # Outer loop: C_taxa, then B_taxa, then A_taxa, then out_taxa.
    for c in taxa_C:
        for b in taxa_B:
            for a in taxa_A:
                for out in taxa_out:
                    # Ensure all four taxa are distinct
                    if len({c, b, a, out}) == 4:
                        # The bash prints the line as: a (from C_taxa), b (from B_taxa), c (from A_taxa), d (from out_taxa)
                        # Later, the variables are re-assigned so that:
                        #   out = 4th field, c = 1st, b = 2nd, a = 3rd.
                        combinations.append((c, b, a, out))
    return combinations

def run_command(command_list, capture_output=False):
    try:
        if capture_output:
            result = subprocess.check_output(command_list)
            return result.decode("utf-8")
        else:
            subprocess.check_call(command_list)
            return ""
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running command: {}\n".format(" ".join(command_list)))
        sys.exit(1)

def process_combination(line_tuple, args, script_dir):
    # Unpack the combination tuple (as generated above):
    # tuple order: (c, b, a, out) where:
    #   c from C_taxa, b from B_taxa, a from A_taxa, out from out_taxa.
    c_taxon, b_taxon, a_taxon, out_taxon = line_tuple
    # For output, the original bash saves the combination line as: "c b a out"
    combination_line = "{} {} {} {}".format(c_taxon, b_taxon, a_taxon, out_taxon)

    # Run nw_prune -v input_trees out a b c
    nw_command = ["nw_prune", "-v", args.input_trees, out_taxon, a_taxon, b_taxon, c_taxon]
    nw_output = run_command(nw_command, capture_output=True)

    # Perform the sed-like replacements in the order given:
    # replace a_taxon -> "A", b_taxon -> "B", c_taxon -> "C", out_taxon -> "Out"
    triplet_content = nw_output.replace(a_taxon, "A")
    triplet_content = triplet_content.replace(b_taxon, "B")
    triplet_content = triplet_content.replace(c_taxon, "C")
    triplet_content = triplet_content.replace(out_taxon, "Out")

    # Write the output to a unique temporary file instead of a fixed "triplet.txt"
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, prefix="triplet_", suffix=".txt")
    tmp.write(triplet_content)
    tmp.close()
    triplet_filename = tmp.name

    try:
        # The tree stats file is shared, so acquire a lock to prevent concurrent access.
        with lock:
            tree_stats_file = args.input_trees + ".tree_stats.txt"
            process_trees_cmd = ["python", os.path.join(script_dir, "process_trees.py"), triplet_filename, tree_stats_file]
            run_command(process_trees_cmd)
    
            branch_stats_cmd = ["Rscript", os.path.join(script_dir, "branch_stats2.0.r"), tree_stats_file]
            bl_stats = run_command(branch_stats_cmd, capture_output=True).strip().replace("\n", " ")
            
            # Remove the intermediate tree_stats file.
            try:
                os.remove(tree_stats_file)
            except OSError:
                pass
    
        # Process the output from the R script
        # Expected format: <chisq>~<second_field>
        fields = bl_stats.split("~")
        chisq = fields[0].strip() if len(fields) >= 1 else ""
        second_field = fields[1].strip() if len(fields) >= 2 else ""
        # Determine if the first field contains "are significantly different"
        result_signif = ("are significantly different" in chisq)
    
        # Process second_field to extract out_pval and out_higher
        out_pval = ""
        out_higher = ""
        if second_field:
            colon_parts = second_field.split(":")
            if len(colon_parts) >= 4:
                out_pval = colon_parts[3].strip()
            # For out_higher, take the first semicolon-delimited segment then split by colon
            semicolon_parts = second_field.split(";")
            if len(semicolon_parts) >= 1:
                colon_parts_first = semicolon_parts[0].split(":")
                if len(colon_parts_first) >= 2:
                    out_higher = colon_parts_first[1].strip()
        # The entire second_field is used as blt output
        blt = second_field
    
        # Determine the conclusion based on the parsed output
        if result_signif:
            try:
                out_pval_decimal = float(out_pval)
            except ValueError:
                out_pval_decimal = 1.0  # if conversion fails, default to non-significant
            if out_pval_decimal < 0.05 and out_higher == "AC":
                conclusion = "Evidence of unsampled introgression"
            else:
                conclusion = "Evidence of sampled introgression"
        else:
            conclusion = "No evidence of introgression"
    
        # Return the combination line and results in the desired format
        # Note: The bash script outputs: combination_line, conclusion, chisq, blt (tab-separated)
        return "{}\t{}\t{}\t{}".format(combination_line, conclusion, chisq, blt)
    finally:
        # Remove the intermediate temporary file
        try:
            os.remove(triplet_filename)
        except OSError:
            pass

def main():
    args = parse_args()

    # Check that all taxa in the provided taxa files appear in the input_trees file.
    taxa_files = [args.out_taxa, args.A_taxa, args.B_taxa, args.C_taxa]
    check_taxa_in_trees(taxa_files, args.input_trees)

    # Display parsed arguments
    sys.stdout.write("Parsed arguments:\n")
    sys.stdout.write("  out_taxa: {}\n".format(args.out_taxa))
    sys.stdout.write("  A_taxa: {}\n".format(args.A_taxa))
    sys.stdout.write("  B_taxa: {}\n".format(args.B_taxa))
    sys.stdout.write("  C_taxa: {}\n".format(args.C_taxa))
    sys.stdout.write("  input_trees: {}\n".format(args.input_trees))
    sys.stdout.write("  threads: {}\n".format(args.threads))

    # Write header to the output file
    with open(args.output_file, "w") as out_f:
        header = "Triplet Tested\tConclusion\tDCT result\tBLT result\n"
        out_f.write(header)

    # Determine the directory where this script is located (for calling other scripts)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Generate every valid combination from the four taxa files.
    combinations = generate_combinations(args.C_taxa, args.B_taxa, args.A_taxa, args.out_taxa)

    # Process each combination in parallel and collect the results.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = list(executor.map(lambda comb: process_combination(comb, args, script_dir), combinations))

    # Append the results to the output file.
    with open(args.output_file, "a") as out_f:
        for result_line in results:
            out_f.write(result_line + "\n")

if __name__ == "__main__":
    main()
