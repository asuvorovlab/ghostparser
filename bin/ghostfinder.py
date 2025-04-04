#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="GhostFinder: Build taxa lists from a GhostBuster output, then run ghostbuster.py to analyze gene trees.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--out_taxa", required=True,
                        help="File with outgroup taxa, one per line.")
    parser.add_argument("--ghostbuster_output", required=True,
                        help="GhostBuster output file.")
    parser.add_argument("--input_trees", required=True,
                        help="File containing the gene trees, one per line.")
    parser.add_argument("--output_file", required=True,
                        help="Path to the final output file.")
    parser.add_argument("--threads", type=int, default=os.cpu_count(),
                        help="Number of threads to use for ghostbuster.py")
    return parser.parse_args()

def generate_GhostFinder_A_taxa(ghostbuster_output_file, output_filename="GhostFinder_A_taxa.txt"):
    taxa_set = set()
    with open(ghostbuster_output_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            if fields[1].strip() == "Evidence of unsampled introgression":
                tokens = fields[0].split()
                if len(tokens) >= 2:
                    taxa_set.add(tokens[1])
    print(f"[DEBUG] GhostFinder_A_taxa.txt: found {len(taxa_set)} taxa")
    with open(output_filename, "w") as out:
        for taxon in sorted(taxa_set):
            out.write(taxon + "\n")
    return output_filename

def generate_GhostFinder_B_taxa(ghostbuster_output_file, output_filename="GhostFinder_B_taxa.txt"):
    taxa_set = set()
    with open(ghostbuster_output_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            if fields[1].strip() == "Evidence of unsampled introgression":
                tokens = fields[0].split()
                if len(tokens) >= 3:
                    taxa_set.add(tokens[2])
    print(f"[DEBUG] GhostFinder_B_taxa.txt: found {len(taxa_set)} taxa")
    with open(output_filename, "w") as out:
        for taxon in sorted(taxa_set):
            out.write(taxon + "\n")
    return output_filename

def generate_additional_taxa_to_exclude(ghostbuster_output_file, output_filename="additional_taxa_to_exclude.txt"):
    taxa_set = set()
    with open(ghostbuster_output_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            if fields[1].strip() == "Evidence of unsampled introgression":
                tokens = fields[0].split()
                if tokens:
                    taxa_set.add(tokens[0])
    print(f"[DEBUG] additional_taxa_to_exclude.txt: found {len(taxa_set)} taxa")
    with open(output_filename, "w") as out:
        for taxon in sorted(taxa_set):
            out.write(taxon + "\n")
    return output_filename

def generate_GhostFinder_C_taxa(out_taxa_file, input_trees_file,
                                gf_B_taxa_file="GhostFinder_B_taxa.txt",
                                additional_exclude_file="additional_taxa_to_exclude.txt",
                                gf_A_taxa_file="GhostFinder_A_taxa.txt",
                                output_filename="GhostFinder_C_taxa.txt"):
    # Build the exclusion set from the four files.
    exclude_set = set()
    def add_file_to_exclusion(filename):
        try:
            with open(filename, "r") as f:
                for line in f:
                    line = line.strip()
                    if line:
                        token = line.split()[0]
                        exclude_set.add(token)
        except IOError:
            sys.stderr.write("Error reading file: {}\n".format(filename))
            sys.exit(1)
    add_file_to_exclusion(gf_B_taxa_file)
    add_file_to_exclusion(additional_exclude_file)
    add_file_to_exclusion(gf_A_taxa_file)
    add_file_to_exclusion(out_taxa_file)  # Adjust if out_taxa file isn’t exactly named Out.txt

    output_set = set()
    with open(input_trees_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            # Replace characters ( ) , ; with a space
            line = re.sub(r'[(),;]', ' ', line)
            fields = line.split()
            for field in fields:
                # Remove branch lengths (remove colon and everything after)
                field = re.sub(r':.*', '', field)
                if field and field not in exclude_set and not field.isdigit():
                    output_set.add(field)
    print(f"[DEBUG] GhostFinder_C_taxa.txt: found {len(output_set)} taxa after exclusion")
    with open(output_filename, "w") as out:
        for taxon in sorted(output_set):
            out.write(taxon + "\n")
    return output_filename

def run_ghostbuster(args, script_dir, ghostbuster_py="ghostbuster.py"):
    ghostbuster_path = os.path.join(script_dir, ghostbuster_py)
    command = [
        "python", ghostbuster_path,
        "--out_taxa", args.out_taxa,
        "--A_taxa", "GhostFinder_A_taxa.txt",
        "--B_taxa", "GhostFinder_B_taxa.txt",
        "--C_taxa", "GhostFinder_C_taxa.txt",
        "--input_trees", args.input_trees,
        "--output_file", "Full_GhostFinder_output.txt",
        "--threads", str(args.threads)
    ]
    print("[DEBUG] Running ghostbuster.py with command:")
    print(" ".join(command))
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running ghostbuster.py: {}\n".format(e))
        sys.exit(1)

def filter_ghostbuster_output(full_output_file, final_output_file):
    try:
        with open(full_output_file, "r") as fin, open(final_output_file, "w") as fout:
            for line in fin:
                # Filter to only keep lines containing the required phrase.
                if "Evidence of sampled introgression" not in line:
                    continue

                # Extract ac_count and bc_count values.
                match_ac = re.search(r"ac_count:\s*(\d+)", line)
                match_bc = re.search(r"bc_count:\s*(\d+)", line)
                if match_ac and match_bc:
                    ac_val = int(match_ac.group(1))
                    bc_val = int(match_bc.group(1))
                    if ac_val > bc_val:
                        fout.write(line)
    except IOError as e:
        sys.stderr.write("Error processing output files: {}\n".format(e))
        sys.exit(1)

def cleanup_intermediate_files(file_list):
    for filename in file_list:
        try:
            if os.path.exists(filename):
                os.remove(filename)
                print(f"[DEBUG] Removed intermediate file: {filename}")
        except Exception as e:
            sys.stderr.write(f"Error removing file {filename}: {e}\n")

def main():
    args = parse_args()

    # Generate intermediate taxa files from ghostbuster_output.
    gf_A = generate_GhostFinder_A_taxa(args.ghostbuster_output, "GhostFinder_A_taxa.txt")
    gf_B = generate_GhostFinder_B_taxa(args.ghostbuster_output, "GhostFinder_B_taxa.txt")
    additional_exclude = generate_additional_taxa_to_exclude(args.ghostbuster_output, "additional_taxa_to_exclude.txt")

    generate_GhostFinder_C_taxa(args.out_taxa, args.input_trees,
                                gf_B_taxa_file=gf_B,
                                additional_exclude_file=additional_exclude,
                                gf_A_taxa_file=gf_A,
                                output_filename="GhostFinder_C_taxa.txt")

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Run ghostbuster.py with the constructed taxa files and threads option.
    run_ghostbuster(args, script_dir)

    # Filter the ghostbuster output for the desired lines.
    filter_ghostbuster_output("Full_GhostFinder_output.txt", args.output_file)
    print(f"[DEBUG] Final output written to {args.output_file}")

    try:
        if os.path.getsize(args.output_file) == 0:
            with open(args.output_file, "w") as f:
                f.write("No putative ghost lineages found.\n")
            print("[DEBUG] Final output was empty; wrote default message.")
    except Exception as e:
        sys.stderr.write(f"Error checking or writing final output file: {e}\n")

    # Cleanup intermediate files.
    cleanup_intermediate_files([
        "GhostFinder_A_taxa.txt",
        "GhostFinder_B_taxa.txt",
        "GhostFinder_C_taxa.txt",
        "additional_taxa_to_exclude.txt",
        "Full_GhostFinder_output.txt"
    ])

if __name__ == "__main__":
    main()
