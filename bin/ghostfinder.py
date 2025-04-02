#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Ghostbuster wrapper that runs ghostbuster.py in the current directory."
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
                        help="Number of parallel threads (default: number of CPU cores).")
    return parser.parse_args()

def generate_taxa_files(ghostbuster_output):
    taxaA = set()
    taxaB = set()
    additional_exclude = set()

    with open(ghostbuster_output, "r") as fin:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            # Compare after stripping extra spaces.
            if parts[1].strip() == "Evidence of unsampled introgression":
                fields = parts[0].split()
                if len(fields) >= 2:
                    taxaA.add(fields[1])
                if len(fields) >= 3:
                    taxaB.add(fields[2])
                if len(fields) >= 1:
                    additional_exclude.add(fields[0])
    with open("GhostFinder_A_taxa.txt", "w") as fout:
        for t in sorted(taxaA):
            fout.write(t + "\n")
    with open("GhostFinder_B_taxa.txt", "w") as fout:
        for t in sorted(taxaB):
            fout.write(t + "\n")
    with open("additional_taxa_to_exclude.txt", "w") as fout:
        for t in sorted(additional_exclude):
            fout.write(t + "\n")

def generate_GhostFinder_C_taxa(out_taxa, input_trees):
    ex = set()
    for fname in ["GhostFinder_B_taxa.txt", "additional_taxa_to_exclude.txt", "GhostFinder_A_taxa.txt"]:
        if os.path.exists(fname):
            with open(fname, "r") as fin:
                for line in fin:
                    token = line.strip().split()[0] if line.strip() else ""
                    if token:
                        ex.add(token)
    if os.path.basename(out_taxa) == "Out.txt" and os.path.exists(out_taxa):
        with open(out_taxa, "r") as fin:
            for line in fin:
                token = line.strip().split()[0] if line.strip() else ""
                if token:
                    ex.add(token)
    C_taxa = set()
    for fname in [out_taxa, input_trees]:
        if os.path.exists(fname):
            with open(fname, "r") as fin:
                for line in fin:
                    newline = re.sub(r"[(),;]", " ", line)
                    tokens = newline.strip().split()
                    for token in tokens:
                        token = token.split(":", 1)[0]
                        if token and token not in ex and not re.fullmatch(r"[0-9]+", token):
                            C_taxa.add(token)
    with open("GhostFinder_C_taxa.txt", "w") as fout:
        for t in sorted(C_taxa):
            fout.write(t + "\n")

def run_ghostbuster_instance(instance_id, args, ghostbuster_py_path,
                             taxa_A_path, taxa_B_path, taxa_C_path):
    """
    Runs ghostbuster.py in the current directory using the provided input files.
    Reads the output file (Full_GhostFinder_output.txt) and returns the filtered output.
    Note: The output file is NOT removed.
    """
    cmd = [
        sys.executable, ghostbuster_py_path,
        "--out_taxa", os.path.abspath(args.out_taxa),
        "--A_taxa", os.path.abspath(taxa_A_path),
        "--B_taxa", os.path.abspath(taxa_B_path),
        "--C_taxa", os.path.abspath(taxa_C_path),
        "--input_trees", os.path.abspath(args.input_trees),
        "--output_file", "Full_GhostFinder_output.txt"
    ]
    sys.stderr.write(f"[Instance {instance_id}] Running command: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    curr_listing = os.listdir(os.getcwd())
    sys.stderr.write(f"[Instance {instance_id}] Current directory files: {curr_listing}\n")
    
    output_path = os.path.join(os.getcwd(), "Full_GhostFinder_output.txt")
    content = ""
    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        with open(output_path, "r") as fin:
            content = fin.read()
        sys.stderr.write(f"[Instance {instance_id}] Read {len(content)} bytes from output file.\n")
        sys.stderr.write(f"[Instance {instance_id}] Content of output file:\n{content}\n")
    else:
        content = result.stdout.decode()
        sys.stderr.write(f"[Instance {instance_id}] Output file empty; using stdout ({len(content)} bytes).\n")
        sys.stderr.write(f"[Instance {instance_id}] Content from stdout:\n{content}\n")
    
    # Do NOT remove the Full_GhostFinder_output.txt file.
    
    # Filter lines containing " sampled"
    filtered_lines = "\n".join([line for line in content.splitlines() if " sampled" in line])
    if not filtered_lines:
        sys.stderr.write(f"[Instance {instance_id}] No lines with ' sampled' found in the output.\n")
    else:
        sys.stderr.write(f"[Instance {instance_id}] Found {len(filtered_lines.splitlines())} matching lines.\n")
    return filtered_lines + "\n" if filtered_lines else ""

def main():
    args = parse_args()

    # Check that required files exist.
    for f in [args.out_taxa, args.ghostbuster_output, args.input_trees]:
        if not os.path.exists(f):
            sys.stderr.write(f"Error: file {f} does not exist.\n")
            sys.exit(1)

    # Generate the taxa files.
    generate_taxa_files(args.ghostbuster_output)
    generate_GhostFinder_C_taxa(args.out_taxa, args.input_trees)

    # Determine ghostbuster.py path (assumed to be in the same directory as this script).
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ghostbuster_py_path = os.path.join(script_dir, "ghostbuster.py")
    if not os.path.exists(ghostbuster_py_path):
        sys.stderr.write(f"Error: ghostbuster.py not found at {ghostbuster_py_path}\n")
        sys.exit(1)

    # Use absolute paths for the taxa files.
    taxa_A_path = os.path.abspath("GhostFinder_A_taxa.txt")
    taxa_B_path = os.path.abspath("GhostFinder_B_taxa.txt")
    taxa_C_path = os.path.abspath("GhostFinder_C_taxa.txt")

    # Run ghostbuster.py without parallelization.
    result = run_ghostbuster_instance(0, args, ghostbuster_py_path,
                                      taxa_A_path, taxa_B_path, taxa_C_path)
    with open(args.output_file, "w") as fout:
        fout.write(result)
    print(f"Combined output written to {args.output_file}")

if __name__ == "__main__":
    main()
