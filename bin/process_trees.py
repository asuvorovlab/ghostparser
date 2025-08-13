import argparse
import os
from ete3 import Tree
import csv
import sys

def main():
    parser = argparse.ArgumentParser(description="Process and modify trees in a Newick file.")
    parser.add_argument("file_path", help="Path to the file to be processed")
    parser.add_argument("output_path", help="Path to the output TSV file")
    args = parser.parse_args()

    file_path = args.file_path
    output_path = args.output_path

    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' does not exist.", file=sys.stderr)
        # create an empty output with header so caller can detect "no data"
        with open(output_path, 'w', newline='') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t')
            tsv_writer.writerow(["Topology", "tree_height"])
        return

    try:
        with open(file_path, 'r') as infile, open(output_path, 'w', newline='') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t')
            tsv_writer.writerow(["Topology", "tree_height"])

            for line_number, line in enumerate(infile, start=1):
                s = line.strip()
                if not s:
                    continue

                # Require all 4 taxa
                required_taxa = {"A", "B", "C", "Out"}
                if not all(("(" + t + ":" in s) or ("," + t + ":" in s) or (")" + t + ":" in s) or (t + ":" in s) for t in required_taxa):
                    # fall back to substring if branch lengths missing; still skip if absent
                    if not all(t in s for t in required_taxa):
                        continue

                try:
                    tree = Tree(s, format=1)
                except Exception as e:
                    # bad line; skip
                    continue

                # Root at Out if possible
                try:
                    tree.set_outgroup("Out")
                except Exception:
                    # If Out is missing or rooting fails, skip
                    continue

                # Get nodes
                try:
                    a_node = tree & "A"
                    b_node = tree & "B"
                    c_node = tree & "C"
                except Exception:
                    continue

                sister_pairs = []
                try:
                    if a_node.up == b_node.up:
                        sister_pairs.append("AB")
                    if a_node.up == c_node.up:
                        sister_pairs.append("AC")
                    if b_node.up == c_node.up:
                        sister_pairs.append("BC")
                except Exception:
                    continue

                if len(sister_pairs) != 1:
                    # ambiguous or unresolved; skip
                    continue

                # compute distances safely
                pair = sister_pairs[0]
                if pair == "AB":
                    sis1, sis2, out = a_node, b_node, c_node
                elif pair == "AC":
                    sis1, sis2, out = a_node, c_node, b_node
                else:
                    sis1, sis2, out = b_node, c_node, a_node

                try:
                    mrca_sister = tree.get_common_ancestor(sis1, sis2)
                    mrca_abc = tree.get_common_ancestor(a_node, b_node, c_node)

                    b1 = sis1.get_distance(mrca_sister)
                    b2 = sis2.get_distance(mrca_sister)
                    b3 = mrca_sister.get_distance(mrca_abc)
                    b4 = out.get_distance(mrca_abc)

                    tree_height = (b3 + b4 + (b1 + b2) / 2.0) / 2.0
                except Exception:
                    # numeric failure; skip
                    continue

                tsv_writer.writerow([pair, f"{tree_height:.10f}"])

    except Exception as e:
        # Last-ditch: emit header-only file so caller can detect emptiness
        with open(output_path, 'w', newline='') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t')
            tsv_writer.writerow(["Topology", "tree_height"])
        print(f"Error processing file: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()

