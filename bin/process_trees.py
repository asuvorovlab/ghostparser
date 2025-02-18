import argparse
import os
from ete3 import Tree
import csv

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process and modify trees in a Newick file.")
    parser.add_argument("file_path", help="Path to the file to be processed")
    parser.add_argument("output_path", help="Path to the output TSV file")

    # Parse arguments
    args = parser.parse_args()

    # Get file paths from arguments
    file_path = args.file_path
    output_path = args.output_path

    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        return

    # Open and process the file
    try:
        with open(file_path, 'r') as file, open(output_path, 'w', newline='') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t')
            # Updated header to include the new A_length column
            tsv_writer.writerow(["Topology", "Distance", "A_length", "Out", "Internal", "Internal_ABCO"])

            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                # Skip lines that do not contain all required taxa
                required_taxa = {"A", "B", "C", "Out"}
                if not all(taxon in line for taxon in required_taxa):
                    continue

                if line:
                    tree = Tree(line, format=1)
                    # Root the tree at "Out"
                    tree.set_outgroup("Out")
                    # Prune "Ghost" from the tree if present
                    tree.prune([node.name for node in tree if node.name != "Ghost"], preserve_branch_length=True)

                    # Get nodes for the taxa
                    a_node = tree & "A"
                    b_node = tree & "B"
                    c_node = tree & "C"

                    # Determine sister relationships among A, B, and C and calculate distances
                    sister_pairs = []
                    distances = {}
                    if a_node.up == b_node.up:
                        sister_pairs.append("AB")
                        distances["AB"] = a_node.get_distance(b_node) / 2
                    if a_node.up == c_node.up:
                        sister_pairs.append("AC")
                        distances["AC"] = a_node.get_distance(c_node) / 2
                    if b_node.up == c_node.up:
                        sister_pairs.append("BC")
                        distances["BC"] = b_node.get_distance(c_node) / 2

                    # Initialize values for additional metrics
                    distance_to_mrca = "None"
                    internal_branch_length = "None"
                    internal_abco = "None"
                    a_length = "None"  # New metric: distance from A to MRCA(A, B, C)

                    # Calculate distances if only one sister pair is found
                    if len(sister_pairs) == 1:
                        non_sister = None
                        if "AB" in sister_pairs:
                            non_sister = c_node
                        elif "AC" in sister_pairs:
                            non_sister = b_node
                        elif "BC" in sister_pairs:
                            non_sister = a_node

                        if non_sister:
                            # MRCA of A, B, and C
                            mrca_abc = tree.get_common_ancestor(a_node, b_node, c_node)
                            # New A_length: distance from A to the MRCA of A, B, and C
                            a_length = a_node.get_distance(mrca_abc)
                            distance_to_mrca = non_sister.get_distance(mrca_abc)

                            # Calculate internal branch length from MRCA of A, B, C to MRCA of the sister pair
                            if "AB" in sister_pairs:
                                mrca_sister_pair = tree.get_common_ancestor(a_node, b_node)
                            elif "AC" in sister_pairs:
                                mrca_sister_pair = tree.get_common_ancestor(a_node, c_node)
                            else:
                                mrca_sister_pair = tree.get_common_ancestor(b_node, c_node)
                            internal_branch_length = mrca_abc.get_distance(mrca_sister_pair)

                            # Calculate internal branch length from MRCA of A, B, C to MRCA of all four taxa (A, B, C, Out)
                            mrca_abco = tree.get_common_ancestor(a_node, b_node, c_node, tree & "Out")
                            internal_abco = mrca_abco.get_distance(mrca_abc)

                    # Write data to TSV file; each row corresponds to one sister pair
                    for pair, distance in distances.items():
                        tsv_writer.writerow([pair, distance, a_length, distance_to_mrca, internal_branch_length, internal_abco])
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()

