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
            # Write header to TSV file
            tsv_writer.writerow(["Topology", "Distance", "Out", "Internal"])

            for line_number, line in enumerate(file, start=1):
                line = line.strip()
                if line:
                    tree = Tree(line, format=1)
                    # Root the tree at "Out"
                    tree.set_outgroup("Out")
                    # Prune "Ghost" from the tree
                    tree.prune([node.name for node in tree if node.name != "Ghost"], preserve_branch_length=True)
                    
                    # Determine sister relationships among A, B, C and calculate distances
                    a_node = tree & "A"
                    b_node = tree & "B"
                    c_node = tree & "C"
                    
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
                    
                    distance_to_mrca = "None"
                    internal_branch_length = "None"
                    
                    # Determine which of A, B, C is not sister and calculate distance to MRCA
                    if len(sister_pairs) == 1:
                        non_sister = None
                        if "AB" in sister_pairs:
                            non_sister = c_node
                        elif "AC" in sister_pairs:
                            non_sister = b_node
                        elif "BC" in sister_pairs:
                            non_sister = a_node
                        
                        if non_sister:
                            mrca_abc = tree.get_common_ancestor(a_node, b_node, c_node)
                            distance_to_mrca = non_sister.get_distance(mrca_abc)
                            
                            # Calculate internal branch length from MRCA of A, B, C to MRCA of sister pair
                            mrca_sister_pair = tree.get_common_ancestor(a_node, b_node) if "AB" in sister_pairs else (
                                tree.get_common_ancestor(a_node, c_node) if "AC" in sister_pairs else tree.get_common_ancestor(b_node, c_node))
                            internal_branch_length = mrca_abc.get_distance(mrca_sister_pair)
                    
                    # Write data to TSV file
                    for pair, distance in distances.items():
                        tsv_writer.writerow([pair, distance, distance_to_mrca, internal_branch_length])
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()
