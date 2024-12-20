import argparse
from ete3 import Tree

# Set up argument parsing
parser = argparse.ArgumentParser(description="Rescale branch lengths of a phylogenetic tree.")
parser.add_argument("tree", type=str, help="Path to the input tree file")
parser.add_argument("outfile", type=str, help="Path to save the rescaled tree")

# Parse the arguments
args = parser.parse_args()

# Load the tree from the provided file
tree = Tree(args.tree)

# Define the mutation rate (e.g., 1e-9 substitutions per generation)
mutation_rate = 1e-9

# Rescale branch lengths
for node in tree.traverse():
    node.dist *= mutation_rate

# Save the rescaled tree to the specified output file
tree.write(outfile=args.outfile)

