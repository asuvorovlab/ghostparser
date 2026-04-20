"""Entry point for the ghostparser package.

This module allows the package to be executed with `python -m ghostparser`.
For tree parsing functionality, use `python -m ghostparser.tree_parser`.
"""

import sys


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Available modules:")
    print("  tree_parser    - Parse and standardize Newick format trees")
    print("  triplet_processor - Run Fig. 6 triplet introgression classification pipeline")
    print("  orchestrator   - Run end-to-end tree parsing + per-triplet inference pipeline")
    print()
    print("Usage:")
    print("  python -m ghostparser.tree_parser -st <species_tree> -gt <gene_trees> -og <outgroup>")
    print("  python -m ghostparser.triplet_processor -i <unique_triplets_gene_trees.txt>")
    print("  python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>")
    print("  python -m ghostparser.orchestrator -c <config.yaml|config.json>")
    print("  (Only orchestrator supports --config-file)")
    print()
    print("Example:")
    print("  python -m ghostparser.tree_parser -st species.nwk -gt genes.nwk -og Outgroup1")
    print("  python -m ghostparser.triplet_processor -i unique_triplets_gene_trees.txt")
    print("  python -m ghostparser.orchestrator -st species.nwk -gt genes.nwk -og Outgroup1")
    print("  python -m ghostparser.orchestrator -c run_config.yaml")
    print()
    print("For more information, run:")
    print("  python -m ghostparser.tree_parser --help")
    print("  python -m ghostparser.triplet_processor --help")
    print("  python -m ghostparser.orchestrator --help")


if __name__ == "__main__":
    main()
