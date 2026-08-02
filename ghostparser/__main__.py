"""Entry point for the ghostparser package.

This module allows the package to be executed with `python -m ghostparser`.
It prints a usage banner; the runnable work lives in the subpackages below.
"""


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Available modules:")
    print("  pipeline             - Run the streaming introgression pipeline (primary entry point)")
    print("  ml                   - Train multi-label classifiers on summary_statistics.tsv")
    print("  introgression_mapper - Build introgression maps from pipeline results")
    print()
    print("Usage:")
    print("  python -m ghostparser.pipeline -st <species_tree> -gt <gene_trees> -og <outgroup>")
    print("  python -m ghostparser.pipeline -c <config.yaml|config.json>")
    print("  python -m ghostparser.ml --model <random_forest|multi_knn> -i <input.tsv> -o <out>")
    print("  python -m ghostparser.introgression_mapper --help")
    print()
    print("Example:")
    print("  python -m ghostparser.pipeline -st species.nwk -gt genes.nwk -og Outgroup1")
    print("  python -m ghostparser.pipeline -c run_config.yaml")
    print()
    print("For more information, run:")
    print("  python -m ghostparser.pipeline --help")
    print("  python -m ghostparser.ml --help")


if __name__ == "__main__":
    main()
