"""Entry point for the ghostparser package.

This module allows the package to be executed with `python -m ghostparser`.
It prints a usage banner; the runnable work lives in the subpackages below.
"""


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Available modules:")
    print(
        "  orchestrator         - Run the introgression inference pipeline (primary entry point)"
    )
    print(
        "  ml                   - Show ML trainer usage; run trainer modules directly on summary_statistics.tsv"
    )
    print(
        "  introgression_mapper - Build introgression maps from orchestrator results"
    )
    print()
    print("Usage:")
    print(
        "  python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>"
    )
    print("  python -m ghostparser.orchestrator -c <config.yaml|config.json>")
    print("  python -m ghostparser.ml.random_forest -i <input.tsv> -o <out>")
    print("  python -m ghostparser.ml.multi_knn -i <input.tsv> -o <out>")
    print("  python -m ghostparser.introgression_mapper --help")
    print()
    print("Example:")
    print(
        "  python -m ghostparser.orchestrator -st species.nwk -gt genes.nwk -og Outgroup1"
    )
    print("  python -m ghostparser.orchestrator -c run_config.yaml")
    print()
    print("For more information, run:")
    print("  python -m ghostparser.orchestrator --help")
    print("  python -m ghostparser.ml")


if __name__ == "__main__":
    main()
