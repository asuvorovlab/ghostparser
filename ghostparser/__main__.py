"""Usage banner for ``python -m ghostparser``; the work lives in the subpackages."""


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Available modules:")
    print(
        "  orchestrator - Run the introgression inference pipeline (primary entry point).\n"
        "                 Consolidation maps are produced as its final stage."
    )
    print(
        "  ml           - Show ML trainer usage; run trainer modules directly on summary_statistics.tsv"
    )
    print()
    print("Usage:")
    print(
        "  python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>"
    )
    print("  python -m ghostparser.orchestrator -c <config.yaml|config.json>")
    print("  python -m ghostparser.ml.random_forest -i <input.tsv> -o <out>")
    print("  python -m ghostparser.ml.multi_knn -i <input.tsv> -o <out>")
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
