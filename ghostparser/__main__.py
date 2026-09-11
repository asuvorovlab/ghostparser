"""Usage banner for ``python -m ghostparser``; the work lives in the subpackages."""


def main():
    """Display usage information for the ghostparser package."""
    print("GhostParser")
    print()
    print("Subpackages:")
    print(
        "  orchestrator - Run the introgression inference pipeline (primary entry point).\n"
        "                 Consolidation maps are produced as its final stage."
    )
    print(
        "  ml           - Train classifiers on the pipeline's summary_statistics.tsv\n"
        "                 and tune their hyperparameters."
    )
    print()
    print("Usage:")
    print(
        "  python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>"
    )
    print("  python -m ghostparser.orchestrator -c <config.yaml|config.json>")
    print("  python -m ghostparser.ml.random_forest -i <summary_statistics.tsv> -o <out>")
    print("  python -m ghostparser.ml.multi_knn -i <summary_statistics.tsv> -o <out>")
    print("  python -m ghostparser.ml.hyper_tune -c <config.yaml|config.json>")
    print()
    print("Example:")
    print(
        "  python -m ghostparser.orchestrator -st species.nwk -gt genes.nwk -og Outgroup1"
    )
    print("  python -m ghostparser.orchestrator -c run_config.yaml")
    print(
        "  python -m ghostparser.ml.random_forest -i results/summary_statistics.tsv -o results/ml_out"
    )
    print()
    print("For more information, run:")
    print("  python -m ghostparser.orchestrator --help")
    print("  python -m ghostparser.ml")


if __name__ == "__main__":
    main()
