"""Package entry point for Ghostparser ML.

Run the individual trainer modules directly instead, for example
``python -m ghostparser.ml.random_forest`` or
``python -m ghostparser.ml.multi_knn``.
"""

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ghostparser.ml package entrypoint. Run a trainer module directly."
    )
    parser.parse_args()
    print("Ghostparser ML")
    print()
    print("This package entrypoint does not dispatch to model modules.")
    print("Run one of:")
    print("  python -m ghostparser.ml.random_forest -i <input.tsv> -o <out>")
    print("  python -m ghostparser.ml.multi_knn -i <input.tsv> -o <out>")
    print("  python -m ghostparser.ml.hyper_tune -c <config.yaml|config.json>")


if __name__ == "__main__":
    main()
