"""Package entry point for Ghostparser ML.

Run with ``python -m ghostparser.ml --model random_forest [args...]`` to dispatch to
the requested model module. By default (no --model) this entrypoint will do nothing
and will print available models to avoid implicit redirection.
"""

from __future__ import annotations

import argparse
import importlib
import sys


def _available_models() -> list[str]:
    # Keep this list in sync with available modules under ghostparser.ml
    return ["multi_knn", "random_forest"]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ghostparser.ml package entrypoint. Use --model to dispatch to a model module."
    )
    parser.add_argument("--model", choices=_available_models(), help="Model module to run (explicit)")
    # parse only the --model flag, leave the rest for the selected module
    args, remaining = parser.parse_known_args()

    if args.model is None:
        print("No model requested. Available models: {}".format(", ".join(_available_models())))
        print("Run with: python -m ghostparser.ml --model random_forest [model args]")
        return

    # Rebuild argv so the selected module's main() sees only its own args
    sys.argv = [sys.argv[0]] + remaining

    module_name = f"ghostparser.ml.{args.model}"
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:  # import errors should surface to the user
        print(f"Failed to import model module '{module_name}': {exc}")
        raise

    # Call the module's main entrypoint
    if not hasattr(module, "main"):
        raise RuntimeError(f"Model module {module_name} does not expose a 'main' function")
    module.main()


if __name__ == "__main__":
    main()