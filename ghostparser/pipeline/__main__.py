"""Command-line entry point for ``python -m ghostparser.pipeline``.

Wires argument parsing and config resolution (from ``config``) to the pipeline
runner, keeping the CLI orchestration separate from config handling.
"""

from __future__ import annotations

from .config import ConfigError, build_argument_parser, resolve_config
from .runner import run_pipeline


def main() -> None:
    """Run the streaming pipeline from the command line.

    Parses process argv, resolves the config, and invokes
    :func:`ghostparser.pipeline.runner.run_pipeline`. On a ``ConfigError`` the
    parser prints the error and exits.
    """
    parser = build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        config = resolve_config(parsed_args)
    except ConfigError as exc:
        parser.error(str(exc))
        return

    run_pipeline(config)


if __name__ == "__main__":
    main()
