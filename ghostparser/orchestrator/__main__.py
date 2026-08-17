"""Command-line entry point for ``python -m ghostparser.orchestrator``.

Wires argument parsing and config resolution (from ``config``) to the orchestrator
runner, keeping the CLI orchestration separate from config handling.
"""

from .config import ConfigError, build_argument_parser, resolve_config
from .runner import run_orchestrator


def main() -> None:
    """Run the orchestrator from the command line.

    Parses process argv, resolves the config, and invokes
    :func:`ghostparser.orchestrator.runner.run_orchestrator`. On a ``ConfigError`` the
    parser prints the error and exits.
    """
    parser = build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        config = resolve_config(parsed_args)
    except ConfigError as exc:
        parser.error(str(exc))
        return

    run_orchestrator(config)


if __name__ == "__main__":
    main()
