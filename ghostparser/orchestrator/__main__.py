"""Command-line entry point wiring argument parsing and config resolution to the runner."""

from .config import ConfigError, build_argument_parser, resolve_config
from .runner import run_orchestrator


def main() -> None:
    """Parse the command line, resolve the config and run the orchestrator."""
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
