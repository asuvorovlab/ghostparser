"""Command-line entry point wiring argument parsing and config resolution to the runner."""

import sys

from ..cli_config import EXIT_CHECK_FAILED, run_cli
from .config import build_argument_parser, resolve_config
from .preflight import PreflightResult
from .runner import run_orchestrator


def _run(parsed_args):
    """Resolve the config and run the orchestrator.

    Args:
        parsed_args: The parsed command line.

    Returns:
        ``EXIT_CHECK_FAILED`` when a preflight check finds defects, else
        ``None``.
    """
    result = run_orchestrator(resolve_config(parsed_args))
    if isinstance(result, PreflightResult) and not result.passed:
        return EXIT_CHECK_FAILED
    return None


def main() -> None:
    """Run the orchestrator from the command line and exit with its status."""
    sys.exit(run_cli(build_argument_parser(), _run))


if __name__ == "__main__":
    main()
