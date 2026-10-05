"""Shared CLI/config resolution helpers for GhostParser commands."""

import argparse
import sys

from .config import ConfigError, GhostParserError, _load_raw_config

# Exit statuses shared by every command-line entry point, so a job script can
# tell a finished run from each kind of failure.
EXIT_OK = 0
EXIT_RUN_FAILED = 1
EXIT_CONFIG_ERROR = 2
EXIT_CHECK_FAILED = 3
EXIT_INTERNAL_ERROR = 70
EXIT_INTERRUPTED = 130


def run_cli(parser: argparse.ArgumentParser, run, argv=None) -> int:
    """Parse the command line, call ``run`` and turn its outcome into an exit
    status, reporting a failure as one line unless ``--debug`` is given.

    A :class:`~ghostparser.config.ConfigError` exits with
    ``EXIT_CONFIG_ERROR``, any other
    :class:`~ghostparser.config.GhostParserError` with ``EXIT_RUN_FAILED``, an
    interrupt with ``EXIT_INTERRUPTED``, and anything else, a bug, with
    ``EXIT_INTERNAL_ERROR``.

    Args:
        parser: The entry point's argument parser; ``--debug`` is added to it.
        run: Called with the parsed namespace; returns an exit status, or
            ``None`` for ``EXIT_OK``.
        argv: Arguments to parse instead of ``sys.argv[1:]``.

    Returns:
        The exit status.

    Raises:
        Exception: Whatever ``run`` raised, when ``--debug`` is given, so the
            full traceback is shown.
    """
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show the full traceback when the run fails, instead of one line",
    )
    args = parser.parse_args(argv)
    try:
        status = run(args)
    except KeyboardInterrupt:
        if args.debug:
            raise
        print("Interrupted.", file=sys.stderr)
        return EXIT_INTERRUPTED
    except ConfigError as exc:
        if args.debug:
            raise
        print(f"Config error: {exc}", file=sys.stderr)
        return EXIT_CONFIG_ERROR
    except GhostParserError as exc:
        if args.debug:
            raise
        print(f"Error: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    except Exception as exc:  # noqa: BLE001
        if args.debug:
            raise
        print(
            f"Internal error ({type(exc).__name__}): {exc}\n"
            "This is a bug; rerun with --debug for the full traceback.",
            file=sys.stderr,
        )
        return EXIT_INTERNAL_ERROR
    return EXIT_OK if status is None else status


def _flag_name(arg_name: str) -> str:
    """Turn an argument name into its ``--flag`` spelling."""
    return f"--{arg_name.replace('_', '-')}"


def _given_cli_args(args: argparse.Namespace, arg_names: list[str]) -> dict:
    """Collect the CLI arguments the user gave, by their config key.

    Every flag the modules expose defaults to ``None``, so a value that is not
    ``None`` (``False`` from a ``--no-<x>`` switch included) was given on
    the command line.
    """
    given = {}
    for arg_name in arg_names:
        value = getattr(args, arg_name, None)
        if value is not None:
            given[arg_name] = value
    return given


def resolve_cli_or_config_args(
    args: argparse.Namespace,
    *,
    normalize_payload,
    payload_arg_names: list[str],
) -> argparse.Namespace:
    """Resolve the runtime config from the config file, if any, with every
    given CLI flag laid over it, then normalize the merged payload.

    Args:
        args: The parsed CLI namespace; ``config_file`` names the file, or is
            ``None``.
        normalize_payload: The module's normalizer, taking the raw payload dict.
        payload_arg_names: The CLI argument names that are also config keys.

    Returns:
        The resolved config as a namespace.
    """
    given = _given_cli_args(args, payload_arg_names)
    if "no_overwrite" in given:
        # The switch is the negation of the config key; hand the normalizer
        # the key itself so the switch beats a file's ``overwrite: true``.
        given["overwrite"] = not given.pop("no_overwrite")
    if args.config_file:
        payload = dict(_load_raw_config(args.config_file))
        overridden = sorted(name for name in given if name in payload)
        payload.update(given)
        if overridden:
            print(
                "Note: command-line flags override these config-file settings: "
                + ", ".join(overridden)
            )
    else:
        payload = {name: None for name in payload_arg_names if name != "no_overwrite"}
        payload.update(given)
    return argparse.Namespace(**normalize_payload(payload))
