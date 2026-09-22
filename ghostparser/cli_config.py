"""Shared CLI/config resolution helpers for GhostParser commands."""

import argparse

from .config import _load_raw_config


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
