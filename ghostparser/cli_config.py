"""Shared CLI/config resolution helpers for GhostParser commands."""

import argparse


def _flag_name(arg_name: str) -> str:
    return f"--{arg_name.replace('_', '-')}"


def _provided_cli_args(args: argparse.Namespace, arg_names: list[str]) -> list[str]:
    provided: list[str] = []
    for arg_name in arg_names:
        value = getattr(args, arg_name, None)
        if isinstance(value, bool):
            if value:
                provided.append(_flag_name(arg_name))
            continue
        if value is not None:
            provided.append(_flag_name(arg_name))
    return provided


def resolve_cli_or_config_args(
    args: argparse.Namespace,
    *,
    load_config,
    normalize_payload,
    payload_arg_names: list[str],
) -> argparse.Namespace:
    """Resolve runtime args from either config-file mode or plain CLI mode."""
    if args.config_file:
        ignored_cli_args = _provided_cli_args(args, payload_arg_names)
        if ignored_cli_args:
            print(
                "Warning: --config-file provided; CLI arguments not in config will be ignored: "
                + ", ".join(ignored_cli_args)
            )
        config = load_config(args.config_file)
        return argparse.Namespace(**config)

    payload = {name: getattr(args, name, None) for name in payload_arg_names}
    config = normalize_payload(payload)
    return argparse.Namespace(**config)
