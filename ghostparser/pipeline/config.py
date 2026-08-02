"""Configuration for the streaming pipeline.

Owns the pipeline's own defaults/choices, its validation rules, the CLI parser,
and the CLI/config resolution. The helpers whose behaviour is shared verbatim
with the ML subpackage — ``ConfigError``, path resolution, raw config-file
loading, required-path validation, overwrite resolution, and output-directory
preparation — are imported from the :mod:`ghostparser.config` trunk.

Config-file mode: ``-c/--config-file`` is CLI-only, and when a config file is
given the other CLI flags are ignored with a warning (config wins). A subset of
runtime knobs is exposed both on the CLI and in the config file; the rest are
config-file-only.
"""

from __future__ import annotations

import argparse

from ..cli_config import resolve_cli_or_config_args
from ..config import (
    DEFAULT_OVERWRITE,
    ConfigError,
    _load_raw_config,
    _resolve_path,
    _validate_overwrite_flag,
    _validate_required_path,
    prepare_output_directory,
)

__all__ = ["ConfigError", "prepare_output_directory", "build_argument_parser",
           "resolve_config", "load_pipeline_config", "normalize_pipeline_payload"]

# Runtime defaults specific to the pipeline.
DEFAULT_OUTPUT_FOLDER = "results"
DEFAULT_PROCESSES = 0
DEFAULT_MIN_SUPPORT_VALUE = 0.5
DEFAULT_DISCORDANT_TEST = "chi-square"
DEFAULT_SUMMARY_STATISTIC = "mean"
DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY = "AVG"
DEFAULT_P_VALUE_CORRECTION = "bfn"
DEFAULT_ALPHA_DCT = 0.05
DEFAULT_ALPHA_KS = 0.05
DEFAULT_BOOTSTRAP = True
DEFAULT_BOOTSTRAP_ITERATIONS = 100
DEFAULT_BOOTSTRAP_DEBUG_MODE = False
DEFAULT_BOOTSTRAP_SUMMARY_ONLY = False
DEFAULT_GENERATE_SUMMARY_STATS = False
DEFAULT_CONSOLIDATION = True

DISCORDANT_TEST_CHOICES = ("chi-square", "z-test")
SUMMARY_STATISTIC_CHOICES = ("mean", "median", "mode")
TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES = ("AVG", "A", "B", "C", "SIS", "INT")
P_VALUE_CORRECTION_CHOICES = ("no", "bfn", "holm", "fdr_bh", "fdr_by", "fdr_tsbh")

# Pipeline-specific parallelization knobs (no orchestrator equivalent to reuse).
PARALLELIZATION_MODE_CHOICES = ("auto", "taxon", "gene")
DEFAULT_PARALLELIZATION_MODE = "auto"
# `auto` picks `gene` for small-taxa/large-gene-tree runs, else `taxon`.
AUTO_TAXA_SMALL_THRESHOLD = 15
AUTO_GENE_TREES_THRESHOLD = 3500

# CLI argument dest names that also map to config-file payload keys. These are
# the config+CLI options; config-file-only keys (discordant_test,
# tree_height_calculation_strategy, min_support_value, bootstrap_iterations,
# bootstrap_seed, generate_summary_stats, bootstrap_debug_mode,
# bootstrap_summary_only) are intentionally absent so they are read only from a
# config file and otherwise take their defaults.
_PIPELINE_PAYLOAD_ARG_NAMES = [
    "species_tree_path",
    "gene_trees_path",
    "outgroups",
    "output_folder",
    "triplet_filter",
    "no_overwrite",
    "processes",
    "parallelization_mode",
    "alpha_dct",
    "alpha_ks",
    "p_value_correction",
    "summary_statistic",
    "consolidation",
    "bootstrap",
]


def _validate_optional_path(payload: dict, key: str) -> str | None:
    """Validate and resolve an optional path field.

    Args:
        payload: The config/CLI payload.
        key: The field name.

    Returns:
        The resolved absolute path, or ``None`` when the field is absent.

    Raises:
        ConfigError: If the field is present but not a non-empty string.
    """
    value = payload.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(
            f"Config field {key} must be a non-empty string when provided"
        )
    return _resolve_path(value.strip())


def _validate_non_negative_int(payload: dict, key: str, default: int) -> int:
    """Validate a non-negative integer field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated integer.

    Raises:
        ConfigError: If the value is not an integer >= 0.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, int) or value < 0:
        raise ConfigError(f"Config field {key} must be an integer >= 0")
    return value


def _validate_optional_bool(payload: dict, key: str, default: bool) -> bool:
    """Validate an optional boolean field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated boolean.

    Raises:
        ConfigError: If the value is present but not a boolean.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, bool):
        raise ConfigError(f"Config field {key} must be a boolean when provided")
    return value


def _validate_optional_float(payload: dict, key: str, default: float) -> float:
    """Validate an optional numeric field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated float.

    Raises:
        ConfigError: If the value cannot be coerced to a float.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"Config field {key} must be a numeric value") from exc


def _validate_choice(
    payload: dict, key: str, default: str, choices: tuple[str, ...]
) -> str:
    """Validate an enumerated string field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.
        choices: The allowed values.

    Returns:
        The validated choice.

    Raises:
        ConfigError: If the value is not one of ``choices``.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, str) or value not in choices:
        raise ConfigError(f"Config field {key} must be one of: {', '.join(choices)}")
    return value


def _validate_bootstrap_options(payload: dict) -> tuple[bool, dict]:
    """Validate the bootstrap toggle and its nested/flat options.

    Accepts the canonical ``bootstrap`` toggle plus either a nested
    ``bootstrap_options`` block or the flat ``bootstrap_iterations``/
    ``bootstrap_seed``/``bootstrap_debug_mode``/``bootstrap_summary_only`` keys.

    Args:
        payload: The config/CLI payload.

    Returns:
        A tuple ``(bootstrap, options)`` where ``options`` has
        ``iterations``/``seed``/``debug_mode``/``summary_only``.

    Raises:
        ConfigError: If any value is malformed.
    """
    bootstrap = _validate_optional_bool(payload, "bootstrap", DEFAULT_BOOTSTRAP)

    raw_options = payload.get("bootstrap_options")
    if raw_options is None:
        raw_options = {}
    if not isinstance(raw_options, dict):
        raise ConfigError(
            "Config field bootstrap_options must be a key/value object when provided"
        )

    iterations = payload.get(
        "bootstrap_iterations",
        raw_options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS),
    )
    if iterations is None:
        iterations = DEFAULT_BOOTSTRAP_ITERATIONS
    if not isinstance(iterations, int) or iterations < 1:
        raise ConfigError(
            "Config field bootstrap_options.iterations must be an integer >= 1"
        )

    seed = payload.get("bootstrap_seed", raw_options.get("seed"))
    if seed is not None and not isinstance(seed, int):
        raise ConfigError(
            "Config field bootstrap_options.seed must be an integer when provided"
        )

    debug_mode = payload.get(
        "bootstrap_debug_mode",
        raw_options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE),
    )
    if debug_mode is None:
        debug_mode = DEFAULT_BOOTSTRAP_DEBUG_MODE
    if not isinstance(debug_mode, bool):
        raise ConfigError(
            "Config field bootstrap_options.debug_mode must be a boolean when provided"
        )

    summary_only = payload.get(
        "bootstrap_summary_only",
        raw_options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY),
    )
    if summary_only is None:
        summary_only = DEFAULT_BOOTSTRAP_SUMMARY_ONLY
    if not isinstance(summary_only, bool):
        raise ConfigError(
            "Config field bootstrap_options.summary_only must be a boolean when provided"
        )

    return bootstrap, {
        "iterations": iterations,
        "seed": seed,
        "debug_mode": debug_mode,
        "summary_only": summary_only,
    }


def _parse_outgroups(value) -> list[str]:
    """Parse the outgroup(s) value into a list of taxon labels.

    Accepts a comma-separated string or a list/tuple/set of labels.

    Args:
        value: The raw outgroup(s) value.

    Returns:
        A non-empty list of outgroup labels.

    Raises:
        ConfigError: If no outgroup labels can be parsed.
    """
    if isinstance(value, str):
        outgroups = [part.strip() for part in value.split(",") if part.strip()]
    elif isinstance(value, (list, tuple, set)):
        outgroups = [str(part).strip() for part in value if str(part).strip()]
    else:
        outgroups = []

    if not outgroups:
        raise ConfigError("Missing required config field: outgroup(s)")
    return outgroups


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the streaming-pipeline CLI parser.

    Returns:
        A configured ``argparse.ArgumentParser`` for the pipeline entry point.
    """
    parser = argparse.ArgumentParser(
        prog="python -m ghostparser.pipeline",
        description=(
            "GhostParser streaming pipeline: fuse triplet extraction and inference "
            "in a single pass without materializing the intermediate triplet file."
        ),
    )
    parser.add_argument(
        "-c",
        "--config-file",
        default=None,
        help="Path to a JSON or YAML config file (config-file mode; other CLI flags are ignored)",
    )
    parser.add_argument(
        "-st",
        "--species-tree-path",
        default=None,
        help="Path to the species tree file in Newick format",
    )
    parser.add_argument(
        "-gt",
        "--gene-trees-path",
        default=None,
        help="Path to the gene trees file in Newick format",
    )
    parser.add_argument(
        "-og",
        "--outgroups",
        default=None,
        help="Outgroup species identifier(s), comma-separated",
    )
    parser.add_argument(
        "--output-folder",
        default=None,
        help=f"Output folder (default: ./{DEFAULT_OUTPUT_FOLDER})",
    )
    parser.add_argument(
        "--triplet-filter",
        default=None,
        help="Path to triplet filter file (comma-separated taxa per line)",
    )
    parser.add_argument(
        "--no-overwrite",
        dest="no_overwrite",
        action="store_true",
        default=None,
        help="Append a numeric suffix when the output folder already exists",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes (0 = all cores)",
    )
    parser.add_argument(
        "--parallelization-mode",
        choices=PARALLELIZATION_MODE_CHOICES,
        default=None,
        help=(
            "Parallelization strategy: 'taxon' dispatches triplet chunks across "
            "workers, 'gene' parallelizes subtree extraction within a triplet, "
            f"'auto' selects one from the input size (default: {DEFAULT_PARALLELIZATION_MODE})"
        ),
    )
    parser.add_argument(
        "--alpha-dct",
        type=float,
        default=None,
        help=f"DCT significance threshold (default: {DEFAULT_ALPHA_DCT})",
    )
    parser.add_argument(
        "--alpha-ks",
        type=float,
        default=None,
        help=f"KS significance threshold (default: {DEFAULT_ALPHA_KS})",
    )
    parser.add_argument(
        "--p-value-correction",
        choices=P_VALUE_CORRECTION_CHOICES,
        default=None,
        help=f"Multiple-testing correction for triplet p-values (default: {DEFAULT_P_VALUE_CORRECTION})",
    )
    parser.add_argument(
        "--summary-statistic",
        choices=SUMMARY_STATISTIC_CHOICES,
        default=None,
        help=f"Summary statistic after the KS test (default: {DEFAULT_SUMMARY_STATISTIC})",
    )
    parser.add_argument(
        "--no-consolidation",
        dest="consolidation",
        action="store_false",
        default=None,
        help="Disable introgression map consolidation outputs (default: enabled)",
    )
    parser.add_argument(
        "--no-bootstrap",
        dest="bootstrap",
        action="store_false",
        default=None,
        help="Disable bootstrap sampling-with-replacement during triplet inference",
    )
    return parser


def normalize_pipeline_payload(payload: dict) -> dict:
    """Normalize a pipeline config/CLI payload into the runtime config dict.

    Validates and fills every runtime key from the payload (a parsed config file
    or a CLI-derived dict), applying pipeline defaults for anything omitted.
    Config-file-only keys default when absent, which is what happens in CLI mode
    since they have no corresponding flag.

    Args:
        payload: Raw config/CLI key-value mapping.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.pipeline.runner.run_pipeline`.

    Raises:
        ConfigError: If a required field is missing or a value is invalid.
    """
    species_tree = _validate_required_path(payload, "species_tree_path")
    gene_trees = _validate_required_path(payload, "gene_trees_path")

    outgroups_source = payload.get("outgroups")
    if outgroups_source is None:
        outgroups_source = payload.get("outgroup")
    outgroups = _parse_outgroups(outgroups_source)

    output = _validate_optional_path(payload, "output_folder")
    if output is None:
        output = _resolve_path(DEFAULT_OUTPUT_FOLDER)

    bootstrap, bootstrap_options = _validate_bootstrap_options(payload)

    return {
        "species_tree": species_tree,
        "gene_trees": gene_trees,
        "outgroup": outgroups,
        "triplet_filter": _validate_optional_path(payload, "triplet_filter"),
        "output": output,
        "overwrite": _validate_overwrite_flag(payload),
        "processes": _validate_non_negative_int(
            payload, "processes", DEFAULT_PROCESSES
        ),
        "parallelization_mode": _validate_choice(
            payload,
            "parallelization_mode",
            DEFAULT_PARALLELIZATION_MODE,
            PARALLELIZATION_MODE_CHOICES,
        ),
        "consolidation": _validate_optional_bool(
            payload, "consolidation", DEFAULT_CONSOLIDATION
        ),
        "generate_summary_stats": _validate_optional_bool(
            payload, "generate_summary_stats", DEFAULT_GENERATE_SUMMARY_STATS
        ),
        "min_support_value": _validate_optional_float(
            payload, "min_support_value", DEFAULT_MIN_SUPPORT_VALUE
        ),
        "discordant_test": _validate_choice(
            payload, "discordant_test", DEFAULT_DISCORDANT_TEST, DISCORDANT_TEST_CHOICES
        ),
        "summary_statistic": _validate_choice(
            payload,
            "summary_statistic",
            DEFAULT_SUMMARY_STATISTIC,
            SUMMARY_STATISTIC_CHOICES,
        ),
        "tree_height_calculation_strategy": _validate_choice(
            payload,
            "tree_height_calculation_strategy",
            DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
            TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
        ),
        "p_value_correction": _validate_choice(
            payload,
            "p_value_correction",
            DEFAULT_P_VALUE_CORRECTION,
            P_VALUE_CORRECTION_CHOICES,
        ),
        "alpha_dct": _validate_optional_float(payload, "alpha_dct", DEFAULT_ALPHA_DCT),
        "alpha_ks": _validate_optional_float(payload, "alpha_ks", DEFAULT_ALPHA_KS),
        "bootstrap": bootstrap,
        "bootstrap_iterations": bootstrap_options["iterations"],
        "bootstrap_seed": bootstrap_options["seed"],
        "bootstrap_debug_mode": bootstrap_options["debug_mode"],
        "bootstrap_summary_only": bootstrap_options["summary_only"],
    }


def load_pipeline_config(config_file: str) -> dict:
    """Load and normalize a pipeline config from a JSON/YAML file.

    Args:
        config_file: Path to the JSON or YAML config file.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.pipeline.runner.run_pipeline`.

    Raises:
        ConfigError: If the file is malformed or a value is invalid.
    """
    payload = _load_raw_config(config_file)
    return normalize_pipeline_payload(payload)


def resolve_config(args: argparse.Namespace) -> dict:
    """Resolve a parsed CLI namespace into the pipeline runtime config.

    In config-file mode (``-c/--config-file``) the file is loaded and the other
    CLI flags are ignored with a warning; otherwise the exposed CLI flags are
    normalized directly. Config-file-only keys take their defaults in CLI mode.

    Args:
        args: Parsed ``argparse.Namespace`` from :func:`build_argument_parser`.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.pipeline.runner.run_pipeline`.

    Raises:
        ConfigError: If a required argument is missing or a value is invalid.
    """
    resolved = resolve_cli_or_config_args(
        args,
        load_config=load_pipeline_config,
        normalize_payload=normalize_pipeline_payload,
        payload_arg_names=_PIPELINE_PAYLOAD_ARG_NAMES,
    )
    return vars(resolved)
