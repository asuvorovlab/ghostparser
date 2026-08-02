"""CLI parsing and config resolution for the streaming pipeline.

Turns process argv (via an argparse ``Namespace``) or a JSON/YAML config file
into the resolved config dict consumed by
:func:`ghostparser.pipeline.runner.run_pipeline`.

Config-file mode mirrors the orchestrator: ``-c/--config-file`` is CLI-only, and
when a config file is given the other CLI flags are ignored with a warning
(config wins). A subset of runtime knobs is exposed both on the CLI and in the
config file; the rest are config-file-only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ghostparser.cli_config import resolve_cli_or_config_args
from ghostparser.config import (
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_CONSOLIDATION,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_GENERATE_SUMMARY_STATS,
    DEFAULT_MIN_SUPPORT_VALUE,
    DEFAULT_OUTPUT_FOLDER,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_PROCESSES,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    DISCORDANT_TEST_CHOICES,
    P_VALUE_CORRECTION_CHOICES,
    SUMMARY_STATISTIC_CHOICES,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
    _load_raw_config,
    _parse_outgroups,
    _validate_bootstrap_options,
    _validate_choice,
    _validate_non_negative_int,
    _validate_optional_bool,
    _validate_optional_float,
    _validate_optional_path,
    _validate_overwrite_flag,
    _validate_required_path,
)

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
    "consolidation",
    "alpha_dct",
    "alpha_ks",
    "p_value_correction",
    "summary_statistic",
    "bootstrap",
]


def _resolve_path(path_str: str) -> str:
    """Resolve a path string to an absolute path.

    Handles ``~`` expansion and relative paths (resolved from the current
    working directory).

    Args:
        path_str: Path string to resolve.

    Returns:
        The absolute path as a string.
    """
    return str(Path(path_str).expanduser().resolve())


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
    or a CLI-derived dict), applying shared defaults for anything omitted.
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
