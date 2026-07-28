"""CLI parsing and default resolution for the streaming pipeline.

Turns process argv (via an argparse ``Namespace``) into the resolved config dict
consumed by :func:`ghostparser.pipeline.runner.run_pipeline`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ghostparser.config import (
    ConfigError,
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_BOOTSTRAP,
    DEFAULT_BOOTSTRAP_ITERATIONS,
    DEFAULT_CONSOLIDATION,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_MIN_SUPPORT_VALUE,
    DEFAULT_OUTPUT_FOLDER,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_PROCESSES,
    DEFAULT_STATS_BACKEND,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
)

# Pipeline-specific parallelization knobs (no orchestrator equivalent to reuse).
PARALLELIZATION_MODE_CHOICES = ("auto", "taxon", "gene")
DEFAULT_PARALLELIZATION_MODE = "auto"
# `auto` picks `gene` for small-taxa/large-gene-tree runs, else `taxon`.
AUTO_TAXA_SMALL_THRESHOLD = 15
AUTO_GENE_TREES_THRESHOLD = 3500


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
        "-st",
        "--species-tree-path",
        required=True,
        help="Path to the species tree file in Newick format",
    )
    parser.add_argument(
        "-gt",
        "--gene-trees-path",
        required=True,
        help="Path to the gene trees file in Newick format",
    )
    parser.add_argument(
        "-og",
        "--outgroups",
        required=True,
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


def resolve_config(args: argparse.Namespace) -> dict:
    """Resolve a parsed CLI namespace into the pipeline runtime config.

    Fills in every runtime key: the resolved required paths and outgroups, the
    exposed optional flags, and the pinned shared defaults that are not exposed
    as flags in v1.

    Args:
        args: Parsed ``argparse.Namespace`` from :func:`build_argument_parser`.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.pipeline.runner.run_pipeline`.

    Raises:
        ConfigError: If a required argument is missing or ``--processes`` is
            negative.
    """
    species_tree = getattr(args, "species_tree_path", None)
    gene_trees = getattr(args, "gene_trees_path", None)
    outgroups_raw = getattr(args, "outgroups", None)

    if not species_tree:
        raise ConfigError("Missing required argument: --species-tree-path")
    if not gene_trees:
        raise ConfigError("Missing required argument: --gene-trees-path")
    if not outgroups_raw:
        raise ConfigError("Missing required argument: --outgroups")

    outgroups = [part.strip() for part in str(outgroups_raw).split(",") if part.strip()]
    if not outgroups:
        raise ConfigError("Missing required argument: --outgroups")

    output_folder = getattr(args, "output_folder", None)
    output = _resolve_path(output_folder or DEFAULT_OUTPUT_FOLDER)

    triplet_filter = getattr(args, "triplet_filter", None)
    triplet_filter = _resolve_path(triplet_filter) if triplet_filter else None

    processes = getattr(args, "processes", None)
    if processes is None:
        processes = DEFAULT_PROCESSES
    if not isinstance(processes, int) or processes < 0:
        raise ConfigError("--processes must be an integer >= 0")

    parallelization_mode = getattr(args, "parallelization_mode", None)
    if parallelization_mode is None:
        parallelization_mode = DEFAULT_PARALLELIZATION_MODE

    consolidation = getattr(args, "consolidation", None)
    if consolidation is None:
        consolidation = DEFAULT_CONSOLIDATION

    bootstrap = getattr(args, "bootstrap", None)
    if bootstrap is None:
        bootstrap = DEFAULT_BOOTSTRAP

    return {
        "species_tree": _resolve_path(species_tree),
        "gene_trees": _resolve_path(gene_trees),
        "outgroup": outgroups,
        "triplet_filter": triplet_filter,
        "output": output,
        "processes": processes,
        "parallelization_mode": parallelization_mode,
        "consolidation": consolidation,
        "bootstrap": bootstrap,
        # Pinned defaults: not exposed as flags in v1.
        "overwrite": True,
        "min_support_value": DEFAULT_MIN_SUPPORT_VALUE,
        "discordant_test": DEFAULT_DISCORDANT_TEST,
        "summary_statistic": DEFAULT_SUMMARY_STATISTIC,
        "stats_backend": DEFAULT_STATS_BACKEND,
        "tree_height_calculation_strategy": DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
        "p_value_correction": DEFAULT_P_VALUE_CORRECTION,
        "alpha_dct": DEFAULT_ALPHA_DCT,
        "alpha_ks": DEFAULT_ALPHA_KS,
        "bootstrap_iterations": DEFAULT_BOOTSTRAP_ITERATIONS,
        "bootstrap_seed": None,
        "generate_summary_stats": False,
    }
