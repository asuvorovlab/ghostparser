"""Run the GhostParser orchestrator over each subdirectory of a parent directory.

This mirrors the folder-based orchestration pattern from the
Slurm shell script,
while executing all discovered subdirectories in a single Python process.

Each run is driven by a generated JSON config file rather than CLI flags,
because ``generate_summary_stats`` is a config-file-only orchestrator setting.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_SPECIES_TREE = (
    "/projects/asuvorov/ghostparser_shared/simulations_64way/species.tree"
)
DEFAULT_OUTGROUPS = "O"
DEFAULT_GENE_TREES_FILE = "trees_subst_1.tre"
DEFAULT_OUTPUT_SUBDIR = "bfn_correction_results"
DEFAULT_P_VALUE_CORRECTION = "bfn"


def build_parser() -> argparse.ArgumentParser:
    """Create CLI parser for parent-directory orchestrator runs.

    Returns:
        The configured ``argparse.ArgumentParser``.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run ghostparser.orchestrator once per subdirectory under "
            "--parent-dir, "
            "deriving gene tree and output paths from each folder."
        )
    )
    parser.add_argument(
        "--parent-dir",
        required=True,
        help="Directory that contains per-simulation subdirectories.",
    )
    parser.add_argument(
        "--species-tree",
        default=DEFAULT_SPECIES_TREE,
        help=f"Species tree path (default: {DEFAULT_SPECIES_TREE}).",
    )
    parser.add_argument(
        "--outgroups",
        default=DEFAULT_OUTGROUPS,
        help=f"Outgroups value passed to the orchestrator (default: {DEFAULT_OUTGROUPS}).",
    )
    parser.add_argument(
        "--gene-trees-file",
        default=DEFAULT_GENE_TREES_FILE,
        help=(
            "Gene tree filename expected in each subdirectory "
            f"(default: {DEFAULT_GENE_TREES_FILE})."
        ),
    )
    parser.add_argument(
        "--output-subdir",
        default=DEFAULT_OUTPUT_SUBDIR,
        help=(
            "Output subdirectory created inside each simulation subdirectory "
            f"(default: {DEFAULT_OUTPUT_SUBDIR})."
        ),
    )
    parser.add_argument(
        "--p-value-correction",
        default=DEFAULT_P_VALUE_CORRECTION,
        help=(
            "P-value correction method passed to the orchestrator "
            f"(default: {DEFAULT_P_VALUE_CORRECTION})."
        ),
    )
    parser.add_argument(
        "--no-generate-summary-stats",
        action="store_true",
        help="Disable generate_summary_stats when invoking the orchestrator.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print derived configs without executing them.",
    )
    return parser


def discover_subdirectories(parent_dir: Path) -> list[Path]:
    """Return sorted immediate subdirectories under parent_dir.

    Args:
        parent_dir: The parent directory to scan.

    Returns:
        The sorted list of immediate subdirectories.
    """
    return sorted(path for path in parent_dir.iterdir() if path.is_dir())


def build_orchestrator_config(
    species_tree: Path,
    gene_trees: Path,
    out_dir: Path,
    outgroups: str,
    p_value_correction: str,
    generate_summary_stats: bool,
) -> dict:
    """Build the orchestrator config payload for one run.

    Args:
        species_tree: Path to the species tree.
        gene_trees: Path to this folder's gene trees.
        out_dir: Output directory for this run.
        outgroups: Comma-separated outgroup taxa.
        p_value_correction: Correction method to apply run-wide.
        generate_summary_stats: Whether to also write summary_statistics.tsv.

    Returns:
        The config mapping to serialize as JSON.
    """
    return {
        "species_tree_path": str(species_tree),
        "gene_trees_path": str(gene_trees),
        "outgroup": outgroups,
        "output_folder": str(out_dir),
        "p_value_correction": p_value_correction,
        "generate_summary_stats": generate_summary_stats,
    }


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    parent_dir = Path(args.parent_dir).expanduser().resolve()
    species_tree = Path(args.species_tree).expanduser().resolve()
    generate_summary_stats = not args.no_generate_summary_stats

    if not parent_dir.exists() or not parent_dir.is_dir():
        raise FileNotFoundError(
            f"Parent directory does not exist or is not a directory: {parent_dir}"
        )
    if not species_tree.exists() or not species_tree.is_file():
        raise FileNotFoundError(
            f"Species tree does not exist or is not a file: {species_tree}"
        )

    folders = discover_subdirectories(parent_dir)
    if not folders:
        raise FileNotFoundError(
            f"No subdirectories found under parent directory: {parent_dir}"
        )

    print(f"Parent directory: {parent_dir}")
    print(f"Discovered folders: {len(folders)}")

    failures: list[tuple[str, int]] = []
    skipped_missing_gene_trees: list[str] = []

    for index, folder in enumerate(folders, start=1):
        gene_trees = folder / args.gene_trees_file
        out_dir = folder / args.output_subdir

        if not gene_trees.exists() or not gene_trees.is_file():
            skipped_missing_gene_trees.append(str(folder.name))
            print(
                f"[{index}/{len(folders)}] Skipping {folder.name}: missing {gene_trees}"
            )
            continue

        config = build_orchestrator_config(
            species_tree=species_tree,
            gene_trees=gene_trees,
            out_dir=out_dir,
            outgroups=args.outgroups,
            p_value_correction=args.p_value_correction,
            generate_summary_stats=generate_summary_stats,
        )

        print(f"[{index}/{len(folders)}] Processing folder: {folder.name}")
        print(f"Config: {json.dumps(config)}")

        if args.dry_run:
            continue

        # The orchestrator reads every setting from the config file, so write one
        # per folder and hand it to the CLI.
        with tempfile.NamedTemporaryFile(
            "w", suffix=".json", prefix=f"gp_{folder.name}_", delete=False
        ) as handle:
            json.dump(config, handle)
            config_path = handle.name

        command = [sys.executable, "-m", "ghostparser.orchestrator", "-c", config_path]
        print(f"Command: {shlex.join(command)}")
        try:
            completed = subprocess.run(command, check=False)
        finally:
            Path(config_path).unlink(missing_ok=True)

        if completed.returncode != 0:
            failures.append((folder.name, completed.returncode))

    print("Run summary:")
    print(f"  Total folders discovered: {len(folders)}")
    print(f"  Missing gene tree file: {len(skipped_missing_gene_trees)}")
    print(f"  Failed orchestrator runs: {len(failures)}")

    if skipped_missing_gene_trees:
        print("  Folders skipped (missing gene trees):")
        for folder_name in skipped_missing_gene_trees:
            print(f"    - {folder_name}")

    if failures:
        print("  Folders with non-zero exit codes:")
        for folder_name, return_code in failures:
            print(f"    - {folder_name}: exit code {return_code}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
