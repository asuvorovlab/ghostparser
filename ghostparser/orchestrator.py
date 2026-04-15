"""End-to-end orchestration for GhostParser.

This module orchestrates both existing stages:
1) tree preprocessing / triplet extraction (from ``tree_parser``)
2) per-triplet introgression inference (from ``triplet_processor``)

It accepts the same CLI options as ``tree_parser`` and runs per-triplet
processing in parallel. Per-triplet statistics dictionaries are accumulated and
written to a final TSV report.
"""

from __future__ import annotations

import argparse
from multiprocessing import cpu_count
from pathlib import Path
import time

import dendropy

from .cli_config import resolve_cli_or_config_args

from .tree_parser import (
    _build_species_triplet_metadata,
    _parse_outgroup_arg,
    _root_tree_on_outgroup,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    MetricsLogger,
    read_triplet_filter_file,
    read_tree_file,
    write_clean_trees,
    write_triplet_gene_trees_multiprocess,
)
from .triplet_processor import (
    analyze_triplet_gene_tree_file,
    write_summary_statistics_tsv,
    write_pipeline_results,
)
from .config import (
    ConfigError,
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_MIN_SUPPORT_VALUE,
    DEFAULT_OUTPUT_FOLDER,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_STATS_BACKEND,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    P_VALUE_CORRECTION_CHOICES,
    SUMMARY_STATISTIC_CHOICES,
    STATS_BACKEND_CHOICES,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
    load_orchestrator_config,
    normalize_orchestrator_payload,
)


def _default_output_path():
    """Default output directory path for CLI mode."""
    return str(Path.cwd() / DEFAULT_OUTPUT_FOLDER)


def _resolve_processes(processes):
    """Resolve process count where 0 means all available CPU cores."""
    if processes == 0:
        return cpu_count()
    return processes


def _resolve_parallel_mode(processes):
    """Resolve worker count and whether multiprocessing should be enabled."""
    resolved_processes = _resolve_processes(processes)
    use_multiprocessing = resolved_processes is not None and resolved_processes > 1
    return resolved_processes, use_multiprocessing


def _now_times():
    """Return current wall-clock and CPU-process times."""
    return time.time(), time.process_time()


def _elapsed_times(start_wall, start_cpu):
    """Return elapsed wall-clock and CPU-process times."""
    return time.time() - start_wall, time.process_time() - start_cpu


def _log_stage_timing(metrics, wall_seconds, cpu_seconds):
    """Log wall-clock and CPU time for a stage."""
    metrics.log(f"  Time taken (wall): {wall_seconds:.2f}s")
    metrics.log(f"  Time taken (CPU): {cpu_seconds:.2f}s")


ORCHESTRATOR_PAYLOAD_ARG_NAMES = [
    "species_tree_path",
    "gene_trees_path",
    "outgroups",
    "triplet_filter",
    "output_folder",
    "processes",
    "generate_summary_stats",
    "min_support_value",
    "discordant_test",
    "summary_statistic",
    "stats_backend",
    "tree_height_calculation_strategy",
    "p_value_correction",
    "alpha_dct",
    "alpha_ks",
    "bootstrap",
    "bootstrap_iterations",
    "bootstrap_seed",
    "bootstrap_debug_mode",
    "bootstrap_summary_only",
]


def _build_argument_parser():
    """Build the orchestrator CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="GhostParser orchestrator: run tree_parser and triplet_processor end-to-end."
    )

    parser.add_argument("-c", "--config-file", type=str, default=None, help="Path to a JSON or YAML config file")
    parser.add_argument("-st", "--species-tree-path", default=None, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene-trees-path", default=None, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroups", default=None, help="Outgroup species identifier(s), comma-separated")
    parser.add_argument(
        "--triplet-filter",
        type=str,
        default=None,
        help="Path to triplet filter file (comma-separated taxa per line)",
    )
    parser.add_argument(
        "--output-folder",
        type=str,
        default=None,
        help=f"Output folder (default: ./{DEFAULT_OUTPUT_FOLDER})",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes for triplet extraction/processing (0 = all cores)",
    )
    parser.add_argument(
        "--generate-summary-stats",
        dest="generate_summary_stats",
        action="store_true",
        default=None,
        help="Generate summary_statistics.tsv output (default: False)",
    )
    parser.add_argument(
        "--min-support-value",
        type=float,
        default=None,
        help=f"Support threshold for tree cleaning (default: {DEFAULT_MIN_SUPPORT_VALUE})",
    )
    parser.add_argument(
        "--discordant-test",
        choices=("chi-square", "z-test"),
        default=None,
        help=f"Discordant count test (default: {DEFAULT_DISCORDANT_TEST})",
    )
    parser.add_argument(
        "--summary-statistic",
        choices=SUMMARY_STATISTIC_CHOICES,
        default=None,
        help=f"Summary statistic after KS test (default: {DEFAULT_SUMMARY_STATISTIC})",
    )
    parser.add_argument(
        "--stats-backend",
        choices=STATS_BACKEND_CHOICES,
        default=None,
        help=f"Statistical backend for DCT/KS (default: {DEFAULT_STATS_BACKEND})",
    )
    parser.add_argument(
        "--tree-height-calculation-strategy",
        choices=TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
        default=None,
        help=(
            "Tree-height strategy: AVG uses mean root-to-tip distance, "
            "A/B/C use the selected taxon's root-to-tip distance, "
            "SIS uses sister-taxon distance, and INT uses internal branch length "
            f"(default: {DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY})"
        ),
    )
    parser.add_argument(
        "--p-value-correction",
        choices=P_VALUE_CORRECTION_CHOICES,
        default=None,
        help=f"Multiple-testing correction for triplet p-values (default: {DEFAULT_P_VALUE_CORRECTION})",
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
        "--bootstrap",
        action="store_true",
        help="Enable bootstrap sampling-with-replacement during triplet inference",
    )
    parser.add_argument(
        "--bootstrap-iterations",
        type=int,
        default=None,
        help="Number of bootstrap iterations per triplet (default: 100)",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=None,
        help="Optional bootstrap random seed for reproducibility",
    )
    parser.add_argument(
        "--bootstrap-debug-mode",
        dest="bootstrap_debug_mode",
        action="store_true",
        default=None,
        help="Enable bootstrap debug outputs (detailed bootstrap metric columns)",
    )
    parser.add_argument(
        "--bootstrap-summary-only",
        dest="bootstrap_summary_only",
        action="store_true",
        default=None,
        help="When bootstrap debug mode is enabled, emit compact summaries instead of full per-iteration lists",
    )
    return parser


def _resolve_runtime_args(args):
    """Resolve runtime arguments from either config-file mode or pure CLI mode."""
    return resolve_cli_or_config_args(
        args,
        load_config=load_orchestrator_config,
        normalize_payload=normalize_orchestrator_payload,
        payload_arg_names=ORCHESTRATOR_PAYLOAD_ARG_NAMES,
    )


def main():
    """CLI entry point for orchestrated end-to-end run."""
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except (ValueError, ConfigError) as exc:
        print(f"Error: {exc}")
        return

    species_tree_path = Path(args.species_tree)
    gene_trees_path = Path(args.gene_trees)

    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {args.species_tree}")
        return

    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {args.gene_trees}")
        return

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    species_tree_clean = str(output_dir / f"processed_{species_tree_path.name}")
    gene_trees_clean = str(output_dir / f"processed_{gene_trees_path.name}")
    metrics_filepath = str(output_dir / "metrics.txt")

    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        outgroup_taxa = _parse_outgroup_arg(args.outgroup)
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log(f"Discordant count test: {args.discordant_test}")
        metrics.log(f"Summary statistic after KS: {args.summary_statistic}")
        metrics.log(f"Statistical backend: {args.stats_backend}")
        metrics.log(f"Tree height strategy: {args.tree_height_calculation_strategy}")
        metrics.log(f"P-value correction: {args.p_value_correction}")
        metrics.log(f"DCT alpha: {args.alpha_dct}")
        metrics.log(f"KS alpha: {args.alpha_ks}")
        metrics.log(f"Bootstrap enabled: {args.bootstrap}")
        metrics.log(f"Bootstrap iterations: {args.bootstrap_options['iterations']}")
        metrics.log(f"Bootstrap seed: {args.bootstrap_options['seed']}")
        metrics.log(f"Bootstrap debug mode: {args.bootstrap_options['debug_mode']}")
        metrics.log(f"Bootstrap summary-only: {args.bootstrap_options['summary_only']}")
        metrics.log(f"Generate summary statistics TSV: {args.generate_summary_stats}")
        support_threshold = (
            args.min_support_value
            if args.min_support_value is not None
            else DEFAULT_MIN_SUPPORT_VALUE
        )
        metrics.log(f"Support threshold: {support_threshold}")
        metrics.log("")

        try:
            species_start_wall, species_start_cpu = _now_times()
            species_trees, dropped_species = clean_and_save_trees(
                str(species_tree_path),
                species_tree_clean,
                min_avg_support=support_threshold,
            )
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < {support_threshold}"
                )

            species_trees = read_tree_file(species_tree_clean)
            species_wall_time, species_cpu_time = _elapsed_times(species_start_wall, species_start_cpu)
            _log_stage_timing(metrics, species_wall_time, species_cpu_time)
        except Exception as exc:
            metrics.log(f"✗ Error processing species tree: {exc}")
            return

        try:
            triplets = []
            species_triplet_trees = {}
            if species_trees:
                taxa = get_taxa_from_tree(species_trees[0])
                metrics.log(f"\n✓ Found {len(taxa)} taxa in species tree")

                pruned_tree, _excluded_taxa, missing_taxa, ingroup_taxa = _root_tree_on_outgroup(
                    species_trees[0], outgroup_taxa
                )

                if missing_taxa:
                    metrics.log(
                        f"⚠ Warning: Outgroup taxa not found in species tree: {', '.join(sorted(missing_taxa))}"
                    )

                if pruned_tree is None or not ingroup_taxa:
                    metrics.log("⚠ Warning: Unable to root and prune species tree on outgroups")
                    return

                species_trees = [pruned_tree]
                write_clean_trees(species_trees, species_tree_clean)

                if args.triplet_filter:
                    filter_path = Path(args.triplet_filter)
                    if not filter_path.exists():
                        metrics.log(f"✗ Error: Triplet filter file not found: {args.triplet_filter}")
                        return

                    raw_triplets, invalid_lines = read_triplet_filter_file(str(filter_path))
                    for line_number, raw in invalid_lines:
                        metrics.log(f"⚠ Warning: Skipping invalid triplet line {line_number} in {filter_path}: {raw}")

                    triplets, skipped_triplets = filter_triplets_by_taxa(raw_triplets, set(ingroup_taxa))
                    for triplet, missing in skipped_triplets:
                        metrics.log(
                            "⚠ Warning: Skipping triplet with missing taxa: "
                            f"{','.join(triplet)} (missing: {', '.join(missing)})"
                        )
                else:
                    triplets = generate_triplets(sorted(ingroup_taxa), [])
                    metrics.log(f"✓ Generated {len(triplets)} unique triplets")

                species_tree_newick = format_newick_with_precision(pruned_tree)
                species_dendro_tree = dendropy.Tree.get(
                    data=species_tree_newick,
                    schema="newick",
                    preserve_underscores=True,
                )

                triplets, species_triplet_trees, skipped_species_triplets = _build_species_triplet_metadata(
                    species_dendro_tree,
                    triplets,
                )

                if skipped_species_triplets:
                    metrics.log(
                        "⚠ Warning: Skipping triplets that could not be mapped on species tree: "
                        + "; ".join(",".join(triplet) for triplet in skipped_species_triplets)
                    )

                metrics.log(f"✓ Normalized {len(triplets)} triplets to A,B,C (A and B are sisters)")
        except Exception as exc:
            metrics.log(f"✗ Error generating triplets: {exc}")
            return

        try:
            genes_start_wall, genes_start_cpu = _now_times()
            gene_trees, dropped_genes, rooted_count, missing_root_indices = clean_and_save_gene_trees(
                str(gene_trees_path),
                gene_trees_clean,
                outgroup_taxa,
                min_avg_support=support_threshold,
            )
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            metrics.log(f"  Rooted {rooted_count} tree(s) on outgroup taxa")
            if missing_root_indices:
                metrics.log(f"  Discarded {len(missing_root_indices)} gene tree(s) without outgroup taxa")
            if dropped_genes:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < {support_threshold}"
                )
            genes_wall_time, genes_cpu_time = _elapsed_times(genes_start_wall, genes_start_cpu)
            _log_stage_timing(metrics, genes_wall_time, genes_cpu_time)
        except Exception as exc:
            metrics.log(f"✗ Error processing gene trees: {exc}")
            return

        processes, use_multiprocessing = _resolve_parallel_mode(args.processes)
        metrics.log(f"\nNumber of parallel cores being utilized: {processes}")
        metrics.log("")

        try:
            triplet_output_path = str(output_dir / "unique_triplets_gene_trees.txt")
            extraction_wall_time = 0.0
            extraction_cpu_time = 0.0
            inference_wall_time = 0.0
            inference_cpu_time = 0.0

            # Stage 1: Triplet extraction to mapping file
            metrics.log("✓ Starting triplet extraction stage...")
            extraction_start_wall, extraction_start_cpu = _now_times()

            total_subtrees, triplets_with_trees, _ = write_triplet_gene_trees_multiprocess(
                triplets,
                gene_trees_clean,
                triplet_output_path,
                species_triplet_trees=species_triplet_trees,
                use_multiprocessing=use_multiprocessing,
                processes=processes,
            )

            extraction_wall_time, extraction_cpu_time = _elapsed_times(
                extraction_start_wall,
                extraction_start_cpu,
            )
            metrics.log("✓ Triplet extraction complete")
            metrics.log(f"  Output: {triplet_output_path}")
            metrics.log(f"  Triplets with gene trees: {triplets_with_trees}")
            metrics.log(f"  Total subtrees extracted: {total_subtrees}")
            _log_stage_timing(metrics, extraction_wall_time, extraction_cpu_time)
            metrics.log("")

            # Stage 2: Introgression inference from mapping file
            metrics.log("✓ Starting introgression inference stage...")
            inference_start_wall, inference_start_cpu = _now_times()

            results = analyze_triplet_gene_tree_file(
                triplet_output_path,
                alpha_dct=args.alpha_dct,
                alpha_ks=args.alpha_ks,
                discordant_test=args.discordant_test,
                summary_statistic=args.summary_statistic,
                stats_backend=args.stats_backend,
                tree_height_calculation_strategy=args.tree_height_calculation_strategy,
                p_value_correction=args.p_value_correction,
                bootstrap=args.bootstrap,
                bootstrap_options=args.bootstrap_options,
                use_multiprocessing=use_multiprocessing,
                processes=processes,
            )
            inference_wall_time, inference_cpu_time = _elapsed_times(
                inference_start_wall,
                inference_start_cpu,
            )

            final_tsv = str(output_dir / "orchestrator_triplet_results.tsv")
            write_pipeline_results(
                results,
                final_tsv,
                dct_method=args.discordant_test,
                summary_statistic=args.summary_statistic,
                p_value_correction=args.p_value_correction,
                bootstrap=args.bootstrap,
                bootstrap_debug_mode=args.bootstrap_options["debug_mode"],
            )
            summary_tsv = str(output_dir / "summary_statistics.tsv")
            if args.generate_summary_stats:
                write_summary_statistics_tsv(
                    results,
                    summary_tsv,
                    bootstrap=args.bootstrap,
                )

            metrics.log("✓ Introgression inference complete")
            metrics.log(f"  Output: {final_tsv}")
            if args.generate_summary_stats:
                metrics.log(f"  Summary statistics output: {summary_tsv}")
            else:
                metrics.log("  Summary statistics output: skipped")
            metrics.log(f"  Triplets analyzed: {len(results)}")
            if use_multiprocessing:
                metrics.log("  Parallelization: enabled (multiprocessing)")
            else:
                metrics.log("  Parallelization: disabled (single-worker mode)")
            _log_stage_timing(metrics, inference_wall_time, inference_cpu_time)
            metrics.log("")

            # Summary
            total_wall_time = species_wall_time + genes_wall_time + extraction_wall_time + inference_wall_time
            total_cpu_time = species_cpu_time + genes_cpu_time + extraction_cpu_time + inference_cpu_time
            metrics.log("✓ End-to-end orchestration complete")
            metrics.log(f"  Total time (wall): {total_wall_time:.2f}s")
            metrics.log(f"  Total time (CPU): {total_cpu_time:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error in triplet inference or introgression inference stage: {exc}")
            return

        metrics.log("\nProcessing complete!")
        metrics.log(f"\nMetrics saved to: {metrics_filepath}")


if __name__ == "__main__":
    main()
