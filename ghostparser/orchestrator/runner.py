"""Orchestrator coordinator: run_orchestrator drives cleaning, triplet setup, the fused streaming engine, correction, TSV writing, and consolidation end-to-end with per-stage timing."""

from __future__ import annotations

import time
from pathlib import Path

import dendropy

from ghostparser.introgression_mapper import generate_introgression_maps

from .config import prepare_output_directory, resolve_config
from .inference import write_pipeline_results, write_summary_statistics_tsv
from .preflight import run_preflight_data_check
from .stream import resolve_parallelization_mode, stream_triplet_results
from .trees import (
    MetricsLogger,
    _build_species_triplet_metadata,
    _parse_outgroup_arg,
    _read_gene_trees_file,
    _root_tree_on_outgroup,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    read_tree_file,
    read_triplet_filter_file,
    write_clean_trees,
)


def _now_times():
    """Capture the current wall and CPU time.

    Returns:
        A tuple ``(wall, cpu)`` of floats.
    """
    return time.time(), time.process_time()


def _elapsed_times(start_wall, start_cpu):
    """Compute elapsed wall and CPU time since a start.

    Args:
        start_wall: Wall-clock start time.
        start_cpu: CPU start time.

    Returns:
        A tuple ``(wall_elapsed, cpu_elapsed)`` of floats.
    """
    return time.time() - start_wall, time.process_time() - start_cpu


def _log_stage_timing(metrics, wall_time, cpu_time):
    """Log wall and CPU timing for one stage.

    Args:
        metrics: The ``MetricsLogger`` to write to.
        wall_time: Elapsed wall-clock seconds.
        cpu_time: Elapsed CPU seconds.
    """
    metrics.log(f"  Time taken (wall): {wall_time:.2f}s")
    metrics.log(f"  Time taken (CPU): {cpu_time:.2f}s")


def _format_metric_float(value):
    """Format an optional float for the metrics log.

    Args:
        value: A float, or ``None``.

    Returns:
        The value at 6 significant digits, or ``n/a`` for ``None``.
    """
    return "n/a" if value is None else f"{value:.6g}"


def _log_permutation_diagnostics(metrics, results):
    """Report permutation-test convergence and decision-rule agreement.

    Triplets whose adaptive run exhausted the resample budget still produce a
    decision, but that decision sits closer to the Monte Carlo error than the
    confidence-interval criterion is willing to certify, so they are named
    individually. The same is done for any triplet where the two one-tailed rule
    and the two-tailed-gate-then-sign rule disagree, which should not happen and
    signals that the inputs or the accumulators need looking at.

    Args:
        metrics: The ``MetricsLogger`` to write to.
        results: List of ``TripletPipelineResult`` objects.
    """
    # A guarded test reports zero resamples and no p-values, so "ran" means the
    # test actually resampled rather than merely having been attempted.
    ran = [result for result in results if result.perm_n_resamples]
    guarded = [
        result
        for result in results
        if not result.perm_n_resamples and result.perm_note
    ]
    no_comparison = [
        result
        for result in results
        if result.perm_n_resamples is None and result.perm_note is None
    ]
    unconverged = [result for result in ran if result.perm_converged is False]
    inconsistent = [result for result in ran if result.perm_consistent is False]
    skewed = [result for result in ran if result.perm_null_skewed]

    metrics.log(f"  Permutation tests run: {len(ran)}")
    metrics.log(
        f"  Permutation resamples drawn: {sum(result.perm_n_resamples for result in ran)}"
    )

    if no_comparison:
        metrics.log(
            f"  ⚠ No concordant/discordant1 heights to compare for "
            f"{len(no_comparison)} triplet(s):"
        )
        for result in no_comparison:
            metrics.log(f"      {','.join(result.triplet)}")
    if guarded:
        metrics.log(
            f"  ⚠ Permutation test guarded (not resampled) for {len(guarded)} triplet(s):"
        )
        for result in guarded:
            metrics.log(f"      {','.join(result.triplet)}: {result.perm_note}")
    if unconverged:
        metrics.log(
            f"  ⚠ Permutation test hit max_resamples without the confidence "
            f"interval excluding alpha for {len(unconverged)} triplet(s):"
        )
        for result in unconverged:
            metrics.log(
                f"      {','.join(result.triplet)}: "
                f"p_greater={_format_metric_float(result.perm_p_greater)}, "
                f"p_less={_format_metric_float(result.perm_p_less)}, "
                f"resamples={result.perm_n_resamples}"
            )
    if skewed:
        # Reported as a count rather than per triplet: tree heights are bounded
        # below by zero and routinely right-skewed, so most triplets in a real
        # dataset land here. It is context for reading the two-tailed column,
        # not a warning to act on -- the per-triplet flag is in the results TSV.
        metrics.log(
            f"  Asymmetric permutation null for {len(skewed)} of {len(ran)} "
            f"triplet(s); directional decisions are unaffected "
            f"(see the perm_null_skewed column)."
        )
    if inconsistent:
        metrics.log(
            f"  ⚠ Permutation decision rules disagree for {len(inconsistent)} triplet(s):"
        )
        for result in inconsistent:
            metrics.log(
                f"      {','.join(result.triplet)}: "
                f"one-tailed={result.perm_decision}, "
                f"p_two_sided={_format_metric_float(result.perm_p_two_sided)}"
            )


def _run_preflight_only(config, output_dir):
    """Run the structural preflight check and stop before any analysis.

    Args:
        config: A resolved config dict.
        output_dir: Prepared output directory the report is written into.

    Returns:
        The :class:`~ghostparser.orchestrator.preflight.PreflightResult`, or
        ``None`` when the check could not run at all (unrootable species tree,
        missing outgroups) — the reason is printed in that case.
    """
    outgroup_taxa = _parse_outgroup_arg(config["outgroup"])
    print("Preflight data check enabled — no analysis will be run.")
    try:
        result = run_preflight_data_check(
            species_tree_path=config["species_tree"],
            gene_trees_path=config["gene_trees"],
            outgroups=outgroup_taxa,
            output_dir=output_dir,
            triplet_filter=config["triplet_filter"],
        )
    except ValueError as exc:
        print(f"✗ Error: Preflight data check could not run: {exc}")
        return None

    print(result.report_text, end="")
    print(f"Saved preflight report: {result.report_path}")
    return result


def run_orchestrator(config):
    """Run the orchestrator end-to-end for one config.

    Cleans and roots the trees, builds and normalizes triplets, runs the fused
    streaming engine with run-wide p-value correction, writes the results TSV,
    and (unless disabled) runs consolidation into a dedicated subfolder, logging
    per-stage timing throughout.

    Args:
        config: A resolved config dict, or an ``argparse.Namespace`` which is
            resolved via :func:`ghostparser.orchestrator.config.resolve_config`.

    Returns:
        The list of ``TripletPipelineResult`` objects, or ``None`` if the run
        exits early on an input or preprocessing error.
    """
    if not isinstance(config, dict):
        config = resolve_config(config)

    species_tree_path = Path(config["species_tree"])
    gene_trees_path = Path(config["gene_trees"])

    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {config['species_tree']}")
        return None
    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {config['gene_trees']}")
        return None

    output_dir = Path(
        prepare_output_directory(config["output"], overwrite=config["overwrite"])
    )

    if config.get("preflight_data_check"):
        return _run_preflight_only(config, output_dir)

    species_tree_clean = str(output_dir / f"processed_{species_tree_path.name}")
    gene_trees_clean = str(output_dir / f"processed_{gene_trees_path.name}")
    metrics_filepath = str(output_dir / "metrics.txt")

    support_threshold = config["min_support_value"]

    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {config['species_tree']}")
        metrics.log(f"Processing gene trees: {config['gene_trees']}")
        outgroup_taxa = _parse_outgroup_arg(config["outgroup"])
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log(f"Discordant count test: {config['discordant_test']}")
        metrics.log(f"Permutation test enabled: {config['permutation_test']}")
        if config["permutation_test"]:
            metrics.log(
                "Permutation resamples: "
                f"{config['permutation_min_resamples']}-"
                f"{config['permutation_max_resamples']}"
            )
            metrics.log(f"Permutation CI method: {config['permutation_ci_method']}")
        metrics.log(
            f"Tree height strategy: {config['tree_height_calculation_strategy']}"
        )
        metrics.log(f"P-value correction: {config['p_value_correction']}")
        metrics.log(f"DCT alpha: {config['alpha_dct']}")
        metrics.log(f"KS alpha: {config['alpha_ks']}")
        metrics.log(f"Permutation alpha: {config['alpha_perm']}")
        metrics.log(f"Bootstrap enabled: {config['bootstrap']}")
        metrics.log(f"Bootstrap iterations: {config['bootstrap_iterations']}")
        metrics.log(f"Bootstrap debug mode: {config['bootstrap_debug_mode']}")
        metrics.log(f"Generate summary statistics TSV: {config['generate_summary_stats']}")
        metrics.log(f"Parallelization mode: {config['parallelization_mode']}")
        metrics.log(f"Consolidation enabled: {config['consolidation']}")
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
            species_wall_time, species_cpu_time = _elapsed_times(
                species_start_wall, species_start_cpu
            )
            _log_stage_timing(metrics, species_wall_time, species_cpu_time)
        except Exception as exc:
            metrics.log(f"✗ Error processing species tree: {exc}")
            return None

        triplets = []
        species_triplet_trees = {}
        plot_taxa = None
        try:
            if species_trees:
                taxa = get_taxa_from_tree(species_trees[0])
                metrics.log(f"\n✓ Found {len(taxa)} taxa in species tree")

                pruned_tree, _excluded_taxa, missing_taxa, ingroup_taxa = (
                    _root_tree_on_outgroup(species_trees[0], outgroup_taxa)
                )

                if missing_taxa:
                    metrics.log(
                        f"⚠ Warning: Outgroup taxa not found in species tree: {', '.join(sorted(missing_taxa))}"
                    )

                if pruned_tree is None or not ingroup_taxa:
                    metrics.log(
                        "⚠ Warning: Unable to root and prune species tree on outgroups"
                    )
                    return None

                species_trees = [pruned_tree]
                write_clean_trees(species_trees, species_tree_clean)

                if config["triplet_filter"]:
                    filter_path = Path(config["triplet_filter"])
                    if not filter_path.exists():
                        metrics.log(
                            f"✗ Error: Triplet filter file not found: {config['triplet_filter']}"
                        )
                        return None

                    raw_triplets, invalid_lines = read_triplet_filter_file(
                        str(filter_path)
                    )
                    for line_number, raw in invalid_lines:
                        metrics.log(
                            f"⚠ Warning: Skipping invalid triplet line {line_number} in {filter_path}: {raw}"
                        )

                    triplets, skipped_triplets = filter_triplets_by_taxa(
                        raw_triplets, set(ingroup_taxa)
                    )
                    plot_taxa = sorted(
                        {taxon for triplet in triplets for taxon in triplet}
                    )
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

                triplets, species_triplet_trees, skipped_species_triplets = (
                    _build_species_triplet_metadata(
                        species_dendro_tree,
                        triplets,
                    )
                )

                if skipped_species_triplets:
                    metrics.log(
                        "⚠ Warning: Skipping triplets that could not be mapped on species tree: "
                        + "; ".join(
                            ",".join(triplet) for triplet in skipped_species_triplets
                        )
                    )

                metrics.log(
                    f"✓ Normalized {len(triplets)} triplets to A,B,C (A and B are sisters)"
                )
        except Exception as exc:
            metrics.log(f"✗ Error generating triplets: {exc}")
            return None

        try:
            genes_start_wall, genes_start_cpu = _now_times()
            gene_trees, dropped_genes, rooted_count, missing_root_indices = (
                clean_and_save_gene_trees(
                    str(gene_trees_path),
                    gene_trees_clean,
                    outgroup_taxa,
                    min_avg_support=support_threshold,
                )
            )
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            metrics.log(f"  Rooted {rooted_count} tree(s) on outgroup taxa")
            if missing_root_indices:
                metrics.log(
                    f"  Discarded {len(missing_root_indices)} gene tree(s) without outgroup taxa"
                )
            if dropped_genes:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < {support_threshold}"
                )
            genes_wall_time, genes_cpu_time = _elapsed_times(
                genes_start_wall, genes_start_cpu
            )
            _log_stage_timing(metrics, genes_wall_time, genes_cpu_time)
        except Exception as exc:
            metrics.log(f"✗ Error processing gene trees: {exc}")
            return None

        try:
            gene_trees_newick = _read_gene_trees_file(gene_trees_clean)
            n_taxa = len(ingroup_taxa)
            resolved_mode = resolve_parallelization_mode(
                config["parallelization_mode"], n_taxa, len(gene_trees_newick)
            )
            metrics.log(f"\n✓ Resolved parallelization mode: {resolved_mode}")

            inference_kwargs = {
                "alpha_dct": config["alpha_dct"],
                "alpha_ks": config["alpha_ks"],
                "discordant_test": config["discordant_test"],
                "permutation_test": config["permutation_test"],
                "permutation_kwargs": {
                    "alpha": config["alpha_perm"],
                    "min_resamples": config["permutation_min_resamples"],
                    "max_resamples": config["permutation_max_resamples"],
                    "ci_method": config["permutation_ci_method"],
                    "correction": config["p_value_correction"],
                },
                "tree_height_calculation_strategy": config[
                    "tree_height_calculation_strategy"
                ],
                "collect_summary_statistics": config["generate_summary_stats"],
                "bootstrap_options": {
                    "iterations": config["bootstrap_iterations"],
                    "debug_mode": config["bootstrap_debug_mode"],
                    "summary_only": config["bootstrap_summary_only"],
                },
                "triplet_seed": config["bootstrap_seed"],
            }

            metrics.log("✓ Starting fused extraction + inference stage...")
            stream_start_wall, stream_start_cpu = _now_times()
            results, stream_worker_cpu = stream_triplet_results(
                triplets,
                species_triplet_trees,
                gene_trees_newick,
                mode=config["parallelization_mode"],
                processes=config["processes"],
                inference_kwargs=inference_kwargs,
                p_value_correction=config["p_value_correction"],
                return_worker_cpu=True,
            )
            stream_wall_time, stream_cpu_time = _elapsed_times(
                stream_start_wall, stream_start_cpu
            )
            # process_time() misses CPU spent in pool workers; add it back so the
            # reported CPU total reflects all work, not just the parent process.
            stream_cpu_time += stream_worker_cpu

            final_tsv = str(output_dir / "orchestrator_triplet_results.tsv")
            write_pipeline_results(
                results,
                final_tsv,
                dct_method=config["discordant_test"],
                p_value_correction=config["p_value_correction"],
                bootstrap=config["bootstrap"],
                bootstrap_debug_mode=config["bootstrap_debug_mode"],
                permutation_test=config["permutation_test"],
            )
            metrics.log("✓ Fused extraction + inference complete")
            metrics.log(f"  Output: {final_tsv}")
            metrics.log(f"  Triplets analyzed: {len(results)}")
            if config["permutation_test"]:
                _log_permutation_diagnostics(metrics, results)

            if config["generate_summary_stats"]:
                summary_tsv = str(output_dir / "summary_statistics.tsv")
                write_summary_statistics_tsv(
                    results,
                    summary_tsv,
                    bootstrap=config["bootstrap"],
                )
                metrics.log(f"  Summary statistics: {summary_tsv}")
            _log_stage_timing(metrics, stream_wall_time, stream_cpu_time)
            metrics.log("")

            map_wall_time = 0.0
            map_cpu_time = 0.0
            if config["consolidation"]:
                metrics.log("✓ Starting introgression map generation stage...")
                map_start_wall, map_start_cpu = _now_times()
                # Consolidation writes into a dedicated subfolder so its own
                # output-directory reset (rmtree) never touches the run folder's
                # results TSV, processed trees, or the still-open metrics.txt.
                consolidation_dir = output_dir / "consolidation"
                map_artifacts = generate_introgression_maps(
                    results,
                    species_tree_path=species_tree_clean,
                    output_dir=str(consolidation_dir),
                    plot_taxa=plot_taxa,
                    outgroups=outgroup_taxa,
                    overwrite=config["overwrite"],
                )
                map_wall_time, map_cpu_time = _elapsed_times(
                    map_start_wall, map_start_cpu
                )
                metrics.log("✓ Introgression map generation complete")
                metrics.log(f"  Output folder: {consolidation_dir}")
                metrics.log(f"  Taxa represented: {map_artifacts.taxa_count}")
                metrics.log(f"  Combined plot: {map_artifacts.plot_path}")
                _log_stage_timing(metrics, map_wall_time, map_cpu_time)
                metrics.log("")
            else:
                metrics.log(
                    "✓ Introgression map generation stage skipped (consolidation disabled)"
                )
                metrics.log("")

            total_wall_time = (
                species_wall_time
                + genes_wall_time
                + stream_wall_time
                + map_wall_time
            )
            total_cpu_time = (
                species_cpu_time + genes_cpu_time + stream_cpu_time + map_cpu_time
            )
            metrics.log("✓ End-to-end pipeline complete")
            metrics.log(f"  Total time (wall): {total_wall_time:.2f}s")
            metrics.log(f"  Total time (CPU): {total_cpu_time:.2f}s")
        except Exception as exc:
            metrics.log(f"✗ Error in fused extraction/inference stage: {exc}")
            return None

        metrics.log("\nProcessing complete!")
        metrics.log(f"\nMetrics saved to: {metrics_filepath}")

    return results
