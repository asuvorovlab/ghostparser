"""Orchestrator coordinator: run_orchestrator drives cleaning, triplet setup, the fused streaming engine, correction, TSV writing, and consolidation end-to-end with per-stage timing."""

import os
import time
from collections import Counter
from dataclasses import replace
from pathlib import Path

import dendropy

from .config import prepare_output_directory, resolve_config
from .consolidation import generate_introgression_maps
from .correction import is_inline_correction
from .inference import (
    PERM_NOTE_NOT_CONSULTED,
    write_pipeline_results,
    write_summary_statistics_tsv,
)
from .preflight import run_preflight_data_check
from .stream import available_cpu_count, resolve_worker_count, stream_triplet_results
from .trees import (
    MetricsLogger,
    _build_species_triplet_metadata,
    _parse_outgroup_arg,
    _read_gene_trees_file,
    _root_tree_on_outgroup,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    filter_triplets_by_taxa,
    load_species_rename_map,
    format_newick_with_precision,
    generate_triplets,
    get_taxa_from_tree,
    read_tree_file,
    rename_newick_labels,
    rename_taxon_labels,
    read_species_filter_file,
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


def _describe_diagnostic(config):
    """Spell out what the point estimate skips under this run's correction.

    Args:
        config: The resolved orchestrator config.

    Returns:
        The ``diagnostic`` setting, followed when it is off by which tests the
        point estimate declines and a reminder that the results are unchanged.
    """
    if config["diagnostic"]:
        return "True (every test is measured for every triplet)"
    if is_inline_correction(config["p_value_correction"]):
        skipped = "tree-height and direction tests below a settled gate"
    else:
        skipped = (
            "direction test below a settled gate; the tree-height test is "
            f"measured for every triplet because {config['p_value_correction']} "
            "corrects it as a whole-run family"
        )
    return f"False (skips the {skipped}; results are identical to a diagnostic run)"


def _log_permutation_diagnostics(metrics, results):
    """Report permutation-test coverage and convergence.

    Every category is reported as a count rather than named triplet by triplet:
    on a real run these lists reach the hundreds, and the per-triplet detail is
    already in the results TSV under ``n_con``/``n_dis1``, ``perm_note``,
    ``perm_converged``, and the ``perm_p_*`` columns.

    Args:
        metrics: The ``MetricsLogger`` to write to.
        results: List of ``TripletPipelineResult`` objects.
    """
    # A guarded test reports zero resamples and no p-values, so "ran" means the
    # test actually resampled rather than merely having been attempted.
    ran = [result for result in results if result.perm_n_resamples]
    # A skipped test and a guarded one both leave the block empty, but they mean
    # opposite things: the first is a non-diagnostic run declining work the
    # cascade could not consult, the second is a test that could not be run.
    skipped = [
        result for result in results if result.perm_note == PERM_NOTE_NOT_CONSULTED
    ]
    ks_skipped = [result for result in results if result.ks_p_value is None]
    guarded = [
        result
        for result in results
        if not result.perm_n_resamples
        and result.perm_note
        and result.perm_note != PERM_NOTE_NOT_CONSULTED
    ]
    no_comparison = [
        result
        for result in results
        if result.perm_n_resamples is None and result.perm_note is None
    ]
    unconverged = [result for result in ran if result.perm_converged is False]
    skews = [
        abs(result.perm_null_skew)
        for result in ran
        if result.perm_null_skew is not None
    ]

    # How far the cascade let each triplet go, which is also what bounds how
    # much of the direction test's cost a non-diagnostic run can decline.
    dct_cleared = [result for result in results if result.dct_significant]
    tht_cleared = [result for result in dct_cleared if result.ks_significant]
    metrics.log(
        f"  Triplets clearing the discordant count gate: {len(dct_cleared)} "
        f"of {len(results)}"
    )
    metrics.log(
        f"  Triplets clearing the tree-height gate: {len(tht_cleared)} "
        f"of {len(results)}"
    )
    metrics.log(f"  Permutation tests run: {len(ran)}")
    metrics.log(
        f"  Permutation resamples drawn: {sum(result.perm_n_resamples for result in ran)}"
    )

    if ks_skipped:
        metrics.log(
            f"  Tree-height tests skipped as already settled: {len(ks_skipped)} "
            f"of {len(results)} triplet(s)"
        )
    if skipped:
        metrics.log(
            f"  Direction tests skipped as already settled: {len(skipped)} "
            f"of {len(results)} triplet(s)"
        )
    if no_comparison:
        metrics.log(
            f"  ⚠ No concordant/discordant1 heights to compare for "
            f"{len(no_comparison)} triplet(s); see the n_con and n_dis1 columns."
        )
    if guarded:
        # Which guard fired is worth keeping, since the four mean different
        # things; which triplet fired it is in the perm_note column.
        by_note = Counter(result.perm_note for result in guarded)
        reasons = ", ".join(
            f"{note}={count}" for note, count in sorted(by_note.items())
        )
        metrics.log(
            f"  ⚠ Permutation test guarded (not resampled) for {len(guarded)} "
            f"triplet(s): {reasons}"
        )
    if unconverged:
        metrics.log(
            f"  ⚠ Permutation test hit max_resamples without the confidence "
            f"interval excluding alpha for {len(unconverged)} triplet(s); "
            f"see the perm_converged and perm_p_* columns."
        )
    if skews:
        # A distribution summary rather than a count: tree heights are bounded
        # below by zero and routinely right-skewed, so some asymmetry is normal
        # and only its magnitude is informative. Per-triplet values are in the
        # perm_null_skew column.
        strong = sum(1 for value in skews if value >= 0.5)
        metrics.log(
            f"  Permutation null skewness: median |skew| "
            f"{sorted(skews)[len(skews) // 2]:.2f}, max {max(skews):.2f}; "
            f"{strong} of {len(skews)} triplet(s) at or above 0.5 "
            f"(see the perm_null_skew column)."
        )


def _rename_result_taxa(results, rename_map):
    """Rebuild every result under its display names.

    The run reads, roots, filters and measures in the trees' own labels; the
    display names enter here, once the decision pass is done and before any
    output is written. Only two fields name taxa: ``triplet`` and the
    ``species_tree`` Newick, whose labels are quoted as the names require.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        rename_map: Mapping of tree label to display name.

    Returns:
        The same list when the map is empty, otherwise a new list of renamed
        results in the same order.
    """
    if not rename_map:
        return results
    return [
        replace(
            result,
            triplet=tuple(rename_taxon_labels(result.triplet, rename_map)),
            species_tree=rename_newick_labels(result.species_tree, rename_map),
        )
        for result in results
    ]


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
            species_filter=config["species_filter"],
            max_triplets=config["preflight_triplet_cap"],
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
        # The whole run -- trees, outgroup, filter, log -- works in the trees'
        # own labels; the map is applied to the results just before they are
        # written, so no renamed label is ever read back out of a Newick. The
        # config layer has already read and validated the file.
        rename_map = {}
        if config["species_rename_map"]:
            rename_map = load_species_rename_map(config["species_rename_map"])
            metrics.log(
                f"Species rename map: {config['species_rename_map']} "
                f"({len(rename_map)} taxa)"
            )
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log(f"Discordant count test: {config['discordant_test']}")
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
        # One seed drives every random draw in the run. When none is
        # configured a fresh one is drawn and reported, so an exploratory run
        # stays reproducible from its own metrics file.
        run_seed = config["seed"]
        seed_origin = "configured"
        if run_seed is None:
            run_seed = int.from_bytes(os.urandom(8), "little")
            seed_origin = "generated"
        metrics.log(f"Seed: {run_seed} ({seed_origin})")
        metrics.log(f"Bootstrap enabled: {config['bootstrap']}")
        metrics.log(f"Bootstrap iterations: {config['bootstrap_iterations']}")
        metrics.log(f"Bootstrap debug mode: {config['bootstrap_debug_mode']}")
        metrics.log(f"Generate summary statistics TSV: {config['generate_summary_stats']}")
        metrics.log(f"Shape diagnostics: {config['shape_diagnostics']}")
        metrics.log(f"Diagnostic: {_describe_diagnostic(config)}")
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
                elif config["species_filter"]:
                    filter_path = Path(config["species_filter"])
                    if not filter_path.exists():
                        metrics.log(
                            f"✗ Error: Species filter file not found: {config['species_filter']}"
                        )
                        return None

                    requested_species = read_species_filter_file(str(filter_path))
                    ingroup_set = set(ingroup_taxa)
                    kept_species = [
                        taxon for taxon in requested_species if taxon in ingroup_set
                    ]
                    for taxon in requested_species:
                        if taxon in outgroup_taxa:
                            metrics.log(
                                f"⚠ Warning: Skipping species filter entry {taxon}: "
                                "it is an outgroup taxon"
                            )
                        elif taxon not in ingroup_set:
                            metrics.log(
                                f"⚠ Warning: Skipping species filter entry {taxon}: "
                                "not found in the species tree"
                            )
                    if len(kept_species) < 3:
                        metrics.log(
                            f"✗ Error: Species filter names {len(kept_species)} "
                            "ingroup taxa; at least 3 are needed to form a triplet"
                        )
                        return None

                    plot_taxa = sorted(kept_species)
                    triplets = generate_triplets(plot_taxa, [])
                    metrics.log(
                        f"✓ Generated {len(triplets)} unique triplets among "
                        f"{len(plot_taxa)} filtered species"
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

            inference_kwargs = {
                "alpha_dct": config["alpha_dct"],
                "alpha_ks": config["alpha_ks"],
                "discordant_test": config["discordant_test"],
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
                "triplet_seed": run_seed,
                "shape_diagnostics": config["shape_diagnostics"],
                "diagnostic": config["diagnostic"],
                "bootstrap": config["bootstrap"],
            }

            # The resolved count is what the run actually does with
            # ``processes: 0``, and on a scheduler allocation it is the
            # allocation, not the node.
            worker_count = resolve_worker_count(len(triplets), config["processes"])
            metrics.log(
                f"Worker processes: {worker_count} "
                f"(requested {config['processes'] or 'all'}; "
                f"{available_cpu_count()} CPU(s) available to this process)"
            )
            metrics.log("✓ Starting fused extraction + inference stage...")
            stream_start_wall, stream_start_cpu = _now_times()
            results, stream_worker_cpu = stream_triplet_results(
                triplets,
                species_triplet_trees,
                gene_trees_newick,
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

            results = _rename_result_taxa(results, rename_map)

            final_tsv = str(output_dir / "orchestrator_triplet_results.tsv")
            write_pipeline_results(
                results,
                final_tsv,
                dct_method=config["discordant_test"],
                p_value_correction=config["p_value_correction"],
                bootstrap=config["bootstrap"],
                bootstrap_debug_mode=config["bootstrap_debug_mode"],
            )
            metrics.log("✓ Fused extraction + inference complete")
            metrics.log(f"  Output: {final_tsv}")
            metrics.log(f"  Triplets analyzed: {len(results)}")
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
                # The results now carry display names while the processed
                # species tree on disk keeps the tree labels; the map lets
                # consolidation line the two up.
                map_artifacts = generate_introgression_maps(
                    results,
                    species_tree_path=species_tree_clean,
                    output_dir=str(consolidation_dir),
                    plot_taxa=(
                        None
                        if plot_taxa is None
                        else rename_taxon_labels(plot_taxa, rename_map)
                    ),
                    outgroups=rename_taxon_labels(outgroup_taxa, rename_map),
                    rename_map=rename_map,
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
