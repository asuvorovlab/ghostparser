"""End-to-end pipeline tests.

``run_pipeline`` reproduces the orchestrator's per-triplet inference (built from
tree_parser + triplet_processor's parquet-style observation path) on every
non-bootstrap field, and all parallelization modes agree bit-for-bit (bootstrap
included, since it is deterministic per triplet under a fixed seed).
"""

import argparse

import dendropy
import pytest

from ghostparser import triplet_processor as tp
from ghostparser import tree_parser as tpz
from ghostparser.pipeline.config import resolve_config
from ghostparser.pipeline.runner import run_pipeline

_OUTGROUP = ["OUT"]
_SEED = 20240724
_ITERATIONS = 40

# Fields compared between pipeline runs (bootstrap is deterministic under a seed).
_ALL_FIELDS = (
    "species_tree",
    "most_frequent_matches_concordant",
    "n_con",
    "n_dis1",
    "n_dis2",
    "dis1_topology",
    "dct_statistic",
    "dct_p_value",
    "dct_p_value_corrected",
    "dct_significant",
    "ks_statistic",
    "ks_p_value",
    "ks_p_value_corrected",
    "ks_significant",
    "summary_con",
    "summary_dis",
    "classification",
    "analyzed_trees",
    "bootstrap_value",
    "all_bootstrap",
)
# Fields compared against the orchestrator reference: bootstrap is excluded
# because the pipeline resamples with NumPy (statistically equivalent, not
# identical to triplet_processor's random.Random).
_NON_BOOTSTRAP_FIELDS = _ALL_FIELDS[:-2]


def _make_config(
    species_path, genes_path, output_folder, *, mode, processes, consolidation=False
):
    """Build a resolved pipeline config for a run with a fixed bootstrap seed.

    Args:
        species_path: Path to the species tree.
        genes_path: Path to the gene trees.
        output_folder: Output directory for the run.
        mode: Parallelization mode.
        processes: Worker process count.
        consolidation: Whether to enable consolidation.

    Returns:
        The resolved config dict with a fixed bootstrap seed and iterations.
    """
    args = argparse.Namespace(
        config_file=None,
        species_tree_path=str(species_path),
        gene_trees_path=str(genes_path),
        outgroups="OUT",
        output_folder=str(output_folder),
        triplet_filter=None,
        no_overwrite=None,
        processes=processes,
        parallelization_mode=mode,
        alpha_dct=None,
        alpha_ks=None,
        p_value_correction=None,
        summary_statistic=None,
        consolidation=consolidation,
        bootstrap=True,
    )
    config = resolve_config(args)
    config["bootstrap_seed"] = _SEED
    config["bootstrap_iterations"] = _ITERATIONS
    return config


def _reference_results(species_path, genes_path, tmp_path):
    """Reproduce per-triplet results using tree_parser and triplet_processor.

    Args:
        species_path: Path to the species tree.
        genes_path: Path to the gene trees.
        tmp_path: Pytest temporary directory for intermediate files.

    Returns:
        A dict mapping each triplet to its reference ``TripletPipelineResult``.
    """
    out_s = tmp_path / "ref_species.tree"
    out_g = tmp_path / "ref_genes.tree"

    tpz.clean_and_save_trees(str(species_path), str(out_s), min_avg_support=0.5)
    species_trees = tpz.read_tree_file(str(out_s))
    pruned, _excluded, _missing, ingroup = tpz._root_tree_on_outgroup(
        species_trees[0], _OUTGROUP
    )
    sp_newick = tpz.format_newick_with_precision(pruned)
    dendro = dendropy.Tree.get(
        data=sp_newick, schema="newick", preserve_underscores=True
    )
    raw_triplets = tpz.generate_triplets(sorted(ingroup), [])
    triplets, species_map, _skipped = tpz._build_species_triplet_metadata(
        dendro, raw_triplets
    )

    tpz.clean_and_save_gene_trees(
        str(genes_path), str(out_g), _OUTGROUP, min_avg_support=0.5
    )
    gene_newick = tpz._read_gene_trees_file(str(out_g))

    reference = []
    for triplet in triplets:
        # Build observation rows the way the orchestrator's parquet path does:
        # metrics computed from the extracted subtree object (unrounded), which
        # is what the pipeline now uses instead of a serialize-then-reparse.
        rows = []
        for newick in gene_newick:
            tree = dendropy.Tree.get(
                data=newick, schema="newick", preserve_underscores=True
            )
            taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
            if set(triplet).issubset(taxa):
                subtree = tpz.extract_triplet_subtree(tree, triplet)
                if subtree:
                    rows.append(tpz._build_observation_metrics(subtree, triplet))
        entry = {"species_tree": species_map[triplet], "observation_rows": rows}
        reference.append(
            tp.analyze_triplet_entry(
                triplet,
                entry,
                generate_summary_stats=False,
                bootstrap=True,
                bootstrap_options={"iterations": _ITERATIONS},
                triplet_seed=_SEED,
            )
        )
    reference = tp._apply_triplet_result_p_value_correction(
        reference,
        alpha_dct=0.05,
        alpha_ks=0.05,
        method="no",
        stats_backend="standard",
    )
    return {result.triplet: result for result in reference}


def _assert_result_matches(pipeline_result, reference_result, fields):
    """Assert one result agrees with a reference on the given fields.

    Args:
        pipeline_result: The pipeline's ``TripletPipelineResult``.
        reference_result: The reference ``TripletPipelineResult``.
        fields: Iterable of field names to compare.
    """
    for field in fields:
        pipeline_value = getattr(pipeline_result, field)
        reference_value = getattr(reference_result, field)
        if isinstance(reference_value, float):
            assert pipeline_value == pytest.approx(reference_value), field
        else:
            assert pipeline_value == reference_value, field


def test_run_pipeline_matches_reference(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """run_pipeline (serial) reproduces the tree_parser + triplet_processor reference."""
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        tmp_path / "out",
        mode="taxon",
        processes=1,
    )
    results = run_pipeline(config)
    reference = _reference_results(pipeline_species_tree, pipeline_gene_trees, tmp_path)

    assert results
    assert len(results) == len(reference)
    for result in results:
        assert result.triplet in reference
        _assert_result_matches(
            result, reference[result.triplet], _NON_BOOTSTRAP_FIELDS
        )
        assert 0.0 <= result.bootstrap_value <= 1.0
        assert sum(result.all_bootstrap.values()) == pytest.approx(1.0)


def test_run_pipeline_writes_results_tsv(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """run_pipeline writes pipeline_triplet_results.tsv with the expected header and row count."""
    output_folder = tmp_path / "out"
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        output_folder,
        mode="taxon",
        processes=1,
    )
    results = run_pipeline(config)

    tsv_path = output_folder / "pipeline_triplet_results.tsv"
    assert tsv_path.exists()
    lines = tsv_path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "triplet"
    assert "classification" in header
    assert "bootstrap_value" in header
    assert len(lines) - 1 == len(results)


def test_consolidation_preserves_run_outputs(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """Consolidation preserves the run outputs and writes its artifacts to a subfolder."""
    output_folder = tmp_path / "out"
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        output_folder,
        mode="taxon",
        processes=1,
        consolidation=True,
    )
    results = run_pipeline(config)
    assert results

    # Consolidation must not wipe the run folder's primary outputs.
    assert (output_folder / "pipeline_triplet_results.tsv").exists()
    assert (output_folder / "metrics.txt").exists()
    assert (output_folder / "processed_species.tree").exists()
    assert (output_folder / "processed_genes.tree").exists()
    # Consolidation artifacts land in their own subfolder.
    assert (output_folder / "consolidation").is_dir()
    assert any((output_folder / "consolidation").iterdir())


def test_generate_summary_stats_writes_tsv(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """generate_summary_stats writes summary_statistics.tsv with the 63 metric columns."""
    output_folder = tmp_path / "out"
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        output_folder,
        mode="taxon",
        processes=1,
    )
    config["generate_summary_stats"] = True
    results = run_pipeline(config)

    summary_tsv = output_folder / "summary_statistics.tsv"
    assert summary_tsv.exists()
    header = summary_tsv.read_text().splitlines()[0].split("\t")
    assert "concordant_avg_tree_height_mean" in header
    assert "discordant2_sister_distance_max" in header
    metric_columns = [
        column
        for column in header
        if column.startswith(("concordant_", "discordant1_", "discordant2_"))
    ]
    assert len(metric_columns) == 63
    # The per-triplet metric statistics are populated on the results too.
    assert any(result.topology_metric_statistics for result in results)


def test_bootstrap_debug_mode_writes_debug_columns(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """bootstrap_debug_mode adds the bootstrap-debug columns and populates them."""
    output_folder = tmp_path / "out"
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        output_folder,
        mode="taxon",
        processes=1,
    )
    config["bootstrap_debug_mode"] = True
    results = run_pipeline(config)

    tsv_path = output_folder / "pipeline_triplet_results.tsv"
    header = tsv_path.read_text().splitlines()[0].split("\t")
    for column in (
        "bootstrap_dct_stats",
        "bootstrap_dct_p_value",
        "bootstrap_ks_stats",
        "bootstrap_ks_p_value",
        "bootstrap_gene_tree_heights",
    ):
        assert column in header
    assert any(result.bootstrap_dct_stats is not None for result in results)


@pytest.mark.parametrize(
    "mode,processes",
    [("taxon", 2), ("gene", 2)],
)
def test_parallel_modes_match_serial(
    pipeline_species_tree, pipeline_gene_trees, tmp_path, mode, processes
):
    """taxon and gene parallel modes yield identical results to the serial run."""
    serial_config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        tmp_path / "serial",
        mode="taxon",
        processes=1,
    )
    serial = {r.triplet: r for r in run_pipeline(serial_config)}

    parallel_config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        tmp_path / "parallel",
        mode=mode,
        processes=processes,
    )
    parallel = run_pipeline(parallel_config)

    assert len(parallel) == len(serial)
    for result in parallel:
        _assert_result_matches(result, serial[result.triplet], _ALL_FIELDS)
