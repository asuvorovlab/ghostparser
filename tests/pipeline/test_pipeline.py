"""End-to-end pipeline tests.

``run_pipeline`` produces per-triplet results matching values derived by hand
from the shared fixture (topology counts read off the gene-tree Newick strings,
test statistics recomputed with SciPy, Bonferroni correction applied by
definition), and all parallelization modes agree bit-for-bit (bootstrap
included, since it is deterministic per triplet under a fixed seed). See
``tests/TEST_IO.md`` for the derivation of every expected value.
"""

import argparse

import pytest
from scipy import stats

from ghostparser.pipeline.config import resolve_config
from ghostparser.pipeline.runner import run_pipeline

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


def _make_config(
    species_path,
    genes_path,
    output_folder,
    *,
    mode,
    processes,
    consolidation=False,
    bootstrap=True,
):
    """Build a resolved pipeline config for a run with a fixed bootstrap seed.

    Args:
        species_path: Path to the species tree.
        genes_path: Path to the gene trees.
        output_folder: Output directory for the run.
        mode: Parallelization mode.
        processes: Worker process count.
        consolidation: Whether to enable consolidation.
        bootstrap: Whether to enable bootstrap resampling.

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
        bootstrap=bootstrap,
    )
    config = resolve_config(args)
    config["bootstrap_seed"] = _SEED
    config["bootstrap_iterations"] = _ITERATIONS
    return config


# Hand-derived topology counts for the shared fixture.
#
# The pruned species tree is (((A,B),C),D), giving 4 triplets. Every gene tree
# has the shape ((((X,Y),Z),D),OUT), so for any triplet containing D the two
# non-D taxa always sit inside the ((X,Y),Z) clade with D outside -- the
# concordant topology is the only one observed (12/0/0). Only triplet (A,B,C)
# varies: reading the innermost sister pair off each of the 12 gene trees gives
# (A,B) in trees 0,1,4,6,8,10,11 -> 7 concordant; (B,C) in trees 3,5,9 -> 3; and
# (A,C) in trees 2,7 -> 2. Ranking the discordant pair puts BC first (3 >= 2).
_EXPECTED_COUNTS = {
    ("A", "B", "C"): (7, 3, 2, "((A,B),C);"),
    ("A", "B", "D"): (12, 0, 0, "((A,B),D);"),
    ("A", "C", "D"): (12, 0, 0, "((A,C),D);"),
    ("B", "C", "D"): (12, 0, 0, "((B,C),D);"),
}
_N_TRIPLETS = len(_EXPECTED_COUNTS)


def _expected_dct(n_dis1, n_dis2):
    """Compute the expected chi-square DCT statistic and p-value.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.

    Returns:
        A tuple ``(statistic, p_value)``; an all-zero pair is a no-op test.
    """
    if n_dis1 + n_dis2 == 0:
        return 0.0, 1.0
    result = stats.chisquare([n_dis1, n_dis2])
    return float(result.statistic), float(result.pvalue)


def _bonferroni(p_value):
    """Apply the pipeline's default Bonferroni correction across all triplets.

    Args:
        p_value: The uncorrected p-value.

    Returns:
        ``min(1.0, p_value * number_of_triplets)``.
    """
    return min(1.0, p_value * _N_TRIPLETS)


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


def test_run_pipeline_matches_derived_expectation(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """run_pipeline (serial) reproduces the hand-derived per-triplet expectation."""
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        tmp_path / "out",
        mode="taxon",
        processes=1,
    )
    results = run_pipeline(config)

    assert len(results) == _N_TRIPLETS
    assert {result.triplet for result in results} == set(_EXPECTED_COUNTS)

    for result in results:
        n_con, n_dis1, n_dis2, species_topology = _EXPECTED_COUNTS[result.triplet]
        assert result.n_con == n_con
        assert result.n_dis1 == n_dis1
        assert result.n_dis2 == n_dis2
        assert result.species_tree == species_topology
        assert result.analyzed_trees == 12
        assert result.most_frequent_matches_concordant is True

        dct_statistic, dct_p_value = _expected_dct(n_dis1, n_dis2)
        assert result.dct_statistic == pytest.approx(dct_statistic)
        assert result.dct_p_value == pytest.approx(dct_p_value)
        assert result.dct_p_value_corrected == pytest.approx(_bonferroni(dct_p_value))
        # Every triplet's corrected DCT p-value stays well above alpha (0.05),
        # so the decision logic stops at the first gate for all of them.
        assert result.dct_significant is False
        assert result.ks_p_value_corrected == pytest.approx(
            _bonferroni(result.ks_p_value)
        )
        assert result.classification == "no_introgression"

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


def test_no_bootstrap_omits_the_bootstrap_columns(
    pipeline_species_tree, pipeline_gene_trees, tmp_path
):
    """Disabling bootstrap drops its TSV columns while keeping the inference fields."""
    output_folder = tmp_path / "out"
    config = _make_config(
        pipeline_species_tree,
        pipeline_gene_trees,
        output_folder,
        mode="taxon",
        processes=1,
        bootstrap=False,
    )
    results = run_pipeline(config)

    header = (
        (output_folder / "pipeline_triplet_results.tsv")
        .read_text()
        .splitlines()[0]
        .split("\t")
    )
    assert "bootstrap_value" not in header
    assert "all_bootstrap" not in header
    # The inference columns are still present and the triplets still resolved.
    assert "classification" in header
    assert len(results) == _N_TRIPLETS


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
