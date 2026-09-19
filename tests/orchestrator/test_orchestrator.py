"""End-to-end orchestrator tests.

``run_orchestrator`` produces per-triplet results matching values derived by hand
from the shared fixture (topology counts read off the gene-tree Newick strings,
test statistics recomputed with SciPy, Bonferroni correction applied by
definition), and all parallelization modes agree bit-for-bit (bootstrap
included, since it is deterministic per triplet under a fixed seed). See
``tests/TEST_IO.md`` for the derivation of every expected value.
"""

import argparse

import pytest
from scipy import stats

from ghostparser.orchestrator.config import resolve_config
from ghostparser.orchestrator.runner import run_orchestrator

pytestmark = pytest.mark.integration

_SEED = 20240724
_ITERATIONS = 40

# Fields compared between orchestrator runs (bootstrap is deterministic under a seed).
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
    "perm_decision",
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
    processes,
    consolidation=False,
    bootstrap=True,
    pipeline_mode=None,
    species_rename_map=None,
):
    """Build a resolved orchestrator config for a run with a fixed bootstrap seed.

    Args:
        species_path: Path to the species tree.
        genes_path: Path to the gene trees.
        output_folder: Output directory for the run.
        processes: Worker process count.
        consolidation: Whether to enable consolidation.
        bootstrap: Whether to enable bootstrap resampling.
        pipeline_mode: ``efficient``/``detailed``, or ``None`` for the default.
        species_rename_map: Path to a rename map, or ``None``.

    Returns:
        The resolved config dict with a fixed bootstrap seed and iterations.
    """
    args = argparse.Namespace(
        config_file=None,
        species_tree_path=str(species_path),
        gene_trees_path=str(genes_path),
        outgroup="OUT",
        output_folder=str(output_folder),
        triplet_filter=None,
        species_rename_map=species_rename_map,
        no_overwrite=None,
        processes=processes,
        alpha_dct=None,
        alpha_ks=None,
        alpha_perm=None,
        p_value_correction=None,
        pipeline_mode=pipeline_mode,
        consolidation=consolidation,
        bootstrap=bootstrap,
    )
    config = resolve_config(args)
    config["seed"] = _SEED
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
    """Apply the orchestrator's default Bonferroni correction across all triplets.

    Args:
        p_value: The uncorrected p-value.

    Returns:
        ``min(1.0, p_value * number_of_triplets)``.
    """
    return min(1.0, p_value * _N_TRIPLETS)


def _assert_result_matches(orchestrator_result, reference_result, fields):
    """Assert one result agrees with a reference on the given fields.

    Args:
        orchestrator_result: The orchestrator's ``TripletPipelineResult``.
        reference_result: The reference ``TripletPipelineResult``.
        fields: Iterable of field names to compare.
    """
    for field in fields:
        orchestrator_value = getattr(orchestrator_result, field)
        reference_value = getattr(reference_result, field)
        if isinstance(reference_value, float):
            assert orchestrator_value == pytest.approx(reference_value), field
        else:
            assert orchestrator_value == reference_value, field


@pytest.mark.parametrize("pipeline_mode", ["efficient", "detailed"])
def test_run_orchestrator_matches_derived_expectation(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path, pipeline_mode
):
    """run_orchestrator (serial) reproduces the hand-derived per-triplet expectation.

    Every triplet here stops at the count gate, so the two pipeline modes differ
    only in what they measure below it: the detailed mode still reports the
    tree-height and direction tests, the efficient mode leaves both empty and
    says so in ``perm_note``. Everything the cascade reads is the same.
    """
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "out",
        processes=1,
        pipeline_mode=pipeline_mode,
    )
    results = run_orchestrator(config)

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
        # Set by the run-wide correction pass, which recomputes the gate from the
        # corrected significance alongside the classification.
        assert result.decision_gate == "DCT"
        assert result.classification == "no_introgression"
        if pipeline_mode == "detailed":
            # The count gate settles every triplet, but the tree-height test is
            # still measured and corrected across the same whole-run family,
            # and the direction test is reported even though nothing read it.
            assert 0.0 <= result.ks_p_value <= 1.0
            assert result.ks_p_value_corrected == pytest.approx(
                _bonferroni(result.ks_p_value)
            )
            assert result.ks_significant is (result.ks_p_value_corrected <= 0.05)
            assert result.perm_decision is not None
        else:
            # Under the default bfn the family is fixed, so a value below a
            # settled count gate would go unread and is never measured.
            assert result.ks_p_value is None
            assert result.ks_p_value_corrected is None
            assert result.ks_significant is None
            assert result.perm_decision is None
            assert result.perm_note == "direction_test_not_consulted"

        assert 0.0 <= result.bootstrap_value <= 1.0
        assert sum(result.all_bootstrap.values()) == pytest.approx(1.0)
        # The interval is reported only when some resample yielded two
        # observations in both groups; either way its two bounds agree on
        # whether they exist, and an existing pair is ordered.
        assert (result.bootstrap_stat_ci_low is None) == (result.bootstrap_stat_ci_high is None)
        if result.bootstrap_stat_ci_low is not None:
            assert result.bootstrap_stat_ci_low <= result.bootstrap_stat_ci_high


@pytest.mark.output
def test_run_orchestrator_writes_results_tsv(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """run_orchestrator writes orchestrator_triplet_results.tsv with the expected header and row count."""
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
    )
    results = run_orchestrator(config)

    tsv_path = output_folder / "orchestrator_triplet_results.tsv"
    assert tsv_path.exists()
    lines = tsv_path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "triplet"
    assert "classification" in header
    assert "bootstrap_value" in header
    # The permutation columns that carry the direction call.
    assert "perm_p_greater" in header
    assert "perm_p_less" in header
    # decision_gate is what tells a reader whether perm_decision was consulted.
    assert "decision_gate" in header
    # The equivalence p-value and the interval on the studentized difference.
    assert "perm_p_tost" in header
    assert "bootstrap_stat_ci_low" in header
    assert "bootstrap_stat_ci_high" in header
    assert len(lines) - 1 == len(results)


@pytest.mark.output
def test_no_bootstrap_omits_the_bootstrap_columns(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """Disabling bootstrap drops its TSV columns while keeping the inference fields."""
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        bootstrap=False,
    )
    results = run_orchestrator(config)

    header = (
        (output_folder / "orchestrator_triplet_results.tsv")
        .read_text()
        .splitlines()[0]
        .split("\t")
    )
    assert "bootstrap_value" not in header
    assert "all_bootstrap" not in header
    # The inference columns are still present and the triplets still resolved.
    assert "classification" in header
    assert len(results) == _N_TRIPLETS


@pytest.mark.parametrize("pipeline_mode", ["efficient", "detailed"])
def test_no_bootstrap_skips_the_bootstrap_itself(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path, pipeline_mode
):
    """Disabling the bootstrap stops the work, not just the columns.

    The studentized interval is produced by the bootstrap loop, so its absence
    is the observable proof no iterations ran. The detailed mode declines to
    skip work the cascade cannot consult, which is not a licence to reinstate
    work the user switched off, so the skip holds in both modes.
    """
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "out",
        processes=1,
        bootstrap=False,
        pipeline_mode=pipeline_mode,
    )
    results = run_orchestrator(config)

    assert len(results) == _N_TRIPLETS
    for result in results:
        assert result.bootstrap_value is None
        assert result.all_bootstrap is None
        assert result.bootstrap_stat_ci_low is None
        assert result.bootstrap_stat_ci_high is None
        # The point estimate is unaffected by the bootstrap being off.
        assert result.classification == "no_introgression"


@pytest.mark.output
def test_species_rename_map_reaches_every_output(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """Mapped taxa appear under their display names in every named output.

    The rename is applied as the trees are read, so it has to reach the
    in-memory results, the results TSV, the per-triplet species subtree, the
    processed tree files, and the consolidation artifacts alike.
    """
    rename_path = tmp_path / "names.tsv"
    rename_path.write_text("A\tHomo\nB\tPan\n")
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        consolidation=True,
        species_rename_map=str(rename_path),
    )

    results = run_orchestrator(config)
    assert len(results) == _N_TRIPLETS

    # Mapped taxa are renamed; unmapped ones (C, D) are untouched.
    seen = {taxon for result in results for taxon in result.triplet}
    assert seen == {"Homo", "Pan", "C", "D"}

    # The per-triplet species subtree is rebuilt from the renamed tree.
    assert any("Homo" in (result.species_tree or "") for result in results)

    results_tsv = (output_folder / "orchestrator_triplet_results.tsv").read_text()
    assert "Homo" in results_tsv and "Pan" in results_tsv

    # The cleaned trees written alongside the results carry the display names,
    # which is what consolidation reads back for its taxon ordering.
    assert "Homo" in (output_folder / "processed_species.tree").read_text()
    assert "Homo" in (output_folder / "processed_genes.tree").read_text()

    # The taxa-order file is what labels the heatmap axes and the bar chart, so
    # it standing in display names is the evidence the plots do too.
    taxa_order = (
        output_folder
        / "consolidation"
        / "consolidation_data"
        / "introgression_taxa_order.tsv"
    ).read_text()
    assert "Homo" in taxa_order and "Pan" in taxa_order
    assert (output_folder / "consolidation" / "introgression_combined.png").exists()

    matrix = (
        output_folder
        / "consolidation"
        / "consolidation_data"
        / "introgression_matrix_inflow_outflow.tsv"
    ).read_text()
    assert "Homo" in matrix and "Pan" in matrix


@pytest.mark.output
def test_consolidation_preserves_run_outputs(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """Consolidation preserves the run outputs and writes its artifacts to a subfolder."""
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        consolidation=True,
    )
    results = run_orchestrator(config)
    assert results

    # Consolidation must not wipe the run folder's primary outputs.
    assert (output_folder / "orchestrator_triplet_results.tsv").exists()
    assert (output_folder / "metrics.txt").exists()
    assert (output_folder / "processed_species.tree").exists()
    assert (output_folder / "processed_genes.tree").exists()
    # Consolidation artifacts land in their own subfolder.
    assert (output_folder / "consolidation").is_dir()
    assert any((output_folder / "consolidation").iterdir())


@pytest.mark.output
def test_generate_summary_stats_writes_tsv(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """generate_summary_stats writes summary_statistics.tsv with the 63 metric columns."""
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
    )
    config["generate_summary_stats"] = True
    results = run_orchestrator(config)

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


@pytest.mark.output
def test_bootstrap_debug_mode_writes_debug_columns(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """bootstrap_debug_mode adds the bootstrap-debug columns and populates them."""
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
    )
    config["bootstrap_debug_mode"] = True
    results = run_orchestrator(config)

    tsv_path = output_folder / "orchestrator_triplet_results.tsv"
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


@pytest.mark.parametrize("processes", [2, 4])
def test_parallel_runs_match_serial(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path, processes
):
    """Parallel runs yield identical results to the serial run."""
    serial_config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "serial",
        processes=1,
    )
    serial = {r.triplet: r for r in run_orchestrator(serial_config)}

    parallel_config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "parallel",
        processes=processes,
    )
    parallel = run_orchestrator(parallel_config)

    assert len(parallel) == len(serial)
    for result in parallel:
        _assert_result_matches(result, serial[result.triplet], _ALL_FIELDS)
