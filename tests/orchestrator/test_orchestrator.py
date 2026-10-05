"""End-to-end orchestrator tests.

``run_orchestrator`` produces per-triplet results matching values derived by hand
from the shared fixture (topology counts read off the gene-tree Newick strings,
test statistics recomputed with SciPy, Bonferroni correction applied by
definition), and all parallelization modes agree bit-for-bit (bootstrap
included, since it is deterministic per triplet under a fixed seed). See
``tests/TEST_IO.md`` for the derivation of every expected value.
"""

import argparse
from io import StringIO

import pytest
from Bio import Phylo
from scipy import stats

from ghostparser.config import InputError
from ghostparser.orchestrator.config import resolve_config
from ghostparser.orchestrator.runner import run_orchestrator
from ghostparser.orchestrator.trees import read_tree_file

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
    diagnostic=None,
    species_rename_map=None,
    species_filter=None,
):
    """Build a resolved orchestrator config for a run with a fixed bootstrap seed.

    Args:
        species_path: Path to the species tree.
        genes_path: Path to the gene trees.
        output_folder: Output directory for the run.
        processes: Worker process count.
        consolidation: Whether to enable consolidation.
        bootstrap: Whether to enable bootstrap resampling.
        diagnostic: Whether to measure every test, or ``None`` for the default.
        species_rename_map: Path to a rename map, or ``None``.
        species_filter: Path to a species filter, or ``None``.

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
        species_filter=species_filter,
        species_rename_map=species_rename_map,
        no_overwrite=None,
        processes=processes,
        alpha_dct=None,
        alpha_ks=None,
        alpha_perm=None,
        p_value_correction=None,
        diagnostic=diagnostic,
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
# non-D taxa always sit inside the ((X,Y),Z) clade with D outside, the
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


@pytest.mark.output
@pytest.mark.parametrize("diagnostic", [False, True])
def test_run_orchestrator_matches_derived_expectation(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path, diagnostic
):
    """run_orchestrator (serial) reproduces the hand-derived per-triplet expectation.

    Every triplet here stops at the count gate, so the ``diagnostic`` setting
    changes only what is measured below it: a diagnostic run still reports the
    tree-height and direction tests, the default leaves both empty and says so
    in ``perm_note``. Everything the cascade reads is the same.
    """
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        diagnostic=diagnostic,
    )
    results = run_orchestrator(config)

    assert len(results) == _N_TRIPLETS
    assert {result.triplet for result in results} == set(_EXPECTED_COUNTS)

    # The results TSV carries one row per triplet under the documented header:
    # the direction call, the gate that tells a reader whether it was
    # consulted, the equivalence p-value and the bootstrap interval.
    lines = (output_folder / "orchestrator_triplet_results.tsv").read_text().strip().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "triplet"
    for column in (
        "classification", "bootstrap_value", "perm_p_greater", "perm_p_less",
        "decision_gate", "perm_p_tost", "bootstrap_perm_stat_ci_low",
        "bootstrap_perm_stat_ci_high",
    ):
        assert column in header
    assert len(lines) - 1 == len(results)

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
        if diagnostic:
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
        assert (result.bootstrap_perm_stat_ci_low is None) == (result.bootstrap_perm_stat_ci_high is None)
        if result.bootstrap_perm_stat_ci_low is not None:
            assert result.bootstrap_perm_stat_ci_low <= result.bootstrap_perm_stat_ci_high


@pytest.mark.output
@pytest.mark.parametrize("diagnostic", [False, True])
def test_no_bootstrap_skips_the_bootstrap_and_its_columns(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path, diagnostic
):
    """Disabling the bootstrap stops the work and drops its columns.

    The studentized interval is produced by the bootstrap loop, so its absence
    is the observable proof no iterations ran. A diagnostic run measures every
    test the cascade cannot consult, which is not a licence to reinstate work
    the user switched off, so the skip holds either way. The inference columns
    and the point estimate are untouched.
    """
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        bootstrap=False,
        diagnostic=diagnostic,
    )
    results = run_orchestrator(config)

    assert len(results) == _N_TRIPLETS
    for result in results:
        assert result.bootstrap_value is None
        assert result.all_bootstrap is None
        assert result.bootstrap_perm_stat_ci_low is None
        assert result.bootstrap_perm_stat_ci_high is None
        assert result.classification == "no_introgression"

    header = (
        (output_folder / "orchestrator_triplet_results.tsv")
        .read_text()
        .splitlines()[0]
        .split("\t")
    )
    assert "bootstrap_value" not in header
    assert "all_bootstrap" not in header
    assert "classification" in header


@pytest.mark.parametrize(
    "filter_text, expected_triplets",
    [
        # Names on their own lines and sharing one, a repeat, the outgroup and
        # an unknown name: A, B and D survive, and C is left out of the run.
        ("A, B\nD\n\nOUT\nNOPE\nB\n", {("A", "B", "D")}),
        # Two usable species cannot form a triplet, so the run stops with an
        # input error.
        ("A,B,NOPE\n", None),
    ],
    ids=["three_species", "too_few"],
)
def test_species_filter_runs_every_triplet_among_the_named_species(
    orchestrator_species_tree,
    orchestrator_gene_trees,
    tmp_path,
    filter_text,
    expected_triplets,
):
    """A species filter runs exactly the triplets its usable species can form.

    The names are matched against the pruned species tree's ingroup: the
    outgroup and an unknown name are skipped rather than failing the run, a
    repeated name counts once, and fewer than three survivors is an input
    error, since the run could produce no triplet. The triplets that do run carry the same counts
    as in the unfiltered run, since the filter changes which triplets are
    measured and nothing about how.
    """
    species_filter = tmp_path / "species.txt"
    species_filter.write_text(filter_text)
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "out",
        processes=1,
        species_filter=str(species_filter),
    )
    if expected_triplets is None:
        with pytest.raises(InputError, match="at least 3 are needed"):
            run_orchestrator(config)
        return
    results = run_orchestrator(config)
    assert {result.triplet for result in results} == expected_triplets
    for result in results:
        n_con, n_dis1, n_dis2, species_topology = _EXPECTED_COUNTS[result.triplet]
        assert (result.n_con, result.n_dis1, result.n_dis2) == (n_con, n_dis1, n_dis2)
        assert result.species_tree == species_topology


@pytest.mark.output
def test_species_rename_map_reaches_every_output(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """Mapped taxa appear under their display names in every output, and nowhere else.

    The run works in the trees' own labels and the map is applied to the
    results just before they are written, so the display names (chosen here
    to hold spaces and a dot, which a bare Newick label cannot) must reach
    the returned results, the results TSV, the ``species_tree`` column and the
    consolidation artifacts, while the processed trees keep the tree labels.
    """
    rename_path = tmp_path / "names.tsv"
    rename_path.write_text("A\tHomo sapiens\nB\tPan sp.\n")
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
    assert seen == {"Homo sapiens", "Pan sp.", "C", "D"}

    # The species subtree column is renamed as Newick: the names that need it
    # are quoted, so the column still parses to the display names.
    by_triplet = {result.triplet: result for result in results}
    assert by_triplet[("Homo sapiens", "Pan sp.", "C")].species_tree == (
        "(('Homo sapiens','Pan sp.'),C);"
    )
    for result in results:
        parsed = Phylo.read(StringIO(result.species_tree), "newick")
        assert {leaf.name for leaf in parsed.get_terminals()} == set(result.triplet)

    results_tsv = (output_folder / "orchestrator_triplet_results.tsv").read_text()
    assert "Homo sapiens" in results_tsv and "Pan sp." in results_tsv

    # The processed trees are read back by the run itself, so they stay in the
    # labels the input trees use.
    for name in ("processed_species.tree", "processed_genes.tree"):
        for tree in read_tree_file(str(output_folder / name)):
            assert {t.name for t in tree.get_terminals()} <= {"A", "B", "C", "D", "OUT"}

    # The taxa-order file is what labels the heatmap axes and the bar chart, so
    # it standing in display names is the evidence the plots do too.
    consolidation_data = output_folder / "consolidation" / "consolidation_data"
    taxa_order = (consolidation_data / "introgression_taxa_order.tsv").read_text()
    assert "Homo sapiens" in taxa_order and "Pan sp." in taxa_order
    assert (output_folder / "consolidation" / "introgression_combined.png").exists()

    matrix = (consolidation_data / "introgression_matrix_inflow_outflow.tsv").read_text()
    assert "Homo sapiens" in matrix and "Pan sp." in matrix


@pytest.mark.output
def test_run_outputs_follow_the_settings(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """One run with every optional output on writes each where documented.

    Consolidation must not wipe the run folder's primary outputs and lands in
    its own subfolder; ``generate_summary_stats`` writes the 63 metric columns
    and fills the per-triplet metric statistics; ``bootstrap_diagnostic`` adds
    the per-iteration columns and populates them.
    """
    output_folder = tmp_path / "out"
    config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        output_folder,
        processes=1,
        consolidation=True,
    )
    config["generate_summary_stats"] = True
    config["bootstrap_diagnostic"] = True
    results = run_orchestrator(config)
    assert results

    for name in (
        "orchestrator_triplet_results.tsv",
        "metrics.txt",
        "processed_species.tree",
        "processed_genes.tree",
    ):
        assert (output_folder / name).exists(), name
    assert (output_folder / "consolidation").is_dir()
    assert any((output_folder / "consolidation").iterdir())

    summary_header = (
        (output_folder / "summary_statistics.tsv").read_text().splitlines()[0].split("\t")
    )
    assert "concordant_avg_tree_height_mean" in summary_header
    assert "discordant2_sister_distance_max" in summary_header
    metric_columns = [
        column
        for column in summary_header
        if column.startswith(("concordant_", "discordant1_", "discordant2_"))
    ]
    assert len(metric_columns) == 63
    assert any(result.topology_metric_statistics for result in results)

    results_header = (
        (output_folder / "orchestrator_triplet_results.tsv")
        .read_text()
        .splitlines()[0]
        .split("\t")
    )
    for column in (
        "bootstrap_dct_stats",
        "bootstrap_dct_p_value",
        "bootstrap_ks_stats",
        "bootstrap_ks_p_value",
        "bootstrap_perm_stats",
        "bootstrap_perm_p_greater",
        "bootstrap_perm_p_less",
        "bootstrap_perm_decisions",
        "bootstrap_gene_tree_heights",
    ):
        assert column in results_header
    assert any(result.bootstrap_dct_stats is not None for result in results)
    assert any(result.bootstrap_perm_decisions is not None for result in results)


def test_parallel_runs_match_serial(
    orchestrator_species_tree, orchestrator_gene_trees, tmp_path
):
    """Parallel runs yield results identical to the serial run, at any worker count."""
    serial_config = _make_config(
        orchestrator_species_tree,
        orchestrator_gene_trees,
        tmp_path / "serial",
        processes=1,
    )
    serial = {r.triplet: r for r in run_orchestrator(serial_config)}

    for processes in (2, 4):
        parallel_config = _make_config(
            orchestrator_species_tree,
            orchestrator_gene_trees,
            tmp_path / f"parallel_{processes}",
            processes=processes,
        )
        parallel = run_orchestrator(parallel_config)

        assert len(parallel) == len(serial)
        for result in parallel:
            _assert_result_matches(result, serial[result.triplet], _ALL_FIELDS)
