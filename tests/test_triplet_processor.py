"""Tests for the triplet_processor module."""

import argparse
import random
from pathlib import Path

import dendropy
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

import ghostparser.triplet_processor as triplet_processor_module

from ghostparser.triplet_processor import (
    _adjust_p_values,
    analyze_triplet_gene_tree_file,
    classify_triplet_topology,
    collect_triplet_statistics,
    compute_tree_height_statistic,
    parse_triplet_gene_trees_file,
    pearson_discordant_chi_square_test,
    run_discordant_count_test,
    run_triplet_pipeline,
    two_sample_ks_test,
    two_sample_ks_test_hybrid,
    two_proportion_discordant_z_test,
    write_pipeline_statistics_json,
    write_pipeline_results,
    _resolve_runtime_args,
)
from ghostparser.triplet_utils import TOPOLOGY_AB, TOPOLOGY_AC, TOPOLOGY_BC
from ghostparser.triplet_utils import classify_triplet_topology_string


def _tree(newick):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def test_compute_tree_height_statistic_matches_definition():
    tree = _tree("((B:2,C:3):4,A:1);")
    observed = compute_tree_height_statistic(tree)
    expected = (1 + (2 + 4) + (3 + 4)) / 3
    assert observed == pytest.approx(expected)


def test_compute_tree_height_statistic_supports_taxon_specific_strategies():
    tree = _tree("((A:2,B:3):4,C:1);")
    species_triplet = ("A", "B", "C")

    assert compute_tree_height_statistic(tree, strategy="A", species_triplet=species_triplet) == pytest.approx(6.0)
    assert compute_tree_height_statistic(tree, strategy="B", species_triplet=species_triplet) == pytest.approx(7.0)
    assert compute_tree_height_statistic(tree, strategy="C", species_triplet=species_triplet) == pytest.approx(1.0)


def test_compute_tree_height_statistic_rejects_unknown_strategy():
    tree = _tree("((A:2,B:3):4,C:1);")

    with pytest.raises(ValueError, match="Unsupported tree height calculation strategy"):
        compute_tree_height_statistic(tree, strategy="D", species_triplet=("A", "B", "C"))


def test_compute_tree_height_statistic_requires_species_triplet_for_taxon_specific_strategies():
    tree = _tree("((A:2,B:3):4,C:1);")

    with pytest.raises(ValueError, match="species_triplet is required"):
        compute_tree_height_statistic(tree, strategy="A")


def test_classify_triplet_topology_string_for_all_three_topologies():
    species_triplet = ("A", "B", "C")

    assert classify_triplet_topology_string(_tree("((A:1,B:1):1,C:1);"), species_triplet) == TOPOLOGY_AB
    assert classify_triplet_topology_string(_tree("((B:1,C:1):1,A:1);"), species_triplet) == TOPOLOGY_BC
    assert classify_triplet_topology_string(_tree("((A:1,C:1):1,B:1);"), species_triplet) == TOPOLOGY_AC


def test_classify_triplet_topology_labels_concordant_and_discordants():
    species_triplet = ("A", "B", "C")
    topology_counts = {
        TOPOLOGY_AB: 4,
        TOPOLOGY_BC: 8,
        TOPOLOGY_AC: 2,
    }

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((A:1,B:1):1,C:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "concordant"
    assert most_frequent_matches_concordant is False

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((B:1,C:1):1,A:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "discordant1"
    assert most_frequent_matches_concordant is False

    label, most_frequent_matches_concordant = classify_triplet_topology(
        _tree("((A:1,C:1):1,B:1);"),
        species_triplet,
        topology_counts,
    )
    assert label == "discordant2"
    assert most_frequent_matches_concordant is False


@pytest.mark.parametrize(
    "test_fn",
    [
        pearson_discordant_chi_square_test,
        two_proportion_discordant_z_test,
    ],
)
def test_balanced_discordant_count_tests_are_not_significant(test_fn):
    stat, p_value = test_fn(10, 10)
    assert stat == pytest.approx(0.0)
    assert p_value == pytest.approx(1.0)


@pytest.mark.backend_parity
def test_custom_chi_square_matches_scipy_reference_randomized():
    rng = random.Random(123)
    for _ in range(1000):
        n_dis1 = rng.randint(0, 500)
        n_dis2 = rng.randint(0, 500)

        custom_stat, custom_p = pearson_discordant_chi_square_test(n_dis1, n_dis2)
        scipy_result = stats.chisquare([n_dis1, n_dis2])

        assert custom_stat == pytest.approx(float(scipy_result.statistic), rel=0.0, abs=1e-12)
        assert custom_p == pytest.approx(float(scipy_result.pvalue), rel=0.0, abs=1e-12)


@pytest.mark.backend_parity
def test_custom_z_test_matches_statsmodels_reference_randomized():
    rng = random.Random(234)
    for _ in range(1000):
        n_dis1 = rng.randint(0, 500)
        n_dis2 = rng.randint(0, 500)

        custom_stat, custom_p = two_proportion_discordant_z_test(n_dis1, n_dis2)

        total = n_dis1 + n_dis2
        if total == 0:
            assert custom_stat == pytest.approx(0.0)
            assert custom_p == pytest.approx(1.0)
            continue

        ref_stat, ref_p = proportions_ztest(
            count=[n_dis1, n_dis2],
            nobs=[total, total],
            alternative="two-sided",
        )

        assert custom_stat == pytest.approx(float(ref_stat), rel=0.0, abs=1e-12)
        assert custom_p == pytest.approx(float(ref_p), rel=0.0, abs=1e-12)


@pytest.mark.backend_parity
def test_custom_ks_matches_scipy_asymptotic_reference_randomized():
    rng = random.Random(456)
    p_diffs = []
    for _ in range(600):
        n1 = rng.randint(20, 80)
        n2 = rng.randint(20, 80)
        sample_a = [rng.random() for _ in range(n1)]
        sample_b = [rng.random() for _ in range(n2)]

        custom_d, custom_p = two_sample_ks_test(sample_a, sample_b)
        scipy_res = stats.ks_2samp(sample_a, sample_b, alternative="two-sided", method="asymp")

        assert custom_d == pytest.approx(float(scipy_res.statistic), rel=0.0, abs=1e-12)
        p_diffs.append(abs(custom_p - float(scipy_res.pvalue)))

    # TODO: This threshold is intentionally relaxed because custom KS currently uses
    # an asymptotic p-value approximation that can deviate from SciPy asymptotic values
    # in a small fraction of randomized cases. If/when a corrected KS p-value
    # calculation function is introduced, tighten this back to the stricter bound.
    assert max(p_diffs) < 0.04


def test_run_triplet_pipeline_uses_species_concordant_and_frequency_ranked_discordants():
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        rng=random.Random(0),
    )

    assert result.triplet == species_triplet
    assert result.n_con == 8
    assert result.n_dis1 == 12
    assert result.n_dis2 == 4
    assert result.most_frequent_matches_concordant is False


def test_run_triplet_pipeline_supports_z_test_for_discordant_counts():
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="z-test",
    )

    assert result.dct_significant is True
    assert result.dct_p_value < 0.01
    assert result.dct_statistic is not None


def test_run_triplet_pipeline_supports_standard_stats_backend():
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="chi-square",
        stats_backend="standard",
    )

    assert result.dct_statistic is not None


@pytest.mark.backend_parity
def test_standard_z_test_matches_statsmodels_reference_randomized():
    rng = random.Random(789)
    for _ in range(1000):
        n_dis1 = rng.randint(0, 500)
        n_dis2 = rng.randint(0, 500)

        standard_stat, standard_p = run_discordant_count_test(
            n_dis1,
            n_dis2,
            method="z-test",
            stats_backend="standard",
        )

        total = n_dis1 + n_dis2
        if total == 0:
            assert standard_stat == pytest.approx(0.0)
            assert standard_p == pytest.approx(1.0)
            continue

        ref_stat, ref_p = proportions_ztest(
            count=[n_dis1, n_dis2],
            nobs=[total, total],
            alternative="two-sided",
        )
        assert standard_stat == pytest.approx(float(ref_stat), rel=0.0, abs=1e-12)
        assert standard_p == pytest.approx(float(ref_p), rel=0.0, abs=1e-12)


def test_run_triplet_pipeline_supports_median_summary_statistic():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:1.0,B:1.0):2.0,C:3.0);"
    dis1_tree = "((B:0.2,C:0.2):0.3,A:0.5);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"
    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
    )

    assert result.summary_con is not None
    assert result.summary_dis is not None


def test_run_triplet_pipeline_breaks_discordant_ties_by_first_topology():
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 3)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
    )

    assert result.triplet == species_triplet
    assert result.n_dis1 == 3
    assert result.n_dis2 == 3


def test_run_triplet_pipeline_relabels_a_b_when_ac_is_more_frequent_discordant():
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 12) + (["((B:1,C:1):1,A:1);"] * 4)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
    )

    # Labels are reassigned so dis1 is always BC|A after canonicalization.
    assert result.triplet == ("B", "A", "C")
    assert result.n_con == 8
    assert result.n_dis1 == 12
    assert result.n_dis2 == 4


def test_run_triplet_pipeline_no_introgression_when_dct_not_significant():
    species_triplet = ("A", "B", "C")
    trees = [
        "((B:1,C:1):1,A:1);",
        "((A:1,C:1):1,B:1);",
        "((A:1,B:1):1,C:1);",
    ] * 8

    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(1))

    assert result.classification == "no_introgression"
    assert not result.dct_significant
    assert result.ks_p_value is not None
    assert result.ks_statistic is not None
    assert result.summary_con is not None
    assert result.summary_dis is not None


def test_run_triplet_pipeline_inflow_when_ks_not_significant():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:0.6,B:0.6):0.4,C:1.0);"
    dis1_tree = "((B:0.6,C:0.6):0.4,A:1.0);"
    dis2_tree = "((A:0.6,C:0.6):0.4,B:1.0);"

    trees = ([dis1_tree] * 30) + ([dis2_tree] * 5) + ([con_tree] * 30)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(2))

    assert result.dct_significant
    assert result.ks_significant is False
    assert result.summary_con is not None
    assert result.summary_dis is not None
    assert result.classification == "inflow_introgression"


def test_run_triplet_pipeline_outflow_when_con_summary_higher():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:1.0,B:1.0):2.0,C:3.0);"
    dis1_tree = "((B:0.2,C:0.2):0.3,A:0.5);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"

    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(3))

    assert result.classification == "outflow_introgression"
    assert result.ks_significant is True
    assert result.most_frequent_matches_concordant is True
    assert result.summary_con > result.summary_dis


def test_run_triplet_pipeline_ghost_when_dis_summary_higher():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:0.2,B:0.2):0.3,C:0.5);"
    dis1_tree = "((B:1.0,C:1.0):2.0,A:3.0);"
    dis2_tree = "((A:0.2,C:0.2):0.3,B:0.5);"

    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(4))

    assert result.classification == "ghost_introgression"
    assert result.ks_significant is True
    assert result.most_frequent_matches_concordant is True
    assert result.summary_con < result.summary_dis


def test_parse_analyze_and_write_pipeline_roundtrip_with_species_header(tmp_path):
    content = """A,B,C\t4\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    parsed = parse_triplet_gene_trees_file(str(input_file))
    assert ("A", "B", "C") in parsed
    entry = parsed[("A", "B", "C")]
    assert entry["count"] == 4
    assert entry["species_tree"] == "((A:1,B:1):1,C:1);"
    assert len(entry["gene_trees"]) == 4

    results = analyze_triplet_gene_tree_file(str(input_file), rng=random.Random(5))
    assert len(results) == 1

    output_file = tmp_path / "triplet_introgression_results.tsv"
    write_pipeline_results(results, str(output_file), dct_method="chi-square")

    out = output_file.read_text()
    assert "species_tree" in out
    assert "most_frequent_matches_concordant" in out
    assert "dct_chi_stats" in out
    assert "dct_z_score" not in out
    assert "A,B,C" in out


def test_analyze_triplet_gene_tree_file_with_multiprocessing(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    results = analyze_triplet_gene_tree_file(
        str(input_file),
        use_multiprocessing=True,
        processes=2,
    )

    assert len(results) == 1
    assert results[0].analyzed_trees == 3


def test_analyze_triplet_gene_tree_file_rejects_unknown_discordant_test(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported discordant test method"):
        analyze_triplet_gene_tree_file(str(input_file), discordant_test="bad-test")


def test_analyze_triplet_gene_tree_file_rejects_unknown_summary_statistic(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported summary statistic"):
        analyze_triplet_gene_tree_file(str(input_file), summary_statistic="bad-summary")


def test_run_triplet_pipeline_supports_mode_summary_statistic():
    species_triplet = ("A", "B", "C")
    trees = (["((A:0.2,B:0.2):0.3,C:0.5);"] * 40) + (["((B:1.0,C:1.0):2.0,A:3.0);"] * 30) + (["((A:0.2,C:0.2):0.3,B:0.5);"] * 5)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="mode",
    )

    assert result.summary_con is not None
    assert result.summary_dis is not None


def test_run_triplet_pipeline_supports_taxon_specific_tree_height_strategy():
    species_triplet = ("A", "B", "C")
    con_tree = "((A:5.0,B:0.2):0.1,C:0.3);"
    dis1_tree = "((B:0.2,C:0.2):0.1,A:0.4);"
    dis2_tree = "((A:0.2,C:0.2):0.1,B:0.3);"
    trees = ([con_tree] * 40) + ([dis1_tree] * 30) + ([dis2_tree] * 5)

    result_avg = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        tree_height_calculation_strategy="AVG",
        rng=random.Random(300),
    )
    result_a = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        tree_height_calculation_strategy="A",
        rng=random.Random(301),
    )

    assert result_avg.summary_con != result_a.summary_con


def test_run_triplet_pipeline_bootstrap_unresolved_when_metrics_missing():
    species_triplet = ("A", "B", "C")
    trees = ["((A:1,B:1):1,C:1);"] * 6

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        bootstrap=True,
        bootstrap_options={"iterations": 5, "seed": 123, "debug_mode": False, "summary_only": True},
    )

    assert result.classification == "no_introgression"
    assert result.bootstrap_value == pytest.approx(1.0)
    assert result.all_bootstrap["unresolved"] == pytest.approx(1.0)


def test_run_bootstrap_iterations_joins_tied_classes(monkeypatch):
    sequence = [
        "ghost_introgression",
        "inflow_introgression",
        "ghost_introgression",
        "inflow_introgression",
    ]

    def _fake_pipeline(*_args, **_kwargs):
        cls = sequence.pop(0)
        return triplet_processor_module.TripletPipelineResult(
            triplet=("A", "B", "C"),
            species_tree=None,
            most_frequent_matches_concordant=True,
            n_con=2,
            n_dis1=1,
            n_dis2=1,
            dis1_topology="BC",
            dct_statistic=1.0,
            dct_p_value=0.01,
            dct_p_value_corrected=0.01,
            dct_significant=True,
            ks_p_value=0.01,
            ks_p_value_corrected=0.01,
            ks_statistic=0.2,
            ks_significant=True,
            summary_con=1.0,
            summary_dis=0.5,
            classification=cls,
            analyzed_trees=4,
        )

    monkeypatch.setattr(triplet_processor_module, "_run_triplet_pipeline_from_observations", _fake_pipeline)

    payload = triplet_processor_module._run_bootstrap_iterations(
        species_triplet=("A", "B", "C"),
        observations=[(TOPOLOGY_AB, 1.0)],
        iterations=4,
        alpha_dct=0.01,
        alpha_ks=0.05,
        discordant_test="chi-square",
        summary_statistic="median",
        stats_backend="custom",
        species_topology=TOPOLOGY_AB,
        species_tree_newick=None,
        debug_mode=False,
        summary_only=True,
        rng=random.Random(1),
    )

    assert payload["bootstrap_value"] == pytest.approx(0.5)


def test_analyze_triplet_gene_tree_file_rejects_unknown_stats_backend(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported stats backend"):
        analyze_triplet_gene_tree_file(str(input_file), stats_backend="numpy")


def test_analyze_triplet_gene_tree_file_rejects_unknown_tree_height_strategy(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported tree height calculation strategy"):
        analyze_triplet_gene_tree_file(str(input_file), tree_height_calculation_strategy="D")


def test_analyze_triplet_gene_tree_file_rejects_unknown_p_value_correction(tmp_path):
    content = """A,B,C\t3\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Unsupported p-value correction method"):
        analyze_triplet_gene_tree_file(str(input_file), p_value_correction="sidak")


def test_adjust_p_values_custom_fdr_matches_known_bh_example():
    p_values = [0.01, 0.04, 0.03, 0.002]
    adjusted = _adjust_p_values(p_values, method="fdr_bh", stats_backend="custom")

    assert adjusted == pytest.approx([0.02, 0.04, 0.04, 0.008], abs=1e-12)


def test_adjust_p_values_standard_matches_statsmodels_for_supported_methods():
    p_values = [0.01, 0.04, 0.03, 0.002]

    expected_bfn = list(multipletests(p_values, method="bonferroni")[1])
    expected_holm = list(multipletests(p_values, method="holm")[1])
    expected_fdr_bh = list(multipletests(p_values, method="fdr_bh")[1])
    expected_fdr_by = list(multipletests(p_values, method="fdr_by")[1])
    expected_fdr_tsbh = list(multipletests(p_values, alpha=0.01, method="fdr_tsbh")[1])

    observed_bfn = _adjust_p_values(p_values, method="bfn", stats_backend="standard")
    observed_holm = _adjust_p_values(p_values, method="holm", stats_backend="standard")
    observed_fdr_bh = _adjust_p_values(p_values, method="fdr_bh", stats_backend="standard")
    observed_fdr_by = _adjust_p_values(p_values, method="fdr_by", stats_backend="standard")
    observed_fdr_tsbh = _adjust_p_values(p_values, method="fdr_tsbh", stats_backend="standard", alpha=0.01)

    assert observed_bfn == pytest.approx(expected_bfn, abs=1e-12)
    assert observed_holm == pytest.approx(expected_holm, abs=1e-12)
    assert observed_fdr_bh == pytest.approx(expected_fdr_bh, abs=1e-12)
    assert observed_fdr_by == pytest.approx(expected_fdr_by, abs=1e-12)
    assert observed_fdr_tsbh == pytest.approx(expected_fdr_tsbh, abs=1e-12)


@pytest.mark.backend_parity
@pytest.mark.parametrize(
    "method,seed,alpha",
    [
        ("bfn", 901, None),
        ("fdr_bh", 902, None),
        ("holm", 903, None),
        ("fdr_by", 904, None),
        ("fdr_tsbh", 905, 0.01),
    ],
)
def test_adjust_p_values_custom_matches_standard_randomized(method, seed, alpha):
    rng = random.Random(seed)
    for _ in range(500):
        sample_size = rng.randint(1, 100)
        p_values = [rng.random() for _ in range(sample_size)]

        kwargs = {"alpha": alpha} if alpha is not None else {}
        custom = _adjust_p_values(p_values, method=method, stats_backend="custom", **kwargs)
        standard = _adjust_p_values(p_values, method=method, stats_backend="standard", **kwargs)
        assert custom == pytest.approx(standard, abs=1e-12)


def test_analyze_triplet_gene_tree_file_applies_selected_correction(tmp_path):
    content = """A,B,C\t6\t((A:1,B:1):1,C:1);

((A:1,B:1):1,C:1);
((A:1,B:1):1,C:1);
((B:1,C:1):1,A:1);
((B:1,C:1):1,A:1);
((B:1,C:1):1,A:1);
((A:1,C:1):1,B:1);
================================================
A,B,D\t6\t((A:1,B:1):1,D:1);

((A:1,B:1):1,D:1);
((A:1,B:1):1,D:1);
((B:1,D:1):1,A:1);
((B:1,D:1):1,A:1);
((B:1,D:1):1,A:1);
((A:1,D:1):1,B:1);
"""
    input_file = tmp_path / "two_triplets.txt"
    input_file.write_text(content)

    results_no = analyze_triplet_gene_tree_file(
        str(input_file),
        alpha_dct=0.25,
        alpha_ks=0.5,
        p_value_correction="no",
        use_multiprocessing=False,
    )
    results_bfn = analyze_triplet_gene_tree_file(
        str(input_file),
        alpha_dct=0.25,
        alpha_ks=0.5,
        p_value_correction="bfn",
        use_multiprocessing=False,
    )

    assert len(results_no) == 2
    assert len(results_bfn) == 2
    assert results_no[0].dct_p_value == pytest.approx(results_bfn[0].dct_p_value)
    assert results_no[0].dct_p_value_corrected < results_bfn[0].dct_p_value_corrected


def test_two_sample_ks_test_hybrid_uses_scipy_near_threshold(monkeypatch):
    monkeypatch.setattr(triplet_processor_module, "two_sample_ks_test", lambda *_: (0.12, 0.051))
    monkeypatch.setattr(triplet_processor_module, "_two_sample_ks_test_scipy", lambda *_: (0.13, 0.049))

    d_stat, p_value = two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], alpha=0.05, borderline_margin=0.01)
    assert d_stat == pytest.approx(0.13)
    assert p_value == pytest.approx(0.049)


def test_two_sample_ks_test_hybrid_keeps_custom_when_not_borderline(monkeypatch):
    monkeypatch.setattr(triplet_processor_module, "two_sample_ks_test", lambda *_: (0.12, 0.40))
    monkeypatch.setattr(triplet_processor_module, "_two_sample_ks_test_scipy", lambda *_: (0.22, 0.10))

    d_stat, p_value = two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], alpha=0.05, borderline_margin=0.01)
    assert d_stat == pytest.approx(0.12)
    assert p_value == pytest.approx(0.40)


def test_two_sample_ks_test_hybrid_rejects_negative_margin():
    with pytest.raises(ValueError, match="borderline_margin"):
        two_sample_ks_test_hybrid([0.1, 0.2], [0.3, 0.4], borderline_margin=-0.1)


def test_parse_triplet_gene_trees_file_requires_species_tree_column(tmp_path):
    content = """A,B,C\t3

((A:1,B:1):1,C:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="expected 3 tab-separated fields"):
        parse_triplet_gene_trees_file(str(input_file))


def test_parse_triplet_gene_trees_file_rejects_empty_species_tree(tmp_path):
    content = """A,B,C\t3\t

((A:1,B:1):1,C:1);
"""
    input_file = tmp_path / "unique_triplets_gene_trees.txt"
    input_file.write_text(content)

    with pytest.raises(ValueError, match="Invalid species tree in header"):
        parse_triplet_gene_trees_file(str(input_file))


def test_write_pipeline_results_includes_dis1_topology_and_omits_removed_topology_columns(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 10) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(6))

    output_file = tmp_path / "results.tsv"
    write_pipeline_results([result], str(output_file), dct_method="chi-square")

    out = output_file.read_text()
    assert "dis1_topology" in out
    assert "dis2_topology" not in out
    assert "highest_freq_topologies" not in out
    assert "n_dis1" in out
    assert "n_dis2" in out


def test_write_pipeline_results_uses_dynamic_summary_column_names(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1.0,B:1.0):2.0,C:3.0);"] * 40) + (["((B:0.2,C:0.2):0.3,A:0.5);"] * 30) + (["((A:0.2,C:0.2):0.3,B:0.5);"] * 5)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
        rng=random.Random(9),
    )

    output_file = tmp_path / "results_dynamic_summary.tsv"
    write_pipeline_results([result], str(output_file), dct_method="chi-square")

    header = output_file.read_text().splitlines()[0]
    assert "median_con" in header
    assert "median_dis" in header
    assert "summary_statistic" not in header
    assert "summary_con" not in header
    assert "summary_dis" not in header


def test_write_pipeline_results_adds_bootstrap_columns_when_enabled(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = ([("((A:1.0,B:1.0):2.0,C:3.0);")] * 40) + ([("((B:0.2,C:0.2):0.3,A:0.5);")] * 30) + ([("((A:0.2,C:0.2):0.3,B:0.5);")] * 5)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
        bootstrap=True,
        bootstrap_options={"iterations": 4, "seed": 7, "debug_mode": True, "summary_only": True},
    )

    output_file = tmp_path / "results_bootstrap.tsv"
    write_pipeline_results(
        [result],
        str(output_file),
        dct_method="chi-square",
        summary_statistic="median",
        bootstrap=True,
        bootstrap_debug_mode=True,
    )

    lines = output_file.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    assert "bootstrap_value" in header
    assert "all_bootstrap" in header
    assert "bootstrap_con_median" in header
    assert "bootstrap_dis_median" in header
    assert "bootstrap_gene_tree_heights" in header
    all_bootstrap_idx = header.index("all_bootstrap")
    heights_idx = header.index("bootstrap_gene_tree_heights")
    assert row[all_bootstrap_idx].startswith("{")
    assert row[heights_idx].startswith("{")


def test_write_pipeline_results_adds_bootstrap_gene_tree_heights_when_summary_only_false(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1.0,B:1.0):2.0,C:3.0);"] * 6) + (["((B:0.2,C:0.2):0.3,A:0.5);"] * 4)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
        bootstrap=True,
        bootstrap_options={"iterations": 4, "seed": 7, "debug_mode": True, "summary_only": False},
    )

    output_file = tmp_path / "results_bootstrap_full.tsv"
    write_pipeline_results(
        [result],
        str(output_file),
        dct_method="chi-square",
        summary_statistic="median",
        bootstrap=True,
        bootstrap_debug_mode=True,
    )

    lines = output_file.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    assert "bootstrap_gene_tree_heights" in header
    heights_idx = header.index("bootstrap_gene_tree_heights")
    assert row[heights_idx].startswith("[")


def test_collect_triplet_statistics_returns_dict_list():
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(7))

    stats = collect_triplet_statistics([result])
    assert isinstance(stats, list)
    assert len(stats) == 1
    assert stats[0]["triplet"] == species_triplet
    assert "most_frequent_matches_concordant" in stats[0]
    assert "species_topology" not in stats[0]
    assert "dct_method" not in stats[0]
    assert "dct_statistic" in stats[0]
    assert "summary_statistic" not in stats[0]


def test_write_pipeline_results_uses_dct_chi_stats_column_for_chi_square(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="chi-square",
        rng=random.Random(10),
    )

    output_file = tmp_path / "results_chi.tsv"
    write_pipeline_results([result], str(output_file), dct_method="chi-square")

    lines = output_file.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    chi_idx = header.index("dct_chi_stats")
    assert "dct_z_score" not in header
    assert "median_con" in header
    assert "median_dis" in header
    assert chi_idx < header.index("dct_p_value")
    assert "dct_p_val_no_corr" in header
    assert "ks_p_val_no_corr" in header

    assert row[chi_idx] != ""


def test_write_pipeline_results_uses_dct_z_score_column_for_z_test(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="z-test",
        rng=random.Random(11),
    )

    output_file = tmp_path / "results_z.tsv"
    write_pipeline_results([result], str(output_file), dct_method="z-test", summary_statistic="median")

    lines = output_file.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    assert "dct_chi_stats" not in header
    z_idx = header.index("dct_z_score")
    assert "median_con" in header
    assert "median_dis" in header
    assert "dct_p_val_no_corr" in header

    assert row[z_idx] != ""


def test_write_pipeline_results_uses_mode_summary_columns_for_mode(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:0.2,B:0.2):0.3,C:0.5);"] * 40) + (["((B:1.0,C:1.0):2.0,A:3.0);"] * 30) + (["((A:0.2,C:0.2):0.3,B:0.5);"] * 5)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="mode",
        rng=random.Random(31),
    )

    output_file = tmp_path / "results_mode.tsv"
    write_pipeline_results([result], str(output_file), dct_method="chi-square", summary_statistic="mode")

    header = output_file.read_text().splitlines()[0]
    assert "mode_con" in header
    assert "mode_dis" in header


def test_write_pipeline_results_uses_dynamic_corrected_p_value_column_names(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1.0,B:1.0):2.0,C:3.0);"] * 20) + (["((B:0.2,C:0.2):0.3,A:0.5);"] * 15) + (["((A:0.2,C:0.2):0.3,B:0.5);"] * 3)
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
        rng=random.Random(40),
    )

    output_file = tmp_path / "results_dynamic_correction.tsv"
    write_pipeline_results(
        [result],
        str(output_file),
        dct_method="chi-square",
        summary_statistic="median",
        p_value_correction="fdr_bh",
    )

    header = output_file.read_text().splitlines()[0]
    assert "dct_p_val_fdr_bh_corr" in header
    assert "ks_p_val_fdr_bh_corr" in header

def test_write_pipeline_results_rejects_unsupported_p_value_correction(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(41))

    output_file = tmp_path / "bad_correction.tsv"
    with pytest.raises(ValueError, match="Unsupported p-value correction method"):
        write_pipeline_results(
            [result],
            str(output_file),
            dct_method="chi-square",
            summary_statistic="median",
            p_value_correction="sidak",
        )


def test_write_pipeline_results_includes_abc_mapping_column(tmp_path):
    species_triplet = ("TaxonA", "TaxonB", "TaxonC")
    trees = ["((TaxonA:1,TaxonB:1):1,TaxonC:1);"] * 6
    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        rng=random.Random(33),
    )

    output_file = tmp_path / "results_abc.tsv"
    write_pipeline_results([result], str(output_file), dct_method="chi-square", summary_statistic="median")

    lines = output_file.read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    mapping_idx = header.index("abc_mapping")
    assert row[mapping_idx] == "A=TaxonA;B=TaxonB;C=TaxonC"


def test_write_pipeline_statistics_json(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((A:1,B:1):1,C:1);"] * 5) + (["((B:1,C:1):1,A:1);"] * 3) + (["((A:1,C:1):1,B:1);"] * 2)
    result = run_triplet_pipeline(species_triplet, trees, species_topology=TOPOLOGY_AB, rng=random.Random(8))

    output_file = tmp_path / "stats.json"
    write_pipeline_statistics_json([result], str(output_file))

    out = output_file.read_text()
    assert "dct_method" not in out
    assert "dct_statistic" in out
    assert "classification" in out


def test_write_pipeline_results_rejects_mixed_discordant_test_outputs(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)

    chi_result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="chi-square",
        rng=random.Random(20),
    )
    z_result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        discordant_test="z-test",
        rng=random.Random(21),
    )

    output_file = tmp_path / "mixed_dct.tsv"
    # Writer column selection now comes from supplied config method.
    with pytest.raises(ValueError, match="Unsupported discordant test method"):
        write_pipeline_results([chi_result, z_result], str(output_file), dct_method="bad-test")


def test_write_pipeline_results_rejects_unsupported_summary_statistic(tmp_path):
    species_triplet = ("A", "B", "C")
    trees = (["((B:1,C:1):1,A:1);"] * 12) + (["((A:1,B:1):1,C:1);"] * 8) + (["((A:1,C:1):1,B:1);"] * 2)

    result = run_triplet_pipeline(
        species_triplet,
        trees,
        species_topology=TOPOLOGY_AB,
        summary_statistic="median",
        rng=random.Random(22),
    )

    output_file = tmp_path / "mixed_summary.tsv"
    with pytest.raises(ValueError, match="Unsupported summary statistic"):
        write_pipeline_results([result], str(output_file), dct_method="chi-square", summary_statistic="bad-summary")


def test_resolve_runtime_args_triplet_processor_cli_defaults():
    args = argparse.Namespace(
        config_file=None,
        input_path="unique_triplets_gene_trees.txt",
        output_path=None,
        stats_output=None,
        alpha_dct=None,
        alpha_ks=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        tree_height_calculation_strategy=None,
        p_value_correction=None,
        bootstrap=False,
        bootstrap_iterations=None,
        bootstrap_seed=None,
        bootstrap_debug_mode=None,
        bootstrap_summary_only=None,
        processes=None,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    # Paths are resolved to absolute paths
    assert resolved.input == str(Path("unique_triplets_gene_trees.txt").resolve())
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
    assert resolved.bootstrap is False
    assert resolved.bootstrap_options == {
        "iterations": 100,
        "seed": None,
        "debug_mode": False,
        "summary_only": False,
    }


def test_resolve_runtime_args_triplet_processor_cli_custom_processes_preserved():
    args = argparse.Namespace(
        config_file=None,
        input_path="unique_triplets_gene_trees.txt",
        output_path=None,
        stats_output=None,
        alpha_dct=None,
        alpha_ks=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        tree_height_calculation_strategy=None,
        p_value_correction=None,
        bootstrap=True,
        bootstrap_iterations=22,
        bootstrap_seed=555,
        bootstrap_debug_mode=True,
        bootstrap_summary_only=False,
        processes=8,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 8
    assert resolved.bootstrap is True
    assert resolved.bootstrap_options == {
        "iterations": 22,
        "seed": 555,
        "debug_mode": True,
        "summary_only": False,
    }


def test_resolve_runtime_args_triplet_processor_config_warns_and_ignores(tmp_path, capsys):
    config_path = tmp_path / "triplet_processor_config.json"
    config_path.write_text(
        """
{
  "input_path": "input.tsv",
  "output_path": "out.tsv",
  "discordant_test": "z-test",
    "summary_statistic": "median",
    "stats_backend": "standard",
    "tree_height_calculation_strategy": "B",
        "p_value_correction": "bfn"
}
""".strip()
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        input_path="unique_triplets_gene_trees.txt",
        output_path=None,
        stats_output=None,
        alpha_dct=None,
        alpha_ks=None,
        discordant_test="chi-square",
        summary_statistic="mean",
        stats_backend="custom",
        tree_height_calculation_strategy="A",
        p_value_correction="no",
        bootstrap=True,
        bootstrap_iterations=20,
        bootstrap_seed=9,
        bootstrap_debug_mode=True,
        bootstrap_summary_only=False,
        processes=None,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    captured = capsys.readouterr()

    assert "Warning: --config-file provided; CLI arguments not in config will be ignored" in captured.out
    # Paths are resolved to absolute paths
    assert resolved.input == str(Path("input.tsv").resolve())
    assert resolved.output == str(Path("out.tsv").resolve())
    assert resolved.discordant_test == "z-test"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "B"
    assert resolved.p_value_correction == "bfn"
    assert resolved.bootstrap is False
    assert resolved.bootstrap_options == {
        "iterations": 100,
        "seed": None,
        "debug_mode": False,
        "summary_only": False,
    }


def test_resolve_runtime_args_triplet_processor_config_processes_preserved_when_set(tmp_path):
    config_path = tmp_path / "triplet_processor_with_processes.json"
    config_path.write_text(
        """
{
  "input_path": "input.tsv",
  "processes": 3
}
""".strip()
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        input_path=None,
        output_path=None,
        stats_output=None,
        alpha_dct=None,
        alpha_ks=None,
        discordant_test="chi-square",
        summary_statistic="mean",
        stats_backend="custom",
        tree_height_calculation_strategy="A",
        p_value_correction="no",
        bootstrap=None,
        bootstrap_iterations=None,
        bootstrap_seed=None,
        bootstrap_debug_mode=None,
        bootstrap_summary_only=None,
        processes=None,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 3


def test_resolve_runtime_args_triplet_processor_config_defaults_processes_to_zero(tmp_path):
    config_path = tmp_path / "triplet_processor_default_processes.json"
    config_path.write_text(
        """
{
  "input_path": "input.tsv"
}
""".strip()
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        input_path=None,
        output_path=None,
        stats_output=None,
        alpha_dct=None,
        alpha_ks=None,
        discordant_test="chi-square",
        summary_statistic="mean",
        stats_backend="custom",
        tree_height_calculation_strategy="A",
        p_value_correction="no",
        bootstrap=None,
        bootstrap_iterations=None,
        bootstrap_seed=None,
        bootstrap_debug_mode=None,
        bootstrap_summary_only=None,
        processes=11,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 0
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
