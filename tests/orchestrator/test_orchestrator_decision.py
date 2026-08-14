"""Cover the decision logic, p-value correction, and degenerate inference inputs.

The shared end-to-end fixture never produces a significant discordant count
test, so these tests drive the classifier with crafted observation sets chosen
to land on each branch of the decision tree. Expected values are derived from
the definitions (see ``tests/TEST_IO.md``); the correction tests compare against
``statsmodels.multipletests`` called directly.
"""

import pytest
from statsmodels.stats.multitest import multipletests

from ghostparser.orchestrator import inference as pinf

_TRIPLET = ("A", "B", "C")
_SPECIES_SUBTREE = "((A:1.0,B:1.0):1.0,C:2.0);"

_CON = "((A,B),C)"
_DIS1 = "((B,C),A)"
_DIS2 = "((A,C),B)"


def _observations(con_heights, dis1_heights, dis2_heights):
    """Build an observation list with explicit per-topology heights.

    Args:
        con_heights: Heights assigned to the concordant topology.
        dis1_heights: Heights assigned to the ``((B,C),A)`` topology.
        dis2_heights: Heights assigned to the ``((A,C),B)`` topology.

    Returns:
        A list of ``(topology, height, metrics)`` observation tuples.
    """
    observations = []
    for topology, heights in (
        (_CON, con_heights),
        (_DIS1, dis1_heights),
        (_DIS2, dis2_heights),
    ):
        observations.extend((topology, height, None) for height in heights)
    return observations


def _analyze(observations, **kwargs):
    """Run the per-triplet analysis with bootstrap disabled for determinism.

    Args:
        observations: The observation list.
        **kwargs: Forwarded to ``analyze_triplet_from_observations``.

    Returns:
        The resulting ``TripletPipelineResult``.
    """
    return pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": 0},
        **kwargs,
    )


def test_classify_no_introgression_when_dct_not_significant():
    """An even discordant split leaves the DCT non-significant -> no introgression."""
    # 10 vs 10 -> chi-square statistic 0, p = 1.0 > alpha, so the first gate stops.
    result = _analyze(_observations([0.1] * 20, [0.9] * 10, [0.9] * 10))
    assert result.dct_statistic == pytest.approx(0.0)
    assert result.dct_p_value == pytest.approx(1.0)
    assert result.dct_significant is False
    assert result.classification == "no_introgression"


def test_classify_inflow_when_tree_height_test_not_significant():
    """A significant DCT with an inseparable height distribution -> inflow."""
    # 30 vs 2 -> chi-square (30-16)^2/16 * 2 = 24.5, p ~ 7.4e-07 < 0.05.
    # Identical con/dis1 heights make the KS statistic 0 (p = 1.0), so the
    # tree-height test is not significant.
    result = _analyze(_observations([0.5] * 10, [0.5] * 30, [0.5] * 2))
    assert result.dct_significant is True
    assert result.ks_statistic == pytest.approx(0.0)
    assert result.ks_significant is False
    assert result.classification == "inflow_introgression"


# Spread-out, fully separated samples. Every concordant height sits above every
# discordant1 height, so the KS statistic is 1.0; the within-group spread gives
# the studentized statistic a finite standard error to divide by, and the
# combined sample admits C(40, 10) = 8.5e8 assignments, well past the support
# floor. Both permutation guards therefore stay clear and the test resamples.
_HIGH = [0.85 + 0.01 * i for i in range(10)]
_LOW = [0.05 + 0.01 * i for i in range(30)]


def test_classify_outflow_when_concordant_heights_exceed_discordant():
    """Significant DCT + significant KS with con > dis1 -> outflow."""
    result = _analyze(_observations(_HIGH, _LOW, [0.1] * 2))
    assert result.dct_significant is True
    assert result.ks_statistic == pytest.approx(1.0)
    assert result.ks_significant is True
    assert result.mean_con > result.mean_dis
    assert result.perm_note is None
    assert result.perm_decision == "greater"
    assert result.classification == "outflow_introgression"


def test_classify_ghost_when_discordant_heights_exceed_concordant():
    """Significant DCT + significant KS with con < dis1 -> ghost."""
    result = _analyze(_observations(_LOW[:10], _HIGH * 3, [0.9] * 2))
    assert result.dct_significant is True
    assert result.ks_significant is True
    assert result.mean_con < result.mean_dis
    assert result.perm_decision == "less"
    assert result.classification == "ghost_introgression"


def test_classify_ambiguous_when_direction_is_undetectable():
    """Significant DCT + significant KS with no mean separation -> ambiguous.

    The two samples share a mean but differ in spread, so the KS test separates
    the distributions while the permutation test finds no directional evidence.
    """
    con_heights = [0.5 + 0.30 * (1 if i % 2 else -1) for i in range(30)]
    dis1_heights = [0.5 + 0.02 * (1 if i % 2 else -1) for i in range(30)]
    result = _analyze(_observations(con_heights, dis1_heights, [0.5] * 2))
    assert result.dct_significant is True
    assert result.ks_significant is True
    assert result.perm_decision == "ambiguous"
    assert result.classification == "ambiguous"


def test_permutation_guard_reports_insufficient_support():
    """A sample too small to resolve alpha is guarded, not silently decided."""
    # C(6, 2) = 15 distinct assignments, far below the 2500-resample minimum.
    result = _analyze(_observations([0.9, 0.8, 0.7, 0.6], [0.1, 0.2], [0.1] * 2))
    assert result.perm_note == "insufficient_permutation_support"
    assert result.perm_n_resamples == 0
    assert result.perm_decision == "ambiguous"


def test_permutation_guard_reports_degenerate_scale():
    """Internally constant groups with different means give no scale to studentize."""
    result = _analyze(_observations([0.9] * 10, [0.1] * 30, [0.1] * 2))
    assert result.perm_note == "degenerate_observed_scale"
    assert result.perm_decision == "ambiguous"
    assert result.classification == "ambiguous"


def test_median_fallback_decides_direction_when_permutation_disabled():
    """With the permutation test off, the median sign comparison drives direction."""
    result = _analyze(
        _observations([0.9] * 10, [0.1] * 30, [0.1] * 2), permutation_test=False
    )
    assert result.perm_decision == "greater"
    assert result.perm_statistic is None
    assert result.classification == "outflow_introgression"


def _subtrees(con_heights, bc_heights, ac_heights):
    """Build gene subtrees placing given heights on each topology.

    Args:
        con_heights: Heights for the concordant ``((A,B),C)`` topology.
        bc_heights: Heights for the ``((B,C),A)`` topology.
        ac_heights: Heights for the ``((A,C),B)`` topology.

    Returns:
        A list of Newick strings.
    """
    shapes = {
        _CON: "((A:{h:.4f},B:{h:.4f}):0.10,C:{o:.4f});",
        _DIS1: "((B:{h:.4f},C:{h:.4f}):0.10,A:{o:.4f});",
        _DIS2: "((A:{h:.4f},C:{h:.4f}):0.10,B:{o:.4f});",
    }
    subtrees = []
    for topology, heights in (
        (_CON, con_heights),
        (_DIS1, bc_heights),
        (_DIS2, ac_heights),
    ):
        subtrees.extend(shapes[topology].format(h=h, o=h + 0.2) for h in heights)
    return subtrees


def test_summary_statistics_discordant1_follows_frequency_not_topology_name():
    """`discordant1_*` describes the more frequent discordant, not always BC|A.

    Here ``AC|B`` is the more frequent discordant (9 against 3), so the
    ``discordant1_*`` columns must describe the AC group — the same gene trees
    the DCT, KS, and permutation tests use and the same group ``dis1_topology``
    names — rather than the BC group.
    """
    bc_heights = [0.10, 0.12, 0.14]
    ac_heights = [0.50 + 0.01 * i for i in range(9)]
    result = pinf.analyze_triplet(
        _TRIPLET,
        _subtrees([0.30] * 10, bc_heights, ac_heights),
        species_subtree=_SPECIES_SUBTREE,
        collect_summary_statistics=True,
        bootstrap_options={"iterations": 0},
        triplet_seed=1,
    )

    assert result.dis1_topology == "AC"
    assert (result.n_dis1, result.n_dis2) == (9, 3)

    statistics = result.topology_metric_statistics
    # Each subtree is ((X:h,Y:h):0.10, Z:h+0.2), so the sisters sit at h + 0.10
    # and the outlier at h + 0.20: avg tree height = (3h + 0.4) / 3 = h + 0.4/3.
    # Group means therefore shift the input means by exactly 0.4/3.
    expected_dis1 = sum(ac_heights) / len(ac_heights) + 0.4 / 3.0
    expected_dis2 = sum(bc_heights) / len(bc_heights) + 0.4 / 3.0
    assert statistics["discordant1_avg_tree_height_mean"] == pytest.approx(
        expected_dis1
    )
    assert statistics["discordant2_avg_tree_height_mean"] == pytest.approx(
        expected_dis2
    )
    # And it agrees with the group the decision logic compared.
    assert result.mean_dis == pytest.approx(expected_dis1)


def test_summary_statistics_discordant_roles_swap_with_the_counts():
    """Swapping which discordant is more frequent swaps the summary columns."""
    bc_heights = [0.50 + 0.01 * i for i in range(9)]
    ac_heights = [0.10, 0.12, 0.14]
    result = pinf.analyze_triplet(
        _TRIPLET,
        _subtrees([0.30] * 10, bc_heights, ac_heights),
        species_subtree=_SPECIES_SUBTREE,
        collect_summary_statistics=True,
        bootstrap_options={"iterations": 0},
        triplet_seed=1,
    )

    assert result.dis1_topology == "BC"
    statistics = result.topology_metric_statistics
    assert statistics["discordant1_avg_tree_height_mean"] == pytest.approx(
        sum(bc_heights) / len(bc_heights) + 0.4 / 3.0
    )
    assert statistics["discordant2_avg_tree_height_mean"] == pytest.approx(
        sum(ac_heights) / len(ac_heights) + 0.4 / 3.0
    )


@pytest.mark.parametrize(
    "dct_significant,ks_significant,direction,expected",
    [
        (False, True, "greater", "no_introgression"),
        (False, False, "less", "no_introgression"),
        (True, False, "greater", "inflow_introgression"),
        (True, True, "greater", "outflow_introgression"),
        (True, True, "less", "ghost_introgression"),
        (True, True, "ambiguous", "ambiguous"),
        (True, True, None, "ambiguous"),
    ],
)
def test_classify_introgression_truth_table(
    dct_significant, ks_significant, direction, expected
):
    """Every branch of the decision logic maps to its documented classification."""
    assert (
        pinf._classify_introgression(dct_significant, ks_significant, direction)
        == expected
    )


@pytest.mark.parametrize(
    "method", ["no", "bfn", "holm", "fdr_bh", "fdr_by", "fdr_tsbh"]
)
def test_adjust_p_values_matches_statsmodels(method):
    """Each correction method reproduces statsmodels' multipletests output."""
    p_values = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6]

    adjusted = pinf._adjust_p_values(p_values, method=method, alpha=0.05)

    if method == "no":
        assert adjusted == pytest.approx(p_values)
        return

    mapped = {
        "bfn": "bonferroni",
        "holm": "holm",
        "fdr_bh": "fdr_bh",
        "fdr_by": "fdr_by",
        "fdr_tsbh": "fdr_tsbh",
    }[method]
    _, expected, _, _ = multipletests(p_values, alpha=0.05, method=mapped)
    assert adjusted == pytest.approx(list(expected))


def test_adjust_p_values_bonferroni_by_definition():
    """Bonferroni multiplies by the number of tests and clamps at 1.0."""
    p_values = [0.01, 0.2, 0.5]
    adjusted = pinf._adjust_p_values(p_values, method="bfn", alpha=0.05)
    assert adjusted == pytest.approx([0.03, 0.6, 1.0])


def test_adjust_p_values_rejects_unknown_method():
    """An unsupported correction method raises ValueError."""
    with pytest.raises(ValueError, match="Unsupported p-value correction method"):
        pinf._adjust_p_values([0.1], method="bogus")


@pytest.mark.parametrize("method", ["chi-square", "z-test"])
def test_discordant_count_test_with_no_discordant_observations(method):
    """A zero/zero discordant split short-circuits to a non-significant result."""
    assert pinf.run_discordant_count_test(0, 0, method=method) == (0.0, 1.0)


def test_discordant_count_test_rejects_unknown_method():
    """An unsupported discordant test raises ValueError."""
    with pytest.raises(ValueError, match="Unsupported discordant test method"):
        pinf.run_discordant_count_test(3, 2, method="bogus")


@pytest.mark.parametrize(
    "sample_a,sample_b", [([], [1.0]), ([1.0], []), ([], [])]
)
def test_ks_test_with_an_empty_sample(sample_a, sample_b):
    """An empty sample makes the KS test a non-significant no-op."""
    assert pinf.run_two_sample_ks_test(sample_a, sample_b) == (0.0, 1.0)
