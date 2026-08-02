"""Cover the decision logic, p-value correction, and degenerate inference inputs.

The shared end-to-end fixture never produces a significant discordant count
test, so these tests drive the classifier with crafted observation sets chosen
to land on each branch of the decision tree. Expected values are derived from
the definitions (see ``tests/TEST_IO.md``); the correction tests compare against
``statsmodels.multipletests`` called directly.
"""

import pytest
from statsmodels.stats.multitest import multipletests

from ghostparser.pipeline import inference as pinf

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


def test_classify_outflow_when_concordant_heights_exceed_discordant():
    """Significant DCT + significant KS with con > dis1 -> outflow."""
    # Fully separated samples: every concordant height exceeds every dis1 height,
    # so the KS statistic is 1.0 and the summary comparison is con > dis1.
    result = _analyze(_observations([0.9] * 10, [0.1] * 30, [0.1] * 2))
    assert result.dct_significant is True
    assert result.ks_statistic == pytest.approx(1.0)
    assert result.ks_significant is True
    assert result.summary_con > result.summary_dis
    assert result.classification == "outflow_introgression"


def test_classify_ghost_when_discordant_heights_exceed_concordant():
    """Significant DCT + significant KS with con < dis1 -> ghost."""
    result = _analyze(_observations([0.1] * 10, [0.9] * 30, [0.9] * 2))
    assert result.dct_significant is True
    assert result.ks_significant is True
    assert result.summary_con < result.summary_dis
    assert result.classification == "ghost_introgression"


@pytest.mark.parametrize(
    "dct_significant,ks_significant,summary_con,summary_dis,expected",
    [
        (False, True, 1.0, 2.0, "no_introgression"),
        (False, False, 1.0, 2.0, "no_introgression"),
        (True, False, 1.0, 2.0, "inflow_introgression"),
        (True, True, 2.0, 1.0, "outflow_introgression"),
        (True, True, 1.0, 2.0, "ghost_introgression"),
        (True, True, 1.0, 1.0, "unresolved"),
        (True, True, None, 1.0, "unresolved"),
        (True, True, 1.0, None, "unresolved"),
    ],
)
def test_classify_introgression_truth_table(
    dct_significant, ks_significant, summary_con, summary_dis, expected
):
    """Every branch of the decision logic maps to its documented classification."""
    assert (
        pinf._classify_introgression(
            dct_significant, ks_significant, summary_con, summary_dis
        )
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
