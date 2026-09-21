"""Verify orchestrator.inference derives topology counts, tree heights, the DCT/KS tests, and classification from their definitions.

Expected values are recomputed inside each test from the primitive quantities of
the input Newick strings (root-to-tip distances and the sister-pair internal
branch, hand-derived and tabulated in ``_LEAF_GEOMETRY``) plus direct calls to
the reference statistical libraries (SciPy / statsmodels). Nothing here compares
against another GhostParser module. See ``tests/TEST_IO.md`` for the full
derivation of every literal below.
"""

import pytest
from scipy import stats
from statsmodels.stats.proportion import proportions_ztest

from ghostparser.orchestrator import inference as pinf

_TRIPLET = ("A", "B", "C")
_SPECIES_SUBTREE = "((A:1.0,B:1.0):1.0,C:2.0);"
_GENE_SUBTREES = [
    "((A:0.10,B:0.10):0.10,C:0.30);",
    "((A:0.12,B:0.11):0.09,C:0.32);",
    "((A:0.20,B:0.15):0.10,C:0.40);",
    "((B:0.20,C:0.20):0.10,A:0.50);",
    "((B:0.25,C:0.22):0.10,A:0.55);",
    "((A:0.30,C:0.30):0.10,B:0.70);",
    "((A:0.15,B:0.15):0.12,C:0.35);",
    "((B:0.35,C:0.35):0.10,A:0.65);",
    "((A:0.18,C:0.16):0.10,B:0.60);",
    "((A:0.11,B:0.13):0.10,C:0.31);",
]

# Hand-derived geometry of each gene subtree above, in input order:
#   (topology, root-to-tip A, root-to-tip B, root-to-tip C, sister-pair internal
#    branch length).
# The topology label names the sister pair; the sister pair's root-to-tip
# distance is (its own edge + the internal branch), while the third taxon hangs
# directly off the root. Example for index 0, "((A:0.10,B:0.10):0.10,C:0.30);":
# A = 0.10 + 0.10 = 0.20, B = 0.10 + 0.10 = 0.20, C = 0.30, internal = 0.10.
_LEAF_GEOMETRY = [
    ("((A,B),C)", 0.20, 0.20, 0.30, 0.10),
    ("((A,B),C)", 0.21, 0.20, 0.32, 0.09),
    ("((A,B),C)", 0.30, 0.25, 0.40, 0.10),
    ("((B,C),A)", 0.50, 0.30, 0.30, 0.10),
    ("((B,C),A)", 0.55, 0.35, 0.32, 0.10),
    ("((A,C),B)", 0.40, 0.70, 0.40, 0.10),
    ("((A,B),C)", 0.27, 0.27, 0.35, 0.12),
    ("((B,C),A)", 0.65, 0.45, 0.45, 0.10),
    ("((A,C),B)", 0.28, 0.60, 0.26, 0.10),
    ("((A,B),C)", 0.21, 0.23, 0.31, 0.10),
]

_SISTER_TAXA = {
    "((A,B),C)": ("A", "B"),
    "((B,C),A)": ("B", "C"),
    "((A,C),B)": ("A", "C"),
}
_CONCORDANT = "((A,B),C)"

_SEED = 20240724
_ITERATIONS = 40


def _expected_height(entry, strategy):
    """Compute H(T) for one gene subtree straight from the strategy definition.

    Args:
        entry: A ``_LEAF_GEOMETRY`` row.
        strategy: One of ``AVG``/``A``/``B``/``C``/``SIS``/``INT``.

    Returns:
        The expected tree-height value as a float.
    """
    topology, dist_a, dist_b, dist_c, internal = entry
    by_taxon = {"A": dist_a, "B": dist_b, "C": dist_c}

    if strategy == "AVG":
        return (dist_a + dist_b + dist_c) / 3.0
    if strategy in {"A", "B", "C"}:
        return by_taxon[strategy]
    left, right = _SISTER_TAXA[topology]
    if strategy == "SIS":
        return by_taxon[left] + by_taxon[right] - 2.0 * internal
    return internal  # INT


def _expected_dct(n_dis1, n_dis2, discordant_test):
    """Compute the expected DCT statistic/p-value with the reference libraries.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.
        discordant_test: ``chi-square`` or ``z-test``.

    Returns:
        A tuple ``(statistic, p_value)``.
    """
    if discordant_test == "chi-square":
        result = stats.chisquare([n_dis1, n_dis2])
        return float(result.statistic), float(result.pvalue)
    total = n_dis1 + n_dis2
    statistic, p_value = proportions_ztest(
        count=[n_dis1, n_dis2], nobs=[total, total], alternative="two-sided"
    )
    return float(statistic), float(p_value)


def _expected_result(strategy, discordant_test, diagnostic=True, alpha=0.05):
    """Derive every asserted inference field for the shared fixture.

    Groups the hand-derived heights by topology, ranks the two discordant
    topologies, runs the DCT and the con-vs-dis1 KS test through SciPy /
    statsmodels, and applies the GhostParser decision logic. A non-diagnostic
    run stops measuring once the cascade is settled, so the expectation carries
    ``None`` for whatever it never computes.

    Args:
        strategy: The tree-height strategy.
        discordant_test: ``chi-square``/``z-test``.
        diagnostic: ``True`` measures all three tests; ``False`` skips the
            tests the cascade cannot consult.
        alpha: Significance threshold shared by both tests.

    Returns:
        A dict of the expected field values.
    """
    heights = {topology: [] for topology in _SISTER_TAXA}
    for entry in _LEAF_GEOMETRY:
        heights[entry[0]].append(_expected_height(entry, strategy))

    # Rank the two non-concordant topologies; ties fall to the first listed.
    discordant = [t for t in ("((B,C),A)", "((A,C),B)") if t != _CONCORDANT]
    first, second = discordant[0], discordant[1]
    if len(heights[first]) >= len(heights[second]):
        dis1, dis2 = first, second
    else:
        dis1, dis2 = second, first

    n_con, n_dis1, n_dis2 = (
        len(heights[_CONCORDANT]),
        len(heights[dis1]),
        len(heights[dis2]),
    )

    dct_statistic, dct_p_value = _expected_dct(n_dis1, n_dis2, discordant_test)
    dct_significant = dct_p_value < alpha

    # A failed count gate settles the call on its own. The decision pass here
    # runs under ``no``, an inline correction, so a non-diagnostic run leaves
    # the tree-height test below it unmeasured; a diagnostic run measures and
    # corrects it like any other member of the family.
    ks_statistic = ks_p_value = ks_significant = None
    if dct_significant or diagnostic:
        ks_result = stats.ks_2samp(
            heights[_CONCORDANT], heights[dis1], alternative="two-sided", method="auto"
        )
        ks_statistic, ks_p_value = float(ks_result.statistic), float(ks_result.pvalue)
        ks_significant = ks_p_value < alpha

    # This fixture has 5 concordant and 3 discordant1 trees, so the pooled
    # sample admits only C(8, 3) = 56 distinct group assignments -- far fewer
    # than the 2500-resample minimum. The support guard fires and reports no
    # conclusion without any resampling. A non-diagnostic run does not even get
    # that far unless both earlier gates cleared.
    direction = "inconclusive"
    if not diagnostic and not (dct_significant and ks_significant):
        direction = None

    # GhostParser decision logic: DCT gate, then the tree-height test, then the
    # concordant-vs-discordant1 direction test.
    if not dct_significant:
        classification = "no_introgression"
    elif not ks_significant:
        classification = "inflow_introgression"
    elif direction == "greater":
        classification = "outflow_introgression"
    elif direction == "less":
        classification = "ghost_introgression"
    else:
        classification = "ambiguous"

    return {
        "n_con": n_con,
        "n_dis1": n_dis1,
        "n_dis2": n_dis2,
        "dis1_topology": "BC" if dis1 == "((B,C),A)" else "AC",
        "most_frequent_matches_concordant": n_con >= n_dis1 and n_con >= n_dis2,
        "dct_statistic": dct_statistic,
        "dct_p_value": dct_p_value,
        "dct_significant": dct_significant,
        "ks_statistic": ks_statistic,
        "ks_p_value": ks_p_value,
        "ks_significant": ks_significant,
        "perm_decision": direction,
        "classification": classification,
        "analyzed_trees": len(_LEAF_GEOMETRY),
    }


def _assert_matches_expected(result, expected):
    """Assert a result equals a derived expectation field by field.

    Args:
        result: The ``TripletPipelineResult`` under test.
        expected: The mapping returned by :func:`_expected_result`.
    """
    for field, expected_value in expected.items():
        actual = getattr(result, field)
        if isinstance(expected_value, float):
            assert actual == pytest.approx(expected_value), field
        else:
            assert actual == expected_value, field


def _decided(observations, **kwargs):
    """Measure a triplet and run the decision pass over it as a family of one.

    Args:
        observations: The observation list.
        **kwargs: Forwarded to ``analyze_triplet_from_observations``.

    Returns:
        The decided ``TripletPipelineResult``.
    """
    result = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
        **kwargs,
    )
    return pinf._apply_triplet_result_p_value_correction(
        [result], alpha_dct=0.05, alpha_ks=0.05, method="no"
    )[0]


@pytest.mark.parametrize("diagnostic", [False, True])
@pytest.mark.parametrize("discordant_test", ["chi-square", "z-test"])
def test_inference_matches_derived_expectation(discordant_test, diagnostic):
    """The tests and the decision reproduce values derived from the definitions.

    Runs the shared 10-gene-subtree fixture through the reference serializer,
    the per-triplet measurement and the decision pass, and checks every
    asserted field against SciPy / statsmodels and the cascade written out by
    hand, with ``diagnostic`` both off and on: the derived statistics are the
    same wherever a run measures them, and a non-diagnostic run leaves the rest
    unmeasured.
    The heights themselves are pinned per strategy by
    ``test_observation_heights_match_derived_geometry``, so one strategy is
    enough here.
    """
    observations = pinf._serialize_triplet_gene_trees(
        _TRIPLET, _GENE_SUBTREES, tree_height_calculation_strategy="AVG"
    )
    result = _decided(
        observations, discordant_test=discordant_test, diagnostic=diagnostic
    )

    _assert_matches_expected(
        result, _expected_result("AVG", discordant_test, diagnostic)
    )
    assert tuple(result.triplet) == _TRIPLET
    assert result.species_tree == "((A,B),C);"
    assert 0.0 <= result.bootstrap_value <= 1.0
    assert sum(result.all_bootstrap.values()) == pytest.approx(1.0)


def test_observation_heights_match_derived_geometry():
    """Serialized observations carry the hand-derived topology and H(T) under every strategy."""
    for strategy in ("AVG", "A", "B", "C", "SIS", "INT"):
        observations = pinf._serialize_triplet_gene_trees(
            _TRIPLET, _GENE_SUBTREES, tree_height_calculation_strategy=strategy
        )
        assert len(observations) == len(_LEAF_GEOMETRY)
        for observation, entry in zip(observations, _LEAF_GEOMETRY):
            assert observation[0] == entry[0], strategy
            assert observation[1] == pytest.approx(_expected_height(entry, strategy)), strategy


def test_empty_observations_decide_no_introgression():
    """Zero observations measure as empty groups and decide ``no_introgression``."""
    result = _decided([])
    assert result.analyzed_trees == 0
    assert result.n_con == 0
    assert result.classification == "no_introgression"
