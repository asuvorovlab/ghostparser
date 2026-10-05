"""Cover the decision logic, p-value correction, and degenerate inference inputs.

The shared end-to-end fixture never produces a significant discordant count
test, so these tests drive the classifier with crafted observation sets chosen
to land on each branch of the decision tree. Expected values are derived from
the definitions (see ``tests/TEST_IO.md``); the correction tests compare against
``statsmodels.multipletests`` called directly.
"""

from dataclasses import replace

import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

from ghostparser.orchestrator import inference as pinf
from tests.orchestrator.tree_references import serialize_triplet_gene_trees
from ghostparser.orchestrator.config import P_VALUE_CORRECTION_CHOICES
from ghostparser.orchestrator.correction import is_inline_correction

_TRIPLET = ("A", "B", "C")
_SPECIES_SUBTREE = "((A,B),C);"

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


def _measure(observations, method, family_size=1, iterations=40, seed=3, **kwargs):
    """Measure one triplet without deciding it.

    Args:
        observations: The observation list.
        method: The p-value correction method.
        family_size: Triplet count the correction should work against.
        iterations: Bootstrap iterations to run.
        seed: Base seed for deterministic resampling.
        **kwargs: Forwarded to ``analyze_triplet_from_observations``.

    Returns:
        The undecided ``TripletPipelineResult``.
    """
    return pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": iterations},
        p_value_correction=method,
        family_size=family_size,
        triplet_seed=seed,
        **kwargs,
    )


def _decide(results, method):
    """Run the run-wide decision pass over measured results."""
    return pinf._apply_triplet_result_p_value_correction(
        results, alpha_dct=0.05, alpha_ks=0.05, method=method
    )


def _corrected(observations, method, family_size=1, iterations=40, seed=3, **kwargs):
    """Measure one triplet and run the run-wide decision pass over it.

    Both correction tiers finish in this pass: the inline methods have already
    tallied their bootstrap votes, the rank-based ones are resolved here. When
    ``family_size`` exceeds 1 the triplet is padded out to that many results,
    each padding entry carrying a p-value of 1.0 on both tests, so the point
    estimate is corrected against the same family the bootstrap used.

    Args:
        observations: The observation list.
        method: The p-value correction method.
        family_size: Triplet count the correction should work against.
        iterations: Bootstrap iterations to run.
        seed: Base seed for deterministic resampling.
        **kwargs: Forwarded to ``analyze_triplet_from_observations``.

    Returns:
        The corrected ``TripletPipelineResult`` for the real triplet.
    """
    result = _measure(
        observations, method, family_size=family_size, iterations=iterations,
        seed=seed, **kwargs,
    )
    # Padding stands in for triplets that never clear a gate. Under a rank-based
    # method they need their own deferred records, not an empty one: every
    # triplet in a real run carries a record, and the per-iteration family is
    # sized from the records present, so dropping them would silently
    # under-correct the bootstrap. An inline method has already voted.
    deferred = None
    if result.bootstrap_deferred is not None:
        deferred = pinf.DeferredBootstrapRecord(
            dct_p_values=np.ones(iterations, dtype=np.float64),
            ks_p_values=np.ones(iterations, dtype=np.float64),
            directions=np.zeros(iterations, dtype=np.int8),
        )
    padding = [
        replace(
            result,
            dct_p_value=1.0,
            ks_p_value=1.0,
            all_bootstrap=None,
            bootstrap_deferred=deferred,
        )
        for _ in range(family_size - 1)
    ]
    return _decide([result] + padding, method)[0]


# Spread-out, fully separated samples. Every concordant height sits above every
# discordant1 height, so the KS statistic is 1.0; the within-group spread gives
# the studentized statistic a finite standard error to divide by, and the
# combined sample admits C(40, 10) = 8.5e8 assignments, well past the support
# floor. Both permutation guards therefore stay clear and the test resamples.
_HIGH = [0.85 + 0.01 * i for i in range(10)]
_LOW = [0.05 + 0.01 * i for i in range(30)]

# Samples sharing a mean but differing in spread: the KS test separates the
# distributions while the permutation test finds no direction.
_WIDE = [0.5 + 0.30 * (1 if i % 2 else -1) for i in range(30)]
_NARROW = [0.5 + 0.02 * (1 if i % 2 else -1) for i in range(30)]


@pytest.mark.parametrize(
    "con, dis1, dis2, dct_significant, ks_significant, decisions, gate, expected",
    [
        # 10 vs 10 -> chi-square statistic 0, p = 1.0, so the first gate stops.
        # KS is significant here too, so this row also shows the DCT gate
        # stopping the cascade before a later gate can be consulted; a
        # non-diagnostic run never measures it, so its flag is undefined there.
        ([0.1] * 20, [0.9] * 10, [0.9] * 10, False, True, None, "DCT",
         "no_introgression"),
        # 30 vs 2 -> chi-square 24.5, p ~ 7.4e-07. Identical con/dis1 heights
        # make the KS statistic 0 (p = 1.0), so gate 2 stops.
        ([0.5] * 10, [0.5] * 30, [0.5] * 2, True, False, None, "THT",
         "inflow_introgression"),
        (_HIGH, _LOW, [0.1] * 2, True, True, {"greater"}, "PERM",
         "outflow_introgression"),
        (_LOW[:10], _HIGH * 3, [0.9] * 2, True, True, {"less"}, "PERM",
         "ghost_introgression"),
        (_WIDE, _NARROW, [0.5] * 2, True, True, {"equivalent", "inconclusive"},
         "PERM", "ambiguous"),
    ],
    ids=["no_introgression", "inflow", "outflow", "ghost", "ambiguous"],
)
@pytest.mark.parametrize("diagnostic", [False, True])
def test_decision_cascade_lands_on_each_classification(
    con, dis1, dis2, dct_significant, ks_significant, decisions, gate, expected,
    diagnostic,
):
    """Crafted observation sets drive the cascade onto each of its five outcomes.

    The gate asserts which test settled the call, so a case that reaches its
    classification by the wrong route fails rather than passing by coincidence.
    The ``diagnostic`` setting must not move the classification or the gate: a
    diagnostic run measures and reports every test on every row, while the
    default declines a test whose result the cascade would have ignored: it
    reports ``perm_decision`` exactly on the rows the permutation gate settled,
    and under the inline ``no`` correction leaves the tree-height flag
    undefined below a settled count gate.
    """
    result = _corrected(
        _observations(con, dis1, dis2), "no", iterations=0, diagnostic=diagnostic
    )
    assert result.dct_significant is dct_significant
    assert result.classification == expected
    assert result.decision_gate == gate

    if diagnostic or gate != "DCT":
        assert result.ks_significant is ks_significant
    else:
        assert result.ks_p_value is None
        assert result.ks_significant is None

    if diagnostic or gate == "PERM":
        assert result.perm_decision is not None
        if decisions is not None:
            assert result.perm_decision in decisions
    else:
        assert result.perm_decision is None
        assert result.perm_note == pinf.PERM_NOTE_NOT_CONSULTED


@pytest.mark.parametrize(
    "con, dis1, note",
    [
        # C(6, 2) = 15 distinct assignments, far below the 2500-resample minimum.
        ([0.9, 0.8, 0.7, 0.6], [0.1, 0.2], "insufficient_permutation_support"),
        # Internally constant groups leave no scale to studentize by.
        ([0.9] * 10, [0.1] * 30, "degenerate_observed_scale"),
    ],
)
def test_permutation_guards_surface_on_the_triplet_result(con, dis1, note):
    """A guarded direction test reports its reason instead of a direction.

    Runs diagnostic so the test is reached regardless of what the earlier
    gates decided; a guard is a property of the samples, not of the cascade
    position.
    """
    result = _analyze(_observations(con, dis1, [0.1] * 2), diagnostic=True)
    assert result.perm_note == note
    assert result.perm_n_resamples == 0
    assert result.perm_decision == "inconclusive"


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


@pytest.mark.parametrize(
    "bc_heights, ac_heights, expected_dis1_topology",
    [
        ([0.10, 0.12, 0.14], [0.50 + 0.01 * i for i in range(9)], "AC"),
        ([0.50 + 0.01 * i for i in range(9)], [0.10, 0.12, 0.14], "BC"),
    ],
    ids=["AC_more_frequent", "BC_more_frequent"],
)
def test_summary_statistics_discordant_roles_follow_the_counts(
    bc_heights, ac_heights, expected_dis1_topology
):
    """`discordant1_*` describes the more frequent discordant, not always BC|A.

    The summary columns must name the same gene trees ``dis1_topology`` names and
    the three tests operate on, in either direction of the count.
    """
    observations = serialize_triplet_gene_trees(
        _TRIPLET,
        _subtrees([0.30] * 10, bc_heights, ac_heights),
        collect_summary_statistics=True,
    )
    result = _analyze(observations)

    assert result.dis1_topology == expected_dis1_topology
    assert (result.n_dis1, result.n_dis2) == (9, 3)

    dis1, dis2 = (
        (ac_heights, bc_heights)
        if expected_dis1_topology == "AC"
        else (bc_heights, ac_heights)
    )
    statistics = result.topology_metric_statistics
    # Each subtree is ((X:h,Y:h):0.10, Z:h+0.2), so the sisters sit at h + 0.10
    # and the outlier at h + 0.20: avg tree height = (3h + 0.4) / 3 = h + 0.4/3.
    assert statistics["discordant1_avg_tree_height_mean"] == pytest.approx(
        sum(dis1) / len(dis1) + 0.4 / 3.0
    )
    assert statistics["discordant2_avg_tree_height_mean"] == pytest.approx(
        sum(dis2) / len(dis2) + 0.4 / 3.0
    )


@pytest.mark.parametrize(
    "dct_significant,ks_significant,direction,expected,expected_gate",
    [
        (False, True, "greater", "no_introgression", "DCT"),
        (False, False, "less", "no_introgression", "DCT"),
        (True, False, "greater", "inflow_introgression", "THT"),
        (True, None, "greater", "inflow_introgression", "THT"),
        (True, True, "greater", "outflow_introgression", "PERM"),
        (True, True, "less", "ghost_introgression", "PERM"),
        (True, True, "equivalent", "ambiguous", "PERM"),
        (True, True, "inconclusive", "ambiguous", "PERM"),
        (True, True, None, "ambiguous", "PERM"),
    ],
)
def test_classify_introgression_truth_table(
    dct_significant, ks_significant, direction, expected, expected_gate
):
    """Every branch of the decision logic maps to its classification and its gate.

    ``_classify_introgression`` returns both in one pass, so the pair is
    asserted over the same truth table: the gate names the test that settled the
    call, and only ``Permutation`` means the direction decision was consulted.
    ``_resolve_deferred_bootstrap`` classifies whole arrays at once rather than
    calling the scalar form per iteration, so the array form is held to the
    same table, an unmeasured tree-height flag reads as not significant
    there, as the scalar form documents.
    """
    assert pinf._classify_introgression(
        dct_significant, ks_significant, direction
    ) == (expected, expected_gate)

    code = pinf._classification_codes(
        np.array([dct_significant]),
        np.array([bool(ks_significant)]),
        np.array([pinf._DIRECTION_CODES.get(direction, pinf._DIRECTION_SKIPPED)]),
    )[0]
    assert pinf._BOOTSTRAP_CLASSES[code] == expected


@pytest.mark.parametrize("method", P_VALUE_CORRECTION_CHOICES)
def test_adjust_p_values_matches_statsmodels_and_never_lowers_a_value(method):
    """Each correction reproduces statsmodels and never pulls a p-value down.

    The first half compares against ``multipletests`` called directly (``no``
    is the identity). The second is what licenses skipping a settled gate:
    both short-circuits (the point estimate skipping the direction test and
    a bootstrap iteration skipping it) rest on ``raw > alpha`` implying
    ``adjusted > alpha``, so the property is asserted over the whole choice
    list rather than a fixed set of names, on the hostile family of eight
    strong signals against two weak ones, where a true-null-count estimator
    would drive its multiplier below 1.
    """
    p_values = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6]
    adjusted = pinf._adjust_p_values(p_values, method=method, alpha=0.05)
    if method == "no":
        assert adjusted == pytest.approx(p_values)
    else:
        mapped = {"bfn": "bonferroni"}.get(method, method)
        _, expected, _, _ = multipletests(p_values, alpha=0.05, method=mapped)
        assert adjusted == pytest.approx(list(expected))

    hostile = [0.001] * 8 + [0.4, 0.9]
    adjusted = pinf._adjust_p_values(hostile, method=method, alpha=0.05)
    assert all(a >= p - 1e-12 for a, p in zip(adjusted, hostile))


@pytest.mark.parametrize("method", P_VALUE_CORRECTION_CHOICES)
def test_inline_correction_matches_the_family_pass_or_refuses(method):
    """A method needing only the family size reproduces the full pass; the rest refuse.

    ``no`` and ``bfn`` vote inline during the stream, correcting one p-value
    from the triplet count alone, so the value must equal what the run-wide
    pass gives that p-value inside a family of that size. A rank-based method
    cannot be applied without the rest of its family and says so.
    """
    p_value = 0.004
    if not is_inline_correction(method):
        with pytest.raises(ValueError, match="needs the whole family"):
            pinf._adjust_p_value_inline(p_value, method, 10)
        return

    for family_size in (1, 7, 250):
        family = [p_value] + [0.5] * (family_size - 1)
        expected = pinf._adjust_p_values(family, method=method, alpha=0.05)[0]
        assert pinf._adjust_p_value_inline(p_value, method, family_size) == pytest.approx(
            expected
        ), family_size


@pytest.mark.parametrize(
    "call, match",
    [
        (lambda: pinf._adjust_p_values([0.1], method="bogus"), "Unsupported p-value correction method"),
        (lambda: pinf.run_discordant_count_test(3, 2, method="bogus"), "Unsupported discordant test method"),
    ],
    ids=["correction", "discordant_test"],
)
def test_unknown_methods_are_rejected(call, match):
    """An unsupported correction or count-test method raises ValueError."""
    with pytest.raises(ValueError, match=match):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: pinf.run_discordant_count_test(0, 0, method="chi-square"),
        lambda: pinf.run_discordant_count_test(0, 0, method="z-test"),
        lambda: pinf.run_two_sample_ks_test([], [1.0]),
        lambda: pinf.run_two_sample_ks_test([], []),
    ],
    ids=["dct_chi_square_no_discordants", "dct_z_test_no_discordants", "ks_one_empty", "ks_both_empty"],
)
def test_degenerate_samples_are_non_significant(call):
    """A test with nothing to compare reports a zero statistic and p = 1."""
    assert call() == (0.0, 1.0)


@pytest.mark.parametrize("method", P_VALUE_CORRECTION_CHOICES)
def test_inline_and_deferred_correction_agree_on_a_single_triplet(method):
    """Every method votes the same way on a family of one, whichever tier it uses.

    A family of one leaves each correction as the identity, so all five methods
    must produce the same bootstrap tally. ``no`` and ``bfn`` reach it inline
    during the stream while the rank-based three park their raw p-values in a
    deferred record (one entry per iteration, whatever the family size, since
    every rank depends on the other triplets' values) which the decision pass
    resolves and clears. Agreeing here is what shows the two code paths
    implement one decision rule rather than two.
    """
    observations = _observations([0.9] * 25, [0.2] * 3 + [0.35] * 22, [0.3] * 4)
    measured = _measure(observations, method)
    if is_inline_correction(method):
        assert measured.bootstrap_deferred is None
        assert measured.all_bootstrap is not None
    else:
        assert measured.all_bootstrap is None
        assert len(measured.bootstrap_deferred.dct_p_values) == 40

    result = _decide([measured], method)[0]
    baseline = _corrected(observations, "no")
    assert result.bootstrap_deferred is None
    assert result.all_bootstrap == baseline.all_bootstrap
    assert result.bootstrap_value == baseline.bootstrap_value


def test_bootstrap_votes_answer_to_the_corrected_threshold():
    """Bootstrap iterations are judged against the same corrected alpha as the point estimate.

    The raw DCT p-value here clears alpha, but Bonferroni over a family of 5000
    pushes it well past. The point estimate therefore reports
    ``no_introgression``, and the bootstrap must agree: judged on raw p-values
    the same iterations would vote for an introgression class instead, which is
    the mismatch that made ``bootstrap_value`` unreadable.
    """
    observations = _observations([0.9] * 40, [0.55] * 18, [0.55] * 6)
    raw = _corrected(observations, "no")
    corrected = _corrected(observations, "bfn", family_size=5000)

    assert raw.dct_p_value <= 0.05
    assert corrected.dct_p_value_corrected > 0.05
    assert corrected.classification == "no_introgression"
    assert corrected.all_bootstrap["no_introgression"] == 1.0
    assert corrected.bootstrap_value == 1.0
    # The same iterations vote differently when nothing corrects them, which is
    # what makes the corrected agreement above meaningful rather than vacuous.
    assert raw.all_bootstrap["no_introgression"] < 1.0


@pytest.mark.parametrize(
    "con, dis1, dis2",
    [
        # An even 10/10 discordant split fails the count gate in most
        # resamples, so most of the direction tests are the record's own.
        (_HIGH, _LOW[:10], [0.5] * 10),
        # Every gate clears in most resamples, so most direction tests are the
        # vote's own and the record must carry exactly those.
        (_HIGH, _LOW, [0.1] * 2),
    ],
    ids=["settled_at_count_gate", "reaches_direction_gate"],
)
def test_diagnostic_bootstrap_records_every_test_without_moving_a_vote(con, dis1, dis2):
    """A diagnostic bootstrap measures all three tests per iteration and changes no vote.

    The record must be complete: one entry per iteration for every test,
    the direction test included even where a failed gate meant the vote never
    read it, and it must be the vote's own numbers: replaying the cascade
    over the recorded p-values and decisions rebuilds ``all_bootstrap``
    exactly. The extra direction tests draw from their own stream, so the
    votes and the studentized interval are identical with the record on or
    off. In summary form the decision counts sum to the iteration count.
    """
    iterations = 30
    runs = {}
    for diagnostic in (False, True):
        result = pinf.analyze_triplet_from_observations(
            _TRIPLET,
            _observations(con, dis1, dis2),
            species_subtree=_SPECIES_SUBTREE,
            bootstrap_options={"iterations": iterations, "diagnostic": diagnostic},
            p_value_correction="bfn",
            family_size=1,
            triplet_seed=11,
        )
        runs[diagnostic] = pinf._apply_triplet_result_p_value_correction(
            [result], alpha_dct=0.05, alpha_ks=0.05, method="bfn"
        )[0]
    lean, full = runs[False], runs[True]

    assert full.all_bootstrap == lean.all_bootstrap
    assert full.bootstrap_perm_stat_ci_low == lean.bootstrap_perm_stat_ci_low
    assert full.bootstrap_perm_stat_ci_high == lean.bootstrap_perm_stat_ci_high
    assert lean.bootstrap_perm_decisions is None

    record = {
        "dct_p": full.bootstrap_dct_p_value,
        "ks_p": full.bootstrap_ks_p_value,
        "perm_stat": full.bootstrap_perm_stats,
        "perm_p_greater": full.bootstrap_perm_p_greater,
        "perm_p_less": full.bootstrap_perm_p_less,
        "decision": full.bootstrap_perm_decisions,
    }
    for values in record.values():
        assert len(values) == iterations
    assert set(record["decision"]) <= {"greater", "less", "inconclusive"}
    # Every iteration resampled with spread in both groups, so the direction
    # test ran to a statistic and a p-value pair everywhere.
    assert all(value is not None for value in record["perm_stat"])
    assert all(value is not None for value in record["perm_p_greater"])

    # Under bfn as a family of one, a gate is judged on the raw p-value; the
    # cascade over the recorded values must rebuild the votes exactly.
    tally = {}
    for dct_p, ks_p, decision in zip(record["dct_p"], record["ks_p"], record["decision"]):
        classification, _ = pinf._classify_introgression(
            dct_p <= 0.05, ks_p <= 0.05, decision
        )
        tally[classification] = tally.get(classification, 0) + 1
    rebuilt = {label: tally.get(label, 0) / iterations for label in full.all_bootstrap}
    assert rebuilt == full.all_bootstrap

    summarized = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        _observations(con, dis1, dis2),
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={
            "iterations": iterations,
            "diagnostic": True,
            "summary_only": True,
        },
        p_value_correction="bfn",
        family_size=1,
        triplet_seed=11,
    )
    decisions = summarized.bootstrap_perm_decisions
    assert decisions["count"] == iterations
    assert (
        decisions["greater"] + decisions["less"] + decisions["inconclusive"]
        == iterations
    )
    assert decisions["greater"] == record["decision"].count("greater")
    assert summarized.bootstrap_perm_p_greater["count"] == iterations
    assert summarized.bootstrap_perm_p_greater["non_null_count"] == iterations


def test_studentized_interval_brackets_the_observed_statistic():
    """The percentile interval is ordered and covers the statistic it describes.

    Each iteration recomputes the studentized difference on its own resample, so
    the interval is centred on the observed value. With too few observations for
    a variance in both groups it is left unreported rather than invented.
    """
    con_heights = [0.50 + 0.01 * i for i in range(60)]
    dis1_heights = [0.20 + 0.01 * i for i in range(40)]
    result = _corrected(
        _observations(con_heights, dis1_heights, [0.3] * 5), "no", iterations=200
    )
    assert result.bootstrap_perm_stat_ci_low < result.bootstrap_perm_stat_ci_high
    assert result.bootstrap_perm_stat_ci_low <= result.perm_statistic <= result.bootstrap_perm_stat_ci_high

    degenerate = _corrected(_observations([0.9] * 2, [0.2], [0.3]), "no", iterations=10)
    assert degenerate.bootstrap_perm_stat_ci_low is None
    assert degenerate.bootstrap_perm_stat_ci_high is None


@pytest.mark.output
@pytest.mark.parametrize(
    "method, expected",
    [
        ("no", False),
        ("bfn", True),
        ("fdr_bh", True),
    ],
)
def test_results_tsv_carries_corrected_columns_only_when_correcting(
    method, expected, tmp_path
):
    """The corrected p-value columns appear only when a correction is in effect.

    Under ``no`` the corrected value equals the raw one by definition, so the
    four columns would repeat their neighbours and suggest an adjustment that
    never happened. The raw columns and the significance flags are written
    either way, since those are what the classification rests on.
    """
    result = _corrected(_observations([0.9] * 25, [0.35] * 20, [0.3] * 4), method)
    path = tmp_path / "results.tsv"
    pinf.write_pipeline_results([result], str(path), p_value_correction=method)

    lines = path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    corrected_columns = [
        f"dct_p_val_{method}_corr",
        f"ks_p_val_{method}_corr",
        f"perm_p_greater_{method}_corr",
        f"perm_p_less_{method}_corr",
    ]
    assert all((column in header) is expected for column in corrected_columns)
    assert not any(column.endswith("_no_corr") for column in header)

    # Whatever the method, the header and the row stay aligned.
    assert len(lines[1].split("\t")) == len(header)
    for column in ("dct_p_value", "ks_p_value", "dct_significant", "ks_significant"):
        assert column in header


# One observation set per cascade outcome, so a whole-run family spans triplets
# settled at each of the three gates rather than only the cheap ones.
_CASCADE_FAMILY = [
    ([0.1] * 20, [0.9] * 10, [0.9] * 10),
    ([0.5] * 10, [0.5] * 30, [0.5] * 2),
    (_HIGH, _LOW, [0.1] * 2),
    (_LOW[:10], _HIGH * 3, [0.9] * 2),
    (_WIDE, _NARROW, [0.5] * 2),
]


# Everything the cascade reads, which the ``diagnostic`` setting must not move.
# The raw ``ks_statistic``/``ks_p_value`` are deliberately absent: under an
# inline correction a non-diagnostic run does not measure them below a settled
# count gate, and nothing reads them there.
_DIAGNOSTIC_INVARIANT_FIELDS = (
    "n_con",
    "n_dis1",
    "n_dis2",
    "dis1_topology",
    "dct_statistic",
    "dct_p_value",
    "dct_p_value_corrected",
    "dct_significant",
    "classification",
    "decision_gate",
    "bootstrap_value",
    "all_bootstrap",
)


@pytest.mark.parametrize("method", ["bfn", "holm", "fdr_bh"])
def test_diagnostic_changes_what_is_measured_and_nothing_concluded(method):
    """Skipping the tests the cascade cannot consult changes no conclusion.

    Skipping the direction test is licensed by every correction being monotone
    and by the permutation p-values being corrected inside the test rather than
    across triplets: a gate that failed raw cannot clear once corrected, and an
    unrun test moves no other triplet's numbers. A run with ``diagnostic`` off
    therefore has to agree field for field with a diagnostic one on everything
    the cascade reads, across an inline correction and two deferred ones, with
    the bootstrap votes included, those are judged against the corrected
    threshold, so a skip that shifted a gate would move them too.

    The diagnostic family, where every raw value exists, also shows what a
    corrected column is: each test's p-values form one family over all
    triplets, whatever the other test said, so a corrected column is the plain
    whole-run correction of its raw column, and blanking every count p-value
    leaves the tree-height column untouched. Under a rank-based method the
    tree-height column must agree everywhere between the two runs, since every
    triplet is a member the default still measures. Under ``bfn`` the default
    leaves it unmeasured below a settled count gate, and the survivors must
    still be corrected by the triplet count rather than by the number measured:
    a family shrunk to the survivors would lower their corrected values and
    move classifications past the tree-height gate.
    """
    measured = {}
    families = {}
    for diagnostic in (False, True):
        measured[diagnostic] = [
            _measure(
                _observations(*case), method, family_size=len(_CASCADE_FAMILY),
                iterations=30, seed=11, diagnostic=diagnostic,
            )
            for case in _CASCADE_FAMILY
        ]
        families[diagnostic] = _decide(measured[diagnostic], method)

    # The family has to exercise every gate, or the comparison would only prove
    # the two runs agree where neither of them resamples.
    assert {result.decision_gate for result in families[False]} == {
        "DCT", "THT", "PERM"
    }

    for lean, full in zip(families[False], families[True]):
        for field in _DIAGNOSTIC_INVARIANT_FIELDS:
            assert getattr(lean, field) == getattr(full, field), field

        if lean.ks_p_value is None:
            assert method == "bfn" and lean.decision_gate == "DCT"
            assert lean.ks_p_value_corrected is None
            assert lean.ks_significant is None
        else:
            assert lean.ks_p_value == full.ks_p_value
            assert lean.ks_p_value_corrected == full.ks_p_value_corrected
            assert lean.ks_significant is full.ks_significant
            if method == "bfn":
                assert lean.ks_p_value_corrected == pytest.approx(
                    min(1.0, lean.ks_p_value * len(_CASCADE_FAMILY))
                )

        if lean.decision_gate == "PERM":
            assert lean.perm_decision == full.perm_decision
            assert lean.perm_note == full.perm_note
        else:
            assert lean.perm_decision is None
            assert lean.perm_statistic is None
            assert lean.perm_note == pinf.PERM_NOTE_NOT_CONSULTED
            assert full.perm_decision is not None

    for test in ("dct", "ks"):
        raw = [getattr(result, f"{test}_p_value") for result in measured[True]]
        expected = pinf._adjust_p_values(raw, method=method, alpha=0.05)
        corrected = [
            getattr(result, f"{test}_p_value_corrected") for result in families[True]
        ]
        assert corrected == pytest.approx(expected), test

    blanked = _decide(
        [replace(result, dct_p_value=1.0) for result in measured[True]], method
    )
    assert [result.ks_p_value_corrected for result in blanked] == [
        result.ks_p_value_corrected for result in families[True]
    ]
    assert {result.classification for result in blanked} == {"no_introgression"}


@pytest.mark.parametrize("diagnostic", [False, True])
@pytest.mark.parametrize("method", ["bfn", "holm"])
def test_bootstrap_measures_the_tree_height_test_its_correction_reads(
    monkeypatch, method, diagnostic
):
    """The bootstrap measures the tree-height test exactly where its correction reads it.

    Under ``bfn`` a count gate that failed on the exactly corrected value
    classifies the iteration before the tree-height flag is read, so the test
    goes unmeasured there, exactly as often as the count gate fails, and never
    elsewhere. Under ``holm`` every iteration's value is a member of a family
    corrected by rank, so it is measured in every iteration. ``diagnostic``
    reaches only the point estimate: on, it adds the one measurement the point
    estimate would otherwise decline under ``bfn``, and the bootstrap's count is
    the same either way.
    """
    calls = []
    real_ks_test = pinf.run_two_sample_ks_test

    def counting_ks_test(*args, **kwargs):
        calls.append(None)
        return real_ks_test(*args, **kwargs)

    monkeypatch.setattr(pinf, "run_two_sample_ks_test", counting_ks_test)

    iterations = 20
    # 10 vs 10 discordant trees: the point estimate's count gate fails, and so
    # do most resamples'.
    result = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        _observations(*_CASCADE_FAMILY[0]),
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": iterations},
        p_value_correction=method,
        family_size=1,
        triplet_seed=11,
        diagnostic=diagnostic,
    )
    decided = pinf._apply_triplet_result_p_value_correction(
        [result], alpha_dct=0.05, alpha_ks=0.05, method=method
    )[0]
    cleared = round(iterations * (1.0 - decided.all_bootstrap.get("no_introgression", 0.0)))
    assert cleared < iterations

    point_estimate = 1 if diagnostic or method == "holm" else 0
    assert (result.ks_p_value is None) == (point_estimate == 0)
    if method == "bfn":
        assert len(calls) == point_estimate + cleared
    else:
        assert len(calls) == point_estimate + iterations
