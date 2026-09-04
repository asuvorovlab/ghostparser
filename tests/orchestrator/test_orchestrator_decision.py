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
from ghostparser.orchestrator.config import P_VALUE_CORRECTION_CHOICES

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


def _corrected(observations, method, family_size=1, iterations=40, seed=3):
    """Analyze one triplet and run the run-wide correction pass over it.

    Both correction tiers finish in this pass: the inline methods have already
    tallied their bootstrap votes, the rank-based ones are resolved here. When
    ``family_size`` exceeds 1 the triplet is padded out to that many results,
    each padding entry carrying a p-value of 1.0, so the point estimate is
    corrected against the same family the bootstrap used.

    Args:
        observations: The observation list.
        method: The p-value correction method.
        family_size: Triplet count the correction should work against.
        iterations: Bootstrap iterations to run.
        seed: Base seed for deterministic resampling.

    Returns:
        The corrected ``TripletPipelineResult`` for the real triplet.
    """
    result = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": iterations},
        p_value_correction=method,
        family_size=family_size,
        triplet_seed=seed,
    )
    # Padding stands in for triplets that never clear a gate. They need their own
    # deferred records, not an empty one: every triplet in a real run carries a
    # record, and the per-iteration family is sized from the records present, so
    # dropping them would silently under-correct the bootstrap.
    padding = [
        replace(
            result,
            dct_p_value=1.0,
            ks_p_value=1.0,
            all_bootstrap=None,
            bootstrap_deferred=pinf.DeferredBootstrapRecord(
                dct_p_values=np.ones(iterations, dtype=np.float64),
                ks_p_values=np.full(iterations, np.nan, dtype=np.float64),
                directions=np.zeros(iterations, dtype=np.int8),
            ),
        )
        for _ in range(family_size - 1)
    ]
    return pinf._apply_triplet_result_p_value_correction(
        [result] + padding, alpha_dct=0.05, alpha_ks=0.05, method=method
    )[0]


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
        # stopping the cascade before a later gate can be consulted.
        # The count gate fails, so the tree-height test joins no correction
        # family and its significance is undefined however the raw value fell.
        ([0.1] * 20, [0.9] * 10, [0.9] * 10, False, None, None, "DCT",
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
@pytest.mark.parametrize("pipeline_mode", ["efficient", "detailed"])
def test_decision_cascade_lands_on_each_classification(
    con, dis1, dis2, dct_significant, ks_significant, decisions, gate, expected,
    pipeline_mode,
):
    """Crafted observation sets drive the cascade onto each of its five outcomes.

    The gate asserts which test settled the call, so a case that reaches its
    classification by the wrong route fails rather than passing by coincidence.
    Both modes must agree on the classification and the gate: the efficient mode
    only declines to run a test whose result the cascade would have ignored, so
    it reports ``perm_decision`` exactly on the rows the permutation gate
    settled, while the detailed mode reports one everywhere.
    """
    result = _analyze(
        _observations(con, dis1, dis2),
        pipeline_mode=pipeline_mode,
    )
    assert result.dct_significant is dct_significant
    assert result.ks_significant is ks_significant
    assert result.classification == expected
    assert result.decision_gate == gate

    if pipeline_mode == "detailed" or gate == "PERM":
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

    Runs in the detailed mode so the test is reached regardless of what the
    earlier gates decided; a guard is a property of the samples, not of the
    cascade position.
    """
    result = _analyze(
        _observations(con, dis1, [0.1] * 2), pipeline_mode="detailed"
    )
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
    result = pinf.analyze_triplet(
        _TRIPLET,
        _subtrees([0.30] * 10, bc_heights, ac_heights),
        species_subtree=_SPECIES_SUBTREE,
        collect_summary_statistics=True,
        bootstrap_options={"iterations": 0},
        triplet_seed=1,
    )

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
    """
    assert pinf._classify_introgression(
        dct_significant, ks_significant, direction
    ) == (expected, expected_gate)


@pytest.mark.parametrize("method", ["no", "bfn", "holm", "fdr_bh", "fdr_by"])
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


@pytest.mark.parametrize("method", ["no", "bfn", "holm", "fdr_bh", "fdr_by"])
def test_inline_and_deferred_correction_agree_on_a_single_triplet(method):
    """Every short-circuiting method votes the same way on a family of one.

    A family of one leaves each correction as the identity, so all five methods
    must produce the same bootstrap tally. ``no`` and ``bfn`` reach it inline
    during the stream while the rank-based three park their raw p-values and are
    corrected afterwards, so agreeing here is what shows the two code paths
    implement one decision rule rather than two.
    """
    observations = _observations([0.9] * 25, [0.2] * 3 + [0.35] * 22, [0.3] * 4)
    result = _corrected(observations, method)
    baseline = _corrected(observations, "no")
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


def test_deferred_bootstrap_record_is_cleared_after_correction():
    """A rank-based method defers its votes and the correction pass resolves them.

    A family size above one is what makes the deferral observable: the
    tree-height family spans triplets, so a lone triplet is fully determined
    where it is analyzed and comes back already resolved.
    """
    observations = _observations([0.9] * 25, [0.2] * 3 + [0.35] * 22, [0.3] * 4)
    deferred = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": 20},
        p_value_correction="holm",
        family_size=2,
        triplet_seed=5,
    )
    assert deferred.all_bootstrap is None
    assert deferred.bootstrap_deferred is not None
    assert len(deferred.bootstrap_deferred.dct_p_values) == 20

    resolved = pinf._apply_triplet_result_p_value_correction(
        [deferred], alpha_dct=0.05, alpha_ks=0.05, method="holm"
    )[0]
    assert resolved.bootstrap_deferred is None
    assert sum(resolved.all_bootstrap.values()) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "dct_significant,ks_significant,direction",
    [
        (dct, ks, direction)
        for dct in (True, False)
        for ks in (True, False)
        for direction in ("greater", "less", "equivalent", "inconclusive")
    ],
)
def test_vectorized_bootstrap_codes_match_classify_introgression(
    dct_significant, ks_significant, direction
):
    """The array-form cascade agrees with the scalar one on every branch.

    ``_resolve_deferred_bootstrap`` classifies whole arrays at once rather than
    calling ``_classify_introgression`` per iteration, so the two must be pinned
    against each other or they can silently drift apart.
    """
    import numpy as np

    code = pinf._classification_codes(
        np.array([dct_significant]),
        np.array([ks_significant]),
        np.array([pinf._DIRECTION_CODES.get(direction, pinf._DIRECTION_SKIPPED)]),
    )[0]
    expected, _ = pinf._classify_introgression(
        dct_significant, ks_significant, direction
    )
    assert pinf._BOOTSTRAP_CLASSES[code] == expected


@pytest.mark.parametrize("method", P_VALUE_CORRECTION_CHOICES)
def test_every_supported_correction_is_monotone(method):
    """No supported correction may lower a p-value below its raw value.

    Both short-circuits — the point estimate skipping the direction test and a
    bootstrap iteration skipping it — rest on ``raw > alpha`` implying
    ``adjusted > alpha``. A method that can pull a p-value back under alpha
    would break both silently, so the property is asserted over the whole
    choice list rather than a fixed set of names. The family is built to be the
    hostile case: eight strong signals against two weak ones, which is where a
    true-null-count estimator would drive its multiplier below 1.
    """
    p_values = [0.001] * 8 + [0.4, 0.9]
    adjusted = pinf._adjust_p_values(p_values, method=method, alpha=0.05)

    assert all(a >= p - 1e-12 for a, p in zip(adjusted, p_values))


@pytest.mark.parametrize("family_size", [1, 7, 250])
def test_inline_bonferroni_matches_the_family_correction(family_size):
    """Correcting one p-value from the family size alone reproduces the full pass."""
    p_value = 0.004
    family = [p_value] + [0.5] * (family_size - 1)
    expected = pinf._adjust_p_values(family, method="bfn", alpha=0.05)[0]
    assert pinf._adjust_p_value_inline(p_value, "bfn", family_size) == pytest.approx(
        expected
    )


def test_inline_correction_rejects_a_rank_based_method():
    """A rank-based method cannot be applied without the rest of its family."""
    with pytest.raises(ValueError, match="needs the whole family"):
        pinf._adjust_p_value_inline(0.01, "holm", 10)


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
    assert result.bootstrap_stat_ci_low < result.bootstrap_stat_ci_high
    assert result.bootstrap_stat_ci_low <= result.perm_statistic <= result.bootstrap_stat_ci_high

    degenerate = _corrected(_observations([0.9] * 2, [0.2], [0.3]), "no", iterations=10)
    assert degenerate.bootstrap_stat_ci_low is None
    assert degenerate.bootstrap_stat_ci_high is None


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


# One observation set per cascade outcome, so a mode comparison spans triplets
# settled at each of the three gates rather than only the cheap ones.
_CASCADE_FAMILY = [
    ([0.1] * 20, [0.9] * 10, [0.9] * 10),
    ([0.5] * 10, [0.5] * 30, [0.5] * 2),
    (_HIGH, _LOW, [0.1] * 2),
    (_LOW[:10], _HIGH * 3, [0.9] * 2),
    (_WIDE, _NARROW, [0.5] * 2),
]

# Everything the cascade reads, which both modes must agree on exactly. The raw
# ``ks_statistic``/``ks_p_value`` are deliberately absent: the efficient mode
# does not measure them below a settled count gate, and they take no part in the
# correction family there.
_MODE_INVARIANT_FIELDS = (
    "n_con",
    "n_dis1",
    "n_dis2",
    "dis1_topology",
    "dct_statistic",
    "dct_p_value",
    "dct_p_value_corrected",
    "dct_significant",
    "ks_p_value_corrected",
    "ks_significant",
    "classification",
    "decision_gate",
    "bootstrap_value",
    "all_bootstrap",
)


@pytest.mark.parametrize("method", ["bfn", "holm"])
def test_efficient_and_detailed_modes_agree_on_every_classification(method):
    """The efficient mode changes what is computed, never what is concluded.

    Skipping the direction test is licensed by every correction being monotone,
    so a gate that failed raw cannot clear once corrected. Both modes therefore
    have to agree field for field on everything the cascade reads, across an
    inline correction and a deferred one, with the bootstrap votes included --
    those are judged against the corrected threshold, so a mode that shifted a
    gate would move them too. The tree-height correction family is the triplets
    the count gate cleared, which is a property of the results rather than of the
    mode, so both runs correct it over the same members.

    The only permitted differences are the ``perm_*`` block and a raw
    ``ks_p_value`` below a settled count gate, both of which the efficient mode
    declines to compute because nothing reads them.
    """
    families = {}
    for mode in ("efficient", "detailed"):
        results = [
            pinf.analyze_triplet_from_observations(
                _TRIPLET,
                _observations(*case),
                species_subtree=_SPECIES_SUBTREE,
                bootstrap_options={"iterations": 30},
                p_value_correction=method,
                family_size=len(_CASCADE_FAMILY),
                triplet_seed=11,
                pipeline_mode=mode,
            )
            for case in _CASCADE_FAMILY
        ]
        families[mode] = pinf._apply_triplet_result_p_value_correction(
            results, alpha_dct=0.05, alpha_ks=0.05, method=method
        )

    # The family has to exercise the permutation gate, or the comparison would
    # only prove the two modes agree where neither of them resamples.
    assert "PERM" in {result.decision_gate for result in families["efficient"]}

    for efficient, detailed in zip(families["efficient"], families["detailed"]):
        for field in _MODE_INVARIANT_FIELDS:
            assert getattr(efficient, field) == getattr(detailed, field), field

        if efficient.decision_gate == "PERM":
            assert efficient.perm_decision == detailed.perm_decision
            assert efficient.perm_note == detailed.perm_note
        else:
            assert efficient.perm_decision is None
            assert efficient.perm_statistic is None
            assert efficient.perm_note == pinf.PERM_NOTE_NOT_CONSULTED
            assert detailed.perm_decision is not None


@pytest.mark.parametrize("pipeline_mode", ["efficient", "detailed"])
def test_tree_height_family_holds_only_the_count_gate_survivors(pipeline_mode):
    """The tree-height correction family is the triplets that cleared gate one.

    A triplet the count test already settled contributes nothing the cascade
    reads, so enrolling its tree-height p-value would inflate the family and push
    the corrected values of the triplets that do decide something towards
    non-significance. The family is therefore the count-gate survivors, and it is
    the same set in both pipeline modes: the detailed mode still measures the
    others, and they still take no part in the correction.
    """
    # Two triplets the count gate settles (10 vs 10 gives chi-square p = 1.0),
    # one it does not.
    cases = [
        ([0.1] * 20, [0.9] * 10, [0.9] * 10),
        ([0.1] * 20, [0.9] * 10, [0.9] * 10),
        (_HIGH, _LOW, [0.1] * 2),
    ]
    results = [
        pinf.analyze_triplet_from_observations(
            _TRIPLET,
            _observations(*case),
            species_subtree=_SPECIES_SUBTREE,
            bootstrap_options={"iterations": 0},
            p_value_correction="bfn",
            family_size=len(cases),
            triplet_seed=5,
            pipeline_mode=pipeline_mode,
        )
        for case in cases
    ]
    corrected = pinf._apply_triplet_result_p_value_correction(
        results, alpha_dct=0.05, alpha_ks=0.05, method="bfn"
    )

    for settled in corrected[:2]:
        assert settled.dct_significant is False
        assert settled.ks_p_value_corrected is None
        assert settled.ks_significant is None
        assert settled.classification == "no_introgression"

    # One family member, so Bonferroni multiplies by 1 rather than by 3.
    survivor = corrected[2]
    assert survivor.dct_significant is True
    assert survivor.ks_p_value_corrected == pytest.approx(survivor.ks_p_value)
