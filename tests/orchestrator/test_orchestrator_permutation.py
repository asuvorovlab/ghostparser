"""Verify the adaptive studentized permutation test against SciPy and by definition.

Three kinds of check live here:

- **Parity.** The observed statistic and both one-tailed p-values are compared
  against ``scipy.stats.permutation_test`` driven with the same studentized
  statistic. SciPy tests one alternative per call, so parity needs two calls
  where GhostParser accumulates both tails from a single set of resamples.
  p-values are Monte Carlo estimates from two independent resampling streams, so
  they are compared within a tolerance derived from the binomial standard error
  rather than for exact equality.
- **Sampler correctness.** The vectorized sampler skips materializing the larger
  group and reconstructs it by subtracting from the pooled totals; a small case
  is checked against exhaustive enumeration of every group assignment, computed
  with the scalar reference statistic.
- **Randomized inputs.** Group sizes, distributions, and separations are drawn
  at random so the invariants are exercised on shapes no fixed fixture covers.

See ``tests/TEST_IO.md`` for the derivation of each expected value.
"""

import math
from itertools import combinations

import numpy as np
import pytest
from scipy import stats

from ghostparser.orchestrator import permutation as pperm

_ALPHA = 0.05


def _studentized(x, y, axis=-1):
    """Welch-studentized mean difference, vectorized along ``axis`` for SciPy.

    Args:
        x: First sample.
        y: Second sample.
        axis: Axis along which the samples lie.

    Returns:
        The studentized statistic, broadcast over the remaining axes.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx = x.shape[axis]
    ny = y.shape[axis]
    numerator = np.mean(x, axis=axis) - np.mean(y, axis=axis)
    denominator = np.sqrt(
        np.var(x, axis=axis, ddof=1) / nx + np.var(y, axis=axis, ddof=1) / ny
    )
    return numerator / denominator


def _scipy_p_value(x, y, alternative, n_resamples, seed):
    """Run SciPy's permutation test with the same statistic and alternative.

    Args:
        x: First sample.
        y: Second sample.
        alternative: ``greater`` or ``less``.
        n_resamples: Number of permutations.
        seed: Seed for SciPy's resampling stream.

    Returns:
        A tuple ``(statistic, p_value)``.
    """
    result = stats.permutation_test(
        (np.asarray(x, dtype=float), np.asarray(y, dtype=float)),
        _studentized,
        permutation_type="independent",
        alternative=alternative,
        n_resamples=n_resamples,
        vectorized=True,
        rng=np.random.default_rng(seed),
    )
    return float(result.statistic), float(result.pvalue)


def _monte_carlo_tolerance(p_value, n_resamples, sigmas=5.0):
    """Bound the gap between two independent Monte Carlo estimates of one p-value.

    Each estimate has binomial standard error ``sqrt(p(1-p)/n)``; the difference
    of two independent estimates has ``sqrt(2)`` times that. A 5-sigma band keeps
    the test from flaking while still failing on any systematic disagreement.

    Args:
        p_value: The reference p-value.
        n_resamples: Permutations behind each estimate.
        sigmas: Width of the band in standard errors.

    Returns:
        The absolute tolerance.
    """
    standard_error = math.sqrt(max(p_value * (1.0 - p_value), 1e-9) / n_resamples)
    return sigmas * math.sqrt(2.0) * standard_error


def _random_samples(rng):
    """Draw a random pair of samples with random sizes, spreads, and separation.

    Args:
        rng: A ``numpy.random.Generator``.

    Returns:
        A tuple ``(x, y)`` of float arrays.
    """
    nx = int(rng.integers(8, 120))
    ny = int(rng.integers(8, 120))
    shift = float(rng.uniform(-1.0, 1.0))
    scale_x = float(rng.uniform(0.2, 2.0))
    scale_y = float(rng.uniform(0.2, 2.0))

    family = int(rng.integers(0, 3))
    if family == 0:
        x = rng.normal(1.0, scale_x, nx)
        y = rng.normal(1.0 + shift, scale_y, ny)
    elif family == 1:
        x = rng.lognormal(0.0, scale_x, nx)
        y = rng.lognormal(shift, scale_y, ny)
    else:
        x = rng.exponential(scale_x, nx)
        y = rng.exponential(scale_y, ny) + shift
    return x, y


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_observed_statistic_matches_scipy(seed):
    """The observed statistic equals SciPy's for the same statistic function."""
    rng = np.random.default_rng(seed)
    x, y = _random_samples(rng)

    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=500, max_resamples=500, rng=np.random.default_rng(seed + 100)
    )
    scipy_statistic, _ = _scipy_p_value(x, y, "greater", 500, seed + 200)

    assert result.statistic == pytest.approx(scipy_statistic)
    assert result.statistic == pytest.approx(float(_studentized(x, y)))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5, 6, 7])
def test_p_values_match_scipy_within_monte_carlo_error(seed):
    """Both one-tailed p-values agree with SciPy's single-alternative runs.

    Correction is disabled so the comparison is against SciPy's raw p-values;
    the adaptive stopping rule is pinned off by making the minimum and maximum
    resample counts equal, so both implementations draw the same number of
    permutations.
    """
    resamples = 4000
    rng = np.random.default_rng(seed)
    x, y = _random_samples(rng)

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        min_resamples=resamples,
        max_resamples=resamples,
        correction="no",
        rng=np.random.default_rng(seed + 100),
    )
    assert result.n_resamples == resamples

    for alternative, ours in (
        ("greater", result.p_greater),
        ("less", result.p_less),
    ):
        _, theirs = _scipy_p_value(x, y, alternative, resamples, seed + 300)
        assert ours == pytest.approx(
            theirs, abs=_monte_carlo_tolerance(theirs, resamples)
        ), alternative


@pytest.mark.parametrize("seed", [11, 12, 13])
def test_decision_matches_scipy_directional_verdict(seed):
    """The directional decision agrees with SciPy's two one-tailed verdicts.

    With Bonferroni over the one-tailed family, GhostParser compares ``2p`` to
    ``alpha``, which is the same threshold as SciPy's raw ``p <= alpha / 2``.
    """
    resamples = 4000
    rng = np.random.default_rng(seed)
    x = rng.normal(1.0, 0.4, 60)
    y = rng.normal(0.4, 0.6, 45)

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        alpha=_ALPHA,
        min_resamples=resamples,
        max_resamples=resamples,
        correction="bfn",
        rng=np.random.default_rng(seed + 100),
    )

    _, scipy_greater = _scipy_p_value(x, y, "greater", resamples, seed + 300)
    _, scipy_less = _scipy_p_value(x, y, "less", resamples, seed + 400)
    if scipy_greater <= _ALPHA / 2:
        expected = "greater"
    elif scipy_less <= _ALPHA / 2:
        expected = "less"
    else:
        expected = None

    if expected is None:
        assert result.decision in {"equivalent", "inconclusive"}
    else:
        assert result.decision == expected


def test_permutation_statistics_match_exhaustive_enumeration():
    """The optimized sampler only ever emits statistics from the exact null set.

    The sampler draws the smaller group and reconstructs the larger one by
    subtracting from the pooled sum and sum-of-squares. Enumerating all
    ``C(9, 4) = 126`` assignments and computing each with the scalar reference
    statistic gives the exact support; every sampled value must fall in it, and
    over enough draws every support value must appear.
    """
    x = np.array([0.11, 0.24, 0.37, 0.52])
    y = np.array([0.63, 0.71, 0.88, 0.95, 1.10])
    nx, ny = x.size, y.size
    pooled = np.concatenate([x, y])
    centered = pooled - pooled.mean()

    exact = set()
    for picked in combinations(range(pooled.size), nx):
        mask = np.zeros(pooled.size, dtype=bool)
        mask[list(picked)] = True
        exact.add(round(pperm.studentized_mean_diff(pooled[mask], pooled[~mask]), 9))

    assert len(exact) == math.comb(nx + ny, nx)

    # Row 0 is the unshifted draw; no equivalence shifts are requested here.
    sampled = pperm._permutation_statistics(
        centered, nx, ny, 5000, np.random.default_rng(7)
    )[0]
    observed = {round(float(value), 9) for value in sampled}

    assert observed <= exact
    # 5000 draws over 126 equally likely assignments leaves a vanishing chance
    # of missing any one of them.
    assert observed == exact


@pytest.mark.parametrize("seed", range(12))
def test_random_inputs_preserve_test_invariants(seed):
    """Structural invariants hold for randomly shaped and distributed inputs."""
    rng = np.random.default_rng(1000 + seed)
    x, y = _random_samples(rng)

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        alpha=_ALPHA,
        min_resamples=600,
        max_resamples=3000,
        rng=np.random.default_rng(seed),
    )

    assert result.note != "both_tails_significant"
    assert result.decision in {"greater", "less", "equivalent", "inconclusive"}
    assert 0.0 < result.p_greater <= 1.0
    assert 0.0 < result.p_less <= 1.0
    # The two one-tailed counts both include ties, so together they cover every
    # resample at least once and their p-values must sum past 1.
    assert result.p_greater + result.p_less > 1.0
    assert 600 <= result.n_resamples <= 3000
    # A significant direction must agree with the sign of the observed statistic.
    if result.decision == "greater":
        assert result.statistic > 0
    elif result.decision == "less":
        assert result.statistic < 0


def test_equal_samples_give_a_zero_statistic_and_no_direction():
    """Two samples holding the same values cannot separate in either direction."""
    values = [0.10, 0.22, 0.31, 0.44, 0.55, 0.61, 0.78, 0.83]
    result = pperm.run_studentized_permutation_test(
        values, list(values), min_resamples=600, max_resamples=600,
        rng=np.random.default_rng(6),
    )
    assert result.statistic == pytest.approx(0.0)
    assert result.decision in {"equivalent", "inconclusive"}


@pytest.mark.parametrize(
    "n, expected",
    [(8, "inconclusive"), (400, "equivalent")],
)
def test_equivalence_needs_enough_data_to_conclude(n, expected):
    """TOST separates "shown to be close" from "nothing shown" as data accrues.

    Both cases feed the test two samples drawn from the same distribution, so
    neither direction can be significant and the equivalence step decides. The
    margin is an effect size (0.5 pooled standard deviations), so it shrinks
    relative to the standard error as n grows: 8 observations per group cannot
    rule out a medium effect, while 400 can. A margin expressed in standard-error
    units would report ``inconclusive`` at every n, since the studentized
    statistic is a pivot whose null spread does not shrink with sample size.
    """
    rng = np.random.default_rng(11)
    x = rng.normal(1.0, 0.2, n)
    y = rng.normal(1.0, 0.2, n)
    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=1000, max_resamples=1000,
        rng=np.random.default_rng(12),
    )
    assert result.decision == expected
    assert (result.p_tost <= 0.05) is (expected == "equivalent")


def test_type_one_error_rate_tracks_alpha_under_unequal_variance():
    """Under the null the test rejects at about alpha, not more.

    Both samples share a mean but differ in size and spread (n=60 at sd 1.0
    against n=180 at sd 0.3). Unequal sizes paired with unequal variances are
    exactly where a permutation test of the raw mean difference loses its
    nominal level; holding level here is what the Welch studentization buys.

    Over 300 null replicates at alpha=0.05 the expected count is 15 with a
    standard deviation of 3.8. The assertion band spans 1% to 10% -- wide enough
    that fixed seeds make it stable, tight enough to catch a test that has
    stopped controlling its error rate. The level does degrade once the smaller
    group falls below roughly 30 observations *and* carries the larger spread;
    that limitation is documented in ORCHESTRATOR.md rather than asserted here.
    """
    replicates = 300
    rejections = 0
    for replicate in range(replicates):
        rng = np.random.default_rng(7000 + replicate)
        x = rng.normal(1.0, 1.0, 60)
        y = rng.normal(1.0, 0.3, 180)
        result = pperm.run_studentized_permutation_test(
            x,
            y,
            alpha=_ALPHA,
            min_resamples=1000,
            max_resamples=1000,
            correction="bfn",
            rng=np.random.default_rng(replicate),
        )
        if result.decision in {"greater", "less"}:
            rejections += 1

    assert 3 <= rejections <= 30, rejections


def test_seeded_runs_are_reproducible():
    """The same inputs and seed give byte-identical results."""
    rng = np.random.default_rng(9)
    x, y = _random_samples(rng)
    kwargs = dict(min_resamples=800, max_resamples=2000)

    first = pperm.run_studentized_permutation_test(
        x, y, rng=np.random.default_rng(42), **kwargs
    )
    second = pperm.run_studentized_permutation_test(
        x, y, rng=np.random.default_rng(42), **kwargs
    )
    assert first == second


@pytest.mark.parametrize(
    "x,y,expected_note",
    [
        ([0.4], [0.1, 0.2, 0.3], "insufficient_group_size"),
        ([0.3, 0.3, 0.3], [0.3, 0.3, 0.3], "zero_pooled_variance"),
        ([0.9] * 10, [0.1] * 30, "degenerate_observed_scale"),
        ([0.9, 0.8, 0.7, 0.6], [0.1, 0.2], "insufficient_permutation_support"),
    ],
)
def test_guards_short_circuit_without_resampling(x, y, expected_note):
    """Each guard returns an inconclusive result and draws no permutations."""
    result = pperm.run_studentized_permutation_test(
        x, y, rng=np.random.default_rng(0)
    )
    assert result.note == expected_note
    assert result.decision == "inconclusive"
    assert result.n_resamples == 0
    assert result.statistic is None
    assert result.converged is False


def test_null_skewness_is_measured_and_matches_scipy():
    """The reported null skewness equals the skewness of the drawn statistics.

    ``perm_null_skew`` is accumulated from running power sums so batches can be
    discarded, so it is checked against `scipy.stats.skew` over the same draws.
    The fixture is a shape seen on real data — a large concordant sample against
    a tiny discordant1 sample carrying a few extreme heights — which splits the
    null into clusters by how many extremes land in the small group and leaves
    it strongly asymmetric.
    """
    rng = np.random.default_rng(21)
    con = rng.normal(0.42, 0.10, 700)
    dis1 = np.concatenate([rng.normal(0.5, 0.1, 15), rng.normal(25.0, 5.0, 4)])
    resamples = 2500

    result = pperm.run_studentized_permutation_test(
        con,
        dis1,
        alpha=_ALPHA,
        min_resamples=resamples,
        max_resamples=resamples,
        correction="bfn",
        rng=np.random.default_rng(22),
    )

    pooled = np.concatenate([con, dis1])
    pooled = pooled - pooled.mean()
    drawn = pperm._permutation_statistics(
        pooled, con.size, dis1.size, resamples, np.random.default_rng(22)
    )[0]

    assert result.statistic < 0
    assert result.decision == "less"
    assert result.n_resamples_skew == resamples
    assert result.null_skew == pytest.approx(float(stats.skew(drawn)), rel=1e-9)
    # This fixture's null is strongly asymmetric; a symmetric one sits near 0.
    assert abs(result.null_skew) > 1.0
    assert result.note is None


def test_null_skewness_is_near_zero_for_a_symmetric_null():
    """Balanced samples from one symmetric family leave the null unskewed."""
    rng = np.random.default_rng(4)
    x = rng.normal(1.0, 1.0, 150)
    y = rng.normal(1.0, 1.0, 150)
    result = pperm.run_studentized_permutation_test(
        x, y, alpha=_ALPHA, min_resamples=4000,
        max_resamples=4000, rng=np.random.default_rng(5),
    )
    assert abs(result.null_skew) < 0.15


def test_adaptive_run_grows_batches_until_it_converges():
    """A clearly separated pair converges on the first batch."""
    rng = np.random.default_rng(8)
    x = rng.normal(2.0, 0.2, 80)
    y = rng.normal(0.5, 0.2, 80)
    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=1000, max_resamples=20000, rng=np.random.default_rng(9)
    )
    assert result.converged is True
    assert result.batches == 1
    assert result.n_resamples == 1000
    assert result.decision == "greater"


def test_undecided_runs_grow_their_batches_until_the_budget_is_reached():
    """An undecided run escalates its batches and stops once the budget is met.

    A marginal shift keeps the corrected p-value close enough to alpha that the
    interval never excludes it, so the run draws every batch it is allowed and
    ends unconverged. ``max_resamples`` is the point at which it stops asking
    for more, not a hard cap: the batch that crosses the line is drawn at full
    size, so the total lands at or above the budget by at most one batch. The
    growth factor itself is a performance knob and is deliberately not pinned.
    """
    rng = np.random.default_rng(3)
    x = rng.normal(1.4, 1.0, 30)
    y = rng.normal(1.0, 1.0, 30)

    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=100, max_resamples=1000, rng=np.random.default_rng(103)
    )

    assert result.batches > 1
    # Growth: a schedule that repeated the opening batch would total exactly
    # batches * 100, so a larger total shows the batches grew.
    assert result.n_resamples > result.batches * 100
    # The run stops only once the budget is met, and overshoots by at most the
    # final (largest) batch rather than being trimmed to land on it exactly.
    assert result.n_resamples >= 1000
    assert result.n_resamples < 2 * 1000
    assert result.converged is False
    assert result.note == "max_resamples_reached"


def test_bootstrap_resample_budget_scales_by_one_fifth():
    """Bootstrap iterations run at a fifth of the configured budget."""
    assert pperm.bootstrap_resample_budget(2500, 25000) == (500, 5000)
    # The floor keeps a usable budget when the configured numbers are tiny.
    assert pperm.bootstrap_resample_budget(2, 3) == (1, 1)


@pytest.mark.parametrize("nx, ny", [(10, 30), (30, 10), (7, 7), (4, 25), (2, 60)])
def test_shifted_statistics_match_an_explicit_shift(nx, ny):
    """A shifted row equals re-running the draw on explicitly shifted data.

    The equivalence test rides on the directional test's permutations: rather
    than resampling shifted data, the kernel folds the shift into the power sums
    it already computed. Requesting the same draw count under the same seed
    yields the same permutations either way, so the two must agree to floating
    point, not merely in distribution. Both group orderings and several lopsided
    splits are covered because the kernel samples whichever group is smaller and
    recovers the other by subtraction.
    """
    rng = np.random.default_rng(4)
    x = rng.gamma(2.0, 1.0, nx)
    y = rng.gamma(2.5, 1.2, ny)
    pooled = np.concatenate([x, y])
    pooled = pooled - pooled.mean()
    shift = 0.75

    fused = pperm._permutation_statistics(
        pooled, nx, ny, 2000, np.random.default_rng(31), shifts=(shift, -shift)
    )

    for row, value in enumerate((0.0, shift, -shift)):
        shifted = pooled.copy()
        shifted[:nx] += value
        shifted = shifted - shifted.mean()
        expected = pperm._permutation_statistics(
            shifted, nx, ny, 2000, np.random.default_rng(31)
        )[0]
        assert fused[row] == pytest.approx(expected, rel=1e-9, abs=1e-12)


def test_equivalence_reuses_the_directional_resamples():
    """TOST answers at the resample count the directional test settled on.

    The equivalence p-value is an add-one estimator over the same draws, so it
    is a multiple of ``1 / (n_resamples + 1)`` and can never sit below that
    floor. Fusing the two also means no extra resampling happens for it.
    """
    rng = np.random.default_rng(8)
    # Same mean, so no direction is significant and the equivalence step runs.
    x = rng.normal(0.0, 1.0, 40)
    y = rng.normal(0.0, 1.0, 40)

    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=2000, max_resamples=2000, rng=np.random.default_rng(2)
    )

    assert result.decision in {"equivalent", "inconclusive"}
    assert result.p_tost is not None
    resolution = 1.0 / (result.n_resamples + 1)
    assert result.p_tost >= resolution
    assert result.p_tost * (result.n_resamples + 1) == pytest.approx(
        round(result.p_tost * (result.n_resamples + 1)), abs=1e-9
    )


def test_equivalence_test_disabled_reports_no_tost():
    """With the equivalence step off the decision falls through to inconclusive."""
    rng = np.random.default_rng(8)
    x = rng.normal(0.0, 1.0, 40)
    y = rng.normal(0.0, 1.0, 40)

    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=2000, max_resamples=2000,
        equivalence_test=False, rng=np.random.default_rng(2),
    )

    assert result.p_tost is None
    assert result.decision == "inconclusive"
