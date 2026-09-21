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
def test_statistic_p_values_and_verdict_match_scipy(seed):
    """The statistic, both one-tailed p-values and the verdict agree with SciPy.

    The statistic must equal SciPy's for the same statistic function. The raw
    p-values are Monte Carlo estimates from two independent resampling
    streams, so they are compared within the binomial tolerance; the adaptive
    stopping rule is pinned off by making the minimum and maximum resample
    counts equal, so both implementations draw the same number of
    permutations. With Bonferroni over the one-tailed family GhostParser
    compares ``2p`` to ``alpha``, the same threshold as SciPy's raw
    ``p <= alpha / 2``, so the verdicts must agree wherever SciPy's p-value
    sits clear of that line by more than the Monte Carlo tolerance.
    """
    resamples = 4000
    rng = np.random.default_rng(seed)
    x, y = _random_samples(rng)

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        alpha=_ALPHA,
        min_resamples=resamples,
        max_resamples=resamples,
        correction="bfn",
        rng=np.random.default_rng(seed + 100),
    )
    assert result.n_resamples == resamples
    assert result.statistic == pytest.approx(float(_studentized(x, y)))

    theirs = {}
    for alternative, ours in (("greater", result.p_greater), ("less", result.p_less)):
        scipy_statistic, theirs[alternative] = _scipy_p_value(
            x, y, alternative, resamples, seed + 300
        )
        assert result.statistic == pytest.approx(scipy_statistic)
        assert ours == pytest.approx(
            theirs[alternative],
            abs=_monte_carlo_tolerance(theirs[alternative], resamples),
        ), alternative

    if all(
        abs(p - _ALPHA / 2) > _monte_carlo_tolerance(p, resamples)
        for p in theirs.values()
    ):
        if theirs["greater"] <= _ALPHA / 2:
            assert result.decision == "greater"
        elif theirs["less"] <= _ALPHA / 2:
            assert result.decision == "less"
        else:
            assert result.decision in {"equivalent", "inconclusive"}


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


def test_random_inputs_preserve_test_invariants():
    """Structural invariants hold for randomly shaped and distributed inputs.

    Twelve random pairs cover group sizes, spreads, families and separations no
    fixed fixture does. The same inputs under the same seed give byte-identical
    results, and two samples holding the same values cannot separate in either
    direction.
    """
    for seed in range(12):
        rng = np.random.default_rng(1000 + seed)
        x, y = _random_samples(rng)
        kwargs = dict(alpha=_ALPHA, min_resamples=600, max_resamples=3000)

        result = pperm.run_studentized_permutation_test(
            x, y, rng=np.random.default_rng(seed), **kwargs
        )

        assert result.note != "both_tails_significant", seed
        assert result.decision in {"greater", "less", "equivalent", "inconclusive"}, seed
        assert 0.0 < result.p_greater <= 1.0, seed
        assert 0.0 < result.p_less <= 1.0, seed
        # The two one-tailed counts both include ties, so together they cover
        # every resample at least once and their p-values must sum past 1.
        assert result.p_greater + result.p_less > 1.0, seed
        assert 600 <= result.n_resamples <= 3000, seed
        # A significant direction must agree with the sign of the statistic.
        if result.decision == "greater":
            assert result.statistic > 0, seed
        elif result.decision == "less":
            assert result.statistic < 0, seed

        assert result == pperm.run_studentized_permutation_test(
            x, y, rng=np.random.default_rng(seed), **kwargs
        ), seed

    values = [0.10, 0.22, 0.31, 0.44, 0.55, 0.61, 0.78, 0.83]
    result = pperm.run_studentized_permutation_test(
        values, list(values), min_resamples=600, max_resamples=600,
        rng=np.random.default_rng(6),
    )
    assert result.statistic == pytest.approx(0.0)
    assert result.decision in {"equivalent", "inconclusive"}


@pytest.mark.parametrize(
    "n, equivalence_test, expected",
    [(8, True, "inconclusive"), (400, True, "equivalent"), (400, False, "inconclusive")],
    ids=["too_little_data", "enough_data", "disabled"],
)
def test_equivalence_step_decides_from_the_directional_resamples(
    n, equivalence_test, expected
):
    """TOST separates "shown to be close" from "nothing shown", on the draws already made.

    Every case feeds the test two samples drawn from the same distribution, so
    neither direction can be significant and the equivalence step decides. The
    margin is an effect size (0.5 pooled standard deviations), so it shrinks
    relative to the standard error as n grows: 8 observations per group cannot
    rule out a medium effect, while 400 can. A margin expressed in
    standard-error units would report ``inconclusive`` at every n, since the
    studentized statistic is a pivot whose null spread does not shrink with
    sample size. The equivalence p-value is an add-one estimator over the
    directional test's own resamples, so it is a multiple of
    ``1 / (n_resamples + 1)`` and never sits below that floor; with the step
    switched off there is no p-value and the decision falls through.
    """
    rng = np.random.default_rng(11)
    x = rng.normal(1.0, 0.2, n)
    y = rng.normal(1.0, 0.2, n)
    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=1000, max_resamples=1000,
        equivalence_test=equivalence_test, rng=np.random.default_rng(12),
    )
    assert result.decision == expected

    if not equivalence_test:
        assert result.p_tost is None
        return
    assert (result.p_tost <= 0.05) is (expected == "equivalent")
    resolution = 1.0 / (result.n_resamples + 1)
    assert result.p_tost >= resolution
    assert result.p_tost / resolution == pytest.approx(
        round(result.p_tost / resolution), abs=1e-9
    )


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


@pytest.mark.parametrize(
    "seed, sizes, bound",
    [
        # A shape seen on real data: a large concordant sample against a tiny
        # discordant1 sample carrying a few extreme heights, which splits the
        # null into clusters by how many extremes land in the small group and
        # leaves it strongly asymmetric.
        (21, "asymmetric", ("gt", 1.0)),
        # Balanced samples from one symmetric family leave the null unskewed.
        (4, "symmetric", ("lt", 0.15)),
    ],
    ids=["asymmetric", "symmetric"],
)
def test_null_skewness_matches_scipy_and_the_null_shape(seed, sizes, bound):
    """The reported null skewness is SciPy's skewness of the drawn statistics.

    ``perm_null_skew`` is accumulated from running power sums so batches can
    be discarded, so it is checked against ``scipy.stats.skew`` over the same
    draws, and its magnitude against what the null's shape implies.
    """
    rng = np.random.default_rng(seed)
    if sizes == "asymmetric":
        x = rng.normal(0.42, 0.10, 700)
        y = np.concatenate([rng.normal(0.5, 0.1, 15), rng.normal(25.0, 5.0, 4)])
        resamples = 2500
    else:
        x = rng.normal(1.0, 1.0, 150)
        y = rng.normal(1.0, 1.0, 150)
        resamples = 4000

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        alpha=_ALPHA,
        min_resamples=resamples,
        max_resamples=resamples,
        correction="bfn",
        rng=np.random.default_rng(seed + 1),
    )

    pooled = np.concatenate([x, y])
    pooled = pooled - pooled.mean()
    drawn = pperm._permutation_statistics(
        pooled, x.size, y.size, resamples, np.random.default_rng(seed + 1)
    )[0]

    assert result.note is None
    assert result.n_resamples_skew == resamples
    assert result.null_skew == pytest.approx(float(stats.skew(drawn)), rel=1e-9)
    if bound[0] == "gt":
        assert abs(result.null_skew) > bound[1]
        assert result.statistic < 0
        assert result.decision == "less"
    else:
        assert abs(result.null_skew) < bound[1]


@pytest.mark.parametrize("separated", [True, False], ids=["converges", "exhausts_budget"])
def test_adaptive_run_converges_or_exhausts_its_budget(separated):
    """A clear separation converges on the first batch; a marginal one grows to the budget.

    A marginal shift keeps the corrected p-value close enough to alpha that the
    interval never excludes it, so the run draws every batch it is allowed and
    ends unconverged. ``max_resamples`` is the point at which it stops asking
    for more, not a hard cap: the batch that crosses the line is drawn at full
    size, so the total lands at or above the budget by at most one batch. The
    growth factor itself is a performance knob and is deliberately not pinned.
    """
    if separated:
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
        return

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


def test_bootstrap_resample_budget_shrinks_but_stays_usable_and_monotone():
    """The bootstrap budget shrinks the configured one without going unusable.

    Asserts the properties rather than the divisor: a bootstrap iteration must
    cost no more than the point estimate, must still draw at least one
    resample, must keep its bounds in order -- otherwise an adaptive run
    inside an iteration has no valid range to grow through -- and a larger
    configured budget must never yield a smaller one.
    """
    for min_resamples, max_resamples in (
        (2500, 25000), (2, 3), (1, 1), (100, 100), (7, 1000), (999, 1001),
    ):
        scaled_min, scaled_max = pperm.bootstrap_resample_budget(
            min_resamples, max_resamples
        )
        assert scaled_max <= max(max_resamples, scaled_min)
        assert scaled_min >= 1
        assert scaled_max >= scaled_min
        # Strictly smaller wherever there is room to divide, so a divisor of
        # 1 -- which would silently restore the full per-iteration cost --
        # fails here.
        if min_resamples >= 2:
            assert scaled_min < min_resamples

    previous = (0, 0)
    for configured in (1, 2, 5, 10, 100, 2500, 25000):
        current = pperm.bootstrap_resample_budget(configured, configured * 10)
        assert current[0] >= previous[0]
        assert current[1] >= previous[1]
        previous = current


def test_shifted_statistics_match_an_explicit_shift():
    """A shifted row equals re-running the draw on explicitly shifted data.

    The equivalence test rides on the directional test's permutations: rather
    than resampling shifted data, the kernel folds the shift into the power sums
    it already computed. Requesting the same draw count under the same seed
    yields the same permutations either way, so the two must agree to floating
    point, not merely in distribution. Both group orderings and several lopsided
    splits are covered because the kernel samples whichever group is smaller and
    recovers the other by subtraction.
    """
    for nx, ny in ((10, 30), (30, 10), (7, 7), (4, 25), (2, 60)):
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
            assert fused[row] == pytest.approx(expected, rel=1e-9, abs=1e-12), (nx, ny, row)
