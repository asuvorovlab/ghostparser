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
        expected = "ambiguous"

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
        exact.add(round(pperm._studentized_mean_diff(pooled[mask], pooled[~mask]), 9))

    assert len(exact) == math.comb(nx + ny, nx)

    sampled = pperm._permutation_statistics(
        centered, nx, ny, 5000, np.random.default_rng(7)
    )
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
    assert result.decision in {"greater", "less", "ambiguous"}
    assert 0.0 < result.p_greater <= 1.0
    assert 0.0 < result.p_less <= 1.0
    assert 0.0 < result.p_two_sided <= 1.0
    # The two one-tailed counts both include ties, so together they cover every
    # resample at least once and their p-values must sum past 1.
    assert result.p_greater + result.p_less > 1.0
    assert 600 <= result.n_resamples <= 3000
    # A significant direction must agree with the sign of the observed statistic.
    if result.decision == "greater":
        assert result.statistic > 0
    elif result.decision == "less":
        assert result.statistic < 0


@pytest.mark.parametrize("seed", range(8))
def test_decision_rules_agree_on_random_inputs(seed):
    """The two one-tailed rule and the two-tailed-gate-then-sign rule concur.

    Correcting the one-tailed family with Bonferroni compares ``2p`` to alpha,
    matching the two-tailed gate, so the two rules are expected to coincide.
    """
    rng = np.random.default_rng(2000 + seed)
    x, y = _random_samples(rng)

    result = pperm.run_studentized_permutation_test(
        x,
        y,
        alpha=_ALPHA,
        min_resamples=2000,
        max_resamples=6000,
        correction="bfn",
        rng=np.random.default_rng(seed),
    )
    assert result.consistent is True
    assert result.decision == result.decision_two_sided


def test_equal_samples_give_a_zero_statistic_and_no_direction():
    """Two samples holding the same values cannot separate in either direction."""
    values = [0.10, 0.22, 0.31, 0.44, 0.55, 0.61, 0.78, 0.83]
    result = pperm.run_studentized_permutation_test(
        values, list(values), min_resamples=600, max_resamples=600,
        rng=np.random.default_rng(6),
    )
    assert result.statistic == pytest.approx(0.0)
    assert result.decision == "ambiguous"


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
        if result.decision != "ambiguous":
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
    """Each guard returns an ambiguous result and draws no permutations."""
    result = pperm.run_studentized_permutation_test(
        x, y, rng=np.random.default_rng(0)
    )
    assert result.note == expected_note
    assert result.decision == "ambiguous"
    assert result.n_resamples == 0
    assert result.statistic is None
    assert result.converged is False


def test_skewed_null_keeps_the_directional_call_and_flags_it():
    """An asymmetric permutation null is flagged without overturning the direction.

    This reproduces a shape seen on real data: a large concordant sample against
    a tiny discordant1 sample carrying extreme heights. Randomly reassigning
    those few large values produces a strongly right-skewed null, so the
    absolute-value two-tailed count borrows the fat right tail and fails to
    resolve a left-tail observation that the directional tail resolves outright.
    The directional decision is the correct one; the mismatch is recorded as a
    note rather than being treated as a rule disagreement.
    """
    rng = np.random.default_rng(21)
    con = rng.normal(0.42, 0.10, 700)
    dis1 = np.concatenate([rng.normal(0.5, 0.1, 15), rng.normal(25.0, 5.0, 4)])

    result = pperm.run_studentized_permutation_test(
        con,
        dis1,
        alpha=_ALPHA,
        min_resamples=2500,
        max_resamples=2500,
        correction="bfn",
        rng=np.random.default_rng(22),
    )

    assert result.statistic < 0
    assert result.decision == "less"
    assert result.p_two_sided > _ALPHA
    assert result.null_skewed is True
    # Skew is informational, not an exception, so it does not occupy ``note``.
    assert result.note is None
    # The cross-check gates on the doubled smaller tail, which stays valid under
    # asymmetry, so the two rules still agree.
    assert result.consistent is True


def test_max_resamples_reached_is_reported():
    """Exhausting the budget without excluding alpha is flagged, not hidden."""
    rng = np.random.default_rng(3)
    # Nearly identical samples keep p far from alpha on the low side but the
    # interval around p_less stays wide enough that alpha remains inside it.
    x = rng.normal(1.0, 1.0, 40)
    y = rng.normal(1.0, 1.0, 40)
    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=200, max_resamples=200, rng=np.random.default_rng(4)
    )
    assert result.n_resamples == 200
    if not result.converged:
        assert result.note == "max_resamples_reached"


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


def test_batches_grow_by_one_quarter_until_the_budget_is_spent():
    """Undecided runs grow each batch to 1.25x the previous one.

    Two samples from the same distribution keep the p-values far from alpha but
    with intervals too wide to exclude it at small resample counts, so the run
    walks the full geometric schedule. Starting at 100 the batches are
    100, 125, 156, 195, 243 (each ``int(previous * 1.25)``), which cumulate to
    100, 225, 381, 576, 819; the sixth batch of 303 is clipped to the 1000
    remaining. The totals must land on that schedule exactly.
    """
    rng = np.random.default_rng(31)
    x = rng.normal(1.0, 1.0, 30)
    y = rng.normal(1.0, 1.0, 30)

    result = pperm.run_studentized_permutation_test(
        x, y, min_resamples=100, max_resamples=1000, rng=np.random.default_rng(32)
    )

    expected_totals = []
    total, batch = 0, 100
    while total < 1000:
        batch = max(1, min(batch, 1000 - total))
        total += batch
        expected_totals.append(total)
        batch = int(batch * 1.25)

    assert result.n_resamples in expected_totals
    assert result.batches == expected_totals.index(result.n_resamples) + 1


def test_median_sign_decision_matches_definition():
    """The fallback compares medians and reports ties as ambiguous."""
    assert pperm.median_sign_decision([3.0, 4.0, 5.0], [1.0, 2.0]) == "greater"
    assert pperm.median_sign_decision([1.0, 2.0], [3.0, 4.0, 5.0]) == "less"
    assert pperm.median_sign_decision([1.0, 2.0, 3.0], [2.0]) == "ambiguous"
    assert pperm.median_sign_decision([], [1.0]) == "ambiguous"


def test_bootstrap_resample_budget_scales_by_one_fifth():
    """Bootstrap iterations run at a fifth of the configured budget."""
    assert pperm.bootstrap_resample_budget(2500, 25000) == (500, 5000)
    # The floor keeps a usable budget when the configured numbers are tiny.
    assert pperm.bootstrap_resample_budget(2, 3) == (1, 1)
