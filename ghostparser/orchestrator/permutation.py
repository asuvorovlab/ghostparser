"""Adaptive studentized permutation test for concordant vs. discordant1 heights.

Decides whether the mean concordant tree height is greater than, less than, or
equivalent to the mean discordant1 height, resampling until the decision is
resolved to within Monte Carlo error.

The method, its citations, and the sampling optimization are written up under
"Gate 3 - Adaptive studentized permutation test" in the orchestrator guide.
"""

import math
from dataclasses import dataclass

import numpy as np
from statsmodels.stats.proportion import proportion_confint

from .config import (
    DEFAULT_ALPHA_PERM,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_PERMUTATION_CI_METHOD,
    DEFAULT_PERMUTATION_MAX_RESAMPLES,
    DEFAULT_PERMUTATION_MIN_RESAMPLES,
)
from .correction import adjust_p_values

__all__ = [
    "EQUIVALENCE_DELTA",
    "PermutationTestResult",
    "bootstrap_resample_budget",
    "run_studentized_permutation_test",
    "studentized_mean_diff",
]

# Decision labels. ``greater``/``less`` refer to the first sample (concordant)
# relative to the second (discordant1). When neither direction is supported the
# equivalence test splits the outcome: ``equivalent`` means the means were shown
# to differ by less than the equivalence margin, ``inconclusive`` means nothing
# was shown either way.
DECISION_GREATER = "greater"
DECISION_LESS = "less"
DECISION_EQUIVALENT = "equivalent"
DECISION_INCONCLUSIVE = "inconclusive"

# TOST equivalence margin, in pooled-standard-deviation units (a Cohen's d of
# 0.5). An effect size, not a multiple of the standard error: the studentized
# statistic is a pivot whose null spread stays near 1 at every sample size, so an
# SE-based margin could never be rejected; see "Gate 3" in the orchestrator
# guide for the measurements behind that.
EQUIVALENCE_DELTA = 0.5

# Each adaptive batch after the first is this multiple of the previous one, so
# a run that is still undecided spends geometrically more effort per check
# instead of paying the confidence-interval overhead on every few resamples.
_BATCH_GROWTH_FACTOR = 1.25

# Confidence level of the stopping rule's interval. An internal precision knob,
# not a statistical choice, so it is fixed rather than configurable.
_CI_LEVEL = 0.95

# Upper bound on the number of random keys held in memory at once. The sampler
# materializes a (batch x n_total) float64 array of sort keys, so this caps a
# single chunk near 16 MB regardless of how large the requested batch is. That
# matters because every pool worker runs its own permutation tests.
_MAX_PERMUTATION_CELLS = 2_000_000

# Above this pooled size the number of distinct group assignments C(n, k) is
# astronomically larger than any resample budget, so the exact support check is
# skipped rather than evaluating a needlessly large binomial coefficient.
_EXACT_SUPPORT_CHECK_MAX_N = 40

# A standard error below this fraction of the data's own magnitude counts as
# zero. Exact-zero tests are not enough: [0.9] * 10 has a float variance near
# 1e-33, which would inflate the statistic to ~1e17.
_DEGENERATE_SCALE_TOLERANCE = 1e-12

# Bootstrap iterations re-run the whole decision, so they use this fraction of
# the configured resample budget. The bootstrap aggregates hundreds of
# iterations into a single support value, which absorbs the extra per-iteration
# Monte Carlo noise that the reduced budget introduces.
BOOTSTRAP_RESAMPLE_DIVISOR = 5


@dataclass(frozen=True)
class PermutationTestResult:
    """Outcome of one adaptive studentized permutation test.

    Attributes:
        statistic: The observed Welch-studentized mean difference, or ``None``
            when a guard short-circuited the test.
        p_greater: Raw one-tailed p-value for "sample x has the larger mean".
        p_less: Raw one-tailed p-value for "sample x has the smaller mean".
        null_skew: Sample skewness of the permutation null of the statistic.
            0 for a symmetric null; the sign gives the direction of the long
            tail. Measured for every test that resampled, whatever it decided.
        p_greater_corrected: ``p_greater`` after correction across the
            one-tailed family.
        p_less_corrected: ``p_less`` after correction across the one-tailed
            family.
        ci_greater: Confidence interval for ``p_greater_corrected``.
        ci_less: Confidence interval for ``p_less_corrected``.
        n_resamples: Total permutations drawn.
        batches: Number of adaptive batches run.
        converged: ``True`` when the interval criterion stopped the run.
        decision: ``greater``, ``less``, ``equivalent``, or ``inconclusive``.
        p_tost: TOST p-value; ``None`` unless the equivalence step ran.
        n_resamples_skew: Permutations behind ``null_skew``.
        note: A short slug naming the guard that short-circuited the test, or
            ``None`` when the test ran.
    """

    statistic: float | None
    p_greater: float | None
    p_less: float | None
    null_skew: float | None
    p_greater_corrected: float | None
    p_less_corrected: float | None
    ci_greater: tuple[float, float] | None
    ci_less: tuple[float, float] | None
    n_resamples: int
    batches: int
    converged: bool
    decision: str
    p_tost: float | None = None
    n_resamples_skew: int = 0
    note: str | None = None


def _guard_result(note):
    """Build an inconclusive result for a test that could not be run.

    Args:
        note: Short slug naming the guard that fired.

    Returns:
        A :class:`PermutationTestResult` with no statistics and an
        ``inconclusive`` decision.
    """
    return PermutationTestResult(
        statistic=None,
        p_greater=None,
        p_less=None,
        null_skew=None,
        p_greater_corrected=None,
        p_less_corrected=None,
        ci_greater=None,
        ci_less=None,
        n_resamples=0,
        batches=0,
        converged=False,
        decision=DECISION_INCONCLUSIVE,
        p_tost=None,
        n_resamples_skew=0,
        note=note,
    )


def studentized_mean_diff(x, y):
    """Compute the Welch-studentized mean difference ``(mx - my) / SE``.

    Args:
        x: First sample of numeric values.
        y: Second sample of numeric values.

    Returns:
        The studentized statistic, or ``nan`` when either sample has fewer than
        two observations or the standard error is not positive.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    nx, ny = x.size, y.size
    if nx < 2 or ny < 2:
        return float("nan")

    standard_error = math.sqrt(
        float(np.var(x, ddof=1)) / nx + float(np.var(y, ddof=1)) / ny
    )
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return float("nan")

    return (float(np.mean(x)) - float(np.mean(y))) / standard_error


def _permutation_statistics(pooled, nx, ny, count, rng):
    """Draw ``count`` random group assignments and studentize each one.

    Samples only the smaller group and recovers the larger by subtracting from
    the pooled totals, using a partial partition instead of a full shuffle.
    ``pooled`` must be mean-centered, which is what keeps that subtraction safe
    in floating point.

    The algebra is derived under "Gate 3 - Adaptive studentized permutation
    test" in the orchestrator guide.

    Args:
        pooled: Mean-centered concatenation of both samples.
        nx: Size of the first group.
        ny: Size of the second group.
        count: Number of permutations to draw.
        rng: A ``numpy.random.Generator``.

    Returns:
        A float array of ``count`` studentized statistics.
    """
    n_total = nx + ny
    # Sample whichever group is smaller; the complement is recovered by
    # subtraction, so the gather cost scales with min(nx, ny) rather than n.
    small_is_x = nx <= ny
    k = nx if small_is_x else ny
    m = n_total - k

    total_sum = float(pooled.sum())
    total_sumsq = float(np.dot(pooled, pooled))

    statistics = np.empty(count, dtype=np.float64)
    chunk_size = max(1, min(count, _MAX_PERMUTATION_CELLS // n_total))

    filled = 0
    while filled < count:
        size = min(chunk_size, count - filled)

        keys = rng.random((size, n_total))
        indices = np.argpartition(keys, k - 1, axis=1)[:, :k]
        values = pooled[indices]

        sum_small = values.sum(axis=1)
        sumsq_small = np.einsum("ij,ij->i", values, values)

        mean_small = sum_small / k
        mean_large = (total_sum - sum_small) / m
        var_small = (sumsq_small - k * mean_small * mean_small) / (k - 1)
        var_large = ((total_sumsq - sumsq_small) - m * mean_large * mean_large) / (m - 1)
        # Rounding can push a zero-variance group microscopically negative.
        np.maximum(var_small, 0.0, out=var_small)
        np.maximum(var_large, 0.0, out=var_large)

        if small_is_x:
            mean_x, var_x, mean_y, var_y = mean_small, var_small, mean_large, var_large
        else:
            mean_x, var_x, mean_y, var_y = mean_large, var_large, mean_small, var_small

        standard_error = np.sqrt(var_x / nx + var_y / ny)
        with np.errstate(divide="ignore", invalid="ignore"):
            chunk_statistics = (mean_x - mean_y) / standard_error
        # A zero standard error means both permuted groups are internally
        # constant. If their means also match, the permutation carries no
        # evidence either way (0/0) and scores 0; if the means differ, there is
        # no scale to divide by and +-inf is the correct extreme value, which
        # numpy already produces.
        np.nan_to_num(chunk_statistics, copy=False, nan=0.0)

        statistics[filled : filled + size] = chunk_statistics
        filled += size

    return statistics


def _scaled_interval(raw_p, corrected_p, interval):
    """Rescale a confidence interval from the raw p-value onto the corrected one.

    Applies the observed correction factor to both bounds, which preserves their
    ordering relative to the threshold because every supported correction is
    monotone non-decreasing.

    Args:
        raw_p: The uncorrected p-value.
        corrected_p: The same p-value after family correction.
        interval: The ``(low, high)`` interval on the raw scale.

    Returns:
        The ``(low, high)`` interval on the corrected scale, clipped to [0, 1].
    """
    low, high = interval
    if raw_p <= 0.0:
        return (float(low), float(high))
    factor = corrected_p / raw_p
    return (min(1.0, float(low) * factor), min(1.0, float(high) * factor))


def _running_skewness(count, sum_t, sum_t2, sum_t3):
    """Compute sample skewness from running power sums.

    Batches are discarded as they are drawn, so the null's shape is accumulated
    as the first three power sums rather than by keeping every statistic.

    Args:
        count: Number of values summed.
        sum_t: Sum of the values.
        sum_t2: Sum of their squares.
        sum_t3: Sum of their cubes.

    Returns:
        The population skewness, or ``None`` when it is undefined (fewer than
        two values, or no spread).
    """
    if count < 2:
        return None
    mean = sum_t / count
    m2 = sum_t2 / count - mean * mean
    if m2 <= 0.0:
        return None
    m3 = sum_t3 / count - 3.0 * mean * sum_t2 / count + 2.0 * mean**3
    return float(m3 / m2**1.5)


def _shifted_null_tail_p(x, y, shift, tail, count, rng):
    """Run a one-sided permutation test against a shifted null.

    Subtracting ``shift`` from ``x`` makes the samples exchangeable under
    ``H0: mean(x) - mean(y) == shift``.

    Args:
        x: First sample.
        y: Second sample.
        shift: The mean difference asserted by the null hypothesis.
        tail: ``greater`` or ``less``, the direction that rejects the null.
        count: Number of permutations to draw.
        rng: A ``numpy.random.Generator``.

    Returns:
        The add-one one-sided p-value.
    """
    shifted = x - shift
    observed = studentized_mean_diff(shifted, y)
    if not math.isfinite(observed):
        return 1.0

    pooled = np.concatenate([shifted, y])
    pooled = pooled - pooled.mean()
    statistics = _permutation_statistics(pooled, x.size, y.size, count, rng)

    if tail == "greater":
        extreme = int(np.count_nonzero(statistics >= observed))
    else:
        extreme = int(np.count_nonzero(statistics <= observed))
    return (1 + extreme) / (count + 1)


def _pooled_standard_deviation(x, y):
    """Compute the pooled standard deviation of two samples.

    Args:
        x: First sample.
        y: Second sample.

    Returns:
        The root-mean-square of the two unbiased sample standard deviations.
    """
    return math.sqrt(
        (float(np.var(x, ddof=1)) + float(np.var(y, ddof=1))) / 2.0
    )


def _equivalence_p_value(x, y, pooled_sd, delta, count, rng):
    """Run the two one-sided tests (TOST) for equivalence of the two means.

    Equivalence needs both nulls rejected, which makes this an
    intersection-union test and is why the pair needs no correction.

    Args:
        x: First sample.
        y: Second sample.
        pooled_sd: The pooled standard deviation, which sets the scale of the
            margin.
        delta: Equivalence margin as an effect size.
        count: Number of permutations to draw per one-sided test.
        rng: A ``numpy.random.Generator``.

    Returns:
        The TOST p-value ``max(p_lower, p_upper)``.
    """
    margin = delta * pooled_sd
    p_lower = _shifted_null_tail_p(x, y, -margin, "greater", count, rng)
    p_upper = _shifted_null_tail_p(x, y, margin, "less", count, rng)
    return max(p_lower, p_upper)


def run_studentized_permutation_test(
    x,
    y,
    *,
    alpha=DEFAULT_ALPHA_PERM,
    min_resamples=DEFAULT_PERMUTATION_MIN_RESAMPLES,
    max_resamples=DEFAULT_PERMUTATION_MAX_RESAMPLES,
    ci_method=DEFAULT_PERMUTATION_CI_METHOD,
    correction=DEFAULT_P_VALUE_CORRECTION,
    equivalence_test=True,
    rng=None,
):
    """Run the adaptive studentized permutation test on two samples.

    Resamples in growing batches until a confidence interval around each
    one-tailed p-value excludes ``alpha``, or the budget is spent. The tail pair
    is corrected against itself; when neither is significant the TOST step at
    :data:`EQUIVALENCE_DELTA` splits ``equivalent`` from ``inconclusive``. The
    skewness of the permutation null is reported alongside, whatever the
    outcome.

    Args:
        x: First sample (concordant tree heights).
        y: Second sample (discordant1 tree heights).
        alpha: Significance threshold, applied to each one-tailed test and to
            the TOST p-value.
        min_resamples: Size of the first batch and the minimum total.
        max_resamples: Resample budget. The run stops once the total reaches it;
            the final batch is drawn whole, so the total may overshoot it by up
            to one batch.
        ci_method: Binomial interval method passed to statsmodels.
        correction: Multiple-testing correction applied across the one-tailed
            family.
        equivalence_test: When ``False``, skip the TOST step and report
            ``inconclusive`` for any non-directional outcome.
        rng: A ``numpy.random.Generator``. A fresh default generator is used
            when omitted.

    Returns:
        A :class:`PermutationTestResult`.
    """
    if rng is None:
        rng = np.random.default_rng()

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx, ny = x.size, y.size
    n_total = nx + ny

    # Guard: a one-observation group has no unbiased variance, so the statistic
    # is undefined. Returning nan here would silently poison the permutation
    # comparisons (every nan comparison is False, which drives both one-tailed
    # counts to zero and makes both directions look maximally significant).
    if nx < 2 or ny < 2:
        return _guard_result("insufficient_group_size")

    mean_x = float(np.mean(x))
    mean_y = float(np.mean(y))
    var_x = float(np.var(x, ddof=1))
    var_y = float(np.var(y, ddof=1))
    standard_error = math.sqrt(var_x / nx + var_y / ny)

    # Both remaining scale guards are relative to the magnitude of the data
    # rather than to exact zero; see _DEGENERATE_SCALE_TOLERANCE.
    data_scale = float(np.max(np.abs(np.concatenate([x, y]))))
    scale_floor = _DEGENERATE_SCALE_TOLERANCE * max(data_scale, 1.0)

    # Guard: every pooled observation is effectively identical, so there is no
    # difference to detect and no scale on which to measure one.
    if standard_error <= scale_floor and abs(mean_x - mean_y) <= scale_floor:
        return _guard_result("zero_pooled_variance")

    # Guard: both groups are internally constant but their means differ. The
    # studentized statistic then divides a real difference by numerical noise,
    # while every permutation that mixes the two constants has genuine spread
    # and stays finite -- so the test would report the p-value floor no matter
    # how few observations back the difference up.
    if not math.isfinite(standard_error) or standard_error <= scale_floor:
        return _guard_result("degenerate_observed_scale")

    # Guard: with very small samples the permutation distribution has fewer
    # distinct values than the requested batch, so its resolution is capped well
    # short of alpha and the interval can never exclude the threshold.
    if n_total <= _EXACT_SUPPORT_CHECK_MAX_N:
        support = math.comb(n_total, min(nx, ny))
        if support < min_resamples:
            return _guard_result("insufficient_permutation_support")

    observed = (mean_x - mean_y) / standard_error
    pooled = np.concatenate([x, y])
    pooled = pooled - pooled.mean()

    max_resamples = max(int(max_resamples), int(min_resamples))
    ci_alpha = 1.0 - _CI_LEVEL

    count_greater = 0
    count_less = 0
    # Running power sums over every permuted statistic, so the null's skewness
    # can be reported without holding the draws from earlier batches.
    sum_t = sum_t2 = sum_t3 = 0.0
    n_done = 0
    batches = 0
    batch = int(min_resamples)
    converged = False

    p_greater = p_less = 1.0
    p_greater_corrected = p_less_corrected = 1.0
    ci_greater = ci_less = (0.0, 1.0)

    while n_done < max_resamples:
        # The final batch is drawn at its full grown size rather than trimmed to
        # the remaining budget. Sampling is vectorized, so a batch costs the same
        # per permutation however large it is, and the extra draws sharpen the
        # interval that decides convergence instead of being spent on a stub that
        # can only narrow it a little. ``max_resamples`` is therefore the point
        # at which the run stops asking for more, not a hard cap on the total.
        # The overshoot is at most one batch. Across a budget wide enough to span
        # several batches that is about a quarter of the total (each batch is
        # 1.25x the previous, so the last one is roughly a quarter of the sum);
        # when ``max_resamples`` sits just above ``min_resamples`` a single grown
        # batch is comparable to the whole budget, and the total can approach
        # 2.25x it. ``_permutation_statistics`` chunks internally against
        # ``_MAX_PERMUTATION_CELLS``, so an oversized batch stays memory-safe.
        batch = max(1, batch)
        statistics = _permutation_statistics(pooled, nx, ny, batch, rng)

        count_greater += int(np.count_nonzero(statistics >= observed))
        count_less += int(np.count_nonzero(statistics <= observed))
        sum_t += float(statistics.sum())
        sum_t2 += float(np.dot(statistics, statistics))
        sum_t3 += float(np.sum(statistics**3))
        n_done += batch
        batches += 1

        # Add-one estimator (Phipson & Smyth 2010): counting the observed
        # arrangement itself keeps the p-value from ever reaching zero and
        # preserves the test's nominal level.
        p_greater = (1 + count_greater) / (n_done + 1)
        p_less = (1 + count_less) / (n_done + 1)

        p_greater_corrected, p_less_corrected = adjust_p_values(
            [p_greater, p_less], method=correction, alpha=alpha
        )

        # The p-value and its interval are two summaries of the *same* pair of
        # numbers, which is what keeps them consistent. The randomness lives
        # entirely in the resampling count, which is Binomial(n_done, p_true)
        # for the exact permutation p-value p_true. ``proportion_confint``
        # takes exactly that count-and-total pair, and its point estimate is
        # ``count / nobs`` -- so passing ``(count + 1, n_done + 1)`` makes the
        # proportion it works from identical to the add-one p-value reported
        # above, rather than an estimator differing from it by ~1/n_done.
        # The interval is mildly conservative because one of those n_done + 1
        # arrangements (the observed one) is fixed rather than drawn.
        raw_ci_greater = proportion_confint(
            count_greater + 1, n_done + 1, alpha=ci_alpha, method=ci_method
        )
        raw_ci_less = proportion_confint(
            count_less + 1, n_done + 1, alpha=ci_alpha, method=ci_method
        )
        ci_greater = _scaled_interval(p_greater, p_greater_corrected, raw_ci_greater)
        ci_less = _scaled_interval(p_less, p_less_corrected, raw_ci_less)

        alpha_inside_greater = ci_greater[0] <= alpha <= ci_greater[1]
        alpha_inside_less = ci_less[0] <= alpha <= ci_less[1]

        if not alpha_inside_greater and not alpha_inside_less:
            converged = True
            break

        batch = int(batch * _BATCH_GROWTH_FACTOR)

    greater_significant = p_greater_corrected <= alpha
    less_significant = p_less_corrected <= alpha

    p_tost = None
    if greater_significant and not less_significant:
        decision = DECISION_GREATER
    elif less_significant and not greater_significant:
        decision = DECISION_LESS
    elif not equivalence_test:
        decision = DECISION_INCONCLUSIVE
    else:
        # The equivalence step reuses the resample count the directional test
        # settled on, so both questions are answered at the same Monte Carlo
        # resolution.
        p_tost = _equivalence_p_value(
            x,
            y,
            _pooled_standard_deviation(x, y),
            EQUIVALENCE_DELTA,
            int(n_done),
            rng,
        )
        decision = (
            DECISION_EQUIVALENT if p_tost <= alpha else DECISION_INCONCLUSIVE
        )

    null_skew = _running_skewness(n_done, sum_t, sum_t2, sum_t3)

    note = None
    if greater_significant and less_significant:
        # Both tails significant at once is not reachable from a coherent
        # permutation distribution, since the two one-tailed counts overlap on
        # ties and therefore sum to more than the total. Reaching it means the
        # inputs or the accumulators are inconsistent, so it is surfaced rather
        # than silently collapsed into a direction.
        note = "both_tails_significant"
    elif not converged:
        note = "max_resamples_reached"

    return PermutationTestResult(
        statistic=float(observed),
        p_greater=float(p_greater),
        p_less=float(p_less),
        null_skew=null_skew,
        p_greater_corrected=float(p_greater_corrected),
        p_less_corrected=float(p_less_corrected),
        ci_greater=ci_greater,
        ci_less=ci_less,
        n_resamples=int(n_done),
        batches=int(batches),
        converged=converged,
        decision=decision,
        p_tost=None if p_tost is None else float(p_tost),
        n_resamples_skew=int(n_done),
        note=note,
    )


def bootstrap_resample_budget(min_resamples, max_resamples):
    """Scale the resample budget down for bootstrap iterations.

    Args:
        min_resamples: The configured minimum resample count.
        max_resamples: The configured maximum resample count.

    Returns:
        A ``(min_resamples, max_resamples)`` tuple reduced by
        :data:`BOOTSTRAP_RESAMPLE_DIVISOR`, with a floor of 1.
    """
    scaled_min = max(1, int(min_resamples) // BOOTSTRAP_RESAMPLE_DIVISOR)
    scaled_max = max(scaled_min, int(max_resamples) // BOOTSTRAP_RESAMPLE_DIVISOR)
    return scaled_min, scaled_max
