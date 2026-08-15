"""Adaptive studentized permutation test for concordant vs. discordant1 heights.

GhostParser's final inference step asks a directional question: is the mean
tree height of the concordant topology greater than, less than, or
indistinguishable from that of the more frequent discordant topology? This
module answers it with a two-sample permutation test on the Welch-studentized
mean difference, resampling adaptively until the decision is resolved to within
Monte Carlo error.

The studentized statistic is what makes the test valid here. The concordant
sample is usually far larger than the discordant1 sample and the two have
different variances; under that combination a permutation test of the raw mean
difference does not hold its nominal level, while the studentized version
remains asymptotically valid (Janssen 1997, *Statistics & Probability Letters*
36(1), 9-21, https://doi.org/10.1016/S0167-7152(97)00043-6).

p-values use the add-one estimator ``(1 + count) / (1 + resamples)``, which is
the unbiased and correctly-sized estimator for a Monte Carlo permutation
p-value; the naive ``count / resamples`` ratio can report zero and understates
the true type-I error rate (Phipson & Smyth 2010, *Statistical Applications in
Genetics and Molecular Biology* 9(1), Article 39,
https://doi.org/10.2202/1544-6115.1585).

See :doc:`ORCHESTRATOR.md <ORCHESTRATOR>` for the full statistical write-up,
including the sampling optimization used by :func:`_permutation_statistics`.
"""

from __future__ import annotations

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
    "PermutationTestResult",
    "run_studentized_permutation_test",
    "median_sign_decision",
    "bootstrap_resample_budget",
]

# Decision labels. ``greater``/``less`` refer to the first sample (concordant)
# relative to the second (discordant1); ``ambiguous`` means the data do not
# support a directional call.
DECISION_GREATER = "greater"
DECISION_LESS = "less"
DECISION_AMBIGUOUS = "ambiguous"

# Each adaptive batch after the first is this multiple of the previous one, so
# a run that is still undecided spends geometrically more effort per check
# instead of paying the confidence-interval overhead on every few resamples.
_BATCH_GROWTH_FACTOR = 1.25

# Confidence level of the interval placed around each p-value by the stopping
# rule. Fixed rather than configurable: it governs how sure the stopping rule
# must be before it commits, which is an internal precision knob rather than a
# statistical choice the analysis depends on. 0.95 matches the statsmodels
# default (``proportion_confint(alpha=0.05)``).
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
# zero. Testing against exact zero is not enough: a sample of nominally
# identical values such as [0.9] * 10 has a floating-point variance around
# 1e-33 rather than 0, which would divide a real mean difference by almost
# nothing and report a studentized statistic of ~1e17 as though it were
# overwhelming evidence. The threshold sits far above double-precision noise
# (~1e-16 relative) and far below any real difference in tree heights.
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
        p_two_sided: Raw two-tailed p-value counting permutations at least as
            extreme in absolute value. Reported as a diagnostic only: it assumes
            a symmetric permutation null, so it disagrees with the directional
            p-values when the null is skewed (see ``null_skewed``).
        p_greater_corrected: ``p_greater`` after correction across the
            one-tailed family.
        p_less_corrected: ``p_less`` after correction across the one-tailed
            family.
        ci_greater: Confidence interval for ``p_greater_corrected``.
        ci_less: Confidence interval for ``p_less_corrected``.
        n_resamples: Total permutations drawn.
        batches: Number of adaptive batches run.
        converged: ``True`` when the confidence-interval criterion stopped the
            run, ``False`` when the resample budget was exhausted or a guard
            fired.
        decision: The directional decision driving classification.
        decision_two_sided: The decision from the two-tailed-gate-then-sign
            rule, kept as a cross-check.
        consistent: Whether the two decision rules agree.
        null_skewed: Whether the absolute-value two-tailed p-value contradicts
            the directional decision, which happens when the permutation null is
            asymmetric. Informational: tree heights are bounded below by zero
            and routinely right-skewed, so this is common rather than
            exceptional, and it never overturns the directional call.
        note: A short slug naming the guard that short-circuited the test, or
            ``None`` when the test ran.
    """

    statistic: float | None
    p_greater: float | None
    p_less: float | None
    p_two_sided: float | None
    p_greater_corrected: float | None
    p_less_corrected: float | None
    ci_greater: tuple[float, float] | None
    ci_less: tuple[float, float] | None
    n_resamples: int
    batches: int
    converged: bool
    decision: str
    decision_two_sided: str
    consistent: bool
    null_skewed: bool = False
    note: str | None = None


def _guard_result(note):
    """Build an ambiguous result for a test that could not be run.

    Args:
        note: Short slug naming the guard that fired.

    Returns:
        A :class:`PermutationTestResult` with no statistics and an ambiguous
        decision.
    """
    return PermutationTestResult(
        statistic=None,
        p_greater=None,
        p_less=None,
        p_two_sided=None,
        p_greater_corrected=None,
        p_less_corrected=None,
        ci_greater=None,
        ci_less=None,
        n_resamples=0,
        batches=0,
        converged=False,
        decision=DECISION_AMBIGUOUS,
        decision_two_sided=DECISION_AMBIGUOUS,
        consistent=True,
        null_skewed=False,
        note=note,
    )


def _studentized_mean_diff(x, y):
    """Compute the Welch-studentized mean difference between two samples.

    The statistic is::

        T = (mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)

    where ``var`` is the unbiased sample variance (``ddof=1``). Dividing by the
    Welch standard error rather than a pooled one is what keeps the permutation
    test valid when the two samples differ in both size and variance, which is
    the normal case for concordant vs. discordant1 tree heights.

    Args:
        x: First sample of numeric values.
        y: Second sample of numeric values.

    Returns:
        The studentized statistic as a float, or ``nan`` when either sample has
        fewer than two observations or the standard error is not positive. The
        callers guard these cases before resampling; the ``nan`` is a
        belt-and-braces return, never a value the decision logic consumes.
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

    Two optimizations keep this cheap enough to run inside every bootstrap
    iteration:

    1. **Sample only the smaller group; derive the larger one exactly.** A
       permutation partitions the pooled values into two groups, so the larger
       group is precisely the complement of the smaller one -- nothing about it
       is unknown or estimated. Both of its sufficient statistics follow by
       subtraction. Writing ``S`` and ``Q`` for the pooled totals (computed once,
       outside the loop) and ``s`` and ``q`` for the drawn group's::

           S = sum(pooled)            Q = sum(pooled^2)
           s = sum(drawn)             q = sum(drawn^2)

       the complement of size ``m = n - k`` has::

           sum      = S - s
           sum_sq   = Q - q
           mean     = (S - s) / m
           var      = ((Q - q) - m * mean^2) / (m - 1)

       That last line is the standard ``E[X^2] - E[X]^2`` identity in unbiased
       (``ddof=1``) form, so the recovered variance is algebraically identical to
       calling ``np.var(complement, ddof=1)`` on the values themselves -- it is
       an exact rearrangement, not an approximation. The drawn group's mean and
       variance come from the same two formulas applied to ``s`` and ``q``
       directly. Both groups therefore have a mean and a variance, which is all
       the Welch statistic needs, and the larger group's values are never
       gathered. Since the discordant1 sample is normally much smaller than the
       concordant one, this cuts the gathered data by the size ratio.
       :func:`test_permutation_statistics_match_exhaustive_enumeration` pins this
       against every one of the 126 group assignments of a 4-vs-5 case computed
       the direct way.
    2. **Partial partition instead of a full shuffle.** Taking the ``k``
       smallest of ``n`` uniform random keys yields a uniformly random size-``k``
       subset, and ``np.argpartition`` finds them in one linear pass rather than
       sorting or Fisher-Yates shuffling the whole row.

    ``pooled`` must already be mean-centered. Centering leaves both the mean
    difference and both variances unchanged, but it is what makes the
    subtraction above safe in floating point. ``Q - q`` is a difference of two
    positive sums of squares, and ``m * mean^2`` is near zero once the pooled
    mean is zero, so neither step cancels significant digits. Without centering,
    tree heights of similar magnitude make ``sum_sq - m * mean^2`` a difference
    of two nearly equal large numbers, which destroys most of the precision.

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

    Every supported correction is monotone non-decreasing in each member of the
    family, so applying the observed correction factor to the interval bounds
    preserves their ordering relative to the threshold. This is exact for
    Bonferroni (a constant factor) and holds for the step-up methods as long as
    the family's rank order is stable across the interval, which it is whenever
    the interval is narrow enough for the stopping rule to fire.

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


def run_studentized_permutation_test(
    x,
    y,
    *,
    alpha=DEFAULT_ALPHA_PERM,
    min_resamples=DEFAULT_PERMUTATION_MIN_RESAMPLES,
    max_resamples=DEFAULT_PERMUTATION_MAX_RESAMPLES,
    ci_method=DEFAULT_PERMUTATION_CI_METHOD,
    correction=DEFAULT_P_VALUE_CORRECTION,
    rng=None,
):
    """Run the adaptive studentized permutation test on two samples.

    The first batch draws ``min_resamples`` permutations. After each batch the
    confidence intervals around both one-tailed p-values are checked against
    ``alpha``: once ``alpha`` lies outside both, the decision cannot change with
    more resampling and the run stops. Otherwise the batch size grows by 25% and
    the run continues until ``max_resamples`` is reached.

    The two one-tailed p-values form a testing family and are corrected against
    each other with ``correction`` before being compared to ``alpha``. The
    two-tailed p-value is reported raw as a cross-check and is not corrected,
    since it is a single test rather than a family.

    Args:
        x: First sample (concordant tree heights).
        y: Second sample (discordant1 tree heights).
        alpha: Significance threshold, applied to each one-tailed test and to
            the two-tailed cross-check.
        min_resamples: Size of the first batch and the minimum total.
        max_resamples: Hard ceiling on total permutations.
        ci_method: Binomial interval method passed to statsmodels.
        correction: Multiple-testing correction applied across the one-tailed
            family.
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
    count_two_sided = 0
    n_done = 0
    batches = 0
    batch = int(min_resamples)
    converged = False

    p_greater = p_less = p_two_sided = 1.0
    p_greater_corrected = p_less_corrected = 1.0
    ci_greater = ci_less = (0.0, 1.0)

    while n_done < max_resamples:
        batch = max(1, min(batch, max_resamples - n_done))
        statistics = _permutation_statistics(pooled, nx, ny, batch, rng)

        count_greater += int(np.count_nonzero(statistics >= observed))
        count_less += int(np.count_nonzero(statistics <= observed))
        count_two_sided += int(np.count_nonzero(np.abs(statistics) >= abs(observed)))
        n_done += batch
        batches += 1

        # Add-one estimator (Phipson & Smyth 2010): counting the observed
        # arrangement itself keeps the p-value from ever reaching zero and
        # preserves the test's nominal level.
        p_greater = (1 + count_greater) / (n_done + 1)
        p_less = (1 + count_less) / (n_done + 1)
        p_two_sided = (1 + count_two_sided) / (n_done + 1)

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
    if greater_significant and not less_significant:
        decision = DECISION_GREATER
    elif less_significant and not greater_significant:
        decision = DECISION_LESS
    else:
        decision = DECISION_AMBIGUOUS

    # Cross-check rule: gate on a two-tailed p-value, then take the sign. The
    # gate uses the doubled smaller tail rather than the absolute-value count,
    # because the latter is only valid when the permutation null is symmetric.
    # A tiny discordant1 group carrying a few extreme tree heights produces a
    # strongly right-skewed null, where the absolute-value count borrows the fat
    # right tail to judge a left-tail observation and hides a real difference.
    p_two_sided_gate = min(1.0, 2.0 * min(p_greater, p_less))
    if p_two_sided_gate <= alpha and mean_x > mean_y:
        decision_two_sided = DECISION_GREATER
    elif p_two_sided_gate <= alpha and mean_x < mean_y:
        decision_two_sided = DECISION_LESS
    else:
        decision_two_sided = DECISION_AMBIGUOUS

    # The directional tails resolving a comparison that the absolute-value count
    # leaves open is the signature of an asymmetric null. It is recorded on its
    # own field rather than in ``note`` because it is the common case for tree
    # heights, and ``note`` is reserved for conditions that need attention.
    null_skewed = decision != DECISION_AMBIGUOUS and p_two_sided > alpha

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
        p_two_sided=float(p_two_sided),
        p_greater_corrected=float(p_greater_corrected),
        p_less_corrected=float(p_less_corrected),
        ci_greater=ci_greater,
        ci_less=ci_less,
        n_resamples=int(n_done),
        batches=int(batches),
        converged=converged,
        decision=decision,
        decision_two_sided=decision_two_sided,
        consistent=decision == decision_two_sided,
        null_skewed=null_skewed,
        note=note,
    )


def median_sign_decision(con_heights, dis1_heights):
    """Decide direction by a sign test on the two sample medians.

    PROVISIONAL. This is the pre-permutation-test behaviour, retained only so
    the two inference paths can be run against each other on data with known
    ground truth, and reached solely through ``permutation_test: false``. It
    attaches no p-value and no notion of significance, so any numerical
    difference between the medians -- however small -- yields a confident
    direction. See the removal checklist at the fallback branch in
    ``inference._decide_direction``.

    Args:
        con_heights: Concordant tree heights.
        dis1_heights: Discordant1 tree heights.

    Returns:
        ``greater``, ``less``, or ``ambiguous``.
    """
    if not len(con_heights) or not len(dis1_heights):
        return DECISION_AMBIGUOUS

    median_con = float(np.median(np.asarray(con_heights, dtype=float)))
    median_dis = float(np.median(np.asarray(dis1_heights, dtype=float)))
    if median_con > median_dis:
        return DECISION_GREATER
    if median_con < median_dis:
        return DECISION_LESS
    return DECISION_AMBIGUOUS


def bootstrap_resample_budget(min_resamples, max_resamples):
    """Scale the resample budget down for use inside bootstrap iterations.

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
