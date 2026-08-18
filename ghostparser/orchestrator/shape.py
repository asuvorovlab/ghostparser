"""Distribution-shape diagnostics for the per-topology tree-height groups.

Modality, skewness and tail weight, described under "Shape diagnostics" in the
orchestrator guide.
"""

import math

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import genpareto, kurtosis, skew

__all__ = [
    "SHAPE_FIELD_NAMES",
    "SHAPE_MIN_OBSERVATIONS",
    "SHAPE_MIN_TAIL_EXCEEDANCES",
    "SHAPE_TAIL_QUANTILE",
    "count_kde_modes",
    "critical_bandwidth",
    "describe_shape",
    "silverman_modality_p_value",
    "tail_shape_index",
]

# Per-group fields written for every topology group when the diagnostics run.
SHAPE_FIELD_NAMES = (
    "n_modes",
    "modes_p",
    "skew",
    "excess_kurtosis",
    "tail_xi",
)

# Below this many observations a KDE has nothing to resolve and the moment
# estimates are dominated by sampling noise.
SHAPE_MIN_OBSERVATIONS = 20

# The upper decile is the peaks-over-threshold region; the generalized Pareto
# fit needs a workable number of points inside it.
SHAPE_TAIL_QUANTILE = 0.90
SHAPE_MIN_TAIL_EXCEEDANCES = 10

# Smoothed-bootstrap replicates for the modality test. Its p-value resolution is
# 1 / (n + 1), which this floors well under a conventional alpha. Unrelated to
# the run's gene-tree bootstrap, which never reaches this module.
MODALITY_BOOTSTRAP_RESAMPLES = 200

# Bisection bracket for the critical bandwidth, in units of the sample's own
# standard deviation, and the relative width at which bisection stops.
_BANDWIDTH_LOW_FACTOR = 1e-3
_BANDWIDTH_HIGH_FACTOR = 10.0
_BANDWIDTH_TOLERANCE = 1e-3

# Evaluation grid for mode counting, padded by three bandwidths on each side.
# Doubling it leaves the modality p-values unchanged and triples the cost, since
# the critical bandwidth is far wider than the grid spacing either way.
_GRID_POINTS = 256
_GRID_PADDING = 3.0

# Sample block size when evaluating the density, so the grid-by-sample matrix
# stays small for large groups.
_KERNEL_BLOCK = 2048


def _as_clean_array(values):
    """Coerce a height group to a finite float array.

    Args:
        values: Sequence of numeric values.

    Returns:
        A 1-D float array holding only the finite entries.
    """
    array = np.asarray(values, dtype=float).ravel()
    return array[np.isfinite(array)]


def count_kde_modes(values, bandwidth):
    """Count local maxima of a Gaussian KDE at one absolute bandwidth.

    Args:
        values: Sample to smooth.
        bandwidth: Kernel standard deviation, in the units of ``values``.

    Returns:
        The number of peaks in the smoothed density.
    """
    array = np.asarray(values, dtype=float)
    if bandwidth <= 0.0:
        return 0

    pad = _GRID_PADDING * bandwidth
    grid = np.linspace(array.min() - pad, array.max() + pad, _GRID_POINTS)

    # Normalizing constants are dropped: only the location of the maxima
    # matters, so the unnormalized kernel sum is enough.
    density = np.zeros(_GRID_POINTS, dtype=float)
    for start in range(0, array.size, _KERNEL_BLOCK):
        block = array[start : start + _KERNEL_BLOCK]
        z = (grid[:, None] - block[None, :]) / bandwidth
        density += np.exp(-0.5 * z * z).sum(axis=1)

    peaks, _ = find_peaks(density)
    return int(peaks.size)


def critical_bandwidth(values, max_modes=1):
    """Find the smallest bandwidth whose KDE has at most ``max_modes`` modes.

    Args:
        values: Sample to smooth.
        max_modes: Mode count the bandwidth must bring the density down to.

    Returns:
        The critical bandwidth, in the units of ``values``.
    """
    array = np.asarray(values, dtype=float)
    spread = float(np.std(array))
    if spread <= 0.0:
        return 0.0

    # Mode count is monotone non-increasing in bandwidth for a Gaussian kernel
    # in one dimension (Silverman 1981), which is what makes bisection valid.
    low = _BANDWIDTH_LOW_FACTOR * spread
    high = _BANDWIDTH_HIGH_FACTOR * spread
    if count_kde_modes(array, high) > max_modes:
        return high
    while high - low > _BANDWIDTH_TOLERANCE * spread:
        middle = 0.5 * (low + high)
        if count_kde_modes(array, middle) > max_modes:
            low = middle
        else:
            high = middle
    return high


def _smoothed_bootstrap_samples(array, bandwidth, n_resamples, rng):
    """Draw variance-corrected smoothed-bootstrap resamples.

    Args:
        array: The observed sample.
        bandwidth: Critical bandwidth defining the sampling density.
        n_resamples: Number of resamples to draw.
        rng: A ``numpy.random.Generator``.

    Returns:
        An ``(n_resamples, len(array))`` array of resamples from the
        critical-bandwidth density.
    """
    n = array.size
    mean = float(np.mean(array))
    variance = float(np.var(array))

    drawn = rng.choice(array, size=(n_resamples, n), replace=True)
    noise = rng.standard_normal((n_resamples, n))
    # Rescaling about the mean keeps each resample's variance equal to the
    # data's, so the test compares smoothing and not spread.
    shrink = math.sqrt(1.0 + bandwidth**2 / variance)
    return mean + (drawn + bandwidth * noise - mean) / shrink


def silverman_modality_p_value(
    values,
    max_modes=1,
    n_resamples=MODALITY_BOOTSTRAP_RESAMPLES,
    rng=None,
):
    """Test "the density has at most ``max_modes`` modes" by critical bandwidth.

    Args:
        values: Sample to test.
        max_modes: Mode count under the null.
        n_resamples: Smoothed-bootstrap replicates.
        rng: A ``numpy.random.Generator``; one is created when omitted.

    Returns:
        The add-one p-value, or ``None`` when the sample is too small or flat.
    """
    array = _as_clean_array(values)
    if array.size < SHAPE_MIN_OBSERVATIONS or float(np.var(array)) <= 0.0:
        return None

    generator = np.random.default_rng() if rng is None else rng
    observed = critical_bandwidth(array, max_modes)

    # A resample needs more smoothing than the data exactly when it still has
    # too many modes at the observed bandwidth, so each replicate costs one
    # mode count instead of a second bisection.
    resamples = _smoothed_bootstrap_samples(array, observed, n_resamples, generator)
    exceedances = sum(
        count_kde_modes(resample, observed) > max_modes for resample in resamples
    )

    return (1.0 + exceedances) / (1.0 + n_resamples)


def tail_shape_index(values, quantile=SHAPE_TAIL_QUANTILE):
    """Estimate the generalized-Pareto shape of the upper tail.

    Args:
        values: Sample whose upper tail is fitted.
        quantile: Threshold quantile defining the exceedances.

    Returns:
        The shape parameter, or ``None`` when too few points clear the
        threshold. Positive is a power-law tail, ``0`` exponential, negative a
        bounded one.
    """
    array = _as_clean_array(values)
    if array.size < SHAPE_MIN_OBSERVATIONS:
        return None

    threshold = float(np.quantile(array, quantile))
    exceedances = array[array > threshold] - threshold
    if exceedances.size < SHAPE_MIN_TAIL_EXCEEDANCES:
        return None

    shape, _, _ = genpareto.fit(exceedances, floc=0.0)
    return float(shape) if math.isfinite(shape) else None


def describe_shape(values, rng=None, n_resamples=MODALITY_BOOTSTRAP_RESAMPLES):
    """Measure modality, skewness and tail weight for one height group.

    Args:
        values: The group's tree heights.
        rng: A ``numpy.random.Generator`` for the modality bootstrap.
        n_resamples: Smoothed-bootstrap replicates for the modality test.

    Returns:
        A dict over :data:`SHAPE_FIELD_NAMES`; every value is ``None`` when the
        group is smaller than :data:`SHAPE_MIN_OBSERVATIONS` or has no spread.
    """
    array = _as_clean_array(values)
    empty = dict.fromkeys(SHAPE_FIELD_NAMES)
    if array.size < SHAPE_MIN_OBSERVATIONS or float(np.var(array)) <= 0.0:
        return empty

    # Counted at Scott's bandwidth, the conventional default, so the number is
    # comparable across groups but stays a bandwidth-dependent description.
    scott_bandwidth = float(np.std(array)) * array.size ** (-1.0 / 5.0)

    return {
        "n_modes": count_kde_modes(array, scott_bandwidth),
        "modes_p": silverman_modality_p_value(
            array,
            max_modes=1,
            n_resamples=n_resamples,
            rng=rng,
        ),
        "skew": float(skew(array)),
        "excess_kurtosis": float(kurtosis(array, fisher=True)),
        "tail_xi": tail_shape_index(array),
    }
