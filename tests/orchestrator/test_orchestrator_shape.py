"""Cover the optional distribution-shape diagnostics.

Three kinds of check live here:

- **Parity.** Skewness and excess kurtosis are compared against
  ``scipy.stats.skew``/``kurtosis`` called directly, and the tail index against
  the generalized-Pareto shape the sample was drawn with.
- **Calibration.** The modality test is driven with samples of known modality:
  five unimodal families must not be rejected, and a well-separated mixture
  must be.
- **Contract.** The results TSV carries the shape columns only when the
  diagnostics ran, and the row stays aligned with the header either way.

See ``tests/TEST_IO.md`` for the derivation of each expected value.
"""

import numpy as np
import pytest
from scipy import stats

from ghostparser.orchestrator import inference as pinf
from ghostparser.orchestrator import shape as pshape

_TRIPLET = ("A", "B", "C")
_SPECIES_SUBTREE = "((A:1.0,B:1.0):1.0,C:2.0);"

_CON = "((A,B),C)"
_DIS1 = "((B,C),A)"
_DIS2 = "((A,C),B)"

_ALPHA = 0.05


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
        observations.extend((topology, float(height), None) for height in heights)
    return observations


def test_shape_moments_match_scipy():
    """Skewness and excess kurtosis reproduce the SciPy reference exactly."""
    values = np.random.default_rng(4).lognormal(0.0, 0.7, 500)

    described = pshape.describe_shape(values, rng=np.random.default_rng(1))

    assert described["skew"] == pytest.approx(float(stats.skew(values)))
    assert described["excess_kurtosis"] == pytest.approx(
        float(stats.kurtosis(values, fisher=True))
    )
    # A lognormal is right-skewed with heavy shoulders; the signs are the point.
    assert described["skew"] > 0
    assert described["excess_kurtosis"] > 0


@pytest.mark.parametrize(
    "name, multimodal",
    [
        ("normal", False),
        ("lognormal", False),
        ("exponential", False),
        ("gamma", False),
        ("mixture", True),
    ],
)
def test_modality_test_rejects_only_a_well_separated_mixture(name, multimodal):
    """Silverman's test holds its level on unimodal shapes and detects a mixture.

    The four unimodal families include strongly skewed ones, where a raw KDE
    peak count at a default bandwidth reports spurious modes. The mixture's two
    components sit four standard deviations apart, which is inside the test's
    resolution.
    """
    rng = np.random.default_rng(17)
    samples = {
        "normal": rng.normal(0.0, 1.0, 400),
        "lognormal": rng.lognormal(0.0, 0.6, 400),
        "exponential": rng.exponential(1.0, 400),
        "gamma": rng.gamma(2.0, 1.0, 400),
        "mixture": np.concatenate(
            [rng.normal(0.0, 1.0, 200), rng.normal(4.0, 1.0, 200)]
        ),
    }

    p_value = pshape.silverman_modality_p_value(
        samples[name], rng=np.random.default_rng(2)
    )

    assert (p_value <= _ALPHA) is multimodal


@pytest.mark.parametrize(
    "name, expected_xi",
    [
        ("pareto", 1.0 / 3.0),
        ("exponential", 0.0),
        ("uniform", -1.0),
    ],
)
def test_tail_index_recovers_known_tail_shapes(name, expected_xi):
    """The generalized-Pareto shape recovers the tail each sample was drawn with.

    A Pareto with index ``a`` has ``xi = 1 / a``, an exponential tail has
    ``xi = 0``, and a uniform has a bounded ``xi = -1``. The peaks-over-threshold
    fit sees only the upper decile, so the tolerance reflects the few hundred
    points it has left.
    """
    rng = np.random.default_rng(23)
    samples = {
        "pareto": rng.pareto(3.0, 4000) + 1.0,
        "exponential": rng.exponential(1.0, 4000),
        "uniform": rng.uniform(0.0, 1.0, 4000),
    }

    xi = pshape.tail_shape_index(samples[name])

    assert xi == pytest.approx(expected_xi, abs=0.25)


def test_shape_is_not_described_for_small_or_flat_groups():
    """Groups too small or with no spread report every field empty."""
    assert all(value is None for value in pshape.describe_shape([1.0] * 5).values())
    assert all(value is None for value in pshape.describe_shape([2.0] * 50).values())

    # A group large enough to describe but with too thin an upper decile for
    # the tail fit keeps its moments and drops only the tail index.
    described = pshape.describe_shape(
        np.random.default_rng(5).normal(0.0, 1.0, 25), rng=np.random.default_rng(1)
    )
    assert described["skew"] is not None
    assert described["tail_xi"] is None


def _shape_result(seed=9, enabled=True):
    """Analyze one lognormal triplet with the diagnostics on or off.

    Args:
        seed: Seed for the drawn heights.
        enabled: Whether to measure the shape diagnostics.

    Returns:
        The decided ``TripletPipelineResult``.
    """
    rng = np.random.default_rng(seed)
    result = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        _observations(
            rng.lognormal(0.0, 0.4, 60),
            rng.lognormal(-0.5, 0.4, 40),
            rng.lognormal(-0.5, 0.4, 30),
        ),
        species_subtree=_SPECIES_SUBTREE,
        triplet_seed=3,
        bootstrap_options={"iterations": 5},
        shape_diagnostics=enabled,
    )
    # The writers take decided results, as they do from the runner.
    return pinf._apply_triplet_result_p_value_correction(
        [result], alpha_dct=0.05, alpha_ks=0.05, method="no"
    )[0]


def test_shape_is_measured_once_and_not_per_bootstrap_iteration():
    """The diagnostics come from the point estimate, so iterations cannot move them.

    They describe the observed height groups; a bootstrap resample is a
    different sample, and measuring one would cost the modality bootstrap on
    every iteration.
    """
    few = _shape_result()
    many = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        _observations(
            np.random.default_rng(9).lognormal(0.0, 0.4, 60),
            np.random.default_rng(9).lognormal(-0.5, 0.4, 40),
            np.random.default_rng(9).lognormal(-0.5, 0.4, 30),
        ),
        species_subtree=_SPECIES_SUBTREE,
        triplet_seed=3,
        bootstrap_options={"iterations": 60},
        shape_diagnostics=True,
    )

    assert few.shape_statistics["con_skew"] == many.shape_statistics["con_skew"]
    assert few.shape_statistics["con_modes_p"] == many.shape_statistics["con_modes_p"]


@pytest.mark.output
@pytest.mark.parametrize("enabled", [False, True])
def test_summary_statistics_tsv_never_carries_shape_columns(enabled, tmp_path):
    """The summary TSV holds no shape diagnostics, measured or not.

    That file is a feature matrix, and the diagnostics are undefined for groups
    below their observation floors — a group of 12 has no modality p-value and a
    thin upper tail has no tail index. Carrying them there would punch holes in
    every row that hit one, so they live in the results TSV alone. Asserted with
    the diagnostics both off and on, since the "on" case is the one that would
    otherwise leak columns.
    """
    result = _shape_result(enabled=enabled)
    path = tmp_path / "summary.tsv"
    pinf.write_summary_statistics_tsv([result], str(path))

    lines = path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")

    assert not any(column.endswith(f"_{field}") for column in header
                   for field in pshape.SHAPE_FIELD_NAMES)
    assert len(row) == len(header)
    # The descriptive per-topology columns the file does carry are unaffected.
    assert "concordant_avg_tree_height_mean" in header
    assert "classification" in header
    # Even with the diagnostics measured, nothing of theirs reaches this file.
    if enabled:
        assert result.shape_statistics is not None


@pytest.mark.output
@pytest.mark.parametrize("enabled", [False, True])
def test_results_tsv_carries_shape_columns_only_when_enabled(enabled, tmp_path):
    """The fifteen shape columns follow the setting and never shift the row.

    They are read off the results rather than passed to the writer, so a run
    that did not measure them cannot emit empty columns claiming it did.
    """
    result = _shape_result(enabled=enabled)
    path = tmp_path / "results.tsv"
    pinf.write_pipeline_results([result], str(path), p_value_correction="no")

    lines = path.read_text().strip().splitlines()
    header = lines[0].split("\t")
    expected = [
        f"{label}_{field}"
        for label in pinf.SHAPE_GROUP_LABELS
        for field in pshape.SHAPE_FIELD_NAMES
    ]
    assert len(expected) == 15
    assert all((column in header) is enabled for column in expected)
    assert len(lines[1].split("\t")) == len(header)

    if enabled:
        row = dict(zip(header, lines[1].split("\t")))
        # 60 concordant heights clear the minimum; 30 discordant2 heights clear
        # it for the moments but leave only 3 in the upper decile.
        assert float(row["con_skew"]) > 0
        assert int(row["con_n_modes"]) >= 1
        assert 0.0 < float(row["con_modes_p"]) <= 1.0
        assert row["dis2_tail_xi"] == ""
    else:
        assert result.shape_statistics is None
