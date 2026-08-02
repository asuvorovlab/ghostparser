"""Assert pipeline.inference.analyze_triplet reproduces triplet_processor.analyze_triplet_entry on identical inputs (the ported inference must stay bit-for-bit equal)."""

import pytest

from ghostparser import triplet_processor as tp
from ghostparser.pipeline import inference as pinf

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

_SEED = 20240724
_ITERATIONS = 40

_COMPARED_FIELDS = (
    "triplet",
    "species_tree",
    "most_frequent_matches_concordant",
    "n_con",
    "n_dis1",
    "n_dis2",
    "dis1_topology",
    "dct_statistic",
    "dct_p_value",
    "dct_p_value_corrected",
    "dct_significant",
    "ks_statistic",
    "ks_p_value",
    "ks_p_value_corrected",
    "ks_significant",
    "summary_con",
    "summary_dis",
    "classification",
    "analyzed_trees",
)
# bootstrap_value / all_bootstrap are intentionally excluded: the pipeline
# resamples with NumPy's RNG, so its bootstrap aggregates are statistically
# equivalent to triplet_processor's random.Random ones but not identical.


def _assert_results_equal(pipeline_result, reference_result):
    """Assert two results agree on every non-bootstrap compared field.

    Args:
        pipeline_result: The pipeline's ``TripletPipelineResult``.
        reference_result: The reference ``TripletPipelineResult``.
    """
    for field in _COMPARED_FIELDS:
        pipeline_value = getattr(pipeline_result, field)
        reference_value = getattr(reference_result, field)
        if isinstance(reference_value, float):
            assert pipeline_value == pytest.approx(reference_value), field
        else:
            assert pipeline_value == reference_value, field


def _assert_valid_bootstrap(result):
    """Assert a result's bootstrap aggregates are well formed.

    Args:
        result: A ``TripletPipelineResult``.
    """
    assert 0.0 <= result.bootstrap_value <= 1.0
    assert result.all_bootstrap is not None
    assert sum(result.all_bootstrap.values()) == pytest.approx(1.0)


@pytest.mark.parametrize("discordant_test", ["chi-square", "z-test"])
@pytest.mark.parametrize("summary_statistic", ["mean", "median", "mode"])
@pytest.mark.parametrize("strategy", ["AVG", "A", "B", "C", "SIS", "INT"])
def test_analyze_triplet_matches_analyze_triplet_entry(
    discordant_test, summary_statistic, strategy
):
    """analyze_triplet equals analyze_triplet_entry across every parameter combination.

    The pipeline uses only the scipy/statsmodels backend, so the orchestrator
    reference is pinned to ``stats_backend="standard"``.
    """
    pipeline_result = pinf.analyze_triplet(
        _TRIPLET,
        _GENE_SUBTREES,
        species_subtree=_SPECIES_SUBTREE,
        alpha_dct=0.05,
        alpha_ks=0.05,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        tree_height_calculation_strategy=strategy,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )

    reference_result = tp.analyze_triplet_entry(
        _TRIPLET,
        {"species_tree": _SPECIES_SUBTREE, "gene_trees": _GENE_SUBTREES},
        alpha_dct=0.05,
        alpha_ks=0.05,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend="standard",
        tree_height_calculation_strategy=strategy,
        generate_summary_stats=False,
        bootstrap=True,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )

    _assert_results_equal(pipeline_result, reference_result)
    _assert_valid_bootstrap(pipeline_result)


def test_analyze_triplet_from_observations_matches_newick_path():
    """analyze_triplet_from_observations equals analyze_triplet on equivalent inputs."""
    observations = pinf._serialize_triplet_gene_trees(
        _TRIPLET, _GENE_SUBTREES, tree_height_calculation_strategy="AVG"
    )
    from_obs = pinf.analyze_triplet_from_observations(
        _TRIPLET,
        observations,
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )
    from_newick = pinf.analyze_triplet(
        _TRIPLET,
        _GENE_SUBTREES,
        species_subtree=_SPECIES_SUBTREE,
        tree_height_calculation_strategy="AVG",
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )
    _assert_results_equal(from_obs, from_newick)
    # Same observations + same seed -> identical NumPy bootstrap.
    assert from_obs.bootstrap_value == from_newick.bootstrap_value
    assert from_obs.all_bootstrap == from_newick.all_bootstrap


def test_bootstrap_is_deterministic_under_seed():
    """A fixed seed yields identical bootstrap aggregates across runs."""
    kwargs = dict(
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )
    first = pinf.analyze_triplet(_TRIPLET, _GENE_SUBTREES, **kwargs)
    second = pinf.analyze_triplet(_TRIPLET, _GENE_SUBTREES, **kwargs)
    assert first.all_bootstrap == second.all_bootstrap
    assert first.bootstrap_value == second.bootstrap_value


def test_analyze_triplet_empty_observations():
    """analyze_triplet on zero gene subtrees returns a no-introgression result."""
    result = pinf.analyze_triplet(
        _TRIPLET,
        [],
        species_subtree=_SPECIES_SUBTREE,
        bootstrap_options={"iterations": _ITERATIONS},
        triplet_seed=_SEED,
    )
    assert result.analyzed_trees == 0
    assert result.n_con == 0
    assert result.classification == "no_introgression"
