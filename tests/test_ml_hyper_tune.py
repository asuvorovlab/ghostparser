import argparse
import json
from pathlib import Path

import pytest

from ghostparser.config import ConfigError
from ghostparser.ml import hyper_tune, tuning_report
from ghostparser.cli_config import resolve_cli_or_config_args
from ghostparser.ml.hyper_tune import (
    load_hyper_tune_config,
    normalize_hyper_tune_payload,
    tune_hyperparameters,
)


_VALID_TUNING = {"model": "random_forest", "search_space": {"n_estimators": [5]}}


def _write_config(tmp_path, name, tuning, **top_level):
    payload = {
        "input_path": "./results/summary_statistics.tsv",
        "output_dir": "./results/hyper_tune_out",
        **top_level,
    }
    if tuning is not None:
        payload["hyperparameter_tuning"] = tuning
    config_path = tmp_path / name
    config_path.write_text(json.dumps(payload))
    return config_path


def _tuning_namespace(input_path, output_dir, **overrides):
    config = argparse.Namespace(
        input_path=str(input_path),
        output_dir=str(output_dir),
        target_column="class",
        test_size=0.25,
        cv_folds=2,
        rare_class_policy="warn_reduce_cv",
        seed=7,
        n_jobs=1,
        overwrite=True,
        model_name="random_forest",
        search_method="grid",
        objective_metric="exact_match_accuracy",
        objective_direction="max",
        top_k=3,
        n_iter=2,
        max_candidates=10,
        use_wandb=False,
        search_space={
            "n_estimators": [5, 10],
            "max_depth": [None],
            "min_samples_split": [2],
            "min_samples_leaf": [1],
            "max_features": ["sqrt"],
            "class_weight": [None],
        },
    )
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


@pytest.mark.config
@pytest.mark.parametrize(
    "max_features", [["sqrt", None, 4], "log2"], ids=["list", "lone_value"]
)
def test_load_hyper_tune_config_accepts_hyperparameter_tuning_section(tmp_path, max_features):
    """A full tuning section resolves, with its keys renamed to their config names.

    The section's `model`/`method`/`objective` keys surface as `model_name`,
    `search_method` and `objective_metric`, so this pins the rename as well as
    the acceptance, a caller reading the resolved config uses the latter
    names. Valid search-space candidates pass the `model`-block rules and keep
    their shape, list or lone value alike.
    """
    config_path = _write_config(
        tmp_path,
        "hyper_tune.json",
        {
            "model": "random_forest",
            "method": "grid",
            "objective": "exact_match_accuracy",
            "top_k": 5,
            "max_candidates": 100,
            "use_wandb": False,
            "search_space": {
                "n_estimators": [5, 10],
                "max_depth": [None],
                "min_samples_split": [2],
                "min_samples_leaf": [1],
                "max_features": max_features,
                "class_weight": [None],
            },
        },
    )

    config = load_hyper_tune_config(str(config_path))

    assert config["model_name"] == "random_forest"
    assert config["search_method"] == "grid"
    assert config["objective_metric"] == "exact_match_accuracy"
    assert config["search_space"]["n_estimators"] == [5, 10]
    assert config["search_space"]["max_features"] == max_features
    assert config["use_wandb"] is False


@pytest.mark.config
def test_tuner_cli_flags_override_the_config_file(tmp_path):
    """`--seed` and `--no-overwrite` beside the tuner's config file replace its values."""
    config_path = _write_config(
        tmp_path, "tune.json", _VALID_TUNING, seed=7, overwrite=True
    )
    args = argparse.Namespace(config_file=str(config_path), seed=9, no_overwrite=True)

    config = resolve_cli_or_config_args(
        args,
        normalize_payload=normalize_hyper_tune_payload,
        payload_arg_names=["seed", "no_overwrite"],
    )

    assert config.seed == 9
    assert config.overwrite is False
    assert config.search_space["n_estimators"] == [5]


@pytest.mark.config
@pytest.mark.parametrize(
    "tuning, top_level, match",
    [
        pytest.param(None, {}, "hyperparameter_tuning", id="missing_section"),
        pytest.param(
            {**_VALID_TUNING, "use_wandb": "yes"},
            {},
            "must be a boolean",
            id="non_boolean_use_wandb",
        ),
        pytest.param(
            {**_VALID_TUNING, "use_wandb": False, "wandb_detailed_payloads": True},
            {},
            "wandb_detailed_payloads requires",
            id="detailed_payloads_without_wandb",
        ),
        pytest.param(
            _VALID_TUNING,
            {"evaluation": {"metrics": "all"}},
            "evaluation",
            id="evaluation_section",
        ),
        pytest.param(
            {**_VALID_TUNING, "search_space": {"max_features": ["sqrt", "auto"]}},
            {},
            r"search_space\.max_features.*'auto'",
            id="invalid_search_space_value",
        ),
    ],
)
def test_load_hyper_tune_config_rejects_malformed_configs(
    tuning, top_level, match, tmp_path
):
    """Each structural rule of the tuning config fails by name.

    One case per rule: the section itself is required (nothing about a search
    can be defaulted, there is no search space to infer); `use_wandb` may be
    omitted but not mistyped, and it gates the detailed-payload flag rather than
    letting it be silently ignored; an `evaluation` block is refused rather
    than ignored, since accepting it would let a user believe their evaluation
    settings applied to every candidate; and a bad search-space candidate is
    named against its `search_space` key before any fit.
    """
    config_path = _write_config(tmp_path, "hyper_tune_bad.json", tuning, **top_level)

    with pytest.raises(ConfigError, match=match):
        load_hyper_tune_config(str(config_path))


@pytest.mark.config
@pytest.mark.parametrize(
    "sample_name, expected_model",
    [
        ("hyperparameter_tuning_random_forest.yaml", "random_forest"),
        ("hyperparameter_tuning_multi_knn.json", "multi_knn"),
    ],
)


def test_shipped_tuning_sample_configs_resolve(sample_name, expected_model):
    """The shipped tuner samples load and resolve to the model each one names.

    Guards against a sample drifting out of step with the validator, which would
    leave users copying a config the loader rejects, in both config formats.
    """
    sample_dir = Path(__file__).resolve().parents[1] / "sample_configs"
    config = load_hyper_tune_config(str(sample_dir / sample_name))

    assert config["model_name"] == expected_model
    assert config["use_wandb"] is False
    assert config["wandb_detailed_payloads"] is False
    assert set(config["search_space"]) <= set(
        hyper_tune.MODEL_SEARCH_KEYS[expected_model]
    )


@pytest.mark.integration
@pytest.mark.output
def test_tune_hyperparameters_grid_search_runs_end_to_end(
    summary_statistics_tsv_tuning, tmp_path
):
    """Grid search enumerates the whole space and writes every local artifact.

    The search space crosses one parameter over two values, so grid search
    must evaluate exactly 2 candidates: the count is what distinguishes
    exhaustive enumeration from sampling. The namespace carries no
    ``use_wandb``: the programmatic entry point takes the same default as the
    config loader, so the run stays local and still emits the ranked TSVs,
    the marginals, the plot and the report.
    """
    output_dir = tmp_path / "hyper_tune_out"
    config = _tuning_namespace(summary_statistics_tsv_tuning, output_dir)
    del config.use_wandb

    result = tune_hyperparameters(config)

    assert result["results"]["model_name"] == "random_forest"
    assert result["results"]["search_method"] == "grid"
    assert result["results"]["use_wandb"] is False
    assert len(result["candidates"]) == 2
    assert (output_dir / "hyper_tune_best_model.pkl").exists()
    assert not (output_dir / "wandb").exists()

    marginals_path = output_dir / "hyper_tune_parameter_marginals.tsv"
    plot_path = output_dir / "figures" / "hyper_tune_search_report.png"
    assert marginals_path.exists()
    assert plot_path.exists()
    assert result["plot_path"] == str(plot_path)

    candidate_header = (
        (output_dir / "hyper_tune_results.tsv").read_text().splitlines()[0].split("\t")
    )
    assert {"rank", "is_best", "cv_score", "elapsed_seconds"} <= set(candidate_header)

    marginal_header = marginals_path.read_text().splitlines()[0].split("\t")
    assert marginal_header == [
        "parameter",
        "value",
        "candidate_count",
        "best_score",
        "mean_score",
        "std_score",
        "worst_score",
        "best_rank",
        "value_rank",
        "at_search_bound",
    ]

    report = (output_dir / "hyper_tune_results.txt").read_text()
    for section in (
        "Search space:",
        "Top 3 candidates:",
        "Per-parameter value summary:",
        "Search-space guidance:",
        "Artifacts:",
        "Weights & Biases logging: disabled",
    ):
        assert section in report

    payload = json.loads((output_dir / "hyper_tune_results.json").read_text())
    assert payload["parameter_marginals"]
    assert payload["artifact_paths"]["parameter_marginals_tsv"] == str(marginals_path)
    assert {entry["parameter"] for entry in payload["parameter_influence"]} == set(
        config.search_space
    )


@pytest.mark.integration
@pytest.mark.output
def test_tune_hyperparameters_random_search_smoke(
    summary_statistics_tsv_tuning, tmp_path
):
    """Random search samples `n_iter` candidates instead of enumerating.

    The space here holds more combinations than `n_iter=1`, so a single
    candidate proves the sampling budget is honoured rather than the grid being
    walked; the best candidate is that one sample.
    """
    config = _tuning_namespace(
        summary_statistics_tsv_tuning,
        tmp_path / "hyper_tune_random_out",
        search_method="random",
        n_iter=1,
        search_space={
            "n_estimators": [5, 10],
            "max_depth": [None, 3],
            "min_samples_split": [2],
            "min_samples_leaf": [1],
            "max_features": ["sqrt"],
            "class_weight": [None],
        },
    )

    result = tune_hyperparameters(config)

    assert result["results"]["search_method"] == "random"
    assert len(result["candidates"]) == 1
    assert result["results"]["best_candidate"]["candidate_index"] == 1


@pytest.mark.parametrize(
    "objective_direction, best_value, best_score, bound_flag",
    [("max", "10", 0.9, "upper_bound"), ("min", "5", 0.4, "lower_bound")],
)


def test_compute_parameter_marginals_summarizes_each_value(
    objective_direction, best_value, best_score, bound_flag
):
    """Values are scored in the objective's direction and ranked against each other."""
    scored = [
        (1, 5, 0.7),
        (2, 5, 0.4),
        (3, 10, 0.9),
        (4, 10, 0.6),
    ]
    ordered = sorted(
        scored, key=lambda item: item[2], reverse=(objective_direction == "max")
    )
    ranked = [
        {
            "candidate_index": candidate_index,
            "n_estimators": n_estimators,
            "cv_score": score,
            "rank": rank,
        }
        for rank, (candidate_index, n_estimators, score) in enumerate(ordered, start=1)
    ]

    marginals = tuning_report.compute_parameter_marginals(
        ranked, ["n_estimators"], objective_direction
    )

    assert [row["value"] for row in marginals] == sorted(
        {"5", "10"}, key=lambda value: 0 if value == best_value else 1
    )
    top = marginals[0]
    assert top["parameter"] == "n_estimators"
    assert top["value"] == best_value
    assert top["candidate_count"] == 2
    assert top["best_score"] == pytest.approx(best_score)
    assert top["best_rank"] == 1
    assert top["value_rank"] == 1
    assert top["at_search_bound"] == bound_flag
    assert marginals[1]["value_rank"] == 2
    assert marginals[1]["at_search_bound"] == ""


def test_parameter_influence_ranks_by_best_score_spread():
    """The parameter whose values separate the objective most is reported first."""
    ranked = [
        {"cv_score": 0.9, "rank": 1, "n_estimators": 10, "max_depth": 5},
        {"cv_score": 0.7, "rank": 2, "n_estimators": 5, "max_depth": 5},
        {"cv_score": 0.6, "rank": 3, "n_estimators": 10, "max_depth": 3},
        {"cv_score": 0.4, "rank": 4, "n_estimators": 5, "max_depth": 3},
    ]

    marginals = tuning_report.compute_parameter_marginals(
        ranked, ["n_estimators", "max_depth"], "max"
    )
    influence = tuning_report.parameter_influence(marginals)

    assert [row["parameter"] for row in influence] == ["max_depth", "n_estimators"]
    assert influence[0]["best_score_spread"] == pytest.approx(0.3)
    assert influence[0]["mean_score_spread"] == pytest.approx(0.3)
    assert influence[1]["best_score_spread"] == pytest.approx(0.2)
    assert influence[1]["mean_score_spread"] == pytest.approx(0.2)


@pytest.mark.output
def test_search_space_guidance_flags_a_dimension_with_no_effect():
    """A parameter whose values all score alike is called out as inert."""
    ranked = [
        {"cv_score": 0.8, "rank": 1, "max_features": "sqrt"},
        {"cv_score": 0.8, "rank": 2, "max_features": "log2"},
    ]
    marginals = tuning_report.compute_parameter_marginals(ranked, ["max_features"], "max")

    lines = "\n".join(tuning_report.format_search_space_guidance(marginals, ranked))

    assert "spread 0.000000" in lines
    assert (
        "max_features: best 0.000000, mean 0.000000 - no effect on the objective"
        in lines
    )


class _StubWandbRun:
    def __init__(self):
        self.summary = {}
        self.finished = False

    def finish(self):
        self.finished = True


class _StubWandb:
    def __init__(self, run):
        self._run = run
        self.logged = []
        self.defined_metrics = []

    def init(self, **_kwargs):
        return self._run

    def log(self, payload):
        self.logged.append(payload)

    def define_metric(self, name, **_kwargs):
        self.defined_metrics.append(name)

    def Table(self, columns, data):  # noqa: N802 - mirrors the wandb API
        return {"columns": columns, "data": data}


@pytest.mark.output
@pytest.mark.parametrize("use_wandb", [False, True])
def test_create_run_logger_routes_on_the_use_wandb_choice(
    use_wandb, tmp_path, monkeypatch
):
    """`use_wandb: false` yields an inert logger; `true` forwards to the W&B run."""
    run = _StubWandbRun()
    stub = _StubWandb(run)
    monkeypatch.setattr(hyper_tune, "_import_wandb", lambda: stub)
    config = _tuning_namespace(
        tmp_path / "input.tsv", tmp_path / "out", use_wandb=use_wandb
    )

    logger = hyper_tune._create_run_logger(config, tmp_path, 2, 2, use_wandb)
    logger.log({"candidate/cv_score": 0.5})
    logger.set_summary("artifact_dir", str(tmp_path))
    logger.finish()

    assert logger.enabled is use_wandb
    assert stub.logged == ([{"candidate/cv_score": 0.5}] if use_wandb else [])
    assert run.summary == ({"artifact_dir": str(tmp_path)} if use_wandb else {})
    assert run.finished is use_wandb
    assert (tmp_path / "wandb").exists() is use_wandb


@pytest.mark.config
def test_import_wandb_raises_config_error_when_unavailable(monkeypatch):
    """A missing wandb install surfaces as an actionable config error, not ImportError."""
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "wandb":
            raise ModuleNotFoundError("No module named 'wandb'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ConfigError, match="use_wandb: false"):
        hyper_tune._import_wandb()


@pytest.mark.integration
@pytest.mark.output
def test_tune_hyperparameters_routes_bulk_artifacts_to_wandb(
    summary_statistics_tsv_tuning, tmp_path, monkeypatch
):
    """With W&B on, the bulk tables and results JSON go to the run, not to disk."""
    run = _StubWandbRun()
    stub = _StubWandb(run)
    monkeypatch.setattr(hyper_tune, "_import_wandb", lambda: stub)
    output_dir = tmp_path / "hyper_tune_wandb_out"
    config = _tuning_namespace(
        summary_statistics_tsv_tuning, output_dir, use_wandb=True
    )

    result = tune_hyperparameters(config)

    for kept in (
        "hyper_tune_best_model.pkl",
        "hyper_tune_results.txt",
        "figures/hyper_tune_search_report.png",
    ):
        assert (output_dir / kept).exists()
    for skipped in (
        "hyper_tune_results.json",
        "hyper_tune_results.tsv",
        "hyper_tune_parameter_marginals.tsv",
        "predictions.tsv",
    ):
        assert not (output_dir / skipped).exists()

    assert result["results_json_path"] is None
    assert result["results_tsv_path"] is None
    assert result["marginals_tsv_path"] is None
    assert result["predictions_path"] is None
    assert result["results"]["artifact_paths"] == {
        "best_model": str(output_dir / "hyper_tune_best_model.pkl"),
        "results_txt": str(output_dir / "hyper_tune_results.txt"),
        "search_report_plot": str(
            output_dir / "figures" / "hyper_tune_search_report.png"
        ),
    }

    logged_keys = {key for payload in stub.logged for key in payload}
    assert {
        "tables/ranked_candidates",
        "tables/parameter_marginals",
        "tables/predictions",
    } <= logged_keys

    assert run.summary["bulk_artifacts_written_locally"] is False
    payload = json.loads(run.summary["results_json"])
    assert payload["parameter_marginals"]
    assert payload["artifact_paths"] == result["results"]["artifact_paths"]

    report = (output_dir / "hyper_tune_results.txt").read_text()
    assert "Weights & Biases logging: enabled" in report
    assert "logged to the Weights & Biases run" in report
