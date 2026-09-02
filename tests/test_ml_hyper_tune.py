import argparse
import json
from pathlib import Path

import pytest

from ghostparser.config import ConfigError
from ghostparser.ml import hyper_tune, tuning_report
from ghostparser.ml.hyper_tune import load_hyper_tune_config, tune_hyperparameters


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
        random_state=7,
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


def test_load_hyper_tune_config_accepts_hyperparameter_tuning_section(tmp_path):
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
                "max_features": ["sqrt"],
                "class_weight": [None],
            },
        },
    )

    config = load_hyper_tune_config(str(config_path))

    assert config["model_name"] == "random_forest"
    assert config["search_method"] == "grid"
    assert config["objective_metric"] == "exact_match_accuracy"
    assert config["search_space"]["n_estimators"] == [5, 10]
    assert config["use_wandb"] is False


@pytest.mark.parametrize(
    "use_wandb, detailed_payloads",
    [(False, None), (True, None), (True, True)],
)
def test_load_hyper_tune_config_fills_model_defaults(
    use_wandb, detailed_payloads, tmp_path
):
    """Keys absent from the tuning block take their defaults; present ones win."""
    tuning = {
        "model": "random_forest",
        "use_wandb": use_wandb,
        "search_space": {"n_estimators": [5, 10]},
    }
    if detailed_payloads is not None:
        tuning["wandb_detailed_payloads"] = detailed_payloads

    config_path = _write_config(tmp_path, "hyper_tune_defaults.json", tuning)

    config = load_hyper_tune_config(str(config_path))

    assert config["class_weight"] is None
    assert config["max_features"] == "sqrt"
    assert config["min_samples_split"] == 2
    assert config["overwrite"] is True
    assert config["use_wandb"] is use_wandb
    assert config["wandb_detailed_payloads"] is bool(detailed_payloads)


@pytest.mark.parametrize(
    "tuning, expected_message",
    [
        pytest.param(
            {"model": "random_forest", "search_space": {"n_estimators": [5]}},
            "use_wandb",
            id="missing_use_wandb",
        ),
        pytest.param(
            {
                "model": "random_forest",
                "use_wandb": "yes",
                "search_space": {"n_estimators": [5]},
            },
            "must be a boolean",
            id="non_boolean_use_wandb",
        ),
        pytest.param(
            {
                "model": "random_forest",
                "use_wandb": False,
                "wandb_detailed_payloads": True,
                "search_space": {"n_estimators": [5]},
            },
            "wandb_detailed_payloads requires",
            id="detailed_payloads_without_wandb",
        ),
    ],
)
def test_load_hyper_tune_config_rejects_invalid_wandb_choice(
    tuning, expected_message, tmp_path
):
    """`use_wandb` must be spelled out as a boolean and must gate the extra payloads."""
    config_path = _write_config(tmp_path, "hyper_tune_wandb.json", tuning)

    with pytest.raises(ConfigError, match=expected_message):
        load_hyper_tune_config(str(config_path))


@pytest.mark.parametrize(
    "sample_name, expected_model",
    [
        ("hyperparameter_tuning_random_forest.yaml", "random_forest"),
        ("hyperparameter_tuning_multi_knn.json", "multi_knn"),
    ],
)
def test_shipped_tuning_sample_configs_resolve(sample_name, expected_model):
    """The shipped tuner samples load and spell out every required key.

    Guards against a sample drifting out of step with the validator, which would
    leave users copying a config the loader rejects. `use_wandb` has no default,
    so a successful load proves the key is present in both config formats.
    """
    sample_dir = Path(__file__).resolve().parents[1] / "sample_configs"
    config = load_hyper_tune_config(str(sample_dir / sample_name))

    assert config["model_name"] == expected_model
    assert config["use_wandb"] is False
    assert config["wandb_detailed_payloads"] is False
    assert set(config["search_space"]) <= set(
        hyper_tune.MODEL_SEARCH_KEYS[expected_model]
    )


def test_load_hyper_tune_config_rejects_evaluation_section(tmp_path):
    config_path = _write_config(
        tmp_path,
        "hyper_tune_with_evaluation.json",
        {
            "model": "random_forest",
            "use_wandb": False,
            "search_space": {"n_estimators": [5, 10]},
        },
        evaluation={"metrics": "all"},
    )

    with pytest.raises(ConfigError, match="evaluation"):
        load_hyper_tune_config(str(config_path))


def test_load_hyper_tune_config_requires_hyperparameter_tuning_section(tmp_path):
    config_path = _write_config(tmp_path, "hyper_tune_missing_section.json", None)

    with pytest.raises(ConfigError, match="hyperparameter_tuning"):
        load_hyper_tune_config(str(config_path))


def test_tune_hyperparameters_requires_explicit_use_wandb(
    summary_statistics_tsv_tuning, tmp_path
):
    """The programmatic entry point enforces the same deliberate choice."""
    config = _tuning_namespace(summary_statistics_tsv_tuning, tmp_path / "no_choice")
    del config.use_wandb

    with pytest.raises(ConfigError, match="use_wandb"):
        tune_hyperparameters(config)


def test_tune_hyperparameters_grid_search_smoke(summary_statistics_tsv_tuning, tmp_path):
    output_dir = tmp_path / "hyper_tune_out"
    config = _tuning_namespace(summary_statistics_tsv_tuning, output_dir)

    result = tune_hyperparameters(config)

    assert result["results"]["model_name"] == "random_forest"
    assert result["results"]["search_method"] == "grid"
    assert result["results"]["use_wandb"] is False
    assert len(result["candidates"]) == 2
    assert (output_dir / "hyper_tune_best_model.pkl").exists()
    assert (output_dir / "hyper_tune_results.json").exists()
    assert not (output_dir / "wandb").exists()


def test_tune_hyperparameters_random_search_smoke(
    summary_statistics_tsv_tuning, tmp_path
):
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


def test_tune_hyperparameters_writes_local_navigation_artifacts(
    summary_statistics_tsv_tuning, tmp_path
):
    """Without W&B the run still emits ranked TSVs, marginals, plots and a report."""
    output_dir = tmp_path / "hyper_tune_local_out"
    config = _tuning_namespace(summary_statistics_tsv_tuning, output_dir)

    result = tune_hyperparameters(config)

    marginals_path = output_dir / "hyper_tune_parameter_marginals.tsv"
    plot_path = output_dir / "hyper_tune_search_report.png"
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
        "hyper_tune_search_report.png",
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
        "search_report_plot": str(output_dir / "hyper_tune_search_report.png"),
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


@pytest.mark.parametrize(
    "values, expected",
    [
        (["sqrt", "none", None, 4], ["sqrt", None, None, 4]),
        ("None", None),
    ],
)
def test_load_hyper_tune_config_normalizes_search_space_values(
    values, expected, tmp_path
):
    """Search-space candidates get the same per-value rules as the `model` block."""
    config_path = _write_config(
        tmp_path,
        "hyper_tune_search_values.json",
        {
            "model": "random_forest",
            "use_wandb": False,
            "search_space": {"max_features": values},
        },
    )

    config = load_hyper_tune_config(str(config_path))

    assert config["search_space"]["max_features"] == expected


def test_load_hyper_tune_config_rejects_invalid_search_space_value(tmp_path):
    """An invalid candidate is named against its search_space key, before any fit."""
    config_path = _write_config(
        tmp_path,
        "hyper_tune_bad_search_value.json",
        {
            "model": "random_forest",
            "use_wandb": False,
            "search_space": {"max_features": ["sqrt", "auto"]},
        },
    )

    with pytest.raises(ConfigError) as excinfo:
        load_hyper_tune_config(str(config_path))

    message = str(excinfo.value)
    assert "hyperparameter_tuning.search_space.max_features" in message
    assert "'auto'" in message
