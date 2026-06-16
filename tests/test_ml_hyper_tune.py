from __future__ import annotations

import argparse
import json

import pytest

from ghostparser.config import ConfigError
from ghostparser.ml.hyper_tune import load_hyper_tune_config, tune_hyperparameters


def test_load_hyper_tune_config_accepts_hyperparameter_tuning_section(tmp_path):
    config_path = tmp_path / "hyper_tune.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/hyper_tune_out",
                "hyperparameter_tuning": {
                    "model": "random_forest",
                    "method": "grid",
                    "objective": "exact_match_accuracy",
                    "top_k": 5,
                    "max_candidates": 100,
                    "search_space": {
                        "n_estimators": [5, 10],
                        "max_depth": [None],
                        "min_samples_split": [2],
                        "min_samples_leaf": [1],
                        "max_features": ["sqrt"],
                        "class_weight": [None],
                    },
                },
            }
        )
    )

    config = load_hyper_tune_config(str(config_path))

    assert config["model_name"] == "random_forest"
    assert config["search_method"] == "grid"
    assert config["objective_metric"] == "exact_match_accuracy"
    assert config["search_space"]["n_estimators"] == [5, 10]


def test_load_hyper_tune_config_fills_model_defaults(tmp_path):
    config_path = tmp_path / "hyper_tune_defaults.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/hyper_tune_out",
                "hyperparameter_tuning": {
                    "model": "random_forest",
                    "search_space": {
                        "n_estimators": [5, 10],
                    },
                },
            }
        )
    )

    config = load_hyper_tune_config(str(config_path))

    assert config["class_weight"] is None
    assert config["max_features"] == "sqrt"
    assert config["min_samples_split"] == 2


def test_load_hyper_tune_config_rejects_evaluation_section(tmp_path):
    config_path = tmp_path / "hyper_tune_with_evaluation.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/hyper_tune_out",
                "evaluation": {"metrics": "all"},
                "hyperparameter_tuning": {
                    "model": "random_forest",
                    "search_space": {"n_estimators": [5, 10]},
                },
            }
        )
    )

    with pytest.raises(ConfigError, match="evaluation"):
        load_hyper_tune_config(str(config_path))


def test_load_hyper_tune_config_requires_hyperparameter_tuning_section(tmp_path):
    config_path = tmp_path / "hyper_tune_missing_section.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "./results/summary_statistics.tsv",
                "output_dir": "./results/hyper_tune_out",
            }
        )
    )

    with pytest.raises(ConfigError, match="hyperparameter_tuning"):
        load_hyper_tune_config(str(config_path))


def test_tune_hyperparameters_grid_search_smoke(
    summary_statistics_tsv_tuning, tmp_path, capsys
):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv_tuning),
        output_dir=str(tmp_path / "hyper_tune_out"),
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
        search_space={
            "n_estimators": [5, 10],
            "max_depth": [None],
            "min_samples_split": [2],
            "min_samples_leaf": [1],
            "max_features": ["sqrt"],
            "class_weight": [None],
        },
    )

    result = tune_hyperparameters(config)
    captured = capsys.readouterr().out

    assert result["results"]["model_name"] == "random_forest"
    assert result["results"]["search_method"] == "grid"
    assert len(result["candidates"]) == 2
    assert (tmp_path / "hyper_tune_out" / "hyper_tune_best_model.pkl").exists()
    assert (tmp_path / "hyper_tune_out" / "hyper_tune_results.json").exists()
    assert "[hyper_tune] Starting hyperparameter tuning" in captured
    assert "Will evaluate 2 candidate cases across 2 CV folds" in captured
    assert "[1/2] evaluating" in captured
    assert "Finished in" in captured


def test_tune_hyperparameters_random_search_smoke(
    summary_statistics_tsv_tuning, tmp_path
):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv_tuning),
        output_dir=str(tmp_path / "hyper_tune_random_out"),
        target_column="class",
        test_size=0.25,
        cv_folds=2,
        rare_class_policy="warn_reduce_cv",
        random_state=7,
        n_jobs=1,
        overwrite=True,
        model_name="random_forest",
        search_method="random",
        objective_metric="exact_match_accuracy",
        objective_direction="max",
        top_k=3,
        n_iter=1,
        max_candidates=10,
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
