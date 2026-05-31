from __future__ import annotations

import argparse
import json

import numpy as np

from ghostparser.ml.random_forest import _parse_classes, train_random_forest


def test_parse_classes_returns_binary_matrix():
    targets, labels = _parse_classes(["101001", "010010"])

    assert labels == ["101001", "010010"]
    assert targets.shape == (2, 6)
    assert targets.tolist() == [[1, 0, 1, 0, 0, 1], [0, 1, 0, 0, 1, 0]]


def test_train_random_forest_smoke(summary_statistics_tsv, tmp_path):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="classes",
        test_size=0.25,
        cv_folds=3,
        random_state=7,
        rare_class_policy="warn_reduce_cv",
        n_estimators=25,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        class_weight=None,
        n_jobs=1,
    )

    result = train_random_forest(config)

    metrics = result["metrics"]
    assert metrics["objective"] == "multi-label classification"
    assert "exact_match_accuracy" in metrics["diagnostic_metrics"]
    assert "hamming_loss" in metrics["primary_metrics"]
    assert metrics["class_distribution"]["overall"]
    assert metrics["bit_distribution"]["overall"]
    assert metrics["cv"] is not None

    metrics_json = json.loads((tmp_path / "ml_out" / "random_forest_metrics.json").read_text())
    assert metrics_json["exact_match_is_diagnostic"] is True

    predictions_path = tmp_path / "ml_out" / "predictions.tsv"
    assert predictions_path.exists()

    model_path = tmp_path / "ml_out" / "random_forest_model.pkl"
    assert model_path.exists()

    label_map = json.loads((tmp_path / "ml_out" / "label_map.json").read_text())
    assert label_map["bit_labels"][0] == "ghost_into_A"


def test_train_random_forest_creates_bitwise_metrics_report(summary_statistics_tsv, tmp_path):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="classes",
        test_size=0.25,
        cv_folds=2,
        random_state=11,
        rare_class_policy="warn_skip_cv",
        n_estimators=15,
        max_depth=4,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        class_weight=None,
        n_jobs=1,
    )

    result = train_random_forest(config)
    metrics_text = (tmp_path / "ml_out" / "random_forest_metrics.txt").read_text()

    assert "Primary objective: multi-label classification" in metrics_text
    assert "Exact-match accuracy is diagnostic" in metrics_text
    assert "Per-bit metrics:" in metrics_text
    assert result["metrics"]["primary_metrics"]["bitwise_accuracy"] <= 1.0
    assert np.isfinite(result["metrics"]["primary_metrics"]["hamming_loss"])