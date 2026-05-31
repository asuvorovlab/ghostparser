from __future__ import annotations

import argparse
import json

import numpy as np

from ghostparser.ml.multi_knn import _build_model, train_multi_knn


def test_multi_knn_train_smoke(summary_statistics_tsv, tmp_path):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="classes",
        test_size=0.25,
        cv_folds=3,
        random_state=7,
        rare_class_policy="warn_reduce_cv",
        n_neighbors=5,
        weights="uniform",
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        n_jobs=1,
        n_estimators=25,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        class_weight=None,
    )

    result = train_multi_knn(config)

    metrics = result["metrics"]
    assert metrics["classifier"] == "multi_knn"
    assert "exact_match_accuracy" in metrics["diagnostic_metrics"]
    assert "hamming_loss" in metrics["primary_metrics"]
    assert metrics["class_distribution"]["overall"]
    assert metrics["bit_distribution"]["overall"]
    assert metrics["cv"] is not None
    assert metrics["knn"]["configured_n_neighbors"] == 5

    metrics_json = json.loads((tmp_path / "ml_out" / "multi_knn_metrics.json").read_text())
    assert metrics_json["exact_match_is_diagnostic"] is True
    assert metrics_json["knn"]["effective_n_neighbors"] >= 1

    predictions_path = tmp_path / "ml_out" / "predictions.tsv"
    assert predictions_path.exists()

    model_path = tmp_path / "ml_out" / "multi_knn_model.pkl"
    assert model_path.exists()

    label_map = json.loads((tmp_path / "ml_out" / "label_map.json").read_text())
    assert label_map["bit_labels"][0] == "ghost_into_A"


def test_multi_knn_build_model_caps_neighbors_to_training_size():
    config = argparse.Namespace(
        n_neighbors=20,
        weights="uniform",
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        n_jobs=1,
    )

    model, effective_n_neighbors = _build_model(config, train_size=2)

    assert effective_n_neighbors == 2
    assert model is not None


def test_multi_knn_metrics_report_mentions_effective_neighbors(summary_statistics_tsv, tmp_path):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="classes",
        test_size=0.25,
        cv_folds=2,
        random_state=11,
        rare_class_policy="warn_skip_cv",
        n_neighbors=20,
        weights="distance",
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        n_jobs=1,
        n_estimators=25,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        class_weight=None,
    )

    result = train_multi_knn(config)
    metrics_text = (tmp_path / "ml_out" / "multi_knn_metrics.txt").read_text()

    assert "Ghostparser ML multi-label KNN baseline" in metrics_text
    assert "Configured n_neighbors:" in metrics_text
    assert "Effective n_neighbors:" in metrics_text
    assert result["metrics"]["knn"]["configured_n_neighbors"] == 20
    assert result["metrics"]["knn"]["effective_n_neighbors"] <= 20
    assert np.isfinite(result["metrics"]["primary_metrics"]["hamming_loss"])
    assert (tmp_path / "ml_out" / "feature_importances.tsv").exists()