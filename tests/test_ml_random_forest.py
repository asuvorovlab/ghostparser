from __future__ import annotations

import argparse
import csv
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
        target_column="class",
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
        overwrite=True,
    )

    result = train_random_forest(config)

    metrics = result["metrics"]
    assert metrics["objective"] == "multi-label classification"
    assert "exact_match_accuracy" in metrics["diagnostic_metrics"]
    assert "hamming_loss" in metrics["primary_metrics"]
    assert metrics["dataset_summary"]["class_distribution"]["overall"]
    assert metrics["dataset_summary"]["bit_distribution"]["overall"]
    assert metrics["cv"] is not None
    assert "confusion_matrix_64_classes" in metrics

    metrics_json = json.loads(
        (tmp_path / "ml_out" / "random_forest_overall_metrics.json").read_text()
    )
    assert metrics_json["exact_match_is_diagnostic"] is True
    assert (
        metrics_json["dataset_summary"]["label_map"]["bit_labels"][0] == "ghost_into_A"
    )
    assert metrics_json["timings_seconds"]["total"] >= 0

    predictions_path = tmp_path / "ml_out" / "predictions.tsv"
    assert predictions_path.exists()
    with open(predictions_path, "r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    assert "matched_label_count" in header

    model_path = tmp_path / "ml_out" / "random_forest_model.pkl"
    assert model_path.exists()
    assert (tmp_path / "ml_out" / "random_forest_confusion_matrices.png").exists()
    assert (
        tmp_path / "ml_out" / "random_forest_confusion_matrix_64_classes.png"
    ).exists()
    assert not (tmp_path / "ml_out" / "label_map.json").exists()
    assert not (tmp_path / "ml_out" / "class_distribution.tsv").exists()
    assert not (tmp_path / "ml_out" / "bit_distribution.tsv").exists()


def test_train_random_forest_creates_bitwise_metrics_report(
    summary_statistics_tsv, tmp_path
):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="class",
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
        overwrite=True,
    )

    result = train_random_forest(config)
    metrics_text = (tmp_path / "ml_out" / "random_forest_metrics.txt").read_text()

    assert "Primary objective: multi-label classification" in metrics_text
    assert "Exact-match accuracy is diagnostic" in metrics_text
    assert "Per-bit metrics:" in metrics_text
    assert "Timings (seconds):" in metrics_text
    assert result["metrics"]["primary_metrics"]["bitwise_accuracy"] <= 1.0
    assert np.isfinite(result["metrics"]["primary_metrics"]["hamming_loss"])


def test_metrics_txt_leads_with_hyperparameters(summary_statistics_tsv, tmp_path):
    """Every hyperparameter used for the run is printed above the metrics."""
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="class",
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
        overwrite=True,
    )

    result = train_random_forest(config)
    text = open(result["metrics_txt_path"], encoding="utf-8").read()

    # The block precedes the metrics it describes.
    assert text.index("Hyperparameters:") < text.index("Test metrics:")

    # Parse the block back into a mapping; names are padded to a common width,
    # so split on whitespace rather than asserting exact column positions.
    block = text.split("Hyperparameters:\n", 1)[1].split("\n\n", 1)[0]
    printed = dict(line.split(maxsplit=1) for line in block.splitlines())

    # Values come from the config, not from estimator defaults.
    assert printed["n_estimators"] == "25"
    assert printed["random_state"] == "7"
    assert printed["test_size"] == "0.25"
    # An unset knob renders as `none`, distinct from an empty value.
    assert printed["max_depth"] == "none"
    # The same mapping is machine-readable in the JSON payload.
    hyperparameters = result["metrics"]["hyperparameters"]
    assert hyperparameters["n_estimators"] == 25
    assert hyperparameters["max_depth"] is None
    assert hyperparameters["cv_folds_requested"] == 3
