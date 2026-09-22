import argparse
import csv
import json

import pytest

from ghostparser.ml.random_forest import train_random_forest

pytestmark = [pytest.mark.integration, pytest.mark.output]


def test_train_random_forest_smoke(summary_statistics_tsv, tmp_path):
    """Training runs end to end and writes every artifact it promises.

    The one smoke test for this entry point: it proves the forest path fits,
    scores and persists, and pins the metrics-field and output-file contracts --
    both metric tiers, the dataset summary, the 64-class matrix, one timing per
    stage, and the four written artifacts.
    """
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="class",
        test_size=0.25,
        cv_folds=3,
        seed=7,
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
    assert set(metrics_json["timings_seconds"]) == {
        "load", "split", "cross_validation", "fit_and_predict",
        "feature_importance", "prediction_row_build", "artifact_write",
        "total",
    }

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
