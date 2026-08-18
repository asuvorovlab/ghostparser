import argparse
import csv
import json

from ghostparser.ml.multi_knn import _build_model, train_multi_knn


def test_multi_knn_train_smoke(summary_statistics_tsv, tmp_path):
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(tmp_path / "ml_out"),
        target_column="class",
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
        report_feature_importance=False,
        overwrite=True,
    )

    result = train_multi_knn(config)

    metrics = result["metrics"]
    assert metrics["classifier"] == "multi_knn"
    assert "exact_match_accuracy" in metrics["diagnostic_metrics"]
    assert "hamming_loss" in metrics["primary_metrics"]
    assert metrics["dataset_summary"]["class_distribution"]["overall"]
    assert metrics["dataset_summary"]["bit_distribution"]["overall"]
    assert metrics["cv"] is not None
    assert metrics["hyperparameters"]["n_neighbors_requested"] == 5
    assert "confusion_matrix_64_classes" in metrics

    metrics_json = json.loads(
        (tmp_path / "ml_out" / "multi_knn_overall_metrics.json").read_text()
    )
    assert metrics_json["exact_match_is_diagnostic"] is True
    assert metrics_json["hyperparameters"]["n_neighbors_effective"] >= 1
    assert metrics_json["dataset_summary"]["split"]["train_rows"] > 0

    predictions_path = tmp_path / "ml_out" / "predictions.tsv"
    assert predictions_path.exists()
    with open(predictions_path, "r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    assert "matched_label_count" in header

    model_path = tmp_path / "ml_out" / "multi_knn_model.pkl"
    assert model_path.exists()
    assert (tmp_path / "ml_out" / "multi_knn_confusion_matrices.png").exists()
    assert (tmp_path / "ml_out" / "multi_knn_confusion_matrix_64_classes.png").exists()
    assert not (tmp_path / "ml_out" / "label_map.json").exists()
    assert not (tmp_path / "ml_out" / "class_distribution.tsv").exists()
    assert not (tmp_path / "ml_out" / "bit_distribution.tsv").exists()


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

    # The clamp must reach the estimator, not only the reported value: asking
    # KNN for more neighbors than training samples raises at fit time.
    assert effective_n_neighbors == 2
    assert model.estimator.n_neighbors == 2

    # Below the clamp the request passes through untouched.
    _, unclamped = _build_model(config, train_size=50)
    assert unclamped == 20
