import argparse
import csv
import json

import pytest

from ghostparser.ml.multi_knn import _build_model, train_multi_knn
from ghostparser.ml.random_forest import train_random_forest

_MODELS = {
    "random_forest": (
        train_random_forest,
        dict(
            n_estimators=25, max_depth=None, min_samples_split=2,
            min_samples_leaf=1, max_features="sqrt", class_weight=None,
        ),
    ),
    "multi_knn": (
        train_multi_knn,
        dict(
            n_neighbors=5, weights="uniform", algorithm="auto", leaf_size=30,
            metric="minkowski", p=2,
        ),
    ),
}


@pytest.mark.integration
@pytest.mark.output
@pytest.mark.parametrize("model", sorted(_MODELS))
def test_trainer_runs_end_to_end_and_writes_every_artifact(
    model, summary_statistics_tsv, tmp_path
):
    """Each trainer fits, scores and persists, writing every artifact it promises.

    The one smoke test per trainer: both metric tiers, the dataset summary,
    cross-validation, the 64-class matrix, one timing per stage, the
    predictions TSV, the model pickle and the four figures under ``figures/``.
    """
    train, hyperparameters = _MODELS[model]
    output_dir = tmp_path / "ml_out"
    config = argparse.Namespace(
        input_path=str(summary_statistics_tsv),
        output_dir=str(output_dir),
        target_column="class",
        test_size=0.25,
        cv_folds=3,
        seed=7,
        rare_class_policy="warn_reduce_cv",
        n_jobs=1,
        overwrite=True,
        **hyperparameters,
    )

    metrics = train(config)["metrics"]

    assert "exact_match_accuracy" in metrics["diagnostic_metrics"]
    assert "hamming_loss" in metrics["primary_metrics"]
    assert metrics["dataset_summary"]["class_distribution"]["overall"]
    assert metrics["dataset_summary"]["bit_distribution"]["overall"]
    assert metrics["cv"] is not None
    assert "confusion_matrix_64_classes" in metrics

    metrics_json = json.loads((output_dir / f"{model}_overall_metrics.json").read_text())
    assert metrics_json["exact_match_is_diagnostic"] is True
    assert metrics_json["dataset_summary"]["label_map"]["bit_labels"][0] == "ghost_into_A"
    assert {"load", "split", "fit_and_predict", "artifact_write", "total"} <= set(
        metrics_json["timings_seconds"]
    )

    with open(output_dir / "predictions.tsv", encoding="utf-8", newline="") as handle:
        assert "matched_label_count" in next(csv.reader(handle, delimiter="\t"))
    assert (output_dir / f"{model}_model.pkl").exists()
    for figure in (
        "confusion_matrices",
        "confusion_matrix_64_classes",
        "per_bit_accuracy",
        "per_class_accuracy",
    ):
        assert (output_dir / "figures" / f"{model}_{figure}.png").exists()


def test_multi_knn_caps_neighbors_to_the_training_size():
    """`n_neighbors` is clamped to the training-set size, and only when needed.

    The clamp must reach the estimator, not only the reported value: asking
    KNN for more neighbors than training samples raises at fit time.
    """
    config = argparse.Namespace(
        n_neighbors=20, weights="uniform", algorithm="auto", leaf_size=30,
        metric="minkowski", p=2, n_jobs=1,
    )

    model, effective_n_neighbors = _build_model(config, train_size=2)
    assert effective_n_neighbors == 2
    assert model.estimator.n_neighbors == 2

    assert _build_model(config, train_size=50)[1] == 20
