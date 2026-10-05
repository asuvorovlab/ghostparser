"""Multi-label KNN baseline for Ghostparser summary statistics."""

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.multioutput import MultiOutputClassifier
from sklearn.neighbors import KNeighborsClassifier

from ..cli_config import run_cli
from ..config import prepare_output_directory
from . import ml_utils as shared
from .config import (
    DEFAULT_FEATURE_IMPORTANCE_CORRELATION_THRESHOLD,
    DEFAULT_N_JOBS,
    build_trainer_argument_parser,
    resolve_trainer_runtime_args,
)

BIT_LABELS = shared.BIT_LABELS
BIT_COUNT = shared.BIT_COUNT


def _build_model(
    config: argparse.Namespace, train_size: int
) -> tuple[MultiOutputClassifier, int]:
    """Build the one-vs-rest KNN, capping ``n_neighbors`` at the training size.

    Args:
        config: The resolved trainer config.
        train_size: Number of training rows.

    Returns:
        ``(model, effective_n_neighbors)``.
    """
    effective_n_neighbors = max(1, min(int(config.n_neighbors), int(train_size)))
    classifier = KNeighborsClassifier(
        n_neighbors=effective_n_neighbors,
        weights=config.weights,
        algorithm=config.algorithm,
        leaf_size=config.leaf_size,
        metric=config.metric,
        p=config.p,
        n_jobs=config.n_jobs if config.n_jobs is not None else DEFAULT_N_JOBS,
    )
    model = MultiOutputClassifier(
        classifier,
        n_jobs=config.n_jobs if config.n_jobs is not None else DEFAULT_N_JOBS,
    )
    return model, effective_n_neighbors


def _cross_validate(
    x_train: np.ndarray,
    y_train: np.ndarray,
    labels_train: np.ndarray,
    folds: int,
    random_state: int | None,
    config: argparse.Namespace,
) -> dict:
    """Cross-validate the KNN over stratified folds of the training partition.

    Args:
        x_train: The training features.
        y_train: The training binary targets.
        labels_train: The training combination labels, for stratification.
        folds: Number of folds.
        random_state: Seed for the fold shuffle.
        config: The resolved trainer config.

    Returns:
        A dict with the per-fold metrics and their ``aggregate`` means and
        standard deviations.
    """
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random_state)
    fold_metrics: list[dict] = []
    for fold_index, (fit_index, val_index) in enumerate(
        splitter.split(x_train, labels_train), start=1
    ):
        model, effective_n_neighbors = _build_model(config, len(fit_index))
        model.fit(x_train[fit_index], y_train[fit_index])
        fold_result = shared.evaluate_predictions(
            y_train[val_index], model.predict(x_train[val_index])
        )
        fold_metrics.append(
            {
                "fold": fold_index,
                "effective_n_neighbors": effective_n_neighbors,
                **{
                    key: value
                    for key, value in fold_result.items()
                    if key not in {"per_bit", "classification_report"}
                },
            }
        )

    aggregate: dict[str, float] = {}
    for key in [
        "exact_match_accuracy",
        "hamming_loss",
        "bitwise_accuracy",
        "micro_f1",
        "macro_f1",
        "weighted_f1",
    ]:
        values = np.asarray([fold[key] for fold in fold_metrics], dtype=float)
        aggregate[f"{key}_mean"] = float(values.mean())
        aggregate[f"{key}_std"] = float(values.std(ddof=0))
    return {"folds": fold_metrics, "aggregate": aggregate}


def train_multi_knn(config: argparse.Namespace) -> dict:
    """Train and evaluate the multi-label KNN baseline end to end.

    Args:
        config: The resolved trainer config.

    Returns:
        A dict of the metrics and the paths of the artifacts written.
    """
    run_start = time.perf_counter()
    load_start = time.perf_counter()
    rows = shared.read_tsv_rows(config.input_path)
    matrix = shared.rows_to_matrix(rows, config.target_column)
    labels = shared.combination_labels(matrix.train_labels)
    load_seconds = time.perf_counter() - load_start

    split_start = time.perf_counter()
    x_train, x_test, y_train, y_test, labels_train, labels_test, split_notes = (
        shared.split_dataset(
            matrix.train_features,
            matrix.train_targets,
            labels,
            config.test_size,
            config.seed,
        )
    )
    split_seconds = time.perf_counter() - split_start

    def model_factory(train_size: int) -> tuple[MultiOutputClassifier, int]:
        return _build_model(config, train_size)

    cv_start = time.perf_counter()
    cv_folds, cv_warnings = shared.auto_cv_folds(
        labels_train, config.cv_folds, config.rare_class_policy
    )
    cv_results = {"folds": [], "aggregate": {}}
    if cv_folds is not None:
        cv_results = _cross_validate(
            x_train, y_train, labels_train, cv_folds, config.seed, config
        )
    cv_seconds = time.perf_counter() - cv_start

    fit_start = time.perf_counter()
    model, effective_n_neighbors = model_factory(len(x_train))
    model.fit(x_train, y_train)
    test_predictions = model.predict(x_test)
    test_metrics = shared.evaluate_predictions(y_test, test_predictions)
    fit_predict_seconds = time.perf_counter() - fit_start

    output_dir = Path(
        prepare_output_directory(config.output_dir, overwrite=config.overwrite)
    )

    metric_set = getattr(config, "evaluation_metrics", "all")
    report_class_distribution = getattr(config, "report_class_distribution", True)
    report_confusion_matrix = getattr(config, "report_confusion_matrix", True)
    report_feature_importance = getattr(config, "report_feature_importance", True)
    # A neighbours classifier measures no impurity, so it has to permute; mdi
    # is rejected rather than quietly swapped for something else.
    feature_importance_method = (
        getattr(config, "feature_importance_method", None) or "permutation"
    )
    correlation_threshold = getattr(
        config,
        "feature_importance_correlation_threshold",
        DEFAULT_FEATURE_IMPORTANCE_CORRELATION_THRESHOLD,
    )
    save_predictions = getattr(config, "save_predictions", True)
    include_diagnostic = metric_set in {"diagnostic", "all"}
    include_per_bit = metric_set in {"per_bit", "all"}

    dataset_summary = shared.build_dataset_summary(
        matrix.train_labels,
        list(labels_train),
        list(labels_test),
        y_train,
        y_test,
        split_notes,
        config.test_size,
        config.seed,
    )

    feature_start = time.perf_counter()
    feature_rows = None
    if report_feature_importance:
        feature_rows = shared.feature_importance_rows(
            model,
            method=feature_importance_method,
            feature_names=matrix.feature_names,
            features=x_test,
            targets=y_test,
            seed=config.seed,
            correlation_threshold=correlation_threshold,
        )
    feature_importance_seconds = time.perf_counter() - feature_start

    prediction_start = time.perf_counter()
    prediction_rows = None
    if save_predictions:
        prediction_rows = shared.build_prediction_rows(y_test, test_predictions)
    prediction_build_seconds = time.perf_counter() - prediction_start

    confusion_matrices = None
    confusion_matrix_64_classes = None
    if report_confusion_matrix and include_diagnostic:
        confusion_matrices = shared.build_confusion_matrices(y_test, test_predictions)
        confusion_matrix_64_classes = shared.build_64_class_confusion_matrix(
            y_test, test_predictions
        )

    figure_paths = shared.save_evaluation_figures(
        output_dir,
        "multi_knn",
        confusion_matrices,
        confusion_matrix_64_classes,
        shared.classes_are_balanced(labels_train, labels_test),
    )

    # Every knob that affects the fitted model or the split it was fitted on.
    # ``n_neighbors`` is reported twice: the request is capped at the training
    # set size, so the effective value is what the estimator actually used.
    hyperparameters = {
        "n_neighbors_requested": int(config.n_neighbors),
        "n_neighbors_effective": int(effective_n_neighbors),
        "weights": config.weights,
        "algorithm": config.algorithm,
        "leaf_size": int(config.leaf_size),
        "metric": config.metric,
        "p": int(config.p),
        "n_jobs": config.n_jobs if config.n_jobs is not None else DEFAULT_N_JOBS,
        "seed": config.seed,
        "test_size": config.test_size,
        "cv_folds_requested": config.cv_folds,
        "cv_folds_effective": cv_folds,
        "rare_class_policy": config.rare_class_policy,
        "target_column": config.target_column,
    }

    metrics_payload = {
        "objective": "multi-label classification",
        "classifier": "multi_knn",
        "hyperparameters": hyperparameters,
        "metric_set": metric_set,
        "exact_match_is_diagnostic": True,
        "split_notes": split_notes,
        "cv": cv_results,
        "cv_notes": cv_warnings,
        "primary_metrics": {
            "hamming_loss": test_metrics["hamming_loss"],
            "micro_f1": test_metrics["micro_f1"],
            "macro_f1": test_metrics["macro_f1"],
            "weighted_f1": test_metrics["weighted_f1"],
            "bitwise_accuracy": test_metrics["bitwise_accuracy"],
        },
        "dataset_summary": dataset_summary,
        "timings_seconds": {
            "load": load_seconds,
            "split": split_seconds,
            "cross_validation": cv_seconds,
            "fit_and_predict": fit_predict_seconds,
            "feature_importance": feature_importance_seconds,
            "prediction_row_build": prediction_build_seconds,
        },
    }
    if include_diagnostic:
        metrics_payload["diagnostic_metrics"] = {
            "exact_match_accuracy": test_metrics["exact_match_accuracy"],
        }
        metrics_payload["classification_report"] = test_metrics["classification_report"]
        if confusion_matrices is not None:
            metrics_payload["confusion_matrices"] = confusion_matrices
        if confusion_matrix_64_classes is not None:
            metrics_payload["confusion_matrix_64_classes"] = confusion_matrix_64_classes
    if include_per_bit:
        metrics_payload["per_bit"] = test_metrics["per_bit"]
    if feature_rows is not None:
        metrics_payload["feature_importance"] = feature_rows
        metrics_payload["feature_importance_method"] = feature_importance_method
    metrics_payload.update(figure_paths)

    model_path = output_dir / "multi_knn_model.pkl"
    metrics_json_path = output_dir / "multi_knn_overall_metrics.json"
    metrics_txt_path = output_dir / "multi_knn_metrics.txt"
    feature_importance_path = output_dir / "feature_importances.tsv"
    predictions_path = output_dir / "predictions.tsv"

    artifact_start = time.perf_counter()
    with open(model_path, "wb") as handle:
        pickle.dump(model, handle)
    if feature_rows is not None:
        shared.write_tsv(feature_importance_path, feature_rows)
    else:
        feature_importance_path = None
    if prediction_rows is not None:
        shared.write_tsv(predictions_path, prediction_rows)
    else:
        predictions_path = None
    artifact_seconds = time.perf_counter() - artifact_start
    metrics_payload["timings_seconds"]["artifact_write"] = artifact_seconds
    metrics_payload["timings_seconds"]["total"] = time.perf_counter() - run_start

    text_lines = [
        "Ghostparser ML multi-label KNN baseline",
        "Primary objective: multi-label classification on the 6-bit classes bitstring.",
        "Exact-match accuracy is diagnostic; bitwise metrics are primary.",
        "",
    ]
    text_lines.extend(shared.format_hyperparameter_section(hyperparameters))
    text_lines.extend([
        "",
        "Test metrics:",
        f"  Hamming loss: {test_metrics['hamming_loss']:.6f}",
        f"  Bitwise accuracy: {test_metrics['bitwise_accuracy']:.6f}",
        f"  Exact-match accuracy: {test_metrics['exact_match_accuracy']:.6f}",
        f"  Micro F1: {test_metrics['micro_f1']:.6f}",
        f"  Macro F1: {test_metrics['macro_f1']:.6f}",
        f"  Weighted F1: {test_metrics['weighted_f1']:.6f}",
        "",
        "Per-bit metrics:",
    ])
    if include_per_bit:
        for bit_label, values in test_metrics["per_bit"].items():
            text_lines.append(
                (
                    f"  {bit_label}: accuracy={values['accuracy']:.6f}, "
                    f"precision={values['precision']:.6f}, "
                    f"recall={values['recall']:.6f}, f1={values['f1']:.6f}, "
                    f"support={values['support']}"
                )
            )
    if include_diagnostic:
        text_lines.extend(
            ["", "Classification report:", test_metrics["classification_report"]]
        )
        if confusion_matrices is not None:
            text_lines.extend(["", "Confusion matrices:"])
            text_lines.extend(
                shared.format_confusion_matrix_section(confusion_matrices)
            )
    text_lines.extend(
        shared.format_evaluation_figure_section(
            confusion_matrix_64_classes, figure_paths
        )
    )
    if report_class_distribution:
        text_lines.extend(["", "Dataset summary:"])
        text_lines.append("  Label map:")
        for key, value in dataset_summary["label_map"].items():
            text_lines.append(f"    {key}: {value}")
        text_lines.append("  Split:")
        for key, value in dataset_summary["split"].items():
            text_lines.append(f"    {key}: {value}")
    if cv_results is not None:
        text_lines.extend(["", f"Cross-validation folds: {cv_folds}", "CV aggregate:"])
        for key, value in cv_results["aggregate"].items():
            text_lines.append(f"  {key}: {value:.6f}")
    if cv_warnings:
        text_lines.extend(["", "Notes:"] + [f"  {note}" for note in cv_warnings])
    if split_notes:
        text_lines.extend(["", "Split notes:"] + [f"  {note}" for note in split_notes])
    text_lines.extend(
        [
            "",
            "Timings (seconds):",
            f"  load: {load_seconds:.6f}",
            f"  split: {split_seconds:.6f}",
            f"  cross_validation: {cv_seconds:.6f}",
            f"  fit_and_predict: {fit_predict_seconds:.6f}",
            f"  feature_importance: {feature_importance_seconds:.6f}",
            f"  prediction_row_build: {prediction_build_seconds:.6f}",
            f"  artifact_write: {artifact_seconds:.6f}",
            f"  total: {metrics_payload['timings_seconds']['total']:.6f}",
        ]
    )

    shared.write_json(metrics_json_path, metrics_payload)
    shared.write_text(metrics_txt_path, "\n".join(text_lines) + "\n")

    return {
        "model_path": str(model_path),
        "metrics_json_path": str(metrics_json_path),
        "metrics_txt_path": str(metrics_txt_path),
        "feature_importance_path": str(feature_importance_path)
        if feature_importance_path is not None
        else None,
        **{
            f"{key}_path": figure_paths.get(key)
            for key in shared.EVALUATION_FIGURES
        },
        "predictions_path": str(predictions_path)
        if predictions_path is not None
        else None,
        "metrics": metrics_payload,
    }


def _run(parsed_args) -> None:
    """Resolve the config, train, and print where the outputs went.

    Args:
        parsed_args: The parsed command line.
    """
    result = train_multi_knn(resolve_trainer_runtime_args(parsed_args))
    print(f"Saved model: {result['model_path']}")
    print(f"Saved metrics: {result['metrics_txt_path']}")
    print(f"Saved JSON metrics: {result['metrics_json_path']}")


def main() -> None:
    """Run the KNN trainer from the command line and exit with its status."""
    parser = build_trainer_argument_parser(
        "Ghostparser ML multi-label KNN baseline for summary statistics."
    )
    sys.exit(run_cli(parser, _run))


if __name__ == "__main__":
    main()
