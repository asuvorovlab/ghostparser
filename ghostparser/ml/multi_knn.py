"""Multi-label KNN baseline for Ghostparser summary statistics."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
from sklearn.inspection import permutation_importance
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.multioutput import MultiOutputClassifier
from sklearn.neighbors import KNeighborsClassifier

from ..cli_config import resolve_cli_or_config_args
from ..config import ConfigError
from . import ml_utils as shared
from .config import DEFAULT_N_JOBS, load_ml_config, normalize_ml_payload

BIT_LABELS = shared.BIT_LABELS
BIT_COUNT = shared.BIT_COUNT


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ghostparser ML multi-label KNN baseline for summary statistics."
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        default=None,
        help="Path to a JSON or YAML config file",
    )
    parser.add_argument(
        "-i",
        "--input-path",
        type=str,
        default=None,
        help="Path to summary_statistics.tsv",
    )
    parser.add_argument(
        "-o", "--output-dir", type=str, default=None, help="Directory for ML outputs"
    )
    return parser


def _resolve_runtime_args(args: argparse.Namespace) -> argparse.Namespace:
    return resolve_cli_or_config_args(
        args,
        load_config=load_ml_config,
        normalize_payload=normalize_ml_payload,
        payload_arg_names=["input_path", "output_dir"],
    )


def _build_model(
    config: argparse.Namespace, train_size: int
) -> tuple[MultiOutputClassifier, int]:
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


def _build_confusion_matrices(
    y_true: np.ndarray, y_pred: np.ndarray
) -> dict[str, list[list[int]]]:
    return {
        bit_label: confusion_matrix(
            y_true[:, bit_index], y_pred[:, bit_index], labels=[0, 1]
        ).tolist()
        for bit_index, bit_label in enumerate(BIT_LABELS)
    }


def _cross_validate(
    x_train: np.ndarray,
    y_train: np.ndarray,
    labels_train: np.ndarray,
    folds: int,
    random_state: int | None,
    config: argparse.Namespace,
) -> dict:
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
    rows = shared.read_tsv_rows(config.input_path)
    matrix = shared.rows_to_matrix(rows, config.target_column)
    labels = shared.combination_labels(matrix.train_labels)
    x_train, x_test, y_train, y_test, labels_train, labels_test, split_notes = (
        shared.split_dataset(
            matrix.train_features,
            matrix.train_targets,
            labels,
            config.test_size,
            config.random_state,
        )
    )

    model_factory = lambda train_size: _build_model(config, train_size)
    cv_folds, cv_warnings = shared.auto_cv_folds(
        labels_train, config.cv_folds, config.rare_class_policy
    )
    cv_results = {"folds": [], "aggregate": {}}
    if cv_folds is not None:
        cv_results = _cross_validate(
            x_train, y_train, labels_train, cv_folds, config.random_state, config
        )

    model, effective_n_neighbors = model_factory(len(x_train))
    model.fit(x_train, y_train)
    test_predictions = model.predict(x_test)
    test_metrics = shared.evaluate_predictions(y_test, test_predictions)

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metric_set = getattr(config, "evaluation_metrics", "all")
    report_class_distribution = getattr(config, "report_class_distribution", True)
    report_confusion_matrix = getattr(config, "report_confusion_matrix", True)
    report_feature_importance = getattr(config, "report_feature_importance", True)
    save_label_map = getattr(config, "save_label_map", True)
    save_predictions = getattr(config, "save_predictions", True)
    include_diagnostic = metric_set in {"diagnostic", "all"}
    include_per_bit = metric_set in {"per_bit", "all"}

    label_map = shared.build_label_map()
    class_distribution = None
    bit_distribution = None
    if report_class_distribution:
        class_distribution = {
            "overall": shared.summarize_distribution(matrix.train_labels),
            "train": shared.summarize_distribution(list(labels_train)),
            "test": shared.summarize_distribution(list(labels_test)),
        }
        bit_distribution = {
            "overall": shared.bit_distribution(matrix.train_targets),
            "train": shared.bit_distribution(y_train),
            "test": shared.bit_distribution(y_test),
        }

    feature_rows = None
    if report_feature_importance:
        permutation = permutation_importance(
            model,
            x_test,
            y_test,
            n_repeats=10,
            random_state=config.random_state,
            scoring="f1_micro",
        )
        feature_rows = shared.build_feature_importance_rows(
            matrix.feature_names, permutation.importances_mean
        )

    prediction_rows = None
    if save_predictions:
        prediction_rows = shared.build_prediction_rows(y_test, test_predictions)

    confusion_matrices = None
    if report_confusion_matrix and include_diagnostic:
        confusion_matrices = _build_confusion_matrices(y_test, test_predictions)

    metrics_payload = {
        "objective": "multi-label classification",
        "classifier": "multi_knn",
        "metric_set": metric_set,
        "exact_match_is_diagnostic": True,
        "split_notes": split_notes,
        "cv": cv_results,
        "cv_notes": cv_warnings,
        "knn": {
            "configured_n_neighbors": int(config.n_neighbors),
            "effective_n_neighbors": int(effective_n_neighbors),
            "weights": config.weights,
            "algorithm": config.algorithm,
            "leaf_size": int(config.leaf_size),
            "metric": config.metric,
            "p": int(config.p),
        },
    }

    metrics_payload["primary_metrics"] = {
        "hamming_loss": test_metrics["hamming_loss"],
        "micro_f1": test_metrics["micro_f1"],
        "macro_f1": test_metrics["macro_f1"],
        "weighted_f1": test_metrics["weighted_f1"],
        "bitwise_accuracy": test_metrics["bitwise_accuracy"],
    }
    if include_diagnostic:
        metrics_payload["diagnostic_metrics"] = {
            "exact_match_accuracy": test_metrics["exact_match_accuracy"],
        }
        metrics_payload["classification_report"] = test_metrics["classification_report"]
        if confusion_matrices is not None:
            metrics_payload["confusion_matrices"] = confusion_matrices
    if include_per_bit:
        metrics_payload["per_bit"] = test_metrics["per_bit"]
    if class_distribution is not None and bit_distribution is not None:
        metrics_payload["class_distribution"] = class_distribution
        metrics_payload["bit_distribution"] = bit_distribution
    if feature_rows is not None:
        metrics_payload["feature_importance"] = feature_rows

    model_path = output_dir / "multi_knn_model.pkl"
    metrics_json_path = output_dir / "multi_knn_metrics.json"
    metrics_txt_path = output_dir / "multi_knn_metrics.txt"
    label_map_path = output_dir / "label_map.json"
    class_distribution_path = output_dir / "class_distribution.tsv"
    bit_distribution_path = output_dir / "bit_distribution.tsv"
    feature_importance_path = output_dir / "feature_importances.tsv"
    predictions_path = output_dir / "predictions.tsv"

    with open(model_path, "wb") as handle:
        pickle.dump(model, handle)
    shared.write_json(metrics_json_path, metrics_payload)
    if save_label_map:
        shared.write_json(label_map_path, label_map)
    else:
        label_map_path = None
    if class_distribution is not None and bit_distribution is not None:
        shared.write_tsv(
            class_distribution_path,
            [
                {
                    "partition": partition,
                    "label": label,
                    "count": values["count"],
                    "fraction": values["fraction"],
                }
                for partition, distribution in class_distribution.items()
                for label, values in distribution.items()
            ],
        )
        shared.write_tsv(
            bit_distribution_path,
            [
                {
                    "partition": partition,
                    "bit_label": bit_label,
                    "positive_count": values["positive_count"],
                    "fraction": values["fraction"],
                }
                for partition, distribution in bit_distribution.items()
                for bit_label, values in distribution.items()
            ],
        )
    else:
        class_distribution_path = None
        bit_distribution_path = None
    if feature_rows is not None:
        shared.write_tsv(feature_importance_path, feature_rows)
    else:
        feature_importance_path = None
    if prediction_rows is not None:
        shared.write_tsv(predictions_path, prediction_rows)
    else:
        predictions_path = None

    text_lines = [
        "Ghostparser ML multi-label KNN baseline",
        "Primary objective: multi-label classification on the 6-bit classes bitstring.",
        "Exact-match accuracy is diagnostic; bitwise metrics are primary.",
        f"Configured n_neighbors: {int(config.n_neighbors)}",
        f"Effective n_neighbors: {int(effective_n_neighbors)}",
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
    ]
    if include_per_bit:
        for bit_label, values in test_metrics["per_bit"].items():
            text_lines.append(
                f"  {bit_label}: accuracy={values['accuracy']:.6f}, precision={values['precision']:.6f}, recall={values['recall']:.6f}, f1={values['f1']:.6f}, support={values['support']}"
            )
    if include_diagnostic:
        text_lines.extend(
            ["", "Classification report:", test_metrics["classification_report"]]
        )
        if confusion_matrices is not None:
            text_lines.extend(["", "Confusion matrices:"])
            for bit_label, matrix_values in confusion_matrices.items():
                text_lines.append(f"  {bit_label}: {matrix_values}")
    if class_distribution is not None and bit_distribution is not None:
        text_lines.extend(["", "Class distribution (bitstrings):"])
        for partition, distribution in class_distribution.items():
            text_lines.append(f"  {partition}:")
            for label, values in distribution.items():
                text_lines.append(
                    f"    {label}: count={values['count']}, fraction={values['fraction']:.6f}"
                )
    if cv_results is not None:
        text_lines.extend(["", f"Cross-validation folds: {cv_folds}", "CV aggregate:"])
        for key, value in cv_results["aggregate"].items():
            text_lines.append(f"  {key}: {value:.6f}")
    if cv_warnings:
        text_lines.extend(["", "Notes:"] + [f"  {note}" for note in cv_warnings])
    if split_notes:
        text_lines.extend(["", "Split notes:"] + [f"  {note}" for note in split_notes])
    shared.write_text(metrics_txt_path, "\n".join(text_lines) + "\n")

    return {
        "model_path": str(model_path),
        "metrics_json_path": str(metrics_json_path),
        "metrics_txt_path": str(metrics_txt_path),
        "label_map_path": str(label_map_path) if label_map_path is not None else None,
        "class_distribution_path": str(class_distribution_path)
        if class_distribution_path is not None
        else None,
        "bit_distribution_path": str(bit_distribution_path)
        if bit_distribution_path is not None
        else None,
        "feature_importance_path": str(feature_importance_path)
        if feature_importance_path is not None
        else None,
        "predictions_path": str(predictions_path)
        if predictions_path is not None
        else None,
        "metrics": metrics_payload,
    }


def main() -> None:
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except (ValueError, ConfigError) as exc:
        print(f"Error: {exc}")
        return

    try:
        result = train_multi_knn(args)
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}")
        return

    print(f"Saved model: {result['model_path']}")
    print(f"Saved metrics: {result['metrics_txt_path']}")
    print(f"Saved JSON metrics: {result['metrics_json_path']}")


if __name__ == "__main__":
    main()
