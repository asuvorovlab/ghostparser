"""Shared helpers for Ghostparser machine-learning baselines."""

from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from sklearn.metrics import accuracy_score, classification_report, f1_score, hamming_loss, precision_recall_fscore_support
from sklearn.model_selection import train_test_split

BIT_LABELS = (
    "ghost_into_A",
    "ghost_into_B",
    "inflow_into_A_from_C",
    "inflow_into_B_from_C",
    "outflow_from_A_to_C",
    "outflow_from_B_to_C",
)
BIT_COUNT = len(BIT_LABELS)
DEFAULT_EXCLUDED_COLUMNS = {
    "triplet",
    "abc_mapping",
    "species_tree",
    "classification",
    "bootstrap_value",
    "source_folder",
    "classes",
}


@dataclass(frozen=True)
class DatasetSplit:
    """Prepared train/test arrays for ML training."""

    feature_names: tuple[str, ...]
    train_features: np.ndarray
    test_features: np.ndarray
    train_targets: np.ndarray
    test_targets: np.ndarray
    train_labels: list[str]
    test_labels: list[str]


def read_tsv_rows(input_path: str) -> list[dict[str, str]]:
    with open(input_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV is missing a header row")
        rows = list(reader)
    if not rows:
        raise ValueError("Input TSV contains no data rows")
    return rows


def is_valid_bitstring(value: str) -> bool:
    return len(value) == BIT_COUNT and set(value) <= {"0", "1"}


def parse_classes(raw_labels: Iterable[str]) -> tuple[np.ndarray, list[str]]:
    labels: list[str] = []
    binary_rows: list[list[int]] = []
    for raw_label in raw_labels:
        label = str(raw_label).strip()
        if not is_valid_bitstring(label):
            raise ValueError(f"Invalid classes label: {label!r}; expected a {BIT_COUNT}-character 0/1 bitstring")
        labels.append(label)
        binary_rows.append([int(bit) for bit in label])
    return np.asarray(binary_rows, dtype=int), labels


def select_feature_names(fieldnames: list[str]) -> tuple[str, ...]:
    excluded = set(DEFAULT_EXCLUDED_COLUMNS)
    return tuple(name for name in fieldnames if name not in excluded)


def rows_to_matrix(rows: list[dict[str, str]], target_column: str) -> DatasetSplit:
    fieldnames = list(rows[0].keys())
    if target_column not in fieldnames:
        raise ValueError(f"Missing required target column: {target_column}")

    feature_names = select_feature_names(fieldnames)
    if not feature_names:
        raise ValueError("No feature columns found after excluding metadata columns")

    feature_rows: list[list[float]] = []
    raw_labels: list[str] = []
    for row in rows:
        raw_labels.append(row[target_column])
        feature_row: list[float] = []
        for feature_name in feature_names:
            raw_value = row.get(feature_name, "")
            if raw_value is None or str(raw_value).strip() == "":
                raise ValueError(f"Missing feature value for {feature_name!r}")
            try:
                feature_row.append(float(raw_value))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Non-numeric feature value for {feature_name!r}: {raw_value!r}") from exc
        feature_rows.append(feature_row)

    targets, labels = parse_classes(raw_labels)
    features = np.asarray(feature_rows, dtype=float)

    return DatasetSplit(
        feature_names=feature_names,
        train_features=features,
        test_features=np.empty((0, len(feature_names))),
        train_targets=targets,
        test_targets=np.empty((0, BIT_COUNT), dtype=int),
        train_labels=labels,
        test_labels=[],
    )


def combination_labels(bit_labels: Iterable[str]) -> np.ndarray:
    return np.asarray(list(bit_labels))


def split_dataset(
    features: np.ndarray,
    targets: np.ndarray,
    labels: np.ndarray,
    test_size: float,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    try:
        x_train, x_test, y_train, y_test, labels_train, labels_test = train_test_split(
            features,
            targets,
            labels,
            test_size=test_size,
            random_state=random_state,
            stratify=labels,
        )
        return x_train, x_test, y_train, y_test, labels_train, labels_test, ["stratified"]
    except ValueError:
        x_train, x_test, y_train, y_test, labels_train, labels_test = train_test_split(
            features,
            targets,
            labels,
            test_size=test_size,
            random_state=random_state,
            stratify=None,
        )
        return x_train, x_test, y_train, y_test, labels_train, labels_test, ["unstratified"]


def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    per_bit_precision, per_bit_recall, per_bit_f1, per_bit_support = precision_recall_fscore_support(
        y_true, y_pred, average=None, zero_division=0
    )
    per_bit_accuracy = (y_true == y_pred).mean(axis=0)
    report = classification_report(y_true, y_pred, target_names=BIT_LABELS, zero_division=0)
    return {
        "exact_match_accuracy": float(accuracy_score(y_true, y_pred)),
        "hamming_loss": float(hamming_loss(y_true, y_pred)),
        "bitwise_accuracy": float((y_true == y_pred).mean()),
        "micro_f1": float(f1_score(y_true, y_pred, average="micro", zero_division=0)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, average="weighted", zero_division=0)),
        "per_bit": {
            bit_label: {
                "precision": float(per_bit_precision[index]),
                "recall": float(per_bit_recall[index]),
                "f1": float(per_bit_f1[index]),
                "support": int(per_bit_support[index]),
                "accuracy": float(per_bit_accuracy[index]),
            }
            for index, bit_label in enumerate(BIT_LABELS)
        },
        "classification_report": report,
    }


def summarize_distribution(labels: list[str]) -> dict[str, dict[str, float]]:
    counts = Counter(labels)
    total = sum(counts.values()) or 1
    return {
        label: {"count": int(count), "fraction": float(count / total)}
        for label, count in sorted(counts.items())
    }


def bit_distribution(targets: np.ndarray) -> dict[str, dict[str, float]]:
    totals = targets.sum(axis=0)
    total_rows = targets.shape[0] or 1
    return {
        bit_label: {"positive_count": int(totals[index]), "fraction": float(totals[index] / total_rows)}
        for index, bit_label in enumerate(BIT_LABELS)
    }


def auto_cv_folds(labels: np.ndarray, requested_folds: int, policy: str) -> tuple[int | None, list[str]]:
    counts = Counter(labels.tolist())
    if not counts:
        return None, ["No labels available for cross-validation"]

    min_count = min(counts.values())
    warnings: list[str] = []
    if min_count < 2:
        if policy == "error":
            raise ValueError("Cannot run stratified cross-validation because at least one class has fewer than 2 samples")
        return None, ["Skipped cross-validation because at least one class has fewer than 2 samples"]

    folds = min(requested_folds, min_count)
    if folds < requested_folds:
        warnings.append(f"Reduced CV folds from {requested_folds} to {folds} because the smallest class has {min_count} samples")

    if policy == "warn_skip_cv" and folds < requested_folds:
        return None, warnings + ["Skipped cross-validation because the requested fold count was not feasible"]

    return folds, warnings


def write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_label_map() -> dict:
    return {
        "source_label_format": "6-bit bitstring",
        "bit_labels": list(BIT_LABELS),
        "bit_index_order": list(range(BIT_COUNT)),
    }


def build_prediction_rows(y_true: np.ndarray, y_pred: np.ndarray) -> list[dict[str, object]]:
    rows = []
    for index, (true_row, pred_row) in enumerate(zip(y_true, y_pred)):
        rows.append(
            {
                "row_index": index,
                "true_label": "".join(str(int(value)) for value in true_row),
                "pred_label": "".join(str(int(value)) for value in pred_row),
                "exact_match": int(np.array_equal(true_row, pred_row)),
                **{f"true_{bit_label}": int(true_row[bit_index]) for bit_index, bit_label in enumerate(BIT_LABELS)},
                **{f"pred_{bit_label}": int(pred_row[bit_index]) for bit_index, bit_label in enumerate(BIT_LABELS)},
            }
        )
    return rows


def build_feature_importance_rows(feature_names: tuple[str, ...], scores: np.ndarray) -> list[dict[str, object]]:
    rows = [
        {"feature": feature_name, "importance": float(scores[index])}
        for index, feature_name in enumerate(feature_names)
    ]
    rows.sort(key=lambda row: row["importance"], reverse=True)
    return rows