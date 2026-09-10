"""Shared helpers for Ghostparser machine-learning baselines."""

import csv
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    hamming_loss,
    precision_recall_fscore_support,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

BIT_LABELS = (
    "ghost_into_A",
    "ghost_into_B",
    "inflow_into_A_from_C",
    "inflow_into_B_from_C",
    "outflow_from_A_to_C",
    "outflow_from_B_to_C",
)
BIT_COUNT = len(BIT_LABELS)
MAX_STRING_CATEGORIES = 7

# All 64 class labels are drawn on each axis of the 64-class figure, so the tick
# font has to fit a 6-character label into one cell of the square grid.
_CLASS_TICK_FONT_SIZE = 6

# Colormap for both confusion-matrix figures. Both plot a sequential quantity,
# and cividis is perceptually uniform and colour-vision-deficiency safe.
CONFUSION_MATRIX_COLORMAP = "cividis"


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
            raise ValueError(
                f"Invalid classes label: {label!r}; expected a {BIT_COUNT}-character 0/1 bitstring"
            )
        labels.append(label)
        binary_rows.append([int(bit) for bit in label])
    return np.asarray(binary_rows, dtype=int), labels


def select_feature_names(fieldnames: list[str], target_column: str) -> tuple[str, ...]:
    return tuple(name for name in fieldnames if name != target_column)


def _parse_numeric_column(values: list[str]) -> np.ndarray | None:
    try:
        return np.asarray([float(value) for value in values], dtype=float)
    except (TypeError, ValueError):
        return None


def _encode_feature_column(
    values: list[str], feature_name: str
) -> tuple[np.ndarray, tuple[str, ...]]:
    parsed_numeric = _parse_numeric_column(values)
    if parsed_numeric is not None:
        return parsed_numeric.reshape(-1, 1), (feature_name,)

    categories = sorted({str(value).strip() for value in values})
    if len(categories) > MAX_STRING_CATEGORIES:
        raise ValueError(
            f"Non-numeric feature value for {feature_name!r}: {values[0]!r}; "
            f"string-valued columns must have at most {MAX_STRING_CATEGORIES} distinct values"
        )

    encoder = OneHotEncoder(
        categories=[categories],
        sparse_output=False,
        handle_unknown="ignore",
        dtype=float,
    )
    encoded = encoder.fit_transform(np.asarray(values, dtype=object).reshape(-1, 1))
    encoded_feature_names = tuple(
        f"{feature_name}={category}" for category in encoder.categories_[0]
    )
    return encoded, encoded_feature_names


def rows_to_matrix(rows: list[dict[str, str]], target_column: str) -> DatasetSplit:
    fieldnames = list(rows[0].keys())
    if target_column not in fieldnames:
        raise ValueError(f"Missing required target column: {target_column}")

    feature_names = select_feature_names(fieldnames, target_column)
    if not feature_names:
        raise ValueError("No feature columns found after excluding the target column")

    raw_labels: list[str] = []
    encoded_blocks: list[np.ndarray] = []
    encoded_feature_names: list[str] = []
    column_values = {feature_name: [] for feature_name in feature_names}
    for row in rows:
        raw_labels.append(row[target_column])
        for feature_name in feature_names:
            raw_value = row.get(feature_name, "")
            if raw_value is None or str(raw_value).strip() == "":
                raise ValueError(f"Missing feature value for {feature_name!r}")
            column_values[feature_name].append(str(raw_value).strip())

    for feature_name in feature_names:
        encoded_block, block_feature_names = _encode_feature_column(
            column_values[feature_name], feature_name
        )
        encoded_blocks.append(encoded_block)
        encoded_feature_names.extend(block_feature_names)

    targets, labels = parse_classes(raw_labels)
    features = np.hstack(encoded_blocks)

    return DatasetSplit(
        feature_names=tuple(encoded_feature_names),
        train_features=features,
        test_features=np.empty((0, features.shape[1])),
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
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]
]:
    try:
        x_train, x_test, y_train, y_test, labels_train, labels_test = train_test_split(
            features,
            targets,
            labels,
            test_size=test_size,
            random_state=random_state,
            stratify=labels,
        )
        return (
            x_train,
            x_test,
            y_train,
            y_test,
            labels_train,
            labels_test,
            ["stratified"],
        )
    except ValueError:
        x_train, x_test, y_train, y_test, labels_train, labels_test = train_test_split(
            features,
            targets,
            labels,
            test_size=test_size,
            random_state=random_state,
            stratify=None,
        )
        return (
            x_train,
            x_test,
            y_train,
            y_test,
            labels_train,
            labels_test,
            ["unstratified"],
        )


def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    per_bit_precision, per_bit_recall, per_bit_f1, per_bit_support = (
        precision_recall_fscore_support(y_true, y_pred, average=None, zero_division=0)
    )
    per_bit_accuracy = (y_true == y_pred).mean(axis=0)
    report = classification_report(
        y_true, y_pred, target_names=BIT_LABELS, zero_division=0
    )
    return {
        "exact_match_accuracy": float(accuracy_score(y_true, y_pred)),
        "hamming_loss": float(hamming_loss(y_true, y_pred)),
        "bitwise_accuracy": float((y_true == y_pred).mean()),
        "micro_f1": float(f1_score(y_true, y_pred, average="micro", zero_division=0)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "weighted_f1": float(
            f1_score(y_true, y_pred, average="weighted", zero_division=0)
        ),
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
        bit_label: {
            "positive_count": int(totals[index]),
            "fraction": float(totals[index] / total_rows),
        }
        for index, bit_label in enumerate(BIT_LABELS)
    }


def auto_cv_folds(
    labels: np.ndarray, requested_folds: int, policy: str
) -> tuple[int | None, list[str]]:
    counts = Counter(labels.tolist())
    if not counts:
        return None, ["No labels available for cross-validation"]

    min_count = min(counts.values())
    warnings: list[str] = []
    if min_count < 2:
        if policy == "error":
            raise ValueError(
                "Cannot run stratified cross-validation because at least one class has fewer than 2 samples"
            )
        return None, [
            "Skipped cross-validation because at least one class has fewer than 2 samples"
        ]

    folds = min(requested_folds, min_count)
    if folds < requested_folds:
        warnings.append(
            f"Reduced CV folds from {requested_folds} to {folds} because the smallest class has {min_count} samples"
        )

    if policy == "warn_skip_cv" and folds < requested_folds:
        return None, warnings + [
            "Skipped cross-validation because the requested fold count was not feasible"
        ]

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


def format_bit_label_title(bit_label: str) -> str:
    """Render a bit label as a plot title.

    Underscores become spaces and the first character is upper-cased. Only the
    first -- ``str.capitalize`` would lower-case the rest and destroy the taxon
    letters, turning ``ghost_into_A`` into ``Ghost into a``.

    Args:
        bit_label: One of :data:`BIT_LABELS`, e.g. ``inflow_into_A_from_C``.

    Returns:
        The title form, e.g. ``Inflow into A from C``.
    """
    spaced = bit_label.replace("_", " ")
    return spaced[:1].upper() + spaced[1:]


def build_confusion_matrices(
    y_true: np.ndarray, y_pred: np.ndarray
) -> dict[str, list[list[int]]]:
    """Build a 2x2 confusion matrix per label bit.

    Args:
        y_true: True multi-label targets, one column per bit.
        y_pred: Predicted multi-label targets, aligned with ``y_true``.

    Returns:
        A mapping of bit label to its ``[[tn, fp], [fn, tp]]`` matrix.
    """
    return {
        bit_label: confusion_matrix(
            y_true[:, bit_index], y_pred[:, bit_index], labels=[0, 1]
        ).tolist()
        for bit_index, bit_label in enumerate(BIT_LABELS)
    }


def build_label_map() -> dict:
    return {
        "source_label_format": "6-bit bitstring",
        "bit_labels": list(BIT_LABELS),
        "bit_index_order": list(range(BIT_COUNT)),
    }


def build_dataset_summary(
    all_labels: list[str],
    train_labels: list[str],
    test_labels: list[str],
    train_targets: np.ndarray,
    test_targets: np.ndarray,
    split_notes: list[str],
    test_size: float,
    random_state: int | None,
) -> dict:
    return {
        "label_map": build_label_map(),
        "split": {
            "strategy": split_notes[0] if split_notes else None,
            "notes": list(split_notes),
            "test_size": float(test_size),
            "random_state": random_state,
            "total_rows": int(len(all_labels)),
            "train_rows": int(len(train_labels)),
            "test_rows": int(len(test_labels)),
        },
        "class_distribution": {
            "overall": summarize_distribution(all_labels),
            "train": summarize_distribution(train_labels),
            "test": summarize_distribution(test_labels),
        },
        "bit_distribution": {
            "overall": bit_distribution(
                train_targets
                if len(test_targets) == 0
                else np.vstack([train_targets, test_targets])
            ),
            "train": bit_distribution(train_targets),
            "test": bit_distribution(test_targets),
        },
    }


def build_prediction_rows(
    y_true: np.ndarray, y_pred: np.ndarray
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index, (true_row, pred_row) in enumerate(zip(y_true, y_pred)):
        matched_label_count = int(np.sum(true_row == pred_row))
        row: dict[str, object] = {
            "row_index": int(index),
            "true_label": "".join(str(int(value)) for value in true_row),
            "pred_label": "".join(str(int(value)) for value in pred_row),
            "exact_match": int(np.array_equal(true_row, pred_row)),
            "matched_label_count": matched_label_count,
        }
        for bit_index, bit_label in enumerate(BIT_LABELS):
            row[f"true_{bit_label}"] = int(true_row[bit_index])
            row[f"pred_{bit_label}"] = int(pred_row[bit_index])
        rows.append(row)
    return rows


def build_feature_importance_rows(
    feature_names: tuple[str, ...], scores: np.ndarray
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index, feature_name in enumerate(feature_names):
        rows.append({"feature": feature_name, "importance": float(scores[index])})
    rows.sort(key=lambda row: row["importance"], reverse=True)
    return rows


def format_hyperparameter_section(
    hyperparameters: dict[str, object],
    title: str = "Hyperparameters:",
) -> list[str]:
    """Render a trainer's hyperparameters as aligned ``key: value`` lines.

    Args:
        hyperparameters: Mapping of hyperparameter name to the value actually
            used for the run. Insertion order is preserved in the output.
        title: Heading placed above the block.

    Returns:
        A list of text lines: the title followed by one indented line per
        hyperparameter, with names padded to a common width. ``None`` values
        render as ``none`` so an unset knob is visibly distinct from an empty
        string.
    """
    lines = [title]
    if not hyperparameters:
        lines.append("  (none)")
        return lines
    width = max(len(str(name)) for name in hyperparameters)
    for name, value in hyperparameters.items():
        rendered = "none" if value is None else str(value)
        lines.append(f"  {str(name):<{width}}  {rendered}")
    return lines


def format_confusion_matrix_section(
    confusion_matrices: dict[str, list[list[int]]],
) -> list[str]:
    lines: list[str] = []
    for bit_label, matrix_values in confusion_matrices.items():
        matrix_array = np.asarray(matrix_values, dtype=float)
        total = float(matrix_array.sum())

        def _cell_text(row_index: int, col_index: int) -> str:
            count = int(matrix_values[row_index][col_index])
            pct = (count / total * 100.0) if total > 0 else 0.0
            return f"{count} ({pct:.1f}%)"

        lines.append(f"  {bit_label}:")
        lines.append("           pred=0  pred=1")
        lines.append(f"    true=0  {_cell_text(0, 0):>11}  {_cell_text(0, 1):>11}")
        lines.append(f"    true=1  {_cell_text(1, 0):>11}  {_cell_text(1, 1):>11}")
    return lines


def save_confusion_matrix_plot(
    confusion_matrices: dict[str, list[list[int]]],
    output_path: Path,
) -> str | None:
    if not confusion_matrices:
        return None

    matrix_items = list(confusion_matrices.items())
    max_value = max(
        int(np.max(np.asarray(matrix_values, dtype=int)))
        for _, matrix_values in matrix_items
    )
    max_value = max(1, max_value)
    n_plots = len(matrix_items)
    n_cols = min(3, n_plots)
    n_rows = math.ceil(n_plots / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5.2 * n_cols, 4.4 * n_rows),
        constrained_layout=True,
    )
    axes_array = np.atleast_1d(axes).ravel()
    cmap = plt.get_cmap(CONFUSION_MATRIX_COLORMAP)

    for axis_index, (bit_label, matrix_values) in enumerate(matrix_items):
        ax = axes_array[axis_index]
        data = np.asarray(matrix_values, dtype=int)
        total = float(data.sum())
        if total > 0:
            label_grid = np.asarray(
                [
                    [f"{value}\n({value / total * 100:.1f}%)" for value in row]
                    for row in data
                ],
                dtype=object,
            )
        else:
            label_grid = np.asarray(
                [[f"{value}\n(0.0%)" for value in row] for row in data],
                dtype=object,
            )
        sns.heatmap(
            data,
            ax=ax,
            cmap=cmap,
            vmin=0,
            vmax=max_value,
            annot=label_grid,
            fmt="",
            square=True,
            cbar=False,
            linewidths=1,
            linecolor="white",
            annot_kws={"size": 11, "weight": "bold"},
        )
        ax.set_title(format_bit_label_title(bit_label), fontsize=12)
        ax.set_xlabel("Predicted label (0 = predicted zero, 1 = predicted one)")
        ax.set_ylabel("True label (0 = true zero, 1 = true one)")
        ax.set_xticklabels(["0", "1"], rotation=0)
        ax.set_yticklabels(["0", "1"], rotation=0)

    for axis_index in range(n_plots, len(axes_array)):
        axes_array[axis_index].axis("off")

    fig.suptitle("Confusion matrices by bit", fontsize=15)
    colorbar_mappable = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_value)
    )
    colorbar_mappable.set_array([])
    fig.colorbar(
        colorbar_mappable,
        ax=axes_array[:n_plots].tolist(),
        shrink=0.85,
        pad=0.02,
        label="Count",
    )
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return str(output_path)


def build_64_class_confusion_matrix(
    y_true: np.ndarray,
    y_pred: np.ndarray,
) -> dict[str, object]:
    powers = (2 ** np.arange(BIT_COUNT - 1, -1, -1)).astype(int)
    true_indices = (np.asarray(y_true, dtype=int) * powers).sum(axis=1)
    pred_indices = (np.asarray(y_pred, dtype=int) * powers).sum(axis=1)

    matrix = np.zeros((2**BIT_COUNT, 2**BIT_COUNT), dtype=int)
    for true_idx, pred_idx in zip(true_indices, pred_indices):
        matrix[int(true_idx), int(pred_idx)] += 1

    class_labels = [format(index, f"0{BIT_COUNT}b") for index in range(2**BIT_COUNT)]
    return {
        "class_labels": class_labels,
        "matrix": matrix.tolist(),
    }


def row_normalize_confusion_matrix(matrix: object) -> np.ndarray:
    """Convert confusion-matrix counts into per-true-class fractions.

    Each cell becomes the share of its true class that landed in that predicted
    column, so every populated row sums to 1 and the diagonal reads as recall.

    Args:
        matrix: Square array-like of confusion-matrix counts.

    Returns:
        A float array on ``[0, 1]`` with the same shape. Rows whose true class
        has no samples stay all-zero instead of dividing by zero.
    """
    counts = np.asarray(matrix, dtype=float)
    row_totals = counts.sum(axis=1, keepdims=True)
    return np.divide(
        counts, row_totals, out=np.zeros_like(counts), where=row_totals > 0
    )


def save_64_class_confusion_matrix_plot(
    class_confusion: dict[str, object],
    output_path: Path,
) -> str:
    labels = [str(label) for label in class_confusion["class_labels"]]
    # Raw counts depend on how many test rows each class happened to draw, which
    # is fixed per run but arbitrary to a reader. Row-normalizing puts every cell
    # on a common 0-1 scale that is comparable across rows and across runs.
    data = row_normalize_confusion_matrix(class_confusion["matrix"])

    fig, ax = plt.subplots(figsize=(18, 16), constrained_layout=True)
    cmap = plt.get_cmap(CONFUSION_MATRIX_COLORMAP)
    # Zero cells are painted at the colormap's low end like any other value, so
    # the grid reads as one continuous surface and the colour scale covers every
    # cell in it.
    sns.heatmap(
        data,
        ax=ax,
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        square=True,
        cbar=True,
        cbar_kws={"label": "Fraction of true class"},
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
    )
    # Every class is labelled on both axes: a reader looking up one specific
    # 6-bit class cannot count rows inwards from a subsampled tick.
    tick_positions = np.arange(len(labels)) + 0.5
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(labels, rotation=90, fontsize=_CLASS_TICK_FONT_SIZE)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(labels, rotation=0, fontsize=_CLASS_TICK_FONT_SIZE)
    ax.tick_params(axis="both", length=2, pad=1.5)
    ax.set_xlabel("Predicted 6 binary class")
    ax.set_ylabel("True 6 binary class")
    ax.set_title(
        "Confusion matrix across all 64 possible 6-bit classes "
        "(row-normalized: fraction of each true class)"
    )

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return str(output_path)
