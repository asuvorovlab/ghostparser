import numpy as np
import pytest

from ghostparser.ml.ml_utils import (
    row_normalize_confusion_matrix,
    rows_to_matrix,
)


def test_rows_to_matrix_uses_numeric_features_and_excludes_target_column():
    rows = [
        {
            "class": "001010",
            "dis1_topology": "BC",
            "feature_a": "10.0",
            "feature_b": "1.5",
            "feature_x": "1.5",
        },
        {
            "class": "001011",
            "dis1_topology": "AC",
            "feature_a": "20.0",
            "feature_b": "2.5",
            "feature_x": "2.5",
        },
    ]

    matrix = rows_to_matrix(rows, target_column="class")

    assert matrix.feature_names == (
        "dis1_topology=AC",
        "dis1_topology=BC",
        "feature_a",
        "feature_b",
        "feature_x",
    )
    assert matrix.train_features.shape == (2, 5)
    assert matrix.train_features[0].tolist() == [0.0, 1.0, 10.0, 1.5, 1.5]
    assert matrix.train_features[1].tolist() == [1.0, 0.0, 20.0, 2.5, 2.5]
    assert matrix.train_labels == ["001010", "001011"]


def test_rows_to_matrix_encodes_multiple_string_columns():
    rows = [
        {
            "class": "001010",
            "topology_a": "BC",
            "topology_b": "XY",
            "feature_x": "1.0",
        },
        {
            "class": "001011",
            "topology_a": "AC",
            "topology_b": "XZ",
            "feature_x": "2.0",
        },
    ]

    matrix = rows_to_matrix(rows, target_column="class")

    assert matrix.feature_names == (
        "topology_a=AC",
        "topology_a=BC",
        "topology_b=XY",
        "topology_b=XZ",
        "feature_x",
    )
    assert matrix.train_features.shape == (2, 5)
    assert matrix.train_features[0].tolist() == [0.0, 1.0, 1.0, 0.0, 1.0]
    assert matrix.train_features[1].tolist() == [1.0, 0.0, 0.0, 1.0, 2.0]


def test_rows_to_matrix_rejects_string_features():
    rows = [
        {
            "class": "001010",
            "dis1_topology": value,
            "feature_x": "1.0",
        }
        for i, value in enumerate(["A", "B", "C", "D", "E", "F", "G", "H"], start=1)
    ]

    with pytest.raises(
        ValueError,
        match="string-valued columns must have at most 7 distinct values",
    ):
        rows_to_matrix(rows, target_column="class")


def test_row_normalize_confusion_matrix_turns_counts_into_per_class_fractions():
    """Every populated row sums to 1; a true class with no samples stays zero."""
    counts = [
        [3, 1, 0],
        [0, 0, 0],
        [1, 1, 2],
    ]

    fractions = row_normalize_confusion_matrix(counts)

    np.testing.assert_allclose(
        fractions,
        [
            [0.75, 0.25, 0.0],
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.5],
        ],
    )
    assert fractions.min() >= 0.0
    assert fractions.max() <= 1.0
    # Rows with samples normalize to 1; the empty row is left at 0 so the plot
    # masks it rather than dividing by zero.
    np.testing.assert_allclose(fractions.sum(axis=1), [1.0, 0.0, 1.0])
