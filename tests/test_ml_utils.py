import numpy as np
import pytest

from ghostparser.config import ConfigError
from ghostparser.ml.ml_utils import (
    BIT_COUNT,
    correlation_feature_groups,
    feature_importance_rows,
    grouped_permutation_importance,
    row_normalize_confusion_matrix,
    rows_to_matrix,
)


def test_rows_to_matrix_uses_numeric_features_and_excludes_target_column():
    """Numeric columns pass through, the target is dropped, strings are one-hot.

    Feature order is asserted because the model pickle is reused across runs:
    one-hot columns come first in sorted-category order, then the numerics in
    input order, and a reader scoring new data must reproduce that layout.
    """
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
    """Each string column is one-hot encoded on its own categories.

    Two string columns must not share a category space: the names stay
    prefixed by their column, and the blocks stay grouped per column rather
    than interleaved.
    """
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
    """A high-cardinality string column is refused rather than one-hot exploded.

    Eight distinct values exceed the 7-category cap, which exists because a
    free-text column would otherwise silently add a column per value.
    """
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


def test_correlation_feature_groups_collects_the_redundant_columns():
    """Features that rank the rows alike land in one group, others stay apart.

    Column 1 is a monotone transform of column 0, so their Spearman
    correlation is exactly 1 whatever the values; column 2 ranks the rows in
    an unrelated order. The threshold is the default 0.7.
    """
    features = np.array(
        [
            [1.0, 1.0, 3.0],
            [2.0, 4.0, 1.0],
            [3.0, 9.0, 4.0],
            [4.0, 16.0, 2.0],
            [5.0, 25.0, 5.0],
        ]
    )

    groups = correlation_feature_groups(
        features, ("height_mean", "height_square", "unrelated"), 0.7
    )

    assert groups == ((0, 1), (2,))


def test_grouped_permutation_importance_scores_the_group_that_carries_the_label():
    """A group is scored by the micro-F1 lost when the whole group is shuffled.

    The stand-in model reads bit 0 off column 0 alone, so it predicts the
    labels exactly: shuffling the group that holds column 0 breaks the
    predictions, while shuffling the group of unread columns leaves every
    prediction and therefore the score untouched, for an importance of
    exactly zero.
    """

    class ColumnZeroModel:
        def predict(self, features):
            positive = (features[:, 0] > 0.5).astype(int)
            return np.column_stack([positive] * BIT_COUNT)

    rows = 40
    column_zero = np.tile([0.0, 1.0], rows // 2)
    features = np.column_stack(
        [column_zero, column_zero * 2.0, np.arange(rows, dtype=float)]
    )
    targets = ColumnZeroModel().predict(features)

    means, deviations = grouped_permutation_importance(
        ColumnZeroModel(), features, targets, ((0, 1), (2,)), seed=3
    )

    assert means[0] > 0.2
    assert means[1] == pytest.approx(0.0)
    assert deviations[1] == pytest.approx(0.0)


def test_feature_importance_rows_refuses_impurity_without_a_tree():
    """`mdi` is rejected, not silently swapped, for a model without impurity."""

    class NoImpurityModel:
        estimators_ = [object()]

        def predict(self, features):
            return np.zeros((len(features), BIT_COUNT), dtype=int)

    with pytest.raises(ConfigError, match="tree-based model"):
        feature_importance_rows(
            NoImpurityModel(),
            method="mdi",
            feature_names=("feature_1",),
            features=np.zeros((4, 1)),
            targets=np.zeros((4, BIT_COUNT), dtype=int),
            seed=1,
            correlation_threshold=0.7,
        )
