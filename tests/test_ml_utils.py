"""ML utilities: the label contract, feature encoding, metrics, distributions,
CV-fold policy and feature importance.

Expected values are computed by hand from small fixed arrays (see
``tests/TEST_IO.md``) rather than by running a model, so each assertion pins a
definition rather than a training outcome.
"""

import numpy as np
import pytest

from ghostparser.config import ConfigError
from ghostparser.ml.ml_utils import (
    BIT_COUNT,
    BIT_LABELS,
    auto_cv_folds,
    bit_distribution,
    build_64_class_confusion_matrix,
    build_feature_importance_rows,
    build_prediction_rows,
    classes_are_balanced,
    correlation_feature_groups,
    evaluate_predictions,
    feature_importance_rows,
    grouped_permutation_importance,
    is_valid_bitstring,
    parse_classes,
    per_class_recall,
    read_tsv_rows,
    row_normalize_confusion_matrix,
    rows_to_matrix,
    save_evaluation_figures,
    select_feature_names,
    summarize_distribution,
)


def test_rows_to_matrix_one_hot_encodes_strings_and_drops_the_target():
    """Numeric columns pass through, the target is dropped, each string column is one-hot.

    Feature order is asserted because the model pickle is reused across runs:
    one-hot blocks come first, one per string column in sorted-category order
    and never sharing a category space, then the numerics in input order, and
    a reader scoring new data must reproduce that layout.
    """
    rows = [
        {"class": "001010", "topology_a": "BC", "topology_b": "XY", "feature_a": "10.0", "feature_x": "1.5"},
        {"class": "001011", "topology_a": "AC", "topology_b": "XZ", "feature_a": "20.0", "feature_x": "2.5"},
    ]

    matrix = rows_to_matrix(rows, target_column="class")

    assert matrix.feature_names == (
        "topology_a=AC",
        "topology_a=BC",
        "topology_b=XY",
        "topology_b=XZ",
        "feature_a",
        "feature_x",
    )
    assert matrix.train_features.tolist() == [
        [0.0, 1.0, 1.0, 0.0, 10.0, 1.5],
        [1.0, 0.0, 0.0, 1.0, 20.0, 2.5],
    ]
    assert matrix.train_labels == ["001010", "001011"]


def test_rows_to_matrix_rejects_string_features():
    """A high-cardinality string column is refused rather than one-hot exploded.

    Eight distinct values exceed the 7-category cap, which exists because a
    free-text column would otherwise silently add a column per value.
    """
    rows = [
        {"class": "001010", "dis1_topology": value, "feature_x": "1.0"}
        for value in "ABCDEFGH"
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


def test_per_class_recall_skips_classes_without_hold_out_rows():
    """Only classes with hold-out rows get a recall; absent ones are not scored 0."""
    y_true = np.asarray([[0] * BIT_COUNT, [0] * BIT_COUNT, [1] * BIT_COUNT])
    y_pred = np.asarray([[0] * BIT_COUNT, [0] * (BIT_COUNT - 1) + [1], [1] * BIT_COUNT])

    set_bit_counts, recalls = per_class_recall(
        build_64_class_confusion_matrix(y_true, y_pred)
    )

    np.testing.assert_array_equal(set_bit_counts, [0, BIT_COUNT])
    np.testing.assert_allclose(recalls, [0.5, 1.0])


_ALL_CLASSES = [format(index, f"0{BIT_COUNT}b") for index in range(2**BIT_COUNT)]


@pytest.mark.parametrize(
    ("train_labels", "test_labels", "expected"),
    [
        (_ALL_CLASSES * 4, _ALL_CLASSES, True),
        (_ALL_CLASSES * 4 + _ALL_CLASSES[:10] * 2, _ALL_CLASSES, True),
        (_ALL_CLASSES * 3 + _ALL_CLASSES[:10] * 2, _ALL_CLASSES, False),
        (_ALL_CLASSES * 4, _ALL_CLASSES + _ALL_CLASSES[:10], False),
        (_ALL_CLASSES * 4, _ALL_CLASSES[1:], False),
        (_ALL_CLASSES[1:] * 4, _ALL_CLASSES, False),
    ],
    ids=[
        "equal",
        "ratio_at_limit",
        "train_ratio_over_limit",
        "test_ratio_over_limit",
        "test_missing_a_class",
        "train_missing_a_class",
    ],
)
def test_classes_are_balanced_bounds_the_count_ratio(
    train_labels, test_labels, expected
):
    """Balanced means all 64 classes in every partition, no count over 1.5x another."""
    assert classes_are_balanced(train_labels, test_labels) is expected


@pytest.mark.output
def test_save_evaluation_figures_creates_no_folder_without_figures(tmp_path):
    """With neither matrix built, nothing is written and no ``figures/`` appears."""
    figure_paths = save_evaluation_figures(tmp_path, "random_forest", None, None, True)

    assert figure_paths == {}
    assert not (tmp_path / "figures").exists()


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


def test_64_class_matrix_orders_classes_by_set_bits():
    """Both axes run from 000000 up to 111111 by the number of set bits.

    The label list must be non-decreasing in set-bit count, lexical within a
    count, and the count matrix must be permuted along with it: a sort that
    reordered the labels but left the counts in binary order would mislabel
    every off-diagonal cell. Three rows with known true/predicted classes pin
    the permutation: each lands where its labels say, and nowhere else.
    """
    y_true = np.array([[1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [1, 1, 1, 1, 1, 1]])
    y_pred = np.array([[0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1], [1, 1, 1, 1, 1, 1]])

    result = build_64_class_confusion_matrix(y_true, y_pred)
    labels = result["class_labels"]
    matrix = np.asarray(result["matrix"])

    assert len(labels) == 64 and len(set(labels)) == 64
    set_bits = [label.count("1") for label in labels]
    assert set_bits == sorted(set_bits)
    for count in range(BIT_COUNT + 1):
        group = [label for label in labels if label.count("1") == count]
        assert group == sorted(group)
    assert labels[0] == "000000" and labels[-1] == "111111"

    assert matrix.sum() == 3
    assert matrix[labels.index("110000"), labels.index("000011")] == 1
    assert matrix[labels.index("000001"), labels.index("000001")] == 1
    assert matrix[labels.index("111111"), labels.index("111111")] == 1


@pytest.mark.parametrize(
    "value, valid",
    [
        ("000000", True),
        ("111111", True),
        ("100001", True),
        (" 010010 ", True),  # surrounding whitespace is stripped
        ("10000", False),  # too short
        ("1000010", False),  # too long
        ("100002", False),  # non-binary digit
        ("10000a", False),  # non-binary character
        ("", False),
    ],
)
def test_bitstring_labels_are_validated_and_parsed(value, valid):
    """Only six-character strings of 0/1 are class labels, and they parse to their bits.

    A valid label parses to one binary row that re-joins to the trimmed
    label, so the matrix is a faithful bit expansion; a malformed one raises
    naming the expected format, even beside a valid neighbour.
    """
    assert is_valid_bitstring(value.strip()) is valid
    if not valid:
        assert is_valid_bitstring(value) is False
        with pytest.raises(ValueError, match="Invalid classes label"):
            parse_classes(["100001", value])
        return

    matrix, labels = parse_classes([value])
    assert matrix.shape == (1, BIT_COUNT)
    assert labels == [value.strip()]
    assert "".join(str(bit) for bit in matrix[0]) == value.strip()
    assert matrix.sum() == value.count("1")


def test_prediction_metrics_and_rows_match_hand_computation():
    """Metrics and per-row records on a 2x6 prediction pair match their definitions."""
    y_true = np.array([[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0]])
    # Row 0 is predicted exactly; row 1 flips the second bit (1 -> 0).
    y_pred = np.array([[1, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0]])

    metrics = evaluate_predictions(y_true, y_pred)

    # 1 of 2 rows matches exactly; 11 of 12 bit positions agree; Hamming loss
    # is the complement of bitwise accuracy.
    assert metrics["exact_match_accuracy"] == pytest.approx(0.5)
    assert metrics["bitwise_accuracy"] == pytest.approx(11 / 12)
    assert metrics["hamming_loss"] == pytest.approx(1 / 12)
    assert set(metrics["per_bit"]) == set(BIT_LABELS)
    # ghost_into_A: predicted correctly in both rows.
    assert metrics["per_bit"]["ghost_into_A"]["accuracy"] == pytest.approx(1.0)
    # ghost_into_B: the one positive was missed -> recall 0, support 1.
    assert metrics["per_bit"]["ghost_into_B"]["recall"] == pytest.approx(0.0)
    assert metrics["per_bit"]["ghost_into_B"]["support"] == 1

    rows = build_prediction_rows(y_true, y_pred)

    assert len(rows) == 2
    # Row 0 agrees on all 6 bits.
    assert rows[0]["matched_label_count"] == 6
    assert rows[0]["exact_match"] == 1
    assert rows[0]["true_label"] == "100001"
    assert rows[0]["pred_label"] == "100001"
    # Row 1 differs on exactly one bit -> 5 of 6 matched, not an exact match.
    assert rows[1]["matched_label_count"] == 5
    assert rows[1]["exact_match"] == 0
    assert rows[1]["true_label"] == "010000"
    assert rows[1]["pred_label"] == "000000"
    # Per-bit columns are emitted for both the true and predicted vectors.
    assert rows[1]["true_ghost_into_B"] == 1
    assert rows[1]["pred_ghost_into_B"] == 0


def test_distributions_count_labels_and_positives_per_bit():
    """Label counts and fractions are reported per distinct label, sorted; positives per bit from the column sums."""
    summary = summarize_distribution(["100001", "000000", "100001", "111111"])

    assert list(summary) == ["000000", "100001", "111111"]
    assert summary["100001"] == {"count": 2, "fraction": pytest.approx(0.5)}
    assert summary["000000"] == {"count": 1, "fraction": pytest.approx(0.25)}

    distribution = bit_distribution(np.array([[1, 0, 0, 0, 0, 1], [1, 1, 0, 0, 0, 0]]))

    # Column 0 is positive in both rows; column 1 in one; columns 2-4 in none.
    assert distribution["ghost_into_A"] == {
        "positive_count": 2,
        "fraction": pytest.approx(1.0),
    }
    assert distribution["ghost_into_B"] == {
        "positive_count": 1,
        "fraction": pytest.approx(0.5),
    }
    assert distribution["inflow_into_A_from_C"]["positive_count"] == 0


def test_build_feature_importance_rows_sorts_descending():
    """Feature importances are paired with names and ranked most-important first."""
    rows = build_feature_importance_rows(
        ("feature_1", "feature_2", "feature_3"), np.array([0.2, 0.5, 0.3])
    )
    assert [row["feature"] for row in rows] == [
        "feature_2",
        "feature_3",
        "feature_1",
    ]
    assert rows[0]["importance"] == pytest.approx(0.5)


@pytest.mark.parametrize(
    "labels, policy, expected_folds, warning",
    [
        # No reduction happens when the smallest class supports the request.
        (["a"] * 5 + ["b"] * 5, "warn_reduce_cv", 5, None),
        # warn_reduce_cv caps folds at the smallest class size and warns.
        (["a"] * 10 + ["b"] * 3, "warn_reduce_cv", 3, "Reduced CV folds from 5 to 3"),
        # warn_skip_cv returns no folds rather than silently reducing them.
        (["a"] * 10 + ["b"] * 3, "warn_skip_cv", None, "Skipped cross-validation"),
        # A class with a single sample makes stratified CV impossible.
        (["a"] * 10 + ["b"], "warn_reduce_cv", None, "fewer than 2 samples"),
        # The error policy turns an infeasible split into an exception.
        (["a"] * 10 + ["b"], "error", ValueError, "fewer than 2 samples"),
        # An empty label array yields no folds and an explanatory warning.
        ([], "error", None, "No labels available for cross-validation"),
    ],
    ids=["enough", "reduce", "skip", "singleton", "error", "empty"],
)
def test_auto_cv_folds_follows_the_rare_class_policy(labels, policy, expected_folds, warning):
    """Fold selection under each rare-class policy."""
    labels = np.array(labels)
    if expected_folds is ValueError:
        with pytest.raises(ValueError, match=warning):
            auto_cv_folds(labels, requested_folds=5, policy=policy)
        return

    folds, warnings = auto_cv_folds(labels, requested_folds=5, policy=policy)
    assert folds == expected_folds
    if warning is None:
        assert warnings == []
    else:
        assert any(warning in text for text in warnings)


def test_read_tsv_rows_reads_records_and_rejects_a_header_only_file(tmp_path):
    """A well-formed TSV reads into one dict per data row, and its features exclude the target."""
    path = tmp_path / "data.tsv"
    path.write_text("class\tfeature_1\tdis1_topology\n100001\t1.5\tBC\n000000\t2.5\tAC\n", encoding="utf-8")
    rows = read_tsv_rows(str(path))
    assert rows == [
        {"class": "100001", "feature_1": "1.5", "dis1_topology": "BC"},
        {"class": "000000", "feature_1": "2.5", "dis1_topology": "AC"},
    ]
    # Feature selection keeps header order and drops the target column.
    assert select_feature_names(list(rows[0]), "class") == ("feature_1", "dis1_topology")

    empty = tmp_path / "empty.tsv"
    empty.write_text("class\tfeature_1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no data rows"):
        read_tsv_rows(str(empty))
