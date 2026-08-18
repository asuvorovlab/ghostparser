"""Cover the ML label contract, evaluation metrics, distribution summaries, and CV-fold policy.

Expected values are computed by hand from small fixed arrays (see
``tests/TEST_IO.md``) rather than by running a model, so each assertion pins a
definition rather than a training outcome.
"""

import numpy as np
import pytest

from ghostparser.ml.ml_utils import (
    BIT_COUNT,
    BIT_LABELS,
    auto_cv_folds,
    bit_distribution,
    build_feature_importance_rows,
    build_prediction_rows,
    evaluate_predictions,
    is_valid_bitstring,
    parse_classes,
    read_tsv_rows,
    select_feature_names,
    summarize_distribution,
)


def test_bit_labels_stay_in_step_with_the_bit_count():
    """Every bit has exactly one name, so a label can be read back positionally.

    The width itself is pinned behaviorally by ``test_is_valid_bitstring``;
    what this guards is the two constants drifting apart.
    """
    assert len(BIT_LABELS) == BIT_COUNT
    assert len(set(BIT_LABELS)) == BIT_COUNT


@pytest.mark.parametrize(
    "value,expected",
    [
        ("000000", True),
        ("111111", True),
        ("100001", True),
        ("10000", False),  # too short
        ("1000010", False),  # too long
        ("100002", False),  # non-binary digit
        ("10000a", False),  # non-binary character
        ("", False),
    ],
)
def test_is_valid_bitstring(value, expected):
    """Only six-character strings of 0/1 are valid class labels."""
    assert is_valid_bitstring(value) is expected


def test_parse_classes_round_trips_bitstrings():
    """Labels parse to a binary matrix whose rows re-join to the input strings."""
    raw = ["100001", "000000", "111111", " 010010 "]
    matrix, labels = parse_classes(raw)

    assert matrix.shape == (4, BIT_COUNT)
    # Whitespace is stripped, so the stored labels are the trimmed originals.
    assert labels == ["100001", "000000", "111111", "010010"]
    # Every row re-joins to its label: the matrix is a faithful bit expansion.
    assert ["".join(str(bit) for bit in row) for row in matrix] == labels
    assert matrix[0].tolist() == [1, 0, 0, 0, 0, 1]
    # Set bits per label: 100001 -> 2, 000000 -> 0, 111111 -> 6, 010010 -> 2.
    assert matrix.sum() == 2 + 0 + 6 + 2


@pytest.mark.parametrize("bad_label", ["10000", "1000010", "10000x", ""])
def test_parse_classes_rejects_malformed_labels(bad_label):
    """A non-bitstring label raises ValueError naming the expected format."""
    with pytest.raises(ValueError, match="Invalid classes label"):
        parse_classes(["100001", bad_label])


def test_select_feature_names_excludes_the_target_column():
    """Feature selection keeps header order and drops the target column."""
    fieldnames = ["class", "feature_1", "feature_2", "dis1_topology"]
    assert select_feature_names(fieldnames, "class") == (
        "feature_1",
        "feature_2",
        "dis1_topology",
    )


def test_evaluate_predictions_matches_hand_computed_metrics():
    """Metrics on a 2x6 prediction pair match their definitions."""
    y_true = np.array([[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0]])
    # Row 0 is predicted exactly; row 1 flips the second bit (1 -> 0).
    y_pred = np.array([[1, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0]])

    metrics = evaluate_predictions(y_true, y_pred)

    # 1 of 2 rows matches exactly.
    assert metrics["exact_match_accuracy"] == pytest.approx(0.5)
    # 11 of 12 bit positions agree.
    assert metrics["bitwise_accuracy"] == pytest.approx(11 / 12)
    # Hamming loss is the complement of bitwise accuracy.
    assert metrics["hamming_loss"] == pytest.approx(1 / 12)
    assert set(metrics["per_bit"]) == set(BIT_LABELS)
    # ghost_into_A: predicted correctly in both rows.
    assert metrics["per_bit"]["ghost_into_A"]["accuracy"] == pytest.approx(1.0)
    # ghost_into_B: the one positive was missed -> recall 0, support 1.
    assert metrics["per_bit"]["ghost_into_B"]["recall"] == pytest.approx(0.0)
    assert metrics["per_bit"]["ghost_into_B"]["support"] == 1


def test_build_prediction_rows_reports_matched_label_count():
    """Each prediction row carries its per-bit agreement count and exact-match flag."""
    y_true = np.array([[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0]])
    y_pred = np.array([[1, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0]])

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


def test_summarize_distribution_counts_and_fractions():
    """Label counts and fractions are reported per distinct label, sorted."""
    summary = summarize_distribution(["100001", "000000", "100001", "111111"])

    assert list(summary) == ["000000", "100001", "111111"]
    assert summary["100001"] == {"count": 2, "fraction": pytest.approx(0.5)}
    assert summary["000000"] == {"count": 1, "fraction": pytest.approx(0.25)}


def test_bit_distribution_counts_positives_per_bit():
    """Per-bit positive counts and fractions come from the target column sums."""
    targets = np.array([[1, 0, 0, 0, 0, 1], [1, 1, 0, 0, 0, 0]])

    distribution = bit_distribution(targets)

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


class TestAutoCvFolds:
    """Fold selection under each rare-class policy."""

    def test_keeps_requested_folds_when_every_class_is_large_enough(self):
        """No reduction happens when the smallest class supports the request."""
        labels = np.array(["a"] * 5 + ["b"] * 5)
        folds, warnings = auto_cv_folds(labels, requested_folds=5, policy="warn_reduce_cv")
        assert folds == 5
        assert warnings == []

    def test_reduces_folds_to_the_smallest_class_count(self):
        """warn_reduce_cv caps folds at the smallest class size and warns."""
        labels = np.array(["a"] * 10 + ["b"] * 3)
        folds, warnings = auto_cv_folds(labels, requested_folds=5, policy="warn_reduce_cv")
        assert folds == 3
        assert any("Reduced CV folds from 5 to 3" in warning for warning in warnings)

    def test_warn_skip_cv_skips_when_reduction_would_be_needed(self):
        """warn_skip_cv returns no folds rather than silently reducing them."""
        labels = np.array(["a"] * 10 + ["b"] * 3)
        folds, warnings = auto_cv_folds(labels, requested_folds=5, policy="warn_skip_cv")
        assert folds is None
        assert any("Skipped cross-validation" in warning for warning in warnings)

    def test_singleton_class_skips_cross_validation(self):
        """A class with a single sample makes stratified CV impossible."""
        labels = np.array(["a"] * 10 + ["b"])
        folds, warnings = auto_cv_folds(labels, requested_folds=5, policy="warn_reduce_cv")
        assert folds is None
        assert any("fewer than 2 samples" in warning for warning in warnings)

    def test_error_policy_raises_on_a_singleton_class(self):
        """The error policy turns an infeasible split into an exception."""
        labels = np.array(["a"] * 10 + ["b"])
        with pytest.raises(ValueError, match="fewer than 2 samples"):
            auto_cv_folds(labels, requested_folds=5, policy="error")

    def test_empty_labels_return_no_folds(self):
        """An empty label array yields no folds and an explanatory warning."""
        folds, warnings = auto_cv_folds(np.array([]), requested_folds=5, policy="error")
        assert folds is None
        assert warnings == ["No labels available for cross-validation"]


def test_read_tsv_rows_rejects_a_header_only_file(tmp_path):
    """A TSV with a header but no data rows is an error."""
    path = tmp_path / "empty.tsv"
    path.write_text("class\tfeature_1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no data rows"):
        read_tsv_rows(str(path))


def test_read_tsv_rows_reads_records(tmp_path):
    """A well-formed TSV reads into one dict per data row."""
    path = tmp_path / "data.tsv"
    path.write_text("class\tfeature_1\n100001\t1.5\n000000\t2.5\n", encoding="utf-8")
    rows = read_tsv_rows(str(path))
    assert rows == [
        {"class": "100001", "feature_1": "1.5"},
        {"class": "000000", "feature_1": "2.5"},
    ]
