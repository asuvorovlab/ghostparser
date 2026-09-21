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
    build_64_class_confusion_matrix,
    build_feature_importance_rows,
    build_prediction_rows,
    evaluate_predictions,
    format_bit_label_title,
    is_valid_bitstring,
    parse_classes,
    read_tsv_rows,
    select_feature_names,
    summarize_distribution,
)


def test_bit_labels_and_their_titles_cover_every_bit():
    """Every bit has exactly one name and one prose title with the taxon letters intact.

    A label is read back positionally, so the name list and the bit count
    must not drift apart. The title trap is `str.capitalize`, which lower-cases
    everything after the first character -- turning `ghost_into_A` into
    `Ghost into a` and renaming the taxon -- so every title is pinned.
    """
    assert len(BIT_LABELS) == BIT_COUNT
    assert len(set(BIT_LABELS)) == BIT_COUNT
    assert {label: format_bit_label_title(label) for label in BIT_LABELS} == {
        "ghost_into_A": "Ghost into A",
        "ghost_into_B": "Ghost into B",
        "inflow_into_A_from_C": "Inflow into A from C",
        "inflow_into_B_from_C": "Inflow into B from C",
        "outflow_from_A_to_C": "Outflow from A to C",
        "outflow_from_B_to_C": "Outflow from B to C",
    }


def test_64_class_matrix_orders_classes_by_set_bits():
    """Both axes run from 000000 up to 111111 by the number of set bits.

    The label list must be non-decreasing in set-bit count, lexical within a
    count, and the count matrix must be permuted along with it -- a sort that
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
