"""Shared fixtures, plus the rule that sorts every test into a category.

Four explicit markers (``config``, ``output``, ``integration`` and
``parity``) are applied in the test files. Anything carrying none of the first three is
the pipeline's statistics and decisions, and is marked ``core`` here so
``pytest -m core`` selects it without every such test naming itself. A test
may carry several: consolidation tests that check a computed average by reading
the TSV it was written to are both ``core`` and ``output``.
"""

import csv

import pytest

_EXPLICIT_CATEGORIES = {"config", "output", "integration"}


def pytest_collection_modifyitems(items):
    """Add ``core`` to every test that carries none of the explicit categories."""
    for item in items:
        if not any(marker.name in _EXPLICIT_CATEGORIES for marker in item.iter_markers()):
            item.add_marker(pytest.mark.core)


@pytest.fixture()
def summary_statistics_tsv(tmp_path):
    """Write a small ``summary_statistics.tsv`` the ML trainers and tuner can fit.

    Four label classes, four rows each, with a string column and four numeric
    features that separate the classes, so a 25% hold-out keeps one row per
    class and three cross-validation folds remain feasible.

    Args:
        tmp_path: pytest temporary directory.

    Returns:
        The path to the written TSV.
    """
    classes = [
        ("100001", "BC", 1.0, 0.6, 0.4, 0.7),
        ("010010", "AC", 1.2, 0.8, 0.6, 0.8),
        ("001100", "BC", 1.4, 0.9, 0.7, 0.85),
        ("110000", "AC", 1.6, 1.0, 0.8, 0.92),
    ]
    rows = [
        {
            "class": label,
            "dis1_topology": topology,
            **{
                f"feature_{index}": f"{value + 0.01 * repeat:.2f}"
                for index, value in enumerate(features, start=1)
            },
        }
        for repeat in range(4)
        for label, topology, *features in classes
    ]
    path = tmp_path / "summary_statistics.tsv"
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path
