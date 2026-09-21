"""Tests for the orchestrator's consolidation stage.

One hand-built result set drives the whole stage: every edge and ghost
target it produces is derived in ``tests/TEST_IO.md``, and the TSVs the stage
writes are read back and checked against those derivations.
"""

from types import SimpleNamespace

import pytest

from ghostparser.orchestrator.consolidation import (
    _collect_counts,
    _collect_non_sister_counts,
    _sampled_introgression_presence,
    generate_introgression_maps,
)

_SPECIES_TREE = "(((A:1,B:1):1,C:1):1,D:1);\n"
_BALANCED_TREE = "((A:1,B:1):1,(C:1,D:1):1);\n"
_OUTGROUP_TREE = "(((A:1,B:1):1,C:1):1,OG:1);\n"

_CONSOLIDATION_FILES = (
    "introgression_matrix_inflow_outflow.tsv",
    "introgression_matrix_inflow_outflow_raw_sum.tsv",
    "introgression_matrix_inflow_outflow_supporting_count.tsv",
    "introgression_ghost_target_strength.tsv",
    "introgression_ghost_target_strength_raw_sum.tsv",
    "introgression_ghost_target_strength_supporting_count.tsv",
    "introgression_taxa_order.tsv",
    "introgression_matrix_sampled_non_sister.tsv",
)


def _row(triplet, classification, dis1_topology, bootstrap_value):
    """Build one result row with the four fields consolidation reads."""
    return SimpleNamespace(
        triplet=triplet,
        classification=classification,
        dis1_topology=dis1_topology,
        bootstrap_value=bootstrap_value,
    )


# Every kind of event the mapping produces, with repeats so the averages
# divide by the supporting count and a no_introgression row that must count
# for nothing. Edges are (source -> target); see TEST_IO.md for the mapping.
_RESULTS = [
    _row(("A", "B", "C"), "inflow_introgression", "BC", 0.6),   # C -> B
    _row(("A", "B", "C"), "inflow_introgression", "BC", 0.2),   # C -> B
    _row(("A", "B", "C"), "no_introgression", "BC", 0.0),       # no edge
    _row(("A", "B", "C"), "inflow_introgression", "AC", 0.5),   # C -> A
    _row(("A", "B", "C"), "outflow_introgression", "BC", 0.3),  # B -> C
    _row(("A", "B", "D"), "ghost_introgression", "BC", 0.8),    # ghost target A
    _row(("A", "B", "D"), "ghost_introgression", "BC", 0.4),    # ghost target A
    _row(("D", "B", "C"), "ghost_introgression", "BC", 0.4),    # ghost target D
    _row(("A", "C", "D"), "ghost_introgression", "AC", 0.5),    # ghost target C
]


def _read_matrix(path):
    """Read a target x source matrix TSV into ``{target: {source: value}}``."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")[1:]
    matrix = {}
    for line in lines[1:]:
        target, *values = line.split("\t")
        matrix[target] = dict(zip(header, (float(value) for value in values)))
    return matrix


def _read_column(path, column):
    """Read one named column of a per-taxon TSV into ``{taxon: value}``."""
    lines = path.read_text().splitlines()
    index = lines[0].split("\t").index(column)
    return {
        line.split("\t")[0]: float(line.split("\t")[index]) for line in lines[1:]
    }


@pytest.mark.core
@pytest.mark.output
def test_generate_introgression_maps_writes_every_sheet_from_the_results(tmp_path):
    """The stage writes the figure and eight TSVs, with the hand-derived values.

    The averages divide each edge's and ghost target's bootstrap weight by the
    number of triplets that produced it -- never by the triplets that merely
    contain the taxa -- so the raw-sum and supporting-count sheets must
    reproduce them. The ghost sheet flags a taxon that is also the target of a
    non-zero sampled edge, the non-sister sheet counts the pairs the normalized
    triplets leave apart, and a second call with overwriting disabled leaves
    an existing folder alone and writes beside it.
    """
    species_tree = tmp_path / "species.tree"
    species_tree.write_text(_SPECIES_TREE)
    output_dir = tmp_path / "out"

    artifacts = generate_introgression_maps(
        _RESULTS, species_tree_path=str(species_tree), output_dir=str(output_dir)
    )

    assert artifacts.taxa_count == 4
    assert artifacts.non_ghost_edge_count == 3
    assert artifacts.ghost_target_count == 3
    assert artifacts.plot_path == str(output_dir / "introgression_combined.png")
    assert (output_dir / "introgression_combined.png").exists()
    data = output_dir / "consolidation_data"
    for name in _CONSOLIDATION_FILES:
        assert (data / name).exists(), name

    # Sampled edges, target x source: C->B carries 0.6 + 0.2 over two rows.
    avg = _read_matrix(data / "introgression_matrix_inflow_outflow.tsv")
    raw = _read_matrix(data / "introgression_matrix_inflow_outflow_raw_sum.tsv")
    count = _read_matrix(
        data / "introgression_matrix_inflow_outflow_supporting_count.tsv"
    )
    assert list(avg) == ["A", "B", "C", "D"]
    assert (avg["B"]["C"], raw["B"]["C"], count["B"]["C"]) == (
        pytest.approx(0.4), pytest.approx(0.8), 2,
    )
    assert (avg["A"]["C"], raw["A"]["C"], count["A"]["C"]) == (0.5, 0.5, 1)
    assert (avg["C"]["B"], raw["C"]["B"], count["C"]["B"]) == (0.3, 0.3, 1)
    assert sum(sum(row.values()) for row in raw.values()) == pytest.approx(1.6)
    assert max(value for row in avg.values() for value in row.values()) == 0.5

    # Ghost targets: A carries 0.8 + 0.4 over two rows, C and D one row each.
    strength = _read_column(data / "introgression_ghost_target_strength.tsv", "raw_strength")
    ghost_raw = _read_column(data / "introgression_ghost_target_strength_raw_sum.tsv", "raw_sum")
    ghost_count = _read_column(
        data / "introgression_ghost_target_strength_supporting_count.tsv", "supporting_count"
    )
    assert strength == pytest.approx({"A": 0.6, "B": 0.0, "C": 0.5, "D": 0.4})
    assert ghost_raw == pytest.approx({"A": 1.2, "B": 0.0, "C": 0.5, "D": 0.4})
    assert ghost_count == {"A": 2, "B": 0, "C": 1, "D": 1}
    # A, B and C are targets of a sampled edge; D has ghost signal only.
    flags = _read_column(
        data / "introgression_ghost_target_strength.tsv", "has_sampled_introgression"
    )
    assert flags == {"A": 1, "B": 1, "C": 1, "D": 0}

    # Non-sister pairs: five A,B,C rows leave {A,C} and {B,C} apart, and so on.
    non_sister = _read_matrix(data / "introgression_matrix_sampled_non_sister.tsv")
    expected_pairs = {("A", "C"): 5, ("B", "C"): 6, ("A", "D"): 3, ("B", "D"): 2, ("C", "D"): 2}
    for first in "ABCD":
        for second in "ABCD":
            expected = expected_pairs.get(tuple(sorted((first, second))), 0)
            assert non_sister[first][second] == expected, (first, second)

    # With overwriting disabled an existing folder is kept and a suffixed one
    # is written beside it.
    existing = tmp_path / "introgression"
    existing.mkdir()
    (existing / "stale.txt").write_text("stale")
    again = generate_introgression_maps(
        _RESULTS,
        species_tree_path=str(species_tree),
        output_dir=str(existing),
        overwrite=False,
    )
    assert again.plot_path == str(tmp_path / "introgression_1" / "introgression_combined.png")
    assert (tmp_path / "introgression_1" / "introgression_combined.png").exists()
    assert (existing / "stale.txt").exists()


@pytest.mark.core
@pytest.mark.output
@pytest.mark.parametrize(
    "tree, results, kwargs, absent, expected_count",
    [
        (_BALANCED_TREE, [_row(("A", "B", "C"), "ghost_introgression", "AC", 1.0)], {}, None, 4),
        (
            _BALANCED_TREE,
            [_row(("A", "B", "C"), "ghost_introgression", "AC", 1.0)],
            {"plot_taxa": ["A", "B", "C"]},
            "D",
            3,
        ),
        (
            _OUTGROUP_TREE,
            [
                _row(("A", "B", "C"), "ghost_introgression", "AC", 1.0),
                _row(("A", "B", "OG"), "ghost_introgression", "BC", 0.4),
            ],
            {"outgroups": ["OG"]},
            "OG",
            3,
        ),
    ],
    ids=["full_tree", "plot_taxa", "outgroups"],
)
def test_generate_introgression_maps_selects_the_requested_taxa(
    tmp_path, tree, results, kwargs, absent, expected_count
):
    """`plot_taxa` and `outgroups` decide which taxa reach every output.

    The matrix header and the ghost rows must agree with each other and with
    ``taxa_count``, so a taxon dropped from one output is dropped from all.
    """
    species_tree = tmp_path / "species.tree"
    species_tree.write_text(tree)
    output_dir = tmp_path / "out"

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
        **kwargs,
    )

    data = output_dir / "consolidation_data"
    taxa_from_matrix = list(_read_matrix(data / "introgression_matrix_inflow_outflow.tsv"))
    taxa_from_ghost = list(
        _read_column(data / "introgression_ghost_target_strength.tsv", "raw_strength")
    )

    assert artifacts.taxa_count == expected_count
    assert taxa_from_matrix == taxa_from_ghost
    assert len(taxa_from_matrix) == expected_count
    if absent is not None:
        assert absent not in taxa_from_matrix


def test_count_helpers_follow_the_classification():
    """The counting helpers count events, not co-occurrence.

    A ``no_introgression`` row contains the same taxa as its neighbours but
    produces no edge and no ghost target, so it raises neither count; an
    unrelated triplet maps to a different edge. Non-sister counts come from
    every row -- for a normalized triplet ``(A, B, C)`` the pairs ``{A, C}``
    and ``{B, C}`` -- so the same ``no_introgression`` row does count there.
    A taxon is flagged as having sampled introgression only when it is the
    target of a non-zero sampled edge.
    """
    results = [
        _row(("A", "B", "C"), "inflow_introgression", "BC", 0.8),   # C -> B
        _row(("A", "B", "C"), "no_introgression", "BC", 0.0),       # nothing
        _row(("A", "B", "D"), "inflow_introgression", "BC", 0.5),   # D -> B
        _row(("A", "B", "C"), "ghost_introgression", "BC", 0.9),    # target A
        _row(("A", "C", "D"), "ghost_introgression", "AC", 0.5),    # target C
    ]
    non_ghost_counts, ghost_counts = _collect_counts(results)
    assert non_ghost_counts == {("C", "B"): 1, ("D", "B"): 1}
    assert ghost_counts == {"A": 1, "C": 1}

    assert _collect_non_sister_counts(results) == {
        ("A", "C"): 3,
        ("B", "C"): 3,
        ("A", "D"): 2,
        ("B", "D"): 1,
        ("C", "D"): 1,
    }

    presence = _sampled_introgression_presence(
        ["A", "B", "C", "D"], {("C", "A"): 0.6, ("A", "B"): 0.0, ("D", "C"): 0.3}
    )
    assert presence == {"A": 1, "B": 0, "C": 1, "D": 0}
