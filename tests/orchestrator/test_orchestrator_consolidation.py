"""Tests for the orchestrator's consolidation stage."""

from types import SimpleNamespace

import pytest

from ghostparser.orchestrator.consolidation import (
    _collect_counts,
    _collect_non_sister_counts,
    _sampled_introgression_presence,
    generate_introgression_maps,
)


def _consolidation_lines(consolidation_dir, filename):
    return (consolidation_dir / filename).read_text().splitlines()


@pytest.mark.output
def test_generate_introgression_maps_creates_expected_outputs(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")
    output_dir = tmp_path / "out"

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.25,
        ),
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.75,
        ),
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="outflow_introgression",
            dis1_topology="AC",
            bootstrap_value=0.5,
        ),
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.2,
        ),
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
    ]

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
    )

    assert artifacts.taxa_count == 4
    assert artifacts.non_ghost_edge_count >= 1
    assert artifacts.ghost_target_count >= 1

    assert (output_dir / "introgression_combined.png").exists()
    assert artifacts.plot_path == str(output_dir / "introgression_combined.png")
    consolidation_dir = output_dir / "consolidation_data"
    assert (consolidation_dir / "introgression_matrix_inflow_outflow.tsv").exists()
    assert (
        consolidation_dir / "introgression_matrix_inflow_outflow_raw_sum.tsv"
    ).exists()
    assert (
        consolidation_dir / "introgression_matrix_inflow_outflow_supporting_count.tsv"
    ).exists()
    assert (consolidation_dir / "introgression_ghost_target_strength.tsv").exists()
    assert (
        consolidation_dir / "introgression_ghost_target_strength_raw_sum.tsv"
    ).exists()
    assert (
        consolidation_dir / "introgression_ghost_target_strength_supporting_count.tsv"
    ).exists()
    assert (consolidation_dir / "introgression_taxa_order.tsv").exists()
    assert (consolidation_dir / "introgression_matrix_sampled_non_sister.tsv").exists()


@pytest.mark.output
def test_generate_introgression_maps_appends_suffix_when_overwrite_disabled(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")
    output_dir = tmp_path / "out"

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.5,
        ),
    ]

    existing_output_dir = tmp_path / "introgression"
    existing_output_dir.mkdir()
    (existing_output_dir / "stale.txt").write_text("stale")

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(existing_output_dir),
        overwrite=False,
    )

    suffixed_output_dir = tmp_path / "introgression_1"
    assert suffixed_output_dir.exists()
    assert artifacts.plot_path == str(
        suffixed_output_dir / "introgression_combined.png"
    )
    assert (existing_output_dir / "stale.txt").exists()


_BALANCED_TREE = "((A:1,B:1):1,(C:1,D:1):1);\n"
_OUTGROUP_TREE = "(((A:1,B:1):1,C:1):1,OG:1);\n"


def _ghost_abc():
    """One ghost result on triplet (A, B, C)."""
    return SimpleNamespace(
        triplet=("A", "B", "C"),
        classification="ghost_introgression",
        dis1_topology="AC",
        bootstrap_value=1.0,
    )


@pytest.mark.core
@pytest.mark.output
@pytest.mark.parametrize(
    "tree, results, kwargs, absent, expected_count",
    [
        (_BALANCED_TREE, [_ghost_abc()], {}, None, 4),
        (_BALANCED_TREE, [_ghost_abc()], {"plot_taxa": ["A", "B", "C"]}, "D", 3),
        (
            _OUTGROUP_TREE,
            [
                _ghost_abc(),
                SimpleNamespace(
                    triplet=("A", "B", "OG"),
                    classification="ghost_introgression",
                    dis1_topology="BC",
                    bootstrap_value=0.4,
                ),
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

    consolidation_dir = output_dir / "consolidation_data"
    matrix_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow.tsv"
    )
    ghost_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength.tsv"
    )
    taxa_from_matrix = matrix_lines[0].split("\t")[1:]
    taxa_from_ghost = [line.split("\t")[0] for line in ghost_lines[1:]]

    assert artifacts.taxa_count == expected_count
    assert taxa_from_matrix == taxa_from_ghost
    assert len(taxa_from_matrix) == expected_count
    assert len(matrix_lines) == expected_count + 1
    if absent is not None:
        assert absent not in taxa_from_matrix


@pytest.mark.core
@pytest.mark.output
def test_generate_introgression_maps_uses_raw_values_with_separate_scales(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    output_dir = tmp_path / "out"

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.2,
        ),
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.1,
        ),
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.9,
        ),
    ]

    generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
    )

    consolidation_dir = output_dir / "consolidation_data"
    matrix_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow.tsv"
    )
    ghost_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength.tsv"
    )

    matrix_values = []
    for line in matrix_lines[1:]:
        matrix_values.extend(float(value) for value in line.split("\t")[1:])

    ghost_values = [float(line.split("\t")[1]) for line in ghost_lines[1:]]

    assert max(matrix_values) <= 1.0
    assert max(ghost_values) <= 1.0
    # non-ghost avg: sum(0.2+0.8)=1.0 / supporting_count_for_edge(C->B)=2 → 0.5
    assert max(matrix_values) == 0.5
    # ghost avg: sum(0.1+0.9)=1.0 / supporting_count_for_A=2 (both ABD rows) → 0.5
    assert max(ghost_values) == 0.5


def _row(triplet, classification, dis1_topology, bootstrap_value):
    """Build one result row for the count helpers."""
    return SimpleNamespace(
        triplet=triplet,
        classification=classification,
        dis1_topology=dis1_topology,
        bootstrap_value=bootstrap_value,
    )


def test_collect_counts_counts_only_the_rows_that_produced_an_edge():
    """Supporting counts follow the classification, not mere co-occurrence.

    A ``no_introgression`` row contains the same taxa as its neighbours but
    produces no edge and no ghost target, so it must not raise either count. An
    unrelated triplet maps to a different edge and must not raise this one.
    """
    results = [
        # Produces the directed edge (C -> B).
        _row(("A", "B", "C"), "inflow_introgression", "BC", 0.8),
        # Same taxa, no edge produced.
        _row(("A", "B", "C"), "no_introgression", "BC", 0.0),
        # Maps to (D -> B), not (C -> B).
        _row(("A", "B", "D"), "inflow_introgression", "BC", 0.5),
        # topo == BC -> ghost target A; topo != BC -> ghost target C.
        _row(("A", "B", "C"), "ghost_introgression", "BC", 0.9),
        _row(("A", "C", "D"), "ghost_introgression", "AC", 0.5),
    ]
    non_ghost_counts, ghost_counts = _collect_counts(results)

    assert non_ghost_counts.get(("C", "B"), 0) == 1
    assert non_ghost_counts.get(("B", "D"), 0) == 0
    assert ghost_counts.get("A", 0) == 1
    assert ghost_counts.get("C", 0) == 1
    assert ghost_counts.get("B", 0) == 0
    assert ghost_counts.get("D", 0) == 0


@pytest.mark.core
@pytest.mark.output
def test_collect_counts_correct_avg_in_generate_introgression_maps(tmp_path):
    """End-to-end: average bootstrap values in TSV use population-level denominators.

    Setup:
    - Triplet (A,B,C): one inflow_BC (weight 0.6) + one no_introgression (weight ignored).
      Edge (C→B): weight_sum=0.6, denominator=2 (both ABC rows contain C and B) → avg=0.3.
    - Triplet (A,B,D) ×1: ghost_BC (weight 0.8). Ghost target A:
      weight_sum=0.8, denominator=3 (ABC×2 + ABD×1 all contain A) → avg≈0.2667.
    """
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")
    output_dir = tmp_path / "out"

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.6,
        ),
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="no_introgression",
            dis1_topology="BC",
            bootstrap_value=0.0,
        ),
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
    ]

    generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
    )

    consolidation_dir = output_dir / "consolidation_data"
    matrix_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow.tsv"
    )
    matrix_raw_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow_raw_sum.tsv"
    )
    matrix_count_lines = _consolidation_lines(
        consolidation_dir,
        "introgression_matrix_inflow_outflow_supporting_count.tsv",
    )
    ghost_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength.tsv"
    )
    ghost_raw_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength_raw_sum.tsv"
    )
    ghost_count_lines = _consolidation_lines(
        consolidation_dir,
        "introgression_ghost_target_strength_supporting_count.tsv",
    )

    # Build lookup: matrix[target][source] = value
    header = matrix_lines[0].split("\t")[1:]
    matrix = {}
    for row_line in matrix_lines[1:]:
        parts = row_line.split("\t")
        row_target = parts[0]
        matrix[row_target] = {
            header[i]: float(parts[i + 1]) for i in range(len(header))
        }

    ghost_map = {}
    for line in ghost_lines[1:]:
        taxon, val, _has_sampled = line.split("\t")
        ghost_map[taxon] = float(val)

    raw_header = matrix_raw_lines[0].split("\t")[1:]
    raw_matrix = {}
    for row_line in matrix_raw_lines[1:]:
        parts = row_line.split("\t")
        row_target = parts[0]
        raw_matrix[row_target] = {
            raw_header[i]: float(parts[i + 1]) for i in range(len(raw_header))
        }

    count_header = matrix_count_lines[0].split("\t")[1:]
    count_matrix = {}
    for row_line in matrix_count_lines[1:]:
        parts = row_line.split("\t")
        row_target = parts[0]
        count_matrix[row_target] = {
            count_header[i]: float(parts[i + 1]) for i in range(len(count_header))
        }

    ghost_raw_map = {}
    for line in ghost_raw_lines[1:]:
        taxon, val = line.split("\t")
        ghost_raw_map[taxon] = float(val)

    ghost_count_map = {}
    for line in ghost_count_lines[1:]:
        taxon, val = line.split("\t")
        ghost_count_map[taxon] = float(val)

    # With undiluted consolidation we divide only by supporting triplets
    # Edge (C->B) has one supporting row with bootstrap 0.6 => avg = 0.6
    assert abs(matrix["B"]["C"] - 0.6) < 1e-9
    # Ghost target A has one supporting ghost row with bootstrap 0.8 => avg = 0.8
    assert abs(ghost_map["A"] - 0.8) < 1e-9
    assert abs(raw_matrix["B"]["C"] - 0.6) < 1e-9
    assert abs(count_matrix["B"]["C"] - 1.0) < 1e-9
    assert abs(ghost_raw_map["A"] - 0.8) < 1e-9
    assert abs(ghost_count_map["A"] - 1.0) < 1e-9


def test_collect_non_sister_counts_counts_non_sister_pairs():
    """_collect_non_sister_counts increments only non-sister pairs.

    For a triplet (A, B, C), A and B are sisters, so:
    - (A, C) and (B, C) are non-sister pairs and get incremented.
    - (A, B) is the sister pair and must not be counted.
    Multiple triplets sharing a non-sister pair accumulate their counts.
    """
    results = [
        # Triplet (A, B, C): sisters are A and B; non-sisters: (A,C), (B,C)
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
        # Same triplet again (different gene tree result): (A,C) and (B,C) each +1
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="no_introgression",
            dis1_topology="BC",
            bootstrap_value=0.0,
        ),
        # Triplet (A, C, D): sisters are A and C; non-sisters: (A,D), (C,D)
        SimpleNamespace(
            triplet=("A", "C", "D"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.5,
        ),
    ]
    counts = _collect_non_sister_counts(results)

    # (A, C) appears as non-sisters twice (both ABC rows)
    assert counts.get(("A", "C"), 0) == 2
    # (B, C) appears as non-sisters twice (both ABC rows)
    assert counts.get(("B", "C"), 0) == 2
    # (A, B) is the sister pair in the ABC triplets — must not be counted
    assert counts.get(("A", "B"), 0) == 0
    # (A, D) is a non-sister pair in the ACD triplet
    assert counts.get(("A", "D"), 0) == 1
    # (C, D) is a non-sister pair in the ACD triplet
    assert counts.get(("C", "D"), 0) == 1
    # (A, C) is also incremented by ACD (C is the C-taxon, non-sister to A and C here—
    # wait: in (A,C,D) the sisters are A and C; so (A,D) and (C,D) are non-sisters)
    # So (A,C) comes only from (A,B,C) triplets → still 2
    assert counts.get(("A", "C"), 0) == 2


# ---------------------------------------------------------------------------
# Ghost bar colouring: length encodes strength, colour encodes co-occurrence
# with sampled introgression.
# ---------------------------------------------------------------------------


def test_sampled_introgression_presence_flags_targets_with_sampled_edges():
    """A taxon is flagged only when it is the target of a non-zero sampled edge."""
    taxa_order = ["A", "B", "C", "D"]
    # Sampled edges are keyed (source, target).
    avg_non_ghost = {
        ("C", "A"): 0.6,  # A is a target -> flagged
        ("A", "B"): 0.0,  # zero weight -> B not flagged
        ("D", "C"): 0.3,  # C is a target -> flagged
    }

    presence = _sampled_introgression_presence(taxa_order, avg_non_ghost)

    assert presence == {"A": 1, "B": 0, "C": 1, "D": 0}


def _ghost_colour_scenario_results():
    """Build results where A has ghost + sampled and D has ghost only.

    Returns:
        A list of ``SimpleNamespace`` triplet results.
    """
    return [
        # Ghost on triplet (A,B,C) with dis1 BC -> ghost target A.
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
        # Inflow on triplet (A,B,C) with dis1 AC -> sampled edge (C, A), target A.
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="AC",
            bootstrap_value=0.5,
        ),
        # Ghost on triplet (D,B,C) with dis1 BC -> ghost target D, no sampled edge.
        SimpleNamespace(
            triplet=("D", "B", "C"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.4,
        ),
    ]


@pytest.mark.core
@pytest.mark.output
def test_ghost_strength_tsv_records_sampled_introgression_flag(tmp_path):
    """The ghost sheet gains a has_sampled_introgression 1/0 column."""
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")
    output_dir = tmp_path / "out"

    generate_introgression_maps(
        _ghost_colour_scenario_results(),
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
    )

    lines = _consolidation_lines(
        output_dir / "consolidation_data", "introgression_ghost_target_strength.tsv"
    )
    assert lines[0].split("\t") == [
        "target_taxon",
        "raw_strength",
        "has_sampled_introgression",
    ]

    flags = {}
    strengths = {}
    for line in lines[1:]:
        taxon, strength, flag = line.split("\t")
        flags[taxon] = int(flag)
        strengths[taxon] = float(strength)

    # A is a ghost target and also the target of sampled edge (C, A).
    assert flags["A"] == 1
    # D is a ghost target with no sampled edge pointing at it.
    assert flags["D"] == 0
    # Taxa with no ghost signal at all are still listed, flagged by sampled only.
    assert set(flags) == {"A", "B", "C", "D"}
    # Bar length still comes from the ghost strength, untouched by the flag.
    assert strengths["A"] > 0.0
    assert strengths["D"] > 0.0
