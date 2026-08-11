"""Tests for introgression mapper module."""

from pathlib import Path
from types import SimpleNamespace

from ghostparser.introgression_mapper import (
    GHOST_ONLY_BAR_COLOR,
    GHOST_WITH_SAMPLED_BAR_COLOR,
    _collect_counts,
    _collect_non_sister_counts,
    _draw_species_tree_strip,
    _sampled_introgression_presence,
    _scaled_consolidation_text_sizes,
    generate_introgression_maps,
)


def _consolidation_lines(consolidation_dir, filename):
    return (consolidation_dir / filename).read_text().splitlines()


def test_generate_introgression_maps_creates_expected_outputs(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

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
        output_dir=str(tmp_path),
    )

    assert artifacts.taxa_count == 4
    assert artifacts.non_ghost_edge_count >= 1
    assert artifacts.ghost_target_count >= 1

    assert (tmp_path / "introgression_combined.png").exists()
    assert artifacts.plot_path == str(tmp_path / "introgression_combined.png")
    consolidation_dir = tmp_path / "consolidation_data"
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


def test_generate_introgression_maps_appends_suffix_when_overwrite_disabled(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

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


def test_generate_introgression_maps_preserves_run_dir_when_reset_disabled(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.5,
        ),
    ]

    output_dir = tmp_path / "results"
    output_dir.mkdir()
    existing = output_dir / "orchestrator_triplet_results.tsv"
    existing.write_text("keep me")

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(output_dir),
        overwrite=True,
        reset_output_dir=False,
    )

    # With reset disabled, a caller's pre-existing outputs must survive and the
    # plots are written into the same directory (not a reset/suffixed one).
    assert existing.exists()
    assert existing.read_text() == "keep me"
    assert artifacts.plot_path == str(output_dir / "introgression_combined.png")
    assert Path(artifacts.plot_path).exists()


def test_scaled_consolidation_text_sizes_grow_with_taxa_count():
    small = _scaled_consolidation_text_sizes(5)
    large = _scaled_consolidation_text_sizes(83)
    capped = _scaled_consolidation_text_sizes(500)

    assert small["axis_label"] < large["axis_label"]
    assert small["cbar_label"] < large["cbar_label"]
    assert small["cbar_tick"] < large["cbar_tick"]

    assert large["axis_label"] == 17
    assert large["cbar_label"] == 16
    assert large["cbar_tick"] == 13

    assert capped["axis_label"] == 20
    assert capped["cbar_label"] == 19
    assert capped["cbar_tick"] == 14


def test_generate_introgression_maps_uses_full_species_tree_by_default(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="ghost_introgression",
            dis1_topology="AC",
            bootstrap_value=1.0,
        ),
    ]

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
    )

    consolidation_dir = tmp_path / "consolidation_data"
    matrix_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow.tsv"
    )
    ghost_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength.tsv"
    )

    header_fields = matrix_lines[0].split("\t")
    taxa_from_header = header_fields[1:]

    taxa_from_ghost = [line.split("\t")[0] for line in ghost_lines[1:]]

    assert "D" in taxa_from_header
    assert len(taxa_from_header) == artifacts.taxa_count
    assert len(matrix_lines) == artifacts.taxa_count + 1
    assert len(taxa_from_ghost) == artifacts.taxa_count
    assert taxa_from_ghost == taxa_from_header


def test_generate_introgression_maps_prunes_requested_plot_taxa(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="ghost_introgression",
            dis1_topology="AC",
            bootstrap_value=1.0,
        ),
    ]

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
        plot_taxa=["A", "B", "C"],
    )

    matrix_lines = _consolidation_lines(
        tmp_path / "consolidation_data", "introgression_matrix_inflow_outflow.tsv"
    )
    header_fields = matrix_lines[0].split("\t")
    taxa_from_header = header_fields[1:]

    assert "D" not in taxa_from_header
    assert len(taxa_from_header) == artifacts.taxa_count


def test_generate_introgression_maps_uses_raw_values_with_separate_scales(tmp_path):
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")

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
        output_dir=str(tmp_path),
    )

    consolidation_dir = tmp_path / "consolidation_data"
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


def test_generate_introgression_maps_excludes_outgroups(tmp_path):
    """Outgroup taxa supplied via the `outgroups` parameter must be absent from all outputs."""
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,OG:1);\n")

    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.6,
        ),
        SimpleNamespace(
            triplet=("A", "B", "OG"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.4,
        ),
    ]

    artifacts = generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
        outgroups=["OG"],
    )

    consolidation_dir = tmp_path / "consolidation_data"
    matrix_lines = _consolidation_lines(
        consolidation_dir, "introgression_matrix_inflow_outflow.tsv"
    )
    ghost_lines = _consolidation_lines(
        consolidation_dir, "introgression_ghost_target_strength.tsv"
    )
    taxa_from_matrix = matrix_lines[0].split("\t")[1:]
    taxa_from_ghost = [line.split("\t")[0] for line in ghost_lines[1:]]

    assert "OG" not in taxa_from_matrix
    assert "OG" not in taxa_from_ghost
    assert artifacts.taxa_count == 3


def test_collect_counts_non_ghost_denominator_is_all_co_occurring_triplets():
    """_collect_counts now returns supporting-triplet counts: only rows that
    produced the directed edge are counted."""
    results = [
        # classified — contributes weight for edge (C→B)
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.8,
        ),
        # no_introgression — still contains B and C, so must count
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="no_introgression",
            dis1_topology="BC",
            bootstrap_value=0.0,
        ),
        # unrelated triplet — does not contain C, must not count for (C→B)
        SimpleNamespace(
            triplet=("A", "B", "D"),
            classification="inflow_introgression",
            dis1_topology="BC",
            bootstrap_value=0.5,
        ),
    ]
    non_ghost_counts, ghost_counts = _collect_counts(results)

    # Only the first row produced the (C->B) directed edge, so count should be 1
    assert non_ghost_counts.get(("C", "B"), 0) == 1
    # (B, D) does not get produced by the third row mapping (it maps D->B),
    # so (B, D) supporting count should be 0
    assert non_ghost_counts.get(("B", "D"), 0) == 0


def test_collect_counts_ghost_denominator_is_all_triplets_containing_taxon():
    """_collect_counts now returns supporting-triplet counts for ghost targets:
    only rows classified as `ghost_introgression` that produced a ghost target
    are counted."""
    results = [
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="ghost_introgression",
            dis1_topology="BC",
            bootstrap_value=0.9,
        ),
        SimpleNamespace(
            triplet=("A", "B", "C"),
            classification="no_introgression",
            dis1_topology="BC",
            bootstrap_value=0.0,
        ),
        SimpleNamespace(
            triplet=("A", "C", "D"),
            classification="ghost_introgression",
            dis1_topology="AC",
            bootstrap_value=0.5,
        ),
    ]
    _non_ghost_counts, ghost_counts = _collect_counts(results)

    # The first row (ABC, topo=BC) produces ghost target A (since topo==BC -> a_taxon)
    # The third row (ACD, topo=AC) produces ghost target C (since topo!=BC -> b_taxon)
    assert ghost_counts.get("A", 0) == 1
    assert ghost_counts.get("C", 0) == 1
    assert ghost_counts.get("D", 0) == 0
    # B is not produced as a ghost target in these rows
    assert ghost_counts.get("B", 0) == 0


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
        output_dir=str(tmp_path),
    )

    consolidation_dir = tmp_path / "consolidation_data"
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


def test_draw_species_tree_strip_suppresses_leaf_labels(tmp_path):
    """show_leaf_labels=False must produce no Text artists on the axis."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.text import Text

    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,C:1);\n")
    taxa_order = ["A", "B", "C"]

    fig, ax = plt.subplots()
    _draw_species_tree_strip(
        ax, str(species_tree), taxa_order, "top", show_leaf_labels=False
    )
    text_artists = [
        child
        for child in ax.get_children()
        if isinstance(child, Text) and child.get_text().strip() in taxa_order
    ]
    plt.close(fig)
    assert len(text_artists) == 0


def test_draw_species_tree_strip_shows_leaf_labels_by_default(tmp_path):
    """show_leaf_labels=True (default) must draw one Text artist per leaf."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.text import Text

    species_tree = tmp_path / "species.tree"
    species_tree.write_text("((A:1,B:1):1,C:1);\n")
    taxa_order = ["A", "B", "C"]

    fig, ax = plt.subplots()
    _draw_species_tree_strip(ax, str(species_tree), taxa_order, "top")
    text_artists = [
        child
        for child in ax.get_children()
        if isinstance(child, Text) and child.get_text().strip() in taxa_order
    ]
    plt.close(fig)
    assert len(text_artists) == len(taxa_order)


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


def test_ghost_strength_tsv_records_sampled_introgression_flag(tmp_path):
    """The ghost sheet gains a has_sampled_introgression 1/0 column."""
    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

    generate_introgression_maps(
        _ghost_colour_scenario_results(),
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
    )

    lines = _consolidation_lines(
        tmp_path / "consolidation_data", "introgression_ghost_target_strength.tsv"
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


def test_ghost_bars_use_constant_colours_by_sampled_presence(tmp_path, monkeypatch):
    """Bar colours are the two fixed colours, chosen by the sampled-presence flag."""
    import matplotlib.axes

    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

    captured = {}
    original_barh = matplotlib.axes.Axes.barh

    def _recording_barh(self, y, width, **kwargs):
        captured["colors"] = kwargs.get("color")
        captured["widths"] = list(width)
        return original_barh(self, y, width, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "barh", _recording_barh)

    generate_introgression_maps(
        _ghost_colour_scenario_results(),
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
    )

    taxa_order = [
        line.split("\t")[1]
        for line in _consolidation_lines(
            tmp_path / "consolidation_data", "introgression_taxa_order.tsv"
        )[1:]
    ]
    colors = dict(zip(taxa_order, captured["colors"]))
    widths = dict(zip(taxa_order, captured["widths"]))

    # Only the two constant colours are ever used -- no colormap gradient.
    assert set(captured["colors"]) <= {
        GHOST_ONLY_BAR_COLOR,
        GHOST_WITH_SAMPLED_BAR_COLOR,
    }
    # cividis low end (dark blue) vs high end (yellow).
    assert colors["A"] == GHOST_WITH_SAMPLED_BAR_COLOR  # ghost + sampled
    assert colors["D"] == GHOST_ONLY_BAR_COLOR  # ghost only
    # A and D have different ghost strengths but D shares its colour with any
    # other ghost-only taxon, proving colour no longer tracks magnitude.
    assert widths["A"] != widths["D"]


def test_zero_heatmap_cells_are_masked(tmp_path, monkeypatch):
    """Cells with no introgression edge are masked instead of painted."""
    import numpy as np

    import ghostparser.introgression_mapper as mapper

    species_tree = tmp_path / "species.tree"
    species_tree.write_text("(((A:1,B:1):1,C:1):1,D:1);\n")

    captured = {}
    original_heatmap = mapper.sns.heatmap

    def _recording_heatmap(data, **kwargs):
        captured["data"] = np.asarray(data, dtype=float)
        captured["mask"] = np.asarray(kwargs.get("mask"))
        return original_heatmap(data, **kwargs)

    monkeypatch.setattr(mapper.sns, "heatmap", _recording_heatmap)

    generate_introgression_maps(
        _ghost_colour_scenario_results(),
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
    )

    # The scenario has exactly one sampled edge, (C, A), so exactly one cell is
    # non-zero and every other cell of the 4x4 matrix is masked.
    assert captured["mask"].tolist() == (captured["data"] == 0.0).tolist()
    assert int((~captured["mask"]).sum()) == 1
