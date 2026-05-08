"""Tests for introgression mapper module."""

from types import SimpleNamespace

from ghostparser.introgression_mapper import (
    _collect_counts,
    _draw_species_tree_strip,
    generate_introgression_maps,
)


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
    assert (tmp_path / "introgression_matrix_inflow_outflow.tsv").exists()
    assert (tmp_path / "introgression_ghost_target_strength.tsv").exists()
    assert (tmp_path / "introgression_taxa_order.tsv").exists()


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

    matrix_lines = (tmp_path / "introgression_matrix_inflow_outflow.tsv").read_text().splitlines()
    ghost_lines = (tmp_path / "introgression_ghost_target_strength.tsv").read_text().splitlines()

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

    matrix_lines = (tmp_path / "introgression_matrix_inflow_outflow.tsv").read_text().splitlines()
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

    matrix_lines = (tmp_path / "introgression_matrix_inflow_outflow.tsv").read_text().splitlines()
    ghost_lines = (tmp_path / "introgression_ghost_target_strength.tsv").read_text().splitlines()

    matrix_values = []
    for line in matrix_lines[1:]:
        matrix_values.extend(float(value) for value in line.split("\t")[1:])

    ghost_values = [float(line.split("\t")[1]) for line in ghost_lines[1:]]

    assert max(matrix_values) <= 1.0
    assert max(ghost_values) <= 1.0
    # non-ghost avg: sum(0.2+0.8)=1.0 / count(triplets containing both C and B)=2 → 0.5
    assert max(matrix_values) == 0.5
    # ghost avg: sum(0.1+0.9)=1.0 / count(triplets containing A)=4 (2 from ABC + 2 from ABD) → 0.25
    assert max(ghost_values) == 0.25


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

    matrix_lines = (tmp_path / "introgression_matrix_inflow_outflow.tsv").read_text().splitlines()
    ghost_lines = (tmp_path / "introgression_ghost_target_strength.tsv").read_text().splitlines()
    taxa_from_matrix = matrix_lines[0].split("\t")[1:]
    taxa_from_ghost = [line.split("\t")[0] for line in ghost_lines[1:]]

    assert "OG" not in taxa_from_matrix
    assert "OG" not in taxa_from_ghost
    assert artifacts.taxa_count == 3


def test_collect_counts_non_ghost_denominator_is_all_co_occurring_triplets():
    """_collect_counts denominator for a directed pair must be all triplets containing
    both taxa, regardless of how those triplets were classified."""
    results = [
        # classified — contributes weight for edge (C→B)
        SimpleNamespace(triplet=("A", "B", "C"), classification="inflow_introgression",
                        dis1_topology="BC", bootstrap_value=0.8),
        # no_introgression — still contains B and C, so must count
        SimpleNamespace(triplet=("A", "B", "C"), classification="no_introgression",
                        dis1_topology="BC", bootstrap_value=0.0),
        # unrelated triplet — does not contain C, must not count for (C→B)
        SimpleNamespace(triplet=("A", "B", "D"), classification="inflow_introgression",
                        dis1_topology="BC", bootstrap_value=0.5),
    ]
    non_ghost_counts, ghost_counts = _collect_counts(results)

    # (C, B) appears in triplet (A,B,C) twice — both rows contain B and C
    assert non_ghost_counts.get(("C", "B"), 0) == 2
    # (B, D) appears only in the ABD triplet
    assert non_ghost_counts.get(("B", "D"), 0) == 1


def test_collect_counts_ghost_denominator_is_all_triplets_containing_taxon():
    """_collect_counts denominator for a ghost target must be all triplets where that
    taxon appears in any position, regardless of classification."""
    results = [
        SimpleNamespace(triplet=("A", "B", "C"), classification="ghost_introgression",
                        dis1_topology="BC", bootstrap_value=0.9),
        SimpleNamespace(triplet=("A", "B", "C"), classification="no_introgression",
                        dis1_topology="BC", bootstrap_value=0.0),
        SimpleNamespace(triplet=("A", "C", "D"), classification="ghost_introgression",
                        dis1_topology="AC", bootstrap_value=0.5),
    ]
    _non_ghost_counts, ghost_counts = _collect_counts(results)

    # A appears in all 3 triplets
    assert ghost_counts.get("A", 0) == 3
    # C appears in all 3 triplets
    assert ghost_counts.get("C", 0) == 3
    # D appears only in triplet (A,C,D)
    assert ghost_counts.get("D", 0) == 1
    # B appears only in triplets (A,B,C) ×2
    assert ghost_counts.get("B", 0) == 2


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
        SimpleNamespace(triplet=("A", "B", "C"), classification="inflow_introgression",
                        dis1_topology="BC", bootstrap_value=0.6),
        SimpleNamespace(triplet=("A", "B", "C"), classification="no_introgression",
                        dis1_topology="BC", bootstrap_value=0.0),
        SimpleNamespace(triplet=("A", "B", "D"), classification="ghost_introgression",
                        dis1_topology="BC", bootstrap_value=0.8),
    ]

    generate_introgression_maps(
        results,
        species_tree_path=str(species_tree),
        output_dir=str(tmp_path),
    )

    matrix_lines = (tmp_path / "introgression_matrix_inflow_outflow.tsv").read_text().splitlines()
    ghost_lines = (tmp_path / "introgression_ghost_target_strength.tsv").read_text().splitlines()

    # Build lookup: matrix[target][source] = value
    header = matrix_lines[0].split("\t")[1:]
    matrix = {}
    for row_line in matrix_lines[1:]:
        parts = row_line.split("\t")
        row_target = parts[0]
        matrix[row_target] = {header[i]: float(parts[i + 1]) for i in range(len(header))}

    ghost_map = {}
    for line in ghost_lines[1:]:
        taxon, val = line.split("\t")
        ghost_map[taxon] = float(val)

    # 0.6 / 2 = 0.3
    assert abs(matrix["B"]["C"] - 0.3) < 1e-9
    # 0.8 / 3 ≈ 0.2667
    assert abs(ghost_map["A"] - 0.8 / 3) < 1e-9


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
    _draw_species_tree_strip(ax, str(species_tree), taxa_order, "top", show_leaf_labels=False)
    text_artists = [child for child in ax.get_children() if isinstance(child, Text)
                    and child.get_text().strip() in taxa_order]
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
    text_artists = [child for child in ax.get_children() if isinstance(child, Text)
                    and child.get_text().strip() in taxa_order]
    plt.close(fig)
    assert len(text_artists) == len(taxa_order)
