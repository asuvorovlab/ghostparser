"""Introgression map generation for GhostParser results.

This module builds a single combined figure from per-triplet introgression calls
containing a directed inflow/outflow heatmap (source x target) and a ghost
target-strength bar chart side by side, plus companion TSV artifacts for raw
bootstrap sums, supporting counts, and undiluted averages.
"""

from __future__ import annotations

import argparse as _argparse
from dataclasses import dataclass
from pathlib import Path

import dendropy
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from .config import prepare_output_directory


@dataclass(frozen=True)
class IntrogressionMapArtifacts:
    """Artifacts produced by introgression map generation."""

    plot_path: str
    non_ghost_matrix_tsv: str
    non_ghost_matrix_raw_sum_tsv: str
    non_ghost_matrix_supporting_count_tsv: str
    ghost_strength_tsv: str
    ghost_strength_raw_sum_tsv: str
    ghost_strength_supporting_count_tsv: str
    taxa_order_tsv: str
    non_sister_count_tsv: str
    taxa_count: int
    non_ghost_edge_count: int
    ghost_target_count: int


def _extract_result_fields(result):
    """Extract core fields from either dataclass results or dict-like rows."""
    if isinstance(result, dict):
        triplet = result.get("triplet")
        if isinstance(triplet, str):
            parts = [part.strip() for part in triplet.split(",")]
            triplet = tuple(parts) if len(parts) == 3 else None
        return (
            triplet,
            result.get("classification"),
            result.get("dis1_topology"),
            result.get("bootstrap_value"),
        )

    return (
        getattr(result, "triplet", None),
        getattr(result, "classification", None),
        getattr(result, "dis1_topology", None),
        getattr(result, "bootstrap_value", None),
    )


def _species_tree_taxa_order(species_tree_path, allowed_taxa=None):
    """Return taxa in species-tree traversal order."""
    tree = dendropy.Tree.get(
        path=str(species_tree_path), schema="newick", preserve_underscores=True
    )
    if allowed_taxa:
        tree.retain_taxa_with_labels(sorted({str(taxon) for taxon in allowed_taxa}))
    order = []
    seen = set()
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None or not leaf.taxon.label:
            continue
        label = str(leaf.taxon.label)
        if label in seen:
            continue
        seen.add(label)
        order.append(label)
    return order


def _load_species_tree(species_tree_path, allowed_taxa=None):
    tree = dendropy.Tree.get(
        path=str(species_tree_path), schema="newick", preserve_underscores=True
    )
    if allowed_taxa:
        tree.retain_taxa_with_labels(sorted({str(taxon) for taxon in allowed_taxa}))
    return tree


def _map_event(a_taxon, b_taxon, c_taxon, classification, dis1_topology):
    """Map a classified triplet to directed non-ghost edge or ghost target."""
    if not dis1_topology:
        return None, None

    topo = str(dis1_topology).strip().upper()
    if topo not in {"BC", "AC"}:
        return None, None

    if classification == "inflow_introgression":
        if topo == "BC":
            return (c_taxon, b_taxon), None
        return (c_taxon, a_taxon), None

    if classification == "outflow_introgression":
        if topo == "BC":
            return (b_taxon, c_taxon), None
        return (a_taxon, c_taxon), None

    if classification == "ghost_introgression":
        if topo == "BC":
            return None, a_taxon
        return None, b_taxon

    return None, None


def _collect_weights(results):
    """Aggregate directed non-ghost edge weights and ghost target weights."""
    non_ghost = {}
    ghost = {}
    taxa_seen = set()

    for result in results:
        triplet, classification, dis1_topology, bootstrap_value = (
            _extract_result_fields(result)
        )
        if not triplet or len(triplet) != 3:
            continue
        a_taxon, b_taxon, c_taxon = triplet
        taxa_seen.update((a_taxon, b_taxon, c_taxon))

        weight = float(bootstrap_value) if bootstrap_value is not None else 1.0
        if weight < 0:
            continue

        edge, ghost_target = _map_event(
            a_taxon, b_taxon, c_taxon, classification, dis1_topology
        )
        if edge is not None:
            non_ghost[edge] = non_ghost.get(edge, 0.0) + weight
        if ghost_target is not None:
            ghost[ghost_target] = ghost.get(ghost_target, 0.0) + weight

    return non_ghost, ghost, taxa_seen


def _build_taxa_order(species_order, taxa_seen):
    """Use species-tree order first, then append any missing taxa deterministically."""
    order = []
    seen = set()
    for taxon in species_order:
        if taxon in seen:
            continue
        seen.add(taxon)
        order.append(taxon)

    for taxon in sorted(taxa_seen):
        if taxon in seen:
            continue
        seen.add(taxon)
        order.append(taxon)

    return order


def _build_non_ghost_matrix(taxa_order, non_ghost_norm):
    """Create target x source matrix in species order."""
    matrix = []
    for target in taxa_order:
        row = []
        for source in taxa_order:
            row.append(non_ghost_norm.get((source, target), 0.0))
        matrix.append(row)
    return matrix


def _write_non_ghost_matrix_tsv(path, taxa_order, matrix):
    with open(path, "w") as out_f:
        out_f.write("target_taxon\t" + "\t".join(taxa_order) + "\n")
        for idx, target in enumerate(taxa_order):
            row_values = [f"{value:.12g}" for value in matrix[idx]]
            out_f.write(target + "\t" + "\t".join(row_values) + "\n")


def _write_ghost_strength_tsv(path, taxa_order, ghost_norm):
    with open(path, "w") as out_f:
        out_f.write("target_taxon\traw_strength\n")
        for target in taxa_order:
            out_f.write(f"{target}\t{ghost_norm.get(target, 0.0):.12g}\n")


def _write_single_value_tsv(path, taxa_order, values, header_name):
    with open(path, "w") as out_f:
        out_f.write(f"target_taxon\t{header_name}\n")
        for target in taxa_order:
            out_f.write(f"{target}\t{values.get(target, 0.0):.12g}\n")


def _write_taxa_order_tsv(path, taxa_order):
    with open(path, "w") as out_f:
        out_f.write("index\ttaxon\n")
        for idx, taxon in enumerate(taxa_order):
            out_f.write(f"{idx}\t{taxon}\n")


def _node_edge_length(node):
    if node.edge_length is None:
        return 0.0
    return float(node.edge_length)


def _species_tree_layout(species_tree_path, taxa_order, orientation):
    tree = _load_species_tree(species_tree_path, taxa_order)
    leaf_positions = {taxon: idx for idx, taxon in enumerate(taxa_order)}

    root_dist = {}

    def _record_root_dist(node, distance):
        root_dist[id(node)] = distance
        for child in node.child_node_iter():
            _record_root_dist(child, distance + _node_edge_length(child))

    _record_root_dist(tree.seed_node, 0.0)

    leaf_depths = []
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None or not leaf.taxon.label:
            continue
        label = str(leaf.taxon.label)
        if label in leaf_positions:
            leaf_depths.append(root_dist[id(leaf)])

    max_depth = max(leaf_depths) if leaf_depths else 1.0

    positions = {}
    leaf_nodes = []

    def _place(node):
        depth_value = root_dist[id(node)]
        if node.is_leaf():
            label = str(node.taxon.label) if node.taxon and node.taxon.label else None
            leaf_index = leaf_positions.get(label, 0)
            if orientation == "top":
                pos = (float(leaf_index), float(max_depth - depth_value))
            else:
                pos = (float(max_depth - depth_value), float(leaf_index))
            positions[id(node)] = pos
            leaf_nodes.append(node)
            return pos

        child_positions = [_place(child) for child in node.child_node_iter()]
        if orientation == "top":
            x = sum(position[0] for position in child_positions) / len(child_positions)
            y = float(max_depth - depth_value)
        else:
            x = float(max_depth - depth_value)
            y = sum(position[1] for position in child_positions) / len(child_positions)
        positions[id(node)] = (x, y)
        return positions[id(node)]

    _place(tree.seed_node)
    return tree, positions, max_depth, leaf_nodes


def _draw_species_tree_strip(
    ax, species_tree_path, taxa_order, orientation, show_leaf_labels=True
):
    tree, positions, max_depth, leaf_nodes = _species_tree_layout(
        species_tree_path, taxa_order, orientation
    )
    segments = []
    for node in tree.preorder_node_iter():
        parent_pos = positions[id(node)]
        for child in node.child_node_iter():
            child_pos = positions[id(child)]
            segments.append([parent_pos, child_pos])

    if segments:
        ax.add_collection(LineCollection(segments, colors="#1f1f1f", linewidths=1.0))

    if orientation == "top":
        for leaf in leaf_nodes:
            x_pos, y_pos = positions[id(leaf)]
            if y_pos > 0.0:
                ax.plot(
                    [x_pos, x_pos],
                    [y_pos, 0.0],
                    linestyle=":",
                    color="#666666",
                    linewidth=0.8,
                )
            if show_leaf_labels:
                label = str(leaf.taxon.label) if leaf.taxon and leaf.taxon.label else ""
                ax.text(
                    x_pos, -0.04, label, rotation=90, ha="center", va="top", fontsize=7
                )

        ax.set_xlim(-0.45, len(taxa_order) - 0.55)
        ax.set_ylim(-0.9 if show_leaf_labels else 0.0, max_depth + 0.95)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        # For left orientation we intentionally do not draw the tree (left strip removed)
        # Keep the axis empty and invisible.
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])

    ax.set_frame_on(False)
    ax.set_facecolor("none")


def _collect_counts(results):
    """Count triplets contributing to each average denominator.

    New behavior (undiluted consolidation):
    - For a directed non-ghost edge (source, target): count only triplets
      where the inference actually produced that directed edge (i.e. supporting
      classifications).
    - For a ghost target taxon: count only triplets classified as
      `ghost_introgression` for that taxon.
    Returns two dicts: `non_ghost_counts` and `ghost_counts` with supporting
    counts for each edge/taxon.
    """
    non_ghost_counts = {}
    ghost_counts = {}
    for result in results:
        triplet, classification, dis1_topology, _bootstrap_value = (
            _extract_result_fields(result)
        )
        if not triplet or len(triplet) != 3:
            continue
        a_taxon, b_taxon, c_taxon = triplet

        # Map the event for this row; only increment counts when the row
        # produced a supporting non-ghost edge or a ghost target.
        edge, ghost_target = _map_event(
            a_taxon, b_taxon, c_taxon, classification, dis1_topology
        )
        if edge is not None:
            non_ghost_counts[edge] = non_ghost_counts.get(edge, 0) + 1
        if ghost_target is not None:
            ghost_counts[ghost_target] = ghost_counts.get(ghost_target, 0) + 1

    return non_ghost_counts, ghost_counts


def _collect_non_sister_counts(results):
    """Count, for each unordered pair of taxa, how many triplets containing
    both have them as non-sister species.

    In every triplet ``(A, B, C)`` the convention is that A and B are sisters
    in the species tree, so:
    - The pair ``{A, B}`` is a *sister* pair in this triplet.
    - The pairs ``{A, C}`` and ``{B, C}`` are *non-sister* pairs.

    Returns a dict mapping canonical (sorted) 2-tuples of taxon names to
    integer counts.
    """
    non_sister_counts = {}
    for result in results:
        triplet, _classification, _dis1_topology, _bootstrap_value = (
            _extract_result_fields(result)
        )
        if not triplet or len(triplet) != 3:
            continue
        a_taxon, b_taxon, c_taxon = triplet

        # (A, B) are sisters — only the non-sister pairs get incremented.
        for sp1, sp2 in ((a_taxon, c_taxon), (b_taxon, c_taxon)):
            key = (sp1, sp2) if sp1 < sp2 else (sp2, sp1)
            non_sister_counts[key] = non_sister_counts.get(key, 0) + 1

    return non_sister_counts


def _build_non_sister_matrix(taxa_order, non_sister_counts):
    """Build a symmetric pairwise non-sister count matrix in taxa order."""
    matrix = []
    for row_taxon in taxa_order:
        row = []
        for col_taxon in taxa_order:
            if row_taxon == col_taxon:
                row.append(0)
            else:
                key = (
                    (row_taxon, col_taxon)
                    if row_taxon < col_taxon
                    else (col_taxon, row_taxon)
                )
                row.append(non_sister_counts.get(key, 0))
        matrix.append(row)
    return matrix


def _write_non_sister_matrix_tsv(path, taxa_order, matrix):
    with open(path, "w") as out_f:
        out_f.write("taxon\t" + "\t".join(taxa_order) + "\n")
        for idx, row_taxon in enumerate(taxa_order):
            row_values = [str(value) for value in matrix[idx]]
            out_f.write(row_taxon + "\t" + "\t".join(row_values) + "\n")


def _scaled_consolidation_text_sizes(n):
    """Return capped annotation font sizes for the consolidation plot."""
    scale_ref = max(min(n, 200), 1)
    frac = (scale_ref - 1) / 199.0 if scale_ref > 1 else 0.0
    eased_frac = frac**0.5

    def _interp(min_size, max_size):
        return int(round(min_size + (max_size - min_size) * eased_frac))

    return {
        "axis_label": _interp(13, 20),
        "cbar_label": _interp(12, 19),
        "cbar_tick": _interp(11, 14),
    }


def _plot_combined(path, species_tree_path, taxa_order, matrix_avg, ghost_avg):
    """Plot combined heatmap (sampled introgressions) with ghost bar chart to the right.

    Layout (left to right):
      heatmap (with species tree on top) | centered target labels | ghost bar | colorbar
    """
    from matplotlib.gridspec import GridSpec

    norm = Normalize(vmin=0.0, vmax=1.0)
    cmap = plt.get_cmap("PuBuGn")
    label_fontsize = 9
    n = len(taxa_order)
    text_sizes = _scaled_consolidation_text_sizes(n)

    # --- dynamic cell size (shrinks as n grows, same logic as reference) ---
    scale_ref = max(min(n, 200), 1)
    frac = (scale_ref - 1) / 199.0 if scale_ref > 1 else 0.0
    cell_size = 0.35 - (0.35 - 0.18) * frac

    heatmap_width = max(n * cell_size, 1.0)
    fig_height_heat = max(4.0, n * cell_size)
    tree_height = max(1.2, fig_height_heat * 0.22)
    fig_height = min(40.0, fig_height_heat + tree_height + 0.4)

    # --- measure longest target label and longest x-tick label in one pass ---
    tick_fontsize = max(5, label_fontsize - 1)
    temp_fig = plt.figure(figsize=(6, 2))
    temp_fig.canvas.draw()
    renderer = temp_fig.canvas.get_renderer()
    max_lw_px = 0  # width at label_fontsize  → used for side label panel
    max_tick_lw_px = (
        0  # width at tick_fontsize   → when rotated 90° this becomes height
    )
    for lbl in taxa_order:
        t = temp_fig.text(0, 0, lbl, fontsize=label_fontsize)
        max_lw_px = max(max_lw_px, t.get_window_extent(renderer=renderer).width)
        t.remove()
        t = temp_fig.text(0, 0, lbl, fontsize=tick_fontsize)
        max_tick_lw_px = max(
            max_tick_lw_px, t.get_window_extent(renderer=renderer).width
        )
        t.remove()
    dpi = temp_fig.dpi
    plt.close(temp_fig)

    pad_inches = (18 + 18) / 72.0
    label_panel_w = max(1.2, min(8.0, (max_lw_px / dpi) + pad_inches))
    bar_panel_w = 3.0
    cbar_w = 0.45
    fig_width = min(
        40.0, max(8.0, heatmap_width + label_panel_w + bar_panel_w + cbar_w + 1.0)
    )

    # --- measured label height drives the dedicated label-strip row ---
    # When rotated 90°, pixel-width of the longest label becomes the required row height.
    # Add explicit top and bottom padding so labels are not flush against adjacent rows.
    src_top_pad = (
        5.0 / 72.0
    )  # inches of whitespace above the text (gap from tree bottom)
    src_bot_pad = 9.0 / 72.0  # inches of whitespace below the text (gap to heatmap top)
    show_src_labels = n <= 120
    if show_src_labels:
        label_strip_h = max(0.4, max_tick_lw_px / dpi + src_top_pad + src_bot_pad)
    else:
        label_strip_h = 0.01  # effectively invisible strip for very large datasets

    # Compact tree strip (no leaf labels needed; label strip handles them)
    tree_h = max(1.0, fig_height_heat * 0.18)
    fig_height = min(40.0, tree_h + label_strip_h + fig_height_heat + 0.4)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=150)
    gs = GridSpec(
        nrows=3,
        ncols=4,
        height_ratios=[tree_h, label_strip_h, fig_height_heat],
        width_ratios=[heatmap_width, label_panel_w, bar_panel_w, cbar_w],
        hspace=0.0,
        wspace=0.02,
        figure=fig,
    )

    ax_tree = fig.add_subplot(gs[0, 0])
    ax_tree_right = fig.add_subplot(gs[0, 1:])
    ax_tree_right.axis("off")

    ax_src_labels = fig.add_subplot(gs[1, 0])  # dedicated source-taxon label strip
    ax_src_right = fig.add_subplot(gs[1, 1:])
    ax_src_right.axis("off")

    ax_heat = fig.add_subplot(gs[2, 0])
    ax_label = fig.add_subplot(gs[2, 1])
    ax_bar = fig.add_subplot(gs[2, 2])
    ax_cbar = fig.add_subplot(gs[2, 3])

    # --- species tree strip (no leaf labels — label strip below handles them) ---
    _draw_species_tree_strip(
        ax_tree, species_tree_path, taxa_order, "top", show_leaf_labels=False
    )

    # --- source-taxon label strip ---
    # x range matches seaborn's heatmap column centres: 0..n in data coords,
    # cell centre at i + 0.5.  y range [0, 1] (axis fraction = data coord here).
    # Text is anchored at y_text < 1.0 so there is src_top_pad inches of whitespace
    # between the tree bottom and the first character, and src_bot_pad inches of
    # whitespace between the last character and the heatmap top edge.
    ax_src_labels.set_xlim(0, n)
    ax_src_labels.set_ylim(0, 1)
    ax_src_labels.axis("off")
    if show_src_labels:
        src_fontsize = max(5, min(tick_fontsize, int(60 / max(n, 1) * 5 + 4)))
        y_text = 1.0 - src_top_pad / label_strip_h  # shift down by top-pad fraction
        for i, taxon in enumerate(taxa_order):
            ax_src_labels.text(
                i + 0.5,
                y_text,
                taxon,
                rotation=90,
                ha="center",
                va="top",
                fontsize=src_fontsize,
                fontstyle="italic",
                clip_on=False,
            )

    # --- heatmap ---
    sns.heatmap(
        matrix_avg,
        ax=ax_heat,
        cmap=cmap,
        norm=norm,
        xticklabels=False,
        yticklabels=False,
        linewidths=0.5,
        linecolor="white",
        cbar=False,
    )
    ax_heat.set_ylabel("Target taxon", fontsize=text_sizes["axis_label"])
    ax_heat.set_xlabel(
        "Sampled Introgression",
        fontweight="bold",
        fontsize=text_sizes["axis_label"],
    )
    ax_heat.tick_params(
        axis="x", bottom=False, labelbottom=False, top=False, labeltop=False
    )

    # --- centered target labels between heatmap and bar ---
    for spine in ax_label.spines.values():
        spine.set_visible(False)
    ax_label.set_xlim(0, 1)
    ax_label.set_xticks([])
    # align y with heatmap: seaborn sets ylim to [n, 0] with ticks at i+0.5
    ax_label.set_ylim(n, 0)
    ax_label.set_yticks([i + 0.5 for i in range(n)])
    ax_label.set_yticklabels(taxa_order, fontsize=label_fontsize)
    ax_label.tick_params(
        axis="y", left=False, right=False, labelleft=True, length=0, pad=0
    )
    for tick in ax_label.get_yticklabels():
        tick.set_horizontalalignment("center")
        tick.set_x(0.5)
        tick.set_fontstyle("italic")

    # --- ghost bar chart ---
    values = [ghost_avg.get(taxon, 0.0) for taxon in taxa_order]
    bar_colors = [cmap(norm(v)) for v in values]
    # bars at seaborn cell centers (i+0.5), height=0.8 to match cell boundaries
    ax_bar.barh(
        [i + 0.5 for i in range(n)],
        values,
        height=0.8,
        color=bar_colors,
        edgecolor="none",
    )
    ax_bar.set_xlim(0.0, 1.0)
    ax_bar.set_ylim(n, 0)
    ax_bar.set_xlabel(
        "Ghost Introgression",
        fontweight="bold",
        fontsize=text_sizes["axis_label"],
    )
    ax_bar.tick_params(
        axis="y", left=False, right=False, labelleft=False, labelright=False
    )
    for spine in ax_bar.spines.values():
        spine.set_visible(False)
    ax_bar.spines["bottom"].set_visible(True)

    # --- shared colorbar ---
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax_cbar)
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.ax.tick_params(labelsize=text_sizes["cbar_tick"])
    cbar.set_label(
        "Average Bootstrap Value",
        labelpad=8,
        fontsize=text_sizes["cbar_label"],
    )

    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def generate_introgression_maps(
    results,
    species_tree_path,
    output_dir,
    plot_taxa=None,
    outgroups=None,
    overwrite=True,
    reset_output_dir=True,
):
    """Generate non-ghost heatmap and ghost target-strength bar plot.

    Args:
        results: Iterable of TripletPipelineResult-like objects or dict rows.
        species_tree_path: Path to processed species tree used for ordering taxa.
        output_dir: Output directory for plots and tabular artifacts.
        plot_taxa: Optional iterable of taxa to retain in the plotted tree.
        outgroups: Optional iterable of taxon names to exclude from the plots
            (e.g. outgroup taxa used for rooting).  When ``None`` or empty no
            taxa are excluded.
        overwrite: When ``reset_output_dir`` is ``True``, whether to overwrite an
            existing output directory or write to an auto-suffixed sibling.
        reset_output_dir: When ``True`` (standalone use), (re)create a clean
            output directory for the plots. When ``False``, write into an
            existing, already-prepared run directory without resetting it, so a
            caller's other outputs (results TSV, processed trees, open metrics
            file) in that directory are never deleted.
    """
    if reset_output_dir:
        # Tests often place the species tree inside the requested output
        # directory. When overwrite=True, the reset would otherwise delete it.
        species_tree_source = Path(species_tree_path).expanduser().resolve()
        output_source = Path(output_dir).expanduser().resolve()
        preserved_species_tree_text: str | None = None
        if (
            overwrite
            and output_source.exists()
            and species_tree_source.exists()
            and species_tree_source.is_file()
            and output_source in species_tree_source.parents
        ):
            preserved_species_tree_text = species_tree_source.read_text(
                encoding="utf-8"
            )

        output_path = Path(prepare_output_directory(output_dir, overwrite=overwrite))

        if preserved_species_tree_text is not None:
            species_tree_source.parent.mkdir(parents=True, exist_ok=True)
            species_tree_source.write_text(
                preserved_species_tree_text, encoding="utf-8"
            )
    else:
        # Write into the caller's already-prepared run directory without
        # resetting it, so its existing outputs are preserved.
        output_path = Path(output_dir).expanduser().resolve()
        output_path.mkdir(parents=True, exist_ok=True)

    outgroup_set = set(outgroups) if outgroups else set()

    non_ghost_weights, ghost_weights, taxa_seen = _collect_weights(results)
    species_order = _species_tree_taxa_order(species_tree_path, plot_taxa)
    species_order = [t for t in species_order if t not in outgroup_set]
    taxa_order = _build_taxa_order(species_order, taxa_seen)
    taxa_order = [t for t in taxa_order if t not in outgroup_set]

    # compute occurrence counts so we can average bootstrap weights
    non_ghost_counts, ghost_counts = _collect_counts(results)

    # build averaged matrices (raw sums, supporting counts, and undiluted averages)
    avg_non_ghost = {}
    for edge, total in non_ghost_weights.items():
        cnt = non_ghost_counts.get(edge, 0)
        avg_non_ghost[edge] = (total / cnt) if cnt > 0 else 0.0

    avg_ghost = {}
    for taxon, total in ghost_weights.items():
        cnt = ghost_counts.get(taxon, 0)
        avg_ghost[taxon] = (total / cnt) if cnt > 0 else 0.0

    raw_non_ghost_matrix = _build_non_ghost_matrix(taxa_order, non_ghost_weights)
    supporting_non_ghost_matrix = _build_non_ghost_matrix(taxa_order, non_ghost_counts)
    avg_non_ghost_matrix = _build_non_ghost_matrix(taxa_order, avg_non_ghost)
    raw_ghost_values = {taxon: ghost_weights.get(taxon, 0.0) for taxon in taxa_order}
    supporting_ghost_values = {
        taxon: ghost_counts.get(taxon, 0) for taxon in taxa_order
    }

    combined_plot = output_path / "introgression_combined.png"

    # ensure consolidation_data subfolder for TSV artifacts
    consolidation_dir = output_path / "consolidation_data"
    consolidation_dir.mkdir(parents=True, exist_ok=True)

    non_ghost_tsv = consolidation_dir / "introgression_matrix_inflow_outflow.tsv"
    non_ghost_raw_sum_tsv = (
        consolidation_dir / "introgression_matrix_inflow_outflow_raw_sum.tsv"
    )
    non_ghost_supporting_count_tsv = (
        consolidation_dir / "introgression_matrix_inflow_outflow_supporting_count.tsv"
    )
    ghost_tsv = consolidation_dir / "introgression_ghost_target_strength.tsv"
    ghost_raw_sum_tsv = (
        consolidation_dir / "introgression_ghost_target_strength_raw_sum.tsv"
    )
    ghost_supporting_count_tsv = (
        consolidation_dir / "introgression_ghost_target_strength_supporting_count.tsv"
    )
    taxa_order_tsv = consolidation_dir / "introgression_taxa_order.tsv"
    non_sister_count_tsv = (
        consolidation_dir / "introgression_matrix_sampled_non_sister.tsv"
    )

    _write_non_ghost_matrix_tsv(non_ghost_raw_sum_tsv, taxa_order, raw_non_ghost_matrix)
    _write_non_ghost_matrix_tsv(
        non_ghost_supporting_count_tsv, taxa_order, supporting_non_ghost_matrix
    )
    non_ghost_tsv = consolidation_dir / "introgression_matrix_inflow_outflow.tsv"
    _write_non_ghost_matrix_tsv(non_ghost_tsv, taxa_order, avg_non_ghost_matrix)
    _write_single_value_tsv(ghost_raw_sum_tsv, taxa_order, raw_ghost_values, "raw_sum")
    _write_single_value_tsv(
        ghost_supporting_count_tsv,
        taxa_order,
        supporting_ghost_values,
        "supporting_count",
    )
    _write_ghost_strength_tsv(ghost_tsv, taxa_order, avg_ghost)
    _write_taxa_order_tsv(taxa_order_tsv, taxa_order)

    non_sister_counts = _collect_non_sister_counts(results)
    non_sister_matrix = _build_non_sister_matrix(taxa_order, non_sister_counts)
    _write_non_sister_matrix_tsv(non_sister_count_tsv, taxa_order, non_sister_matrix)

    _plot_combined(
        combined_plot, species_tree_path, taxa_order, avg_non_ghost_matrix, avg_ghost
    )

    return IntrogressionMapArtifacts(
        plot_path=str(combined_plot),
        non_ghost_matrix_tsv=str(non_ghost_tsv),
        non_ghost_matrix_raw_sum_tsv=str(non_ghost_raw_sum_tsv),
        non_ghost_matrix_supporting_count_tsv=str(non_ghost_supporting_count_tsv),
        ghost_strength_tsv=str(ghost_tsv),
        ghost_strength_raw_sum_tsv=str(ghost_raw_sum_tsv),
        ghost_strength_supporting_count_tsv=str(ghost_supporting_count_tsv),
        taxa_order_tsv=str(taxa_order_tsv),
        non_sister_count_tsv=str(non_sister_count_tsv),
        taxa_count=len(taxa_order),
        non_ghost_edge_count=len(non_ghost_weights),
        ghost_target_count=len(
            [taxon for taxon, value in ghost_weights.items() if value > 0.0]
        ),
    )


# ---------------------------------------------------------------------------
# Standalone CLI entry point
# ---------------------------------------------------------------------------


def _build_standalone_parser():
    parser = _argparse.ArgumentParser(
        description="Generate introgression maps from a GhostParser orchestrator results TSV."
    )
    parser.add_argument(
        "-r",
        "--results-tsv",
        required=True,
        help="Path to orchestrator_triplet_results.tsv",
    )
    parser.add_argument(
        "-st",
        "--species-tree-path",
        required=True,
        help="Path to the processed species tree (Newick)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory to write output plots and TSVs",
    )
    parser.add_argument(
        "--no-overwrite",
        dest="no_overwrite",
        action="store_true",
        default=False,
        help="Append a numeric suffix when the output directory already exists",
    )
    parser.add_argument(
        "-og",
        "--outgroups",
        default=None,
        help="Comma-separated outgroup taxon names to exclude from plots (e.g. 'Taxon1,Taxon2')",
    )
    return parser


def _parse_outgroups_arg(value):
    if not value:
        return None
    return [part.strip() for part in value.split(",") if part.strip()]


def _read_results_tsv(path):
    """Read orchestrator results TSV into a list of dicts."""
    rows = []
    with open(path, newline="") as fh:
        import csv

        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


if __name__ == "__main__":
    _parser = _build_standalone_parser()
    _args = _parser.parse_args()
    _outgroups = _parse_outgroups_arg(_args.outgroups)
    _results = _read_results_tsv(_args.results_tsv)
    _artifacts = generate_introgression_maps(
        _results,
        species_tree_path=_args.species_tree_path,
        output_dir=_args.output_dir,
        outgroups=_outgroups,
        overwrite=not _args.no_overwrite,
    )
    print(f"Taxa represented:          {_artifacts.taxa_count}")
    print(f"Non-ghost directed edges:  {_artifacts.non_ghost_edge_count}")
    print(f"Ghost targets with signal: {_artifacts.ghost_target_count}")
    print(f"Combined plot:             {_artifacts.plot_path}")
    print(f"Non-ghost matrix TSV:      {_artifacts.non_ghost_matrix_tsv}")
    print(f"Non-ghost raw sum TSV:     {_artifacts.non_ghost_matrix_raw_sum_tsv}")
    print(
        f"Non-ghost count TSV:       {_artifacts.non_ghost_matrix_supporting_count_tsv}"
    )
    print(f"Ghost strength TSV:        {_artifacts.ghost_strength_tsv}")
    print(f"Ghost raw sum TSV:         {_artifacts.ghost_strength_raw_sum_tsv}")
    print(
        f"Ghost count TSV:           {_artifacts.ghost_strength_supporting_count_tsv}"
    )
    print(f"Taxa order TSV:            {_artifacts.taxa_order_tsv}")
    print(f"Non-sister count TSV:      {_artifacts.non_sister_count_tsv}")
