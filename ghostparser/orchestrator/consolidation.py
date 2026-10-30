"""Consolidation stage: introgression maps from the orchestrator's results.

Builds two figures, each an inflow/outflow heatmap beside a ghost bar chart
(one shows bootstrap support, the other supporting triplet counts), plus
companion TSV matrices, consuming the in-memory results list.
"""

import textwrap
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import seaborn as sns
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator

from ghostparser.config import prepare_output_directory

from .trees import _copy_clade_for_taxa, _preorder, _terminals

BOOTSTRAP_PLOT_NAME = "introgression_bootstrap.png"
COUNTS_PLOT_NAME = "introgression_counts.png"

# Colormap shared by the bootstrap figure's heatmap, ghost bars and colorbar.
CONSOLIDATION_COLORMAP = "cividis"

# The counts figure carries its numbers as text, so every non-empty cell and
# every ghost bar takes one fill and no colour scale.
COUNT_FILL_COLOR = "black"
COUNT_TEXT_COLOR = "white"

# Ghost bar outlines: dashed when the target has ghost introgression only,
# solid when it is also the target of a sampled edge. The outline sits a small
# white gap outside the fill, so it stays readable on any fill, black included.
GHOST_OUTLINE_COLOR = "#1f1f1f"
GHOST_OUTLINE_WIDTH = 1.0
GHOST_ONLY_LINESTYLE = (0, (3.0, 2.0))
GHOST_WITH_SAMPLED_LINESTYLE = "solid"
GHOST_OUTLINE_GAP_POINTS = 2.0
# Fill height as a fraction of a heatmap row, leaving room for the outline.
GHOST_BAR_HEIGHT = 0.6

# Heatmap cells with no introgression edge. Lighter than every step of
# CONSOLIDATION_COLORMAP, so an empty cell never reads as a low-valued one.
EMPTY_CELL_COLOR = "#f0f0f0"
PANEL_BORDER_COLOR = "#bdbdbd"

# Vertical gap between the ghost bar legend and the uncoloured-cell note above
# it, in points so it holds at every figure size.
NOTE_LEGEND_GAP_POINTS = 18.0

HEATMAP_PANEL_TITLE = "Sampled Introgression"
GHOST_PANEL_TITLE = "Ghost Introgression"
COUNT_PANEL_SUFFIX = "\n(Supporting Triplet Count)"
# Slack around a panel title when it, rather than the data, sets panel width.
PANEL_TITLE_PAD = 0.25

# Count annotations: capped font size, and a bar longer than this fraction of
# its axis carries its count inside the end, where outside it would run past
# the panel edge.
MAX_ANNOTATION_FONTSIZE = 10.0
BAR_LABEL_INSIDE_FRACTION = 0.85
BAR_LABEL_PAD_POINTS = 3.0

# Figure margins reserved for artists that live outside the GridSpec: the
# rotated "Target taxon" label on the left, the colorbar's label on the right,
# and the bold panel titles underneath. GridSpec ratios divide the area left
# after these, so reserving the space explicitly is what makes a row's ratio
# equal its intended height in inches: matplotlib's default subplot margins
# would otherwise silently absorb roughly a fifth of every row.
# Rendering resolution, and the ceiling the figure is allowed to grow to. The
# figure is sized from the taxon count, so a large tree produces a large canvas
# rather than a squeezed one; the cap only stops a pathological input from
# exhausting memory. Agg refuses a dimension above 65536 px, which at
# FIGURE_DPI is about 436 inches, so this leaves ample headroom.
FIGURE_DPI = 150
MAX_FIGURE_IN = 200.0

# Floor on the heatmap panel itself, so a handful of taxa still get a readable
# panel rather than a sliver.
MIN_HEATMAP_IN = 4.0

# Height of the species-tree strip as a fraction of the heatmap panel.
TREE_STRIP_FRACTION = 0.18

MARGIN_LEFT_IN = 0.55
MARGIN_RIGHT_IN = 0.80
MARGIN_TOP_IN = 0.25
MARGIN_BOTTOM_IN = 0.70

TAXON_LABEL_FONTSIZE = 9
TAXON_TICK_FONTSIZE = 8
# Whitespace above and below the rotated source labels: the gap from the tree
# and the gap to the heatmap's top edge.
SOURCE_LABEL_TOP_PAD_IN = 5.0 / 72.0
SOURCE_LABEL_BOTTOM_PAD_IN = 9.0 / 72.0
COLORBAR_W_IN = 0.45


@dataclass(frozen=True)
class IntrogressionMapArtifacts:
    """Artifacts produced by introgression map generation."""

    bootstrap_plot_path: str
    counts_plot_path: str
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
    """Extract the fields the maps are built from off one triplet result.

    Args:
        result: A ``TripletPipelineResult``.

    Returns:
        A ``(triplet, classification, dis1_topology, bootstrap_value)`` tuple.
    """
    return (
        getattr(result, "triplet", None),
        getattr(result, "classification", None),
        getattr(result, "dis1_topology", None),
        getattr(result, "bootstrap_value", None),
    )


def _load_species_tree(species_tree_path, rename_map=None):
    """Read the processed species tree and put the display names on its leaves.

    Args:
        species_tree_path: Path to the processed species tree, in the labels
            the run measured in.
        rename_map: Optional mapping of those labels to the display names the
            results carry; applied to the leaves in memory only.

    Returns:
        The ``Bio.Phylo`` tree.
    """
    tree = Phylo.read(str(species_tree_path), "newick")
    if rename_map:
        for leaf in _terminals(tree.root):
            leaf.name = rename_map.get(leaf.name, leaf.name)
    return tree


def _species_tree_taxa_order(tree, allowed_taxa=None):
    """Return taxa in species-tree leaf order, restricted to ``allowed_taxa``.

    Pruning keeps the surviving leaves in their relative order, so filtering
    the traversal is the order the pruned tree would list.
    """
    allowed = None if not allowed_taxa else {str(taxon) for taxon in allowed_taxa}
    order = []
    seen = set()
    for leaf in _terminals(tree.root):
        if not leaf.name:
            continue
        label = str(leaf.name)
        if label in seen or (allowed is not None and label not in allowed):
            continue
        seen.add(label)
        order.append(label)
    return order


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
    """Write a target x source matrix TSV with the taxa as row and column labels."""
    with open(path, "w") as out_f:
        out_f.write("target_taxon\t" + "\t".join(taxa_order) + "\n")
        for idx, target in enumerate(taxa_order):
            row_values = [f"{value:.12g}" for value in matrix[idx]]
            out_f.write(target + "\t" + "\t".join(row_values) + "\n")


def _sampled_introgression_presence(taxa_order, non_ghost_norm):
    """Flag each taxon that is also the target of a sampled introgression edge.

    Sampled (non-ghost) edges are keyed ``(source, target)``, so a taxon counts
    as having sampled introgression when any source contributes a non-zero
    weight to it as the target, the same row the ghost bar chart draws.

    Args:
        taxa_order: Ordered taxa matching the plot axes.
        non_ghost_norm: Mapping of ``(source, target)`` to a weight/average.

    Returns:
        A dict mapping each taxon to ``1`` when it also has sampled
        introgression and ``0`` when its only signal is ghost introgression.
    """
    targets_with_sampled = {
        target
        for (_source, target), value in non_ghost_norm.items()
        if value  # non-zero weight
    }
    return {taxon: int(taxon in targets_with_sampled) for taxon in taxa_order}


def _write_ghost_strength_tsv(path, taxa_order, ghost_norm, sampled_presence=None):
    """Write per-taxon ghost strength plus the sampled-introgression flag.

    Args:
        path: Output TSV path.
        taxa_order: Ordered taxa matching the plot axes.
        ghost_norm: Mapping of taxon to its ghost strength.
        sampled_presence: Optional mapping of taxon to ``1``/``0`` indicating
            whether that taxon also has sampled introgression. Defaults to all
            zeros when omitted.
    """
    presence = sampled_presence or {}
    with open(path, "w") as out_f:
        out_f.write("target_taxon\traw_strength\thas_sampled_introgression\n")
        for target in taxa_order:
            strength = ghost_norm.get(target, 0.0)
            flag = int(presence.get(target, 0))
            out_f.write(f"{target}\t{strength:.12g}\t{flag}\n")


def _write_single_value_tsv(path, taxa_order, values, header_name):
    """Write a per-taxon TSV with one value column."""
    with open(path, "w") as out_f:
        out_f.write(f"target_taxon\t{header_name}\n")
        for target in taxa_order:
            out_f.write(f"{target}\t{values.get(target, 0.0):.12g}\n")


def _write_taxa_order_tsv(path, taxa_order):
    """Write the taxa order as an index/taxon TSV."""
    with open(path, "w") as out_f:
        out_f.write("index\ttaxon\n")
        for idx, taxon in enumerate(taxa_order):
            out_f.write(f"{idx}\t{taxon}\n")


def _node_edge_length(clade):
    """Return a clade's branch length, ``0.0`` when it has none."""
    if clade.branch_length is None:
        return 0.0
    return float(clade.branch_length)


def _species_tree_layout(tree, taxa_order, orientation):
    """Lay out the species tree along one axis for the plot strip.

    Args:
        tree: The ``Bio.Phylo`` species tree; a pruned copy keeps only
            ``taxa_order``.
        taxa_order: The leaf order along the axis.
        orientation: ``horizontal`` or ``vertical``.

    Returns:
        ``(tree, positions, max_depth, leaf_nodes)``: node positions along and
        across the axis, the deepest root distance, and the leaves.
    """
    if taxa_order:
        pruned_root = _copy_clade_for_taxa(tree.root, set(taxa_order))
        if pruned_root is not None:
            tree = Tree(root=pruned_root, rooted=True)
    leaf_positions = {taxon: idx for idx, taxon in enumerate(taxa_order)}

    # A tree without branch lengths is drawn as a cladogram, one unit per
    # branch, rather than collapsing onto its root.
    has_lengths = any(
        clade.branch_length is not None for clade in _preorder(tree.root)[1:]
    )
    root_dist = {}

    def _record_root_dist(clade, distance):
        root_dist[id(clade)] = distance
        for child in clade.clades:
            step = _node_edge_length(child) if has_lengths else 1.0
            _record_root_dist(child, distance + step)

    _record_root_dist(tree.root, 0.0)

    leaf_depths = []
    for leaf in _terminals(tree.root):
        if not leaf.name:
            continue
        label = str(leaf.name)
        if label in leaf_positions:
            leaf_depths.append(root_dist[id(leaf)])

    max_depth = max(leaf_depths) if leaf_depths else 1.0

    positions = {}
    leaf_nodes = []

    def _place(node):
        depth_value = root_dist[id(node)]
        if not node.clades:
            label = str(node.name) if node.name else None
            leaf_index = leaf_positions.get(label, 0)
            if orientation == "top":
                pos = (float(leaf_index), float(max_depth - depth_value))
            else:
                pos = (float(max_depth - depth_value), float(leaf_index))
            positions[id(node)] = pos
            leaf_nodes.append(node)
            return pos

        child_positions = [_place(child) for child in node.clades]
        if orientation == "top":
            x = sum(position[0] for position in child_positions) / len(child_positions)
            y = float(max_depth - depth_value)
        else:
            x = float(max_depth - depth_value)
            y = sum(position[1] for position in child_positions) / len(child_positions)
        positions[id(node)] = (x, y)
        return positions[id(node)]

    _place(tree.root)
    return tree, positions, max_depth, leaf_nodes


def _draw_species_tree_strip(
    ax, species_tree, taxa_order, orientation, show_leaf_labels=True
):
    """Draw the species tree as a thin strip aligned to the heatmap axis."""
    tree, positions, max_depth, leaf_nodes = _species_tree_layout(
        species_tree, taxa_order, orientation
    )
    segments = []
    for node in _preorder(tree.root):
        parent_pos = positions[id(node)]
        for child in node.clades:
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
                label = str(leaf.name) if leaf.name else ""
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
    """Count the triplets that produced each edge, for the average denominators.

    A directed non-ghost edge counts only the triplets whose inference produced
    that edge; a ghost target counts only the triplets classified
    ``ghost_introgression`` for it.

    Args:
        results: Iterable of ``TripletPipelineResult`` objects.

    Returns:
        A tuple ``(non_ghost_counts, ghost_counts)`` keyed by directed
        ``(source, target)`` pair and by taxon respectively.
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
    """Count the triplets in which each unordered taxon pair is non-sister.

    Triplets are normalized so A and B are the species-tree sisters, leaving
    ``{A, C}`` and ``{B, C}`` as the non-sister pairs.

    Args:
        results: Iterable of ``TripletPipelineResult`` objects.

    Returns:
        A dict mapping sorted 2-tuples of taxon names to counts.
    """
    non_sister_counts = {}
    for result in results:
        triplet, _classification, _dis1_topology, _bootstrap_value = (
            _extract_result_fields(result)
        )
        if not triplet or len(triplet) != 3:
            continue
        a_taxon, b_taxon, c_taxon = triplet

        # (A, B) are sisters, only the non-sister pairs get incremented.
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
    """Write the symmetric non-sister count matrix TSV."""
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


@dataclass(frozen=True)
class _FigureLayout:
    """Panel sizes in inches, plus the font sizes that scale with the taxa count."""

    fig_width: float
    fig_height: float
    heatmap_width: float
    label_panel_w: float
    bar_panel_w: float
    cbar_w: float
    tree_h: float
    label_strip_h: float
    heatmap_h: float
    text_sizes: dict


def _figure_layout(
    taxa_order, sampled_title, ghost_title, min_bar_panel_w, has_colorbar
):
    """Size every panel from the taxa count and the measured label and title widths.

    Args:
        taxa_order: Ordered taxa for both panels.
        sampled_title: The heatmap's panel title.
        ghost_title: The ghost bar chart's panel title.
        min_bar_panel_w: Floor on the ghost panel's width, in inches.
        has_colorbar: Whether to reserve a column for the colorbar.

    Returns:
        A ``_FigureLayout``.
    """
    n = len(taxa_order)
    text_sizes = _scaled_consolidation_text_sizes(n)

    # Cells shrink as the taxa count grows.
    scale_ref = max(min(n, 200), 1)
    frac = (scale_ref - 1) / 199.0 if scale_ref > 1 else 0.0
    cell_size = 0.35 - (0.35 - 0.18) * frac
    heatmap_h = max(MIN_HEATMAP_IN, n * cell_size)

    temp_fig = plt.figure(figsize=(6, 2))
    temp_fig.canvas.draw()
    renderer = temp_fig.canvas.get_renderer()

    def _width_in(text, **kwargs):
        artist = temp_fig.text(0, 0, text, **kwargs)
        width = artist.get_window_extent(renderer=renderer).width
        artist.remove()
        return width / temp_fig.dpi

    longest_label = max(
        (_width_in(t, fontsize=TAXON_LABEL_FONTSIZE) for t in taxa_order), default=0.0
    )
    longest_tick = max(
        (_width_in(t, fontsize=TAXON_TICK_FONTSIZE) for t in taxa_order), default=0.0
    )
    # The bold panel titles sit outside their axes and are not picked up by
    # `bbox_inches="tight"`, so each panel has to be wide enough to hold its own
    # title. This only ever binds for small taxa counts, where the heatmap is
    # narrower than the words underneath it.
    sampled_title_w, ghost_title_w = (
        _width_in(title, fontsize=text_sizes["axis_label"], fontweight="bold")
        for title in (sampled_title, ghost_title)
    )
    plt.close(temp_fig)

    heatmap_width = max(n * cell_size, 1.0, sampled_title_w + PANEL_TITLE_PAD)
    label_panel_w = max(1.2, min(8.0, longest_label + 36 / 72.0))
    bar_panel_w = max(min_bar_panel_w, ghost_title_w + PANEL_TITLE_PAD)
    cbar_w = COLORBAR_W_IN if has_colorbar else 0.0
    grid_width = heatmap_width + label_panel_w + bar_panel_w + cbar_w
    fig_width = min(
        MAX_FIGURE_IN, max(8.0, grid_width + MARGIN_LEFT_IN + MARGIN_RIGHT_IN)
    )

    # Rotated 90 degrees, the longest source label's width becomes the strip's
    # height; the pads keep it clear of the tree above and the heatmap below.
    label_strip_h = max(
        0.4, longest_tick + SOURCE_LABEL_TOP_PAD_IN + SOURCE_LABEL_BOTTOM_PAD_IN
    )
    tree_h = max(1.0, heatmap_h * TREE_STRIP_FRACTION)

    # If the canvas would exceed the cap, the heatmap absorbs the shortfall and
    # the label strip keeps the inches it was measured to need. Scaling every row
    # down together would shrink the strip too and clip the labels it exists to
    # hold, which is the one thing the layout must not do.
    fixed_h = tree_h + label_strip_h + MARGIN_TOP_IN + MARGIN_BOTTOM_IN
    if heatmap_h + fixed_h > MAX_FIGURE_IN:
        heatmap_h = max(MIN_HEATMAP_IN, MAX_FIGURE_IN - fixed_h)
        tree_h = max(1.0, heatmap_h * TREE_STRIP_FRACTION)

    grid_height = tree_h + label_strip_h + heatmap_h
    fig_height = min(MAX_FIGURE_IN, grid_height + MARGIN_TOP_IN + MARGIN_BOTTOM_IN)

    return _FigureLayout(
        fig_width=fig_width,
        fig_height=fig_height,
        heatmap_width=heatmap_width,
        label_panel_w=label_panel_w,
        bar_panel_w=bar_panel_w,
        cbar_w=cbar_w,
        tree_h=tree_h,
        label_strip_h=label_strip_h,
        heatmap_h=heatmap_h,
        text_sizes=text_sizes,
    )


def _create_figure(layout, species_tree, taxa_order):
    """Build the figure grid and draw the panels both figures share.

    The species tree runs along the top, the source labels under it and the
    target labels between the heatmap and the ghost bars; the heatmap, bars and
    colorbar are left to the caller.

    Args:
        layout: The ``_FigureLayout``.
        species_tree: Processed species tree, pruned to ``taxa_order`` on a copy.
        taxa_order: Ordered taxa for both panels.

    Returns:
        ``(fig, axes)``, with ``axes`` a dict keyed by panel name; ``cbar`` is
        ``None`` when the layout reserves no colorbar.
    """
    width_ratios = [layout.heatmap_width, layout.label_panel_w, layout.bar_panel_w]
    if layout.cbar_w:
        width_ratios.append(layout.cbar_w)
    fig = plt.figure(figsize=(layout.fig_width, layout.fig_height), dpi=FIGURE_DPI)
    gs = GridSpec(
        nrows=3,
        ncols=len(width_ratios),
        height_ratios=[layout.tree_h, layout.label_strip_h, layout.heatmap_h],
        width_ratios=width_ratios,
        hspace=0.0,
        wspace=0.02,
        # Pin the grid to the reserved margins so each ratio keeps its intended
        # size in inches instead of being scaled into matplotlib's defaults.
        left=MARGIN_LEFT_IN / layout.fig_width,
        right=1.0 - MARGIN_RIGHT_IN / layout.fig_width,
        top=1.0 - MARGIN_TOP_IN / layout.fig_height,
        bottom=MARGIN_BOTTOM_IN / layout.fig_height,
        figure=fig,
    )
    axes = {
        "tree": fig.add_subplot(gs[0, 0]),
        "tree_right": fig.add_subplot(gs[0, 1:]),
        "src_labels": fig.add_subplot(gs[1, 0]),
        "src_right": fig.add_subplot(gs[1, 1:]),
        "heat": fig.add_subplot(gs[2, 0]),
        "target_labels": fig.add_subplot(gs[2, 1]),
        "bar": fig.add_subplot(gs[2, 2]),
        "cbar": fig.add_subplot(gs[2, 3]) if layout.cbar_w else None,
    }
    axes["tree_right"].axis("off")
    axes["src_right"].axis("off")

    _draw_species_tree_strip(
        axes["tree"], species_tree, taxa_order, "top", show_leaf_labels=False
    )
    _draw_source_labels(axes["src_labels"], taxa_order, layout)
    _draw_target_labels(axes["target_labels"], taxa_order)
    return fig, axes


def _draw_source_labels(ax, taxa_order, layout):
    """Write the source-taxon labels, rotated, between the tree and the heatmap."""
    n = len(taxa_order)
    # Column centres at i + 0.5, matching seaborn's heatmap; y is the axis
    # fraction, with the text hung `SOURCE_LABEL_TOP_PAD_IN` below the tree.
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.axis("off")
    src_fontsize = max(5, min(TAXON_TICK_FONTSIZE, int(60 / max(n, 1) * 5 + 4)))
    y_text = 1.0 - SOURCE_LABEL_TOP_PAD_IN / layout.label_strip_h
    for i, taxon in enumerate(taxa_order):
        ax.text(
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


def _draw_target_labels(ax, taxa_order):
    """Write the target-taxon labels, centred between the heatmap and the bars."""
    n = len(taxa_order)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    # Seaborn puts the heatmap's rows at i + 0.5 on a [n, 0] y range.
    ax.set_ylim(n, 0)
    ax.set_yticks([i + 0.5 for i in range(n)])
    ax.set_yticklabels(taxa_order, fontsize=TAXON_LABEL_FONTSIZE)
    ax.tick_params(axis="y", left=False, right=False, labelleft=True, length=0, pad=0)
    for tick in ax.get_yticklabels():
        tick.set_horizontalalignment("center")
        tick.set_x(0.5)
        tick.set_fontstyle("italic")


def _draw_heatmap(ax, values, cmap, norm, title, layout):
    """Draw the sampled-introgression heatmap, leaving zero cells uncoloured."""
    # Zero cells mean "no introgression edge observed". Masking them keeps the
    # sparse real signal legible: they fall through to the axes facecolor, and
    # the white gridlines plus panel border keep the row/column structure.
    ax.set_facecolor(EMPTY_CELL_COLOR)
    with warnings.catch_warnings():
        # A run with no sampled introgression masks every cell, and seaborn's
        # autoscale then takes nanmin/nanmax of an all-NaN array. The limits it
        # computes are unused because `norm` fixes them, so the warning reports
        # a calculation that cannot affect the output.
        warnings.filterwarnings("ignore", "All-NaN slice encountered", RuntimeWarning)
        sns.heatmap(
            values,
            ax=ax,
            cmap=cmap,
            norm=norm,
            mask=values == 0.0,
            xticklabels=False,
            yticklabels=False,
            linewidths=0.5,
            linecolor="white",
            cbar=False,
        )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color(PANEL_BORDER_COLOR)
    ax.set_ylabel("Target taxon", fontsize=layout.text_sizes["axis_label"])
    ax.set_xlabel(title, fontweight="bold", fontsize=layout.text_sizes["axis_label"])
    ax.tick_params(axis="x", bottom=False, labelbottom=False, top=False, labeltop=False)


def _annotate_heatmap_counts(ax, values, layout):
    """Write each non-empty cell's count on it, sized to fit the cell."""
    rows, cols = np.nonzero(values)
    if rows.size == 0:
        return
    labels = [str(int(values[r, c])) for r, c in zip(rows, cols)]
    n = len(values)
    cell_w_pt = layout.heatmap_width * 72.0 / n
    cell_h_pt = layout.heatmap_h * 72.0 / n
    # Digits run about 0.6 em wide.
    longest = max(len(label) for label in labels)
    fontsize = min(
        MAX_ANNOTATION_FONTSIZE, 0.5 * cell_h_pt, 0.85 * cell_w_pt / (0.6 * longest)
    )
    for r, c, label in zip(rows, cols, labels):
        ax.text(
            c + 0.5,
            r + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=fontsize,
            color=COUNT_TEXT_COLOR,
        )


def _draw_ghost_bars(fig, ax, lengths, fill_colors, has_sampled, x_limits):
    """Draw one ghost bar per target, outlined dashed (ghost only) or solid (+ sampled).

    Args:
        fig: The figure, for converting the outline gap from points.
        ax: The bar panel.
        lengths: Bar length per taxon in ``taxa_order``; ``0`` draws no bar.
        fill_colors: Fill colour per taxon.
        has_sampled: Per taxon, whether it is also the target of a sampled edge.
        x_limits: ``(x_min, x_max)`` of the panel.

    Returns:
        The outline gap in points, so labels can clear the outline.
    """
    n = len(lengths)
    ax.set_xlim(*x_limits)
    ax.set_ylim(n, 0)
    ax.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # The gap is set in points so it reads the same at any figure size, then
    # converted per axis; it shrinks on crowded rows so outlines never touch.
    position = ax.get_position()
    width_in = position.width * fig.get_figwidth()
    height_in = position.height * fig.get_figheight()
    row_h_pt = height_in * 72.0 / max(n, 1)
    gap_pt = min(GHOST_OUTLINE_GAP_POINTS, 0.1 * row_h_pt)
    gap_x = (gap_pt / 72.0) * (x_limits[1] - x_limits[0]) / width_in
    gap_y = (gap_pt / 72.0) * n / height_in

    for i, (length, fill, sampled) in enumerate(zip(lengths, fill_colors, has_sampled)):
        if not length:
            continue
        y0 = i + 0.5 - GHOST_BAR_HEIGHT / 2.0
        ax.add_patch(
            Rectangle((0.0, y0), length, GHOST_BAR_HEIGHT, facecolor=fill, linewidth=0)
        )
        ax.add_patch(
            Rectangle(
                (-gap_x, y0 - gap_y),
                length + 2.0 * gap_x,
                GHOST_BAR_HEIGHT + 2.0 * gap_y,
                fill=False,
                edgecolor=GHOST_OUTLINE_COLOR,
                linewidth=GHOST_OUTLINE_WIDTH,
                linestyle=(
                    GHOST_WITH_SAMPLED_LINESTYLE if sampled else GHOST_ONLY_LINESTYLE
                ),
                clip_on=False,
            )
        )
    return gap_pt


def _annotate_ghost_counts(ax, counts, x_max, gap_pt, layout):
    """Write each ghost bar's count at its end: outside in black, or inside a long bar."""
    row_h_pt = layout.heatmap_h * 72.0 / max(len(counts), 1)
    fontsize = min(MAX_ANNOTATION_FONTSIZE, 0.5 * row_h_pt)
    for i, count in enumerate(counts):
        if not count:
            continue
        inside = count > BAR_LABEL_INSIDE_FRACTION * x_max
        offset = BAR_LABEL_PAD_POINTS if inside else BAR_LABEL_PAD_POINTS + gap_pt
        ax.annotate(
            str(int(count)),
            xy=(count, i + 0.5),
            xytext=(-offset if inside else offset, 0),
            textcoords="offset points",
            ha="right" if inside else "left",
            va="center",
            fontsize=fontsize,
            color=COUNT_TEXT_COLOR if inside else "black",
        )


def _draw_ghost_legend_and_note(fig, axes, layout, note):
    """Draw the outline legend above the ghost bars, and ``note`` above the legend."""
    outline = {
        "facecolor": "white",
        "edgecolor": GHOST_OUTLINE_COLOR,
        "linewidth": GHOST_OUTLINE_WIDTH,
    }
    handles = [
        Patch(linestyle=GHOST_ONLY_LINESTYLE, label="Ghost only", **outline),
        Patch(
            linestyle=GHOST_WITH_SAMPLED_LINESTYLE, label="Ghost + sampled", **outline
        ),
    ]
    # Anchored to the bottom of the source-label strip, the row directly above
    # the bars, it grows upward into the tree-strip row, empty on this side.
    legend = axes["src_right"].legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        frameon=False,
        fontsize=layout.text_sizes["cbar_tick"],
        title="Ghost Introgression\nBar Outline",
        title_fontsize=layout.text_sizes["cbar_label"],
        borderpad=0.8,
    )
    # The title wraps onto two lines; centre them on each other.
    legend.get_title().set_multialignment("center")

    # The note is wrapped to the bar panel's width and sits on the legend's
    # rendered top edge, hence the draw to measure the legend.
    note_font = layout.text_sizes["cbar_tick"]
    wrap_chars = max(18, int((layout.bar_panel_w * 72.0) / (0.5 * note_font)))
    try:
        renderer = fig.canvas.get_renderer()
    except AttributeError:  # backend without a cached renderer
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    legend_top = (
        legend.get_window_extent(renderer).transformed(fig.transFigure.inverted()).y1
    )
    # The gap is in points so the padding stays visually constant; a figure
    # fraction would shrink as `fig_height` grows with the taxa count.
    note_gap = (NOTE_LEGEND_GAP_POINTS / 72.0) / layout.fig_height
    bar_position = axes["bar"].get_position()
    fig.text(
        bar_position.x0 + bar_position.width / 2.0,
        legend_top + note_gap,
        textwrap.fill(note, width=wrap_chars),
        ha="center",
        va="bottom",
        fontsize=note_font,
        fontweight="bold",
        color="#4d4d4d",
    )


def _plot_bootstrap_map(
    path, species_tree, taxa_order, matrix_avg, ghost_avg, ghost_counts, ghost_has_sampled
):
    """Plot average bootstrap support: the sampled heatmap beside equal-length ghost bars.

    Heatmap cells and ghost bars share one colour scale, so the bar colour, not
    its length, carries the ghost support.

    Args:
        path: Output image path.
        species_tree: Processed species tree (``Bio.Phylo``) drawn as the top
            strip, pruned to ``taxa_order`` on a copy.
        taxa_order: Ordered taxa for both panels.
        matrix_avg: Target x source matrix of sampled-introgression averages.
        ghost_avg: Mapping of taxon to its average ghost bootstrap value.
        ghost_counts: Mapping of taxon to its ghost supporting count; a taxon
            gets a bar when it is non-zero.
        ghost_has_sampled: Mapping of taxon to ``1``/``0``, whether it is also
            the target of a sampled edge.
    """
    cmap = plt.get_cmap(CONSOLIDATION_COLORMAP)
    norm = Normalize(vmin=0.0, vmax=1.0)
    layout = _figure_layout(
        taxa_order,
        HEATMAP_PANEL_TITLE,
        GHOST_PANEL_TITLE,
        min_bar_panel_w=1.5,
        has_colorbar=True,
    )
    fig, axes = _create_figure(layout, species_tree, taxa_order)

    _draw_heatmap(
        axes["heat"],
        np.asarray(matrix_avg, dtype=float),
        cmap,
        norm,
        HEATMAP_PANEL_TITLE,
        layout,
    )
    _draw_ghost_bars(
        fig,
        axes["bar"],
        [1.0 if ghost_counts.get(taxon) else 0.0 for taxon in taxa_order],
        [cmap(norm(ghost_avg.get(taxon, 0.0))) for taxon in taxa_order],
        [bool(ghost_has_sampled.get(taxon)) for taxon in taxa_order],
        x_limits=(-0.05, 1.05),
    )
    # Equal-length bars carry no scale, so the panel keeps only its title.
    axes["bar"].set_xticks([])
    axes["bar"].set_xlabel(
        GHOST_PANEL_TITLE,
        fontweight="bold",
        fontsize=layout.text_sizes["axis_label"],
    )
    _draw_ghost_legend_and_note(
        fig,
        axes,
        layout,
        "Heatmap cells with no introgression are left uncoloured and are not "
        "represented on the colour scale.",
    )

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=axes["cbar"])
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0", "0.5", "1"])
    cbar.ax.tick_params(labelsize=layout.text_sizes["cbar_tick"])
    cbar.set_label(
        "Average Bootstrap Value",
        labelpad=8,
        fontsize=layout.text_sizes["cbar_label"],
    )

    fig.savefig(path, bbox_inches="tight", dpi=FIGURE_DPI)
    plt.close(fig)


def _plot_counts_map(
    path, species_tree, taxa_order, matrix_counts, ghost_counts, ghost_has_sampled
):
    """Plot supporting triplet counts: count-labelled heatmap cells beside ghost bars as long as their counts.

    Args:
        path: Output image path.
        species_tree: Processed species tree (``Bio.Phylo``) drawn as the top
            strip, pruned to ``taxa_order`` on a copy.
        taxa_order: Ordered taxa for both panels.
        matrix_counts: Target x source matrix of sampled supporting counts.
        ghost_counts: Mapping of taxon to its ghost supporting count.
        ghost_has_sampled: Mapping of taxon to ``1``/``0``, whether it is also
            the target of a sampled edge.
    """
    sampled_title = HEATMAP_PANEL_TITLE + COUNT_PANEL_SUFFIX
    ghost_title = GHOST_PANEL_TITLE + COUNT_PANEL_SUFFIX
    layout = _figure_layout(
        taxa_order, sampled_title, ghost_title, min_bar_panel_w=3.0, has_colorbar=False
    )
    fig, axes = _create_figure(layout, species_tree, taxa_order)

    heat_values = np.asarray(matrix_counts, dtype=float)
    _draw_heatmap(
        axes["heat"],
        heat_values,
        ListedColormap([COUNT_FILL_COLOR]),
        Normalize(vmin=0.0, vmax=1.0),
        sampled_title,
        layout,
    )
    _annotate_heatmap_counts(axes["heat"], heat_values, layout)

    counts = [ghost_counts.get(taxon, 0) for taxon in taxa_order]
    x_max = max(max(counts, default=0), 1)
    gap_pt = _draw_ghost_bars(
        fig,
        axes["bar"],
        counts,
        [COUNT_FILL_COLOR] * len(taxa_order),
        [bool(ghost_has_sampled.get(taxon)) for taxon in taxa_order],
        x_limits=(0.0, x_max),
    )
    _annotate_ghost_counts(axes["bar"], counts, x_max, gap_pt, layout)
    axes["bar"].spines["bottom"].set_visible(True)
    axes["bar"].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes["bar"].set_xlabel(
        ghost_title, fontweight="bold", fontsize=layout.text_sizes["axis_label"]
    )
    _draw_ghost_legend_and_note(
        fig, axes, layout, "Heatmap cells with no introgression are left uncoloured."
    )

    fig.savefig(path, bbox_inches="tight", dpi=FIGURE_DPI)
    plt.close(fig)


def generate_introgression_maps(
    results,
    species_tree_path,
    output_dir,
    plot_taxa=None,
    outgroups=None,
    rename_map=None,
    overwrite=True,
):
    """Generate the bootstrap and counts introgression maps and their TSV matrices.

    Args:
        results: Iterable of TripletPipelineResult-like objects or dict rows.
        species_tree_path: Path to processed species tree used for ordering taxa.
        output_dir: Output directory for plots and tabular artifacts.
        plot_taxa: Optional iterable of taxa to retain in the plotted tree.
        outgroups: Optional iterable of taxon names to exclude from the plots
            (e.g. outgroup taxa used for rooting).  When ``None`` or empty no
            taxa are excluded.
        rename_map: Optional mapping of the species tree file's labels to the
            names ``results``, ``plot_taxa`` and ``outgroups`` use, for a run
            whose outputs carry display names while its processed trees keep
            the tree labels.
        overwrite: Whether to overwrite an existing output directory or write
            to an auto-suffixed sibling.
    """
    # The orchestrator hands this a dedicated ``consolidation/`` subfolder of a
    # freshly prepared run directory, so the reset only ever creates it.
    output_path = Path(prepare_output_directory(output_dir, overwrite=overwrite))

    outgroup_set = set(outgroups) if outgroups else set()

    non_ghost_weights, ghost_weights, taxa_seen = _collect_weights(results)
    species_tree = _load_species_tree(species_tree_path, rename_map)
    species_order = _species_tree_taxa_order(species_tree, plot_taxa)
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

    bootstrap_plot = output_path / BOOTSTRAP_PLOT_NAME
    counts_plot = output_path / COUNTS_PLOT_NAME

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
    # Cross-reference the sampled edges so each ghost target records whether it
    # also carries sampled introgression; this drives the bar outline.
    ghost_has_sampled = _sampled_introgression_presence(taxa_order, avg_non_ghost)
    _write_ghost_strength_tsv(
        ghost_tsv, taxa_order, avg_ghost, sampled_presence=ghost_has_sampled
    )
    _write_taxa_order_tsv(taxa_order_tsv, taxa_order)

    non_sister_counts = _collect_non_sister_counts(results)
    non_sister_matrix = _build_non_sister_matrix(taxa_order, non_sister_counts)
    _write_non_sister_matrix_tsv(non_sister_count_tsv, taxa_order, non_sister_matrix)

    _plot_bootstrap_map(
        bootstrap_plot,
        species_tree,
        taxa_order,
        avg_non_ghost_matrix,
        avg_ghost,
        supporting_ghost_values,
        ghost_has_sampled,
    )
    _plot_counts_map(
        counts_plot,
        species_tree,
        taxa_order,
        supporting_non_ghost_matrix,
        supporting_ghost_values,
        ghost_has_sampled,
    )

    return IntrogressionMapArtifacts(
        bootstrap_plot_path=str(bootstrap_plot),
        counts_plot_path=str(counts_plot),
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
