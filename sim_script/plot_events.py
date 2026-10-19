"""Draw one replicate's scenario: the full species tree with every realized event on it."""

from __future__ import annotations

import math
from pathlib import Path

from demography import Event, branch_sizes
from sim_config import SimConfig
from species_tree import Node

# Chart surface, ink and the first three categorical slots of the reference
# palette (validated all-pairs for colour-vision deficiency). Size changes and
# extinctions stay in ink, so three hues cover every event kind.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
EVENT_COLOURS = {"pulse": "#2a78d6", "gene_flow": "#eb6834", "bottleneck": "#1baf7a"}

_MIN_WIDTH, _MAX_EXTRA_WIDTH = 1.2, 3.8


def _x_positions(root: Node) -> dict[int, float]:
    x: dict[int, float] = {}
    for i, tip in enumerate(root.tips()):
        x[tip.index] = float(i)
    for node in root.postorder():
        if not node.is_tip:
            x[node.index] = sum(x[c.index] for c in node.children) / len(node.children)
    return x


def _branch_style(node: Node, cfg: SimConfig) -> tuple[str, str]:
    """Colour and dash style: solid ink when a sampled ingroup tip descends, muted
    for the outgroup, dashed grey for a lineage nothing sampled descends from."""
    roles = {cfg.species.role(name) for name in node.tip_names()}
    if "sampled" in roles:
        return INK, "-"
    if roles == {"outgroup"}:
        return INK_MUTED, "-"
    return INK_MUTED, (0, (4, 3))


def _size_segments(node: Node, start: float, events: list[Event]) -> list[tuple[float, float, float]]:
    """``(bottom, top, Ne)`` pieces of a branch, stepping at each size change on it."""
    top = node.parent_age
    changes = sorted((e.time, e.size) for e in events
                     if e.kind == "size_change" and e.branch is node and e.size is not None)
    pieces, low, size = [], node.age, start
    for time, new in changes:
        pieces.append((low, time, size))
        low, size = time, new
    pieces.append((low, top, size))
    return pieces


def plot_scenario(cfg: SimConfig, root: Node, events: list[Event], path: Path, title: str, caption: str) -> None:
    """Write the scenario figure for one replicate.

    Args:
        cfg: The validated config.
        root: The replicate's full species tree, nodes indexed, ages in generations.
        events: The resolved events, as written to ``events.tsv``.
        path: Output file; its suffix picks the format.
        title: Figure title.
        caption: One line of run settings printed under the plot.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch, Patch, Rectangle

    per_unit = cfg.generations_per_unit
    x = _x_positions(root)
    tips = root.tips()
    stem = root.age * 0.06
    top = (root.age + stem) / per_unit

    starts = branch_sizes(cfg, root)
    segments = {n.index: _size_segments(n, starts[n.index], events) for n in root.preorder()}
    biggest = max(s for pieces in segments.values() for *_, s in pieces)

    def width(size: float) -> float:
        return _MIN_WIDTH + _MAX_EXTRA_WIDTH * math.sqrt(size / biggest)

    fig, ax = plt.subplots(figsize=(max(6.5, 0.75 * len(tips) + 3.0), 6.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    for node in root.preorder():
        colour, dash = _branch_style(node, cfg)
        for low, high, size in segments[node.index]:
            high = root.age + stem if math.isinf(high) else high
            ax.plot([x[node.index]] * 2, [low / per_unit, high / per_unit], color=colour,
                    linestyle=dash, linewidth=width(size), solid_capstyle="butt", zorder=2)
        if not node.is_tip:
            xs = [x[c.index] for c in node.children]
            ax.plot([min(xs), max(xs)], [node.age / per_unit] * 2, color=colour, linewidth=1.0, zorder=2)

    span = top
    for tip in tips:
        role = cfg.species.role(tip.name)
        hidden = role in ("ghost", "extinct")
        ax.text(x[tip.index], -0.03 * span, tip.name, ha="center", va="top", fontsize=9,
                color=INK_MUTED if hidden or role == "outgroup" else INK,
                style="italic" if hidden else "normal", rotation=90 if len(tips) > 14 else 0)
        if role == "extinct" or tip.age > 0:
            ax.text(x[tip.index], tip.age / per_unit, "†", ha="center", va="top", fontsize=12,
                    color=INK, zorder=4)

    kinds_seen: set[str] = set()
    pulses_drawn = 0
    for event in events:
        t = event.time / per_unit
        if event.kind == "pulse":
            kinds_seen.add("pulse_" + event.origin)
            colour = EVENT_COLOURS["pulse"]
            x0, x1 = x[event.donor.index], x[event.recipient.index]
            ax.add_patch(FancyArrowPatch(
                (x0, t), (x1, t), arrowstyle="-|>", mutation_scale=14, linewidth=1.6,
                edgecolor=colour, facecolor=SURFACE if event.origin == "random" else colour, zorder=5,
            ))
            # Alternate the label's place along the arrow so close pulses stay legible.
            along = 0.3 if pulses_drawn % 2 == 0 else 0.7
            pulses_drawn += 1
            ax.text(x0 + along * (x1 - x0), t + 0.012 * span, f"{event.proportion:.2f}", ha="center",
                    va="bottom", fontsize=8, color=INK_SECONDARY, zorder=6)
        elif event.kind == "gene_flow":
            kinds_seen.add("gene_flow")
            colour = EVENT_COLOURS["gene_flow"]
            x0, x1 = x[event.donor.index], x[event.recipient.index]
            low, high = t, event.end / per_unit
            ax.add_patch(Rectangle((min(x0, x1), low), abs(x1 - x0), high - low, facecolor=colour,
                                   alpha=0.12, edgecolor="none", zorder=1))
            # A symmetric pair gets one arrow each way, set apart so both read.
            mid = low + (high - low) * (0.45 if x0 < x1 else 0.55)
            ax.add_patch(FancyArrowPatch((x0, mid), (x1, mid), arrowstyle="-|>", mutation_scale=12,
                                         linewidth=1.2, linestyle=(0, (3, 2)), color=colour, zorder=5))
            ax.text((x0 + x1) / 2, mid + 0.012 * span, f"m={event.rate:.2g}", ha="center", va="bottom",
                    fontsize=8, color=INK_SECONDARY, zorder=6)
        elif event.kind == "bottleneck":
            kinds_seen.add("bottleneck")
            ax.plot(x[event.branch.index], t, marker="D", markersize=8, color=EVENT_COLOURS["bottleneck"],
                    markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=6)
        elif event.kind == "size_change":
            kinds_seen.add("size_change")
            ax.plot(x[event.branch.index], t, marker="_", markersize=14, markeredgewidth=2, color=INK, zorder=6)
            note = f"Ne {event.size:g}" if event.size is not None else f"r {event.growth_rate:g}"
            ax.text(x[event.branch.index] + 0.12, t, note, ha="left", va="center", fontsize=7,
                    color=INK_SECONDARY, zorder=6)
        elif event.kind == "extinction":
            kinds_seen.add("extinction")

    handles = [Line2D([], [], color=INK, linewidth=2, label="sampled lineage")]
    if any(cfg.species.role(n) in ("ghost", "extinct") for n in root.tip_names()):
        handles.append(Line2D([], [], color=INK_MUTED, linewidth=2, linestyle=(0, (4, 3)),
                              label="ghost / extinct lineage"))
    handles.append(Line2D([], [], color=INK_MUTED, linewidth=2, label="outgroup"))
    pulse = EVENT_COLOURS["pulse"]
    if "pulse_config" in kinds_seen:
        handles.append(Line2D([], [], color=pulse, marker=">", markersize=8, linewidth=1.6,
                              label="pulse (configured), proportion"))
    if "pulse_random" in kinds_seen:
        handles.append(Line2D([], [], color=pulse, marker=">", markersize=8, linewidth=1.6,
                              markerfacecolor=SURFACE, label="pulse (random), proportion"))
    if "gene_flow" in kinds_seen:
        handles.append(Patch(facecolor=EVENT_COLOURS["gene_flow"], alpha=0.35, label="gene flow, rate per gen."))
    if "bottleneck" in kinds_seen:
        handles.append(Line2D([], [], color=EVENT_COLOURS["bottleneck"], marker="D", markersize=7,
                              linestyle="none", label="bottleneck"))
    if "size_change" in kinds_seen:
        handles.append(Line2D([], [], color=INK, marker="_", markersize=12, markeredgewidth=2,
                              linestyle="none", label="size change"))
    if "extinction" in kinds_seen:
        handles.append(Line2D([], [], color=INK, marker="$†$", markersize=9, linestyle="none",
                              label="extinction"))
    legend = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
                       fontsize=8, labelcolor=INK_SECONDARY, title="branch width ∝ Ne", title_fontsize=8)
    legend.get_title().set_color(INK_MUTED)

    ax.set_xlim(-0.6, len(tips) - 0.4)
    ax.set_ylim(-0.02 * span, top * 1.02)
    ax.set_xticks([])
    ax.set_ylabel(f"time before present ({cfg['time_units']})", color=INK_SECONDARY, fontsize=9)
    ax.tick_params(axis="y", colors=INK_MUTED, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    fig.text(0.01, 0.01, caption, fontsize=7.5, color=INK_MUTED, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
