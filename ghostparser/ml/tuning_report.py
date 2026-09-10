"""Local reporting for hyperparameter tuning runs.

Builds the plaintext report, the per-parameter marginal table, and the
navigation plots that make a tuning run interpretable without Weights & Biases.
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Two-slot categorical palette for the search-report figure: candidate scores and
# the running-best trace are distinct series, not a magnitude ramp. The pair
# clears the lightness-band, chroma, CVD-separation and contrast checks against
# the light chart surface below.
SERIES_CANDIDATE = "#2a78d6"
SERIES_RUNNING_BEST = "#eb6834"
CHART_SURFACE = "#fcfcfb"
PRIMARY_INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED_INK = "#898781"
GRIDLINE = "#e1e0d9"
AXIS_RULE = "#c3c2b7"

PARAM_PLOT_COLUMNS = 3
DEFAULT_PLOT_TOP_K = 10
_NEGLIGIBLE = 1e-12


def _format_value(value: object) -> str:
    """Render a search-space value the way the config file spells it."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _value_group_key(value: object) -> str:
    """Group key that also works for unhashable search-space values."""
    return repr(value)


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _clean(value: float) -> float:
    """Flatten floating-point noise so identical scores report an exact zero spread."""
    return 0.0 if abs(value) < _NEGLIGIBLE else float(value)


def _render_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    """Render an aligned fixed-width text table as a list of lines."""
    widths = [len(header) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    def render_row(cells: list[str]) -> str:
        return "  ".join(cell.ljust(widths[index]) for index, cell in enumerate(cells))

    lines = [render_row(headers).rstrip(), "  ".join("-" * width for width in widths)]
    lines.extend(render_row(row).rstrip() for row in rows)
    return lines


def compute_parameter_marginals(
    ranked_candidates: list[dict[str, object]],
    parameter_keys: list[str],
    objective_direction: str,
) -> list[dict[str, object]]:
    """Summarize candidate scores for each searched parameter value.

    Args:
        ranked_candidates: Candidate rows already sorted best-first and carrying
            a ``rank`` and a ``cv_score``.
        parameter_keys: Search-space parameter names to summarize.
        objective_direction: ``max`` or ``min``.

    Returns:
        One row per (parameter, value) pair with the count of candidates that
        used the value, the best/mean/worst score it reached, the best rank it
        achieved, and its rank among the values of the same parameter.
    """
    minimize = objective_direction == "min"
    marginals: list[dict[str, object]] = []

    for parameter in parameter_keys:
        groups: dict[str, dict[str, object]] = {}
        for row in ranked_candidates:
            if parameter not in row:
                continue
            group = groups.setdefault(
                _value_group_key(row[parameter]),
                {"value": row[parameter], "scores": [], "ranks": []},
            )
            group["scores"].append(float(row["cv_score"]))
            group["ranks"].append(int(row["rank"]))

        parameter_rows: list[dict[str, object]] = []
        for group in groups.values():
            scores = np.asarray(group["scores"], dtype=float)
            parameter_rows.append(
                {
                    "parameter": parameter,
                    "value": _format_value(group["value"]),
                    "raw_value": group["value"],
                    "candidate_count": int(scores.size),
                    "best_score": float(scores.min() if minimize else scores.max()),
                    "mean_score": float(scores.mean()),
                    "std_score": _clean(float(scores.std(ddof=0))),
                    "worst_score": float(scores.max() if minimize else scores.min()),
                    "best_rank": int(min(group["ranks"])),
                }
            )

        parameter_rows.sort(key=lambda row: row["best_rank"])
        for value_rank, row in enumerate(parameter_rows, start=1):
            row["value_rank"] = value_rank

        raw_values = [row["raw_value"] for row in parameter_rows]
        bounds = None
        if len(raw_values) > 1 and all(_is_number(value) for value in raw_values):
            bounds = (min(raw_values), max(raw_values))
        for row in parameter_rows:
            boundary = ""
            if bounds is not None and row["value_rank"] == 1:
                if row["raw_value"] == bounds[0]:
                    boundary = "lower_bound"
                elif row["raw_value"] == bounds[1]:
                    boundary = "upper_bound"
            row["at_search_bound"] = boundary

        marginals.extend(parameter_rows)

    return marginals


def marginal_tsv_rows(marginals: list[dict[str, object]]) -> list[dict[str, object]]:
    """Drop the un-serializable raw value so the marginals can be written to TSV."""
    return [
        {key: value for key, value in row.items() if key != "raw_value"}
        for row in marginals
    ]


def parameter_influence(
    marginals: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Rank parameters by how far the objective moved across the values tried.

    A large spread means the parameter mattered; a spread of zero on both the
    best and the mean score means every value performed alike and the dimension
    can be dropped from the next search.
    """
    grouped: dict[str, dict[str, list[float]]] = {}
    for row in marginals:
        entry = grouped.setdefault(row["parameter"], {"best": [], "mean": []})
        entry["best"].append(float(row["best_score"]))
        entry["mean"].append(float(row["mean_score"]))

    influence = [
        {
            "parameter": parameter,
            "best_score_spread": _clean(max(scores["best"]) - min(scores["best"])),
            "mean_score_spread": _clean(max(scores["mean"]) - min(scores["mean"])),
        }
        for parameter, scores in grouped.items()
    ]
    influence.sort(
        key=lambda row: (
            -row["best_score_spread"],
            -row["mean_score_spread"],
            row["parameter"],
        )
    )
    return influence


def format_candidate_table(
    ranked_candidates: list[dict[str, object]],
    parameter_keys: list[str],
    objective_metric: str,
    limit: int,
) -> list[str]:
    """Aligned table of the best candidates with their parameters and scores."""
    if not ranked_candidates:
        return ["  (no candidates were evaluated)"]

    headers = ["rank", "score", f"{objective_metric}_std", "seconds"]
    headers.extend(parameter_keys)
    rows: list[list[str]] = []
    for row in ranked_candidates[:limit]:
        std_value = row.get(f"{objective_metric}_std")
        elapsed = row.get("elapsed_seconds")
        cells = [
            str(row.get("rank", "")),
            f"{float(row['cv_score']):.6f}",
            "n/a" if std_value is None else f"{float(std_value):.6f}",
            "n/a" if elapsed is None else f"{float(elapsed):.2f}",
        ]
        cells.extend(_format_value(row.get(key)) for key in parameter_keys)
        rows.append(cells)
    return [f"  {line}" for line in _render_table(headers, rows)]


def format_marginal_table(marginals: list[dict[str, object]]) -> list[str]:
    """Aligned table of the per-parameter value summaries."""
    if not marginals:
        return ["  (no search-space parameters to summarize)"]

    headers = [
        "parameter",
        "value",
        "candidates",
        "best_score",
        "mean_score",
        "std_score",
        "worst_score",
        "best_rank",
        "value_rank",
    ]
    rows = [
        [
            row["parameter"],
            row["value"],
            str(row["candidate_count"]),
            f"{row['best_score']:.6f}",
            f"{row['mean_score']:.6f}",
            f"{row['std_score']:.6f}",
            f"{row['worst_score']:.6f}",
            str(row["best_rank"]),
            str(row["value_rank"]),
        ]
        for row in marginals
    ]
    return [f"  {line}" for line in _render_table(headers, rows)]


def format_search_space_guidance(
    marginals: list[dict[str, object]],
    ranked_candidates: list[dict[str, object]],
) -> list[str]:
    """Plain-language notes on where to take the search space next."""
    if not ranked_candidates:
        return ["  (no candidates were evaluated)"]

    lines: list[str] = []
    scores = [float(row["cv_score"]) for row in ranked_candidates]
    spread = max(scores) - min(scores)
    lines.append(
        f"  Score range across {len(scores)} candidates: "
        f"{min(scores):.6f} to {max(scores):.6f} (spread {spread:.6f})"
    )
    if spread == 0.0:
        lines.append(
            "  Every candidate scored identically, so nothing in this search "
            "space moved the objective. Widen the ranges or change the objective."
        )

    influence = parameter_influence(marginals)
    if influence:
        lines.append("  Parameter influence (score spread across the values tried):")
        for row in influence:
            inert = (
                row["best_score_spread"] == 0.0 and row["mean_score_spread"] == 0.0
            )
            note = " - no effect on the objective" if inert else ""
            lines.append(
                f"    {row['parameter']}: best {row['best_score_spread']:.6f}, "
                f"mean {row['mean_score_spread']:.6f}{note}"
            )

    best_by_parameter = [row for row in marginals if row["value_rank"] == 1]
    if best_by_parameter:
        lines.append("  Best value per parameter:")
        for row in best_by_parameter:
            lines.append(
                f"    {row['parameter']}={row['value']} "
                f"(best_score {row['best_score']:.6f})"
            )

    bound_notes = [row for row in best_by_parameter if row["at_search_bound"]]
    if bound_notes:
        lines.append("  Values sitting at the edge of the searched range:")
        for row in bound_notes:
            edge = "smallest" if row["at_search_bound"] == "lower_bound" else "largest"
            lines.append(
                f"    {row['parameter']}={row['value']} is the {edge} value tried; "
                "extend the range in that direction to check for a better optimum."
            )
    return lines


def build_results_text(
    *,
    config,
    ranked_candidates: list[dict[str, object]],
    marginals: list[dict[str, object]],
    parameter_keys: list[str],
    best_candidate_record: dict[str, object],
    effective_n_neighbors: int | None,
    test_metrics: dict,
    dataset_summary: dict,
    cv_warnings: list[str],
    split_notes: list[str],
    timings: dict[str, float],
    total_seconds: float,
    cv_folds: int,
    wandb_enabled: bool,
    artifact_paths: dict[str, str],
) -> str:
    """Assemble the full plaintext tuning report."""
    objective_direction = config.objective_direction
    lines = [
        f"Ghostparser ML hyperparameter tuning ({config.model_name})",
        f"Search method: {config.search_method}",
        f"Objective metric: {config.objective_metric} ({objective_direction})",
        f"Candidates evaluated: {len(ranked_candidates)}",
        f"CV folds: {cv_folds}",
        f"Weights & Biases logging: {'enabled' if wandb_enabled else 'disabled'}",
        "",
        "Search space:",
    ]
    if config.search_space:
        for key in parameter_keys:
            values = config.search_space[key]
            if not isinstance(values, (list, tuple)):
                values = [values]
            rendered = ", ".join(_format_value(value) for value in values)
            lines.append(f"  {key}: [{rendered}]")
    else:
        lines.append("  (empty)")

    lines.extend(
        [
            "",
            "Best candidate:",
            f"  Candidate index: {best_candidate_record['candidate_index']}",
        ]
    )
    if effective_n_neighbors is not None:
        lines.append(f"  Effective n_neighbors: {effective_n_neighbors}")
    for key, value in best_candidate_record["candidate_params"].items():
        lines.append(f"  {key}: {value}")
    lines.append(f"  CV score: {float(best_candidate_record['cv_score']):.6f}")

    lines.extend(
        [
            "",
            "Test metrics:",
            f"  Hamming loss: {test_metrics['hamming_loss']:.6f}",
            f"  Bitwise accuracy: {test_metrics['bitwise_accuracy']:.6f}",
            f"  Exact-match accuracy: {test_metrics['exact_match_accuracy']:.6f}",
            f"  Micro F1: {test_metrics['micro_f1']:.6f}",
            f"  Macro F1: {test_metrics['macro_f1']:.6f}",
            f"  Weighted F1: {test_metrics['weighted_f1']:.6f}",
            "",
            "Dataset summary:",
            "  Label map:",
        ]
    )
    for key, value in dataset_summary["label_map"].items():
        lines.append(f"    {key}: {value}")
    lines.append("  Split:")
    for key, value in dataset_summary["split"].items():
        lines.append(f"    {key}: {value}")

    lines.extend(["", f"Top {config.top_k} candidates:"])
    lines.extend(
        format_candidate_table(
            ranked_candidates,
            parameter_keys,
            config.objective_metric,
            config.top_k,
        )
    )

    lines.extend(["", "Per-parameter value summary:"])
    lines.extend(format_marginal_table(marginals))

    lines.extend(["", "Search-space guidance:"])
    lines.extend(format_search_space_guidance(marginals, ranked_candidates))

    if cv_warnings:
        lines.extend(["", "CV notes:"] + [f"  {note}" for note in cv_warnings])
    if split_notes:
        lines.extend(["", "Split notes:"] + [f"  {note}" for note in split_notes])

    lines.extend(["", "Timings (seconds):"])
    for key, value in timings.items():
        lines.append(f"  {key}: {value:.6f}")
    lines.append(f"  total: {total_seconds:.6f}")

    lines.extend(["", "Artifacts:"])
    for key, value in artifact_paths.items():
        lines.append(f"  {key}: {value}")
    if wandb_enabled:
        lines.append(
            "  The ranked candidates, parameter marginals, predictions and full "
            "results payload are logged to the Weights & Biases run instead of "
            "this directory."
        )

    return "\n".join(lines) + "\n"


def _style_axis(ax, *, grid_axis: str = "y") -> None:
    """Recessive chrome: hairline solid grid, muted ink, no top/right spines."""
    ax.set_facecolor(CHART_SURFACE)
    ax.grid(True, axis=grid_axis, color=GRIDLINE, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS_RULE)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED_INK, labelsize=9, length=3, width=0.8)
    ax.title.set_color(PRIMARY_INK)
    ax.xaxis.label.set_color(SECONDARY_INK)
    ax.yaxis.label.set_color(SECONDARY_INK)


def _padded_limits(values: np.ndarray, pad_fraction: float = 0.08) -> tuple[float, float]:
    """Data-range limits with breathing room, widened when every value is equal."""
    low = float(values.min())
    high = float(values.max())
    span = high - low
    if span < _NEGLIGIBLE:
        pad = max(abs(high), 1.0) * 0.05
        return low - pad, high + pad
    pad = span * pad_fraction
    return low - pad, high + pad


def _plot_search_progress(ax, ranked_candidates, objective_metric, objective_direction):
    ordered = sorted(ranked_candidates, key=lambda row: int(row["candidate_index"]))
    indices = np.asarray([int(row["candidate_index"]) for row in ordered], dtype=int)
    scores = np.asarray([float(row["cv_score"]) for row in ordered], dtype=float)
    running_best = (
        np.minimum.accumulate(scores)
        if objective_direction == "min"
        else np.maximum.accumulate(scores)
    )

    ax.step(
        indices,
        running_best,
        where="post",
        color=SERIES_RUNNING_BEST,
        linewidth=2.0,
        label=f"running best ({objective_direction})",
        zorder=2,
    )
    ax.scatter(
        indices,
        scores,
        color=SERIES_CANDIDATE,
        s=70,
        edgecolor=CHART_SURFACE,
        linewidth=1.0,
        label="candidate score",
        zorder=3,
    )

    best_position = int(
        np.argmin(scores) if objective_direction == "min" else np.argmax(scores)
    )
    ax.annotate(
        f"best {scores[best_position]:.4f} at candidate {indices[best_position]}",
        xy=(indices[best_position], scores[best_position]),
        xytext=(6, 8),
        textcoords="offset points",
        fontsize=9,
        color=SECONDARY_INK,
    )

    ax.set_title("Search progress in evaluation order", fontsize=12)
    ax.set_xlabel("Candidate index (order evaluated)")
    ax.set_ylabel(objective_metric)
    ax.set_ylim(*_padded_limits(scores, 0.12))
    _style_axis(ax)
    legend = ax.legend(loc="lower right", fontsize=9, frameon=False)
    for entry in legend.get_texts():
        entry.set_color(SECONDARY_INK)


def _plot_top_candidates(
    ax, ranked_candidates, parameter_keys, objective_metric, top_k
):
    """Ranked dot plot.

    Objective scores cluster in a narrow band well away from zero, so a bar from
    a zero baseline would render every candidate as the same length; dots carry
    the comparison honestly on a zoomed axis.
    """
    top_rows = list(reversed(ranked_candidates[:top_k]))
    labels = [
        ", ".join(f"{key}={_format_value(row.get(key))}" for key in parameter_keys)
        or f"candidate {row['candidate_index']}"
        for row in top_rows
    ]
    scores = np.asarray([float(row["cv_score"]) for row in top_rows], dtype=float)
    positions = np.arange(len(top_rows))

    ax.scatter(
        scores,
        positions,
        color=SERIES_CANDIDATE,
        s=90,
        edgecolor=CHART_SURFACE,
        linewidth=1.0,
        zorder=3,
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.7, len(top_rows) - 0.3)
    ax.set_xlabel(objective_metric)
    ax.set_title(f"Top {len(top_rows)} candidates by {objective_metric}", fontsize=12)

    # Extra room on the right so the direct labels sit inside the axes.
    low, high = _padded_limits(scores, 0.08)
    high += (high - low) * 0.22
    ax.set_xlim(low, high)
    offset = (high - low) * 0.015
    for position, score in zip(positions, scores, strict=True):
        ax.text(
            score + offset,
            position,
            f"{score:.4f}",
            va="center",
            ha="left",
            fontsize=8,
            color=SECONDARY_INK,
        )
    _style_axis(ax, grid_axis="x")


def _plot_parameter_marginal(ax, parameter, ranked_candidates, objective_metric):
    groups: dict[str, dict[str, object]] = {}
    for row in ranked_candidates:
        if parameter not in row:
            continue
        group = groups.setdefault(
            _value_group_key(row[parameter]),
            {"value": row[parameter], "scores": []},
        )
        group["scores"].append(float(row["cv_score"]))

    ordered = list(groups.values())
    if all(_is_number(group["value"]) for group in ordered):
        ordered.sort(key=lambda group: group["value"])
    else:
        ordered.sort(key=lambda group: _format_value(group["value"]))

    labels = [_format_value(group["value"]) for group in ordered]
    data = [group["scores"] for group in ordered]
    positions = np.arange(1, len(data) + 1)

    # The box is recessive context; the candidate scores themselves are the series.
    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.5,
        showfliers=False,
        medianprops={"color": SECONDARY_INK, "linewidth": 1.4},
        boxprops={"color": AXIS_RULE, "linewidth": 0.9},
        whiskerprops={"color": AXIS_RULE, "linewidth": 0.9},
        capprops={"color": AXIS_RULE, "linewidth": 0.9},
    )
    del box

    rng = np.random.default_rng(0)
    for position, scores in zip(positions, data, strict=True):
        jitter = rng.uniform(-0.11, 0.11, size=len(scores))
        ax.scatter(
            position + jitter,
            scores,
            s=46,
            color=SERIES_CANDIDATE,
            edgecolor=CHART_SURFACE,
            linewidth=1.0,
            zorder=3,
        )

    all_scores = np.asarray([score for scores in data for score in scores], dtype=float)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9, rotation=30, ha="right")
    ax.set_xlim(0.4, len(data) + 0.6)
    ax.set_ylim(*_padded_limits(all_scores, 0.12))
    ax.set_title(parameter, fontsize=12)
    ax.set_ylabel(objective_metric)
    _style_axis(ax)


def save_tuning_plots(
    ranked_candidates: list[dict[str, object]],
    parameter_keys: list[str],
    objective_metric: str,
    objective_direction: str,
    output_path: Path,
    top_k: int = DEFAULT_PLOT_TOP_K,
) -> str | None:
    """Write the combined search-progress, top-candidate, and marginal figure.

    Args:
        ranked_candidates: Candidate rows sorted best-first.
        parameter_keys: Search-space parameter names; each one that actually
            varied gets its own marginal panel.
        objective_metric: Metric name used for the axis labels.
        objective_direction: ``max`` or ``min``, drives the running-best trace.
        output_path: Destination image path.
        top_k: Number of candidates in the top-candidate panel.

    Returns:
        The written path as a string, or ``None`` when there is nothing to plot.
    """
    if not ranked_candidates:
        return None

    panel_keys = [
        key
        for key in parameter_keys
        if len({_value_group_key(row.get(key)) for row in ranked_candidates}) > 1
    ]
    parameter_rows = math.ceil(len(panel_keys) / PARAM_PLOT_COLUMNS)
    shown_candidates = min(top_k, len(ranked_candidates))

    progress_height = 3.6
    top_height = max(2.4, 0.34 * shown_candidates + 1.3)
    parameter_height = 3.4
    figure_height = progress_height + top_height + parameter_height * parameter_rows
    height_ratios = [progress_height, top_height] + [parameter_height] * parameter_rows

    fig = plt.figure(
        figsize=(5.6 * PARAM_PLOT_COLUMNS, figure_height),
        constrained_layout=True,
        facecolor=CHART_SURFACE,
    )
    grid = fig.add_gridspec(
        2 + parameter_rows, PARAM_PLOT_COLUMNS, height_ratios=height_ratios
    )

    _plot_search_progress(
        fig.add_subplot(grid[0, :]),
        ranked_candidates,
        objective_metric,
        objective_direction,
    )
    _plot_top_candidates(
        fig.add_subplot(grid[1, :]),
        ranked_candidates,
        parameter_keys,
        objective_metric,
        top_k,
    )
    for panel_index, parameter in enumerate(panel_keys):
        row = 2 + panel_index // PARAM_PLOT_COLUMNS
        column = panel_index % PARAM_PLOT_COLUMNS
        _plot_parameter_marginal(
            fig.add_subplot(grid[row, column]),
            parameter,
            ranked_candidates,
            objective_metric,
        )

    fig.suptitle(
        f"Hyperparameter search over {objective_metric} ({objective_direction})",
        fontsize=15,
        color=PRIMARY_INK,
    )
    fig.savefig(
        output_path, dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor()
    )
    plt.close(fig)
    return str(output_path)


__all__ = [
    "build_results_text",
    "compute_parameter_marginals",
    "format_candidate_table",
    "format_marginal_table",
    "format_search_space_guidance",
    "marginal_tsv_rows",
    "parameter_influence",
    "save_tuning_plots",
]
