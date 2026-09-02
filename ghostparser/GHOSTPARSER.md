# GhostParser Package Guide

This is the top-level guide to the `ghostparser` package: what each module is
for, what they share, and where to read more. It deliberately stays at the
level of package structure — the per-module mechanics live in the module
guides, and every configuration key lives in [CONFIG.md](../CONFIG.md).

## Package layout

```
ghostparser/
  orchestrator/         the introgression engine (primary entry point)
  ml/                   optional multi-label classifiers over summary_statistics.tsv
  config.py             shared configuration trunk
  cli_config.py         shared CLI/config-file precedence resolver
  triplet_utils.py      shared triplet topology helpers
  __main__.py           usage banner for `python -m ghostparser`
```

| Module | Purpose | Guide |
| --- | --- | --- |
| `ghostparser.orchestrator` | Streaming triplet extraction + introgression inference. The main entry point. | [orchestrator/ORCHESTRATOR.md](orchestrator/ORCHESTRATOR.md) |
| `ghostparser.ml` | Trains multi-label classifiers on an orchestrator run's `summary_statistics.tsv`, and tunes their hyperparameters. Requires `pip install .[ml]`; Weights & Biases logging is a separate opt-in extra (`pip install .[wandb]`). | [ml/ML.md](ml/ML.md) |

`python -m ghostparser` prints a usage banner; the runnable entry points are
`python -m ghostparser.orchestrator` and `python -m ghostparser.ml`. The
introgression maps are produced by the orchestrator's own consolidation stage
(`orchestrator/consolidation.py`), described in
[orchestrator/ORCHESTRATOR.md](orchestrator/ORCHESTRATOR.md).

## What the modules share

`orchestrator` and `ml` are intentionally near-independent: each owns its own
defaults, choices, validation rules, and config loader, so their settings can
diverge. They share
only a thin trunk of helpers whose behaviour is identical for every caller.

### `ghostparser.config` — the configuration trunk

Holds exactly the pieces that behave the same everywhere:

- `ConfigError` — the shared exception type for invalid configuration.
- `_resolve_path` — `~`-expansion plus absolute/relative resolution.
- `_load_raw_config` — JSON/YAML loading with suffix and root-type validation.
- `_validate_required_path` — required-path validation and resolution.
- `DEFAULT_OVERWRITE` / `_validate_overwrite_flag` — overwrite resolution, where
  the canonical `overwrite` key wins over the CLI-style `no_overwrite`.
- `prepare_output_directory` (with `_next_available_suffixed_path`) — resets an
  existing output directory, or picks the smallest free `_<n>` sibling when
  `overwrite` is false.

Anything that differs between modules deliberately does **not** live here.
`orchestrator/config.py` and `ml/config.py` each define their own defaults and
their own validators — for example both have a `_validate_optional_float`, but
the ML one additionally requires a fraction strictly between 0 and 1.

### `ghostparser.cli_config` — CLI/config precedence

`resolve_cli_or_config_args(args, *, load_config, normalize_payload,
payload_arg_names)` implements the one rule both CLIs follow: if
`-c/--config-file` is given, the file supplies every setting and any other
supplied flags are reported as ignored; otherwise the CLI payload is normalized
directly. Each module passes in its own loader, normalizer, and payload key
list.

### `ghostparser.triplet_utils` — topology helpers

Pure functions for triplet topology handling (`find_sister_pair`,
`normalize_abc_from_sister_pair`, `classify_triplet_topology_string`,
`rank_topologies_by_frequency`), used by the orchestrator for both preprocessing and
inference.

## Consolidation Stage (`ghostparser.orchestrator.consolidation`)


### Role

`orchestrator/consolidation.py` turns the per-triplet results into a single combined visualization plus companion TSV artifacts representing introgression signal across the ingroup taxa. It is the orchestrator's final stage rather than an entry point of its own: the runner calls it with the in-memory results list, so no results TSV is re-read between inference and plotting. See [orchestrator/ORCHESTRATOR.md](orchestrator/ORCHESTRATOR.md) for how it is wired in.

### `generate_introgression_maps(results, species_tree_path, output_dir, plot_taxa=None, outgroups=None, overwrite=True)`

Generates the combined consolidation figure and tabular outputs.

**Arguments:**

- `results`: Iterable of `TripletPipelineResult` objects, as produced by the streaming engine.
- `species_tree_path`: Path to the processed species tree used for taxon ordering.
- `output_dir`: Directory to write all output files.
- `plot_taxa`: Optional list of taxa to retain in the plot; defaults to full ingroup.
- `outgroups`: Optional list of taxon names to exclude from all plots and TSVs (e.g. taxa used for rooting).
- `overwrite`: When `true`, existing output directories are cleared before writing; when `false`, a suffix such as `_1` is appended to avoid reusing an existing directory.

**Outputs:**

- `introgression_combined.png` — combined figure with directed inflow/outflow heatmap and ghost target-strength bar chart.
- `introgression_matrix_inflow_outflow.tsv` — target × source matrix of average bootstrap support values.
- `introgression_ghost_target_strength.tsv` — per-taxon average ghost bootstrap support, plus a `has_sampled_introgression` flag (`1` when the taxon is also the target of a sampled introgression edge).
- `introgression_taxa_order.tsv` — ordered taxa list matching the plot axes.

Returns an `IntrogressionMapArtifacts` dataclass with `plot_path`, TSV paths, `taxa_count`, `non_ghost_edge_count`, and `ghost_target_count`.

### Bootstrap averaging denominators

Average bootstrap support values are computed using population-level co-occurrence denominators:

- **Sampled introgression** for a directed pair (source → target):

  `avg = sum(bootstrap_value for classified triplets) / count(all triplets containing both source and target)`

  The denominator counts every triplet in which both taxa co-appear, regardless of whether that triplet was classified as introgression. For `n` ingroup taxa, a directed pair appears in exactly `n − 2` triplets.

- **Ghost introgression** for a target taxon:

  `avg = sum(bootstrap_value for ghost-classified triplets) / count(all triplets containing that taxon)`

  The denominator counts every triplet in which the target taxon appears in any position, regardless of classification. For `n` ingroup taxa, a single taxon appears in `(n−1)C2` triplets.

This normalizes signal strength by the total number of opportunities at which the event could have been detected, giving a population-level confidence estimate that accounts for the full co-occurrence space.

### Plot layout

The combined figure uses a three-row layout above the data panels:

1. **Species tree strip** (top row) — topology-only tree with leaf labels suppressed.
2. **Source taxon label strip** (middle row) — a dedicated thin row containing the source-taxon names, rotated 90°, aligned to heatmap column centres. Row height is computed from the rendered pixel-width of the longest label so labels are never clipped. The labels are always drawn; the canvas is sized from the taxon count, so a larger tree produces a larger figure.
3. **Data panels** (bottom row, left to right):
    - **Inflow/outflow heatmap** — rows are target taxa, columns are source taxa, coloured by average bootstrap support on the `CONSOLIDATION_COLORMAP` (`cividis`) scale.
    - **Target label panel** — centred target taxon names aligned pixel-exactly to heatmap rows.
    - **Ghost bar chart** — horizontal bars per target taxon. Bar *length* is the average ghost bootstrap support. Bar *colour* is constant per bar, drawn from the two extremes of the same colormap, and encodes only whether the taxon also has sampled introgression: the high end, yellow (`GHOST_ONLY_BAR_COLOR`) when the taxon's only signal is ghost introgression, the low end, dark blue (`GHOST_WITH_SAMPLED_BAR_COLOR`) when it is also the target of a sampled introgression edge. Every bar carries a hairline `GHOST_BAR_EDGE_COLOR` outline so its extent stays legible against the panel. The flag is computed by `_sampled_introgression_presence` and recorded in the ghost TSV's `has_sampled_introgression` column.
    - **Colorbar** — applies to the heatmap only. A two-entry legend sits directly above the bar panel explaining the ghost bar colours, and a note above that states that uncoloured heatmap cells carry no introgression and are not on the colour scale.

All three rows share `hspace=0` so they appear flush. Figure and panel widths scale dynamically with taxon count and rendered label widths.
The shared x/y labels and colorbar text scale with taxon count and are capped to stay readable on large figures, while the species-name labels keep their separate sizing.

### How it is invoked

Consolidation runs automatically as the orchestrator's last stage, writing into
a `consolidation/` subfolder of the run's output directory. Disable it with
`--no-consolidation` / `consolidation: false`. It has no CLI of its own: its
input is the in-memory results list, not a file, so there is nothing to point a
command line at.
