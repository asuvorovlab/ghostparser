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

`python -m ghostparser` prints a usage banner, as does `python -m ghostparser.ml`;
the runnable entry points are `python -m ghostparser.orchestrator`,
`python -m ghostparser.ml.random_forest`, `python -m ghostparser.ml.multi_knn`
and `python -m ghostparser.ml.hyper_tune`. The
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
- Path resolution — `~`-expansion plus absolute/relative resolution.
- Raw config loading — JSON/YAML loading with suffix and root-type validation.
- Required-path validation and resolution.
- `DEFAULT_OVERWRITE` and overwrite resolution, where the canonical `overwrite`
  key wins over the CLI-style `no_overwrite`.
- `prepare_output_directory` — resets an existing output directory, or picks
  the smallest free `_<n>` sibling when `overwrite` is false.

Anything that differs between modules deliberately does **not** live here.
`orchestrator/config.py` and `ml/config.py` each define their own defaults and
their own validators — for example both validate an optional float, but the
ML one additionally requires a fraction strictly between 0 and 1.

Within the ML subpackage, `ml/config.py` is the single home for that module's
validators and for the CLI plumbing its two trainers share
(`build_trainer_argument_parser`, `resolve_trainer_runtime_args`);
`ml/hyper_tune.py` imports them rather than restating them.

### `ghostparser.cli_config` — CLI/config precedence

`resolve_cli_or_config_args(args, *, load_config, normalize_payload,
payload_arg_names)` implements the one rule both CLIs follow: if
`-c/--config-file` is given, the file supplies every setting and any other
supplied flags are reported as ignored; otherwise the CLI payload is normalized
directly. Each module passes in its own loader, normalizer, and payload key
list.

### `ghostparser.triplet_utils` — topology helpers

Pure functions for triplet topology handling (`triplet_taxa_labels`,
`find_sister_pair`, `normalize_abc_from_sister_pair`,
`topology_from_sister_pair`, `classify_triplet_topology_string`), used by the
orchestrator for both preprocessing and inference.

## Consolidation Stage (`ghostparser.orchestrator.consolidation`)


### Role

`orchestrator/consolidation.py` turns the per-triplet results into a single combined visualization plus companion TSV artifacts representing introgression signal across the ingroup taxa. It is the orchestrator's final stage rather than an entry point of its own: the runner calls it with the in-memory results list, so no results TSV is re-read between inference and plotting. See [orchestrator/ORCHESTRATOR.md](orchestrator/ORCHESTRATOR.md) for how it is wired in.

### `generate_introgression_maps(results, species_tree_path, output_dir, plot_taxa=None, outgroups=None, rename_map=None, overwrite=True)`

Generates the combined consolidation figure and tabular outputs.

**Arguments:**

- `results`: Iterable of `TripletPipelineResult` objects, as produced by the streaming engine.
- `species_tree_path`: Path to the processed species tree used for taxon ordering.
- `output_dir`: Directory to write all output files.
- `plot_taxa`: Optional list of taxa to retain in the plot; defaults to full ingroup.
- `outgroups`: Optional list of taxon names to exclude from all plots and TSVs (e.g. taxa used for rooting).
- `rename_map`: Optional mapping from the species tree file's labels to the names the results carry, so the tree can be lined up with results written under display names.
- `overwrite`: When `true`, existing output directories are cleared before writing; when `false`, a suffix such as `_1` is appended to avoid reusing an existing directory.

**Outputs:** the figure in `output_dir`, the TSVs in its `consolidation_data/` subfolder.

- `introgression_combined.png` — combined figure with directed inflow/outflow heatmap and ghost target-strength bar chart.
- `introgression_matrix_inflow_outflow.tsv` — target × source matrix of average bootstrap support, with the summed support in `introgression_matrix_inflow_outflow_raw_sum.tsv` and the number of supporting triplets in `introgression_matrix_inflow_outflow_supporting_count.tsv`.
- `introgression_ghost_target_strength.tsv` — per-taxon average ghost bootstrap support, plus a `has_sampled_introgression` flag (`1` when the taxon is also the target of a sampled introgression edge); `introgression_ghost_target_strength_raw_sum.tsv` and `introgression_ghost_target_strength_supporting_count.tsv` hold the sum and the count.
- `introgression_matrix_sampled_non_sister.tsv` — symmetric matrix counting, per taxon pair, the triplets in which the two are not the species-tree sisters.
- `introgression_taxa_order.tsv` — ordered taxa list matching the plot axes.

Returns an `IntrogressionMapArtifacts` dataclass with `plot_path`, the eight TSV paths, `taxa_count`, `non_ghost_edge_count`, and `ghost_target_count`.

### Bootstrap averaging

Each cell of the heatmap and each ghost bar is an undiluted average: the mean
`bootstrap_value` over the triplets whose classification produced that edge or
that ghost target, and nothing else. A triplet classified `no_introgression`
or `ambiguous`, or one whose edge points elsewhere, is not in the denominator.

- **Sampled introgression** for a directed pair (source → target): an
  `inflow_introgression` triplet contributes the edge from `C` to the
  discordant1 sister, and an `outflow_introgression` triplet the edge from
  that sister to `C`.

  `avg = sum(bootstrap_value over triplets that produced the edge) / count(those triplets)`

- **Ghost introgression** for a target taxon: a `ghost_introgression` triplet
  contributes the sister that the discordant1 topology leaves out.

  `avg = sum(bootstrap_value over triplets naming that ghost target) / count(those triplets)`

The `*_raw_sum.tsv` files hold the numerators and the `*_supporting_count.tsv`
files the denominators, so the averages can be reweighted against any other
denominator — such as the number of triplets in which a pair co-occurs, which
`introgression_matrix_sampled_non_sister.tsv` counts for every non-sister
pair. A run without a bootstrap (`bootstrap: false`) has no `bootstrap_value`,
and every classified triplet then enters the sums with a weight of 1, so every
supported edge and ghost target averages to `1` and the raw sums equal the
supporting counts.

### Plot layout

The combined figure uses a three-row layout above the data panels:

1. **Species tree strip** (top row) — topology-only tree with leaf labels suppressed.
2. **Source taxon label strip** (middle row) — a dedicated thin row containing the source-taxon names, rotated 90°, aligned to heatmap column centres. Row height is computed from the rendered pixel-width of the longest label so labels are never clipped. The labels are always drawn; the canvas is sized from the taxon count, so a larger tree produces a larger figure.
3. **Data panels** (bottom row, left to right):
    - **Inflow/outflow heatmap** — rows are target taxa, columns are source taxa, coloured by average bootstrap support on the `CONSOLIDATION_COLORMAP` (`cividis`) scale.
    - **Target label panel** — centred target taxon names aligned pixel-exactly to heatmap rows.
    - **Ghost bar chart** — horizontal bars per target taxon. Bar *length* is the average ghost bootstrap support. Bar *colour* is constant per bar, drawn from the two extremes of the same colormap, and encodes only whether the taxon also has sampled introgression: the high end, yellow (`GHOST_ONLY_BAR_COLOR`) when the taxon's only signal is ghost introgression, the low end, dark blue (`GHOST_WITH_SAMPLED_BAR_COLOR`) when it is also the target of a sampled introgression edge. Every bar carries a hairline `GHOST_BAR_EDGE_COLOR` outline so its extent stays legible against the panel. The flag is recorded in the ghost TSV's `has_sampled_introgression` column.
    - **Colorbar** — applies to the heatmap only. A two-entry legend sits directly above the bar panel explaining the ghost bar colours, and a note above that states that uncoloured heatmap cells carry no introgression and are not on the colour scale.

All three rows share `hspace=0` so they appear flush. Figure and panel widths scale dynamically with taxon count and rendered label widths.
The shared x/y labels and colorbar text scale with taxon count and are capped to stay readable on large figures, while the species-name labels keep their separate sizing.

### How it is invoked

Consolidation runs automatically as the orchestrator's last stage, writing the
figure into a `consolidation/` subfolder of the run's output directory and the
TSVs into `consolidation/consolidation_data/`. Disable it with
`--no-consolidation` / `consolidation: false`. It has no CLI of its own: its
input is the in-memory results list, not a file, so there is nothing to point a
command line at.
