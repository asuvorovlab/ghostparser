# ghostparser.orchestrator

`ghostparser.orchestrator` is GhostParser's introgression engine. It fuses triplet
subtree extraction and per-triplet inference into a single streaming pass, so
the intermediate triplet-gene-trees dataset is never written to disk or reloaded
into memory.

This document explains how the module works. The complete reference for every
flag and config key lives in [CONFIG.md](../../CONFIG.md#orchestrator-primary-module).

## Running it

```bash
# minimal run
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup

# multiple outgroups, custom output folder, all cores
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og Out1,Out2 \
    --output-folder results --processes 0

# a filtered set of triplets, no consolidation plots
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup \
    --triplet-filter triplets.txt --no-consolidation

# config-file mode (JSON or YAML); other CLI flags are ignored
python -m ghostparser.orchestrator -c run_config.yaml

# check the input data and exit, without running any analysis
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup --preflight-data-check
```

## Configuration

Three inputs are required — the species tree (`-st`), the gene trees (`-gt`),
and the outgroup(s) (`-og`) — and everything else has a default.
`-c/--config-file` is the only CLI-only option; when given, the file supplies
every setting and the other CLI flags are ignored with a warning.

A handful of settings are config-file-only (`discordant_test`,
`tree_height_calculation_strategy`, `min_support_value`,
`generate_summary_stats`, and the `bootstrap_options` block). Two defaults are
worth calling out: `p_value_correction` defaults to `bfn` (Bonferroni) and
`summary_statistic` defaults to `mean`.

See [CONFIG.md](../../CONFIG.md#orchestrator-primary-module) for every key, its
default, and its allowed values.

## Preflight data check

`--preflight-data-check` (or `preflight_data_check: true`) turns the run into a
data validation pass. `runner._run_preflight_only` short-circuits
`run_orchestrator` immediately after the output directory is prepared, so no
analysis runs and the only artifact is `preflight_data_check.txt`.

`preflight.run_preflight_data_check` replays the same structural logic the
engine uses, but collects failures instead of raising on the first one:

1. Root the species tree on the outgroup and normalize each triplet to A/B/C
   via `find_sister_pair` and `normalize_abc_from_sister_pair`, recording
   triplets whose rooted topology cannot be resolved.
2. Root every gene tree with `trees._root_tree_on_any_outgroup`, recording the
   ones where no outgroup label is present.
3. For each triplet contained in a gene tree, extract the subtree and replay
   `triplet_taxa_labels` → `find_sister_pair` → `topology_from_sister_pair` →
   sister-pair MRCA lookup, recording whichever step fails.

Every failure becomes an `Issue` with a dotted category. The report groups them
by category with counts, up to 25 examples each naming the gene-tree index and
triplet plus the offending input line, and an attribution summary separating
species-tree causes from gene-tree causes. `PreflightResult.passed` is `True`
only when nothing was detected.

Three conditions make the check itself impossible and raise `ValueError`
instead: no outgroups given, a species-tree file that does not hold exactly one
tree, and a species tree containing none of the outgroups.

## Orchestrator summary

`runner.run_orchestrator(config)` coordinates the run:

1. **Species preprocessing** — `trees.clean_and_save_trees` standardizes the
   species tree and drops trees whose mean internal support is below
   `min_support_value`. `trees._root_tree_on_outgroup` roots on the outgroup
   MRCA and prunes the outgroup, returning the ingroup taxa.
2. **Triplet setup** — `trees.generate_triplets` enumerates every ingroup
   triplet (or `trees.read_triplet_filter_file` plus
   `trees.filter_triplets_by_taxa` restricts them).
   `trees._build_species_triplet_metadata` normalizes each triplet to
   `(A, B, C)` with A and B the species-tree sisters, and builds the triplet's
   species subtree.
3. **Gene-tree preprocessing** — `trees.clean_and_save_gene_trees` cleans each
   gene tree and roots it on the outgroup.
4. **Fused extraction + inference** — `stream.stream_triplet_results` walks the
   triplets, extracts each one's subtree per gene tree, converts it directly to
   an observation (`inference.observation_from_subtree`), and immediately runs
   `inference.analyze_triplet_from_observations`. Only the small result object
   is retained; the subtrees are discarded.
5. **Run-wide correction** —
   `inference._apply_triplet_result_p_value_correction` applies the
   multiple-testing correction once across all triplets, because a global
   correction needs every p-value in a single pass.
6. **Writing** — `inference.write_pipeline_results` emits
   `orchestrator_triplet_results.tsv`; `inference.write_summary_statistics_tsv`
   emits `summary_statistics.tsv` when `generate_summary_stats` is set.
7. **Consolidation** — `introgression_mapper.generate_introgression_maps`
   writes the map artifacts into a `consolidation/` subfolder.

## Per-triplet inference

For each triplet the engine classifies every gene tree's subtree into one of
three topologies — concordant (matching the species tree) plus two discordant
alternatives — and records a tree height H(T) per the configured strategy. It
then applies a three-gate decision:

1. **Discordant count test (DCT)** — compares the two discordant counts
   (`inference.run_discordant_count_test`, chi-square or z-test). If the
   corrected p-value is not below `alpha_dct`, the triplet is
   `no_introgression` and the remaining gates are skipped.
2. **Tree-height test** — a two-sample KS test between the concordant and
   discordant1 height distributions (`inference.run_two_sample_ks_test`). If it
   is *not* significant, the triplet is `inflow_introgression`.
3. **Summary comparison** — otherwise the configured `summary_statistic` over
   the concordant heights is compared with the same statistic over the
   discordant1 heights: con > dis gives `outflow_introgression`, con < dis gives
   `ghost_introgression`, and equal or missing gives `unresolved`.

`inference._classify_introgression` implements this decision table directly.

Bootstrap resampling (on by default) repeats the analysis over resampled
observations and aggregates the per-iteration classifications into
`bootstrap_value`.

## Outputs

Written under the output folder:

- `orchestrator_triplet_results.tsv` — one row per triplet.
- `summary_statistics.tsv` — only when `generate_summary_stats` is set;
  per-triplet topology/metric summary statistics (63 metric columns covering
  mean/median/mode/variance/entropy/min/max over avg-tree-height/internal-branch/
  sister-distance for concordant/discordant1/discordant2).
- `processed_<species tree>` / `processed_<gene trees>` — cleaned, rooted trees.
- `metrics.txt` — per-stage wall/CPU timing and run parameters.
- `consolidation/` — the combined heatmap/bar-chart plot and TSV matrices from
  `introgression_mapper`, unless `--no-consolidation` is given. Consolidation
  writes into this dedicated subfolder so its own output-directory reset never
  removes the run folder's results TSV, processed trees, or the open
  `metrics.txt`.

### How each results column is produced

| Column | Source | Method |
| --- | --- | --- |
| `triplet` | Triplet setup | The normalized `(A, B, C)` labels, A and B being the species-tree sisters. |
| `species_tree` | Triplet setup | The triplet's species subtree, serialized topology-only (branch lengths omitted). |
| `dis1_topology` | Topology ranking | `BC` or `AC` — whichever discordant topology is more frequent; ties resolve to the first listed. |
| `most_frequent_matches_concordant` | Topology counts | True when the concordant count is at least both discordant counts. |
| `n_con` / `n_dis1` / `n_dis2` | Topology counts | Gene trees observed with each topology. |
| `analyzed_trees` | Extraction | Gene trees from which a subtree for this triplet was extracted. |
| `dct_statistic` / `dct_p_value` | DCT | SciPy chi-square or statsmodels z-test over `[n_dis1, n_dis2]`. An all-zero discordant split short-circuits to `(0.0, 1.0)`. |
| `dct_p_value_<method>_corr` | Correction | Run-wide correction over every triplet's DCT p-value. |
| `dct_significant` | Decision gate 1 | Corrected DCT p-value below `alpha_dct`. |
| `ks_statistic` / `ks_p_value` | Tree-height test | Two-sample KS between concordant and discordant1 heights; an empty sample yields `(0.0, 1.0)`. |
| `ks_p_value_<method>_corr` | Correction | Run-wide correction over every triplet's KS p-value. |
| `ks_significant` | Decision gate 2 | Corrected KS p-value below `alpha_ks`. |
| `summary_con` / `summary_dis` | Decision gate 3 | The configured statistic over the concordant / discordant1 heights. |
| `classification` | Decision logic | `no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, or `unresolved`. |
| `inference_description` | Reporting | Human-readable direction naming the actual species. |
| `bootstrap_value` / `all_bootstrap` | Bootstrap | Fraction of iterations agreeing with the final classification, plus the full class-fraction map. Present unless `--no-bootstrap`. |
| `bootstrap_*` debug columns | Bootstrap debug | Per-iteration DCT/KS statistics, con/dis summaries, and gene-tree heights. Present only with `bootstrap_debug_mode`. |

The results TSV is named `orchestrator_triplet_results.tsv`, distinct from the
inputs a standalone `introgression_mapper` run consumes, so both can share an
output folder without colliding.

## Parallelization modes

- `taxon` — triplets are split into chunks and dispatched across workers; each
  worker runs the fused extract-then-infer loop for its chunk over the shared
  gene-tree list.
- `gene` — triplets are processed serially in the parent; within a single
  triplet, per-gene-tree subtree extraction is parallelized across cores.
- `auto` — selects `gene` when the ingroup taxa count is below
  `AUTO_TAXA_SMALL_THRESHOLD` (15) or the gene-tree count exceeds
  `AUTO_GENE_TREES_THRESHOLD` (3500), otherwise `taxon`.

With one worker (or one triplet) the engine runs the fused loop serially in the
parent process regardless of mode.

## Internal design

### File layout

```
ghostparser/orchestrator/
  __init__.py    exports run_orchestrator
  __main__.py    python -m ghostparser.orchestrator entry point: main() wires parsing -> run_orchestrator
  config.py      orchestrator defaults/choices, validation, CLI parser, and CLI/config resolution
  trees.py       tree/triplet preprocessing
  inference.py   per-triplet inference + summary stats + result type + TSV writers
  stream.py      fused extract+infer streaming engine
  preflight.py   structural data check reached via --preflight-data-check
  runner.py      run_orchestrator coordinator (cleaning -> streaming -> correction -> writing -> consolidation)
  ORCHESTRATOR.md    this document
```

### Shared dependencies

The orchestrator owns its tree preprocessing, inference, and configuration. It
imports only four things from the rest of the package:

- `ghostparser.config` — the shared configuration trunk (`ConfigError`, path
  resolution, raw config-file loading, required-path validation, overwrite
  resolution, `prepare_output_directory`). Orchestrator-specific defaults, choices,
  and validators live in `orchestrator/config.py`, which is why the orchestrator can set
  its own `bfn`/`mean` defaults without affecting the ML subpackage.
- `ghostparser.cli_config` — the generic `resolve_cli_or_config_args` resolver
  implementing config-file-wins precedence.
- `ghostparser.triplet_utils` — pure topology helpers.
- `ghostparser.introgression_mapper` — `generate_introgression_maps` for the
  consolidation stage.

### Fused streaming engine

The per-triplet unit is `inference.analyze_triplet_from_observations(triplet,
observations, species_subtree, ...)`, which takes precomputed `(topology,
tree-height, metrics)` observations. The third element carries the per-tree
summary metrics and is `None` unless `generate_summary_stats` is enabled.

The processing unit is a chunk of triplets that share one parse pass over the
gene trees. For each parsed gene tree the engine extracts every in-chunk
triplet's subtree and computes its observation directly from the subtree object
(`inference.observation_from_subtree`) — there is **no** serialize-to-Newick and
reparse round trip. It then runs `analyze_triplet_from_observations` for each
triplet in the chunk and drops the observations. This bounds live memory to one
chunk while amortizing the DendroPy parse cost across the chunk's triplets. Gene
trees are loaded once in the parent and shared read-only to workers via a
fork/forkserver initializer.

Bootstrap resampling is vectorized with NumPy: per-triplet resample indices are
drawn with a seeded `numpy.random.Generator`, and topology counts and
per-topology height groups are computed with array operations. Bootstrap values
are deterministic under a fixed `bootstrap_seed` — the per-triplet seed is
derived from the run seed and the triplet, so every parallelization mode agrees
exactly.

### Memory rationale

- No intermediate triplet-gene-trees file is written or reloaded.
- Peak memory is roughly the shared gene-tree Newick list, plus one chunk of
  transient subtrees, plus the accumulating list of small result objects.
- Result objects must be accumulated because global p-value correction needs all
  p-values in a single pass.

### Consolidation interface contract

`generate_introgression_maps` reads only four fields from each result
(`triplet`, `classification`, `dis1_topology`, `bootstrap_value`) via its
duck-typed `_extract_result_fields`. `TripletPipelineResult` keeps those four
fields; treat field parity on them as a maintenance constraint.
