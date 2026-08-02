# ghostparser.pipeline

`ghostparser.pipeline` is a self-contained streaming introgression pipeline. It
fuses triplet subtree extraction and per-triplet inference into a single pass so
the intermediate triplet-gene-trees dataset is never written to disk or reloaded
into memory.

## When to use it

Use the pipeline for large gene-tree sets where the orchestrator's two-stage
design (extract every triplet's subtrees to one intermediate file, then reload
that whole file for inference) exhausts memory. The pipeline extracts a triplet's
subtrees, runs its inference immediately, keeps only the small per-triplet result
object, and discards the subtrees — bounding live memory to one chunk of subtrees
plus the shared gene-tree list plus the accumulating results.

For parquet intermediate storage use `ghostparser.orchestrator`.

## Running it

```bash
# minimal run
python -m ghostparser.pipeline -st species.tree -gt genes.tree -og OutGroup

# multiple outgroups, custom output folder, all cores
python -m ghostparser.pipeline \
    -st species.tree -gt genes.tree -og Out1,Out2 \
    --output-folder results --processes 0

# a filtered set of triplets, no consolidation plots
python -m ghostparser.pipeline \
    -st species.tree -gt genes.tree -og OutGroup \
    --triplet-filter triplets.txt --no-consolidation

# config-file mode (JSON or YAML); other CLI flags are ignored
python -m ghostparser.pipeline -c run_config.yaml
```

### Config-file mode

`-c/--config-file` (CLI-only) points at a JSON or YAML file. When it is given,
the file supplies every setting and the other CLI flags are ignored with a
warning (config wins). The file may set both the config+CLI options and the
config-file-only options listed below, using the same key names, and accepts the
nested `bootstrap_options` block (`iterations`, `seed`, `debug_mode`,
`summary_only`).

### CLI flags

Required (in CLI mode; supplied by the file in config-file mode):

- `-st, --species-tree-path` — species tree in Newick format.
- `-gt, --gene-trees-path` — gene trees in Newick format.
- `-og, --outgroups` — outgroup taxon identifier(s), comma-separated.

Config+CLI options (settable on the CLI or in a config file):

| Flag | Config key | Default | Meaning |
| --- | --- | --- | --- |
| `--output-folder` | `output_folder` | `results` | Output directory. |
| `--triplet-filter` | `triplet_filter` | (none) | File of comma-separated taxa triplets, one per line. |
| `--no-overwrite` | `overwrite` | overwrite on | Write to an auto-suffixed sibling directory when the output folder exists. |
| `--processes` | `processes` | `0` | Worker processes; `0` uses all cores. |
| `--parallelization-mode` | `parallelization_mode` | `auto` | `taxon`, `gene`, or `auto` (see below). |
| `--alpha-dct` | `alpha_dct` | `0.05` | DCT significance threshold. |
| `--alpha-ks` | `alpha_ks` | `0.05` | KS significance threshold. |
| `--p-value-correction` | `p_value_correction` | `bfn` | Run-wide multiple-testing correction. |
| `--summary-statistic` | `summary_statistic` | `mean` | Statistic for the con-vs-dis1 comparison (`mean`/`median`/`mode`). |
| `--no-consolidation` | `consolidation` | consolidation on | Disable the introgression map/plot stage. |
| `--no-bootstrap` | `bootstrap` | bootstrap columns on | Omit the bootstrap columns from the results TSV. |

Config-file-only options (no CLI flag; default unless set in a config file):

| Config key | Default | Meaning |
| --- | --- | --- |
| `discordant_test` | `chi-square` | Discordant count test (`chi-square` or `z-test`). |
| `tree_height_calculation_strategy` | `AVG` | Tree-height strategy (`AVG`/`A`/`B`/`C`/`SIS`/`INT`). |
| `min_support_value` | `0.5` | Support threshold for tree cleaning. |
| `bootstrap_iterations` | `100` | Bootstrap iterations per triplet. |
| `bootstrap_seed` | (none) | Bootstrap RNG seed for reproducibility. |
| `generate_summary_stats` | `false` | Also write `summary_statistics.tsv`. |
| `bootstrap_debug_mode` | `false` | Add the bootstrap-debug columns to the results TSV. |
| `bootstrap_summary_only` | `false` | With debug mode, emit compact summaries instead of per-iteration lists. |

The statistical tests always use the scipy/statsmodels backend; there is no
selectable stats backend.

### Outputs

Written under the output folder:

- `pipeline_triplet_results.tsv` — one row per triplet with topology counts,
  DCT/KS statistics, classification, inference description, and (unless
  `--no-bootstrap`) the bootstrap columns. With `bootstrap_debug_mode` the
  bootstrap-debug columns (`bootstrap_dct_stats`, `bootstrap_dct_p_value`,
  `bootstrap_ks_stats`, `bootstrap_ks_p_value`, `bootstrap_con_<statistic>`,
  `bootstrap_dis_<statistic>`, `bootstrap_gene_tree_heights`) are appended.
- `summary_statistics.tsv` — written only when `generate_summary_stats` is set;
  per-triplet topology/metric summary statistics (63 metric columns covering
  mean/median/mode/variance/entropy/min/max over avg-tree-height/internal-branch/
  sister-distance for concordant/discordant1/discordant2).
- `processed_<species tree>` / `processed_<gene trees>` — cleaned, rooted trees.
- `metrics.txt` — per-stage wall/CPU timing and run parameters.
- `consolidation/` — a subfolder holding the consolidation artifacts (combined
  heatmap/bar-chart plot and TSV matrices) from `introgression_mapper`, unless
  `--no-consolidation` is given. Consolidation writes into this dedicated
  subfolder so its own output-directory reset never removes the run folder's
  results TSV, processed trees, or the open `metrics.txt`.

The results TSV is named `pipeline_triplet_results.tsv` (distinct from the
orchestrator's `orchestrator_triplet_results.tsv`) so a pipeline run and an
orchestrator run can share an output folder without colliding.

### Parallelization modes

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
ghostparser/pipeline/
  __init__.py    exports run_pipeline
  __main__.py    python -m ghostparser.pipeline entry point: main() wires parsing -> run_pipeline
  config.py      self-contained config foundation: defaults/choices, ConfigError, validation + config-file loading, prepare_output_directory, CLI/config resolution (build_argument_parser, resolve_config, load_pipeline_config, normalize_pipeline_payload)
  trees.py       PORTED tree/triplet preprocessing (from tree_parser)
  inference.py   PORTED per-triplet inference + summary stats + result type + TSV writers (from triplet_processor)
  stream.py      fused extract+infer streaming engine
  runner.py      run_pipeline coordinator (cleaning -> streaming -> correction -> writing -> consolidation)
  PIPELINE.md    this document
```

### Port-vs-import boundary

The pipeline is self-contained: it does **not** import from `orchestrator`,
`tree_parser`, `triplet_processor`, `config`, or `cli_config` (the modules slated
for deletion or restructuring once the pipeline is a proven replacement). The
logic it needs from `tree_parser`/`triplet_processor` is **ported** into
`trees.py` and `inference.py`, and the config/CLI foundation from
`config`/`cli_config` is **ported** into `config.py` (defaults and choices,
`ConfigError`, the config-file loading/validation helpers,
`prepare_output_directory`, and the CLI/config resolver). Ported logic is copied
and cleaned of the deferred/unused paths (intermediate-file machinery, parquet
I/O, the custom statistical backend, standalone CLIs). The computation of every
ported inference function is kept identical so results match the orchestrator's
scipy/statsmodels backend (verified by the parity tests in `tests/pipeline/`).

It **imports** only the two shared foundation modules that are not slated for
deletion:

- `ghostparser.triplet_utils` — pure topology helpers.
- `ghostparser.introgression_mapper` — `generate_introgression_maps` for the
  consolidation stage.

Because the pipeline owns its config defaults, its `p_value_correction` default
(`bfn`) and `summary_statistic` default (`mean`) can differ from the shared
orchestrator defaults without touching any shared module.

### Fused streaming engine

The per-triplet unit is `inference.analyze_triplet_from_observations(triplet,
observations, species_subtree, ...)`, which takes precomputed `(topology,
tree-height, metrics)` observations, runs the discordant count test and the KS
tree-height test, classifies the triplet, and aggregates a bootstrap value. The
third observation element carries the per-tree summary metrics and is `None`
unless `generate_summary_stats` is enabled.

The processing unit is a chunk of triplets that share one parse pass over the
gene trees. For each parsed gene tree the engine extracts every in-chunk
triplet's subtree and computes its observation directly from the subtree object
(`inference.observation_from_subtree`) — there is **no** serialize-to-Newick and
reparse round trip. It then runs `analyze_triplet_from_observations` for each
triplet in the chunk and drops the observations. This bounds live memory to one
chunk while amortizing the DendroPy parse cost across the chunk's triplets. Gene
trees are loaded once in the parent and shared read-only to workers via a
fork/forkserver initializer.

Because observations are computed from the (unrounded) subtree objects, tree
heights match the orchestrator's **parquet** observation path rather than its
Newick round-trip path.

The bootstrap resampling is vectorized with NumPy: per-triplet resample indices
are drawn with a seeded `numpy.random.Generator` and topology counts and
per-topology height groups are computed with array operations. Bootstrap values
are therefore statistically equivalent to, but not bit-for-bit identical with,
the orchestrator's `random.Random`-based bootstrap; they remain deterministic
under a fixed `bootstrap_seed` (the per-triplet seed is derived from the run seed
and the triplet), so parallelization modes agree exactly.

`stream.stream_triplet_results(...)` accumulates the per-triplet results and
applies the p-value correction once across all of them
(`_apply_triplet_result_p_value_correction`) — global correction requires every
p-value in one pass, matching the orchestrator's run-wide correction.

### Memory rationale

- No intermediate triplet-gene-trees file is written or reloaded.
- Peak memory ≈ shared gene-tree Newick list + one chunk of transient subtrees +
  the accumulating list of small result objects.
- Result objects must be accumulated because global p-value correction needs all
  p-values in a single pass.

### Consolidation interface contract

`generate_introgression_maps` reads only four fields from each result (`triplet`,
`classification`, `dis1_topology`, `bootstrap_value`) via its duck-typed
`_extract_result_fields`. The ported `TripletPipelineResult` keeps those four
fields; treat field parity on them as a maintenance constraint.
