# Configuration Guide

This guide is organized around the main pipeline entry point, `ghostparser.orchestrator`, and then the three submodules (`tree_parser`, `triplet_processor`, `introgression_mapper`).
Only `ghostparser.orchestrator` supports config files.

The streaming module `ghostparser.pipeline` has its own CLI and config-file mode; its flags, config keys, and defaults are documented in [ghostparser/pipeline/PIPELINE.md](ghostparser/pipeline/PIPELINE.md).

## Path Resolution

GhostParser automatically resolves all path fields (input files, output directories, filter files) following standard operating system conventions:

### Path Types

#### 1. Absolute paths (start with `/`)

```yaml
species_tree_path: /home/user/data/species.tree
output_folder: /scratch/results
```

- Used as-is without modification
- Platform-independent representation

#### 2. Relative paths (no leading `/`)

```yaml
species_tree_path: data/species.tree
gene_trees_path: ./genes.tree
output_folder: results
```

- Resolved from the **current working directory** where the command is executed
- Example: If you run the command from `/home/user/project/`, then `data/species.tree` resolves to `/home/user/project/data/species.tree`

#### 3. User home directory (starts with `~`)

```yaml
species_tree_path: ~/data/species.tree
output_folder: ~/results
triplet_filter: ~/filters/triplets.txt
```

- `~` expands to your home directory (e.g., `/home/username/`)
- Example: `~/data/species.tree` becomes `/home/username/data/species.tree`

### Important Notes

- **Path resolution happens at runtime** when the config/CLI is parsed
- **All path types work in both CLI and config file modes**
- **Relative paths are NOT relative to the config file location** - they are relative to the directory where you execute the command
- For portability, consider using relative paths in configs and running commands from a consistent location

### Examples

**Config file at** `~/project/configs/run.yaml`:
```yaml
species_tree_path: ../data/species.tree    # Relative to execution directory, not config file
gene_trees_path: ~/data/genes.tree         # User home directory
output_folder: /scratch/results            # Absolute path
```

**Executed from** `/home/user/project/`:
```bash
python -m ghostparser.orchestrator -c configs/run.yaml
```

**Resolved paths:**

- `species_tree_path` -> `/home/user/project/../data/species.tree` -> `/home/user/data/species.tree`
- `gene_trees_path` -> `/home/user/data/genes.tree`
- `output_folder` -> `/scratch/results`

---

## Orchestrator (Primary Module)

`ghostparser.orchestrator` is the end-to-end entrypoint. It runs tree preprocessing/triplet extraction and then triplet inference in one pipeline.

### Run With a Config File

Supported config file formats:

- `.json`
- `.yaml`
- `.yml`

Run with config mode:

```bash
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

If `-c/--config-file` and other CLI arguments are provided together, config mode is used and other CLI arguments are ignored.

### Orchestrator Config Keys and CLI Options

#### Required keys

##### `species_tree_path`

- Type: string path
- Parallel CLI: `--species-tree-path` (alias: `-st`)
- Description: path to the species tree file in Newick format.

##### `gene_trees_path`

- Type: string path
- Parallel CLI: `--gene-trees-path` (alias: `-gt`)
- Description: path to the gene trees file in Newick format.

##### `outgroups` or `outgroup`

- Type:
  - `outgroups`: list of strings
  - `outgroup`: string
- Parallel CLI: `--outgroups` (alias: `-og`)
- Description: outgroup taxa used for rooting/pruning.
- Accepted forms:

```yaml
outgroups:
  - Taxon1
  - Taxon2
```

```yaml
outgroup: Taxon1
```

```yaml
outgroup: Taxon1,Taxon2
```

#### Optional keys

##### `triplet_filter`

- Type: string path
- Parallel CLI: `--triplet-filter`
- Description: optional path to a triplet filter file (comma-separated taxa per line).

##### `output_folder`

- Type: string path
- Parallel CLI: `--output-folder`
- Description: output directory for pipeline artifacts.
- Default: `./results` from the current working directory.

##### `overwrite`

- Type: boolean
- Parallel CLI: `--no-overwrite` disables overwrite when set on the CLI
- Description: controls whether an existing output folder is cleared before the pipeline writes artifacts.
- Default: `true`

##### `triplet_output_format`

- Type: string
- Parallel CLI: `--triplet-output-format`
- Allowed values: `parquet` (default), `txt`
- Description: extraction output format from the tree parser stage.
- Value meaning:
  - `parquet`: writes partitioned parquet output optimized for downstream processing and larger datasets.
  - `txt`: writes plain-text triplet blocks for manual inspection and debugging.

##### `parquet_partitions`

- Type: integer >= 0
- Parallel CLI: `--parquet-partitions`
- Description: number of hash partitions used when triplet output format is parquet.
- Default: `128`.

##### `parquet_compression`

- Type: string
- Parallel CLI: `--parquet-compression`
- Allowed values: `zstd` (default), `snappy`, `gzip`, `brotli`, `none`
- Description: parquet compression codec.
- Value meaning:
  - `zstd`: best general-purpose balance of compression ratio and speed.
  - `snappy`: faster compression/decompression with larger output files.
  - `gzip`: higher compression ratio with slower runtime.
  - `brotli`: high compression ratio, typically slower than `zstd` for this workflow.
  - `none`: no compression, largest files and fastest write path.

##### `processes`

- Type: integer >= 0
- Parallel CLI: `--processes`
- Description: worker count for extraction and inference.
- Notes:
  - `0` uses all available CPU cores.
  - `1` uses single-worker execution.

##### `generate_summary_stats`

- Type: boolean
- Parallel CLI: `--generate-summary-stats`
- Description: writes `summary_statistics.tsv` when enabled.
- Default: `false`.

##### `consolidation`

- Type: boolean
- Parallel CLI: `--no-consolidation` (disable switch)
- Description: controls the introgression consolidation stage that generates heatmap/bar-chart artifacts from in-memory triplet results.
- Default: `true`.
- Consolidation artifacts:
  - `introgression_combined.png` — single combined figure with inflow/outflow heatmap and ghost bar chart.
- The TSV artifacts are written to `consolidation_data/` inside the configured results directory.
  - `introgression_matrix_inflow_outflow.tsv`
  - `introgression_matrix_inflow_outflow_raw_sum.tsv`
  - `introgression_matrix_inflow_outflow_supporting_count.tsv`
  - `introgression_ghost_target_strength.tsv`
  - `introgression_ghost_target_strength_raw_sum.tsv`
  - `introgression_ghost_target_strength_supporting_count.tsv`
  - `introgression_taxa_order.tsv`
- Average bootstrap values in the artifacts use supporting-triplet denominators: the directed-pair average is `sum(bootstrap weights) / count(triplets that produced that directed edge)`, and the ghost-target average is `sum(bootstrap weights) / count(triplets that produced that ghost target)`.
- Outgroup taxa are automatically excluded from all consolidation plots and TSVs.

##### `min_support_value`

- Type: number
- Parallel CLI: `--min-support-value`
- Description: support threshold for species and gene tree cleaning.
- Default: `0.5`.

##### `discordant_test`

- Type: string
- Parallel CLI: `--discordant-test`
- Allowed values: `chi-square` (default), `z-test`
- Description: discordant count test used in triplet inference.
- Value meaning:
  - `chi-square`: Pearson chi-square test on discordant topology counts.
  - `z-test`: two-proportion z-test on discordant topology frequencies.

##### `summary_statistic`

- Type: string
- Parallel CLI: `--summary-statistic`
- Allowed values: `median` (default), `mean`, `mode`
- Description: statistic used for con/dis1 distributions after KS testing.
- Value meaning:
  - `median`: robust to outliers; preferred for skewed distributions.
  - `mean`: arithmetic average; sensitive to outliers.
  - `mode`: most frequent value after binning heights to 3 decimal places.

##### `stats_backend`

- Type: string
- Parallel CLI: `--stats-backend`
- Allowed values: `standard` (default), `custom`
- Description: backend used for DCT and KS calculations.
- Value meaning:
  - `standard`: SciPy/statsmodels implementations for statistical tests.
  - `custom`: GhostParser's internal implementations (experimental).

##### `tree_height_calculation_strategy`

- Type: string
- Parallel CLI: `--tree-height-calculation-strategy`
- Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`
- Description: strategy used to compute tree-height statistics.
- Value meaning:
  - `AVG`: mean root-to-tip distance across taxa A, B, and C.
  - `A`/`B`/`C`: root-to-tip distance for the selected taxon only.
  - `SIS`: sister-pair distance metric from the inferred sister taxa.
  - `INT`: internal branch distance from sister-pair MRCA to triplet root.

##### `p_value_correction`

- Type: string
- Parallel CLI: `--p-value-correction`
- Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`
- Description: multiple-testing correction for DCT and KS p-values.
- Value meaning:
  - `no`: no correction.
  - `bfn`: Bonferroni single-step family-wise error control.
  - `holm`: Holm step-down family-wise error control.
  - `fdr_bh`: Benjamini-Hochberg false discovery rate control.
  - `fdr_by`: Benjamini-Yekutieli FDR control for dependent tests.
  - `fdr_tsbh`: two-stage Benjamini-Hochberg FDR control.

##### `alpha_dct`

- Type: number
- Parallel CLI: `--alpha-dct`
- Description: p-value threshold for DCT.
- Default: `0.05`.

##### `alpha_ks`

- Type: number
- Parallel CLI: `--alpha-ks`
- Description: p-value threshold for KS.
- Default: `0.05`.

##### `bootstrap`

- Type: boolean
- Parallel CLI: `--no-bootstrap` (disable switch)
- Description: enables bootstrap sampling-with-replacement when true.
- Default: `true`.

##### `bootstrap_options`

- Type: object
- Parallel CLI:
  - `iterations` -> `--bootstrap-iterations`
  - `seed` -> `--bootstrap-seed`
  - `debug_mode` -> `--bootstrap-debug-mode`
  - `summary_only` -> `--bootstrap-summary-only`
- Description: nested bootstrap runtime options.
- Supported keys:
  - `iterations` (integer >= 1), default `100`
  - `seed` (integer, optional)
  - `debug_mode` (boolean), default `false`
  - `summary_only` (boolean), default `false`

##### `bootstrap_iterations`

- Type: integer >= 1
- Parallel CLI: `--bootstrap-iterations`
- Description: flat key alternative to `bootstrap_options.iterations`.

##### `bootstrap_seed`

- Type: integer
- Parallel CLI: `--bootstrap-seed`
- Description: flat key alternative to `bootstrap_options.seed`.

##### `bootstrap_debug_mode`

- Type: boolean
- Parallel CLI: `--bootstrap-debug-mode`
- Description: flat key alternative to `bootstrap_options.debug_mode`.

##### `bootstrap_summary_only`

- Type: boolean
- Parallel CLI: `--bootstrap-summary-only`
- Description: flat key alternative to `bootstrap_options.summary_only`.

### Orchestrator CLI Example

```bash
python -m ghostparser.orchestrator \
  --species-tree-path data/asuv21/species.tree \
  --gene-trees-path data/asuv21/gene_trees.tree \
  --outgroups Ephemera_danica,Isonychia_kiangsinensis \
  --triplet-filter data/asuv21/asuv_all/triplet_filter.txt \
  --output-folder results \
  --triplet-output-format parquet \
  --parquet-partitions 128 \
  --parquet-compression zstd \
  --processes 8 \
  --generate-summary-stats \
  --no-consolidation \
  --min-support-value 0.6 \
  --discordant-test z-test \
  --summary-statistic mean \
  --stats-backend custom \
  --tree-height-calculation-strategy B \
  --p-value-correction fdr_bh \
  --alpha-dct 0.02 \
  --alpha-ks 0.1 \
  --bootstrap-iterations 250 \
  --bootstrap-seed 42 \
  --bootstrap-debug-mode
```

## Machine Learning (ghostparser.ml)

The ML subpackage exposes explicit trainer modules. Invoke a trainer directly (for example `python -m ghostparser.ml.random_forest` or `python -m ghostparser.ml.multi_knn`). The random forest baseline is the main example path in this section.

The loaders treat the 6-bit label column as a multi-label target: each bit becomes one binary label, so the trainer can report both per-label scores and the stricter exact-match result for the whole bitstring.

Install the optional ML dependency set with `pip install .[ml]` when you want these trainers available; the core package can be installed without scikit-learn.

The loader now uses a strict layout:

- top-level keys for core run inputs and split/runtime controls
- `model` for trainer hyperparameters
- `evaluation` for metric selection and report/save toggles

The only trainer CLI flags are `-c/--config-file`, `-i/--input-path`, and `-o/--output-dir`. Other settings are config-file keys.

### Top-Level Config Keys

##### `input_path`

- Type: string
- Parallel CLI: `--input-path` (alias `-i`)
- Description: path to `summary_statistics.tsv` produced by the pipeline.

##### `output_dir`

- Type: string
- Parallel CLI: `--output-dir` (alias `-o`)
- Description: directory where trained model and metric artifacts will be written.

##### `overwrite`

- Type: boolean
- Parallel CLI: `--no-overwrite` disables overwrite when set on the trainer CLI
- Description: controls whether an existing trainer output directory is cleared before artifacts are written.
- Default: `true`

##### `target_column`

- Type: string
- Default: `class`
- Description: column name in the TSV containing the fixed-length binary target string. The loader reads that value as a string bitstring, expands it into one binary label per position for multi-label training, and expects six positions.

##### `test_size`

- Type: float in (0,1)
- Default: `0.2`
- Description: fraction of the dataset reserved for the hold-out test set.

##### `cv_folds`

- Type: int >= 1 or null
- Default: `5`
- Description: number of stratified cross-validation folds to run on the training partition.

##### `rare_class_policy`

- Type: string
- Default: `warn_reduce_cv`
- Description: behaviour when stratified CV is not feasible. Choices: `warn_reduce_cv`, `warn_skip_cv`, `error`.

##### `random_state`

- Type: int or null
- Default: `null`
- Description: optional RNG seed for deterministic splits and model behaviour.

##### `n_jobs`

- Type: int or null
- Default: `-1`
- Description: number of CPU worker jobs used by estimators. `-1` means all available CPU cores for operations that support parallelism; it does not enable GPU acceleration.

### Config Layout

The loader expects a top-level layout like this:

```yaml
input_path: ./results/summary_statistics.tsv
output_dir: ./results/ml_out
target_column: class
test_size: 0.2
cv_folds: 5
rare_class_policy: warn_reduce_cv
random_state: 42
n_jobs: -1

model:
  n_estimators: 200
  max_depth: 10
  min_samples_split: 2
  min_samples_leaf: 1
  max_features: sqrt
  class_weight: null
  n_neighbors: 5
  weights: uniform
  algorithm: auto
  leaf_size: 30
  metric: minkowski
  p: 2

evaluation:
  metrics: all
  report_class_distribution: true
  report_confusion_matrix: true
  report_feature_importance: true
  save_label_map: true
  save_predictions: true
```

### Model Parameters

Feature handling:

- The loader uses every non-target column as a feature.
- Numeric feature columns are used directly.
- String-valued feature columns with 7 or fewer distinct values are one-hot encoded automatically.
- String-valued feature columns with more than 7 distinct values are rejected.
- The configured `target_column` is excluded from the feature matrix automatically, and its raw TSV value is still read as the multi-label target bitstring.

If you want to avoid a string column being encoded, remove it from the TSV before calling the trainer.

- RandomForest (`ghostparser.ml.random_forest`):
  - `n_estimators` (int, default: `200`)
  - `max_depth` (int or null, default: `null`) — `null` leaves tree depth unconstrained.
  - `min_samples_split` (int, default: `2`) — controls how many samples are required before a split is allowed. Larger values make the trees more conservative when the data is noisy or small.
  - `min_samples_leaf` (int, default: `1`) — controls how many samples must remain in a leaf. Larger values smooth the model and can reduce noise.
  - `max_features` (string|int or null, default: `sqrt`) — `null` keeps the estimator's default split-feature behavior.
  - `class_weight` (null|dict, default: `null`) — `null` disables class weighting.

  Leave `min_samples_split` and `min_samples_leaf` out of the config if you want the defaults. The loader does not infer them from the dataset, and explicit `null` values are rejected.

- Multi-label KNN (`ghostparser.ml.multi_knn`):
  - `n_neighbors` (int >= 1, default: `5`)
  - `weights` (string, default: `uniform`) — `uniform` or `distance`.
  - `algorithm` (string, default: `auto`) — `auto`, `ball_tree`, `kd_tree`, `brute`.
  - `leaf_size` (int >= 1, default: `30`)
  - `metric` (string, default: `minkowski`)
  - `p` (int >= 1, default: `2`)

The top-level `n_jobs` and `random_state` keys apply to both trainers.

### Evaluation Parameters

- `metrics`: string metric-set selector. Choices: `all`, `primary`, `diagnostic`, `per_bit`. Default: `all`.
- `report_class_distribution`: boolean, default `true`. When enabled, the text report includes the dataset summary block.
- `report_confusion_matrix`: boolean, default `true`.
- `report_feature_importance`: boolean, default `true`.
- `save_label_map`: boolean, default `true`. The label map is embedded in the overall metrics JSON.
- `save_predictions`: boolean, default `true`. The prediction TSV includes the matched-bit count.

Use `metrics: all` when you want both per-label metrics and exact-match accuracy in the same run. The `diagnostic` set is the strict whole-bitstring view, while `primary` and `per_bit` expose narrower slices of the same evaluation.

### Hyperparameter Tuning Parameters

The `hyperparameter_tuning` section configures `python -m ghostparser.ml.hyper_tune`. It is separate from `model` and `evaluation` so tuning stays explicit. The tuner config only accepts the runtime keys listed above plus `hyperparameter_tuning`; do not provide `model`, `evaluation`, or model hyperparameters at the top level. Those sections belong to the trainer config, not the tuner config.

Inside `hyperparameter_tuning`, the following keys are expected:

- `model` (string, default `random_forest`): tuner target. Choices: `random_forest`, `multi_knn`.
- `method` (string, default `grid`): `grid` or `random`.
- `objective` (string, default `exact_match_accuracy`): metric used to rank candidates. Choices: `exact_match_accuracy`, `hamming_loss`, `bitwise_accuracy`, `micro_f1`, `macro_f1`, `weighted_f1`.
- `top_k` (int, default `10`): number of top candidates to include in the text summary.
- `n_iter` (int, default `20`): number of sampled candidates when `method: random`.
- `max_candidates` (int, default `5000`): hard cap for full grid evaluation.
- `wandb_detailed_payloads` (bool, default `false`): when `true`, log additional per-candidate CV payloads (aggregate and fold-level JSON) and extra summary JSON blobs to Weights & Biases. Keep `false` when network/storage overhead matters.
- `search_space` (mapping): model hyperparameter candidates. Each parameter should map to a list of values. Omit a parameter from `search_space` if you want the trainer default to apply during tuning.

Allowed `search_space` keys depend on `model`:

- `random_forest`: `n_estimators`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `max_features`, `class_weight`
- `multi_knn`: `n_neighbors`, `weights`, `algorithm`, `leaf_size`, `metric`, `p`

Do not place runtime fields such as `input_path`, `output_dir`, `target_column`, `test_size`, `cv_folds`, `rare_class_policy`, `random_state`, or `n_jobs` inside `search_space`.

Example:

```yaml
hyperparameter_tuning:
  model: random_forest
  method: random
  objective: exact_match_accuracy
  top_k: 5
  n_iter: 20
  max_candidates: 5000
  search_space:
    n_estimators: [100, 200, 400]
    max_depth: [null, 10, 20]
    min_samples_split: [2, 5]
    min_samples_leaf: [1, 2]
    max_features: [sqrt, log2]
    class_weight: [null]
```

Use `method: grid` to evaluate every combination in the search space. Use `method: random` when you want to sample a fixed number of combinations from a larger space.

While it runs, the tuner prints progress to the console: it announces the search method, reports the total candidate cases it will evaluate, estimates the total model fits implied by CV, and logs per-candidate timing updates.

### CLI examples

Run RandomForest via module entrypoint:

```bash
python -m ghostparser.ml.random_forest -i ./results/summary_statistics.tsv -o ./results/ml_rf_out
```

Run Multi-KNN via module entrypoint:

```bash
python -m ghostparser.ml.multi_knn -i ./results/summary_statistics.tsv -o ./results/ml_knn_out
```

### Sample Config Files

- `sample_configs/random_forest_minimal.yaml`
- `sample_configs/multi_knn_minimal.yaml`
- `sample_configs/hyperparameter_tuning_random_forest.yaml`
- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`

The ML sample configs illustrate the `input_path`, `output_dir`, `model`, `evaluation`, and `hyperparameter_tuning` sections that the ML loaders expect.

## Tree Parser (CLI Submodule)

`ghostparser.tree_parser` performs cleaning and triplet extraction only.

### CLI Options

#### Required

##### `--species-tree-path` (`-st`)

- Parallel orchestrator key/CLI: `species_tree_path` / `--species-tree-path`
- Description: species tree input path.

##### `--gene-trees-path` (`-gt`)

- Parallel orchestrator key/CLI: `gene_trees_path` / `--gene-trees-path`
- Description: gene trees input path.

##### `--outgroups` (`-og`)

- Parallel orchestrator key/CLI: `outgroups` or `outgroup` / `--outgroups`
- Description: comma-separated outgroup taxa.

#### Optional

##### `--triplet-filter`

- Parallel orchestrator key/CLI: `triplet_filter` / `--triplet-filter`
- Description: optional triplet filter path.

##### `--output-folder`

- Parallel orchestrator key/CLI: `output_folder` / `--output-folder`
- Description: output directory for extracted artifacts.

##### `--triplet-output-format`

- Parallel orchestrator key/CLI: `triplet_output_format` / `--triplet-output-format`
- Allowed values: `parquet` (default), `txt`
- Value meaning:
  - `parquet`: partitioned columnar output for scalable downstream processing.
  - `txt`: plain-text output for quick inspection.

##### `--parquet-partitions`

- Parallel orchestrator key/CLI: `parquet_partitions` / `--parquet-partitions`
- Description: parquet partition count.

##### `--parquet-compression`

- Parallel orchestrator key/CLI: `parquet_compression` / `--parquet-compression`
- Allowed values: `zstd` (default), `snappy`, `gzip`, `brotli`, `none`
- Value meaning: same codec semantics as orchestrator `parquet_compression`.

##### `--min-support-value`

- Parallel orchestrator key/CLI: `min_support_value` / `--min-support-value`
- Description: support threshold for tree cleaning.

##### `--processes`

- Parallel orchestrator key/CLI: `processes` / `--processes`
- Description: worker count for extraction (`0` uses all cores).

##### `--no-multiprocessing`

- Parallel orchestrator key/CLI: no direct config key
- Description: disables multiprocessing in tree parser execution.

### Tree Parser CLI Example

```bash
python -m ghostparser.tree_parser \
  --species-tree-path data/asuv21/species.tree \
  --gene-trees-path data/asuv21/gene_trees.tree \
  --outgroups Ephemera_danica,Isonychia_kiangsinensis \
  --triplet-filter data/asuv21/asuv_all/triplet_filter.txt \
  --output-folder results \
  --triplet-output-format parquet \
  --parquet-partitions 128 \
  --parquet-compression zstd \
  --min-support-value 0.6 \
  --processes 8
```

## Triplet Processor (CLI Submodule)

`ghostparser.triplet_processor` runs inference on precomputed triplet extraction output.

### CLI Options

#### Required

##### `--input-path`

- Parallel orchestrator key/CLI: no direct key (orchestrator wires this stage internally)
- Description: path to triplet extraction output (`.parquet` dataset or `.txt`).

#### Optional

##### `--input-format`

- Parallel orchestrator key/CLI: no direct key
- Allowed values: `auto` (default), `txt`, `parquet`
- Description: input parser mode.
- Value meaning:
  - `auto`: infer format from input path.
  - `txt`: force text parser.
  - `parquet`: force parquet parser.

##### `--output-path`

- Parallel orchestrator key/CLI: no direct key
- Description: output TSV path.

##### `--stats-output`

- Parallel orchestrator key/CLI: no direct key
- Description: optional JSON statistics output path.

##### `--alpha-dct`

- Parallel orchestrator key/CLI: `alpha_dct` / `--alpha-dct`
- Description: DCT significance threshold.

##### `--alpha-ks`

- Parallel orchestrator key/CLI: `alpha_ks` / `--alpha-ks`
- Description: KS significance threshold.

##### `--discordant-test`

- Parallel orchestrator key/CLI: `discordant_test` / `--discordant-test`
- Allowed values: `chi-square` (default), `z-test`
- Value meaning: same method semantics as orchestrator `discordant_test`.

##### `--summary-statistic`

- Parallel orchestrator key/CLI: `summary_statistic` / `--summary-statistic`
- Allowed values: `median` (default), `mean`, `mode`
- Value meaning: same statistic semantics as orchestrator `summary_statistic`.

##### `--stats-backend`

- Parallel orchestrator key/CLI: `stats_backend` / `--stats-backend`
- Allowed values: `standard` (default), `custom`
- Value meaning: same backend semantics as orchestrator `stats_backend`.

##### `--tree-height-calculation-strategy`

- Parallel orchestrator key/CLI: `tree_height_calculation_strategy` / `--tree-height-calculation-strategy`
- Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`
- Value meaning: same strategy semantics as orchestrator `tree_height_calculation_strategy`.

##### `--p-value-correction`

- Parallel orchestrator key/CLI: `p_value_correction` / `--p-value-correction`
- Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`
- Value meaning: same correction-method semantics as orchestrator `p_value_correction`.

##### `--no-bootstrap`

- Parallel orchestrator key/CLI: `bootstrap` / `--no-bootstrap`
- Description: disable bootstrap sampling-with-replacement.

##### `--bootstrap-iterations`

- Parallel orchestrator key/CLI: `bootstrap_options.iterations` / `--bootstrap-iterations`
- Description: iterations per triplet.

##### `--bootstrap-seed`

- Parallel orchestrator key/CLI: `bootstrap_options.seed` / `--bootstrap-seed`
- Description: optional random seed.

##### `--bootstrap-debug-mode`

- Parallel orchestrator key/CLI: `bootstrap_options.debug_mode` / `--bootstrap-debug-mode`
- Description: enable debug bootstrap columns.

##### `--bootstrap-summary-only`

- Parallel orchestrator key/CLI: `bootstrap_options.summary_only` / `--bootstrap-summary-only`
- Description: compact bootstrap debug summaries.

##### `--processes`

- Parallel orchestrator key/CLI: `processes` / `--processes`
- Description: worker count for triplet analysis (`0` uses all cores).

##### `--generate-summary-stats`

- Parallel orchestrator key/CLI: `generate_summary_stats` / `--generate-summary-stats`
- Description: write `summary_statistics.tsv`.

##### `--no-multiprocessing`

- Parallel orchestrator key/CLI: no direct config key
- Description: disable multiprocessing in triplet processor execution.

### Triplet Processor CLI Example

```bash
python -m ghostparser.triplet_processor \
  --input-path results/unique_triplets_gene_trees.parquet \
  --input-format auto \
  --output-path results/triplet_introgression_results.tsv \
  --stats-output results/triplet_introgression_results.json \
  --alpha-dct 0.05 \
  --alpha-ks 0.05 \
  --discordant-test chi-square \
  --summary-statistic median \
  --stats-backend standard \
  --tree-height-calculation-strategy AVG \
  --p-value-correction no \
  --bootstrap-iterations 100 \
  --processes 8
```

---

## Introgression Mapper (CLI Submodule)

`ghostparser.introgression_mapper` can be run independently to regenerate consolidation artifacts from an existing `orchestrator_triplet_results.tsv` without re-running the full pipeline.

### CLI Options

#### Required

##### `-r`, `--results-tsv`

- Description: path to `orchestrator_triplet_results.tsv` produced by the orchestrator.

##### `-st`, `--species-tree-path`

- Description: path to the processed species tree file (Newick) used for taxon ordering.

##### `-o`, `--output-dir`

- Description: directory to write output plots and TSVs.

#### Optional

##### `-og`, `--outgroups`

- Description: comma-separated outgroup taxon names to exclude from all plots and TSVs.
- Example: `--outgroups OutGroup1,OutGroup2`

### Introgression Mapper CLI Example

```bash
python -m ghostparser.introgression_mapper \
  --results-tsv results/orchestrator_triplet_results.tsv \
  --species-tree-path results/processed_species.tree \
  --output-dir results/ \
  --outgroups Ephemera_danica,Isonychia_kiangsinensis
```

### Outputs

- `introgression_combined.png` — combined inflow/outflow heatmap and ghost target-strength bar chart.
- `introgression_matrix_inflow_outflow.tsv` — target × source matrix of average bootstrap support values.
- `introgression_ghost_target_strength.tsv` — per-taxon average ghost bootstrap support.
- `introgression_taxa_order.tsv` — ordered taxa list matching the plot axes.

---

## Notes

- Path values are resolved at runtime to absolute paths.
- Config-file mode (`-c/--config-file`) is available in orchestrator only.
- Use `--processes 1` to run single-worker mode.
