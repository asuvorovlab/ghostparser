# Configuration Guide

This guide is the complete reference for every GhostParser configuration key. It is organized around the main entry point, `ghostparser.pipeline`, followed by the machine-learning subpackage (`ghostparser.ml`) and the standalone consolidation CLI (`ghostparser.introgression_mapper`).

Both `ghostparser.pipeline` and the `ghostparser.ml` trainers support config files. For how each module works internally, see [ghostparser/pipeline/PIPELINE.md](ghostparser/pipeline/PIPELINE.md) and [ghostparser/ml/ML.md](ghostparser/ml/ML.md).

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
python -m ghostparser.pipeline -c configs/run.yaml
```

**Resolved paths:**

- `species_tree_path` -> `/home/user/project/../data/species.tree` -> `/home/user/data/species.tree`
- `gene_trees_path` -> `/home/user/data/genes.tree`
- `output_folder` -> `/scratch/results`

---

## Pipeline (Primary Module)

`ghostparser.pipeline` is the end-to-end entry point. It fuses tree preprocessing, triplet subtree extraction, and per-triplet inference into a single streaming pass, then optionally consolidates the results into introgression maps.

For how the module works internally, see [ghostparser/pipeline/PIPELINE.md](ghostparser/pipeline/PIPELINE.md).

### Run With a Config File

`-c/--config-file` is the only CLI-only option. It accepts a JSON or YAML file:

```bash
python -m ghostparser.pipeline -c run_config.yaml
python -m ghostparser.pipeline -c run_config.json
```

When a config file is given, **the file supplies every setting and the other CLI flags are ignored with a warning** (config wins). This is the only way to set the config-file-only keys listed below.

Minimal YAML:

```yaml
species_tree_path: data/species.tree
gene_trees_path: data/genes.tree
outgroups: OutGroup
output_folder: results
```

Fuller YAML showing the config-file-only keys and the nested bootstrap block:

```yaml
species_tree_path: data/species.tree
gene_trees_path: data/genes.tree
outgroups: Out1,Out2
output_folder: results
processes: 0
parallelization_mode: auto
alpha_dct: 0.05
alpha_ks: 0.05
p_value_correction: bfn
summary_statistic: mean
discordant_test: chi-square
tree_height_calculation_strategy: AVG
min_support_value: 0.5
generate_summary_stats: true
consolidation: true
bootstrap: true
bootstrap_options:
  iterations: 100
  seed: 42
  debug_mode: false
  summary_only: false
```

### Required Keys

These must be supplied either on the CLI or in the config file.

##### `species_tree_path`

- CLI: `-st, --species-tree-path`
- Species tree in Newick format.

##### `gene_trees_path`

- CLI: `-gt, --gene-trees-path`
- Gene trees in Newick format, one per line.

##### `outgroups`

- CLI: `-og, --outgroups`
- Outgroup taxon identifier(s). Accepts a comma-separated string (`Out1,Out2`) or a YAML/JSON list. The species tree is rooted and pruned on the outgroup MRCA; gene trees are rooted on the outgroup.
- The alias `outgroup` is also accepted in config files.

### Config + CLI Keys

Settable either on the CLI or in a config file.

##### `output_folder`

- CLI: `--output-folder`
- Default: `results`
- Directory for all run outputs.

##### `overwrite`

- CLI: `--no-overwrite` (sets `overwrite: false`)
- Default: `true`
- When `true`, an existing output directory is reset. When `false`, the run writes to an auto-suffixed sibling (`results_1`, `results_2`, ...) using the smallest missing suffix.

##### `triplet_filter`

- CLI: `--triplet-filter`
- Default: none (all triplets)
- Path to a file of comma-separated taxa triplets, one per line. Only the listed triplets are analyzed.

##### `processes`

- CLI: `--processes`
- Default: `0`
- Worker process count. `0` uses all available cores; `1` runs serially in the parent process.

##### `parallelization_mode`

- CLI: `--parallelization-mode`
- Default: `auto`
- Allowed: `auto`, `taxon`, `gene`
- `taxon` dispatches chunks of triplets across workers; `gene` parallelizes per-gene-tree subtree extraction within one triplet; `auto` selects `gene` when the ingroup taxa count is below 15 or the gene-tree count exceeds 3500, otherwise `taxon`.

##### `alpha_dct`

- CLI: `--alpha-dct`
- Default: `0.05`
- Significance threshold for the discordant count test (the first gate of the decision logic).

##### `alpha_ks`

- CLI: `--alpha-ks`
- Default: `0.05`
- Significance threshold for the KS tree-height test (the second gate).

##### `p_value_correction`

- CLI: `--p-value-correction`
- Default: `bfn`
- Allowed: `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`
- Multiple-testing correction applied once across every triplet in the run. Corrected p-values drive the significance decisions; the uncorrected values are retained in the output for reporting.

##### `summary_statistic`

- CLI: `--summary-statistic`
- Default: `mean`
- Allowed: `mean`, `median`, `mode`
- Statistic used to compare concordant against discordant1 tree heights in the final classification step. `mode` bins values to three decimals and resolves ties to the largest value.

##### `consolidation`

- CLI: `--no-consolidation` (sets `consolidation: false`)
- Default: `true`
- Generates the introgression map artifacts into a `consolidation/` subfolder of the output directory.

##### `bootstrap`

- CLI: `--no-bootstrap` (sets `bootstrap: false`)
- Default: `true`
- Enables bootstrap resampling per triplet and adds the `bootstrap_value` and `all_bootstrap` columns to the results TSV.

### Config-File-Only Keys

These have no CLI flag. They take their default unless set in a config file.

##### `discordant_test`

- Default: `chi-square`
- Allowed: `chi-square`, `z-test`
- The discordant count test: SciPy's Pearson chi-square, or a statsmodels two-proportion z-test.

##### `tree_height_calculation_strategy`

- Default: `AVG`
- Allowed: `AVG`, `A`, `B`, `C`, `SIS`, `INT`
- How H(T) is computed per gene tree: `AVG` averages the three root-to-tip distances; `A`/`B`/`C` take a single taxon's root-to-tip distance; `SIS` is the sister-pair patristic distance; `INT` is the sister-MRCA-to-root internal branch.

##### `min_support_value`

- Default: `0.5`
- Trees whose mean internal-node support falls below this threshold are dropped during cleaning. Trees without support labels are always kept.

##### `generate_summary_stats`

- Default: `false`
- Also writes `summary_statistics.tsv` with 63 metric columns (mean/median/mode/variance/entropy/min/max over avg-tree-height/internal-branch/sister-distance for concordant/discordant1/discordant2).

##### `bootstrap_options`

A nested block; each key may also be given flat as `bootstrap_<key>`.

- `iterations` (flat: `bootstrap_iterations`) — default `100`. Bootstrap iterations per triplet; must be an integer >= 1.
- `seed` (flat: `bootstrap_seed`) — default none. Base RNG seed; each triplet derives a deterministic per-triplet seed from it, so results are reproducible and independent of the parallelization mode.
- `debug_mode` (flat: `bootstrap_debug_mode`) — default `false`. Appends the per-iteration bootstrap-debug columns to the results TSV.
- `summary_only` (flat: `bootstrap_summary_only`) — default `false`. With debug mode on, emits compact summaries instead of full per-iteration lists.

### Pipeline CLI Example

```bash
python -m ghostparser.pipeline \
  -st data/species.tree \
  -gt data/genes.tree \
  -og Out1,Out2 \
  --output-folder results \
  --processes 0 \
  --parallelization-mode auto \
  --alpha-dct 0.05 \
  --alpha-ks 0.05 \
  --p-value-correction bfn \
  --summary-statistic mean
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

The ML sample configs illustrate the `input_path`, `output_dir`, `model`, `evaluation`, and `hyperparameter_tuning` sections that the ML loaders expect.

---

## Introgression Mapper (CLI Submodule)

`ghostparser.introgression_mapper` can be run independently to regenerate consolidation artifacts from an existing `pipeline_triplet_results.tsv` without re-running the full pipeline.

### CLI Options

#### Required

##### `-r`, `--results-tsv`

- Description: path to `pipeline_triplet_results.tsv` produced by the pipeline.

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
  --results-tsv results/pipeline_triplet_results.tsv \
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
- Config-file mode (`-c/--config-file`) is available in `ghostparser.pipeline` and the `ghostparser.ml` trainers.
- Use `--processes 1` to run single-worker mode.
