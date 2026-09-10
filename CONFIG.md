# Configuration Guide

This guide is the complete reference for every GhostParser configuration key. It is organized around the main entry point, `ghostparser.orchestrator`, followed by the machine-learning subpackage (`ghostparser.ml`).

Both `ghostparser.orchestrator` and the `ghostparser.ml` trainers support config files. For how each module works internally, see [ghostparser/orchestrator/ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md) and [ghostparser/ml/ML.md](ghostparser/ml/ML.md).

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
- **Keep `output_folder` separate from the directories holding your input trees.** With `overwrite: true` the output directory is reset before the run, which deletes whatever is already in it — including a species tree or gene-tree file sitting there. Point `output_folder` at a directory of its own.

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

`ghostparser.orchestrator` is the end-to-end entry point. It fuses tree preprocessing, triplet subtree extraction, and per-triplet inference into a single streaming pass, then optionally consolidates the results into introgression maps.

For how the module works internally, see [ghostparser/orchestrator/ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md).

### Run With a Config File

`-c/--config-file` is the only CLI-only option. It accepts a JSON or YAML file:

```bash
python -m ghostparser.orchestrator -c run_config.yaml
python -m ghostparser.orchestrator -c run_config.json

# or start from a shipped sample
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

When a config file is given, **the file supplies every setting and the other CLI flags are ignored with a warning** (config wins). This is the only way to set the config-file-only keys listed below.

Minimal YAML:

```yaml
species_tree_path: data/species.tree
gene_trees_path: data/genes.tree
outgroup: OutGroup
output_folder: results
```

Fuller YAML showing the config-file-only keys and the nested blocks:

```yaml
species_tree_path: data/species.tree
gene_trees_path: data/genes.tree
outgroup: Out1,Out2
output_folder: results
processes: 0
seed: 42
alpha_dct: 0.05
alpha_ks: 0.05
alpha_perm: 0.05
p_value_correction: bfn
discordant_test: chi-square
tree_height_calculation_strategy: AVG
min_support_value: 0.5
generate_summary_stats: true
shape_diagnostics: false
consolidation: true
bootstrap: true
bootstrap_options:
  iterations: 100
  debug_mode: false
  summary_only: false
permutation_options:
  min_resamples: 2500
  max_resamples: 25000
  ci_method: wilson
```

### Required Keys

These must be supplied either on the CLI or in the config file.

##### `species_tree_path`

- CLI: `-st, --species-tree-path`
- Species tree in Newick format.

##### `gene_trees_path`

- CLI: `-gt, --gene-trees-path`
- Gene trees in Newick format, one per line.

##### `outgroup`

- CLI: `-og, --outgroup`
- Outgroup taxon identifier(s). One key covers both the single- and multiple-outgroup cases: give a single label (`OutGroup`), a comma-separated string (`Out1,Out2`), or a YAML/JSON list (`["Out1", "Out2"]`). List entries may themselves be comma-separated. The species tree is rooted and pruned on the outgroup MRCA; gene trees are rooted on the outgroup.

### Config + CLI Keys

Settable either on the CLI or in a config file.

##### `output_folder`

- CLI: `--output-folder`
- Default: `results`
- Directory for all run outputs.
- Must not be a directory containing your input trees: with `overwrite: true` it is reset before the run and its existing contents are removed. Consolidation writes into a `consolidation/` subfolder of it, which is likewise reset.

##### `overwrite`

- CLI: `--no-overwrite` (sets `overwrite: false`)
- Default: `true`
- When `true`, an existing output directory is reset. When `false`, the run writes to an auto-suffixed sibling (`results_1`, `results_2`, ...) using the smallest missing suffix.

##### `triplet_filter`

- CLI: `--triplet-filter`
- Default: none (all triplets)
- Path to a file of comma-separated taxa triplets, one per line. Only the listed triplets are analyzed.

##### `species_rename_map`

- CLI: `--species-rename-map`
- Default: none (taxa keep their tree labels)
- Path to a map giving the name each taxon should appear under in the outputs.
  Accepts a two-column TSV (tree label, then display name, tab-separated; blank
  lines and `#` comments ignored) or a YAML mapping, chosen by file extension:

  ```tsv
  T1	Homo sapiens
  T2	Pan troglodytes
  ```

  ```yaml
  T1: Homo sapiens
  T2: Pan troglodytes
  ```

  The rename is applied to the species and gene trees as they are read, so every
  later stage uses the display names: the results TSV, `summary_statistics.tsv`,
  the processed tree files, and the consolidation matrices and plots. Taxa absent
  from the map keep their tree labels, so a partial map is fine.

  Because the rename happens first, the outgroup and any triplet-filter entries
  are matched against the display names as well. Supply those in the tree's own
  labels and they are mapped for you.

  The map is rejected if it maps a label more than once, or maps two labels onto
  the same display name — the latter would collapse two taxa into duplicate tree
  labels, which the Newick parser refuses further downstream.

##### `seed`

- CLI: `--seed`
- Default: none
- Base RNG seed for the whole run. Every random draw derives from it: the
  permutation direction test, the bootstrap resampling, the bootstrap's own
  permutations, and the modality bootstrap in the shape diagnostics. Each
  triplet derives its own stream from `(seed, triplet)`, so a run is
  reproducible at any worker count and independent of how many resamples any
  single adaptive test happens to draw.
- When omitted, a seed is drawn at run time instead. Either way the value used
  is written to `metrics.txt` as `Seed: <n> (configured|generated)`, so a run
  started without a seed can still be reproduced by passing back the value it
  reports.

##### `processes`

- CLI: `--processes`
- Default: `0`
- Worker process count. `0` uses all available cores; `1` runs serially in the parent process.

##### `alpha_dct`

- CLI: `--alpha-dct`
- Default: `0.05`
- Significance threshold for the discordant count test (the first gate of the decision logic).

##### `alpha_ks`

- CLI: `--alpha-ks`
- Default: `0.05`
- Significance threshold for the KS tree-height test (the second gate).

##### `alpha_perm`

- CLI: `--alpha-perm`
- Default: `0.05`
- Significance threshold for the studentized permutation test (the third gate). Applied to each of the two one-tailed p-values after they are corrected against each other, and to the two-tailed cross-check.

##### `p_value_correction`

- CLI: `--p-value-correction`
- Default: `bfn`
- Allowed: `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`
- Multiple-testing correction, applied in three places. Run-wide, it adjusts every triplet's DCT and KS p-value in a single pass. Inside each permutation test, it adjusts that test's pair of one-tailed p-values against each other — never across triplets, because a Monte Carlo p-value has a resolution floor that across-triplet correction would fall through. Inside the bootstrap, each iteration's DCT and KS p-values are corrected across triplets for that iteration index, so iterations answer to the same thresholds the reported classification does. Corrected p-values drive the significance decisions; the uncorrected values are retained in the output for reporting. Under `no` the results TSV carries the raw p-values and the significance flags only.
- The choice affects run time. `no` and `bfn` are applied inline while triplets stream; the others must hold every triplet's per-iteration p-values until the stream finishes.
- Every supported method is monotone — none can adjust a p-value *below* its raw value — which is what lets both the point estimate and the bootstrap skip a test once an earlier gate has failed. See [`pipeline_mode`](#pipeline_mode).
- In YAML, `p_value_correction: no` may be written with or without quotes. YAML resolves the bare word `no` to a boolean, and enumerated fields map booleans back to the choice they spell (`no`/`off`/`n`/`false`, `yes`/`on`/`y`/`true`), so both forms select the same value.

##### `pipeline_mode`

- CLI: `--pipeline-mode`
- Default: `efficient`
- Allowed: `efficient`, `detailed`
- `efficient` stops measuring a triplet once the decision cascade is settled: a triplet the count gate settled skips the tree-height test, and one either of the first two gates settled skips the permutation direction test, rather than computing a result nothing reads. `detailed` runs all three gates for every triplet.
- **Both modes produce identical results.** No supported correction can lower a p-value, so a gate that failed on the raw value cannot clear on the corrected one — the efficient mode only ever declines a test the cascade could not have consulted. Neither does the extra work change any correction family: the tree-height family is always the triplets that cleared the count gate, so whatever `detailed` measures below a settled gate is never enrolled. `classification`, `decision_gate`, the `dct_*` columns, `ks_p_value_corrected`, `ks_significant`, and every bootstrap column match exactly.
- What differs is which columns are populated, never their values. Triplets the efficient mode settled early leave the `perm_*` block empty and carry `perm_note: direction_test_not_consulted`, which distinguishes a deliberate skip from a test that ran and hit a guard; they also leave the raw `ks_statistic`/`ks_p_value` empty, since the tree-height test below a settled count gate decides nothing and takes no part in its correction family either way.
- The direction test is the most expensive of the three, so not running it where it cannot matter is where the time goes. How much that is worth depends on how many of your triplets stop at an earlier gate.
- Choose `detailed` when you want the direction test's statistics for every triplet regardless of whether they decided anything, which is a debugging need rather than an analysis one.
- The bootstrap already skipped a settled gate per iteration in both modes; this key governs the point estimate.

##### `consolidation`

- CLI: `--no-consolidation` (sets `consolidation: false`)
- Default: `true`
- Generates the introgression map artifacts into a `consolidation/` subfolder of the output directory.

##### `bootstrap`

- CLI: `--no-bootstrap` (sets `bootstrap: false`)
- Default: `true`
- Enables bootstrap resampling per triplet and adds the `bootstrap_value` and `all_bootstrap` columns to the results TSV. Setting it false skips the iterations entirely, so `bootstrap_stat_ci_low`/`bootstrap_stat_ci_high` are empty too — that interval is a bootstrap percentile interval, not a permutation output.
- This is an instruction about what to compute, so it holds under `pipeline_mode: detailed` as well: `detailed` declines to skip work the cascade cannot consult, which is not the same as reinstating work you switched off. The same is true of `generate_summary_stats` and `shape_diagnostics`.

##### `preflight_data_check`

- CLI: `--preflight-data-check` (sets `preflight_data_check: true`)
- Default: `false`
- Runs only the structural sanity check on the species tree, gene trees, and triplets, writes `preflight_data_check.txt` into the output folder, and exits without any analysis. No results TSV, processed trees, `metrics.txt`, or consolidation artifacts are produced. Every other analysis key is ignored for that run.
- The report lists each detected issue by category (for example `gene_tree.rooting_failed`, `triplet.unresolved_rooted_sister_pair`), a count per category, up to 25 example messages naming the offending gene-tree index and triplet, and a species-tree-versus-gene-tree attribution summary.
- The checks are structural: they establish whether the data can be processed, not whether the result will be biologically meaningful.

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
- Also writes `summary_statistics.tsv` with 63 metric columns (mean/median/mode/variance/entropy/min/max over avg-tree-height/internal-branch/sister-distance for concordant/discordant1/discordant2). The `discordant1_*` columns describe whichever discordant topology is more frequent — the same group the `dis1_topology` column names and the statistical tests use — and `discordant2_*` the other one.
- The results TSV carries no per-group mean or median columns regardless of this setting; the introgression direction comes from the permutation p-values.

##### `shape_diagnostics`

- Default: `false`
- Appends fifteen columns describing the shape of each height group: a KDE mode count, a Silverman modality p-value, skewness, excess kurtosis, and a generalized-Pareto tail index. They are descriptive only and never affect a classification.
- They land in the results TSV as `con_*`/`dis1_*`/`dis2_*`, and nowhere else. `summary_statistics.tsv` never carries them: it is a feature matrix, and these columns are undefined for groups below their observation floors, so including them would leave holes in it.
- Off by default because the modality p-value is a smoothed bootstrap: it costs about 0.2s per group, so roughly 0.6s of extra CPU per triplet. On a large taxon set that dominates the run.
- Measured once per triplet from the observed heights. Bootstrap iterations do not recompute them.
- Groups with fewer than 20 observations are left empty, as are the tail indices of groups whose upper decile holds fewer than 10 points. See "Shape diagnostics" in the [orchestrator guide](ghostparser/orchestrator/ORCHESTRATOR.md) for how to read each column.

##### `bootstrap_options`

A nested block; each key may also be given flat as `bootstrap_<key>`.

- `iterations` (flat: `bootstrap_iterations`) — default `100`. Bootstrap iterations per triplet; must be an integer >= 1.
- `debug_mode` (flat: `bootstrap_debug_mode`) — default `false`. Appends the per-iteration bootstrap-debug columns to the results TSV.
- `summary_only` (flat: `bootstrap_summary_only`) — default `false`. With debug mode on, emits compact summaries instead of full per-iteration lists.

##### `permutation_options`

A nested block tuning the permutation test, the third decision gate. There is no `initial_batch` key — the first adaptive batch is always `min_resamples`, and each subsequent batch is 1.25x the previous one — and no `ci_level` key, since the interval is fixed at 95%.

- `min_resamples` — default `2500`. Size of the first batch and the minimum total permutations. Must be an integer >= 1. A triplet whose pooled sample admits fewer than this many distinct group assignments is skipped with an `insufficient_permutation_support` note, because its permutation distribution cannot resolve `alpha_perm`.
- `max_resamples` — default `25000`. Resample budget. Must be an integer >= `min_resamples`. Reaching it without the confidence interval excluding `alpha_perm` sets `perm_converged` to false and records the triplet in `metrics.txt`. It is the point at which the run stops asking for more rather than a hard cap: the batch that crosses it is drawn whole rather than trimmed, so `perm_n_resamples` can exceed it by up to one batch. Keep the two values a few multiples apart — with `max_resamples` close to `min_resamples` a single grown batch is comparable to the whole budget, so the overshoot is proportionally much larger.
- `ci_method` — default `wilson`. Binomial interval method passed to `statsmodels.stats.proportion.proportion_confint`. Allowed: `wilson`, `beta`, `agresti_coull`, `jeffreys`, `binom_test`, `normal`. `wilson` inverts the score test, stays inside [0, 1], and holds close-to-nominal coverage for the very small proportions this test produces; `normal` degrades badly there and `beta` (Clopper–Pearson) is guaranteed-coverage but conservative, so it resamples longer than necessary. An unrecognized value is rejected when the config is parsed. See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md#how-the-interval-is-computed-and-why-it-matches-the-p-value) for how the interval is derived and why it is consistent with the reported p-value.

Bootstrap iterations re-run the direction test at one fifth of `min_resamples` and `max_resamples`, since the bootstrap aggregate absorbs the extra per-iteration Monte Carlo noise.

### Orchestrator CLI Example

```bash
python -m ghostparser.orchestrator \
  -st data/species.tree \
  -gt data/genes.tree \
  -og Out1,Out2 \
  --output-folder results \
  --processes 0 \
  --alpha-dct 0.05 \
  --alpha-ks 0.05 \
  --alpha-perm 0.05 \
  --p-value-correction bfn
```

## Machine Learning (ghostparser.ml)

The ML subpackage exposes explicit trainer modules. Invoke a trainer directly (for example `python -m ghostparser.ml.random_forest` or `python -m ghostparser.ml.multi_knn`). The random forest baseline is the main example path in this section.

The loaders treat the 6-bit label column as a multi-label target: each bit becomes one binary label, so the trainer can report both per-label scores and the stricter exact-match result for the whole bitstring.

Install the optional ML dependency set with `pip install .[ml]` when you want these trainers available; the core package can be installed without scikit-learn.

The loader uses a strict layout:

- top-level keys for core run inputs and split/runtime controls
- `model` for trainer hyperparameters
- `evaluation` for metric selection and report/save toggles

The only trainer CLI flags are `-c/--config-file`, `-i/--input-path`, and `-o/--output-dir`. Other settings are config-file keys.

### Top-Level Config Keys

##### `input_path`

- Type: string
- Parallel CLI: `--input-path` (alias `-i`)
- Description: path to `summary_statistics.tsv` produced by the orchestrator.

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
  - `max_features` (string, int, float or null, default: `sqrt`) — how many features each split may consider. Accepts `sqrt`, `log2`, an integer `>= 1` (that many features per split), a float in `(0.0, 1.0]` (that fraction of the features), or `null` to use **every** feature at each split. `auto` is rejected: scikit-learn removed it in 1.3, and `sqrt` is its classifier equivalent.
  - `class_weight` (string, dict, list or null, default: `null`) — accepts `balanced`, `balanced_subsample`, a mapping of class label to weight, a list of such mappings (one per label), or `null` for no class weighting.

  In YAML, write the null value as `null` or `~`. The bare words `None` and `none` are read as plain strings, not null, so GhostParser maps them (and `null` written as a string) back to null for `max_features` and `class_weight` rather than passing the literal text to the estimator. Any other value is rejected with a message naming what it received and every accepted form.

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
- `use_wandb` (bool, **required**, no default): whether the run logs to Weights & Biases. There is deliberately no default — every tuning config states the choice, and the loader rejects a config that omits it or gives a non-boolean. With `false` the tuner needs neither a W&B account nor the `wandb` package, makes no network calls, writes no `wandb/` directory, and writes the full local artifact set. With `true` the `wandb` package must be installed (`pip install .[wandb]`) and authenticated, and the bulk outputs (`hyper_tune_results.json`, `hyper_tune_results.tsv`, `hyper_tune_parameter_marginals.tsv`, `predictions.tsv`) are logged to the W&B run instead of the output directory, which then keeps only the model pickle, the plaintext report, and the search-report plot.
- `wandb_detailed_payloads` (bool, default `false`): when `true`, log additional per-candidate CV payloads (aggregate and fold-level JSON) and extra summary JSON blobs to Weights & Biases. Requires `use_wandb: true`; combining it with `use_wandb: false` is rejected. Keep `false` when network/storage overhead matters.
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
  use_wandb: false
  search_space:
    n_estimators: [100, 200, 400]
    max_depth: [null, 10, 20]
    min_samples_split: [2, 5]
    min_samples_leaf: [1, 2]
    max_features: [sqrt, log2]
    class_weight: [null]
```

`max_features` and `class_weight` candidates are validated value by value against
the same rules as the `model` block above, so an invalid entry is reported against
its `search_space` key before the search starts rather than at the first fit.

To log the run to Weights & Biases instead, install and authenticate the extra
(`pip install .[wandb]`, then `wandb login`) and swap the two flags:

```yaml
hyperparameter_tuning:
  use_wandb: true
  wandb_detailed_payloads: true
```

Runnable samples for both formats ship as
`sample_configs/hyperparameter_tuning_random_forest.yaml` and
`sample_configs/hyperparameter_tuning_multi_knn.json`.

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

Orchestrator samples live alongside them:

- `sample_configs/orchestrator_minimal.yaml` — the three required inputs plus an output folder; everything else defaults.
- `sample_configs/orchestrator_full.yaml` — every key at its default value, including the config-file-only ones, as a starting point to trim.

---

## Consolidation Outputs

Consolidation is a stage of the orchestrator, not a separate entry point, and is
controlled by the [`consolidation`](#consolidation) key. It writes into a
`consolidation/` subfolder of the run's output folder:

- `introgression_combined.png` — combined inflow/outflow heatmap and ghost target-strength bar chart.
- `introgression_matrix_inflow_outflow.tsv` — target × source matrix of average bootstrap support values.
- `introgression_ghost_target_strength.tsv` — per-taxon average ghost bootstrap support, plus a `has_sampled_introgression` flag (`1` when that taxon is also the target of a sampled introgression edge) that sets the bar colour.
- `introgression_taxa_order.tsv` — ordered taxa list matching the plot axes.

Taxa named in [`outgroup`](#outgroup) are excluded from every plot and TSV here.

---

## Notes

- Path values are resolved at runtime to absolute paths.
- Config-file mode (`-c/--config-file`) is available in `ghostparser.orchestrator` and the `ghostparser.ml` trainers.
- Use `--processes 1` to run single-worker mode.
