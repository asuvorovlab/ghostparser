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
species_tree_path: data/species.tree     # Newick species tree file path
gene_trees_path: data/genes.tree         # one Newick gene tree per line
outgroup: Out1,Out2                      # comma-separated labels, or a YAML list like ["Out1", "Out2"]
output_folder: results                    # output directory for all run artifacts
overwrite: true                           # true = reset and reuse output dir; false = append suffix
triplet_filter: null                      # path to filter file or null (all triplets)
species_rename_map: null                  # path to rename map or null
processes: 0                              # 0 = every CPU available to the process
seed: null                                # null = generate a run seed at runtime
alpha_dct: 0.05                           # DCT significance threshold
alpha_ks: 0.05                            # KS significance threshold
alpha_perm: 0.05                          # permutation significance threshold
p_value_correction: bfn                   # no, bfn, holm, fdr_bh, fdr_by
pipeline_mode: efficient                  # efficient = skip tests the cascade cannot consult; detailed = run all tests
consolidation: true                       # true = write introgression maps and consolidation CSV/PNG outputs
bootstrap: true                           # true = run bootstrap resampling
preflight_data_check: false               # true = run only the structural preflight check and exit
preflight_triplet_cap: 15000              # 0 = check all triplets
discordant_test: chi-square               # chi-square or z
tree_height_calculation_strategy: AVG     # AVG, A, B, C, SIS, INT
min_support_value: 0.5                    # minimum support threshold for a topology to count
generate_summary_stats: false             # true = also write summary_statistics.tsv
shape_diagnostics: false                  # true = include shape-diagnostics columns
bootstrap_options:
  iterations: 100                         # bootstrap resamples per triplet
  debug_mode: false                       # true = emit debug summaries during bootstrap
  summary_only: false                     # true = keep only the summary output
permutation_options:
  min_resamples: 2500                     # minimum permutation resamples to attempt
  max_resamples: 25000                    # hard upper bound for permutation resamples
  ci_method: wilson                       # CI method for the permutation interval
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

  The display names appear in `orchestrator_triplet_results.tsv`,
  `summary_statistics.tsv`, and the `consolidation/` TSV matrices and plot.
  The outgroup, a triplet filter, the processed tree files and `metrics.txt`
  use the labels as they appear in the input trees. Taxa absent from the map
  keep their tree labels.

  The file is read when the config resolves and rejected if it is missing or
  malformed, maps a label more than once, maps two labels onto the same
  display name, or gives a display name holding a tab, line break, comma,
  semicolon or equals sign — the characters the results TSV uses as
  delimiters.

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
- Worker process count. `0` uses every CPU the process is allowed to run on — which under a container limit or a job scheduler can be fewer than the machine has — and `1` runs serially in the parent process. The count actually used is reported in `metrics.txt` as `Worker processes`, with the CPUs the process could see.
- A worker's memory is one batch of observations at a time plus its share of the read-only geometry cache, which is inherited copy-on-write rather than copied. Starting more workers than there are CPUs gains nothing and multiplies that per-worker memory, which is why `0` stops at the CPUs the process can use.
- The workers are processes on the machine the orchestrator was started on, and that is the only machine a run uses. On a cluster, a job spanning several nodes runs the orchestrator on one of them and leaves the rest idle; request the CPUs and memory on a single node instead. No multi-node execution is implemented.

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
- Multiple-testing correction, applied in three places. Run-wide, once every triplet has been measured, it adjusts the DCT p-values as one family over all triplets and the KS p-values as another — every triplet is a member of both, whatever its other test said — and only then is any triplet classified. Inside each permutation test, it adjusts that test's pair of one-tailed p-values against each other — never across triplets, because a Monte Carlo p-value has a resolution floor that across-triplet correction would fall through. Inside the bootstrap, each iteration's DCT and KS p-values are corrected across triplets for that iteration index, so iterations answer to the same thresholds the reported classification does. Corrected p-values drive the significance decisions; the uncorrected values are retained in the output for reporting. Under `no` the results TSV carries the raw p-values and the significance flags only.
- The choice affects run time and memory. `no` and `bfn` depend only on the family size, which is the triplet count, so each bootstrap iteration votes as it runs and only the tally is kept, and under [`pipeline_mode: efficient`](#pipeline_mode) a tree-height test below a settled count gate can be left unmeasured. The rank-based methods must hold every triplet's per-iteration p-values until the stream finishes, and measure the tree-height test for every triplet and in every bootstrap iteration, because every member's value moves the others' ranks.
- Every supported method is monotone — none can adjust a p-value *below* its raw value — which is what lets a direction test be skipped once an earlier gate has failed, in the bootstrap under every mode and in the point estimate under `efficient`: the skipped result could never have been read.
- In YAML, `p_value_correction: no` may be written with or without quotes. YAML resolves the bare word `no` to a boolean, and enumerated fields map booleans back to the choice they spell (`no`/`off`/`n`/`false`, `yes`/`on`/`y`/`true`), so both forms select the same value.

##### `pipeline_mode`

- CLI: `--pipeline-mode`
- Default: `efficient`
- Allowed: `efficient`, `detailed`
- `efficient` stops measuring a triplet once the decision cascade has settled it, rather than computing a result nothing reads. The permutation direction test is skipped for every triplet that one of the first two gates settled, under every correction method. The tree-height test is skipped below a settled count gate only under `no` and `bfn`; under `holm`, `fdr_bh` and `fdr_by` it is measured for every triplet, because those corrections rank every triplet's value against the others' and leaving one out would change the rest. `detailed` runs all three tests for every triplet.
- The bootstrap follows the same rule. Under `efficient` with `no`/`bfn` an iteration whose count gate failed leaves its tree-height test unmeasured; under `detailed` every iteration measures it, so the bootstrap replicates the procedure the point estimate runs on the real data. An iteration's direction test is skipped below a failed gate in both modes.
- **Both modes produce identical results under every correction method** — not merely the classifications but every value the cascade reads: `classification`, `decision_gate`, the `dct_*` columns, `ks_p_value_corrected`, `ks_significant`, and every bootstrap column match exactly. Two facts make that so. First, every supported correction is monotone, so a gate that failed on the raw p-value cannot clear on the corrected one, and the efficient mode only ever declines a test below a gate that has already failed — under `no`/`bfn` it judges the gate on the exactly corrected value (the family size is the triplet count, known before the run starts), under the rank-based methods on the raw value, which is the conservative side. Second, the permutation p-values are corrected inside each test, across its own pair of one-tailed p-values, never across triplets, so leaving one triplet's direction test unrun changes nothing for any other triplet. The tree-height test is different — its family spans the run — which is exactly why it is skipped only where the family size alone corrects it. The full argument is under "Skipping a settled gate" in the orchestrator guide.
- What differs is which columns are populated, never their values. Triplets the efficient mode settled early leave the `perm_*` block empty and carry `perm_note: direction_test_not_consulted`, which distinguishes a deliberate skip from a test that ran and hit a guard; under `no`/`bfn` a triplet the count gate settled also leaves the raw `ks_statistic`/`ks_p_value` and the corrected `ks_*` columns empty. `metrics.txt` names the mode, says what it skips under the run's correction, and counts the skipped tests.
- Expect `efficient` to be somewhat faster than `detailed`, not dramatically so. What it declines is the point estimate's direction test and, under `no`/`bfn`, the tree-height test, while most of a run's time is the bootstrap's permutation tests: every iteration of a triplet that reaches the third gate runs one, in both modes. The wall time moves less than the CPU time, because the triplets that reach the third gate set the length of the run either way.
- Choose `detailed` when you want the direction test's statistics for every triplet regardless of whether they decided anything, which is a debugging need rather than an analysis one.

##### `consolidation`

- CLI: `--no-consolidation` (sets `consolidation: false`)
- Default: `true`
- Generates the introgression map artifacts into a `consolidation/` subfolder of the output directory.

##### `bootstrap`

- CLI: `--no-bootstrap` (sets `bootstrap: false`)
- Default: `true`
- Enables bootstrap resampling per triplet and adds the `bootstrap_value` and `all_bootstrap` columns to the results TSV. Setting it false skips the iterations entirely, so `bootstrap_perm_stat_ci_low`/`bootstrap_perm_stat_ci_high` are empty too — that interval is a bootstrap percentile interval, not a permutation output. Consolidation then weights every classified triplet as 1 where it would have used `bootstrap_value`, so the introgression maps still build; their cell values become plain counts over the same co-occurrence denominators.
- This is an instruction about what to compute, so it holds under `pipeline_mode: detailed` as well: `detailed` declines to skip work the cascade cannot consult, which is not the same as reinstating work you switched off. The same is true of `generate_summary_stats` and `shape_diagnostics`.

##### `preflight_data_check`

- CLI: `--preflight-data-check` (sets `preflight_data_check: true`)
- Default: `false`
- Runs only the structural sanity check on the species tree, gene trees, and triplets, writes `preflight_data_check.txt` into the output folder, and exits without any analysis. No results TSV, processed trees, `metrics.txt`, or consolidation artifacts are produced. Every other analysis key is ignored for that run.
- The report lists each detected issue by category (for example `gene_tree.rooting_failed`, `triplet.unresolved_rooted_sister_pair`), a count per category, up to 25 example messages naming the offending gene-tree index and triplet, and a species-tree-versus-gene-tree attribution summary.
- The checks are structural: they establish whether the data can be processed, not whether the result will be biologically meaningful.

##### `preflight_triplet_cap`

- CLI: `--preflight-triplet-cap <n>`
- Default: `15000`
- Type: integer `>= 0`; `0` lifts the cap.
- The most ingroup triplets the preflight data check walks. Above it the check takes the first `n` in generation order and reports `analysis.triplet_cap_applied` with the total it skipped, so the report says when it is partial. Parsing the gene trees and building their geometry is a fixed cost, so raising the cap costs less than proportionally: on 1,475 gene trees, 15,000 triplets take about 25 s against 11 s for 1,500.
- Applies only to the check. A real run always processes every triplet, and the key is ignored unless `preflight_data_check` is set. A `triplet_filter` is never capped -- the check walks every triplet the filter names.

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
# required: -st/--species-tree-path, -gt/--gene-trees-path, -og/--outgroup
python -m ghostparser.orchestrator \
  -st data/species.tree \
  -gt data/genes.tree \
  -og Out1,Out2 \
  --output-folder results \
  --processes 0 \
  --alpha-dct 0.05 \
  --alpha-ks 0.05 \
  --alpha-perm 0.05 \
  --p-value-correction bfn \
  --pipeline-mode efficient

# optional CLI switches: -c/--config-file, --triplet-filter, --species-rename-map,
# --seed, --no-overwrite, --no-consolidation, --no-bootstrap,
# --preflight-data-check, --preflight-triplet-cap
# default behaviors: config-file mode ignores other flags; no triplet filter means all triplets;
# overwrite is enabled by default; consolidation and bootstrap are enabled by default;
# preflight checks are disabled by default; 0 checks all triplets.
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
overwrite: true                       # true = clear the existing output directory before writing
target_column: class                 # multi-label target bitstring column
test_size: 0.2                       # fraction reserved for the hold-out split
cv_folds: 5                          # integer >= 1 or null; null disables CV
rare_class_policy: warn_reduce_cv    # warn_reduce_cv, warn_skip_cv, error
random_state: null                   # null = use nondeterministic defaults
n_jobs: -1                           # -1 = all available cores, null = no parallelism

model:
  n_estimators: 200                  # 200 trees by default
  max_depth: 10                      # null = unlimited depth
  min_samples_split: 2               # minimum samples required to split a node
  min_samples_leaf: 1                # minimum samples required in a leaf node
  max_features: null                 # null = use all features; allowed: sqrt, log2, int, float
  class_weight: null                 # null, balanced, balanced_subsample, mapping, or list
  n_neighbors: 5                     # KNN neighbors
  weights: uniform                   # uniform or distance
  algorithm: auto                    # auto, ball_tree, kd_tree, brute
  leaf_size: 30                      # KNN leaf size
  metric: minkowski                  # KNN distance metric
  p: 2                               # Minkowski exponent

evaluation:
  metrics: all                       # all, primary, diagnostic, per_bit
  report_class_distribution: true    # include the dataset summary block
  report_confusion_matrix: true      # include confusion-matrix metrics
  report_feature_importance: true     # include feature importance output
  save_label_map: true               # store the label map in the metrics JSON
  save_predictions: true             # save per-row predictions to TSV
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
  - `max_features` (string, int, float or null, default: `null`) — how many features each split may consider. `null` uses **every** feature at each split; `sqrt` and `log2` take that function of the feature count, an integer `>= 1` that many features, and a float in `(0.0, 1.0]` that fraction of them. `auto` is rejected: scikit-learn removed it in 1.3, and `sqrt` is its classifier equivalent.
  - `class_weight` (string, dict, list or null, default: `null`) — accepts `balanced`, `balanced_subsample`, a mapping of class label to weight, a list of such mappings (one per label), or `null` for no class weighting.

  **Writing null.** For either key, null is the default, so the simplest way to get it is to omit the key. To write it out, use `null` in YAML or JSON (YAML also accepts `~`). Do **not** write the bare words `None` or `none`: YAML reads them as the strings `"None"` and `"none"`, which are not null and are rejected like any other unrecognized value. In a `hyperparameter_tuning.search_space` list the same applies per candidate — `max_features: [sqrt, null]`, not `[sqrt, None]`.

  Any other value for either key is rejected with a message naming what it received, every accepted form, and how to write null.

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
- `use_wandb` (bool, default `false`): whether the run logs to Weights & Biases. Logging is opt-in, so a config that omits the key runs locally; a non-boolean is rejected. With `false` the tuner needs neither a W&B account nor the `wandb` package, makes no network calls, writes no `wandb/` directory, and writes the full local artifact set. With `true` the `wandb` package must be installed (`pip install .[wandb]`) and authenticated, and the bulk outputs (`hyper_tune_results.json`, `hyper_tune_results.tsv`, `hyper_tune_parameter_marginals.tsv`, `predictions.tsv`) are logged to the W&B run instead of the output directory, which then keeps only the model pickle, the plaintext report, and the search-report plot.
- `wandb_detailed_payloads` (bool, default `false`): when `true`, log additional per-candidate CV payloads (aggregate and fold-level JSON) and extra summary JSON blobs to Weights & Biases. Requires `use_wandb: true`; combining it with `use_wandb: false` is rejected. Keep `false` when network/storage overhead matters.
- `search_space` (mapping): model hyperparameter candidates. Each parameter should map to a list of values. Omit a parameter from `search_space` if you want the trainer default to apply during tuning.

Allowed `search_space` keys depend on `model`:

- `random_forest`: `n_estimators`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `max_features`, `class_weight`
- `multi_knn`: `n_neighbors`, `weights`, `algorithm`, `leaf_size`, `metric`, `p`

Do not place runtime fields such as `input_path`, `output_dir`, `target_column`, `test_size`, `cv_folds`, `rare_class_policy`, `random_state`, or `n_jobs` inside `search_space`.

Example:

```yaml
input_path: ./results/summary_statistics.tsv
output_dir: ./results/hyper_tune_out
overwrite: true
target_column: class
test_size: 0.2
cv_folds: 5
rare_class_policy: warn_reduce_cv
random_state: null
n_jobs: -1

hyperparameter_tuning:
  model: random_forest             # random_forest or multi_knn
  method: grid                    # grid or random
  objective: exact_match_accuracy  # exact_match_accuracy, hamming_loss, bitwise_accuracy, micro_f1, macro_f1, weighted_f1
  top_k: 10                       # keep the top K candidates in the summary report
  n_iter: 20                      # sampled candidates when method: random
  max_candidates: 5000            # hard cap for full-grid evaluation
  use_wandb: false                # true = log to Weights & Biases
  wandb_detailed_payloads: false  # true = log extra per-candidate payloads when use_wandb: true
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

The trainer CLIs accept only the options below. Omitted options use the shown
defaults; `-c/--config-file` takes precedence over the other trainer options.

```bash
# Random forest and Multi-KNN trainers use the same four options.
python -m ghostparser.ml.random_forest \
  -i ./results/summary_statistics.tsv \
  -o ./results/ml_out \
  --no-overwrite
python -m ghostparser.ml.multi_knn \
  -i ./results/summary_statistics.tsv \
  -o ./results/ml_out

# optional CLI switches: -c/--config-file, -i/--input-path, -o/--output-dir,
# --no-overwrite
# default behaviors: config-file mode ignores other flags; input/output are required unless
# supplied by config; overwrite is enabled by default.

# The hyperparameter tuner accepts these two options:
python -m ghostparser.ml.hyper_tune \
  -c ./sample_configs/hyperparameter_tuning_random_forest.yaml
# optional CLI switches: -c/--config-file, --no-overwrite
# default behaviors: config-file mode is required; overwrite is enabled by default.
```

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
