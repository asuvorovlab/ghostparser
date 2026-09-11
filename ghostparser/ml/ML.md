
# Machine Learning Module

This document describes the machine-learning baselines included with Ghostparser under `ghostparser.ml`. It focuses on how the trainers work — the data contract, the training and evaluation flow, and how to read the outputs. The random forest baseline is the primary focus; multi-label KNN remains available as a secondary baseline. The full configuration reference lives in [CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml).

Install the optional ML dependency set with `pip install .[ml]` if you want to run these baselines. The ML extra includes `scikit-learn`. Weights & Biases is a separate extra (`pip install .[wandb]`) that only `hyper_tune` can use, and only when a run asks for it — see [Weights & Biases logging](#weights--biases-logging-optional).

## Contents
- [Overview](#overview)
- [Data contract](#data-contract)
- [Usage](#usage)
- [Configuration](#configuration)
- [Choosing an evaluation mode](#choosing-an-evaluation-mode)
- [Hyperparameter tuning](#hyperparameter-tuning)
- [Reading a tuning run](#reading-a-tuning-run)
- [Weights & Biases logging (optional)](#weights--biases-logging-optional)
- [Tuning the estimators](#tuning-the-estimators)
- [Feature Importance](#feature-importance)
- [What the model outputs mean](#what-the-model-outputs-mean)
- [Cross-validation and reproducibility](#cross-validation-and-reproducibility)
- [Outputs](#outputs)
- [Notes and diagnostics](#notes-and-diagnostics)
- [See also](#see-also)

## Overview

Purpose: provide a reproducible multi-label baseline that uses summary statistics from `summary_statistics.tsv` to predict a fixed-length binary label string.

The current implementation expects a 6-bit target string and treats it as six separate binary labels, so a prediction like `001010` versus `001011` is only wrong on the last label under per-label metrics. If you want the strict “all bits must match” view, that is reported as exact-match accuracy. The per-label and exact-match ideas are general, but the current code is hard-wired to a 6-label target and would need code changes to support a different label count.

Models provided:
- `random_forest` — RandomForestClassifier wrapped in MultiOutputClassifier
- `multi_knn` — KNeighborsClassifier wrapped in MultiOutputClassifier

## Data contract

- `class`: canonical 6-bit bitstring. Each bit corresponds to the labels in the following order:
  1. ghost_into_A
  2. ghost_into_B
  3. inflow_into_A_from_C
  4. inflow_into_B_from_C
  5. outflow_from_A_to_C
  6. outflow_from_B_to_C
- Feature columns: provide the columns you want the model to use, plus the target column. Numeric columns are used directly. String-valued feature columns with 7 or fewer distinct values are one-hot encoded automatically. Columns with more than 7 distinct values are rejected so the model does not invent an arbitrary ordering.

The `class` column is stored as a string bitstring in the TSV and that is expected. The loader reads it as a label string, validates that it is a 6-character `0`/`1` pattern, and then expands it into six binary targets for training. That does not cause a problem.

## Usage

Command line (quick):

Run the Random Forest trainer explicitly as the module:

```bash
python -m ghostparser.ml.random_forest -i /path/to/summary_statistics.tsv -o ml_results/
```

There is no package-level model dispatcher. Call `python -m ghostparser.ml.random_forest` or `python -m ghostparser.ml.multi_knn` directly.

To run the KNN baseline explicitly:

```bash
python -m ghostparser.ml.multi_knn -i /path/to/summary_statistics.tsv -o ml_results/
```

To run the hyperparameter tuner explicitly:

```bash
python -m ghostparser.ml.hyper_tune -c sample_configs/hyperparameter_tuning_random_forest.yaml
```

The tuner is config-file only. Two runnable samples ship with the package, one per
format and per model, and both set every required key:

- `sample_configs/hyperparameter_tuning_random_forest.yaml` — random-forest grid search
- `sample_configs/hyperparameter_tuning_multi_knn.json` — multi-label KNN random search

Both set `use_wandb: false`, so they run without a Weights & Biases account.

## Configuration

The ML loaders enforce a strict nested layout. For the full schema and examples, see [CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml), but the most important user-facing idea is simple: the top level tells the trainer where the data is and how to split it, `model` controls the estimator itself, and `evaluation` controls what you want reported or saved. The loader handles numeric features directly and one-hot encodes low-cardinality string features for you.

The top-level keys point the trainer at the data and control the split
(`input_path`, `output_dir`, `overwrite`, `target_column`, `test_size`,
`cv_folds`, `random_state`, `n_jobs`, `rare_class_policy`); `model` holds the
estimator hyperparameters; `evaluation` selects what gets reported and saved.
Every key, its default, and its allowed values are documented in
[CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml).

### Choosing an evaluation mode

`evaluation.metrics` decides which view of correctness you get, and the right
choice depends on how you want a near-miss to count:

- `all` (the default) keeps everything together — the label-level scores, the
  exact-match accuracy, and the per-bit breakdown. Use it when you want to know
  both "how many labels were right?" and "did the whole 6-bit pattern match?"
- `primary` keeps the core aggregate/per-label metrics and hamming loss: a
  compact view for comparing models.
- `diagnostic` adds the strict exact-match accuracy, the full classification
  report, and (when `report_confusion_matrix` is on) per-bit confusion matrices.
- `per_bit` returns the label-by-label breakdown only.

If one mistaken bit should still count as mostly correct, read the per-label
metrics; if the full 6-bit pattern must match, read exact-match accuracy
alongside them.

The `report_*` and `save_*` keys trade runtime and disk for diagnostic depth.
`report_feature_importance` is cheap for the random forest (built-in
importances) but slow for KNN (permutation importance). The label map and the
class/bit distributions are always embedded in the overall metrics JSON.

### Hyperparameter tuning

The `hyperparameter_tuning` section configures the standalone tuner in `ghostparser.ml.hyper_tune`. It is separate from `model` and `evaluation` so the search strategy stays explicit and easy to read. For the tuner, only the runtime keys (`input_path`, `output_dir`, `overwrite`, `target_column`, `test_size`, `cv_folds`, `rare_class_policy`, `random_state`, `n_jobs`) plus `hyperparameter_tuning` are allowed at the top level. Do not provide `evaluation` or `model` sections in a tuning config; the tuner does not read them.

`method: grid` evaluates every combination in the search space and is rejected
if the full grid would exceed `max_candidates`; `method: random` samples
`n_iter` combinations instead, which is what you want once the space is large.
Candidates are ranked by `objective` (any trainer metric, such as
`exact_match_accuracy` or `micro_f1`), and the best `top_k` are kept in the
human-readable summary. Any parameter omitted from `search_space` falls back to
the trainer default. See
[CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml) for the full key
reference.

Supported `search_space` keys are:

- `random_forest`: `n_estimators`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `max_features`, `class_weight`
- `multi_knn`: `n_neighbors`, `weights`, `algorithm`, `leaf_size`, `metric`, `p`

Do not put runtime keys such as `input_path`, `output_dir`, `target_column`, `test_size`, `cv_folds`, `rare_class_policy`, `random_state`, or `n_jobs` inside `search_space`.

Example:

```yaml
hyperparameter_tuning:
  model: random_forest
  method: grid
  objective: exact_match_accuracy
  top_k: 5
  max_candidates: 1000
  search_space:
    n_estimators: [100, 200, 400]
    max_depth: [null, 10, 20]
    min_samples_split: [2, 5]
    min_samples_leaf: [1, 2]
    max_features: [sqrt, log2]
    class_weight: [null]
```

Use `method: random` when the space is large and you want a sampled search instead of checking every combination.

The tuner prints console progress while it runs, including the number of candidate cases it plans to evaluate, the approximate number of model fits implied by CV, and per-candidate timing updates.

### Reading a tuning run

Every tuning run with `use_wandb: false` writes a self-contained set of local
artifacts, so a search can be navigated end to end without any external service. With
`use_wandb: true` the bulk outputs move into the Weights & Biases run and the output
directory keeps only `hyper_tune_best_model.pkl`, `hyper_tune_results.txt`, and
`hyper_tune_search_report.png` — see [Weights & Biases logging](#weights--biases-logging-optional).

`hyper_tune_results.txt` is the place to start. It opens with the run header (model,
method, objective and its direction, candidate count, CV folds, whether W&B logging
was on), echoes the `search_space` that was actually explored, then reports:

- **Top _k_ candidates** — an aligned table of the best candidates, one row each, with
  the CV score, its standard deviation across folds, the wall-clock seconds the
  candidate took, and every search-space parameter as its own column. `top_k` sets the
  row count.
- **Per-parameter value summary** — one row per (parameter, value) pair: how many
  candidates used the value, the best/mean/std/worst score it reached, the best rank
  it achieved, and its rank among the values of the same parameter. This is the marginal
  view: it answers "was `n_estimators: 400` ever worth it?" without re-reading the full
  candidate list.
- **Search-space guidance** — the observed score range; a parameter-influence ordering
  (how far the objective moved across each parameter's values, by both best and mean
  score, so a dimension that changed nothing is named outright); the best value for each
  parameter; and a warning for any winning value that sits at the edge of the range that
  was searched, which is the signal to extend the range in that direction and search again.
- **Timings** and the **artifact paths** for the rest of the run.

The same content is available in machine-readable form, written to the output
directory when `use_wandb: false` and logged to the W&B run when `use_wandb: true`:

- `hyper_tune_results.tsv` — every candidate, ranked best-first, with `rank`, `is_best`,
  `cv_score`, the objective's mean and standard deviation, `elapsed_seconds`, and one
  column per search-space parameter.
- `hyper_tune_parameter_marginals.tsv` — the per-parameter value summary, plus an
  `at_search_bound` column marking a winning value that is the smallest or largest
  numeric value tried.
- `hyper_tune_results.json` — the full run payload, including `parameter_marginals`,
  `parameter_influence`, `artifact_paths`, and `use_wandb`.
- `predictions.tsv` — one row per test sample, with the true and predicted bits.

`hyper_tune_search_report.png` is the visual counterpart, stacked in one figure:

1. **Search progress in evaluation order** — every candidate's score against the order
   it was evaluated, with the running-best trace over it. On a random search this shows
   whether the search had converged or was still improving when it stopped.
2. **Top candidates** — a ranked dot plot of the best candidates, each labelled with its
   parameter combination and score. Objective scores sit in a narrow band away from zero,
   so this is a dot plot on a zoomed axis rather than bars from a zero baseline.
3. **One panel per parameter that actually varied** — the score distribution for each
   value, box plus the individual candidates. A panel with clearly separated levels is a
   parameter that matters; a flat panel is a dimension you can drop from the next search.

### Weights & Biases logging (optional)

`hyperparameter_tuning.use_wandb` is an optional boolean defaulting to `false`, so a
config that never mentions it logs nothing to Weights & Biases. Runs with
`use_wandb: false` need neither a W&B account nor the `wandb` package installed, and
they still produce the full local report described above.

With `use_wandb: false` the tuner writes no `wandb/` directory and makes no network
calls; the report header records `Weights & Biases logging: disabled` and
`hyper_tune_results.json` carries `"use_wandb": false`.

`use_wandb: true` also changes where the bulk outputs land. The W&B run becomes their
home, so the tuner does not duplicate them on disk:

| Output | `use_wandb: false` | `use_wandb: true` |
| --- | --- | --- |
| `hyper_tune_best_model.pkl` | output directory | output directory |
| `hyper_tune_results.txt` | output directory | output directory |
| `hyper_tune_search_report.png` | output directory | output directory |
| `hyper_tune_results.tsv` | output directory | W&B table `tables/ranked_candidates` |
| `hyper_tune_parameter_marginals.tsv` | output directory | W&B table `tables/parameter_marginals` |
| `predictions.tsv` | output directory | W&B table `tables/predictions` |
| `hyper_tune_results.json` | output directory | W&B run summary key `results_json` |

The plaintext report still carries the top candidates, the per-parameter value summary,
and the search-space guidance in full, so the output directory stays readable on its own.
Its `Artifacts:` block lists only the files that were actually written, and names W&B as
the home of the rest. The `artifact_paths` key of the run payload follows the same rule,
and the run summary records `local_artifact_paths` and `bulk_artifacts_written_locally`.

To use `use_wandb: true`:

1. Install the extras:

```bash
pip install .[ml,wandb]
```

2. Authenticate once:

```bash
wandb login
```

3. (Optional) set project/account defaults:

```bash
export WANDB_PROJECT=ghostparser-hyper-tune
export WANDB_ENTITY=<your-wandb-entity>
```

4. (Optional) use offline mode when needed:

```bash
export WANDB_MODE=offline
```

What the tuner logs to W&B:

- run metadata (model, method, objective, candidate count, CV folds)
- one lightweight record per candidate (score, params, elapsed time, best-so-far flag)
- final metrics and timing summaries
- the ranked candidates, parameter marginals, and per-row predictions as W&B tables,
  plus the full results payload as the `results_json` run-summary field

`hyperparameter_tuning.wandb_detailed_payloads: true` adds per-candidate CV payloads
(aggregate and fold-level JSON) and richer run-summary JSON fields. It requires
`use_wandb: true`; setting it alongside `use_wandb: false` is rejected by the config
loader.

To keep network and memory overhead low on long runs, the integration logs scalar
summaries only (no per-fold raw prediction payload uploads and no large artifact uploads
to W&B by default).

### Tuning the estimators

Both estimators are wrapped in a `MultiOutputClassifier`, so each hyperparameter
applies to all six one-vs-rest models. The keys and their defaults are listed in
[CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml); the practical
intuition is:

- **Random forest** — `n_estimators` trades stability for runtime and memory.
  `max_depth`, `min_samples_split`, and `min_samples_leaf` all constrain how far
  the trees can chase noise, so raising them makes the model more conservative
  (useful on small or noisy training sets). `max_features` controls per-split
  feature sampling: smaller values make the trees more diverse, larger values
  make each tree greedier, and `null` lets every split see every feature. `class_weight` can favour rare outcomes.
- **Multi-label KNN** — `n_neighbors` sets how local a prediction is, and
  `weights: distance` lets closer neighbours dominate. `algorithm` and
  `leaf_size` affect search performance rather than the model's meaning, while
  `metric` and `p` define the distance function (`p: 1` is Manhattan, `p: 2` is
  Euclidean).

`min_samples_split` and `min_samples_leaf` are not inferred from the data and
the loader rejects `null` for them — omit the keys to accept the defaults.
Note that KNN caps `n_neighbors` at the training-set size when the configured
value would exceed it.

### Feature Importance

The feature-importance score depends on the trainer:

- Random Forest averages the fitted `feature_importances_` values from the six one-vs-rest estimators. Those values are mean decrease in impurity scores, so larger values mean the feature helped the trees split the data more effectively. We compute `feature_importances_` on each of the six binary estimators and average them to produce a single importance per input feature. The underlying score is the mean decrease in impurity (MDI); it is fast but can be biased toward features with many possible split points and can be difficult to interpret when features are strongly correlated.
- Multi-label KNN uses permutation importance on the held-out test split with `scoring="f1_micro"`. Each feature is shuffled repeatedly, and the score measures how much the micro-F1 drops on average. Larger values mean the model relied on that feature more strongly. We use permutation importance on the held-out test split (the implementation uses `n_repeats=10` by default). Each feature is shuffled `n_repeats` times and we report the mean drop in `f1_micro`; this directly measures the impact on the chosen evaluation metric but is slower and can show negative values when shuffling by chance improves the metric on small test sets.

These are model-relative importance measures, not causal explanations. Scores near zero mean the feature had little effect under the chosen trainer and evaluation metric; treat them as heuristic indicators of influence rather than proofs of causality.

### What the model outputs mean

Both baselines convert the 6-bit target string into a multi-label problem with six binary outputs. That means the model can be “partly right” on a row: if only one bit is wrong, the per-label metrics reflect that single label error, while exact-match accuracy marks the entire 6-bit prediction as incorrect.

If you change the number of labels in the future, the inference interpretation changes only in the obvious way: per-label metrics still measure each label independently, but exact-match accuracy becomes stricter as the label count increases. The current code does not infer a new label count automatically, so changing that target shape requires code changes rather than a config tweak.

The sample config files in [random_forest_minimal.yaml](https://github.com/asif256000/ghostparser/blob/main/sample_configs/random_forest_minimal.yaml) and [multi_knn_minimal.yaml](https://github.com/asif256000/ghostparser/blob/main/sample_configs/multi_knn_minimal.yaml) are good references for the respective trainers.


## Cross-validation and reproducibility

- When `cv_folds` is set, the trainer attempts stratified-like splitting based on label combinations. If exact stratification is impossible due to rare labels, `rare_class_policy` controls fold reduction or CV skipping.
- All splitting uses `random_state` for determinism.
- `n_jobs` only affects CPU-side parallel work in scikit-learn; it does not change the model into a GPU-backed implementation.

Stratified here means we try to preserve the frequency of each 6-bit label combination across folds so each fold has a similar class distribution. If a particular combination is too rare to appear in every fold, the trainer follows `rare_class_policy` (reduce folds or skip CV) to avoid invalid splits.


## Outputs

- `*_model.pkl` — Pickled trained `MultiOutputClassifier` (prefix: `random_forest_` or `multi_knn_`)
- `*_overall_metrics.json` — Structured metrics, timings, and the consolidated dataset summary
- `*_metrics.txt` — Human-readable summary. Opens with a `Hyperparameters:` block listing every knob that shaped the run — the estimator settings, `test_size`, the requested and effective `cv_folds`, `rare_class_policy`, and `target_column` — followed by the metrics, dataset summary, and diagnostic notes. The same mapping is available under the `hyperparameters` key of `*_overall_metrics.json`.
- `*_confusion_matrices.png` — Heatmap grid of the six per-bit confusion matrices on the `cividis` scale, dark for smaller counts and bright for larger ones. Each panel is titled with its bit in prose (`Inflow into A from C`), keeping the taxon letters upper-cased
- `*_confusion_matrix_64_classes.png` — Heatmap of the full 64-class confusion matrix across all possible 6-bit labels, on the same colormap. Cells are row-normalized: each one is the fraction of its true class that was predicted into that column, so the scale is a fixed 0-1 regardless of how many test rows a class drew, every populated row sums to 1, and the diagonal reads as per-class recall. Cells at 0 are painted at the low end of the scale like any other value, so the grid reads as one continuous surface and the colorbar covers every cell in it. Both axes label all 64 classes — `Predicted 6 binary class` across and `True 6 binary class` down — so a specific class can be looked up directly rather than counted inwards from a sparser tick. The classes run in the same order on both axes, by the number of set bits: `000000` first, then the six single-bit classes, up to `111111` last, lexically within a count. The `confusion_matrix_64_classes` entry in the metrics JSON carries its `class_labels` in that same order, so its `matrix` rows and columns read positionally against them.
- `feature_importances.tsv` — Ranked features and importance scores (tree-based for RF, permutation for KNN)
- `predictions.tsv` — Per-row predictions with true/pred bit flags, exact-match indicator, and matched-bit count

`hyper_tune` writes its own set instead — `hyper_tune_best_model.pkl`,
`hyper_tune_results.json`, `hyper_tune_results.txt`, `hyper_tune_results.tsv`,
`hyper_tune_parameter_marginals.tsv`, `hyper_tune_search_report.png`, and
`predictions.tsv` — all described in [Reading a tuning run](#reading-a-tuning-run).
With `use_wandb: true` the bulk members of that set are logged to the W&B run rather
than written to the output directory.

### Using the model pickle

The model pickle stores the fitted `MultiOutputClassifier`, so you can reload it later and score new data with the same feature layout.

```python
import pickle

from ghostparser.ml import ml_utils

with open("results/random_forest_model.pkl", "rb") as handle:
  model = pickle.load(handle)

rows = ml_utils.read_tsv_rows("new_summary_statistics.tsv")
matrix = ml_utils.rows_to_matrix(rows, target_column="class")
predicted_bits = model.predict(matrix.train_features)
predicted_labels = ["".join(str(int(bit)) for bit in row) for row in predicted_bits]
```

The loaded model expects the encoded feature matrix that Ghostparser builds during training, so the input data should use the same columns and string-encoding rules as the training TSV.


## Notes and diagnostics

- The ML modules prioritise multi-label measures; exact-match accuracy is a strict diagnostic useful for end-to-end checks but not optimised as a primary loss.
- The `evaluation` section in the config lets you toggle both computation and persistence of diagnostic artifacts; disabling an artifact saves time and disk space for large runs.

## See also

- Configuration: [CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml)

