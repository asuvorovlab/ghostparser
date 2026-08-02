
# Machine Learning Module

This document describes the machine-learning baselines included with Ghostparser under `ghostparser.ml`. It focuses on what the configuration keys control, what they change for a training run, and how the evaluation metrics behave. The random forest baseline is the primary focus; multi-label KNN remains available as a secondary baseline.

Install the optional ML dependency set with `pip install .[ml]` if you want to run these baselines. The ML extra includes `scikit-learn` and `wandb`.

## Contents
- [Overview](#overview)
- [Data contract](#data-contract)
- [Usage](#usage)
- [Configuration](#configuration)
- [Evaluation config keys](#evaluation-config-keys)
- [Evaluation argument](#evaluation-argument)
- [Random Forest settings](#random-forest-settings)
- [Multi-label KNN settings](#multi-label-knn-settings)
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

Note: there is no package-level model dispatcher here; call `python -m ghostparser.ml.random_forest` or `python -m ghostparser.ml.multi_knn` directly.

To run the KNN baseline explicitly:

```bash
python -m ghostparser.ml.multi_knn -i /path/to/summary_statistics.tsv -o ml_results/
```

To run the hyperparameter tuner explicitly:

```bash
python -m ghostparser.ml.hyper_tune -c sample_configs/hyperparameter_tuning_random_forest.yaml
```

## Configuration

The ML loaders enforce a strict nested layout. For the full schema and examples, see [CONFIG.md](../../CONFIG.md#machine-learning-ghostparserml), but the most important user-facing idea is simple: the top level tells the trainer where the data is and how to split it, `model` controls the estimator itself, and `evaluation` controls what you want reported or saved. The loader handles numeric features directly and one-hot encodes low-cardinality string features for you.

### Core run controls

- `input_path`: path to the TSV input file.
- `output_dir`: directory where artifacts are written.
- `overwrite`: when `true`, existing output directories are cleared before trainer artifacts are written; when `false`, a suffix such as `_1` is appended.
- `target_column`: the label column containing the fixed-length binary target string, which defaults to `class`.
- `test_size`: fraction reserved for hold-out evaluation.
- `cv_folds`: how many cross-validation folds to attempt.
- `random_state`: seed for reproducible splits and model randomness.
- `n_jobs`: CPU parallelism control; `-1` uses all available CPU cores for operations that support parallelism, and it does not use the GPU.
- `rare_class_policy`: what to do if a requested CV split is not feasible because some label combinations are too rare.

### Evaluation argument

Use `evaluation.metrics: all` if you want both per-label metrics and the strict exact-match view. That setting keeps the label-level scores, the exact-match accuracy, and the per-bit breakdown together, which is the most useful mode when you want to know both “how many labels were right?” and “did the whole 6-bit pattern match exactly?”

The other options are narrower slices: `primary` keeps the core label-level metrics, `diagnostic` adds the strict exact-match and classification report, and `per_bit` adds the label-by-label breakdown. If your goal is one mistaken bit should still count as mostly correct, use the per-label metrics; if your goal is the full 6-bit pattern must match exactly, use the exact-match metric alongside the per-label ones.

### Evaluation config keys

The trainers accept an `evaluation` mapping in the config that controls which metrics are computed and which artifacts are written to disk. Keys and effects:

- `metrics` (string, default: `all`): one of `all`, `primary`, `diagnostic`, `per_bit`.
  - `all`: compute and include per-bit metrics, aggregate/per-label metrics (micro/macro/weighted F1, precision/recall), and diagnostic outputs (exact-match accuracy, classification report).
  - `primary`: compute the core aggregate/per-label metrics and hamming loss (a compact view for model comparison).
  - `diagnostic`: compute strict diagnostics such as exact-match accuracy, the full classification report, and confusion matrices (confusion matrices are produced only when `report_confusion_matrix` is enabled).
  - `per_bit`: compute and return the per-label breakdown for each bit separately.

- `report_class_distribution` (bool, default: `true`): when enabled, the text report repeats the dataset summary. The overall metrics JSON always contains the label map plus the class and bit distributions for the overall file and the train/test split.

- `report_confusion_matrix` (bool, default: `true`): when enabled and `metrics` includes diagnostic outputs (`diagnostic` or `all`), compute confusion matrices for each bit and include them in the metrics payload and as a persisted artifact. The trainers also render a single heatmap-style PNG with one subplot per bit so the true/false and predicted 0/1 counts are easy to compare visually.

- `report_feature_importance` (bool, default: `true`): when enabled, compute and persist `feature_importances.tsv`. For tree-based models this is fast (built-in feature importances); for non-tree models (KNN) this uses permutation importance and can be slow.

- `save_label_map` (bool, default: `true`): the label map is embedded in the overall metrics JSON so downstream parsing can read it alongside the split summary.

- `save_predictions` (bool, default: `true`): persist `predictions.tsv` containing per-row true/predicted bit flags, an exact-match indicator, and the number of matched bits.

These keys let you trade computation and storage cost for diagnostic depth: enabling `report_feature_importance` and the `diagnostic` metrics gives the richest outputs but increases runtime and disk usage.

### Hyperparameter tuning

The `hyperparameter_tuning` section configures the standalone tuner in `ghostparser.ml.hyper_tune`. It is separate from `model` and `evaluation` so the search strategy stays explicit and easy to read. For the tuner, only the runtime keys (`input_path`, `output_dir`, `overwrite`, `target_column`, `test_size`, `cv_folds`, `rare_class_policy`, `random_state`, `n_jobs`) plus `hyperparameter_tuning` are allowed at the top level. Do not provide `evaluation` or `model` sections in a tuning config; the tuner does not read them.

Suggested keys:

- `model` (string, default: `random_forest`): which ML module to tune. Choices are `random_forest` and `multi_knn`.
- `method` (string, default: `grid`): `grid` evaluates every combination in the search space, while `random` samples `n_iter` combinations from that space.
- `objective` (string, default: `exact_match_accuracy`): metric used to rank candidates. Choices match the trainer metrics such as `exact_match_accuracy`, `hamming_loss`, `bitwise_accuracy`, `micro_f1`, `macro_f1`, and `weighted_f1`.
- `top_k` (int, default: `10`): how many of the best candidates to keep in the human-readable summary.
- `n_iter` (int, default: `20`): how many candidates to sample when `method: random` is selected.
- `max_candidates` (int, default: `5000`): safety limit for `method: grid`; if the full grid would exceed this value, the tuner rejects the config.
- `wandb_detailed_payloads` (bool, default: `false`): when `true`, send additional detailed candidate payloads (aggregate and fold-level CV JSON) to Weights & Biases; use this only when network/storage overhead is acceptable.
- `search_space` (mapping): parameter grid for the selected model. Each key should be one supported hyperparameter and each value should be a list of candidate values.
- `search_space` (mapping): parameter grid for the selected model. Each key should be one supported hyperparameter and each value should be a list of candidate values. If you omit a parameter from `search_space`, the tuner uses the trainer default for that parameter.

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

#### Weights & Biases setup (required for tuner)

`ghostparser.ml.hyper_tune` initializes Weights & Biases for every tuning run.

Initial setup:

1. Install ML dependencies:

```bash
pip install .[ml]
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

What the tuner logs to WandB:

- run metadata (model, method, objective, candidate count, CV folds)
- one lightweight record per candidate (score, params, elapsed time, best-so-far flag)
- final metrics and timing summaries

If you set `hyperparameter_tuning.wandb_detailed_payloads: true`, the tuner also logs additional JSON payloads per candidate (CV aggregate and fold-level details) and richer run-summary JSON fields.

To keep network and memory overhead low on long runs, the integration logs scalar summaries only (no per-fold raw prediction payload uploads and no large artifact uploads to WandB by default).

### Random Forest settings

- `n_estimators`: more trees usually make the model steadier, but training takes longer and uses more memory.
- `max_depth`: limits how deep each tree can grow; smaller values usually reduce overfitting and make the model faster. A `null` value leaves tree depth unconstrained, so each tree can expand until the split rules stop it.
- `min_samples_split`: requires more samples before a node can split, which makes the trees less sensitive to noise. Larger values make the trees more conservative, which can help when the training set is smaller, noisier, or the model is overfitting.
- `min_samples_leaf`: forces each leaf to contain more samples, which smooths predictions and can help generalization. Larger values make the model smoother and less sensitive to tiny, unstable groups.
- `max_features`: controls how many features each split considers; smaller values increase tree diversity, larger values make each tree more greedy. A `null` value means the estimator uses its default feature-selection behavior for each split.
- `class_weight`: lets you weight label classes differently if you want to favor rare outcomes. A `null` value means no class weighting is applied, so all classes are treated equally.

For `min_samples_split` and `min_samples_leaf`, leave the keys out if you want the configured defaults. The model does not infer these values from the data, and `null` is rejected by the loader.

If your TSV contains string-valued columns, the loader will one-hot encode them only when they have 7 or fewer distinct values.

### Multi-label KNN settings

- `n_neighbors`: controls how many nearby training examples vote for a prediction; smaller values are more local, larger values are smoother.
- `weights`: `uniform` treats all neighbors equally, while `distance` gives closer neighbors more influence.
- `algorithm`: chooses the neighbor-search strategy; this mostly changes performance, not the final meaning of the model.
- `leaf_size`: tuning knob for tree-based neighbor search performance and memory usage.
- `metric`: distance function used to compare samples.
- `p`: the Minkowski distance power, where `1` behaves like Manhattan distance and `2` behaves like Euclidean distance.

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
- `*_metrics.txt` — Human-readable summary of metrics, dataset summary, and diagnostic notes
- `*_confusion_matrices.png` — Heatmap grid of all confusion matrices, colored from red (smaller counts) to green (larger counts)
- `*_confusion_matrix_64_classes.png` — Heatmap of the full 64-class confusion matrix across all possible 6-bit labels
- `feature_importances.tsv` — Ranked features and importance scores (tree-based for RF, permutation for KNN)
- `predictions.tsv` — Per-row predictions with true/pred bit flags, exact-match indicator, and matched-bit count

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

