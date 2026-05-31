
# Machine Learning Module

This document describes the machine-learning baselines included with Ghostparser under `ghostparser.ml`. It focuses on what the configuration keys control, how each trainer is implemented, and where to find more exhaustive configuration examples (see [CONFIG.md](CONFIG.md)).

Install the optional ML dependency set with `pip install .[ml]` if you want to run these baselines.

## Contents
- [Overview](#overview)
- [Data contract](#data-contract)
- [Usage](#usage)
- [Configuration](#configuration)
- [Effects of configuration keys](#effects-of-configuration-keys)
- [Implementation](#implementation)
  - [Random Forest](#random-forest)
  - [Multi-label KNN](#multi-label-knn)
- [Cross-validation and reproducibility](#cross-validation-and-reproducibility)
- [Outputs](#outputs)
- [Notes and diagnostics](#notes-and-diagnostics)
- [See also](#see-also)

## Overview

Purpose: provide a reproducible multi-label baseline that uses summary statistics from `summary_statistics.tsv` to predict a 6-bit `classes` bitstring.

Models provided:
- `random_forest` — RandomForestClassifier wrapped in MultiOutputClassifier
- `multi_knn` — KNeighborsClassifier wrapped in MultiOutputClassifier

## Data contract

- `classes`: canonical 6-bit bitstring. Each bit corresponds to the labels in the following order:
  1. ghost_into_A
  2. ghost_into_B
  3. inflow_into_A_from_C
  4. inflow_into_B_from_C
  5. outflow_from_A_to_C
  6. outflow_from_B_to_C
- Feature columns: all non-metadata columns in the TSV. The trainers exclude these metadata columns: `triplet`, `abc_mapping`, `species_tree`, `classification`, `bootstrap_value`, `source_folder`, and `classes`. All feature values must be numeric.

## Usage

Command line (quick):

Run the Random Forest trainer explicitly as the module:

```bash
python -m ghostparser.ml.random_forest -i /path/to/summary_statistics.tsv -o ml_results/
```

Or dispatch via the package entrypoint (explicit model required):

```bash
python -m ghostparser.ml --model random_forest -i /path/to/summary_statistics.tsv -o ml_results/
```

Note: `python -m ghostparser.ml` without `--model` does not redirect — you must explicitly request the model to run.

To run the KNN baseline explicitly:

```bash
python -m ghostparser.ml.multi_knn -i /path/to/summary_statistics.tsv -o ml_results/
```

## Configuration

The ML loaders enforce a strict, nested configuration layout. For examples and full schema notes, see [CONFIG.md](CONFIG.md#machine-learning-ghostparserml).

Top-level runtime keys (always required or recommended):
- `input_path` (string): path to the TSV input file (required)
- `output_dir` (string): directory where artifacts are written (required)
- `target_column` (string, default: `classes`): column containing the 6-bit labels
- `test_size` (float, 0-1): holdout fraction used for the final test split
- `cv_folds` (int): number of cross-validation folds to attempt
- `random_state` (int): RNG seed used across splitting and model reproducibility
- `n_jobs` (int): parallel jobs for model training/evaluation
- `rare_class_policy` (string): controls CV behavior when labels are rare; values include `warn_reduce_cv`, `reduce_cv`, and `skip_cv`

`model` (object): model-specific hyperparameters.
- For `random_forest`, supported keys include `n_estimators`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `max_features`, and `class_weight`.
- For `multi_knn`, supported keys include `n_neighbors`, `weights`, `algorithm`, `leaf_size`, `metric`, and `p`.
- Any key under `model` is passed to the underlying scikit-learn estimator constructor (subject to supported list in docs).

`evaluation` (object): controls metrics, reports, and saved artifacts. Keys include:
- `metrics` (string or list, default: `all`)
- `report_class_distribution` (bool, default: true)
- `report_confusion_matrix` (bool, default: true)
- `report_feature_importance` (bool, default: true)
- `save_label_map` (bool, default: true)
- `save_predictions` (bool, default: true)

For full examples and a sample config file, see multi_knn_minimal.yaml and the `Machine Learning` section in [CONFIG.md](CONFIG.md).

The sample config file in [sample_configs/multi_knn_minimal.yaml](sample_configs/multi_knn_minimal.yaml) is a good starting point for both trainers.

## Effects of configuration keys

- `test_size` and `cv_folds` determine holdout and CV behavior; `rare_class_policy` influences fold reduction or skipping when labels are rare.
- `random_state` makes splitting and training deterministic.
- `n_jobs` enables parallel computation for training, permutation importance, and CV utilities.
- `model` keys only affect the chosen estimator and are intentionally namespaced under `model`.
- `evaluation` flags gate both computation and artifact writing.

## Implementation

This section gives per-model implementation details and pointers to the code.

### Random Forest

**Summary:** builds a scikit-learn `RandomForestClassifier` and wraps it in `MultiOutputClassifier` to support the 6 independent binary targets.

**Hyperparameters:** supported via `model` in the config: `n_estimators`, `max_depth`, `min_samples_split`, `min_samples_leaf`, `max_features`, `class_weight`.

**Feature importance:** uses the fitted forest's `feature_importances_` (tree-based importances). When enabled the trainer writes `feature_importances.tsv`. Permutation-based importance may also be computed for diagnostics.

**Evaluation:** computes multi-label metrics (Hamming loss, per-bit precision/recall/F1, micro/macro/weighted F1) and exact-match accuracy as a diagnostic. Per-bit confusion matrices are computed when `report_confusion_matrix` is enabled.

### Multi-label KNN

**Summary:** builds a `KNeighborsClassifier` and wraps it in `MultiOutputClassifier` to expose the same 6-bit interface.

**Hyperparameters:** supported via `model` in the config: `n_neighbors`, `weights`, `algorithm`, `leaf_size`, `metric`, `p`.

**Effective neighbors:** configured `n_neighbors` is clamped to available training size per split/CV fold to avoid invalid neighbor counts.

**Feature importance:** uses permutation importance (`sklearn.inspection.permutation_importance`) when `report_feature_importance` is enabled; results are written to `feature_importances.tsv`.


## Cross-validation and reproducibility

- When `cv_folds` is set, the trainer attempts stratified-like splitting based on label combinations. If exact stratification is impossible due to rare labels, `rare_class_policy` controls fold reduction or CV skipping.
- All splitting uses `random_state` for determinism.


## Outputs

- `*_model.pkl` — Pickled trained `MultiOutputClassifier` (prefix: `random_forest_` or `multi_knn_`)
- `*_metrics.json` — Structured metrics including per-bit and aggregate scores
- `*_metrics.txt` — Human-readable summary of metrics and diagnostic notes
- `label_map.json` — Mapping of bit indices to label names and metadata
- `class_distribution.tsv` — Counts and fractions for each 6-bit label across train/test partitions
- `bit_distribution.tsv` — Positive-counts and fractions for each individual bit across partitions
- `feature_importances.tsv` — Ranked features and importance scores (tree-based for RF, permutation for KNN)
- `predictions.tsv` — Per-row predictions with true/pred bit flags and exact-match indicator


## Notes and diagnostics

- The ML modules prioritise multi-label measures; exact-match accuracy is a strict diagnostic useful for end-to-end checks but not optimised as a primary loss.
- The `evaluation` section in the config lets you toggle both computation and persistence of diagnostic artifacts; disabling an artifact saves time and disk space for large runs.

## See also

- Configuration: [CONFIG.md](CONFIG.md#machine-learning-ghostparserml)

