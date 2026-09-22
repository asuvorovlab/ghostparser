# Machine Learning Module

`ghostparser.ml` trains multi-label classifiers on the summary statistics of
orchestrator runs (typically simulations whose true introgression events are
known), so that a run's statistics can be mapped onto the events that produced
them. Two baselines ship, a random forest and a multi-label k-nearest
neighbours, plus a hyperparameter tuner. Install with `pip install .[ml]`;
every key is in [CONFIG.md](../../CONFIG.md#machine-learning).

## Data contract

The input is `summary_statistics.tsv` (`generate_summary_stats: true`) with a
`class` column holding a 6-bit label per row, one bit per event:

| Bit | Event |
| --- | --- |
| 1 | `ghost_into_A` |
| 2 | `ghost_into_B` |
| 3 | `inflow_into_A_from_C` |
| 4 | `inflow_into_B_from_C` |
| 5 | `outflow_from_A_to_C` |
| 6 | `outflow_from_B_to_C` |

The label is treated as six binary targets, so a prediction can be partly
right: per-label metrics score each bit, and exact-match accuracy scores the
whole string. Every other column is a feature: numeric columns as they are,
string columns with at most seven distinct values one-hot encoded, more are
rejected. Empty cells are rejected. Keep the label as a string in the TSV;
a spreadsheet that strips leading zeros corrupts it.

Bit 1 is the leftmost character, bit 6 the rightmost, so a `class` of `101101`
reads as four events in one triplet: ghost introgression into A (bit 1 set),
none into B (bit 2 clear), inflow into A from C and inflow into B from C (bits
3 and 4 set), no outflow from A to C (bit 5 clear), and outflow from B to C
(bit 6 set). `000000` is a triplet with no introgression at all.

## Usage

```bash
python -m ghostparser.ml.random_forest -i results/summary_statistics.tsv -o ml_out
python -m ghostparser.ml.multi_knn -i results/summary_statistics.tsv -o ml_out --seed 7
python -m ghostparser.ml.hyper_tune -c sample_configs/hyperparameter_tuning_random_forest.yaml
```

The trainers take `-c/--config-file`, `-i`, `-o`, `--seed` and
`--no-overwrite`; the tuner takes `-c`, `--seed` and `--no-overwrite`. A flag
given beside a config file overrides the file's value. Sample configs ship as
`sample_configs/random_forest_minimal.yaml`, `multi_knn_minimal.yaml`,
`hyperparameter_tuning_random_forest.yaml` and
`hyperparameter_tuning_multi_knn.json`;
`orchestrator_species_filter.yaml` shows a run that writes the
`summary_statistics.tsv` they read.

## Training and evaluation

The rows are split into a training and a hold-out partition (`test_size`),
stratified on the label combination. Cross-validation with `cv_folds`
stratified folds runs on the training partition; when a label combination is
too rare to appear in every fold, `rare_class_policy` reduces the fold count
(`warn_reduce_cv`), skips cross-validation (`warn_skip_cv`) or stops
(`error`). Both estimators are wrapped one-vs-rest over the six bits, so each
hyperparameter applies to all six models. `seed` fixes the split, the folds
and the estimators.

`evaluation.metrics` chooses the view: `all` (per-label scores, exact-match
accuracy and the per-bit breakdown), `primary` (aggregate and per-label scores
with Hamming loss), `diagnostic` (exact-match accuracy, the full
classification report and per-bit confusion matrices) or `per_bit`. The
`report_*` and `save_*` keys switch the diagnostic outputs on and off;
feature importance is cheap for the forest and slow for KNN.

**Feature importance.** The random forest reports mean decrease in impurity,
averaged over the six estimators (Breiman 2001, *Machine Learning* 45(1),
5-32, https://doi.org/10.1023/A:1010933404324); it is fast, but biased toward
features with many split points and hard to read for correlated features
(Strobl et al. 2007, *BMC Bioinformatics* 8, 25,
https://doi.org/10.1186/1471-2105-8-25). KNN reports permutation importance
on the hold-out split: each feature is shuffled ten times and the mean drop in
micro-F1 is recorded, which measures the feature's effect on the reported
metric directly but can go negative on small test sets. Both are
model-relative measures of influence, not causal statements.

## Hyperparameter tuning

`hyper_tune` searches the `search_space` of the chosen model by full grid
(`method: grid`, refused above `max_candidates`) or by `n_iter` random draws
(`method: random`), scores each candidate by cross-validated `objective`, and
refits the best on the training partition. The plaintext
`hyper_tune_results.txt` reports the top `top_k` candidates, a per-parameter
marginal table (how each value of each parameter scored across the
candidates that used it), a parameter-influence ranking by how far the
objective moved across each parameter's values, the best value per parameter,
and a warning when a winning value sits at the edge of the range searched,
the cue to extend the range. `hyper_tune_search_report.png` shows the search
progress in evaluation order, the top candidates, and one panel per parameter
that varied.

With `use_wandb: false` (the default) everything is written locally:
`hyper_tune_best_model.pkl`, `hyper_tune_results.{txt,json,tsv}`,
`hyper_tune_parameter_marginals.tsv`, `hyper_tune_search_report.png` and
`predictions.tsv`. With `use_wandb: true` (`pip install .[wandb]`, then
`wandb login`; `WANDB_PROJECT`, `WANDB_ENTITY` and `WANDB_MODE=offline` are
honoured) the TSVs and the JSON are logged to the run as tables and the
output directory keeps the pickle, the text report and the figure.
`wandb_detailed_payloads: true` adds per-candidate fold-level payloads.

## Outputs

| File | Contents |
| --- | --- |
| `<model>_model.pkl` | The fitted one-vs-rest classifier. |
| `<model>_overall_metrics.json`, `<model>_metrics.txt` | Metrics, timings, the hyperparameters used, the dataset summary and the label map. |
| `<model>_confusion_matrices.png` | The six per-bit confusion matrices. |
| `<model>_confusion_matrix_64_classes.png` | The row-normalized 64-class confusion matrix over every 6-bit label, ordered by number of set bits; the diagonal reads as per-class recall. The metrics JSON carries the same matrix and class order. |
| `feature_importances.tsv` | Features ranked by importance. |
| `predictions.tsv` | Per-row true and predicted bits, exact-match flag and matched-bit count. |

To score new data with a saved model:

```python
import pickle
from ghostparser.ml import ml_utils

with open("ml_out/random_forest_model.pkl", "rb") as handle:
    model = pickle.load(handle)
rows = ml_utils.read_tsv_rows("new_summary_statistics.tsv")
matrix = ml_utils.rows_to_matrix(rows, target_column="class")
bits = model.predict(matrix.train_features)
labels = ["".join(str(int(bit)) for bit in row) for row in bits]
```

The new TSV must carry the same feature columns, under the same encoding
rules, as the training TSV.
