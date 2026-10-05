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
`report_*` and `save_*` keys switch the diagnostic outputs on and off.

## Feature importance

`evaluation.feature_importance_method` chooses between three measures, which
disagree whenever features are correlated. They are all model-relative
measures of influence, not causal statements.

**`mdi`**, the forest's default, is mean decrease in impurity (Breiman 2001,
*Machine Learning* 45(1), 5-32,
https://doi.org/10.1023/A:1010933404324): every node that split on a feature
contributes its impurity drop weighted by the rows reaching it, summed over
the forest, averaged over the six estimators and normalized to sum to 1. It
comes free with the fitted trees, but it is measured in-sample on in-bag
rows, so it rewards a feature for the noise it fits, in proportion to how
many candidate split points it offers (Strobl et al. 2007, *BMC
Bioinformatics* 8, 25, https://doi.org/10.1186/1471-2105-8-25).

Its other failure is the one that matters for these features. At each node
the forest draws a random subset of features (`max_features`) and uses one of
them, so when several features carry the same information each node books the
whole impurity drop to whichever was drawn, and the credit for one signal is
divided among its carriers roughly in proportion to how often each is drawn
and chosen. Masking compounds it: a feature used high in a tree leaves
nothing for its twin further down the same path. With `k` redundant carriers
each scores about `1/k` of what the signal scores alone, so a real signal can
rank below a unique but weaker feature, and the order within the redundant
set moves with the seed. The input invites exactly this: the summary
statistics are three topologies by three metrics by seven statistics, and the
mean, median, mode, minimum and maximum of one height distribution are five
views of one quantity.

**`permutation`**, the neighbours classifier's default, shuffles one feature's
column in the hold-out partition and records the drop in micro-F1, ten times
per feature (Breiman 2001; Fisher, Rudin & Dominici 2019, *Journal of Machine
Learning Research* 20(177), 1-81, https://jmlr.org/papers/v20/18-760.html).
Scoring held-out rows against the metric the run reports removes both
in-sample biases, and the score is stated in the units of that metric. It
does not remove the redundancy problem, and arguably sharpens it: with the
twin column left intact the model reads the signal off it, so both twins look
unimportant.

**`grouped_permutation`** answers the question the redundancy makes
ambiguous, by scoring a group of features rather than a column. Features are
clustered by average-linkage hierarchical clustering on the distance
`1 - abs(Spearman rho)`, cut at
`evaluation.feature_importance_correlation_threshold`, and every column of a
group is shuffled by the same row order, which breaks the group's link to the
label while leaving the correlations inside the group intact. One signal is
then scored once, whichever of its carriers the model happened to use. The
output names each feature's group and that group's size beside the score the
group earned.

Neither permutation measure is unbiased for correlated features in the strict
sense, and none can be: when two features carry one signal, "what would I
lose without this column" and "what share of the signal is this column" are
different questions with different right answers. Both permutation measures
also evaluate the model on shuffled rows that the original feature
distribution would rarely produce. The scores carry a standard deviation over
the ten shuffles, and on a small hold-out partition an uninformative feature
can come back slightly negative.

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
the cue to extend the range. `figures/hyper_tune_search_report.png` shows the search
progress in evaluation order, the top candidates, and one panel per parameter
that varied.

With `use_wandb: false` (the default) everything is written locally:
`hyper_tune_best_model.pkl`, `hyper_tune_results.{txt,json,tsv}`,
`hyper_tune_parameter_marginals.tsv`, `figures/hyper_tune_search_report.png`
and `predictions.tsv`. With `use_wandb: true` (`pip install .[wandb]`, then
`wandb login`; `WANDB_PROJECT`, `WANDB_ENTITY` and `WANDB_MODE=offline` are
honoured) the TSVs and the JSON are logged to the run as tables and the
output directory keeps the pickle, the text report and the figure.
`wandb_detailed_payloads: true` adds per-candidate fold-level payloads.

## Outputs

| File | Contents |
| --- | --- |
| `<model>_model.pkl` | The fitted one-vs-rest classifier. |
| `<model>_overall_metrics.json`, `<model>_metrics.txt` | Metrics, timings, the hyperparameters used, the dataset summary and the label map. |
| `figures/<model>_confusion_matrices.png` | The six per-bit confusion matrices, each row normalized to the fraction of its true bit, with counts in brackets. |
| `figures/<model>_confusion_matrix_64_classes.png` | The row-normalized 64-class confusion matrix over every 6-bit label, ordered by number of set bits and divided into blocks by that number; the colour scale is square-root, so small off-diagonal fractions stay visible. The diagonal reads as per-class recall. The metrics JSON carries the same matrix and class order. |
| `figures/<model>_per_bit_accuracy.png` | Accuracy of each bit on the hold-out partition. |
| `figures/<model>_per_class_accuracy.png` | Per-class recall (the 64-class diagonal) grouped by the number of set bits in the true class, with each group's mean. A class with no hold-out rows is left out rather than scored 0. The 1/64 chance line is drawn only when all 64 classes occur in both partitions and no class count exceeds 1.5 times another, since it assumes equal class weights. |
| `feature_importances.tsv` | Features ranked by importance: `feature` and `importance` under every measure, `importance_std` over the shuffles under the permutation measures, and each feature's `group` and `group_size` under `grouped_permutation`. The metrics JSON names the measure used. |
| `predictions.tsv` | Per-row true and predicted bits, exact-match flag and matched-bit count. |

The figures need `evaluation.report_confusion_matrix` and an `evaluation.metrics` of `diagnostic` or `all`; without them no `figures/` folder is created. `<model>_metrics.txt` lists the path of each figure written.

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
