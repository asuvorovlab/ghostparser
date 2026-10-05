# Configuration Guide

Every GhostParser key, for `ghostparser.orchestrator` first and the
`ghostparser.ml` trainers and tuner after. The method behind the orchestrator
keys is in [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md); the
trainers are described in [ML.md](ghostparser/ml/ML.md).

## Config files and the command line

Every module accepts `-c/--config-file` with a JSON or YAML file. The file
supplies the settings, and any flag given beside it overrides the file's value
for that setting; the run prints which settings it overrode. Keys without a
flag can only be set in a file. [sample_configs/](sample_configs/) holds a
loadable file for each common scenario to copy from (listed at the end).

To run many datasets with the same settings, keep those in one file and pass
each dataset's paths as flags:

```bash
for folder in /path/to/datasets/*/; do
  [ -f "$folder/genes.tre" ] || continue
  python -m ghostparser.orchestrator -c shared.yaml \
    -gt "$folder/genes.tre" --output-folder "$folder/results"
done
```

Paths are resolved when the config is read: absolute paths as given, `~` to
the home directory, and relative paths from the directory the command runs
in, not from the config file's location. The output folder is reset before
a run under `overwrite: true`, so keep it separate from the directories
holding the input trees.

Every command also takes `--debug`, on the command line only and off by
default: a failed run then shows the full traceback instead of one line. The
exit statuses are listed under Errors in the [README](README.md#errors).

## Orchestrator

### Arguments at a glance

| Argument | Config key | Meaning |
| --- | --- | --- |
| `-st`, `--species-tree-path` | `species_tree_path` | Species tree, Newick. Required. |
| `-gt`, `--gene-trees-path` | `gene_trees_path` | Gene trees, Newick, one per line. Required. |
| `-og`, `--outgroup` | `outgroup` | Outgroup label(s), comma-separated or a list. Required. |
| `--output-folder`, `--no-overwrite` | `output_folder`, `overwrite` | Output directory, reset before the run unless `overwrite` is false. |
| `--triplet-filter` | `triplet_filter` | Analyze only the listed triplets. |
| `--species-filter` | `species_filter` | Analyze every triplet among the listed species. |
| `--species-rename-map` | `species_rename_map` | Names to show in the outputs instead of the tree labels. |
| `--seed` | `seed` | RNG seed for every random draw. |
| `--processes` | `processes` | Worker processes on one machine; `0` = every available CPU. |
| `--alpha-dct`, `--alpha-ks`, `--alpha-perm` | `alpha_dct`, `alpha_ks`, `alpha_perm` | The three gates' thresholds. |
| `--p-value-correction` | `p_value_correction` | `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`. |
| `--diagnostic` | `diagnostic` | Measure every test for every triplet. |
| `--no-consolidation`, `--no-bootstrap` | `consolidation`, `bootstrap` | Skip the maps; skip the bootstrap. |
| `--preflight-data-check`, `--preflight-triplet-cap` | `preflight_data_check`, `preflight_triplet_cap` | Check the inputs and exit; how many triplets the check walks. |
| `--debug` | (command line only) | Show the full traceback when the run fails. |
| (no flag) | `discordant_test`, `tree_height_calculation_strategy`, `min_support_value`, `generate_summary_stats`, `shape_diagnostics`, `bootstrap_options`, `permutation_options` | Config-file only. |

### Input file formats

The trees are Newick. Support values on internal nodes are read for
`min_support_value` and stripped from the processed trees; labels a bare
Newick token cannot hold are single-quoted.

A **triplet filter** names one triplet per line, comma-separated, in the
trees' labels; order within a line does not matter:

```text
A,B,C
A,B,D
C,D,E
```

A **species filter** names taxa, comma-separated, any number per line; every
triplet among them is analyzed:

```text
A, B
C
D, E
```

A **rename map** gives, per tree label, the name to show in the outputs, as a
two-column TSV (`#` comments and blank lines ignored) or a YAML mapping,
chosen by the file extension:

```text
A	Homo sapiens
B	Pan troglodytes
```

```yaml
A: Homo sapiens
B: Pan troglodytes
```

Minimal file:

```yaml
species_tree_path: data/species.tree
gene_trees_path: data/genes.tree
outgroup: OutGroup
output_folder: results
```

Every key at its default, as `sample_configs/orchestrator_full.yaml` ships it:

```yaml
species_tree_path: data/species.tree     # Newick species tree
gene_trees_path: data/genes.tree         # one Newick gene tree per line
outgroup: Out1,Out2                      # comma-separated, or a list
output_folder: results
overwrite: true                          # false = write to a suffixed sibling
triplet_filter: null                     # not with species_filter
species_filter: null
species_rename_map: null
processes: 0                             # 0 = every CPU available to the process
seed: null                               # null = drawn and reported
alpha_dct: 0.05
alpha_ks: 0.05
alpha_perm: 0.05
p_value_correction: bfn                  # no, bfn, holm, fdr_bh, fdr_by
diagnostic: false
consolidation: true
bootstrap: true
preflight_data_check: false
preflight_triplet_cap: 15000             # 0 = no cap
discordant_test: chi-square              # chi-square or z-test
tree_height_calculation_strategy: AVG    # AVG, A, B, C, SIS, INT
min_support_value: 0.5
generate_summary_stats: false
shape_diagnostics: false
bootstrap_options:
  iterations: 100
  diagnostic: false
  summary_only: false
permutation_options:
  min_resamples: 2500
  max_resamples: 25000
  ci_method: wilson
```

### Required

##### `species_tree_path` (`-st`)

Species tree in Newick format. The file must hold exactly one tree. Branch
lengths are optional: the topology alone places the triplets, and lengths
only rank the outgroups (see `outgroup`).

##### `gene_trees_path` (`-gt`)

Gene trees in Newick format, one per line. Every branch should carry a length,
since the tests compare tree heights: a missing length is read as 0, and
`metrics.txt` counts the trees lacking any.

##### `outgroup` (`-og`)

One label, a comma-separated string, or a list. The species tree is rooted
where the outgroups branch off and pruned of them; if other taxa sit between
the outgroups the run stops and names them.

Each gene tree is rooted from its farthest outgroup: the one with the longest
mean path to the ingroup taxa in that gene tree. The outgroups that sit among
the ingroup taxa once the tree is rooted there are pruned without being used,
and the tree is rooted at the common ancestor of the farthest and the others
outside the ingroup.

When two outgroups tie in a gene tree (as they all do in a tree without
branch lengths), the species tree decides: it ranks the
outgroups by the summed branch lengths from the ingroup root, the listed order
breaking exact ties. A species tree lacking any branch length keeps the listed
order. An outgroup the species tree lacks ranks last. `metrics.txt` gives the ranking and counts, per outgroup,
the trees in which it was the farthest and the trees in which it was pruned
unused.

### Config + CLI

##### `output_folder` (`--output-folder`)

Default `results`. Reset before the run under `overwrite: true`; consolidation
writes into its `consolidation/` subfolder.

##### `overwrite` (`--no-overwrite` sets `false`)

Default `true`. With `false` the run writes to the smallest free suffixed
sibling (`results_1`, `results_2`, ...).

##### `triplet_filter` (`--triplet-filter`)

Path to a file of comma-separated triplets, one per line (see *Input file
formats*); only those are analyzed. A triplet naming a taxon that is not an ingroup taxon is skipped
with a warning. Not with `species_filter`.

##### `species_filter` (`--species-filter`)

Path to a file of taxon names, comma-separated, any number per line, spelled
as in the trees (see *Input file formats*). Every triplet among the listed species is analyzed. A name
that is not an ingroup taxon is skipped with a warning; fewer than three left
stops the run. Not with `triplet_filter`.

##### `species_rename_map` (`--species-rename-map`)

Path to a two-column TSV (tree label, display name; `#` comments and blank
lines ignored) or a YAML mapping, chosen by extension (see *Input file
formats*). Read when the config loads, so the file must exist. Display names appear in
the results TSV, `summary_statistics.tsv` and the consolidation outputs; the
outgroup, the filters, the processed trees and `metrics.txt` use the tree
labels. The file is rejected if it maps a label twice, gives two labels one
name, or a name holds a tab, line break, comma, semicolon or `=`.

##### `seed` (`--seed`)

Base seed for every random draw: the direction test, the bootstrap, its own
permutations and the modality bootstrap. Each triplet derives its stream from
`(seed, triplet)`, so results are identical at any worker count. When omitted
a seed is drawn and written to `metrics.txt` as `Seed: <n> (generated)`.

##### `processes` (`--processes`)

Default `0`: every CPU the process may run on, which under a scheduler or a
container is the allocation, not the machine. `1` runs serially. The workers
clean and root the gene trees, build the gene-tree cache, run the triplet
inference and, under `preflight_data_check`, walk the gene trees. Workers are
processes on the machine the run starts on; a job spanning several machines
uses one.

##### `alpha_dct`, `alpha_ks`, `alpha_perm` (`--alpha-dct`, `--alpha-ks`, `--alpha-perm`)

Default `0.05` each: the thresholds of the three gates, compared against the
corrected p-values. `alpha_perm` also bounds the equivalence p-value.

##### `p_value_correction` (`--p-value-correction`)

Default `bfn`; `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`. Applied once across
every triplet to the count-test p-values and again to the tree-height
p-values, and inside each direction test to its two one-tailed p-values.
Bootstrap iterations are judged against the same corrected thresholds. Under
`no` the results carry raw p-values and the flags only. The rank-based
methods hold every triplet's per-iteration p-values until the run finishes
and measure the tree-height test for every triplet, so they cost more memory
and time than `no`/`bfn`. In YAML, `no` may be written bare or quoted.

##### `diagnostic` (`--diagnostic` sets `true`)

Default `false`: a triplet is measured only as far as the cascade reads, so
rows settled by an earlier gate leave the later tests' columns empty
(`perm_note: direction_test_not_consulted`). `true` measures every test for
every triplet and changes no result. It does not reach the bootstrap, whose
own switch is `bootstrap_options.diagnostic`.

##### `consolidation` (`--no-consolidation` sets `false`)

Default `true`: write the introgression maps.

##### `bootstrap` (`--no-bootstrap` sets `false`)

Default `true`. `false` skips the iterations, so `bootstrap_value`,
`all_bootstrap` and the `bootstrap_perm_stat_ci_*` interval are absent and
consolidation weighs every classified triplet as 1.

##### `preflight_data_check` (`--preflight-data-check` sets `true`)

Default `false`. `true` runs only the structural check, writes
`preflight_data_check.txt` and exits; no other output is produced.

##### `preflight_triplet_cap` (`--preflight-triplet-cap`)

Default `15000`; `0` lifts it. The most triplets the check walks; a bound cap
is reported in the check. A `triplet_filter` is never capped.

### Config-file only

##### `discordant_test`

Default `chi-square`; or `z-test`. Pearson's chi-square, or a two-proportion
z-test for which `z^2 = 2 chi^2` on the same counts, so it rejects more
readily.

##### `tree_height_calculation_strategy`

Default `AVG`; `A`, `B`, `C`, `SIS`, `INT`. `AVG` averages the three
root-to-tip distances of the triplet, `A`/`B`/`C` take one taxon's, `SIS` is
the patristic distance between the sisters and `INT` the internal branch from
the triplet's root to the sisters' node.

##### `min_support_value`

Default `0.5`. Trees whose mean internal-node support is below it are
dropped; trees without support labels are kept.

##### `generate_summary_stats`

Default `false`. `true` also writes `summary_statistics.tsv`: mean, median,
mode, variance, entropy, minimum and maximum of the average tree height, the
internal branch and the sister distance, per topology group (63 columns).

##### `shape_diagnostics`

Default `false`. `true` adds the mode count, Silverman modality p-value,
skewness, excess kurtosis and generalized-Pareto tail index of each height
group to the results TSV. Descriptive only; the modality bootstrap makes it
expensive.

##### `bootstrap_options`

Nested, or flat as `bootstrap_<key>`.

- `iterations`: default `100`, integer `>= 1`.
- `diagnostic`: default `false`. Measures all three tests in every iteration
  and writes them per iteration (`bootstrap_dct_*`, `bootstrap_ks_*`,
  `bootstrap_perm_*`, `bootstrap_con_mean`, `bootstrap_dis_mean`,
  `bootstrap_gene_tree_heights`) without moving any vote. Costly; pair it
  with a filter.
- `summary_only`: default `false`. With `diagnostic`, write per-column
  summaries (`count`, `non_null_count`, `mean`, `median`, `min`, `max`)
  instead of full lists.

##### `permutation_options`

- `min_resamples`: default `2500`. The first batch and the minimum total.
- `max_resamples`: default `25000`, at least `min_resamples`. The budget;
  the batch that crosses it is drawn whole, so `perm_n_resamples` can exceed
  it by up to one batch. Keep the two a few multiples apart.
- `ci_method`: default `wilson`; `wilson`, `beta`, `agresti_coull`,
  `jeffreys`, `binom_test`, `normal`. The binomial interval of the stopping
  rule; Wilson keeps its coverage near the small p-values the test produces,
  `normal` does not, and `beta` (Clopper-Pearson) is conservative and
  resamples longer.

Bootstrap iterations run the direction test at a fifth of both resample
counts.

## Machine learning

The trainers read `summary_statistics.tsv`, treat the 6-bit `class` column
as six binary labels, use every other column as a feature (numeric columns
directly, string columns with at most seven distinct values one-hot encoded,
more are rejected), and report per-label and exact-match metrics. Install with
`pip install .[ml]`.

Flags: `-c/--config-file`, `-i/--input-path`, `-o/--output-dir`, `--seed`,
`--no-overwrite`, `--debug`. The tuner takes `-c`, `--seed`, `--no-overwrite`
and `--debug`.

```yaml
input_path: ./results/summary_statistics.tsv
output_dir: ./results/ml_out
overwrite: true
target_column: class                 # the 6-bit label column
test_size: 0.2                       # hold-out fraction, in (0, 1)
cv_folds: 5                          # integer >= 1, or null for no CV
rare_class_policy: warn_reduce_cv    # warn_reduce_cv, warn_skip_cv, error
seed: null
n_jobs: -1                           # -1 = every core

model:                               # random forest; the KNN keys are below
  n_estimators: 200
  max_depth: null                    # null = unlimited
  min_samples_split: 2
  min_samples_leaf: 1
  max_features: null                 # null = every feature; sqrt, log2, int, float in (0, 1]
  class_weight: null                 # null, balanced, balanced_subsample, mapping, or list

evaluation:
  metrics: all                       # all, primary, diagnostic, per_bit
  report_class_distribution: true
  report_confusion_matrix: true
  report_feature_importance: true
  feature_importance_method: null    # null = mdi on the forest, permutation on KNN
  feature_importance_correlation_threshold: 0.7   # grouped_permutation only
  save_label_map: true
  save_predictions: true
```

For `multi_knn` the `model` block holds `n_neighbors` (`5`), `weights`
(`uniform` or `distance`), `algorithm` (`auto`, `ball_tree`, `kd_tree`,
`brute`), `leaf_size` (`30`), `metric` (`minkowski`) and `p` (`2`).

Runtime keys sit at the top level, hyperparameters under `model`, reporting
controls under `evaluation`; a key in the wrong section is rejected.
`rare_class_policy` says what to do when a label combination is too rare to
stratify: reduce the fold count, skip cross-validation, or stop. Write null
as `null` or `~`, never the bare word `None`; `max_features: auto` is
rejected (scikit-learn removed it; `sqrt` is the equivalent). Omit
`min_samples_split` and `min_samples_leaf` to take their defaults.

`feature_importance_method` chooses how influence is measured, and the three
answers differ when features are correlated, as the summary statistics of one
height distribution are:

| Value | Measure |
| --- | --- |
| `mdi` | Mean decrease in impurity, summed over the splits a feature makes and averaged over the six estimators. Fast, read off the fitted trees, but scored in-sample: it favours features with many split points, and correlated features divide the credit for one signal between them. Needs a tree-based model, so the neighbours classifier rejects it. |
| `permutation` | The drop in hold-out micro-F1 when a feature's column is shuffled. Scores the metric actually reported, on rows the model never saw. Correlated features still score low, because shuffling one leaves its twin to carry the signal. |
| `grouped_permutation` | The same shuffle applied to a whole group of correlated features at once, so one signal is scored once. Features are grouped by average-linkage clustering on `1 - abs(Spearman rho)`, cut at `feature_importance_correlation_threshold`; a lower threshold builds larger groups. |

Null takes each trainer's own measure: `mdi` for `random_forest`, which reads
impurity off its trees, and `permutation` for `multi_knn`, which has no
impurity to read. The permutation measures shuffle ten times and report the
mean drop with its standard deviation; on a small hold-out partition the
scores are noisy, and an uninformative feature can score slightly negative.
No measure is unbiased for correlated features in the strict sense, because
two features carrying one signal have no unique split of the credit between
them; `grouped_permutation` sidesteps the question by scoring the signal
rather than the columns.

### Hyperparameter tuning

`python -m ghostparser.ml.hyper_tune -c file` reads the runtime keys above
plus a `hyperparameter_tuning` block, and no `model` or `evaluation` section:

```yaml
hyperparameter_tuning:
  model: random_forest             # random_forest or multi_knn
  method: grid                     # grid or random
  objective: exact_match_accuracy  # exact_match_accuracy, hamming_loss, bitwise_accuracy, micro_f1, macro_f1, weighted_f1
  top_k: 10                        # candidates in the text summary
  n_iter: 20                       # candidates sampled by random
  max_candidates: 5000             # cap on a full grid
  use_wandb: false                 # true = log to Weights & Biases
  wandb_detailed_payloads: false   # per-candidate CV payloads; needs use_wandb
  search_space:                    # the model's hyperparameters, each a list
    n_estimators: [100, 200, 400]
    max_depth: [null, 10, 20]
    max_features: [sqrt, log2]
```

A parameter omitted from `search_space` keeps the trainer default; runtime
keys inside it are rejected. With `use_wandb: true` (`pip install .[wandb]`,
`wandb login`) the ranked candidates, parameter marginals, predictions and the
results JSON are logged to the run instead of written to disk.

## Sample configs

Every file under `sample_configs/` loads as shipped; change the paths and the
outgroup labels to your data.

| File | Scenario |
| --- | --- |
| `orchestrator_minimal.yaml` | The three required inputs and an output folder; everything else at its default. |
| `orchestrator_full.yaml` | Every orchestrator key at its default, with a comment on each; trim it to the keys you change. |
| `orchestrator_preflight.yaml` | Check the trees and the outgroup order before a run; nothing else is produced. |
| `orchestrator_species_filter.yaml` | Every triplet among a chosen set of species, a fixed seed, and `summary_statistics.tsv` for the trainers. |
| `orchestrator_triplet_filter_diagnostic.yaml` | A few named triplets in depth: every test for every triplet, every bootstrap iteration recorded, shape diagnostics, a different height strategy. |
| `orchestrator_screen.json` | A fast screen of a large taxon set, in JSON: no bootstrap, no maps, `fdr_bh`, the z-test, a smaller permutation budget. |
| `random_forest_minimal.yaml`, `multi_knn_minimal.yaml` | One trainer run each on a run's `summary_statistics.tsv`. |
| `hyperparameter_tuning_random_forest.yaml` | A grid search over the random forest, reported locally. |
| `hyperparameter_tuning_multi_knn.json` | A random search over the KNN, in JSON. |
