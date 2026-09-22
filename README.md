# Ghostparser

<p align="center">
   <img width="120" height="120" alt="GhostParser icon" src="docs/assets/ghostparser_icon.svg" />
</p>

## Overview

**Ghostparser** is a phylogenetic introgression pipeline built around [`ghostparser.orchestrator`](#orchestrator-primary-entry-point), with an optional machine-learning module, [`ghostparser.ml`](#machine-learning-ghostparserml), for training models on simulated summary statistics.

The core pipeline implements the DCT+THT workflow for detecting sampled and ghost introgressions in large trees. Because the method operates on species triplets, the input trees are decomposed into all possible triplets by default, or narrowed with [`triplet-filter`](#configuration) when you only want to analyze specific triplets. Results are consolidated into a heatmap (for sampled directed introgressions) and bar chart (for ghost introgression where only target of introgression is inferred). The pipeline supports both command-line and configuration-file modes.

The optional machine-learning subpackage trains multi-label classifiers on a run's summary statistics and includes utilities for model selection and tuning.

---

## Jump to Sections:

1. [Quick Start](#quick-start)
2. [Modules](#modules)
   - [Orchestrator (Primary Entry Point)](#orchestrator-primary-entry-point)
   - [Machine Learning](#machine-learning-ghostparserml)
3. [Orchestrator Input/Output](#orchestrator-inputoutput)
4. [Configuration](#configuration)
   - [Defaults](#defaults-at-a-glance)
6. [Handled Errors](#handled-errors)
7. [Testing](#testing)
8. [For Maintainers](#for-maintainers)

---

## Quick Start

#### Installation from Wheel (Recommended for Users)

Pre-built wheel distributions are automatically created and published when a release is tagged. Users can download and install them without building locally.

**Downloading the wheel from GitHub:**

1. Visit the [Releases page](https://github.com/asif256000/ghostparser/releases)
2. Find the latest release (e.g., `v0.1.0`)
3. Download the `.whl` file (e.g., `ghostparser-0.1.0-py3-none-any.whl`)

**Installing from wheel:**

```bash
# Install from downloaded wheel file
pip install ghostparser-0.1.0-py3-none-any.whl

# OR install directly from GitHub releases URL
pip install https://github.com/asif256000/ghostparser/releases/download/v0.1.0/ghostparser-0.1.0-py3-none-any.whl

# OR install from source distribution
pip install https://github.com/asif256000/ghostparser/releases/download/v0.1.0/ghostparser-0.1.0.tar.gz
```

**After installation, users can run from any directory:**

```bash
python -m ghostparser.orchestrator -c config.yaml

python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup

python -m ghostparser.ml.random_forest -i summary_statistics.tsv -o ml_out

python -m ghostparser.ml.multi_knn -i summary_statistics.tsv -o ml_out
```

**Note:** The wheel installation installs the package into your Python environment, so you don't need to be in the project directory to run it. All module commands (`python -m ghostparser.*`) work from anywhere.


### Installation via pip (For Developers)

```bash
pip install -r requirements.txt
```

If you want the machine-learning baselines, install the optional ML extra instead of the core-only package:

```bash
pip install .[ml]
```

That extra pulls in the scikit-learn dependency used by `ghostparser.ml.random_forest` and `ghostparser.ml.multi_knn`.

The same extra also enables `ghostparser.ml.hyper_tune`, which reads a `hyperparameter_tuning` config section and can run either grid search or random search across the supported ML trainers.

Run the tuner directly with:

```bash
python -m ghostparser.ml.hyper_tune -c sample_configs/hyperparameter_tuning_random_forest.yaml
```

Weights & Biases logging is opt-in: `hyperparameter_tuning.use_wandb` defaults to `false`, and the tuner then runs entirely locally, writing a plaintext report, ranked-candidate and per-parameter TSVs, and a search-report figure. Set it to `true` to log to W&B instead, which ships as its own extra:

```bash
pip install .[ml,wandb]
```

#### Poetry 2.x+ Alternative

If you prefer Poetry 2.x+ instead of plain pip:

```bash
poetry install
```

Run commands with Poetry:

```bash
poetry run pytest -q
poetry run python -m ghostparser.orchestrator -c run_config.yaml
```

---

## Modules

### Orchestrator (Primary Entry Point)

`ghostparser.orchestrator` is the introgression inference engine.

Run with a config file (recommended for reproducibility):

```bash
python -m ghostparser.orchestrator -c run_config.yaml
```

Run via CLI flags:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

With an explicit worker count:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup \
    --processes 4
```

#### How the Orchestrator Works

1. The species tree is standardized, filtered on mean internal support, rooted where the outgroups branch off, and pruned of them; gene trees are cleaned and each is rooted where the largest set of its outgroups branches off, then pruned of them.
2. Every ingroup triplet is enumerated (or restricted to the triplets named by `--triplet-filter`, or to every triplet among the species named by `--species-filter`) and normalized to `(A, B, C)` with A and B the species-tree sisters.
3. Every gene tree is measured once up front, and for each triplet the engine reads that triplet's rooted shape back out of those measurements: it classifies the topology as concordant or one of two discordant alternatives, and records a tree height H(T).
4. Three tests are measured per triplet: the discordant count test, the KS tree-height test, and a studentized permutation test on the concordant-versus-discordant1 mean heights.
5. Once every triplet is in, multiple-testing correction is applied across the whole run and a three-gate decision reads the corrected results in that order. Each triplet lands on `no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, or `ambiguous`.
6. Results are written, and consolidation renders the introgression maps.

See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md) for the mechanism in detail.

#### Arguments

**Required:**

- `-st, --species-tree-path`
- `-gt, --gene-trees-path`
- `-og, --outgroup`

**Config-file mode (CLI-only):**

- `-c, --config-file` — when given, the file supplies every setting and the other CLI flags are ignored with a warning.

**Config + CLI:**

- `--output-folder`, `--no-overwrite`, `--triplet-filter`, `--species-filter`, `--species-rename-map`
- `--seed`

   - Base RNG seed for the whole run; every random draw derives from it.
   - Omit for a fresh seed each run; the value used is reported in `metrics.txt`.

- `--processes`
- `--alpha-dct`, `--alpha-ks`, `--alpha-perm`, `--p-value-correction`
- `--no-consolidation`, `--no-bootstrap`, `--diagnostic`
- `--preflight-data-check`, `--preflight-triplet-cap`

**Config-file only:** `discordant_test`, `tree_height_calculation_strategy`, `min_support_value`, `generate_summary_stats`, `shape_diagnostics`, the `bootstrap_options` block (`iterations`, `diagnostic`, `summary_only`) and the `permutation_options` block (`min_resamples`, `max_resamples`, `ci_method`).

The statistical tests always use the scipy/statsmodels backend. The full key reference is in the **[Configuration Guide](CONFIG.md#orchestrator-primary-module)**.

#### Checking Your Data First

Before committing to a full run, `--preflight-data-check` validates the input trees and exits without any analysis:

```bash
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup --preflight-data-check
```

It writes `preflight_data_check.txt` into the output folder, listing every structural problem it found — gene trees missing the outgroup, polytomous triplets with no resolvable sister pair, triplet-filter lines naming unknown taxa — with counts per category, examples naming the offending gene tree and triplet, and a summary attributing the issues to the species tree or the gene trees. It also counts how many gene trees rooted on each outgroup and how many carry outgroups that do not all lie on one side of the other taxa. These are the failures that would otherwise surface as errors partway through a long run. The checks are structural only: passing means the data can be processed, not that the result will be biologically meaningful. By default the check walks at most 15,000 triplets and says so in the report when that cap binds; `--preflight-triplet-cap` (or `preflight_triplet_cap`) raises it, and `0` checks every triplet.

#### Primary Outputs

1. `orchestrator_triplet_results.tsv` — per-triplet classification results
2. `summary_statistics.tsv` — only when `generate_summary_stats` is set
3. `processed_<species tree file>` / `processed_<gene trees file>` — cleaned, rooted trees
4. `metrics.txt` — per-stage wall/CPU timing and run parameters
5. `consolidation/` — the combined figure, with the TSV matrices under `consolidation/consolidation_data/

The output folder is reset at the start of a run under the default `overwrite: true`, so point it at a directory of its own — never at a directory holding your input trees, which would be deleted with it. `--no-overwrite` writes to an auto-suffixed sibling instead.

Consolidation details:

- The figure uses the `cividis` colormap throughout. The colorbar applies to the inflow/outflow heatmap. In the ghost bar chart, bar *length* encodes the ghost bootstrap value while colour encodes only whether that taxon also has sampled introgression, using the two extremes of the same colormap: yellow for ghost-only, dark blue for ghost plus sampled. A legend above the bars states the mapping.
- Heatmap cells with no introgression edge are left unpainted, so a sparse matrix shows its real signal rather than a wall of colour. The gridlines and panel border still mark the row and column structure.
- The species tree topology is stitched onto the plot axes so the source and target axes read like tree labels.
- By default the plots use the processed species tree after outgroup pruning; with a triplet filter or a species filter, the plotted tree can be pruned to the filtered taxa.
- Consolidation is enabled by default; disable it with `--no-consolidation`.
- Its artifacts go in a dedicated `consolidation/` subfolder so its own output-directory reset cannot remove the run's results.

#### Where the Decisions Are Made

- Nothing is classified while the triplets are being processed. Once every triplet is in, the count p-values are corrected together across all triplets, the tree-height p-values likewise, and only then does the cascade read the corrected flags to classify each triplet. The direction test's p-values are corrected inside the test itself (across its two one-tailed p-values) and never across triplets, so they are complete as soon as the test finishes. See [CONFIG.md](CONFIG.md#p_value_correction).
- The `decision_gate` column names the test that settled the call. A test reported on the same row below that gate took no part in it.

#### How Much Gets Computed

- By default a triplet is measured only as far as the decision cascade reads. Only the third gate consults the permutation direction test, so a triplet an earlier gate already decided skips it instead of computing a result nothing reads. The tree-height test is skipped below a settled count gate under `no` and `bfn`, whose correction of each p-value depends on nothing but the triplet count; under `holm`, `fdr_bh` and `fdr_by` it is measured for every triplet, because those corrections rank every triplet's value against the others' and need every one of them. See [CONFIG.md](CONFIG.md#diagnostic).
- `--diagnostic` (`diagnostic: true`) measures all three tests for every triplet, so every `ks_*` and `perm_*` column is filled. Use it when you want the direction test's statistics everywhere for debugging. Expect the default to be somewhat faster, not dramatically so: most of a run's time is the bootstrap's permutation tests, and every iteration of a triplet that reaches the third gate runs one either way.
- **This changes no classification, under any correction method.** Two things guarantee it. No supported correction can lower a p-value below its raw value, so a gate that failed raw cannot clear once corrected — the skipped test could never have been reached. And the direction test's p-values are corrected inside the test, never across triplets, so leaving one triplet's test unrun moves no other triplet's numbers; that is what makes it safe to skip triplet by triplet while the run is still streaming, before the run-wide correction has seen the whole family.
- The bootstrap is not affected by `diagnostic`: every iteration measures only what its vote reads, under the same per-correction rule. `bootstrap_options.diagnostic` is its own switch: it measures all three tests in every iteration and writes them per iteration, without moving a single vote.
- Skipped triplets have an empty `perm_*` block and carry `perm_note: direction_test_not_consulted`, so a deliberate skip is distinguishable from a test that ran and hit a guard. Under `no`/`bfn` a triplet the count gate settled also leaves the tree-height columns empty. `metrics.txt` reports the setting, says what the run skips under its correction, and counts the skipped tests and the triplets clearing each gate.

#### Direction Test Behavior

- The third decision gate is an adaptive studentized permutation test on the concordant versus discordant1 mean tree heights. It resamples until a confidence interval around the p-value excludes `alpha_perm`, or until the total reaches `max_resamples` — `metrics.txt` counts the triplets that spend the budget, and the `perm_converged` column identifies them. The final batch is drawn whole rather than trimmed, so the reported resample count can sit just above `max_resamples`.
- Its tuning knobs (`min_resamples`, `max_resamples`, `ci_method`) live in the config file's `permutation_options` block.
- The direction is read off the corrected one-tailed p-values, so the results TSV carries no per-group mean or median columns. Descriptive per-group statistics live in `summary_statistics.tsv` (`generate_summary_stats`).
- Triplets whose samples are too small or too degenerate to support the test are reported as `inconclusive` with the reason in the `perm_note` column, rather than being given a direction the data cannot justify.
- When neither direction is significant, the run does not stop at "no difference found": failing to find a difference is not evidence that there is none, so whether the mean heights were *shown* to be close needs its own test. That test is TOST — two one-sided tests, the standard equivalence procedure. It turns the null around: instead of asking whether the difference is zero, it states two nulls, that the concordant mean is at least a margin *below* the discordant mean and that it is at least a margin *above* it, and tests each one-sidedly; rejecting both confines the difference to within the margin. The margin is a Cohen's *d* of 0.5, half a pooled standard deviation. Rejecting both gives `perm_decision = equivalent` (with the p-value in `perm_p_tost`); otherwise nothing was established and the decision is `inconclusive`. Both classify the triplet as `ambiguous`.
- TOST rides on the direction test's own resamples rather than drawing its own. Its two nulls sit a margin either side of zero, which is a constant added to one group, and a shifted permutation's statistic follows from the unshifted power sums plus two extra reductions — so the equivalence step costs a fraction of a pass instead of two more, and both questions are answered at the same Monte Carlo resolution by construction.
- `bootstrap_perm_stat_ci_low` / `bootstrap_perm_stat_ci_high` give a bootstrap-percentile interval on the studentized mean difference, so a direction can be read as an effect size rather than only as a threshold crossing. It comes from the bootstrap loop, so `--no-bootstrap` leaves it empty.
- See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md#the-statistical-tests) for the method, its citations, and its known small-sample limitation.

#### Bootstrap Behavior

- Bootstrap is enabled by default and can be disabled with `--no-bootstrap`, which skips the iterations rather than merely dropping their columns; the remaining controls (`iterations`, `diagnostic`, `summary_only`) are set through the config file's `bootstrap_options` block.
- Bootstrap iterations re-run the direction test at one fifth of the configured resample budget.
- Iterations are judged against the same *corrected* p-value thresholds as the reported classification, so `bootstrap_value` measures support for the decision actually made. With `no` or `bfn` the correction is applied as each iteration runs; the rank-based methods need every triplet's p-value for the same iteration, so those are corrected after the streaming pass.
- An iteration whose direction test hits a guard votes `ambiguous`, since the guard reports `inconclusive`.
- `bootstrap_value` reports the bootstrap fraction for the final `classification` value after correction.
- The run-wide `--seed` makes results reproducible and identical at any worker count, because each triplet derives its own stream from it. It seeds every random draw in the run, not just the bootstrap. When omitted, a seed is drawn and reported in `metrics.txt`, so any run can be reproduced from its own log.

#### Orchestrator Input/Output

##### Input Expectations

**Required Arguments:**

- `--species-tree-path` (alias `-st`)

   - Path to the species tree file in Newick format.

- `--gene-trees-path` (alias `-gt`)

   - Path to the gene trees file in Newick format.

- `--outgroup` (alias `-og`)

   - Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups; a config file may also give a list.

**Optional Arguments:**

- `--output-folder`

   - Output folder for the run (default: `./results`).

- `--triplet-filter`

   - Path to a triplet filter file (comma-separated taxa per line).
   - When provided, only those triplets are processed. Triplets naming a taxon
     missing from the species tree are skipped with a warning.

- `--species-filter`

   - Path to a species filter file: taxon names as they appear in the trees,
     comma-separated, any number per line.
   - When provided, every triplet among the listed species is processed and
     no other. A name not found among the ingroup taxa — misspelled, absent
     from the species tree, or an outgroup — is skipped with a warning, and
     the run stops if fewer than three names remain.
   - Cannot be combined with `--triplet-filter`; naming both is a config error.

- `--species-rename-map`

   - Path to a two-column TSV or YAML mapping from the taxon labels used in the
     trees to the names that should appear in the outputs.
   - The display names appear in the results TSV, `summary_statistics.tsv`, and
     the consolidation matrices and plot; the outgroup, a triplet or species
     filter and the processed tree files use the tree labels. Taxa absent from the map keep
     their tree labels. See [CONFIG.md](CONFIG.md).

- `--processes`

   - Number of worker processes.
   - Defaults to `0`: every CPU the process is allowed to run on, which under a
     container limit or a job scheduler can be fewer than the machine has. Use
     `--processes 1` to run serially in the parent process, which is useful for
     debugging or constrained systems. `metrics.txt` reports the count actually
     used.
   - Workers are processes on one machine. On a cluster, a job that spans
     several nodes runs the orchestrator on a single node and leaves the
     others idle; no multi-node execution is implemented.

##### Output Files

An orchestrator run generates these output files:

1. **`processed_<species tree file>`** - Processed species tree, named after the input file, with support values removed and the outgroup rooting and pruning applied
2. **`processed_<gene trees file>`** - Processed gene trees, named after the input file, with support values removed and the outgroup rooting and pruning applied
3. **`metrics.txt`** - Metrics log with warnings, timings, and counts
4. **`orchestrator_triplet_results.tsv`** - Final triplet-level classification results (`no_introgression`, `outflow_introgression`, `inflow_introgression`, `ghost_introgression`, or `ambiguous`)
5. **`summary_statistics.tsv`** - Optional per-triplet summary table, written only when `generate_summary_stats` is enabled. It carries no shape columns; those go to the results TSV alone. Its `discordant1_*` columns describe whichever discordant topology is more frequent (matching the `dis1_topology` column) and `discordant2_*` the other. It includes:
   - identity columns (`triplet`, `abc_mapping`, `species_tree`, `dis1_topology`)
   - topology counts (`n_con`, `n_dis1`, `n_dis2`)
   - 63 topology/metric summary columns (7 statistics x 3 topology classes x 3 metric types)
   - final `classification` and `bootstrap_value` (when bootstrap is enabled)
6. **`consolidation/`** - Introgression map figure, with the TSV matrices under `consolidation/consolidation_data/`

Base TSV output includes `dis1_topology` and a topology-only `species_tree` value for each triplet.
Base TSV output also includes a `decision_gate` column naming which test settled the classification (`DCT`, `THT`, or `PERM`). Only `PERM` means `perm_decision` was actually consulted: by default the permutation columns are filled in only on those rows, and under `diagnostic: true` they are filled in for every triplet, so `decision_gate` is what tells you whether they took part.
Base TSV output also includes an `inference` column with human-readable direction text using actual species names.

When bootstrap is enabled, the TSV adds:

- `bootstrap_value`
- `all_bootstrap`

When the bootstrap is diagnostic (`bootstrap_options.diagnostic: true`), every iteration measures all three tests and the TSV also adds, one entry per iteration:

- `bootstrap_dct_stats`, `bootstrap_dct_p_value`
- `bootstrap_ks_stats`, `bootstrap_ks_p_value`
- `bootstrap_perm_stats`, `bootstrap_perm_p_greater`, `bootstrap_perm_p_less` — the studentized statistic and the raw one-tailed p-values of the iteration's direction test
- `bootstrap_perm_decisions` — `greater`, `less` or `inconclusive` per iteration
- `bootstrap_con_mean`, `bootstrap_dis_mean`
- `bootstrap_gene_tree_heights`

The columns are JSON strings. With `summary_only: true` each numeric column holds a `count`/`non_null_count`/`mean`/`median`/`min`/`max` summary and `bootstrap_perm_decisions` a count per decision. A diagnostic bootstrap runs a direction test in every iteration of every triplet, including the ones an earlier gate settled, so it costs substantially more than a plain run; use it on a filtered set of triplets or species.

#### Example Usage

**Species tree** (`species.tree`):

```
(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);
```

**Gene trees** (`genes.tree`):

```
((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

**Run:**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

**Multiple outgroups (comma-separated):**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup1,OutGroup2
```

**With triplet filter:**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --triplet-filter triplets.txt
```

**With species filter** (every triplet among the listed species):

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --species-filter species.txt
```

The outgroups must branch off the species tree at a single point: read unrooted, one branch (or, inside a polytomy, one node) parts every outgroup from every other taxon, whatever orientation the file was written in. The tree is rooted there, the outgroups are pruned, and the remaining taxa are the ingroup. If other taxa sit between the outgroups -- so that rooting on the outgroups would leave those taxa on different sides of them -- the run stops and names the groups those taxa fall into, so that you can add the groups that are outgroups to `--outgroup`, or correct the species tree, and rerun. Each gene tree is rooted the same way, on the outgroups it carries, and pruned of them; a gene tree carrying no outgroup, or nothing but outgroups, is dropped. When the outgroups present all lie on one side of the other taxa, they root the tree together. When they do not -- a single gene often places a distant outgroup on a long branch somewhere inside the ingroup -- the largest set of them that still parts from the other taxa at one point roots the tree, and an outgroup tangled among the ingroup taxa is pruned without being used: the outgroups that still sit together outvote the one that strayed, and rooting on the stray would have moved the root into the ingroup and changed the shape of every triplet spanning the two. When no set holds a majority -- two outgroups apart, for instance -- the one listed earliest wins, so list first the outgroup you trust most: close enough to the ingroup to be placed well by a single gene, far enough to be outside it in every gene, with a distant long-branch outgroup last as the fallback for trees carrying no other. `metrics.txt` reports how many gene trees were rooted using each outgroup, in how many each was tangled and pruned unused, how many trees had a tangled outgroup at all, and in how many of those the listed order chose.

**Output** (`orchestrator_triplet_results.tsv`, abbreviated):

```
triplet	species_tree	n_con	n_dis1	n_dis2	dis1_topology	classification	bootstrap_value
TaxaA,TaxaB,TaxaC	((TaxaA,TaxaB),TaxaC);	7	3	2	BC	no_introgression	0.82
TaxaA,TaxaC,TaxaD	((TaxaA,TaxaC),TaxaD);	12	0	0	BC	no_introgression	1.00
```

One row per triplet; the full column list is documented in [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md).

---


### Machine Learning (ghostparser.ml)

A small machine-learning baseline lives under `ghostparser.ml`. It consumes `summary_statistics.tsv` (the optional summary output from the orchestrator) and provides explicit trainer modules for a multi-label Random Forest and a multi-label KNN baseline. Use them for quick prototyping and diagnostics; see [ML.md](ghostparser/ml/ML.md) for full usage and the data contract.

Run the trainers directly:

Run example:

```bash
python -m ghostparser.ml.random_forest -i results/summary_statistics.tsv -o results/ml_out
```

The KNN baseline is available as:

```bash
python -m ghostparser.ml.multi_knn -i results/summary_statistics.tsv -o results/ml_out
```

The package entrypoint `python -m ghostparser.ml` only prints those direct-run commands; it does not select a trainer itself.

**Configuration Precedence:**

`ghostparser.orchestrator` and the `ghostparser.ml` trainers each support `-c/--config-file`. When a config file is given, it supplies every setting and the other CLI flags are ignored with a warning.

---

## Configuration

See the **[Configuration Guide](CONFIG.md)** for complete details on:
- Orchestrator config keys and structure
- Machine-learning configuration sections
- JSON/YAML configuration formats
- Configuration precedence and CLI override rules


### Defaults at a Glance

Orchestrator defaults are defined in `ghostparser/orchestrator/config.py`:

**Statistical and Processing Defaults:**

- `discordant_test`: `chi-square`
- `tree_height_calculation_strategy`: `AVG`
- `p_value_correction`: `bfn`
- `alpha_dct`: `0.05`
- `alpha_ks`: `0.05`
- `alpha_perm`: `0.05`
- permutation resamples: `2500`–`25000`, `wilson` intervals

**Execution Defaults:**

- `processes`: `0` (every CPU available to the process)
- `diagnostic`: `false`
- `output_folder`: `./results`
- `overwrite`: `true`
- `min_support_value`: `0.5`
- `consolidation`: `true`
- `generate_summary_stats`: `false`
- `shape_diagnostics`: `false`
- `bootstrap`: `true`
- `bootstrap_options.iterations`: `100`
- `seed`: unset (drawn per run and reported)
- `bootstrap_options.diagnostic`: `false`
- `bootstrap_options.summary_only`: `false`

---

## Handled Errors

This section summarizes user-facing errors and validation failures that GhostParser modules can raise or report during execution.

### Orchestrator (`ghostparser.orchestrator`)

- `Error: Species tree file not found: ...` / `Error: Gene trees file not found: ...`
   The run exits early when required input files are missing.
- `✗ Error: Triplet filter file not found: ...` / `✗ Error: Species filter file not found: ...`
   The path given to `--triplet-filter` or `--species-filter` does not exist.
- `✗ Error: Species filter names N ingroup taxa; at least 3 are needed to form a triplet`
   After skipping names not found among the ingroup taxa, the species filter leaves too few to combine.
- `Config fields triplet_filter and species_filter cannot both be set`
   The two filters are alternatives; name the triplets or the species.
- `Config field species_rename_map: ...`
   The rename map is read when the config resolves; a missing or malformed file, a label mapped twice, two labels sharing a display name, or a display name holding a tab, line break, comma, semicolon or equals sign is a config error.
- `✗ Error processing species tree: ...`
   Species-tree cleaning/parsing failed (typically malformed Newick, missing taxa, or filtering issues).
- `✗ Error: Could not root the species tree ...`
   The outgroups do not root the species tree: none of them is in it, every taxon in it is an outgroup, or the outgroups branch off at more than one point. In the last case the message lists the groups the other taxa fall into, largest first, so the groups that are outgroups can be added to `--outgroup`.
- `✗ Error generating triplets: ...`
   Triplet-generation setup failed (for example rooting/pruning/mapping failures).
- `✗ Error processing gene trees: ...`
   Gene-tree cleaning/rooting stage failed before extraction.
- `✗ Error in fused extraction/inference stage: ...`
   An exception occurred while streaming triplet extraction and inference; the appended message is the originating error.

### Config Loading (`ghostparser.config` and the module config loaders)

- `Config file not found: ...`
   The config path does not exist.
- `YAML support requires PyYAML to be installed`
   YAML config was provided but `PyYAML` is unavailable.
- `Config file must be .json, .yaml, or .yml`
   Unsupported config extension.
- `Config root must be a key/value object`
   Top-level config payload is not a mapping.
- `Missing required config field: ...`
   A required field (for example the species or gene tree path) is absent or empty.
- `Missing required config field: outgroup`
   No usable outgroup taxa were provided.
- `Config field ... must be a non-empty string when provided`
   Optional string/path fields were passed as empty or the wrong type.
- `Config field ... must be an integer >= 0`
   Non-negative integer settings (for example `processes`) are invalid.
- `Config field ... must be a boolean when provided` / `Config field ... must be a numeric value`
   Boolean/float-style fields were provided with incompatible types.
- `Config field ... must be one of: ...`
   Choice-constrained fields (discordant test, tree-height strategy, p-value correction) contain unsupported values.
- `Config field overwrite must be a boolean when provided` / `Config field no_overwrite must be a boolean when provided`
   The overwrite flags were given non-boolean values.
- `Config field bootstrap_options.* ...`
   Bootstrap options failed validation (`iterations >= 1`, boolean diagnostic/summary flags).

### Tree Preprocessing (`ghostparser.orchestrator.trees`)

- `Tree file not found: ...`
   Input tree file path is missing.
- `Invalid Newick format in ...`
   Tree parsing failed, or the file is empty/structurally unusable.

### Triplet Topology Utilities (`ghostparser.triplet_utils`)

- `Triplet tree must contain exactly 3 terminal taxa`
   A triplet tree has missing/extra terminal taxa labels.
- `Could not determine rooted sister pair for triplet tree`
   Rooted triplet is unresolved/ambiguous (commonly polytomy or ambiguous rooting).
- `Tree taxa do not match provided ABC triplet`
   Topology classification was requested with an incompatible ABC taxon mapping.

### Triplet Inference (`ghostparser.orchestrator.inference`)

- `Unsupported discordant test method: ...` / `Unsupported tree height calculation strategy: ...` / `Unsupported p-value correction method: ...`
   A selected method is outside the supported choices; the message lists the valid ones.
- `Invalid species topology: ...` / `Resolved topology roles require a valid species topology`
   Internal topology state is inconsistent with the supported canonical topologies.
- `Triplet tree must contain exactly 3 terminal taxa`
   Per-tree triplet metrics require exactly three labeled leaves.
- `species_triplet is required for tree height strategies A, B, and C`
   Taxon-specific height strategies were requested without ABC triplet labels.
- `Selected taxon ... not found in triplet tree`
   The A/B/C-selected taxon is absent from the observed triplet tree.
- `Could not determine sister-pair MRCA for triplet tree`
   Sister-pair branch metrics could not be resolved for a triplet.
- `Invalid dis1_topology '...'. Expected 'BC' or 'AC'.`
   Inference text generation got an invalid discordant topology label.
- `Unsupported summary statistic name: ...`
   An unsupported statistic was requested during summary-metric computation.
- `Bootstrap payload is not JSON-serializable: ...`
   Bootstrap diagnostic output could not be converted to TSV-safe JSON.

### ML Dataset and Feature Validation (`ghostparser.ml.ml_utils`)

- `Input TSV is missing a header row` / `Input TSV contains no data rows`
   Training input file is structurally incomplete.
- `Invalid classes label: ... expected a 6-character 0/1 bitstring`
   Class label format is invalid for multi-label decoding.
- `Missing required target column: ...`
   Target column is absent from the TSV header.
- `No feature columns found after excluding the target column`
   Input has no usable predictors.
- `Missing feature value for '...'`
   At least one required feature cell is empty.
- `Non-numeric feature value for '...': ... string-valued columns must have at most 7 distinct values`
   Categorical feature expansion exceeded supported cardinality cap.
- `Cannot run stratified cross-validation because at least one class has fewer than 2 samples`
   Rare-class policy `error` blocks CV when class support is too low.

### ML Config Validation (`ghostparser.ml.config` and `ghostparser.ml.random_forest`)

- `Config file not found: ...` / `YAML support requires PyYAML to be installed`
   ML config file is missing or YAML dependency is unavailable.
- `Config file must be .json, .yaml, or .yml` / `Config root must be a key/value object`
   Unsupported config format or wrong top-level payload shape.
- `Missing required config field: ...` and typed validation messages (`must be ...`)
   Required runtime/model/evaluation fields are missing or typed incorrectly.
- `Do not place model hyperparameters or evaluation reporting controls at the top level ...`
   ML config shape is invalid; expected nested `model` and `evaluation` sections.
- `Place 'seed' and 'n_jobs' at the top level, not under 'model'. ...`
   Runtime controls were placed in the wrong config section.
- `Config field model.max_features got 'None'. Valid values are ...` (likewise `model.class_weight`)
   The bare word `None` in YAML is a string, not null. Omit the key (null is the default) or write `null` / `~`.
- `ConfigError(str(exc))` wrapping `rows_to_matrix` failures
   Random-forest CLI forwards TSV/feature validation failures as config errors.

### Hyperparameter Tuning (`ghostparser.ml.hyper_tune`)

- `hyperparameter_tuning.use_wandb is true but the 'wandb' package could not be imported ...`
   W&B logging was requested without the optional dependency; install `.[wandb]` or set `use_wandb: false`.
- `Config field hyperparameter_tuning.wandb_detailed_payloads requires hyperparameter_tuning.use_wandb: true`
   Detailed W&B payloads were requested for a run that is not logging to W&B.
- `Failed to initialize Weights & Biases for hyperparameter tuning ...`
   W&B initialization failed (authentication/network/mode setup).
- `Missing required config field: hyperparameter_tuning`
   Tuning config block is mandatory.
- `Unexpected top-level keys in hyperparameter_tuning config: ...`
   Config includes unsupported top-level keys.
- `Do not place ... at the top level in hyperparameter_tuning configs ...`
   Runtime/search/model/evaluation fields are in invalid sections.
- `Config field hyperparameter_tuning.search_space must be a mapping/object ...`
   Search-space block has invalid type.
- `Search space lists must contain at least one value`
   A tunable parameter list is empty.
- `Config field search_space must define at least one parameter`
   No search dimensions were provided.
- `Do not place runtime fields inside hyperparameter_tuning.search_space ...`
   Runtime keys were incorrectly placed inside search-space.
- `Unsupported search_space keys for ...`
   Search parameters do not match the selected model.
- `Cross-validation results do not include objective metric: ...`
   Objective metric name is not present in candidate CV output.
- `Unsupported tuning model: ...`
   Tuner received a model name outside supported model backends.
- `Hyperparameter tuning requires a feasible cross-validation split ...`
   CV fold setup is infeasible for class distribution/policy.
- `Search space expands to ... exceeds max_candidates=...`
   Grid size exceeded configured candidate cap.
- `No candidates were evaluated`
   Search loop finished without any valid candidate evaluation.

---

## Testing

Run all tests:

```bash
pytest
```

Every test carries a category marker, so a change to one part of the code can be
checked with the tests that cover it: `pytest -m core` for the statistics and
decisions, `-m config` for config loading and validation, `-m output` for the
files and columns written, `-m integration` for the entry points run end to end,
`-m parity` for the cached geometry against the DendroPy reference. Markers
combine (`-m "core and not integration"`).

See [TESTS.md](tests/TESTS.md) for the criteria, the test map, and per-test documentation.

---

## For Maintainers

### Triggering a Release Build

Release builds are triggered by pushing a version tag that matches `v*.*.*`.

```bash
# 1) Ensure main is up to date
git checkout main
git pull origin main

# 2) Create an annotated release tag
git tag -a v0.1.1 -m "Release v0.1.1"

# 3) Push the tag to trigger GitHub Actions release workflow
git push origin v0.1.1
```

Optional cleanup if you tagged the wrong commit:

```bash
# 1) Delete the local tag
git tag -d v0.1.1

# 2) Delete the remote tag
git push origin --delete v0.1.1

# 3) Manually delete the GitHub Release (if created)
# Go to: https://github.com/asif256000/ghostparser/releases
# Find the v0.1.1 release and click "Delete"
# OR use GitHub CLI:
gh release delete v0.1.1
```

Note: Deleting the tag removes the git version reference, but the GitHub Release and uploaded assets remain unless explicitly deleted. For a full cleanup, delete both the tag and the release.
