# Changelog

## v0.1.3 - Unreleased

### Cached tree geometry

- Replaced per-triplet subtree extraction with a per-tree geometry cache. Previously every (triplet, gene tree) pair copied the whole gene tree, pruned the copy down to the three taxa and measured it, so a run's extraction cost was the tree size times the number of triplets times the number of gene trees. Now each gene tree is walked once -- pre-order for a parent/edge-length array, post-order for a table of every pair's lowest common ancestor (LCA) -- and each triplet's observation is read out of those tables: the rooted topology is the pair whose LCA differs from the other two pairs' (all three coinciding is a polytomy, and the triplet is skipped in that tree as before), and the tree height is the sum of edge lengths along the short path from the triplet's own LCA to the leaves the strategy names.
    - The cache is built once in the parent and inherited copy-on-write by the workers, so no gene tree is parsed or walked more than once per run.
    - Topology comes from LCA node identity rather than any distance comparison, so a zero-length internal branch cannot be mistaken for a polytomy; heights are summed along one path rather than taken as the difference of two root depths, which avoids the cancellation that form introduces.
    - The species tree goes through the same cache, which also yields each triplet's induced subtree for the `species_tree` column, child order and root edge included.
    - Counts, topologies, test statistics, p-values and classifications are unchanged; tree heights under the `AVG` strategy can differ from the previous path by floating-point rounding only.
    - The parity tests (`pytest -m parity`) hold the new path to agreement with two independent implementations: the DendroPy subtree extraction it replaced (`trees.extract_triplet_subtree` + `inference.observation_from_subtree`, kept for that purpose and called by no run) and a BioPython one written in the test suite from the definitions alone, across tree shapes, every tree-height strategy and every triplet of a nine-taxon tree, skip decisions included.
- `preflight_data_check` replays the same path and so reports the skip decisions a run would make. Its cap is the `preflight_triplet_cap` key (`--preflight-triplet-cap`, default 15000, `0` = no cap). The unreachable categories `triplet.invalid_leaf_count`, `triplet.abc_mapping_mismatch`, `triplet.sister_mrca_missing` and `species_triplet.subtree_missing` are gone; `triplet.subtree_extraction_failed` is now `triplet.geometry_unavailable`. A `usable_pairs` line accounts for every triplet/gene-tree pair, including the pairs skipped because the gene tree lacks a taxon, which previously went uncounted. The check roots and prunes the species tree with the run's own rooting rather than on the first outgroup label present, so its ingroup and triplet set are the run's; a species tree the run could not root is reported as a check that could not run, with the run's message.

### Skipped and reused computation

- The point estimate no longer measures tests the decision cascade cannot read. The direction test is skipped below a settled gate under every correction; the tree-height test is skipped below a settled count gate under `no`/`bfn`, whose correction of a p-value depends on the triplet count alone, and is still measured for every triplet under `holm`/`fdr_bh`/`fdr_by`, which rank every triplet's value against the others'. The bootstrap iterations follow the same rule, so nothing the bootstrap reads is affected. Results, bootstrap values included, are unchanged: every supported correction is monotone, and permutation p-values are corrected within the test only. A skipped test leaves its columns empty, with `perm_note: direction_test_not_consulted` for the direction test; `metrics.txt` reports what the run skips under its correction and counts the skips. The saving is modest, because the bootstrap's permutation tests dominate a run and are not skipped.
- Added `diagnostic` (`--diagnostic`, default `false`). `true` measures all three tests for every triplet in the point estimate, filling every `ks_*` and `perm_*` column, and changes no result. It does not reach the bootstrap, whose own switch is `bootstrap_options.diagnostic`.
- The TOST equivalence step reuses the direction test's own resamples instead of drawing two further permutation passes. `perm_p_tost` values are a different Monte Carlo realization of the same quantity.
- Fixed `--no-bootstrap` running the bootstrap and only dropping its columns. It now skips the iterations, under `diagnostic: true` too.

### Parallelization

- Removed `parallelization_mode` (`--parallelization-mode`). Triplet chunks are the only unit worth distributing; `--processes` alone controls parallelism.
- Fixed `processes: 0` starting one worker per core in the machine even when the process may only use some of them -- under a container's CPU limit or a job scheduler's per-job share, which can be a small fraction of the machine. The surplus workers were time-sliced over the CPUs actually available while each held its own working memory. It now counts the CPUs the process is allowed to run on, and `metrics.txt` reports the resolved count as `Worker processes`.
- Workers extract observations in bounded batches rather than a whole chunk at once, and the per-triplet result objects are slotted, so memory no longer grows with the chunk size. Results are unchanged.
- A run uses one machine. On a cluster, a job spanning several nodes runs on one of them; no multi-node execution is implemented.

### Multiple-testing correction

- Removed `fdr_tsbh` from `p_value_correction`. It was the one method that can lower a p-value, and so the one that could not license skipping a settled gate.
- Every decision is made in the run-wide pass. The stream returns measurements only, and `_apply_triplet_result_p_value_correction` sets the corrected p-values, both significance flags, the classification, its gate and `bootstrap_value` once every triplet's p-values are in. `analyze_triplet` is removed; `analyze_triplet_from_observations` returns an undecided result.

### Result columns and reporting

- Renamed `perm_stat_ci_low`/`perm_stat_ci_high` to `bootstrap_perm_stat_ci_low`/`bootstrap_perm_stat_ci_high`: a bootstrap percentile interval on the permutation test's statistic, each bootstrap iteration contributing the studentized difference of its own resample. It is empty without a bootstrap.
- Removed the shape-diagnostic columns from `summary_statistics.tsv`, where they are undefined for groups below their observation floors and left holes the trainers reject. They remain in the results TSV as `con_*`/`dis1_*`/`dis2_*`.
- `metrics.txt` reports the seed as `Seed: <n> (configured|generated)`, the `diagnostic` setting with what the run skips, and the triplets clearing each gate.
- Renamed `bootstrap_options.debug_mode` (flat `bootstrap_debug_mode`) to `bootstrap_options.diagnostic` (`bootstrap_diagnostic`). A diagnostic bootstrap now runs the direction test in every iteration as well, and adds `bootstrap_perm_stats`, `bootstrap_perm_p_greater`, `bootstrap_perm_p_less` and `bootstrap_perm_decisions` beside the count and tree-height columns; under `summary_only` the decisions are summarized as a count per decision. The direction tests added below a failed gate draw from a separate stream, so `bootstrap_value` and the `bootstrap_perm_stat_ci_*` interval are unchanged by the switch.

### Outgroup rooting

- The species tree is rooted where the outgroups branch off, in whatever orientation the file was written. Previously it was rooted on the outgroups' most recent common ancestor as the file happened to be oriented, so an unrooted tree written with the outgroups on either side of its first node could not be rooted at all (the ancestor was the whole tree), and any taxon that fell under that ancestor was silently pruned with the outgroups. Now the outgroups must branch off the rest of the tree at a single point -- one branch, or one node inside a polytomy, parting every outgroup from every other taxon -- and nothing but the outgroups is pruned. When other taxa sit between the outgroups, the run stops and the message lists the groups those taxa fall into, largest first, so the groups that are outgroups can be added to `outgroup`. `scripts/plot_triplet_tree_heights.py` imports the package's rooting instead of carrying a copy of the old one.
- Gene trees are rooted as before, each on the first listed outgroup it carries, and the run now reports what that rule did: `metrics.txt` counts the gene trees rooted on each outgroup and the gene trees in which the outgroups present do not all lie on one side of the other taxa, where rooting on another of them would give some triplets a different shape. The preflight report carries the same counts as `rooted_on` and `outgroups_apart` lines, and both advise listing first the outgroup whose placement in the gene trees is most reliable. No result changes.

### Configuration keys

- Replaced `bootstrap_options.seed`/`bootstrap_seed` with a run-wide `seed` key (`--seed`), which drives every random draw in a run. When unset, one seed is drawn for the run and reported, so any run can be reproduced from its metrics file.
- Added `species_rename_map` (`--species-rename-map`), a two-column TSV or YAML map from tree labels to the names shown in the results TSV, `summary_statistics.tsv` and the consolidation outputs. The run itself works in the tree labels, so a display name may hold spaces, dots or quotes. The file is read when the config resolves; a missing or malformed file, a label mapped twice, two labels sharing a name, or a name holding a tab, line break, comma, semicolon or equals sign is a config error.
- Newick writers single-quote a label holding whitespace or `()[]{}':;,"\=`, so an input tree with quoted labels round-trips through the processed tree files.
- Added `species_filter` (`--species-filter`), a file of taxon names — comma-separated, any number per line — that restricts the run to every triplet among the named species. Names that are not ingroup taxa are skipped with a warning; fewer than three left stops the run. It excludes `triplet_filter`; setting both is a config error. The preflight check accepts it too and caps the triplets it generates like the full ingroup's.

### Hyperparameter tuning

- Weights & Biases logging is opt-in: `hyperparameter_tuning.use_wandb` (default `false`), with `wandb` in its own `wandb` extra (`pip install .[ml,wandb]`) and imported only when used. `wandb_detailed_payloads` requires `use_wandb: true`. With W&B on, `hyper_tune_results.json`/`.tsv`, `hyper_tune_parameter_marginals.tsv` and `predictions.tsv` are logged to the run as tables instead of written to disk, and `artifact_paths` lists only what was written. Added `sample_configs/hyperparameter_tuning_multi_knn.json`.
- `hyper_tune_results.txt` echoes the search space, tabulates the top candidates and ranks the parameters by influence. Added `hyper_tune_parameter_marginals.tsv`, `hyper_tune_search_report.png`, an `elapsed_seconds` column in `hyper_tune_results.tsv`, and `use_wandb`, `parameter_marginals`, `parameter_influence` and `artifact_paths` in `hyper_tune_results.json`.

### Trainers

- The 64-class confusion-matrix heatmap shows row-normalized fractions on a fixed 0-1 scale, labels every class on both axes, and orders the classes by number of set bits; the metrics JSON's `class_labels` and `matrix` follow the same order. The per-bit panels are titled in prose.
- `model.max_features` defaults to `null` (every feature at each split) instead of `sqrt`. `model.max_features`, `model.class_weight` and their `search_space` candidates are validated in config, with messages naming the accepted forms.
- Removed dead code (`inference.compute_tree_height_statistic`, `inference._build_triplet_np_rng`, `triplet_utils.rank_topologies_by_frequency`, `tuning_report._better`) and deduplicated the trainers' validators, confusion-matrix builder and CLI resolver into `ml/config.py` and `ml_utils`. No behavior change.

### Tests and documentation

- Every test carries a category marker -- `core`, `config`, `output`, `integration`, `parity` -- so `pytest -m` can run one part of the suite.
- Default-value pins are replaced by two invariants (CLI and config-file modes resolve identically; `orchestrator_full.yaml` names every runtime key at its default), one-key rejection tests are folded into one parametrized case per loader, assertions on report wording and on artifacts nothing writes are dropped, and the bootstrap-budget divisor pin is replaced by property tests.
- Added parity coverage for the cached geometry and tests for the rename map, the preflight pair accounting and the `diagnostic` setting, each with its derivation in [tests/TEST_IO.md](tests/TEST_IO.md).

## v0.1.2 - August 18, 2026

### Architecture

- Rewrote `ghostparser.orchestrator` as a package running one fused streaming pass: triplet extraction and inference share a single loop, with no intermediate `unique_triplets_gene_trees` file and no Newick round-trip between them. This removes the memory blowup on large gene-tree sets; a 3000-gene-tree / 35-triplet serial run went from roughly 44s to 25s.
    - Replaced the `tree_parser` and `triplet_processor` submodules and their standalone CLIs with internal modules of the package, making `python -m ghostparser.orchestrator` the single analysis entry point.
    - Moved `ghostparser/introgression_mapper.py` to `ghostparser/orchestrator/consolidation.py` and dropped its standalone CLI. Consolidation is the orchestrator's final stage and its input is the in-memory results list, so there was nothing for a command line to point at that the orchestrator did not already produce; `generate_introgression_maps` takes `TripletPipelineResult` objects only. `tests/test_introgression_mapper.py` moved to `tests/orchestrator/test_orchestrator_consolidation.py`.
    - Removed the custom statistical backend — the discordant count test, KS test, and p-value corrections always use SciPy/statsmodels, and `stats_backend` is gone — along with the parquet/intermediate-file options (`triplet_output_format`, `input_format`, `parquet_partitions`, `parquet_compression`).
- Added `--parallelization-mode {auto,taxon,gene}`: `taxon` dispatches triplet chunks across workers, `gene` parallelizes subtree extraction within one triplet, and `auto` picks from input size (`gene` below 15 ingroup taxa or above 3500 gene trees, else `taxon`). Bootstrap values stay deterministic under a fixed `bootstrap_seed`, so every mode agrees exactly.
- Restructured configuration into a shared trunk plus per-module ownership. `ghostparser/config.py` keeps only the helpers that behave identically for every consumer (460 lines down to 188); `ghostparser/orchestrator/config.py` and `ghostparser/ml/config.py` own their defaults, choices, and validators, so their settings can diverge. Both CLIs share `cli_config.resolve_cli_or_config_args` for config-file-wins precedence.

### Statistical inference

- Replaced the final inference step: the direction of the concordant-versus-discordant1 tree-height difference comes from an adaptive Welch-studentized permutation test (`ghostparser/orchestrator/permutation.py`), which fully replaces the median sign comparison rather than sitting alongside it, so the call carries a p-value and a confidence statement instead of following any numerical difference however small. Resampling continues until a confidence interval around the p-value excludes `alpha_perm` or the budget is spent, and samples too small or too degenerate to support the test are reported with a reason rather than given a direction. See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md#the-statistical-tests) for the method and its citations.
    - Added `alpha_perm` (default `0.05`, CLI `--alpha-perm`) and the config-file-only `permutation_options` block (`min_resamples`, `max_resamples`, `ci_method`).
    - Replaced the `unresolved` classification with `ambiguous`, which also covers the case where the KS test separates the height distributions but the permutation test finds no directional evidence.
- Added TOST equivalence testing, splitting the permutation outcome into `equivalent` (the mean heights were shown to differ by less than the margin) and `inconclusive` (nothing was established either way). Both still classify the triplet as `ambiguous`; the distinction says whether the data supported a claim of similarity or merely failed to support a difference. The margin is a Cohen's *d* of 0.5 — half a pooled standard deviation — tested by two one-sided permutation tests against shifted nulls and reported in `perm_p_tost`.
- Added a `decision_gate` column naming the test that produced the classification (`DCT`, `THT`, or `PERM`). It is derived in the same pass as the classification — `_classify_introgression` returns both — so the two cannot disagree. The point estimate computes all three tests for every triplet, so the `perm_*` columns are populated even where the cascade stopped earlier; previously nothing in the row distinguished a `perm_decision` that produced the classification from one that was merely recorded, which made rows such as a non-directional `perm_decision` beside `classification = inflow_introgression` read as self-contradictory.
- Changed one orchestrator default: `p_value_correction` from `no` to `bfn`.

### Multiple-testing correction

- Fixed bootstrap iterations measuring a different decision rule than the classification they support. Iterations compared *raw* DCT and KS p-values to `alpha_dct`/`alpha_ks` while the reported classification came from *corrected* ones, so `bootstrap_value` answered a question nobody asked — real rows showed `classification = inflow_introgression` alongside `all_bootstrap[inflow_introgression] = 0`. Every iteration is now judged against the same corrected thresholds as the point estimate.
    - `no` and `bfn` are applied inline during the stream; `holm`, `fdr_bh`, and `fdr_by` park their raw per-iteration p-values and are corrected across triplets afterwards, one family per iteration index.
    - `fdr_tsbh` takes the same deferred path but runs every test unconditionally, since two-stage BH can lower a p-value and so cannot license skipping a gate. The other four provably cannot, which is what lets them skip the direction test once a raw gate has failed — the run's largest cost centre. See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md#correction-inside-the-bootstrap) for the per-method argument and citations.
- Changed the results TSV to omit the four `*_corr` columns under `p_value_correction: no`. With no adjustment applied they duplicated the raw p-value sitting beside them, and a header reading `dct_p_val_no_corr` implied a correction that had not been performed. The raw p-values and the `dct_significant`/`ks_significant` flags are written regardless, so nothing the classification rests on is lost; a `no` run writes 32 columns instead of 36.
- Fixed `p_value_correction: no` in a YAML config failing validation. YAML 1.1 resolves the bare word `no` to boolean `False`, so the value never reached the validator as a string and the run aborted with a message listing `no` among the valid choices. Boolean values are now mapped back to the choice they were written as (`no`/`off`/`n`/`false` and `yes`/`on`/`y`/`true`), so the unquoted spelling works for any enumerated field. The error message also names the value it received.

### Result columns and reporting

- Changed the results TSV columns: `summary_con`/`summary_dis` are replaced by the `perm_*` family (statistic, raw and corrected one-tailed p-values, `perm_p_tost`, resample count, convergence, `perm_null_skew`, `perm_note`, `perm_decision`). No per-group mean or median columns are written there at all — the direction comes from the p-values, and descriptive per-group statistics belong to `summary_statistics.tsv`. Removed the `summary_statistic` key and its `--summary-statistic` flag, which no longer selected anything.
    - `perm_null_skew` reports the permutation null's asymmetry as the sample skewness of the null, accumulated from running power sums so batches are still discarded, for every test that resampled whatever it decided. It never drives a decision; a large magnitude says a few extreme heights in the smaller group dominate the resampling.
- Added `perm_stat_ci_low` and `perm_stat_ci_high`: a bootstrap-percentile interval on the studentized mean difference at the `1 - 2 * alpha_perm` level, which is the level at which an interval and a one-sided test at `alpha_perm` agree. Unlike the binomial intervals the stopping rule places around the p-values, this one is an interval on the effect, so it says how large the separation might plausibly be rather than only whether it cleared a threshold.
- Added `shape_diagnostics` (config-file only, default `false`), which appends fifteen columns describing the shape of each height group: a KDE mode count, a Silverman critical-bandwidth modality p-value, skewness, excess kurtosis, and a generalized-Pareto tail index over the upper decile. They are descriptive and never enter a classification.
    - The modality test is used rather than a raw mode count or a BIC-selected Gaussian mixture because both of those mistake skew for structure — on 400 unimodal draws a mode count at the default bandwidth reports 3 modes for an exponential and 4 for a Pareto, and BIC selects two or more components for every skewed unimodal family tried. Silverman's test rejects none of them, at the price of needing roughly three standard deviations of separation to detect a real mixture.
    - The tail index uses a peaks-over-threshold generalized-Pareto fit rather than the Hill estimator, which is non-negative by construction and so cannot represent a bounded tail.
    - Off by default because the modality bootstrap costs about 0.2s per group, roughly 0.6s of extra CPU per triplet. The diagnostics are measured once per triplet from the observed heights; bootstrap iterations never recompute them.
    - When `generate_summary_stats` is also set they are repeated in `summary_statistics.tsv` under that file's full-word group names (`concordant_*`/`discordant1_*`/`discordant2_*`), so they are available as features to the ML trainers.
- Fixed the `discordant1_*` and `discordant2_*` columns of `summary_statistics.tsv` naming the wrong gene trees. They were hardwired to the canonical `BC|A` and `AC|B` topologies while the `dis1_topology` column and every statistical test use whichever discordant topology is more frequent, so on any triplet where `AC|B` outnumbered `BC|A` a `discordant1_*` feature described a different set of gene trees than its own row claimed — which also mislabelled the features `ghostparser.ml` trains on. The roles are now resolved per triplet from the observed counts.
- Changed `metrics.txt` to report guarded and unconverged permutation tests as counts rather than naming every triplet. On a real run those lists ran to hundreds of lines each and buried the rest of the log, while every per-triplet detail they carried is already in the results TSV under `perm_note`, `perm_converged`, and the `perm_p_*` columns. The guarded line still breaks its count down by which guard fired, since the four mean different things.

### Consolidation figure

- Changed the figure to use `cividis` throughout. Ghost bar length alone carries the bootstrap value, while bar colour encodes co-occurrence — the colormap's high end (yellow) for ghost-only targets and its low end (dark blue) for targets that also have sampled introgression — recorded in a new `has_sampled_introgression` column in `introgression_ghost_target_strength.tsv`. Heatmap cells with no introgression edge are left unpainted so sparse signal stays legible, with a note stating they are off the colour scale. Axis species names are italicized.
- Fixed the figure dropping its source-taxon labels above 120 taxa. The labels fit comfortably well past that: at 200 taxa the column pitch is 13 pt against a 5 pt font, and at 400 it is unchanged. The real constraint was a hard-coded 40-inch cap on the figure, which only binds around 150 taxa and squeezed every row — including the label strip — once it did. The cap is now the named `MAX_FIGURE_IN` at 200 inches (Agg's 65536-pixel limit is about 436 inches at the render DPI), the labels are always drawn, and when the cap does bind the heatmap absorbs the shortfall while the tree and label strips keep the inches they were measured to need.
- Fixed clipped text: each panel is at least as wide as its own title, and the figure grid is pinned to reserved margins so a row's height ratio holds in inches. Previously matplotlib's default margins absorbed roughly a fifth of every row, pushing the rotated source labels under the heatmap.
- Fixed the consolidation stage deleting a run's own outputs by writing artifacts into a dedicated `consolidation/` subfolder. Previously the results TSV, processed trees, and open `metrics.txt` were silently wiped, or the run failed with `Directory not empty` on network filesystems.

### Configuration and CLI

- Simplified the orchestrator's outgroup configuration to the single key `outgroup` (CLI `-og/--outgroup`), replacing the `outgroups`/`outgroup` pair. It accepts a single label, a comma-separated string, or a list, and list entries may themselves be comma-separated.
- Added overwrite control across the orchestrator, consolidation, and ML trainers (`--no-overwrite` / `overwrite: false`): disabled overwrite writes to an auto-suffixed sibling directory instead of reusing the configured one.
- Added `--preflight-data-check` / `preflight_data_check` (default `false`), which runs only the structural check on the trees and triplets, writes `preflight_data_check.txt` to the output folder, and exits without any analysis. It collects every failure rather than raising on the first, grouped by category with counts and examples naming the offending gene tree and triplet. Replaces the removed `scripts/preflight_triplet_sanity_check.py`.

### Machine learning

- Added the `ghostparser.ml` subpackage: `random_forest` and `multi_knn` trainers for multi-label prediction over `summary_statistics.tsv`, plus the `hyper_tune` grid/random hyperparameter search with optional Weights & Biases logging (`wandb_detailed_payloads`). Each run writes a pickled model, ranked `feature_importances.tsv`, `predictions.tsv` with per-row bit flags and `matched_label_count`, `*_metrics.txt` opening with a `Hyperparameters:` block (estimator settings, `test_size`, requested and effective `cv_folds`, `rare_class_policy`, `target_column`), `*_overall_metrics.json` carrying the same hyperparameters plus timings and a consolidated `dataset_summary`, and `cividis` confusion matrices exported as both text blocks and plots with percentages.

### Tests and documentation

- Reworked the test suite so expectations derive from definitions rather than from a second implementation: topology counts read off the fixture Newick by hand, tree heights recomputed from root-to-tip geometry, and test statistics recomputed inline with SciPy/statsmodels. The custom-vs-standard backend parity tests are gone, leaving only the DendroPy-vs-BioPython triplet-collapse checks under a new `parity` marker (replacing `backend_parity`). Added suites for the decision logic, the ML label/metric contract, and the config trunk.
    - Added `tests/orchestrator/test_orchestrator_permutation.py`, covering the permutation test with parity against `scipy.stats.permutation_test`, exhaustive-enumeration validation of the vectorized sampler, randomized-input invariants, and a type-I error calibration.
    - Added `tests/orchestrator/test_orchestrator_shape.py`, covering the shape diagnostics with moment parity against SciPy, modality calibration on samples of known modality, tail-index recovery on known tail shapes, and the two TSV column contracts.
    - Merged near-duplicate tests into parametrized ones (the five classification-cascade cases, the two permutation guards, the two summary-statistic role cases, the two correction-monotonicity cases, the three consolidation taxa-selection cases, the two count-aggregation cases, the two output-directory cases, the four ML config default/override cases, and the hyperparameter-tuning payload flag).
    - Dropped tests that exercised no logic: `test_max_resamples_reached_is_reported`, whose only assertion sat behind `if not result.converged` so it could pass without checking anything, and `test_parse_classes_returns_binary_matrix`, which drove a one-line delegation to `ml_utils.parse_classes` already covered in full by `test_parse_classes_round_trips_bitstrings`. The KNN neighbour-cap test now asserts the clamp reaches the estimator — where an over-large `n_neighbors` actually raises — and that it does not fire when unneeded.
- Restructured the documentation so configuration detail lives in exactly one place: `CONFIG.md` is the complete key reference, the per-module guides cover how each module works and link to it, and `GHOSTPARSER.md` is the top-level package guide. Added `tests/TEST_IO.md`, documenting each test's literal inputs and the derivation of every expected value.
- Trimmed docstring prose that restated derivations already written up in ORCHESTRATOR.md, keeping the Google-style `Args:`/`Returns:`/`Raises:` sections; module docstrings are a summary plus a pointer.

### Packaging

- Removed `pyarrow`, which only served the deleted parquet path, added the missing `numpy` entry to `requirements.txt`, and bumped the packaged version to `0.1.2`.

## v0.1.1 - May 17, 2026

- Changed default `alpha_dct` from `0.01` to `0.05` and default `bootstrap` from `false` to `true` across orchestrator/runtime defaults.
- Changed CLI architecture so config-file mode is available only in `ghostparser.orchestrator`; `tree_parser` and `triplet_processor` entry points use CLI-only runtime resolution with required core inputs.
- Changed multiprocessing start-method selection to prefer safe modes (`forkserver`/`spawn`) and updated CPU timing aggregation to include explicit per-worker deltas in extraction and inference pipelines.
- Added `ghostparser.introgression_mapper` consolidation module with a standalone CLI entry point, producing a single combined figure (`introgression_combined.png`) — inflow/outflow heatmap and ghost target-strength bar chart side by side with a shared colorbar and species tree strip above. Controlled via orchestrator `consolidation` config key (default `true`; `--no-consolidation` to disable).
- Added population-level co-occurrence denominators for bootstrap averaging in introgression maps, with automatic outgroup exclusion from all plot and TSV outputs.
- Changed `all_bootstrap` TSV column format from JSON to comma-separated `classification=value` pairs.
- Updated tests and documentation.


## v0.1.0 - April 20, 2026

- Added bootstrap sampling-with-replacement controls for orchestrator and triplet processor configs:

    - `bootstrap` (default `false`)
    - `bootstrap_options.iterations` (default `100`)
    - `bootstrap_options.seed` (optional reproducibility)
    - `bootstrap_options.debug_mode` (default `false`)
    - `bootstrap_options.summary_only` (default `false`, applied when debug mode is enabled)
- Bootstrap processing reuses per-triplet serialized gene-tree observations and runs per-iteration reanalysis over sampled observations.
- Iterations with incomplete required metrics are classified as `unresolved` and do not stop processing.
- `bootstrap_value` is aligned to the final corrected `classification` class fraction.
- Added bootstrap output columns to TSV output when bootstrap is enabled:

    - `bootstrap_value`
    - `all_bootstrap`

- Added bootstrap debug output columns to TSV output when bootstrap debug mode is enabled:

    - `bootstrap_dct_stats`
    - `bootstrap_dct_p_value`
    - `bootstrap_ks_stats`
    - `bootstrap_ks_p_value`
    - `bootstrap_con_<mean|median|mode>`
    - `bootstrap_dis_<mean|median|mode>`
    - `bootstrap_gene_tree_heights`

- Bootstrap payload columns are serialized as JSON strings by default with compact key:value fallback.
- Species-tree triplet output in TSV is now topology-only Newick (branch lengths omitted).
- Added tree-height strategies `SIS` (sister-taxon distance) and `INT` (sister-MRCA to triplet-root internal branch).
- Added `dis1_topology` as a base TSV output column for all runs.
- Added CPU-process timing alongside wall-clock timing in logs and `metrics.txt` for species processing, gene processing, triplet extraction, inference, and total runtime.
- Added `summary_statistics.tsv` output with per-triplet topology summaries: 63 metric-stat columns (mean/median/mode/variance/entropy/min/max × concordant/discordant1/discordant2 × avg-tree-height/internal-branch/sister-distance), plus identity fields, topology counts, classification, and bootstrap value when enabled.
- Added parquet triplet dataset support across extraction and inference (`unique_triplets_gene_trees.parquet`), including partitioned observation storage and cached per-observation tree-height metrics.
- Updated orchestrator pipeline wiring so extraction format settings are passed end-to-end (format-aware writer selection and matching inference input format).
- Added `pyarrow>=15.0` as a required dependency for parquet IO and updated focused tests to cover parquet writing, parquet parsing/analysis, and config/runtime default handling.


## v0.0.4 (alpha) - April 1, 2026

- Corrected versioning of releases for consistency between setup.py and GitHub tags.


## v0.0.3 (alpha)

- Added `--p-value-correction` option for all p-values used in summary statistics across triplets.
- Default p-value correction: `no`.
- Supported correction options: `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.
- Output records include both original and corrected p-values, with dynamic corrected-column names (for example `dct_p_val_fdr_bh_corr`).
- Inference uses corrected p-values; original p-values are retained for reporting.
- Randomized parity tests validate custom correction implementations against the standard backend.