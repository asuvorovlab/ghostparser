# Changelog


## v0.1.2 - Unreleased

- Rewrote `ghostparser.orchestrator` as a package running one fused streaming pass: triplet extraction and inference share a single loop, with no intermediate `unique_triplets_gene_trees` file and no Newick round-trip between them. This removes the memory blowup on large gene-tree sets; a 3000-gene-tree / 35-triplet serial run went from roughly 44s to 25s. Bootstrap values stay deterministic under a fixed `bootstrap_seed`, so every parallelization mode agrees exactly.
- Replaced the `tree_parser` and `triplet_processor` submodules and their standalone CLIs with internal modules of the orchestrator package. `python -m ghostparser.orchestrator` is the single analysis entry point, alongside the unchanged `ghostparser.introgression_mapper` CLI.
- Added `--parallelization-mode {auto,taxon,gene}`: `taxon` dispatches triplet chunks across workers, `gene` parallelizes subtree extraction within one triplet, and `auto` picks from input size (`gene` below 15 ingroup taxa or above 3500 gene trees, else `taxon`).
- Added `--preflight-data-check` / `preflight_data_check` (default `false`), which runs only the structural check on the trees and triplets, writes `preflight_data_check.txt` to the output folder, and exits without any analysis. It collects every failure rather than raising on the first, grouped by category with counts and examples naming the offending gene tree and triplet. Replaces the removed `scripts/preflight_triplet_sanity_check.py`.
- Added overwrite control across the orchestrator, mapper, and ML trainers (`--no-overwrite` / `overwrite: false`): disabled overwrite writes to an auto-suffixed sibling directory instead of reusing the configured one. Mapper cleanup was hardened so it cannot remove species-tree inputs when paths overlap.
- Added the `ghostparser.ml` subpackage: `random_forest` and `multi_knn` trainers for multi-label prediction over `summary_statistics.tsv`, plus the `hyper_tune` grid/random hyperparameter search with optional Weights & Biases logging (`wandb_detailed_payloads`). Each run writes a pickled model, ranked `feature_importances.tsv`, `predictions.tsv` with per-row bit flags and `matched_label_count`, `*_metrics.txt` opening with a `Hyperparameters:` block (estimator settings, `test_size`, requested and effective `cv_folds`, `rare_class_policy`, `target_column`), `*_overall_metrics.json` carrying the same hyperparameters plus timings and a consolidated `dataset_summary`, and `cividis` confusion matrices exported as both text blocks and plots with percentages.
- Replaced the final inference step: the direction of the concordant-versus-discordant1 tree-height difference now comes from an adaptive Welch-studentized permutation test (`ghostparser/orchestrator/permutation.py`) rather than a bare comparison of summary statistics, so the call carries a p-value and a confidence statement instead of following any numerical difference however small. Resampling continues until a confidence interval around the p-value excludes `alpha_perm` or the budget is spent, and samples too small or too degenerate to support the test are reported as `ambiguous` with a reason rather than given a direction. See [ORCHESTRATOR.md](ghostparser/orchestrator/ORCHESTRATOR.md#the-statistical-tests) for the method and its citations.
- Added `alpha_perm` (default `0.05`, CLI `--alpha-perm`), the `permutation_test` toggle (default `true`, CLI `--no-permutation-test`, falling back to a median sign comparison), and the config-file-only `permutation_options` block (`min_resamples`, `max_resamples`, `ci_method`). Replaced the `unresolved` classification with `ambiguous`, which now also covers the case where the KS test separates the height distributions but the permutation test finds no directional evidence.
- Changed the results TSV columns: `summary_con`/`summary_dis` are replaced by the `perm_*` family (statistic, raw and corrected one-tailed p-values, two-tailed, resample count, convergence, `permutation_consistency_flag`, `perm_null_skewed`, `perm_note`, `perm_decision`). No per-group mean or median columns are written there at all — the direction comes from the p-values, and descriptive per-group statistics belong to `summary_statistics.tsv`. Removed the `summary_statistic` key and its `--summary-statistic` flag, which no longer selected anything.
- Simplified the orchestrator's outgroup configuration to the single key `outgroup` (CLI `-og/--outgroup`), replacing the `outgroups`/`outgroup` pair. It accepts a single label, a comma-separated string, or a list, and list entries may themselves be comma-separated.
- Changed one orchestrator default: `p_value_correction` from `no` to `bfn`.
- Changed the consolidation figure to use `cividis` throughout. Ghost bar length alone now carries the bootstrap value, while bar colour encodes co-occurrence — the colormap's high end (yellow) for ghost-only targets and its low end (dark blue) for targets that also have sampled introgression — recorded in a new `has_sampled_introgression` column in `introgression_ghost_target_strength.tsv`. Heatmap cells with no introgression edge are left unpainted so sparse signal stays legible, with a note stating they are off the colour scale. Axis species names are italicized.
- Removed the custom statistical backend — the discordant count test, KS test, and p-value corrections always use SciPy/statsmodels, and `stats_backend` is gone — along with the parquet/intermediate-file options (`triplet_output_format`, `input_format`, `parquet_partitions`, `parquet_compression`).
- Removed `pyarrow`, which only served the deleted parquet path, added the missing `numpy` entry to `requirements.txt`, and bumped the packaged version to `0.1.2`.
- Fixed the `discordant1_*` and `discordant2_*` columns of `summary_statistics.tsv` naming the wrong gene trees. They were hardwired to the canonical `BC|A` and `AC|B` topologies while the `dis1_topology` column and every statistical test use whichever discordant topology is more frequent, so on any triplet where `AC|B` outnumbered `BC|A` a `discordant1_*` feature described a different set of gene trees than its own row claimed — which also mislabelled the features `ghostparser.ml` trains on. The roles are now resolved per triplet from the observed counts.
- Fixed the consolidation stage deleting a run's own outputs: `generate_introgression_maps` gained a `reset_output_dir` flag, and the orchestrator writes artifacts into a dedicated `consolidation/` subfolder. Previously the results TSV, processed trees, and open `metrics.txt` were silently wiped, or the run failed with `Directory not empty` on network filesystems.
- Fixed clipped text in the consolidation figure: each panel is now at least as wide as its own title, and the figure grid is pinned to reserved margins so a row's height ratio holds in inches. Previously matplotlib's default margins absorbed roughly a fifth of every row, pushing the rotated source labels under the heatmap.
- Restructured configuration into a shared trunk plus per-module ownership. `ghostparser/config.py` keeps only the helpers that behave identically for every consumer (460 lines down to 188); `ghostparser/orchestrator/config.py` and `ghostparser/ml/config.py` own their defaults, choices, and validators, so their settings can diverge. Both CLIs share `cli_config.resolve_cli_or_config_args` for config-file-wins precedence.
- Added `tests/orchestrator/test_orchestrator_permutation.py`, covering the permutation test with parity against `scipy.stats.permutation_test`, exhaustive-enumeration validation of the vectorized sampler, randomized-input invariants, and a type-I error calibration.
- Reworked the test suite so expectations derive from definitions rather than from a second implementation: topology counts read off the fixture Newick by hand, tree heights recomputed from root-to-tip geometry, and test statistics recomputed inline with SciPy/statsmodels. The custom-vs-standard backend parity tests are gone, leaving only the DendroPy-vs-BioPython triplet-collapse checks under a new `parity` marker (replacing `backend_parity`). Added suites for the decision logic, the ML label/metric contract, and the config trunk.
- Restructured the documentation so configuration detail lives in exactly one place: `CONFIG.md` is the complete key reference, the per-module guides cover how each module works and link to it, and `GHOSTPARSER.md` is the top-level package guide. Added `tests/TEST_IO.md`, documenting each test's literal inputs and the derivation of every expected value.

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