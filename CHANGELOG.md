# Changelog


## v0.1.2 - Unreleased

- Added a new self-contained streaming subpackage `ghostparser.pipeline` (`python -m ghostparser.pipeline`) that fuses triplet subtree extraction and per-triplet inference into a single pass, so the intermediate `unique_triplets_gene_trees` dataset is never materialized or reloaded — removing the two-stage memory blowup on large gene-tree sets. It ports (copies and cleans) the preprocessing and inference logic it needs from `tree_parser`/`triplet_processor` rather than importing them, and reuses the shared `config`, `triplet_utils`, and `introgression_mapper` foundation. The v1 surface is CLI-only (no config-file mode) with a `--parallelization-mode {auto,taxon,gene}` selector; it pins the remaining knobs to shared defaults, skips summary-statistics gathering, and writes results to `pipeline_triplet_results.tsv`. Per-triplet results match the orchestrator (verified by parity tests under `tests/pipeline/`). Consolidation writes its artifacts into a dedicated `consolidation/` subfolder of the output directory so its output-directory reset cannot delete the run's results TSV, processed trees, or the open `metrics.txt` (or fail with `Directory not empty` on network filesystems). The streaming engine computes each triplet's `(topology, tree-height)` observation directly from the extracted subtree object instead of serializing to Newick and reparsing (removing the dominant redundant DendroPy parse), and vectorizes bootstrap resampling with NumPy; on a 3000-gene-tree / 35-triplet benchmark this cut a serial run from roughly 44s to 25s. Because observations now come from the unrounded subtree objects, tree heights match the orchestrator's parquet observation path, and the NumPy-based bootstrap values are statistically equivalent to (not bit-identical with) the orchestrator's `random.Random` bootstrap while remaining deterministic under a fixed seed.
- Fixed the consolidation stage deleting a run's own outputs: `generate_introgression_maps` gained a `reset_output_dir` flag (default `True` for the standalone CLI). The orchestrator now calls it with `reset_output_dir=False`, so consolidation writes into the already-prepared run directory instead of resetting it — the results TSV, processed trees, intermediate parquet, and open `metrics.txt` are preserved (previously they were silently wiped on local filesystems, or the run failed with `Directory not empty` on network filesystems).
- Added overwrite control for result directories across orchestrator, tree parser, introgression mapper, and ML trainers (`--no-overwrite` / `overwrite: false`). Default behavior reuses/overwrites configured output directories; disabled overwrite writes to an auto-suffixed sibling directory using a one-scan smallest-missing-suffix allocator. Mapper cleanup was also hardened to avoid removing species-tree inputs when paths overlap.
- Added new machine-learning subpackage `ghostparser.ml`, including `random_forest` and `multi_knn` trainers for multi-label prediction workflows, plus a standalone hyperparameter tuner (`ghostparser.ml.hyper_tune`). This release also adds stricter `hyperparameter_tuning` config validation and defaults/observability updates (WandB run naming, optional `wandb_detailed_payloads`, and KNN tuning example).
- Updated ML artifacts and diagnostics: text metrics are now `*_metrics.txt`, JSON remains `*_overall_metrics.json`, `predictions.tsv` includes `matched_label_count`, metrics include timing + consolidated `dataset_summary`, confusion matrices include percentages and export as both text blocks and plots (`*_confusion_matrices.png`, `*_confusion_matrix_64_classes.png`).
- Updated consolidation plotting, tests, and documentation: species names on consolidation axes are now italicized, and docs/tests were refreshed for overwrite behavior, tuner constraints, split strategy behavior, and feature-importance guidance.


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