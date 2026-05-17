# Changelog


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