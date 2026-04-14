# Changelog

## v0.0.3 (alpha)

- Added `--p-value-correction` option for all p-values used in summary statistics across triplets.
- Default p-value correction: `no`.
- Supported correction options: `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.
- Output records include both original and corrected p-values, with dynamic corrected-column names (for example `dct_p_val_fdr_bh_corr`).
- Inference uses corrected p-values; original p-values are retained for reporting.
- Randomized parity tests validate custom correction implementations against the standard backend.

## v0.0.4 (alpha) - April 1, 2026

- Corrected versioning of releases for consistency between setup.py and GitHub tags.

## v0.1.0 - Unreleased

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
- Orchestrator now defaults to in-memory triplet mapping/inference (no `unique_triplets_gene_trees.txt` written by default). Full mapping output is opt-in via `write_full_triplet_gene_trees_mapping` / `--write-full-triplet-gene-trees-mapping` for debugging, with potential large-file I/O overhead on large triplet sets.
- Added CPU-process timing alongside wall-clock timing in logs and `metrics.txt` for species processing, gene processing, triplet extraction, inference, and total runtime.