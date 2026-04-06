# Changelog

## Unreleased

- Added bootstrap sampling-with-replacement controls for orchestrator and triplet processor configs:
	- `bootstrap` (default `false`)
	- `bootstrap_options.iterations` (default `100`)
	- `bootstrap_options.seed` (optional reproducibility)
	- `bootstrap_options.summary_only` (default `true`)
- Bootstrap processing reuses per-triplet serialized gene-tree observations and runs per-iteration reanalysis over sampled observations.
- Iterations with incomplete required metrics are classified as `unresolved` and do not stop processing.
- Final classification uses bootstrap majority class; ties are reported as comma-joined class labels.
- Added bootstrap output columns to TSV output when bootstrap is enabled:
	- `bootstrap_value`
	- `bootstrap_classification`
	- `all_bootstrap`
	- `bootstrap_dct_stats`
	- `bootstrap_dct_p_value`
	- `bootstrap_ks_stats`
	- `bootstrap_ks_p_value`
	- `bootstrap_con_<mean|median|mode>`
	- `bootstrap_dis_<mean|median|mode>`
	- `bootstrap_gene_tree_heights`
- Bootstrap payload columns are serialized as JSON strings by default with compact key:value fallback.

## v0.0.3 (alpha)

- Added `--p-value-correction` option for all p-values used in summary statistics across triplets.
- Default p-value correction: `no`.
- Supported correction options: `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.
- Output records include both original and corrected p-values, with dynamic corrected-column names (for example `dct_p_val_fdr_bh_corr`).
- Inference uses corrected p-values; original p-values are retained for reporting.
- Randomized parity tests validate custom correction implementations against the standard backend.

## v0.0.4 (alpha) - April 1, 2026

- Corrected versioning of releases for consistency between setup.py and GitHub tags.