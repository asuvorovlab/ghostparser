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