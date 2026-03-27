# Orchestrator Module Guide

This document describes the end-to-end orchestrator in [ghostparser/orchestrator.py](https://github.com/asif256000/ghostparser/blob/main/ghostparser/orchestrator.py):

1. tree preprocessing/triplet extraction (from `tree_parser`)
2. per-triplet introgression inference (from `triplet_processor`)
3. final TSV reporting (`orchestrator_triplet_results.tsv`)

## CLI Input Options

Run:

```bash
python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>
```

Or with a config file:

```bash
python -m ghostparser.orchestrator -c <config.yaml>
```

Options:

- `-st`, `--species-tree-path` (required)

    - Path to species tree file (Newick).

- `-gt`, `--gene-trees-path` (required)

    - Path to gene trees file (Newick, one tree per line).

- `-og`, `--outgroups` (required)

    - Outgroup taxon name(s). Multiple outgroups are comma-separated.

- `-c`, `--config-file` (optional)

    - Path to JSON/YAML config file.
    - When provided, other CLI options are ignored and a warning is printed.

- `--triplet-filter` (optional)

    - Path to triplet file (`taxon1,taxon2,taxon3` per line). If omitted, triplets are generated from species-tree ingroup taxa.

- `--output-folder` (optional)

    - Output folder path. Default is `./results` (from the current working directory).

- `--processes` (optional)

    - Worker processes for both triplet extraction (`tree_parser`) and per-triplet inference (`triplet_processor`).
    - Defaults to `0`, which means all available CPU cores (`cpu_count()`).

- `--tree-height-calculation-strategy` (optional)

    - Tree-height strategy used by `triplet_processor`.
    - Allowed values: `AVG` (default), `A`, `B`, `C`.

- `--p-value-correction` (optional)

    - Multiple-testing correction applied across triplets for DCT and KS p-values.
    - Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.

CLI mode inputs are normalized into the same key/value payload used by config files, so defaults and validation are consistent across both modes.

## Pipeline Summary

1. Species tree cleaning/rooting
   - Uses `clean_and_save_trees`, then roots/prunes with `_root_tree_on_outgroup`.
2. Triplet generation and normalization
   - Uses filter file or combinations from ingroup taxa.
   - Converts each triplet to A/B/C orientation where A and B are species sisters (`_build_species_triplet_metadata`).
3. Gene tree cleaning/rooting
   - Uses `clean_and_save_gene_trees`.
4. Triplet extraction file generation
  - Uses `write_triplet_gene_trees_multiprocess` from `tree_parser`.
  - Writes `unique_triplets_gene_trees.txt` in the same section format as `tree_parser`.
5. Per-triplet inference from file
  - Uses `analyze_triplet_gene_tree_file` from `triplet_processor` on `unique_triplets_gene_trees.txt`.
6. Final reporting
   - Uses `write_pipeline_results` to write `orchestrator_triplet_results.tsv`.

## Final TSV Columns (How Each Is Produced)

The final file is `orchestrator_triplet_results.tsv`.

### Identity and Mapping

- `triplet`

    - Source: `row["triplet"]` from `TripletPipelineResult.to_dict()`.
    - Method: tuple `(A, B, C)` joined as `A,B,C`.

### Topology Labels

Canonical topology strings:

- `((A,B),C)`
- `((B,C),A)`
- `((A,C),B)`

- `species_tree`

    - Source: species-triplet subtree header passed into `run_triplet_pipeline`.
    - Method: exact extracted species-triplet Newick string for that triplet.

Triplet labeling is canonicalized after topology-frequency counting so that:

- concordant is always `((A,B),C)`
- discordant1 is always `((B,C),A)` (and represented by `n_dis1`)
- discordant2 is always `((A,C),B)` (and represented by `n_dis2`)

If the two discordant topologies tie in frequency, canonical ordering is kept.

### Frequency Metadata

- `most_frequent_matches_concordant`

    - Source: `run_triplet_pipeline`.
    - Method: `True` when `n_con >= n_dis1` and `n_con >= n_dis2`; otherwise `False`.

### Inference Group Counts

- `n_con`

    - Method: count for species-matching topology.

- `n_dis1`

    - Method: count for more frequent discordant topology.

- `n_dis2`

    - Method: count for less frequent discordant topology.

### DCT-like Test Outputs

- `dct_chi_stats`

    - Method: included only when `--discordant-test chi-square` is used for the run.

- `dct_z_score`

    - Method: included only when `--discordant-test z-test` is used for the run.

- `dct_p_value`

    - Method: original (uncorrected) discordant-test p-value from (`n_dis1`, `n_dis2`), using selected backend (`standard` default, `custom` optional):

        - `chi-square` (default): custom chi-square in `custom`, SciPy chi-square in `standard`
        - `z-test`: custom manual z-test in `custom`, statsmodels two-proportion z-test in `standard`

- `dct_p_val_<correction>_corr`

    - Method: corrected DCT p-value after applying selected `--p-value-correction` across triplets (`no` default; also `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`).

- `dct_significant`

    - Method: `dct_p_val_<correction>_corr <= alpha_dct` (`alpha_dct` default `0.01`).

### Tree-Height Test Outputs

Per-gene-tree height uses:

- `compute_tree_height_statistic` with selected strategy:

    - `AVG`: mean root-to-tip distance across 3 leaves
    - `A`, `B`, or `C`: root-to-tip distance of that taxon
- for `((X:b2,Y:b3):b4,Z:b1)`, `AVG` gives: $H(T) = (b1 + b2 + b3 + 2b4)/3$

- `ks_statistic`

    - Method: two-sample KS statistic between height samples of `dis1_topology` vs `con_topology`.

- `ks_p_value`

    - Method: original (uncorrected) p-value from selected backend (`standard` default, `custom` optional) for two-sided KS.

- `ks_p_val_<correction>_corr`

    - Method: corrected KS p-value after applying selected `--p-value-correction` across triplets (`no` default; also `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`).

- `ks_significant`

    - Method: `ks_p_val_<correction>_corr <= alpha_ks` (`alpha_ks` default `0.05`).

### Height Summary Values and Final Classification

- Dynamic summary columns in TSV:

    - `median_con` / `median_dis` when `summary_statistic=median`
    - `mean_con` / `mean_dis` when `summary_statistic=mean`
    - `mode_con` / `mode_dis` when `summary_statistic=mode`

`mode` uses 3-decimal binning before mode selection; if multiple modes remain, the maximum mode value is used.

- `classification`

    - Decision path:

        - `no_introgression` if DCT is not significant
        - `inflow_introgression` if DCT significant but KS not significant
        - if KS significant:

            - `outflow_introgression` if `<summary>_con > <summary>_dis`
            - `ghost_introgression` if `<summary>_con < <summary>_dis`
            - `unresolved` if equal/undefined.

### Data Coverage Tracking

- `analyzed_trees`

    - Method (orchestrator output):

        - number of extracted triplet gene trees that are actually used in per-triplet inference.

    - This is the direct denominator behind per-triplet topology counts and statistics.

## Example TSV Representation

Below is an example of how one row appears in `orchestrator_triplet_results.tsv`.

Header (truncated for readability):

```tsv
triplet	species_tree	n_con	n_dis1	n_dis2	most_frequent_matches_concordant	dct_chi_stats	dct_p_value	dct_significant	ks_statistic	ks_p_value	ks_significant	median_con	median_dis	classification	analyzed_trees
```

Example data row:

```tsv
TaxaA,TaxaB,TaxaC	((TaxaA:1,TaxaB:1):1,TaxaC:1);	8	12	4	False	8	0.001	True	0.2	0.07	False			inflow_introgression	24
```

How to read this example quickly:

- `most_frequent_matches_concordant=False` because a discordant topology frequency is higher than concordant frequency.
- `classification = inflow_introgression` because DCT is significant, but KS is not significant.

## Additional Outputs from Orchestrator Run

- `processed_<species_tree_filename>`
- `processed_<gene_trees_filename>`
- `unique_triplets_gene_trees.txt`
- `orchestrator_triplet_results.tsv`
- `metrics.txt`

## Notes on Parallelism

- Multiprocessing is applied during triplet extraction via `write_triplet_gene_trees_multiprocess`.
- Multiprocessing is also applied during inference via `triplet_processor.analyze_triplet_gene_tree_file`.
- `--processes 0` explicitly resolves to all cores via `cpu_count()`.
- Use `--processes 1` for effective single-worker execution.
