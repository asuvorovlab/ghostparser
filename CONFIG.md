# Configuration Guide

This guide is organized around the main pipeline entry point, `ghostparser.orchestrator`, and then the two submodules (`tree_parser`, `triplet_processor`).

## Path Resolution

GhostParser automatically resolves all path fields (input files, output directories, filter files) following standard operating system conventions:

### Path Types

#### 1. Absolute paths (start with `/`)

```yaml
species_tree_path: /home/user/data/species.tree
output_folder: /scratch/results
```

- Used as-is without modification
- Platform-independent representation

#### 2. Relative paths (no leading `/`)

```yaml
species_tree_path: data/species.tree
gene_trees_path: ./genes.tree
output_folder: results
```

- Resolved from the **current working directory** where the command is executed
- Example: If you run the command from `/home/user/project/`, then `data/species.tree` resolves to `/home/user/project/data/species.tree`

#### 3. User home directory (starts with `~`)

```yaml
species_tree_path: ~/data/species.tree
output_folder: ~/results
triplet_filter: ~/filters/triplets.txt
```

- `~` expands to your home directory (e.g., `/home/username/`)
- Example: `~/data/species.tree` becomes `/home/username/data/species.tree`

### Important Notes

- **Path resolution happens at runtime** when the config/CLI is parsed
- **All path types work in both CLI and config file modes**
- **Relative paths are NOT relative to the config file location** - they are relative to the directory where you execute the command
- For portability, consider using relative paths in configs and running commands from a consistent location

### Examples

**Config file at** `~/project/configs/run.yaml`:
```yaml
species_tree_path: ../data/species.tree    # Relative to execution directory, not config file
gene_trees_path: ~/data/genes.tree         # User home directory
output_folder: /scratch/results            # Absolute path
```

**Executed from** `/home/user/project/`:
```bash
python -m ghostparser.orchestrator -c configs/run.yaml
```

**Resolved paths:**

- `species_tree_path` → `/home/user/project/../data/species.tree` → `/home/user/data/species.tree`
- `gene_trees_path` → `/home/user/data/genes.tree`
- `output_folder` → `/scratch/results`

---

## Orchestrator Usage

Run orchestrator with a config file:

```bash
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

Run orchestrator with plain CLI arguments:

```bash
python -m ghostparser.orchestrator \
  -st data/asuv21/species.tree \
  -gt data/asuv21/gene_trees.tree \
  -og Ephemera_danica,Isonychia_kiangsinensis
```

If `-c/--config-file` and other CLI args are provided together, config mode is used and the other CLI args are ignored with a warning.

CLI mode internally normalizes provided flags into the same key/value configuration payload used by config files, so validation and defaults are consistent across both modes.

Submodule config-file modes:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

`triplet_processor` config-file mode:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

---

## Orchestrator Configuration (Primary)

### Consolidated CLI Example (Non-Default Values)

```bash
python -m ghostparser.orchestrator \
    --species-tree-path data/asuv21/species.tree \
    --gene-trees-path data/asuv21/gene_trees.tree \
    --outgroups Ephemera_danica,Isonychia_kiangsinensis \
    --triplet-filter data/asuv21/asuv_all/triplet_filter.txt \
    --output-folder results \
    --processes 8 \
    --min-support-value 0.6 \
    --discordant-test z-test \
    --summary-statistic mean \
    --stats-backend custom \
    --tree-height-calculation-strategy B \
    --p-value-correction fdr_bh \
    --alpha-dct 0.02 \
    --alpha-ks 0.1 \
    --bootstrap \
    --bootstrap-iterations 250 \
    --bootstrap-seed 42 \
    --bootstrap-debug-mode
```

### Supported Formats

- `.json`
- `.yaml`
- `.yml`

---

### Keys

### Required

- `species_tree_path` (string)

    - CLI flag: `--species-tree-path` (alias: `-st`)

        - Path to species tree file.

- `gene_trees_path` (string)

    - CLI flag: `--gene-trees-path` (alias: `-gt`)

        - Path to gene trees file.

- `outgroups` (list of strings) or `outgroup` (string)

    - CLI flag: `--outgroups` (alias: `-og`)

        - A single string is parsed as one taxon.
        - A comma-separated string is parsed as multiple taxa.
        - Examples:

                ```yaml
                outgroups:
                    - Taxon1
                    - Taxon2
                ```

                ```yaml
                outgroup: Taxon1
                ```

                ```yaml
                outgroup: Taxon1,Taxon2
                ```

### Optional

- `output_folder` (string)

    - CLI flag: `--output-folder`

    - Output directory path.
    - Default: `./results` from the current working directory.

- `processes` (integer >= 0)

    - CLI flag: `--processes`

    - Worker count for extraction and inference.
    - `0` means all available CPU cores.
    - `1` means single-worker execution (no multiprocessing).

- `write_full_triplet_gene_trees_mapping` (boolean)

    - CLI flag: `--write-full-triplet-gene-trees-mapping`

    - When `true`, writes the full `unique_triplets_gene_trees.txt` mapping before inference.
    - Useful for debugging and auditing extracted triplet-to-gene-tree mappings.
    - For large numbers of triplets, this file can become very large and may add noticeable I/O overhead.
    - Default: `false` (in-memory streaming mode, no intermediate full mapping file).

- `triplet_filter` (string)

    - CLI flag: `--triplet-filter`

    - Path to triplet filter file (comma-separated taxa per line).

- `min_support_value` (number)

    - CLI flag: `--min-support-value`

    - Support filtering threshold for species and gene tree cleaning.
    - Default behavior when omitted is equivalent to `0.5`.

- `discordant_test` (string)

    - CLI flag: `--discordant-test`

    - Discordant count test method used by `triplet_processor` stage.
    - Allowed values: `chi-square` (default), `z-test`.

- `summary_statistic` (string)

    - CLI flag: `--summary-statistic`

    - Statistic used for con/dis1 distributions after KS test.
    - Allowed values: `median` (default), `mean`, `mode`.
    - `mode` bins heights to 3 decimal places before computing the mode; ties keep the maximum mode value.

- `stats_backend` (string)

    - CLI flag: `--stats-backend`

    - Statistical backend used for DCT and KS computations.
    - Allowed values: `standard` (default), `custom`.
    - `custom` uses GhostParser manual statistical implementations.
    - `standard` uses SciPy for chi-square/KS and statsmodels for two-proportion z-test.

- `tree_height_calculation_strategy` (string)

    - CLI flag: `--tree-height-calculation-strategy`

    - Tree-height statistic strategy used in the triplet processor stage.
    - Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`.
    - `AVG` uses mean root-to-tip distance across all three taxa.
    - `A`, `B`, and `C` use only the corresponding taxon's root-to-tip distance.
    - `SIS` uses the pairwise distance between the two sister taxa in each gene-tree topology.
    - `INT` uses the internal branch length from the sister-pair MRCA to the triplet root.

- `p_value_correction` (string)

    - CLI flag: `--p-value-correction`

    - Multiple-testing correction applied across triplets for DCT and KS p-values.
    - Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.
    - Method details:

        - `bfn` (Bonferroni, single-step family-wise error control): for $m$ tests, corrected value is $p_i' = \min(1, m\cdot p_i)$.
        - `holm` (Holm-Bonferroni, step-down family-wise error control): sort p-values ascending and scale each by remaining hypotheses, enforcing monotonicity.
        - `fdr_bh` (single-stage FDR Benjamini-Hochberg): sort p-values, scale by rank using $p_{(i)}' = p_{(i)}\cdot m / i$, then enforce monotonicity from largest to smallest rank.
        - `fdr_by` (Benjamini-Yekutieli): BH-style scaling with an additional harmonic-factor multiplier for dependence robustness.
        - `fdr_tsbh` (two-stage BH): performs BH with an adaptive estimate of true nulls using the selected alpha.
        - `no`: no correction; corrected values are identical to original p-values.

    - `no` disables p-value correction.

- `alpha_dct` (number)

    - CLI flag: `--alpha-dct`

    - P-value threshold for the discordant count test.
    - Default: `0.01`.

- `alpha_ks` (number)

    - CLI flag: `--alpha-ks`

    - P-value threshold for KS tree-height test.
    - Default: `0.05`.

- `bootstrap` (boolean)

    - CLI flag: `--bootstrap`

    - Enables bootstrap sampling-with-replacement per triplet.
    - Default: `false`.

- `bootstrap_options` (object)

    - CLI flags for nested keys:
        - `iterations` -> `--bootstrap-iterations`
        - `seed` -> `--bootstrap-seed`
        - `debug_mode` -> `--bootstrap-debug-mode`
        - `summary_only` -> `--bootstrap-summary-only`

    - Bootstrap runtime options.
    - Supported keys:
        - `iterations` (integer >= 1): number of bootstrap iterations per triplet. Default: `100`.
        - `seed` (integer, optional): enables reproducible bootstrap sampling when provided.
        - `debug_mode` (boolean): when `true`, detailed bootstrap metric columns are written. Default: `false`.
        - `summary_only` (boolean): when `true` and debug mode is enabled, detailed metric columns store compact summaries; otherwise they store full per-iteration lists. Default: `false`.

    - When bootstrap is enabled, the final TSV includes additional columns:
        - `bootstrap_value`
        - `all_bootstrap`

    - When bootstrap debug mode is enabled, the final TSV also includes:
        - `bootstrap_dct_stats`
        - `bootstrap_dct_p_value`
        - `bootstrap_ks_stats`
        - `bootstrap_ks_p_value`
        - `bootstrap_con_<mean|median|mode>`
        - `bootstrap_dis_<mean|median|mode>`
        - `bootstrap_gene_tree_heights`

### Sample Configs

See examples in:

- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`
- `sample_configs/orchestrator_full.json`


## Tree Parser Configuration (Submodule)

### Required

- `species_tree_path` (string)

    - CLI flag: `--species-tree-path` (alias: `-st`)

- `gene_trees_path` (string)

    - CLI flag: `--gene-trees-path` (alias: `-gt`)

- `outgroups` (list of strings) or `outgroup` (string)

    - CLI flag: `--outgroups` (alias: `-og`)

        - A single string is parsed as one taxon.
        - A comma-separated string is parsed as multiple taxa.

### Optional

- `output_folder` (string)

    - CLI flag: `--output-folder`

    - Output folder relative to the input species-tree folder.

- `processes` (integer >= 0)

    - CLI flag: `--processes`

    - Worker count for triplet extraction (`0` = all cores).

- `triplet_filter` (string)

    - CLI flag: `--triplet-filter`

    - Path to optional triplet filter file.

- `min_support_value` (number)

    - CLI flag: `--min-support-value`

    - Support filtering threshold for species and gene tree cleaning.
    - Default: `0.5`.

- `no_multiprocessing` (boolean)

    - CLI flag: `--no-multiprocessing`

    - `true` forces single-worker extraction.

### Sample Configs

- `sample_configs/tree_parser_minimal.yaml`
- `sample_configs/tree_parser_full.yaml`


## Triplet Processor Configuration (Submodule)

### Required

- `input_path` (string)

    - CLI flag: `--input-path`

        - Path to `unique_triplets_gene_trees.txt`.

### Optional

- `output_path` (string)

    - CLI flag: `--output-path`

    - Output TSV path.
    - Default: `<input_dir>/triplet_introgression_results.tsv`.

- `stats_output` (string)

    - CLI flag: `--stats-output`

    - Optional JSON statistics output path.
    - Default: same path as output TSV with `.json` extension.

- `alpha_dct` (number)

    - CLI flag: `--alpha-dct`

    - Default: `0.01`.

- `alpha_ks` (number)

    - CLI flag: `--alpha-ks`

    - Default: `0.05`.

- `discordant_test` (string)

    - CLI flag: `--discordant-test`

    - Allowed values: `chi-square` (default), `z-test`.

- `summary_statistic` (string)

    - CLI flag: `--summary-statistic`

    - Allowed values: `median` (default), `mean`, `mode`.
    - `mode` bins heights to 3 decimal places before computing the mode; ties keep the maximum mode value.

- `stats_backend` (string)

    - CLI flag: `--stats-backend`

    - Statistical backend used for DCT and KS computations.
    - Allowed values: `standard` (default), `custom`.
    - `custom` uses GhostParser manual statistical implementations.
    - `standard` uses SciPy for chi-square/KS and statsmodels for two-proportion z-test.

- `tree_height_calculation_strategy` (string)

    - CLI flag: `--tree-height-calculation-strategy`

    - Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`.
    - `AVG` uses mean root-to-tip distance across all three taxa.
    - `A`, `B`, and `C` use only the corresponding taxon's root-to-tip distance.
    - `SIS` uses the pairwise distance between the two sister taxa in each gene-tree topology.
    - `INT` uses the internal branch length from the sister-pair MRCA to the triplet root.

- `p_value_correction` (string)

    - CLI flag: `--p-value-correction`

    - Multiple-testing correction applied across triplets for DCT and KS p-values.
    - Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.

- `processes` (integer >= 0)

    - CLI flag: `--processes`

    - Worker count for triplet inference (`0` = all cores).

- `no_multiprocessing` (boolean)

    - CLI flag: `--no-multiprocessing`

    - `true` forces single-worker analysis.

- `bootstrap` (boolean)

    - CLI flag: `--bootstrap`

    - Enables bootstrap sampling-with-replacement per triplet.
    - Default: `false`.

- `bootstrap_options` (object)

    - CLI flags for nested keys:
        - `iterations` -> `--bootstrap-iterations`
        - `seed` -> `--bootstrap-seed`
        - `debug_mode` -> `--bootstrap-debug-mode`
        - `summary_only` -> `--bootstrap-summary-only`

    - Bootstrap runtime options.
    - Supported keys:
        - `iterations` (integer >= 1): number of bootstrap iterations per triplet. Default: `100`.
        - `seed` (integer, optional): enables reproducible bootstrap sampling when provided.
        - `debug_mode` (boolean): when `true`, detailed bootstrap metric columns are written. Default: `false`.
        - `summary_only` (boolean): when `true` and debug mode is enabled, detailed metric columns store compact summaries; otherwise they store full per-iteration lists. Default: `false`.

    - When bootstrap is enabled, the TSV writer appends:
        - `bootstrap_value`
        - `all_bootstrap`

    - When bootstrap debug mode is enabled, the TSV writer also appends:
        - `bootstrap_dct_stats`
        - `bootstrap_dct_p_value`
        - `bootstrap_ks_stats`
        - `bootstrap_ks_p_value`
        - `bootstrap_con_<mean|median|mode>`
        - `bootstrap_dis_<mean|median|mode>`
        - `bootstrap_gene_tree_heights`

### Sample Configs

- `sample_configs/triplet_processor_minimal.yaml`
- `sample_configs/triplet_processor_full.yaml`

---

## Notes

- **Path Resolution**: All paths (absolute, relative, or `~`-prefixed) are automatically resolved to absolute paths at runtime. See the [Path Resolution](#path-resolution) section above for details.
- In config mode (`-c/--config-file`), other CLI options are ignored with a warning.
- Use `--processes 1` (or `processes: 1`) to disable multiprocessing while still using the same pipeline.
- Default output folder is `./results` relative to the current working directory when not specified.
