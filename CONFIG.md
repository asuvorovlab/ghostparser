# Configuration Guide

This guide is organized around the main pipeline entry point, `ghostparser.orchestrator`, and then the two submodules (`tree_parser`, `triplet_processor`).
Only `ghostparser.orchestrator` supports config files.

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

- `species_tree_path` -> `/home/user/project/../data/species.tree` -> `/home/user/data/species.tree`
- `gene_trees_path` -> `/home/user/data/genes.tree`
- `output_folder` -> `/scratch/results`

---

## Orchestrator (Primary Module)

`ghostparser.orchestrator` is the end-to-end entrypoint. It runs tree preprocessing/triplet extraction and then triplet inference in one pipeline.

### Run With a Config File

Supported config file formats:

- `.json`
- `.yaml`
- `.yml`

Run with config mode:

```bash
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml
```

If `-c/--config-file` and other CLI arguments are provided together, config mode is used and other CLI arguments are ignored.

### Orchestrator Config Keys and CLI Options

#### Required keys

##### `species_tree_path`

- Type: string path
- Parallel CLI: `--species-tree-path` (alias: `-st`)
- Description: path to the species tree file in Newick format.

##### `gene_trees_path`

- Type: string path
- Parallel CLI: `--gene-trees-path` (alias: `-gt`)
- Description: path to the gene trees file in Newick format.

##### `outgroups` or `outgroup`

- Type:
  - `outgroups`: list of strings
  - `outgroup`: string
- Parallel CLI: `--outgroups` (alias: `-og`)
- Description: outgroup taxa used for rooting/pruning.
- Accepted forms:

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

#### Optional keys

##### `triplet_filter`

- Type: string path
- Parallel CLI: `--triplet-filter`
- Description: optional path to a triplet filter file (comma-separated taxa per line).

##### `output_folder`

- Type: string path
- Parallel CLI: `--output-folder`
- Description: output directory for pipeline artifacts.
- Default: `./results` from the current working directory.

##### `triplet_output_format`

- Type: string
- Parallel CLI: `--triplet-output-format`
- Allowed values: `parquet` (default), `txt`
- Description: extraction output format from the tree parser stage.
- Value meaning:
  - `parquet`: writes partitioned parquet output optimized for downstream processing and larger datasets.
  - `txt`: writes plain-text triplet blocks for manual inspection and debugging.

##### `parquet_partitions`

- Type: integer >= 0
- Parallel CLI: `--parquet-partitions`
- Description: number of hash partitions used when triplet output format is parquet.
- Default: `128`.

##### `parquet_compression`

- Type: string
- Parallel CLI: `--parquet-compression`
- Allowed values: `zstd` (default), `snappy`, `gzip`, `brotli`, `none`
- Description: parquet compression codec.
- Value meaning:
  - `zstd`: best general-purpose balance of compression ratio and speed.
  - `snappy`: faster compression/decompression with larger output files.
  - `gzip`: higher compression ratio with slower runtime.
  - `brotli`: high compression ratio, typically slower than `zstd` for this workflow.
  - `none`: no compression, largest files and fastest write path.

##### `processes`

- Type: integer >= 0
- Parallel CLI: `--processes`
- Description: worker count for extraction and inference.
- Notes:
  - `0` uses all available CPU cores.
  - `1` uses single-worker execution.

##### `generate_summary_stats`

- Type: boolean
- Parallel CLI: `--generate-summary-stats`
- Description: writes `summary_statistics.tsv` when enabled.
- Default: `false`.

##### `min_support_value`

- Type: number
- Parallel CLI: `--min-support-value`
- Description: support threshold for species and gene tree cleaning.
- Default: `0.5`.

##### `discordant_test`

- Type: string
- Parallel CLI: `--discordant-test`
- Allowed values: `chi-square` (default), `z-test`
- Description: discordant count test used in triplet inference.
- Value meaning:
  - `chi-square`: Pearson chi-square test on discordant topology counts.
  - `z-test`: two-proportion z-test on discordant topology frequencies.

##### `summary_statistic`

- Type: string
- Parallel CLI: `--summary-statistic`
- Allowed values: `median` (default), `mean`, `mode`
- Description: statistic used for con/dis1 distributions after KS testing.
- Value meaning:
  - `median`: robust to outliers; preferred for skewed distributions.
  - `mean`: arithmetic average; sensitive to outliers.
  - `mode`: most frequent value after binning heights to 3 decimal places.

##### `stats_backend`

- Type: string
- Parallel CLI: `--stats-backend`
- Allowed values: `standard` (default), `custom`
- Description: backend used for DCT and KS calculations.
- Value meaning:
  - `standard`: SciPy/statsmodels implementations for statistical tests.
  - `custom`: GhostParser's internal implementations (experimental).

##### `tree_height_calculation_strategy`

- Type: string
- Parallel CLI: `--tree-height-calculation-strategy`
- Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`
- Description: strategy used to compute tree-height statistics.
- Value meaning:
  - `AVG`: mean root-to-tip distance across taxa A, B, and C.
  - `A`/`B`/`C`: root-to-tip distance for the selected taxon only.
  - `SIS`: sister-pair distance metric from the inferred sister taxa.
  - `INT`: internal branch distance from sister-pair MRCA to triplet root.

##### `p_value_correction`

- Type: string
- Parallel CLI: `--p-value-correction`
- Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`
- Description: multiple-testing correction for DCT and KS p-values.
- Value meaning:
  - `no`: no correction.
  - `bfn`: Bonferroni single-step family-wise error control.
  - `holm`: Holm step-down family-wise error control.
  - `fdr_bh`: Benjamini-Hochberg false discovery rate control.
  - `fdr_by`: Benjamini-Yekutieli FDR control for dependent tests.
  - `fdr_tsbh`: two-stage Benjamini-Hochberg FDR control.

##### `alpha_dct`

- Type: number
- Parallel CLI: `--alpha-dct`
- Description: p-value threshold for DCT.
- Default: `0.05`.

##### `alpha_ks`

- Type: number
- Parallel CLI: `--alpha-ks`
- Description: p-value threshold for KS.
- Default: `0.05`.

##### `bootstrap`

- Type: boolean
- Parallel CLI: `--no-bootstrap` (disable switch)
- Description: enables bootstrap sampling-with-replacement when true.
- Default: `true`.

##### `bootstrap_options`

- Type: object
- Parallel CLI:
  - `iterations` -> `--bootstrap-iterations`
  - `seed` -> `--bootstrap-seed`
  - `debug_mode` -> `--bootstrap-debug-mode`
  - `summary_only` -> `--bootstrap-summary-only`
- Description: nested bootstrap runtime options.
- Supported keys:
  - `iterations` (integer >= 1), default `100`
  - `seed` (integer, optional)
  - `debug_mode` (boolean), default `false`
  - `summary_only` (boolean), default `false`

##### `bootstrap_iterations`

- Type: integer >= 1
- Parallel CLI: `--bootstrap-iterations`
- Description: flat key alternative to `bootstrap_options.iterations`.

##### `bootstrap_seed`

- Type: integer
- Parallel CLI: `--bootstrap-seed`
- Description: flat key alternative to `bootstrap_options.seed`.

##### `bootstrap_debug_mode`

- Type: boolean
- Parallel CLI: `--bootstrap-debug-mode`
- Description: flat key alternative to `bootstrap_options.debug_mode`.

##### `bootstrap_summary_only`

- Type: boolean
- Parallel CLI: `--bootstrap-summary-only`
- Description: flat key alternative to `bootstrap_options.summary_only`.

### Orchestrator CLI Example

```bash
python -m ghostparser.orchestrator \
  --species-tree-path data/asuv21/species.tree \
  --gene-trees-path data/asuv21/gene_trees.tree \
  --outgroups Ephemera_danica,Isonychia_kiangsinensis \
  --triplet-filter data/asuv21/asuv_all/triplet_filter.txt \
  --output-folder results \
  --triplet-output-format parquet \
  --parquet-partitions 128 \
  --parquet-compression zstd \
  --processes 8 \
  --generate-summary-stats \
  --min-support-value 0.6 \
  --discordant-test z-test \
  --summary-statistic mean \
  --stats-backend custom \
  --tree-height-calculation-strategy B \
  --p-value-correction fdr_bh \
  --alpha-dct 0.02 \
  --alpha-ks 0.1 \
  --bootstrap-iterations 250 \
  --bootstrap-seed 42 \
  --bootstrap-debug-mode
```

### Sample Config Files

- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`
- `sample_configs/orchestrator_full.json`

## Tree Parser (CLI Submodule)

`ghostparser.tree_parser` performs cleaning and triplet extraction only.

### CLI Options

#### Required

##### `--species-tree-path` (`-st`)

- Parallel orchestrator key/CLI: `species_tree_path` / `--species-tree-path`
- Description: species tree input path.

##### `--gene-trees-path` (`-gt`)

- Parallel orchestrator key/CLI: `gene_trees_path` / `--gene-trees-path`
- Description: gene trees input path.

##### `--outgroups` (`-og`)

- Parallel orchestrator key/CLI: `outgroups` or `outgroup` / `--outgroups`
- Description: comma-separated outgroup taxa.

#### Optional

##### `--triplet-filter`

- Parallel orchestrator key/CLI: `triplet_filter` / `--triplet-filter`
- Description: optional triplet filter path.

##### `--output-folder`

- Parallel orchestrator key/CLI: `output_folder` / `--output-folder`
- Description: output directory for extracted artifacts.

##### `--triplet-output-format`

- Parallel orchestrator key/CLI: `triplet_output_format` / `--triplet-output-format`
- Allowed values: `parquet` (default), `txt`
- Value meaning:
  - `parquet`: partitioned columnar output for scalable downstream processing.
  - `txt`: plain-text output for quick inspection.

##### `--parquet-partitions`

- Parallel orchestrator key/CLI: `parquet_partitions` / `--parquet-partitions`
- Description: parquet partition count.

##### `--parquet-compression`

- Parallel orchestrator key/CLI: `parquet_compression` / `--parquet-compression`
- Allowed values: `zstd` (default), `snappy`, `gzip`, `brotli`, `none`
- Value meaning: same codec semantics as orchestrator `parquet_compression`.

##### `--min-support-value`

- Parallel orchestrator key/CLI: `min_support_value` / `--min-support-value`
- Description: support threshold for tree cleaning.

##### `--processes`

- Parallel orchestrator key/CLI: `processes` / `--processes`
- Description: worker count for extraction (`0` uses all cores).

##### `--no-multiprocessing`

- Parallel orchestrator key/CLI: no direct config key
- Description: disables multiprocessing in tree parser execution.

### Tree Parser CLI Example

```bash
python -m ghostparser.tree_parser \
  --species-tree-path data/asuv21/species.tree \
  --gene-trees-path data/asuv21/gene_trees.tree \
  --outgroups Ephemera_danica,Isonychia_kiangsinensis \
  --triplet-filter data/asuv21/asuv_all/triplet_filter.txt \
  --output-folder results \
  --triplet-output-format parquet \
  --parquet-partitions 128 \
  --parquet-compression zstd \
  --min-support-value 0.6 \
  --processes 8
```

## Triplet Processor (CLI Submodule)

`ghostparser.triplet_processor` runs inference on precomputed triplet extraction output.

### CLI Options

#### Required

##### `--input-path`

- Parallel orchestrator key/CLI: no direct key (orchestrator wires this stage internally)
- Description: path to triplet extraction output (`.parquet` dataset or `.txt`).

#### Optional

##### `--input-format`

- Parallel orchestrator key/CLI: no direct key
- Allowed values: `auto` (default), `txt`, `parquet`
- Description: input parser mode.
- Value meaning:
  - `auto`: infer format from input path.
  - `txt`: force text parser.
  - `parquet`: force parquet parser.

##### `--output-path`

- Parallel orchestrator key/CLI: no direct key
- Description: output TSV path.

##### `--stats-output`

- Parallel orchestrator key/CLI: no direct key
- Description: optional JSON statistics output path.

##### `--alpha-dct`

- Parallel orchestrator key/CLI: `alpha_dct` / `--alpha-dct`
- Description: DCT significance threshold.

##### `--alpha-ks`

- Parallel orchestrator key/CLI: `alpha_ks` / `--alpha-ks`
- Description: KS significance threshold.

##### `--discordant-test`

- Parallel orchestrator key/CLI: `discordant_test` / `--discordant-test`
- Allowed values: `chi-square` (default), `z-test`
- Value meaning: same method semantics as orchestrator `discordant_test`.

##### `--summary-statistic`

- Parallel orchestrator key/CLI: `summary_statistic` / `--summary-statistic`
- Allowed values: `median` (default), `mean`, `mode`
- Value meaning: same statistic semantics as orchestrator `summary_statistic`.

##### `--stats-backend`

- Parallel orchestrator key/CLI: `stats_backend` / `--stats-backend`
- Allowed values: `standard` (default), `custom`
- Value meaning: same backend semantics as orchestrator `stats_backend`.

##### `--tree-height-calculation-strategy`

- Parallel orchestrator key/CLI: `tree_height_calculation_strategy` / `--tree-height-calculation-strategy`
- Allowed values: `AVG` (default), `A`, `B`, `C`, `SIS`, `INT`
- Value meaning: same strategy semantics as orchestrator `tree_height_calculation_strategy`.

##### `--p-value-correction`

- Parallel orchestrator key/CLI: `p_value_correction` / `--p-value-correction`
- Allowed values: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`
- Value meaning: same correction-method semantics as orchestrator `p_value_correction`.

##### `--no-bootstrap`

- Parallel orchestrator key/CLI: `bootstrap` / `--no-bootstrap`
- Description: disable bootstrap sampling-with-replacement.

##### `--bootstrap-iterations`

- Parallel orchestrator key/CLI: `bootstrap_options.iterations` / `--bootstrap-iterations`
- Description: iterations per triplet.

##### `--bootstrap-seed`

- Parallel orchestrator key/CLI: `bootstrap_options.seed` / `--bootstrap-seed`
- Description: optional random seed.

##### `--bootstrap-debug-mode`

- Parallel orchestrator key/CLI: `bootstrap_options.debug_mode` / `--bootstrap-debug-mode`
- Description: enable debug bootstrap columns.

##### `--bootstrap-summary-only`

- Parallel orchestrator key/CLI: `bootstrap_options.summary_only` / `--bootstrap-summary-only`
- Description: compact bootstrap debug summaries.

##### `--processes`

- Parallel orchestrator key/CLI: `processes` / `--processes`
- Description: worker count for triplet analysis (`0` uses all cores).

##### `--generate-summary-stats`

- Parallel orchestrator key/CLI: `generate_summary_stats` / `--generate-summary-stats`
- Description: write `summary_statistics.tsv`.

##### `--no-multiprocessing`

- Parallel orchestrator key/CLI: no direct config key
- Description: disable multiprocessing in triplet processor execution.

### Triplet Processor CLI Example

```bash
python -m ghostparser.triplet_processor \
  --input-path results/unique_triplets_gene_trees.parquet \
  --input-format auto \
  --output-path results/triplet_introgression_results.tsv \
  --stats-output results/triplet_introgression_results.json \
  --alpha-dct 0.05 \
  --alpha-ks 0.05 \
  --discordant-test chi-square \
  --summary-statistic median \
  --stats-backend standard \
  --tree-height-calculation-strategy AVG \
  --p-value-correction no \
  --bootstrap-iterations 100 \
  --processes 8
```

---

## Notes

- Path values are resolved at runtime to absolute paths.
- Config-file mode (`-c/--config-file`) is available in orchestrator only.
- Use `--processes 1` to run single-worker mode.
