# Configuration Guide

This guide is organized around the main pipeline entry point, `ghostparser.orchestrator`, and then the two submodules (`tree_parser`, `triplet_processor`).

## Path Resolution

GhostParser automatically resolves all path fields (input files, output directories, filter files) following standard operating system conventions:

### Path Types

1. **Absolute paths** (start with `/`):
   ```yaml
   species_tree_path: /home/user/data/species.tree
   output_folder: /scratch/results
   ```
   - Used as-is without modification
   - Platform-independent representation

2. **Relative paths** (no leading `/`):
   ```yaml
   species_tree_path: data/species.tree
   gene_trees_path: ./genes.tree
   output_folder: results
   ```
   - Resolved from the **current working directory** where the command is executed
   - Example: If you run the command from `/home/user/project/`, then `data/species.tree` resolves to `/home/user/project/data/species.tree`

3. **User home directory** (starts with `~`):
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

## Orchestrator-First Usage

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

### Supported Formats

- `.json`
- `.yaml`
- `.yml`

---

### Keys

### Required

- `species_tree_path` (string)
  - Path to species tree file.
- `gene_trees_path` (string)
  - Path to gene trees file.
- `outgroups` (list of strings)
  - Example:
    ```yaml
    outgroups:
      - Taxon1
      - Taxon2
    ```

### Optional

- `output_folder` (string)
  - Output directory path.
  - Default: `./results` from the current working directory.
- `processes` (integer >= 0)
  - Worker count for extraction and inference.
  - `0` means all available CPU cores.
  - `1` means single-worker execution (no multiprocessing).
- `triplet_filter` (string)
  - Path to triplet filter file (comma-separated taxa per line).
- `min_support_value` (number)
  - Support filtering threshold for species and gene tree cleaning.
  - Default behavior when omitted is equivalent to `0.5`.
- `discordant_test` (string)
  - Discordant count test method used by `triplet_processor` stage.
  - Allowed values: `chi-square` (default), `z-test`.
- `summary_statistic` (string)
  - Statistic used for con/dis1 distributions after KS test.
  - Allowed values: `median` (default), `mean`, `mode`.
  - `mode` bins heights to 3 decimal places before computing the mode; ties keep the maximum mode value.
- `stats_backend` (string)
  - Statistical backend used for DCT and KS computations.
  - Allowed values: `standard` (default), `custom`.
  - `custom` uses GhostParser manual statistical implementations.
  - `standard` uses SciPy for chi-square/KS and statsmodels for two-proportion z-test.
- `tree_height_calculation_strategy` (string)
  - Tree-height statistic strategy used in the triplet processor stage.
  - Allowed values: `AVG` (default), `A`, `B`, `C`.
  - `AVG` uses mean root-to-tip distance across all three taxa.
  - `A`, `B`, and `C` use only the corresponding taxon's root-to-tip distance.
- `alpha_dct` (number)
  - P-value threshold for the discordant count test.
  - Default: `0.01`.
- `alpha_ks` (number)
  - P-value threshold for KS tree-height test.
  - Default: `0.05`.

### Sample Configs

See examples in:

- `sample_configs/orchestrator_minimal.yaml`
- `sample_configs/orchestrator_full.yaml`
- `sample_configs/orchestrator_full.json`


## Tree Parser Configuration (Submodule)

### Required

- `species_tree_path` (string)
- `gene_trees_path` (string)
- `outgroups` (list of strings)
  - Example:
    ```yaml
    outgroups:
      - Taxon1
      - Taxon2
    ```

### Optional

- `output_folder` (string)
  - Output folder relative to the input species-tree folder.
- `processes` (integer >= 0)
  - Worker count for triplet extraction (`0` = all cores).
- `triplet_filter` (string)
  - Path to optional triplet filter file.
- `min_support_value` (number)
  - Support filtering threshold for species and gene tree cleaning.
  - Default: `0.5`.
- `no_multiprocessing` (boolean)
  - `true` forces single-worker extraction.

### Sample Configs

- `sample_configs/tree_parser_minimal.yaml`
- `sample_configs/tree_parser_full.yaml`


## Triplet Processor Configuration (Submodule)

### Required

- `input_path` (string)
  - Path to `unique_triplets_gene_trees.txt`.

### Optional

- `output_path` (string)
  - Output TSV path.
  - Default: `<input_dir>/triplet_introgression_results.tsv`.
- `stats_output` (string)
  - Optional JSON statistics output path.
  - Default: same path as output TSV with `.json` extension.
- `alpha_dct` (number)
  - Default: `0.01`.
- `alpha_ks` (number)
  - Default: `0.05`.
- `discordant_test` (string)
  - Allowed values: `chi-square` (default), `z-test`.
- `summary_statistic` (string)
  - Allowed values: `median` (default), `mean`, `mode`.
  - `mode` bins heights to 3 decimal places before computing the mode; ties keep the maximum mode value.
- `stats_backend` (string)
  - Statistical backend used for DCT and KS computations.
  - Allowed values: `standard` (default), `custom`.
  - `custom` uses GhostParser manual statistical implementations.
  - `standard` uses SciPy for chi-square/KS and statsmodels for two-proportion z-test.
- `tree_height_calculation_strategy` (string)
  - Allowed values: `AVG` (default), `A`, `B`, `C`.
  - `AVG` uses mean root-to-tip distance across all three taxa.
  - `A`, `B`, and `C` use only the corresponding taxon's root-to-tip distance.
- `processes` (integer >= 0)
  - Worker count for triplet inference (`0` = all cores).
- `no_multiprocessing` (boolean)
  - `true` forces single-worker analysis.

### Sample Configs

- `sample_configs/triplet_processor_minimal.yaml`
- `sample_configs/triplet_processor_full.yaml`

---

## Notes

- **Path Resolution**: All paths (absolute, relative, or `~`-prefixed) are automatically resolved to absolute paths at runtime. See the [Path Resolution](#path-resolution) section above for details.
- In config mode (`-c/--config-file`), other CLI options are ignored with a warning.
- Use `--processes 1` (or `processes: 1`) to disable multiprocessing while still using the same pipeline.
- Default output folder is `./results` relative to the current working directory when not specified.
