# Helper Scripts

The scripts in this folder are lightweight utilities for repetitive data-prep and analysis support tasks that are separate from the main GhostParser pipeline.

## Current Helper Scripts

- `consolidate_tabular_files.py`

  - Recursively finds tabular files by subfolder + name pattern + extension.
    - Consolidates all rows into one output file.
    - Adds a source-folder column to track where each row came from.

- `profile_stats_methods.py`

    - Benchmarks GhostParser custom statistical helpers vs SciPy-backed versions.
    - Reports per-call timing summaries for chi-square and KS helper methods.

- `plot_triplet_tree_heights.py`

    - Roots species and gene trees using outgroup logic consistent with `tree_parser`.
    - Computes triplet tree-height arrays for `concordant`, `discordant1`, and `discordant2` topologies.
    - Writes a histogram+KDE overlay plot and a text dump (arrays + exact topology counts + metadata).

Example:

```bash
python3 scripts/plot_triplet_tree_heights.py \
  --species-tree-path /path/to/species.tree \
  --gene-trees-path /path/to/genes.tre \
  --outgroup O \
  --output-dir /path/to/output
```

## Quick Usage

```bash
conda run -n ghostparser python scripts/consolidate_tabular_files.py \
  --input-folder /path/to/root \
  --match-subfolder fdr_corrected \
  --file-starts-with orchestrator_triplet_results \
  --file-ends-with "" \
  --extension tsv \
  --output-folder /path/to/output \
  --output-file-name consolidated_orchestrator_triplet_results.tsv
```

Use `--help` for full options:

```bash
conda run -n ghostparser python scripts/consolidate_tabular_files.py --help
```

Notes:

- Use `--match-subfolder` to restrict matching to files inside specific subfolders.
- `--file-starts-with` and `--file-ends-with` are matched against the file name stem (not full path).
- If `--output-file-name` is omitted, the script writes `consolidated<extension>` (for example, `consolidated.tsv`).
- `--source-column-name` can rename the source tracking column (default: `source_folder`).

## Profiling Helper Usage

```bash
conda run -n ghostparser python scripts/profile_stats_methods.py
```

For custom benchmark settings, run:

```bash
conda run -n ghostparser python scripts/profile_stats_methods.py --help
```
