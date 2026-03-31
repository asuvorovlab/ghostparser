# Helper Scripts

The scripts in this folder are lightweight utilities for repetitive data-prep and analysis support tasks that are separate from the main GhostParser pipeline.

## Current Helper Scripts

- `consolidate_tabular_files.py`
    - Recursively finds tabular files by subfolder + name pattern + extension.
    - Consolidates all rows into one output file.
    - Adds a source-folder column to track where each row came from.
    - Example usage:

      ```bash
      python3 scripts/consolidate_tabular_files.py \
        --input-folder /path/to/root \
        --match-subfolder fdr_corrected \
        --file-starts-with orchestrator_triplet_results \
        --file-ends-with "" \
        --extension tsv \
        --output-folder /path/to/output \
        --output-file-name consolidated_orchestrator_triplet_results.tsv
      ```

- `profile_stats_methods.py`
    - Benchmarks GhostParser custom statistical helpers vs SciPy-backed versions.
    - Reports per-call timing summaries for chi-square and KS helper methods.
    - Example usage:

      ```bash
      python3 scripts/profile_stats_methods.py
      ```

- `plot_triplet_tree_heights.py`
    - Roots species and gene trees using outgroup logic consistent with `tree_parser`.
    - Computes triplet tree-height arrays for `concordant`, `discordant1`, and `discordant2` topologies.
    - Writes a histogram+KDE overlay plot and a text dump (arrays + exact topology counts + metadata).
    - Example usage:

      ```bash
      python3 scripts/plot_triplet_tree_heights.py \
        --species-tree-path /path/to/species.tree \
        --gene-trees-path /path/to/genes.tre \
        --outgroup O \
        --output-dir /path/to/output
      ```

For full CLI options and argument descriptions, run `--help` for any helper script:

```bash
python3 scripts/consolidate_tabular_files.py --help
python3 scripts/profile_stats_methods.py --help
python3 scripts/plot_triplet_tree_heights.py --help
```
