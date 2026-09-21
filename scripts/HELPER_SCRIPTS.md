# Helper Scripts

The scripts in this folder are lightweight utilities for repetitive data-prep and analysis support tasks that are separate from the main GhostParser orchestrator.

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

- `plot_triplet_tree_heights.py`
    - Roots species and gene trees with the orchestrator's own outgroup rooting, so the package must be importable (run from the repository root or install it).
    - Computes branch-metric arrays for `concordant`, `discordant1`, and `discordant2` topologies.
    - Supports `--branch-strategy` with: `avg_height`, `internal_branch`, `sister_distance`, `external_branch`.
    - Writes a histogram+KDE overlay plot.
    - Writes topology height arrays to `triplet_tree_heights_arrays.json`.
    - Writes metadata to `triplet_tree_heights_output.txt`, including exact topology counts and discordant mapping (`disc1_topology`, `disc2_topology`).
    - Example usage:

      ```bash
      python3 scripts/plot_triplet_tree_heights.py \
        --species-tree-path /path/to/species.tree \
        --gene-trees-path /path/to/genes.tre \
        --outgroup O \
        --branch-strategy internal_branch \
        --output-dir /path/to/output \
        --bins 100
      ```

- `run_orchestrator_parent_dir.py`
    - Iterates all immediate subdirectories under a provided parent directory.
    - Reads each folder's gene trees from `<folder>/<--gene-trees-file>` and writes its results to `<folder>/<--output-subdir>` (default `bfn_correction_results`).
    - The species tree, the gene-tree file name and the outgroups are required; nothing about the data is assumed.
    - Invokes `ghostparser.orchestrator` for each folder via a generated JSON config.
    - Example usage:

      ```bash
      python3 scripts/run_orchestrator_parent_dir.py \
        --parent-dir /path/to/simulations \
        --species-tree /path/to/species.tree \
        --outgroups OutGroup \
        --gene-trees-file genes.tre
      ```

For full CLI options and argument descriptions, run `--help` for any helper script:

```bash
python3 scripts/consolidate_tabular_files.py --help
python3 scripts/plot_triplet_tree_heights.py --help
python3 scripts/run_orchestrator_parent_dir.py --help
```
