# Ghostparser Module Documentation

Detailed documentation for the GhostParser pipeline, centered on `ghostparser.orchestrator`, with `ghostparser.tree_parser` and `ghostparser.triplet_processor` as supporting submodules.

## Overview

GhostParser is an orchestrator-first pipeline:

1. `orchestrator` coordinates the full run and output generation.
2. `tree_parser` handles tree cleaning, outgroup rooting/pruning, and triplet extraction.
3. `triplet_processor` performs per-triplet statistical inference and classification.

The pipeline processes phylogenetic trees in Newick format, standardizes trees by removing support values, extracts triplet subtrees from gene trees, and runs the GhostParser decision path for introgression calls.

When multiple outgroup taxa are provided (comma-separated), the species tree is rooted on their most recent common
ancestor (MRCA) and the outgroup clade is pruned. Any additional taxa that fall inside the outgroup clade are
excluded from triplet generation and logged as a warning (including the full list of excluded taxa) in the metrics
file.

For very large datasets, triplet processing parallelizes over triplet chunks and streams gene trees from disk,
keeping memory bounded to the active triplet chunk.

`triplet_processor` consumes `unique_triplets_gene_trees.txt` and applies the GhostParser decision pipeline:

1. Count concordant and discordant topology frequencies, then relabel taxa so concordant is `AB|C` and discordant1 is `BC|A`.
2. Compute `H(T)` using a configurable tree-height strategy (`AVG` default, or taxon-specific `A|B|C`).
3. Run discordant count test (configurable: Pearson chi-square or two-proportion z-test, alpha `alpha_dct`, default `0.01`).
4. If significant, run two-sample KS tree-height test (alpha `alpha_ks`, default `0.05`).
5. Apply selected multiple-testing correction across triplets for DCT and KS p-values (`no` default; also `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`).
6. If significant, compare selected summary statistics (median by default; mean and binned mode optional) to infer outflow vs ghost introgression.

This is a configurable pipeline: users can choose supported statistical methods and thresholds while preserving the same core stage order.

Triplet sections in `unique_triplets_gene_trees.txt` are written as:

- `A,B,C<TAB>gene_tree_count<TAB>species_triplet_tree_newick`

where `A` and `B` are species sisters for that triplet.

## Orchestrator (Primary Entry Point)

### Role

`ghostparser.orchestrator` is the main pipeline interface. It executes preprocessing and inference in sequence and writes final orchestrator outputs.

### CLI usage

```bash
python -m ghostparser.orchestrator -st <species_tree> -gt <gene_trees> -og <outgroup>
```

Config-file mode:

```bash
python -m ghostparser.orchestrator -c <config.yaml>
```

### Why the pipeline works

- `tree_parser` produces rooted/cleaned species and gene trees plus normalized triplet blocks.
- `triplet_processor` analyzes those triplet blocks with DCT + KS + summary-comparison logic.
- `orchestrator` preserves one consistent runtime configuration path and final reporting format.

### Primary orchestrator outputs

- `unique_triplets_gene_trees.txt`
- `orchestrator_triplet_results.tsv`
- `metrics.txt`

## Triplet Processor Module (`ghostparser.triplet_processor`)

### `compute_tree_height_statistic(tree)`

Computes the triplet tree-height statistic with a selected strategy:

`H(T) = mean(root-to-tip distances)` when strategy is `AVG`.

For `((X:b2,Y:b3):b4,Z:b1)`, this is:

`H(T) = (b1 + b2 + b3 + 2*b4) / 3`.

With taxon-specific strategy `A`, `B`, or `C`, `H(T)` is the root-to-tip distance of the selected taxon.

### `classify_triplet_topology(tree, species_triplet, topology_counts, species_topology=TOPOLOGY_AB)`

Classifies rooted triplet topology relative to `(A,B,C)` into:

- `concordant`: species-matching topology
- `discordant1`: more frequent discordant topology label (based on `topology_counts`)
- `discordant2`: less frequent discordant topology label (based on `topology_counts`)
- ties between discordants are resolved deterministically by picking the first discordant topology as `discordant1`

Returns:

- `(label, most_frequent_matches_concordant)` where `most_frequent_matches_concordant` is `True` when concordant count is not lower than either discordant count.

### `run_triplet_pipeline(species_triplet, triplet_gene_trees, alpha_dct=0.01, alpha_ks=0.05, discordant_test='chi-square', summary_statistic='median', stats_backend='standard', tree_height_calculation_strategy='AVG')`

Runs full sequential GhostParser logic and returns counts, p-values, summary values, and final classification.

Discordant test options:

- `chi-square` (default): custom implementation in `custom` backend, `scipy.stats.chisquare` in `standard` backend
- `z-test`: two-proportion z-test (manual in `custom` backend; `statsmodels.stats.proportion.proportions_ztest` in `standard` backend)

Summary statistic options after KS test:

- `median` (default)
- `mean`
- `mode` (values rounded to 3 decimals before computing mode; if multiple modes, maximum is used)

Statistical backend options:

- `standard` (default): uses SciPy reference implementations for chi-square and KS, and statsmodels for two-proportion z-test
- `custom`: uses GhostParser manual chi-square, two-proportion z-test, and KS implementations

Topology convention:

- labels are canonicalized after counting so concordant is always `((A,B),C)`
- `n_dis1` always maps to `((B,C),A)`
- `n_dis2` always maps to `((A,C),B)`
- ties between discordants keep canonical ordering
- `most_frequent_matches_concordant`: `True` when concordant frequency is not lower than either discordant frequency

Output includes `species_tree` (the extracted species-tree Newick for the triplet).

Possible `classification` values:

- `no_introgression`
- `outflow_introgression`
- `inflow_introgression`
- `ghost_introgression`
- `unresolved`

### `parse_triplet_gene_trees_file(filepath)`

Parses `unique_triplets_gene_trees.txt` into a dictionary:

- triplet -> `{count, species_tree, gene_trees}`

### `analyze_triplet_gene_tree_file(filepath, alpha_dct=0.01, alpha_ks=0.05, discordant_test='chi-square', summary_statistic='median', stats_backend='standard', p_value_correction='bfn', rng=None, use_multiprocessing=True, processes=None)`

Runs the pipeline for all triplets in an input file with configurable discordant test and summary statistic.

### `write_pipeline_results(results, output_filepath)`

Writes per-triplet results to a TSV file with counts, DCT/KS statistics, raw and corrected p-values (`dct_p_value`, `ks_p_value`, plus dynamic corrected columns like `dct_p_val_bfn_corr`), dynamic summary columns (`median_con`/`median_dis`, `mean_con`/`mean_dis`, or `mode_con`/`mode_dis`), and classification.

### CLI usage

```bash
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt
```

Config-file mode:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

Optional arguments:

- `--output-path`: output TSV path (default: `triplet_introgression_results.tsv` next to input)
- `--alpha-dct`: DCT threshold (default: `0.01`)
- `--alpha-ks`: KS threshold (default: `0.05`)
- `--summary-statistic`: `median` (default), `mean`, or `mode`
- `--stats-backend`: `standard` (default) or `custom`
- `--tree-height-calculation-strategy`: `AVG` (default), `A`, `B`, or `C`
- `--p-value-correction`: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`
- `--processes`: worker count for triplet inference (`0` = all cores)
- `--no-multiprocessing`: disable multiprocessing for triplet inference

When `-c/--config-file` is provided, other CLI options are ignored with a warning.

CLI mode and config-file mode both resolve through the same normalization and default-validation path.

### Tree parser CLI usage

Quick CLI mode (required inputs only):

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

Config-file mode:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

When `-c/--config-file` is provided, other CLI options are ignored with a warning.

## Tree Parser Module (`ghostparser.tree_parser`)

### Module Functions

### Tree Reading and Validation

#### `read_tree_file(filepath)`

Read Newick trees from a file with validation.

**Arguments:**

- `filepath`: Path to the Newick tree file

**Returns:**

- List of Bio.Phylo tree objects

**Raises:**

- `FileNotFoundError`: If the file does not exist
- `ValueError`: If the file contains invalid Newick format

**Example:**

```python
from ghostparser.tree_parser import read_tree_file

trees = read_tree_file("species.tree")
print(f"Read {len(trees)} trees")
```

### Support Value Processing

#### `calculate_average_support(tree)`

Calculate the average support value for a tree from internal nodes.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The average support value across all internal nodes, or `None` if no support values found

**Example:**

```python
avg_support = calculate_average_support(tree)
if avg_support is not None:
    print(f"Average support: {avg_support:.4f}")
```

#### `remove_support_values(tree)`

Remove support values (bootstrap/posterior probabilities) from internal nodes.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The tree with support values removed from all internal nodes

**Example:**

```python
clean_tree = remove_support_values(tree)
```

### Tree Standardization

#### `standardize_tree(tree)`

Standardize a tree by removing support values while preserving branch lengths.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- The standardized tree with support values removed

**Example:**

```python
standardized = standardize_tree(tree)
```

### Formatting and Writing

#### `format_newick_with_precision(tree, decimal_places=10)`

Format a tree as Newick string with custom branch length precision.

**Arguments:**

- `tree`: A Bio.Phylo tree object
- `decimal_places`: Number of decimal places for branch lengths (default: 10)

**Returns:**

- Newick format string with specified precision and trailing zeros removed

**Example:**

```python
newick = format_newick_with_precision(tree, decimal_places=10)
print(newick)
```

#### `write_clean_trees(trees, output_filepath)`

Write cleaned trees to a file in Newick format.

**Arguments:**

- `trees`: List of Bio.Phylo tree objects
- `output_filepath`: Path where trees will be written

**Example:**

```python
write_clean_trees(cleaned_trees, "output_clean.tree")
```

#### `clean_and_save_trees(input_filepath, output_filepath, support_threshold=0.5)`

Clean trees and save to file, filtering by support threshold.

**Arguments:**

- `input_filepath`: Path to input Newick file
- `output_filepath`: Path for output file
- `support_threshold`: Minimum average support to keep trees (default: 0.5)

**Returns:**

- Tuple of (cleaned_trees, dropped_indices_dict)

**Example:**

```python
trees, dropped = clean_and_save_trees("input.tree", "output.tree")
print(f"Kept {len(trees)} trees, dropped {len(dropped)}")
```

### Taxa and Triplet Operations

#### `get_taxa_from_tree(tree)`

Extract all terminal taxa names from a phylogenetic tree.

**Arguments:**

- `tree`: A Bio.Phylo tree object

**Returns:**

- A sorted list of terminal taxa names

**Example:**

```python
taxa = get_taxa_from_tree(tree)
print(f"Found taxa: {', '.join(taxa)}")
```

#### `generate_triplets(taxa_list, outgroup)`

Generate all unique triplet combinations from taxa excluding the outgroup(s).

**Arguments:**

- `taxa_list`: List of all taxa names
- `outgroup`: Outgroup taxon or iterable of taxa to exclude (comma-separated string supported)

**Returns:**

- A list of tuples, where each tuple contains 3 taxa names

**Formula:** nC3 where n = len(taxa_list) - k (excluding k outgroup taxa)

**Example:**

```python
triplets = generate_triplets(taxa, "OutGroup1,OutGroup2")
print(f"Generated {len(triplets)} triplets")
```

#### `write_triplets_to_file(triplets, output_filepath)`

Write triplets to a file, one triplet per line.

**Arguments:**

- `triplets`: List of triplet tuples
- `output_filepath`: Path where the triplets will be written

**Example:**

```python
write_triplets_to_file(triplets, "triplets.txt")
```

### Triplet Subtree Extraction

#### `extract_triplet_subtree(tree, triplet_taxa)`

Extract a subtree containing only the specified triplet taxa.

**Arguments:**

- `tree`: A Bio.Phylo tree object
- `triplet_taxa`: List or tuple of 3 taxa names to extract

**Returns:**

- A new Bio.Phylo tree containing only the specified taxa, or `None` if not all taxa are present

**Details:**

- Creates a deep copy of the tree to avoid modifying the original
- Prunes away all taxa not in the triplet
- Recalculates branch lengths appropriately
- Branch length aggregation is handled by Bio.Phylo during pruning: when a parent ends up with a single child, Bio.Phylo collapses that edge and folds its length into the child, so root-to-leaf distances for the retained taxa match the original tree (we do not manually sum lengths)

**Example:**

```python
triplet = ("TaxaA", "TaxaB", "TaxaC")
subtree = extract_triplet_subtree(gene_tree, triplet)
if subtree:
    print(f"Extracted subtree with {len(subtree.get_terminals())} taxa")
```

#### `process_gene_trees_for_triplets(gene_trees, triplets)`

Process gene trees to extract subtrees for each triplet.

**Arguments:**

- `gene_trees`: List of Bio.Phylo tree objects
- `triplets`: List of triplet tuples

**Returns:**

- A dictionary mapping triplets to lists of subtree Newick strings

**Details:**

- For each gene tree, attempts to extract a subtree for each triplet
- Skips gene trees that don't contain all 3 taxa from a triplet
- Same gene tree can contribute to multiple triplets

**Example:**

```python
triplet_trees = process_gene_trees_for_triplets(gene_trees, triplets)
for triplet, trees in triplet_trees.items():
    print(f"{triplet}: {len(trees)} gene trees")
```

#### `write_triplet_gene_trees(triplet_gene_trees, output_filepath)`

Write triplet gene trees to a file in the specified format.

**Arguments:**

- `triplet_gene_trees`: Dictionary mapping triplets to lists of Newick strings
- `output_filepath`: Path where the results will be written

**Format:**

```
TaxonA,TaxonB,TaxonC<TAB>count

newick_tree_1
newick_tree_2
...


NextTriplet...
```

**Example:**

```python
write_triplet_gene_trees(triplet_trees, "triplet_gene_trees.txt")
```

### Utility Functions

#### `get_clean_filename(filepath)`

Generate a clean filename with \_clean prefix.

**Arguments:**

- `filepath`: The original file path

**Returns:**

- The path with \_clean inserted before the file extension

**Example:**

```python
clean_name = get_clean_filename("species.tree")
# Returns: "species_clean.tree"
```

### Command Line Interface

### Main Entry Point

#### `main()`

Main entry point for the tree parser module with argparse CLI.

**Usage:**

```bash
python -m ghostparser.tree_parser -st SPECIES_TREE -gt GENE_TREES -og OUTGROUP
```

**Arguments:**

**Required:**
- `-st, --species-tree-path`: Path to species tree file in Newick format
- `-gt, --gene-trees-path`: Path to gene trees file in Newick format
- `-og, --outgroups`: Outgroup species identifier(s), comma-separated when multiple

**Optional:**
- `--processes`: Number of worker processes for multiprocessing (only used when multiprocessing is enabled).
                    Defaults to `0` (all cores). Ignored if `--no-multiprocessing` is set.
- `--no-multiprocessing`: Disable multiprocessing for triplet extraction.
                         Processes triplets sequentially on a single worker.
                         Useful for debugging or systems with limited memory.

**Output Files:**

1. `{species_tree}_clean.tree` - Cleaned species tree
2. `{gene_trees}_clean.tree` - Cleaned gene trees
3. `{gene_trees}_unique_triplet_gene_trees.txt` - Triplet gene trees

**Example:**

```bash
python -m ghostparser.tree_parser \
   -st data/species.tree \
   -gt data/genes.tree \
   -og OutGroup
```

### Workflow Details

### Complete Processing Pipeline

1. **Initialization**

   - Parse command line arguments
   - Validate input files exist
   - Generate output filenames

2. **Species Tree Processing**

   - Read species tree
   - Calculate average support
   - Remove support values
   - Filter by support threshold
   - Write cleaned tree
   - Extract taxa names

3. **Triplet Generation**

   - Remove outgroup from taxa list
   - Generate all nC3 combinations
   - Display triplet count

4. **Gene Tree Processing**

   - Read gene trees
   - Calculate average support for each
   - Remove support values
   - Filter by support threshold
   - Write cleaned trees

5. **Triplet Extraction**

   - For each gene tree and triplet combination:
     - Check if gene tree contains all 3 taxa
     - Extract subtree if all present
     - Store subtree under triplet key
   - Write results to triplet gene trees file

6. **Statistics and Reporting**
   - Total triplets generated
   - Triplets with matching gene trees
   - Total subtrees extracted
   - Average trees per triplet

### Configuration

### Support Threshold

Default: `0.5`

Trees with average support below this threshold are filtered out. This ensures only high-quality trees are used for analysis.

### Branch Length Precision

Default: `10 decimal places`

Branch lengths are preserved with high precision and trailing zeros are removed for cleaner output.

### Dependencies

- **BioPython (>=1.79)**: For Phylo module (tree parsing and manipulation)
- **Python 3.x**: Uses f-strings, pathlib, type hints

### Error Handling

The module provides robust error handling:

- **Invalid Newick Format**: Raises `ValueError` with descriptive message
- **Missing Files**: Raises `FileNotFoundError` with file path
- **Empty Trees**: Validated during reading
- **Missing Taxa**: Returns `None` for invalid triplet extractions
- **Missing Outgroup**: Warns user and shows available taxa
