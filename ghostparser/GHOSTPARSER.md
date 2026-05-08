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

`triplet_processor` consumes triplet datasets from either:

- `unique_triplets_gene_trees.txt` (text sections)
- `unique_triplets_gene_trees.parquet/` (partitioned parquet dataset)

and applies the GhostParser decision pipeline:

1. Count concordant and discordant topology frequencies, set concordant to the species-tree topology, and assign discordant1/discordant2 by discordant counts (higher count -> discordant1; ties keep fixed discordant order).
2. Compute `H(T)` using a configurable tree-height strategy (`AVG` default, taxon-specific `A|B|C`, sister-distance `SIS`, or internal-branch `INT`).
3. Run discordant count test (configurable: Pearson chi-square or two-proportion z-test, alpha `alpha_dct`, default `0.05`).
4. If significant, run two-sample KS tree-height test (alpha `alpha_ks`, default `0.05`).
5. Apply selected multiple-testing correction across triplets for DCT and KS p-values (`no` default; also `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`).
6. If significant, compare selected summary statistics (median by default; mean and binned mode optional) to infer outflow vs ghost introgression.
7. When bootstrap is enabled, resample valid per-triplet observations with replacement and aggregate per-iteration classifications into bootstrap support values.

This is a configurable pipeline: users can choose supported statistical methods and thresholds while preserving the same core stage order.

Text-mode triplet sections in `unique_triplets_gene_trees.txt` are written as:

- `A,B,C<TAB>gene_tree_count<TAB>species_triplet_tree_newick<TAB>[A=speciesA,B=speciesB,C=speciesC]<TAB>AB:n/concordant,BC:n/discordant1|discordant2,AC:n/discordant2|discordant1`

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

- `unique_triplets_gene_trees.txt` (default)
- `unique_triplets_gene_trees.parquet/` (when parquet mode is enabled)
- `orchestrator_triplet_results.tsv`
- `metrics.txt`
- `introgression_combined.png`
- `introgression_matrix_inflow_outflow.tsv`
- `introgression_ghost_target_strength.tsv`
- `introgression_taxa_order.tsv`

Consolidation behavior:

- Consolidation is enabled by default in orchestrator runtime.
- The combined plot shows the inflow/outflow heatmap and the ghost target-strength bar chart side by side, sharing a single colorbar.
- Average bootstrap values written to the TSV and shown in the plot use population-level denominators:
    - Sampled introgression average for a directed pair (source → target): `sum(bootstrap weights) / count(all triplets containing both source and target)`.
    - Ghost introgression average for a target taxon: `sum(bootstrap weights) / count(all triplets containing that taxon)`.
- The species tree topology strip is drawn on top of the heatmap. Outgroup taxa passed via the `outgroups` parameter are excluded from all plots and TSVs.
- Consolidation can be disabled through orchestrator config/CLI when visualization artifacts are not needed.

For full details on the consolidation module, see the [Introgression Mapper Module](#introgression-mapper-module-ghostparserintrogression_mapper) section below.

## Introgression Mapper Module (`ghostparser.introgression_mapper`)

### Role

`ghostparser.introgression_mapper` consumes per-triplet pipeline results and produces a single combined visualization and companion TSV artifacts representing introgression signal across the ingroup taxa. It is called automatically by the orchestrator consolidation stage; see the [Orchestrator Module Guide](ORCHESTRATOR.md#additional-outputs-from-orchestrator-run) for how it is wired into the pipeline.

### `generate_introgression_maps(results, species_tree_path, output_dir, plot_taxa=None, outgroups=None)`

Generates the combined consolidation figure and tabular outputs.

**Arguments:**

- `results`: Iterable of `TripletPipelineResult`-like objects (or dict rows from a TSV).
- `species_tree_path`: Path to the processed species tree used for taxon ordering.
- `output_dir`: Directory to write all output files.
- `plot_taxa`: Optional list of taxa to retain in the plot; defaults to full ingroup.
- `outgroups`: Optional list of taxon names to exclude from all plots and TSVs (e.g. taxa used for rooting).

**Outputs:**

- `introgression_combined.png` — combined figure with directed inflow/outflow heatmap and ghost target-strength bar chart.
- `introgression_matrix_inflow_outflow.tsv` — target × source matrix of average bootstrap support values.
- `introgression_ghost_target_strength.tsv` — per-taxon average ghost bootstrap support.
- `introgression_taxa_order.tsv` — ordered taxa list matching the plot axes.

Returns an `IntrogressionMapArtifacts` dataclass with `plot_path`, TSV paths, `taxa_count`, `non_ghost_edge_count`, and `ghost_target_count`.

### Bootstrap averaging denominators

Average bootstrap support values are computed using population-level co-occurrence denominators:

- **Sampled introgression** for a directed pair (source → target):

  `avg = sum(bootstrap_value for classified triplets) / count(all triplets containing both source and target)`

  The denominator counts every triplet in which both taxa co-appear, regardless of whether that triplet was classified as introgression. For `n` ingroup taxa, a directed pair appears in exactly `n − 2` triplets.

- **Ghost introgression** for a target taxon:

  `avg = sum(bootstrap_value for ghost-classified triplets) / count(all triplets containing that taxon)`

  The denominator counts every triplet in which the target taxon appears in any position, regardless of classification. For `n` ingroup taxa, a single taxon appears in `(n−1)C2` triplets.

This normalizes signal strength by the total number of opportunities at which the event could have been detected, giving a population-level confidence estimate that accounts for the full co-occurrence space.

### Plot layout

The combined figure uses a three-row layout above the data panels:

1. **Species tree strip** (top row) — topology-only tree with leaf labels suppressed.
2. **Source taxon label strip** (middle row) — a dedicated thin row containing the source-taxon names, rotated 90°, aligned to heatmap column centres. Row height is computed from the rendered pixel-width of the longest label so labels are never clipped. Shown for datasets up to 120 taxa; suppressed beyond that.
3. **Data panels** (bottom row, left to right):
    - **Inflow/outflow heatmap** — rows are target taxa, columns are source taxa, coloured by average bootstrap support.
    - **Target label panel** — centred target taxon names aligned pixel-exactly to heatmap rows.
    - **Ghost bar chart** — horizontal bars per target taxon showing average ghost bootstrap support.
    - **Shared colorbar** — single colorbar covering both the heatmap and bar chart.

All three rows share `hspace=0` so they appear flush. Figure and panel widths scale dynamically with taxon count and rendered label widths.

### CLI usage

```bash
python -m ghostparser.introgression_mapper \
    -r orchestrator_triplet_results.tsv \
    -st processed_species.tree \
    -o output_dir/ \
    -og OutGroup1,OutGroup2
```

**Required arguments:**

- `-r`, `--results-tsv`: Path to `orchestrator_triplet_results.tsv`.
- `-st`, `--species-tree-path`: Path to the processed species tree (Newick).
- `-o`, `--output-dir`: Directory to write output plots and TSVs.

**Optional arguments:**

- `-og`, `--outgroups`: Comma-separated outgroup taxon names to exclude from plots.

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

### `run_triplet_pipeline(species_triplet, triplet_gene_trees, alpha_dct=0.05, alpha_ks=0.05, discordant_test='chi-square', summary_statistic='median', stats_backend='standard', tree_height_calculation_strategy='AVG', bootstrap=True, bootstrap_options=None)`

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

- input triplets are emitted as `(A,B,C)` where `A` and `B` are species sisters and `C` is the non-sister
- concordant is always `AB`
- `n_dis1` is whichever discordant topology has higher count
- `n_dis2` is the other discordant topology
- ties between discordants keep fixed discordant ordering (BC before AC)
- `most_frequent_matches_concordant`: `True` when concordant frequency is not lower than either discordant frequency

Output includes `species_tree` as topology-only Newick for the triplet.

Bootstrap behavior in `run_triplet_pipeline`:

- The standard non-bootstrap pipeline is always evaluated once per triplet.
- If bootstrap is enabled, observations are sampled with replacement for each iteration.
- Iterations that cannot compute required metrics are counted as `unresolved`.
- `bootstrap_value` is reported for the final `classification` value.

Possible `classification` values:

- `no_introgression`
- `outflow_introgression`
- `inflow_introgression`
- `ghost_introgression`
- `unresolved`

### `parse_triplet_gene_trees_file(filepath)`

Parses either text sections or parquet datasets into a dictionary.

For text input (`*.txt`), strict section schema is:

- header: `A,B,C<TAB>count<TAB>species_tree_newick<TAB>[A=...,B=...,C=...]<TAB>AB:x/concordant,BC:y/discordant1|discordant2,AC:z/discordant2|discordant1`
- one required blank line after header
- zero or more Newick gene-tree lines (each ending with `;`)
- sections separated by `============================================================`

Validation is strict: malformed sections, duplicate triplets, or header/tree-count mismatches raise `ValueError`.

- triplet -> `{count, species_tree, gene_trees, label_map, header_topology_counts, header_dis1_topology}`

For parquet input (`*.parquet` dataset directory), required layout is:

- `triplets/*.parquet` with one row per triplet (`A`, `B`, `C`, `count`, topology-count metadata)
- `observations/partition_id=*/part-*.parquet` with one row per extracted subtree and cached metrics (`h_a`, `h_b`, `h_c`, `h_avg`, `h_int`, `h_sis`)

### `analyze_triplet_gene_tree_file(filepath, alpha_dct=0.05, alpha_ks=0.05, discordant_test='chi-square', summary_statistic='median', stats_backend='standard', p_value_correction='no', bootstrap=True, bootstrap_options=None, rng=None, use_multiprocessing=True, processes=None)`

Runs the pipeline for all triplets in an input file with configurable discordant test and summary statistic.

### `write_pipeline_results(results, output_filepath, dct_method='chi-square', summary_statistic='median', p_value_correction='no', bootstrap=True)`

Writes per-triplet results to a TSV file with counts, DCT/KS statistics, raw and corrected p-values (`dct_p_value`, `ks_p_value`, plus dynamic corrected columns like `dct_p_val_bfn_corr`), dynamic summary columns (`median_con`/`median_dis`, `mean_con`/`mean_dis`, or `mode_con`/`mode_dis`), final classification, and an `inference` column describing direction with species names.

When `bootstrap=True`, output also includes:

Bootstrap payload columns are serialized as strict JSON; non-serializable payloads raise `ValueError`.

- `bootstrap_value`
- `all_bootstrap`

When bootstrap debug mode is enabled, output also includes:

- `bootstrap_dct_stats`
- `bootstrap_dct_p_value`
- `bootstrap_ks_stats`
- `bootstrap_ks_p_value`
- `bootstrap_con_<mean|median|mode>`
- `bootstrap_dis_<mean|median|mode>`
- `bootstrap_gene_tree_heights`

Bootstrap payload fields are serialized as JSON strings by default.

### CLI usage

```bash
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt

# or parquet dataset mode
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.parquet --input-format parquet
```

Optional arguments:

- `--output-path`: output TSV path (default: `triplet_introgression_results.tsv` next to input)
- `--stats-output`: JSON output path for full per-triplet statistics
- `--alpha-dct`: DCT threshold (default: `0.05`)
- `--alpha-ks`: KS threshold (default: `0.05`)
- `--discordant-test`: `chi-square` (default) or `z-test`
- `--summary-statistic`: `median` (default), `mean`, or `mode`
- `--stats-backend`: `standard` (default) or `custom`
- `--tree-height-calculation-strategy`: `AVG` (default), `A`, `B`, `C`, `SIS`, or `INT`
- `--p-value-correction`: `no` (default), `bfn`, `holm`, `fdr_bh`, `fdr_by`, or `fdr_tsbh`
- `--no-bootstrap`: disable bootstrap sampling-with-replacement
- `--bootstrap-iterations`: number of iterations (default: `100`)
- `--bootstrap-seed`: optional reproducibility seed
- `--bootstrap-debug-mode`: enable detailed bootstrap metric output columns
- `--bootstrap-summary-only`: with debug mode, write compact summaries instead of full per-iteration lists
- `--processes`: worker count for triplet inference (`0` = all cores)
- `--generate-summary-stats`: write `summary_statistics.tsv`
- `--no-multiprocessing`: disable multiprocessing for triplet inference

These CLI parameters can be used for focused inference runs, debugging, and testing.

### Tree parser CLI usage

Quick CLI mode (required inputs only):

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

`tree_parser` supports required arguments `-st/--species-tree-path`, `-gt/--gene-trees-path`, and `-og/--outgroups`, plus optional CLI parameters for focused extraction runs and debugging.

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
