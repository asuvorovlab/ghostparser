# Ghostparser

## Overview

**Ghostparser** is a phylogenetic introgression pipeline centered on `ghostparser.orchestrator`. The orchestrator runs two internal stages (`tree_parser` and `triplet_processor`) and produces final triplet-level inference outputs. Orchestrator supports command-line and configuration-file modes, and the two submodules support CLI parameters that are useful for focused runs, debugging, and testing.

---

## Jump to Sections:

1. [Quick Start](#quick-start)
2. [Modules](#modules)
   - [Orchestrator (Primary Pipeline)](#orchestrator-primary-pipeline)
   - [Tree Parser](#tree-parser-submodule)
   - [Triplet Processor](#triplet-processor-submodule)
3. [Orchestrator Input/Output](#orchestrator-inputoutput)
4. [Configuration](#configuration)
5. [Defaults](#defaults-at-a-glance)
6. [Handled Errors](#handled-errors)
7. [Testing](#testing)
8. [For Maintainers](#for-maintainers)

---

## Quick Start

#### Installation from Wheel (Recommended for Users)

Pre-built wheel distributions are automatically created and published when a release is tagged. Users can download and install them without building locally.

**Downloading the wheel from GitHub:**

1. Visit the [Releases page](https://github.com/asif256000/ghostparser/releases)
2. Find the latest release (e.g., `v0.1.0`)
3. Download the `.whl` file (e.g., `ghostparser-0.1.0-py3-none-any.whl`)

**Installing from wheel:**

```bash
# Install from downloaded wheel file
pip install ghostparser-0.1.0-py3-none-any.whl

# OR install directly from GitHub releases URL
pip install https://github.com/asif256000/ghostparser/releases/download/v0.1.0/ghostparser-0.1.0-py3-none-any.whl

# OR install from source distribution
pip install https://github.com/asif256000/ghostparser/releases/download/v0.1.0/ghostparser-0.1.0.tar.gz
```

**After installation, users can run from any directory:**

```bash
python -m ghostparser.orchestrator -c config.yaml

python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup

python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt
```

**Note:** The wheel installation installs the package into your Python environment, so you don't need to be in the project directory to run it. All module commands (`python -m ghostparser.*`) work from anywhere.


### Installation via pip (For Developers)

```bash
pip install -r requirements.txt
```

If you want the machine-learning baselines, install the optional ML extra instead of the core-only package:

```bash
pip install .[ml]
```

That extra pulls in the scikit-learn dependency used by `ghostparser.ml.random_forest` and `ghostparser.ml.multi_knn`.

The same extra also enables `ghostparser.ml.hyper_tune`, which reads a `hyperparameter_tuning` config section and can run either grid search or random search across the supported ML trainers.

Run the tuner directly with:

```bash
python -m ghostparser.ml.hyper_tune -c sample_configs/hyperparameter_tuning_random_forest.yaml
```

#### Poetry 2.x+ Alternative

If you prefer Poetry 2.x+ instead of plain pip:

```bash
poetry install
```

Run commands with Poetry:

```bash
poetry run pytest -q
poetry run python -m ghostparser.orchestrator -c run_config.yaml
```

---

## Modules

### Orchestrator (Primary Pipeline)

Run with config file (recommended for reproducibility):

```bash
python -m ghostparser.orchestrator -c run_config.yaml
```

Run the full pipeline via CLI flags:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

CLI mode with explicit worker count:

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --processes 4
```

#### How the Pipeline Works

- `tree_parser` standardizes species/gene trees, roots on outgroup(s), and provides triplet extraction utilities.
- `triplet_processor` applies the GhostParser statistical decision pipeline to each triplet payload.
- The orchestrator coordinates both steps and writes the final results table.

GhostParser is configurable (discordant test, backend, thresholds, summary statistic), so execution follows the same core pipeline stages while allowing controlled method choices.

#### Orchestrator Arguments

**Required:**
- `-st, --species-tree-path`
- `-gt, --gene-trees-path`
- `-og, --outgroups`

**Config mode:**
- `-c, --config-file`

**Common optional:**
- `--output-folder`
- `--no-overwrite`
- `--triplet-filter`
- `--processes`
- `--generate-summary-stats`
- `--min-support-value`
- `--discordant-test`
- `--summary-statistic`
- `--stats-backend`
- `--tree-height-calculation-strategy`
- `--p-value-correction`
- `--alpha-dct`, `--alpha-ks`
- `--no-bootstrap`
- `--bootstrap-iterations`
- `--bootstrap-seed`
- `--bootstrap-debug-mode`
- `--bootstrap-summary-only`
- `--triplet-output-format`
- `--parquet-partitions`
- `--parquet-compression`
- `--no-consolidation`

#### Primary Outputs

1. `unique_triplets_gene_trees.parquet` (default; use `--triplet-output-format txt` to write text output)
2. `orchestrator_triplet_results.tsv`
3. `summary_statistics.tsv` (only when `--generate-summary-stats` is enabled)
4. `metrics.txt`
5. `introgression_heatmap_inflow_outflow.png` and `introgression_ghost_target_strength.png` (written when consolidation is enabled)
6. `introgression_matrix_inflow_outflow.tsv`, `introgression_ghost_target_strength.tsv`, and `introgression_taxa_order.tsv`

Heatmap consolidation details:

- The inflow/outflow heatmap and ghost target-strength bar chart use raw values with a shared max/min/mid color legend per plot.
- The species tree topology is stitched onto the top and left edges of the heatmap so the source and target axes read like tree labels.
- By default, the plots use the processed species tree after outgroup pruning; when a triplet filter is supplied, the plotted tree can be pruned to the taxa represented in that filtered set.
- Consolidation is enabled by default and can be disabled with `--no-consolidation`.

### Streaming pipeline (`ghostparser.pipeline`)

`ghostparser.pipeline` fuses triplet subtree extraction and per-triplet inference into a single streaming pass, so the intermediate triplet-gene-trees dataset is never materialized. Use it for large gene-tree sets where the orchestrator's two-stage design exhausts memory.

```bash
python -m ghostparser.pipeline -st species.tree -gt genes.tree -og OutGroup
python -m ghostparser.pipeline -st species.tree -gt genes.tree -og OutGroup --processes 0 --parallelization-mode auto
python -m ghostparser.pipeline -c run_config.yaml
```

It supports both CLI flags and a JSON/YAML config file (`-c/--config-file`, config-file mode; the file wins over other CLI flags). Required inputs are `-st/--species-tree-path`, `-gt/--gene-trees-path`, and `-og/--outgroups`. Config+CLI options include `--output-folder`, `--triplet-filter`, `--no-overwrite`, `--processes`, `--parallelization-mode {auto,taxon,gene}`, `--alpha-dct`, `--alpha-ks`, `--p-value-correction`, `--summary-statistic`, `--no-consolidation`, and `--no-bootstrap`; further knobs (including `generate_summary_stats` and the bootstrap-debug options) are config-file-only. The statistical tests always use the scipy/statsmodels backend. Results are written to `pipeline_triplet_results.tsv`. See [ghostparser/pipeline/PIPELINE.md](ghostparser/pipeline/PIPELINE.md) for the full reference and design.

### Tree Parser (Submodule)

Use this module when you only want preprocessing + triplet extraction.

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

#### Core Behavior

- Removes support labels and preserves branch lengths
- If support values are present, removes trees with average support below `min_support_value`
- Roots on outgroup(s), prunes outgroup clade, and logs excluded taxa
- Writes processed trees and triplet extraction output (`unique_triplets_gene_trees.parquet` by default, `unique_triplets_gene_trees.txt` when `--triplet-output-format txt` is selected)

Useful CLI options for focused runs and debugging:

- `--triplet-filter`
- `--output-folder`
- `--no-overwrite`
- `--min-support-value`
- `--processes`
- `--no-multiprocessing`
- `--triplet-output-format`
- `--parquet-partitions`
- `--parquet-compression`

### Triplet Processor (Submodule)

Use this module when you already have triplet extraction output (`unique_triplets_gene_trees.parquet` by default or `unique_triplets_gene_trees.txt`) and only need inference.

```bash
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt
```

#### Core Behavior

- Runs DCT (`chi-square` or `z-test`)
- Runs KS tree-height test when DCT is significant
- Applies summary-statistic comparison (`median`, `mean`, or binned `mode`) for final classification
- Supports tree-height strategy selection via `tree_height_calculation_strategy`: `AVG` (default), taxon-specific `A|B|C`, sister-distance `SIS`, or internal-branch `INT`
- Optionally runs bootstrap sampling-with-replacement per triplet using reusable per-gene-tree observations

Useful CLI options for focused runs and debugging:

- `--output-path`
- `--stats-output`
- `--input-format`
- `--alpha-dct`, `--alpha-ks`
- `--discordant-test`
- `--summary-statistic`
- `--stats-backend`
- `--tree-height-calculation-strategy`
- `--p-value-correction`
- `--bootstrap-iterations`
- `--bootstrap-seed`
- `--bootstrap-debug-mode`
- `--bootstrap-summary-only`
- `--processes`
- `--generate-summary-stats`
- `--no-multiprocessing`

#### Bootstrap Behavior

- Bootstrap is enabled by default and can be disabled with `--no-bootstrap`; additional bootstrap controls are configured through orchestrator (`--bootstrap-iterations`, `--bootstrap-seed`, `--bootstrap-debug-mode`, `--bootstrap-summary-only`).
- Iterations with incomplete required metrics are counted as `unresolved` and processing continues.
- `bootstrap_value` reports the bootstrap fraction for the final `classification` value after correction.

---

---

## Orchestrator Input/Output

### Input Expectations

**Required Arguments:**

- `--species-tree-path` (alias `-st`)

   - Path to the species tree file in Newick format.

- `--gene-trees-path` (alias `-gt`)

   - Path to the gene trees file in Newick format.

- `--outgroups` (alias `-og`)

   - Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups.

**Optional Arguments:**

- `--output-folder`

   - Output folder relative to the input data folder (default: same folder as input data).

- `--triplet-filter`

   - Path to a triplet filter file (comma-separated taxa per line).
   - When provided, only those triplets are processed.
   - Triplets containing taxa missing from the species tree are skipped with a warning.

- `--processes`

   - Number of worker processes for multiprocessing.
   - Defaults to `0` (all cores).
   - Ignored if `--no-multiprocessing` is set.

- `--no-multiprocessing`

   - Disable multiprocessing.
   - Processes triplets sequentially using a single worker.
   - Useful for debugging or on systems with limited resources.

### Output Files

The orchestrator generates these output files:

1. **`processed_species.tree`** - Processed species tree with support values removed and outgroup rooting applied
2. **`processed_gene_trees.tree`** - Processed gene trees with support values removed and outgroup rooting applied
3. **`unique_triplets_gene_trees.parquet`** - Default triplet-to-gene-tree mapping output used as the input to triplet inference (`.txt` can be selected with `--triplet-output-format txt`)
4. **`metrics.txt`** - Metrics log with warnings, timings, and counts
5. **`orchestrator_triplet_results.tsv`** - Final triplet-level classification results (`no_introgression`, `outflow_introgression`, `inflow_introgression`, `ghost_introgression`, or `unresolved`)
6. **`summary_statistics.tsv`** - Optional per-triplet summary table, written only when `--generate-summary-stats` is enabled, including:
   - identity columns (`triplet`, `abc_mapping`, `species_tree`, `dis1_topology`)
   - topology counts (`n_con`, `n_dis1`, `n_dis2`)
   - 63 topology/metric summary columns (7 statistics × 3 topology classes × 3 metric types)
   - final `classification` and `bootstrap_value` (when bootstrap is enabled)

Base TSV output includes `dis1_topology` and a topology-only `species_tree` value for each triplet.
Base TSV output also includes an `inference` column with human-readable direction text using actual species names.

When bootstrap is enabled, the TSV adds:

- `bootstrap_value`
- `all_bootstrap`

When bootstrap debug mode is enabled, the TSV also adds:

- `bootstrap_dct_stats`
- `bootstrap_dct_p_value`
- `bootstrap_ks_stats`
- `bootstrap_ks_p_value`
- `bootstrap_con_<mean|median|mode>`
- `bootstrap_dis_<mean|median|mode>`
- `bootstrap_gene_tree_heights`

Bootstrap payload columns are serialized as JSON strings by default.

### Example Usage

**Species tree** (`species.tree`):

```
(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);
```

**Gene trees** (`genes.tree`):

```
((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

**Run:**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup
```

**Multiple outgroups (comma-separated):**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup1,OutGroup2
```

**With triplet filter:**

```bash
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --triplet-filter triplets.txt
```

When multiple outgroups are provided, the species tree is rooted on their most recent common ancestor (MRCA) and the outgroup clade is pruned. Any additional taxa that fall inside the outgroup clade are excluded from triplet generation and logged as a warning (including the full list of excluded taxa) in the metrics file.

**Output** (`unique_triplets_gene_trees.txt`):

```
TaxaA,TaxaB,TaxaC	2	((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);

((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.13,TaxaB:0.24):0.35,TaxaC:0.46):0.57;


TaxaA,TaxaC,TaxaD	1	((TaxaA:0.1,TaxaC:0.2):0.3,TaxaD:0.4);

((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
```

---

## Configuration

See the **[Configuration Guide](CONFIG.md)** for complete details on:
- Orchestrator-first config keys and structure
- Module-specific configuration sections
- JSON/YAML configuration formats
- Configuration precedence and CLI override rules

---

## Defaults at a Glance

Core defaults are centralized in orchestrator config/CLI normalization and in module CLI runtime defaults:

**Statistical and Processing Defaults:**

- `discordant_test`: `chi-square`
- `summary_statistic`: `median`
- `stats_backend`: `standard`
- `tree_height_calculation_strategy`: `AVG`
- `p_value_correction`: `no`
- `alpha_dct`: `0.05`
- `alpha_ks`: `0.05`

**Execution Defaults:**

- `processes`: `0` (all available CPU cores)
- `output_folder` (orchestrator/tree_parser): `./results`
- `overwrite` (orchestrator/tree_parser): `true`
- `triplet_output_format` (orchestrator/tree_parser): `parquet`
- `input_format` (triplet_processor): `parquet`
- `parquet_partitions`: `128`
- `parquet_compression`: `zstd`
- `min_support_value`: `0.5`
- `bootstrap`: `true`
- `bootstrap_options.debug_mode`: `false`
- `bootstrap_options.iterations`: `100`
- `bootstrap_options.seed`: unset
- `bootstrap_options.summary_only`: `false`

### Machine Learning (ghostparser.ml)

A small machine-learning baseline lives under `ghostparser.ml`. It consumes `summary_statistics.tsv` (the optional summary output from the pipeline) and provides explicit trainer modules for a multi-label Random Forest and a multi-label KNN baseline. Use them for quick prototyping and diagnostics; see `ML.md` for full usage and the data contract.

Run example:

```bash
python -m ghostparser.ml.random_forest -i results/summary_statistics.tsv -o results/ml_out
```

Or explicitly dispatch via the package entrypoint:

```bash
python -m ghostparser.ml --model random_forest -i results/summary_statistics.tsv -o results/ml_out
```

The KNN baseline is available as:

```bash
python -m ghostparser.ml.multi_knn -i results/summary_statistics.tsv -o results/ml_out
```

Note: `python -m ghostparser.ml` will not redirect to any model by default — you must pass `--model` to dispatch.

**Backend Details:**

- `stats_backend`: `standard` (Uses `scipy.stats` and `statsmodels` for DCT and KS tests)

**Configuration Precedence:**

`ghostparser.orchestrator` supports `-c/--config-file`; `tree_parser` and `triplet_processor` accept CLI parameters.

---

## Handled Errors

This section summarizes user-facing errors and validation failures that GhostParser modules can raise or report during execution.

### Orchestrator (`ghostparser.orchestrator`)

- `Error: Species tree file not found: ...` / `Error: Gene trees file not found: ...`
   Orchestrator exits early when required input files are missing.
- `✗ Error processing species tree: ...`
   Species-tree cleaning/parsing failed (typically malformed Newick, missing taxa, or filtering issues).
- `✗ Error generating triplets: ...`
   Triplet-generation setup failed (for example rooting/pruning/mapping failures).
- `✗ Error processing gene trees: ...`
   Gene-tree cleaning/rooting stage failed before extraction.
- `✗ Error in triplet inference or introgression inference stage: ...`
   A downstream extraction/inference/consolidation exception occurred; the appended message is the originating module error.

### Config Loading (`ghostparser.config`)

- `Config file not found: ...`
   The config path does not exist.
- `YAML support requires PyYAML to be installed`
   YAML config was provided but `PyYAML` is unavailable.
- `Config file must be .json, .yaml, or .yml`
   Unsupported config extension.
- `Config root must be a key/value object`
   Top-level config payload is not a mapping.
- `Missing required config field: ...`
   A required field (for example species path, gene path, or output-critical key) is absent or empty.
- `Missing required config field: outgroup(s)`
   No usable outgroup taxa were provided.
- `Config field ... must be a non-empty string when provided`
   Optional string/path fields were passed as empty or wrong type.
- `Config field ... must be an integer >= 0`
   Non-negative integer settings (for example `processes`, `parquet_partitions`) are invalid.
- `Config field ... must be a boolean when provided` / `Config field ... must be a numeric value`
   Boolean/float-style fields were provided with incompatible types.
- `Config field ... must be one of: ...`
   Choice-constrained fields (test/statistic/backend/format/correction) contain unsupported values.
- `Config field bootstrap_options.* ...`
   Bootstrap options failed validation (`iterations >= 1`, integer seed, boolean debug/summary flags).

### Tree Parsing and Extraction (`ghostparser.tree_parser`)

- `Tree file not found: ...`
   Input tree file path is missing.
- `Invalid Newick format in ...` (including `Tree <idx> has no terminal nodes`)
   Tree parsing failed, file is empty/invalid, or parsed trees are structurally unusable.
- `Missing required CLI argument: --outgroups`
   CLI outgroup argument is empty after parsing.
- `CLI argument --processes must be an integer >= 0`
   Process count is invalid.
- `species_triplet_trees is required and cannot be None`
   Internal extraction writer was called without required species-triplet mapping.
- `Missing species subtree mapping for triplets: ...`
   Some triplets have no mapped species subtree and cannot be serialized.
- `Missing species subtree for triplet header: ...` / `Missing topology summary for triplet header: ...`
   A triplet output header is missing required metadata fields.
- `Parquet output requires pyarrow to be installed`
   Parquet export requested but `pyarrow` is unavailable.
- `parquet_partitions must be >= 1`
   Invalid parquet partition count.
- `Triplet subtree labels do not match expected triplet`
   Extracted subtree taxa do not match the target triplet labels.
- `Could not determine sister-pair MRCA`
   Internal branch metric could not be computed because sister MRCA resolution failed.

### Triplet Topology Utilities (`ghostparser.triplet_utils`)

- `Triplet tree must contain exactly 3 terminal taxa`
   A triplet tree has missing/extra terminal taxa labels.
- `Could not determine rooted sister pair for triplet tree`
   Rooted triplet is unresolved/ambiguous (commonly polytomy or ambiguous rooting).
- `Tree taxa do not match provided ABC triplet`
   Topology classification was requested with an incompatible ABC taxon mapping.

### Triplet Processing (`ghostparser.triplet_processor`)

- `Unsupported ...` for discordant test, stats backend, summary statistic, tree-height strategy, p-value correction, or input format
   A selected method/format is outside supported choices.
- `Unknown triplet topology` / `Resolved topology roles require a valid species topology` / `Invalid species topology: ...`
   Internal topology state is inconsistent with supported canonical topologies.
- `Triplet tree must contain exactly 3 terminal taxa`
   Per-tree triplet metrics require exactly three labeled leaves.
- `species_triplet is required for tree height strategies A, B, and C`
   Taxon-specific height strategies were requested without ABC triplet labels.
- `Selected taxon ... not found in triplet tree`
   A/B/C-selected taxon is absent from the observed triplet tree.
- `Could not determine sister-pair MRCA for triplet tree`
   Sister-pair branch metrics could not be resolved for a triplet.
- `Invalid dis1_topology '...'. Expected 'BC' or 'AC'.`
   Inference text generation got an invalid discordant topology label.
- `borderline_margin must be >= 0`
   Hybrid KS helper margin parameter is invalid.
- `Unsupported summary statistic name: ...`
   Unsupported statistic requested in summary metric computation.
- `Bootstrap payload is not JSON-serializable: ...`
   Bootstrap debug output could not be converted to TSV-safe JSON.
- `alpha must be in (0, 1) for fdr_tsbh`
   Two-stage BH correction requires a strict `(0,1)` alpha.
- `Missing worker triplet entry for: ...`
   Multiprocessing worker was asked to analyze a triplet absent from its entry map.
- `CLI argument --processes must be an integer >= 0` / `--bootstrap-iterations must be an integer >= 1`
   Triplet processor runtime arguments are invalid.

### Triplet Text Input Validation (`parse_triplet_gene_trees_file`)

- `Invalid triplet header format (expected 5 tab-separated fields): ...`
   Triplet section header is malformed.
- `Invalid triplet header: ...` / `Invalid triplet count in header: ...` / `Invalid species tree in header: ...`
   Header fields are missing or cannot be parsed.
- `Invalid ABC label mapping in header: ...` / `Triplet/header label mapping mismatch ...`
   Header ABC mapping is malformed or inconsistent with listed taxa.
- `Invalid topology summary in header: ...` / `Invalid discordant role assignment in header: ...`
   Topology count/role metadata is malformed or logically inconsistent.
- `Triplet count/header mismatch ...`
   Header count disagrees with summary totals or parsed gene-tree rows.
- `Invalid triplet section format: expected blank line after header`
   Required blank separator after header is missing.
- `Invalid gene-tree line (expected Newick ending with ';'): ...`
   Gene-tree line is malformed text in triplet payload.
- `Duplicate triplet header encountered: ...`
   The same triplet appears more than once in a single input payload.

### Triplet Parquet Input Validation (`parse_triplet_gene_trees_parquet`)

- `Parquet input requires pyarrow to be installed`
   Parquet input parsing requested without `pyarrow`.
- `Invalid parquet triplet dataset: expected 'triplets/' and 'observations/' directories`
   Dataset directory layout does not match expected schema.
- `Duplicate triplet header encountered: ...`
   Duplicate triplet rows exist in parquet triplet metadata.
- `Invalid discordant role assignment in parquet row for ...`
   Parquet role metadata is inconsistent.
- `Observation references unknown triplet_id: ...`
   Observation partition references a missing triplet metadata row.
- `Triplet count/header mismatch ...`
   Observed row count does not match declared triplet count.

### ML Package Entrypoint (`ghostparser.ml.__main__`)

- `RuntimeError: Model module ... does not expose a 'main' function`
   Dispatch target module is importable but missing required CLI entry function.

### ML Dataset and Feature Validation (`ghostparser.ml.ml_utils`)

- `Input TSV is missing a header row` / `Input TSV contains no data rows`
   Training input file is structurally incomplete.
- `Invalid classes label: ... expected a 6-character 0/1 bitstring`
   Class label format is invalid for multi-label decoding.
- `Missing required target column: ...`
   Target column is absent from the TSV header.
- `No feature columns found after excluding the target column`
   Input has no usable predictors.
- `Missing feature value for '...'`
   At least one required feature cell is empty.
- `Non-numeric feature value for '...': ... string-valued columns must have at most 7 distinct values`
   Categorical feature expansion exceeded supported cardinality cap.
- `Cannot run stratified cross-validation because at least one class has fewer than 2 samples`
   Rare-class policy `error` blocks CV when class support is too low.

### ML Config Validation (`ghostparser.ml.config` and `ghostparser.ml.random_forest`)

- `Config file not found: ...` / `YAML support requires PyYAML to be installed`
   ML config file is missing or YAML dependency is unavailable.
- `Config file must be .json, .yaml, or .yml` / `Config root must be a key/value object`
   Unsupported config format or wrong top-level payload shape.
- `Missing required config field: ...` and typed validation messages (`must be ...`)
   Required runtime/model/evaluation fields are missing or typed incorrectly.
- `Do not place model hyperparameters or evaluation reporting controls at the top level ...`
   ML config shape is invalid; expected nested `model` and `evaluation` sections.
- `Place 'random_state' and 'n_jobs' at the top level, not under 'model'. ...`
   Runtime controls were placed in the wrong config section.
- `ConfigError(str(exc))` wrapping `rows_to_matrix` failures
   Random-forest CLI forwards TSV/feature validation failures as config errors.

### Hyperparameter Tuning (`ghostparser.ml.hyper_tune`)

- `Failed to initialize Weights & Biases for hyperparameter tuning ...`
   W&B initialization failed (authentication/network/mode setup).
- `Missing required config field: hyperparameter_tuning`
   Tuning config block is mandatory.
- `Unexpected top-level keys in hyperparameter_tuning config: ...`
   Config includes unsupported top-level keys.
- `Do not place ... at the top level in hyperparameter_tuning configs ...`
   Runtime/search/model/evaluation fields are in invalid sections.
- `Config field hyperparameter_tuning.search_space must be a mapping/object ...`
   Search-space block has invalid type.
- `Search space lists must contain at least one value`
   A tunable parameter list is empty.
- `Config field search_space must define at least one parameter`
   No search dimensions were provided.
- `Do not place runtime fields inside hyperparameter_tuning.search_space ...`
   Runtime keys were incorrectly placed inside search-space.
- `Unsupported search_space keys for ...`
   Search parameters do not match the selected model.
- `Cross-validation results do not include objective metric: ...`
   Objective metric name is not present in candidate CV output.
- `Unsupported tuning model: ...`
   Tuner received a model name outside supported model backends.
- `Hyperparameter tuning requires a feasible cross-validation split ...`
   CV fold setup is infeasible for class distribution/policy.
- `Search space expands to ... exceeds max_candidates=...`
   Grid size exceeded configured candidate cap.
- `No candidates were evaluated`
   Search loop finished without any valid candidate evaluation.

---

## Testing

Run all tests:

```bash
pytest
```

See [tests/TESTS.md](tests/TESTS.md) for detailed test documentation and test map.

---

## For Maintainers

### Triggering a Release Build

Release builds are triggered by pushing a version tag that matches `v*.*.*`.

```bash
# 1) Ensure main is up to date
git checkout main
git pull origin main

# 2) Create an annotated release tag
git tag -a v0.1.1 -m "Release v0.1.1"

# 3) Push the tag to trigger GitHub Actions release workflow
git push origin v0.1.1
```

Optional cleanup if you tagged the wrong commit:

```bash
# 1) Delete the local tag
git tag -d v0.1.1

# 2) Delete the remote tag
git push origin --delete v0.1.1

# 3) Manually delete the GitHub Release (if created)
#    Go to: https://github.com/asif256000/ghostparser/releases
#    Find the v0.1.1 release and click "Delete"
#    OR use GitHub CLI:
gh release delete v0.1.1
```

Note: Deleting the tag removes the git version reference, but the GitHub Release and uploaded assets remain unless explicitly deleted. For a full cleanup, delete both the tag and the release.