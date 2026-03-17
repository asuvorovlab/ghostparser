# Ghostparser

## Overview

**Ghostparser** is a phylogenetic introgression pipeline centered on `ghostparser.orchestrator`. The orchestrator runs two internal stages (`tree_parser` and `triplet_processor`) and produces final triplet-level inference outputs. It supports both command-line and configuration-file modes for flexible execution, with support for multiprocessing, multiple statistical backends, and configurable statistical thresholds.

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
6. [Testing](#testing)
7. [For Maintainers](#for-maintainers)

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

- `tree_parser` standardizes species/gene trees, roots on outgroup(s), and writes triplet-specific gene-tree blocks.
- `triplet_processor` applies the GhostParser statistical decision pipeline to each triplet block.
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
- `--triplet-filter`
- `--processes`
- `--min-support-value`
- `--discordant-test`
- `--summary-statistic`
- `--stats-backend`
- `--tree-height-calculation-strategy`
- `--p-value-correction`
- `--alpha-dct`, `--alpha-ks`

#### Primary Outputs

1. `unique_triplets_gene_trees.txt`
2. `orchestrator_triplet_results.tsv`
3. `metrics.txt`

### Tree Parser (Submodule)

Use this module when you only want preprocessing + triplet extraction.

```bash
python -m ghostparser.tree_parser -st species.tree -gt genes.tree -og OutGroup
```

Config mode:

```bash
python -m ghostparser.tree_parser -c sample_configs/tree_parser_minimal.yaml
```

#### Core Behavior

- Removes support labels and preserves branch lengths
- If support values are present, removes trees with average support below `min_support_value`
- Roots on outgroup(s), prunes outgroup clade, and logs excluded taxa
- Writes processed trees and `unique_triplets_gene_trees.txt`

### Triplet Processor (Submodule)

Use this module when you already have `unique_triplets_gene_trees.txt` and only need inference.

```bash
python -m ghostparser.triplet_processor --input-path unique_triplets_gene_trees.txt
```

Config mode:

```bash
python -m ghostparser.triplet_processor -c sample_configs/triplet_processor_minimal.yaml
```

#### Core Behavior

- Runs DCT (`chi-square` or `z-test`)
- Runs KS tree-height test when DCT is significant
- Applies summary-statistic comparison (`median`, `mean`, or binned `mode`) for final classification
- Supports tree-height strategy selection via `tree_height_calculation_strategy`: `AVG` (default) or taxon-specific `A|B|C`

---

---

## Orchestrator Input/Output

### Input Expectations

**Required Arguments:**
- `--species-tree-path` (alias `-st`): Path to the species tree file in Newick format
- `--gene-trees-path` (alias `-gt`): Path to the gene trees file in Newick format
- `--outgroups` (alias `-og`): Outgroup species identifier(s). Use comma-separated taxa for multiple outgroups.

**Optional Arguments:**
- `--output-folder`: Output folder relative to the input data folder (default: same folder as input data).
- `--triplet-filter`: Path to a triplet filter file (comma-separated taxa per line). When provided, only those triplets are processed. Triplets containing taxa missing from the species tree are skipped with a warning.
- `--processes`: Number of worker processes for multiprocessing. Defaults to `0` (all cores). Ignored if `--no-multiprocessing` is set.
- `--no-multiprocessing`: Disable multiprocessing. Processes triplets sequentially using a single worker. Useful for debugging or on systems with limited resources.

### Output Files

The orchestrator generates these output files:

1. **`processed_species.tree`** - Processed species tree with support values removed and outgroup rooting applied
2. **`processed_gene_trees.tree`** - Processed gene trees with support values removed and outgroup rooting applied
3. **`unique_triplets_gene_trees.txt`** - Triplets normalized to `A,B,C` (with `A,B` as species sisters), with required header format `triplet<TAB>count<TAB>species_tree` (non-empty species subtree)
4. **`metrics.txt`** - Metrics log with warnings, timings, and counts
5. **`orchestrator_triplet_results.tsv`** - Final triplet-level classification results (`no_introgression`, `outflow_introgression`, `inflow_introgression`, `ghost_introgression`, or `unresolved`)

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

Core defaults are centralized and applied consistently in both CLI mode and config-file mode:

**Statistical and Processing Defaults:**
- `discordant_test`: `chi-square`
- `summary_statistic`: `median`
- `stats_backend`: `standard`
- `tree_height_calculation_strategy`: `AVG`
- `p_value_correction`: `no`
- `alpha_dct`: `0.01`
- `alpha_ks`: `0.05`

**Execution Defaults:**
- `processes`: `0` (all available CPU cores)
- `output_folder` (orchestrator/tree_parser): `./results`
- `min_support_value`: `0.5`

**Backend Details:**
- `stats_backend`: `standard` (Uses `scipy.stats` and `statsmodels` for DCT and KS tests)

**Configuration Precedence:**
When `-c/--config-file` is provided, other CLI flags are ignored with a warning.

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