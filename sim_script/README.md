# Simulating species and gene trees

`simulate_trees.py` simulates a species tree and gene trees under known
introgression for benchmarking GhostParser. Introgression can come from sampled
species, from ghost lineages (alive today but unsampled) and from extinct
lineages. Gene trees come from the multispecies coalescent in msprime
([Baumdicker et al. 2022](https://doi.org/10.1093/genetics/iyab229)). The script
writes a species tree and a gene-tree file the orchestrator reads directly, the
true events, and a figure of each replicate's scenario.

The simulator's packages are not GhostParser dependencies. Install them into
whichever environment runs the script:

```bash
pip install -r sim_script/requirements.txt
```

## Running

```bash
python sim_script/simulate_trees.py -c sim_script/configs/basic_newick.yaml
python -m ghostparser.orchestrator -c sim_out/basic_newick/rep_000/ghostparser_config.yaml
```

| Flag | Meaning |
|---|---|
| `-c`, `--config` | YAML or JSON config. Required. |
| `-o`, `--output-dir` | Overrides `output_dir`. |
| `--seed` | Overrides `seed`. |
| `--overwrite` | Resets an existing output directory. Without it, an existing directory is left alone and the run writes to `<name>_1`, `<name>_2` and so on. |
| `--debug` | Shows the full traceback on failure. |

A config or input problem prints one line and exits with 2 (config) or 1
(input, or a missing package). Every replicate's tree and events are checked
before anything is written.

| Example | Scenario |
|---|---|
| [configs/basic_newick.yaml](configs/basic_newick.yaml) | A fixed four-taxon tree with one pulse from C into B. |
| [configs/ghost_birth_death.yaml](configs/ghost_birth_death.yaml) | A birth-death tree whose own extinct lineages become ghosts, with random pulses from them into the sampled lineages; three replicates. |
| [configs/deep_newick.yaml](configs/deep_newick.yaml) | A fixed six-taxon tree with pulses between two ancestral lineages, from a ghost into an ancestor, and between two tips. |
| [configs/full.yaml](configs/full.yaml) | Every key at its default, with one example entry for each list. |

## Outputs

Each replicate writes to `rep_000/`, `rep_001/` and so on.

| File | Content |
|---|---|
| `species.tree` | The sampled species and outgroup(s): the orchestrator's species tree. |
| `species_full.tree` | The whole species tree, ghosts and extinct lineages included, each extinct tip ending at its extinction time. |
| `genes.tree` | `gene_tree_count` gene trees, one Newick per line, one tip per sampled species. |
| `events.tsv` | The events that were simulated: pulses, gene flow, size changes, bottlenecks and extinctions, with times in `time_units` and in generations. |
| `scenario.svg` | The species tree drawn against time with every event on it (see [Scenario figure](#scenario-figure)). |
| `demography.txt` | msprime's epoch-by-epoch account of the demography, for checking the model by eye. |
| `resolved_config.yaml` | The config with every default filled in, the species expanded and the master seed recorded. Rerunning it reproduces the replicate. |
| `ghostparser_config.yaml` | `species_tree_path`, `gene_trees_path`, `outgroup` and `output_folder` for the orchestrator. |

## Conventions

**Direction.** Introgression and gene flow are written forward in time: genes
move from `donor` into `recipient`. A pulse with `proportion: 0.3` replaces 30%
of the recipient's ancestry at that moment with donor ancestry.

**Times** run backwards from the present (0) and are in `time_units`:
`generations`, `coalescent` (one unit is ploidy × `population.size`
generations), or `years` (divided by `generation_time`). Rates are always per
generation: `gene_flow.rate`, `size_changes.growth_rate`, `clock.mean_rate`
and `locus.recombination_rate`.

**Branch references.** Wherever a key names a lineage (`donor`, `recipient`,
`branch`), it takes a species name, a prefix, or a list of them. One species
means its terminal branch. Several species mean the branch above their most
recent common ancestor, so any event can sit deep in the tree. This pulse runs
from the ancestor of D and E into the ancestor of A, B and C:

```yaml
introgression:
  - {donor: [D, E], recipient: [A, B, C], time_fraction: 0.5, proportion: 0.2}
```

On a random tree the clades are not known in advance, so `random_introgression`
reaches ancestors through `donor_branches` and `recipient_branches` instead.

**Windows.** A branch exists from its lower node (0 for an extant tip, the
extinction time for an extinct one, the split time for an internal branch) up
to its parent's split. Two branches can exchange genes only while both exist,
and every event must fall inside that shared window; an error names the window
otherwise. `time_fraction` places a pulse at that fraction of the way up the
window. Lineages that die out usually do so before today's terminal branches
begin, so pulses from extinct lineages mostly land on ancestral branches.

**Branch lengths.** Gene-tree lengths are substitutions per site: each branch's
duration in generations times a rate drawn for that branch alone, an
uncorrelated lognormal relaxed clock
([Drummond et al. 2006](https://doi.org/10.1371/journal.pbio.0040088)) whose
mean is `clock.mean_rate`. With `branch_length_units: normalized` (the
default), every tree is divided by the species tree's origin depth in expected
substitutions (its root age times `mean_rate`), so `species.tree` is 1 deep
from origin to present. Gene trees are taller, because their lineages coalesce
above the species root. One constant divides every tree, which keeps the
differences in gene-tree height that GhostParser's tests read. `output.scale`
then multiplies every length, so `scale: 100` puts the origin at depth 100.

## Configuration

| Key | Meaning |
|---|---|
| `seed` | Master seed; `null` draws one and records it in `resolved_config.yaml`. Each replicate gets its own streams for the tree, the events, the coalescent and the clock, so changing one part leaves the others' draws alone. |
| `output_dir`, `overwrite` | Where the replicates go, and whether to reset an existing directory. |
| `replicates` | Independent datasets. Each gets its own random tree (for `yule` and `birth_death`), random pulses and gene trees. |
| `time_units`, `generation_time` | See [Conventions](#conventions). |
| `species` | A list of names (`[A, B, C]`), or a `{prefix: count}` mapping: `{T: 3, E: 2}` gives T0, T1, T2, E0, E1. A group can carry its role inline: `E: {count: 2, extinct: true}` (also `ghost`, `outgroup`). Optional under `model: newick`, which takes the names from the tree. |
| `outgroup` | Names or prefixes of the sampled outgroup(s). Random trees attach them as the ingroup's sister. |
| `ghosts` | Names or prefixes alive today but never sampled. |
| `extinct` | Names or prefixes that die out before the present and are never sampled. |
| `species_tree.model` | `newick`, `yule` (pure birth; alias `pure_birth`) or `birth_death`. See [Species tree models](#species-tree-models). |
| `species_tree.newick`, `newick_file` | The tree for `newick`, with branch lengths in `time_units`. |
| `species_tree.root_age` | Age of the ingroup root for `yule` and `birth_death`. |
| `species_tree.birth_rate`, `death_rate` | Speciation and extinction rates; `death_rate` (default 0.5) belongs to `birth_death` and must stay below `birth_rate`. Only their ratio shapes the tree, because the tree is then scaled to `root_age`. |
| `species_tree.extinction_fraction` | `[low, high]`: an extinct species reaching the present in the tree survives a uniform draw from this share of its terminal branch, then dies. |
| `species_tree.extinct_from_process` | `birth_death` only: take the extinct species from lineages the process itself drove extinct, instead of truncating terminal branches. |
| `species_tree.outgroup_age` | The outgroup split as a multiple of `root_age`. |
| `species_tree.max_attempts` | `birth_death` only: how many simulations to try before giving up. |
| `population.size` | Effective size of every population. |
| `population.sizes` | Per-branch sizes: `[{branch: T0, size: 5000}]`. |
| `population.ploidy` | Ploidy for the coalescent time scale. One genome is sampled per sampled species. |
| `introgression` | Pulses: `{donor, recipient, proportion}` with `time` or `time_fraction`. |
| `random_introgression` | Pulses drawn anew in each replicate: `{count, donors, recipients, proportion, time_fraction, allow_sisters, donor_branches, recipient_branches}`. `donors` and `recipients` are pools of names or prefixes; `proportion` and `time_fraction` take a value or a `[low, high]` range (default `[0.1, 0.9]` for `time_fraction`). `donor_branches` and `recipient_branches` (each `tips`, the default, `ancestral` or `any`) choose which branches of a pool can take part: an ancestral branch belongs to the pool when its clade holds a pool species and otherwise only ghost or extinct lineages. Pulses between sister lineages are skipped unless `allow_sisters` (they leave every gene-tree topology unchanged). A pair is drawn only when its branches coexist, and pairs may repeat. |
| `gene_flow` | Continuous migration: `{donor, recipient, rate, start, end, symmetric}`. `start` and `end` default to the window's limits; `symmetric` adds the reverse direction at the same rate. |
| `size_changes` | `{branch, time, size, growth_rate}`; one of `size` or `growth_rate` at least. |
| `bottlenecks` | `{branch, time, strength}`: an instantaneous bottleneck equivalent to `strength` time units of coalescence. `{branch, time, proportion}`: each lineage joins a single ancestor with that probability. |

Every event takes `time` or `time_fraction`. A fixed `time` suits a fixed
Newick tree; on a random tree a branch may not exist at that time, so
`time_fraction` (that fraction of the way up the branch, or the shared window
for a pair) keeps an event valid in every replicate.
| `gene_tree_count` | Gene trees per replicate, each an independent locus. |
| `locus.sequence_length`, `recombination_rate` | Each locus's length in bases and its recombination rate. The defaults (1 base, no recombination) give one tree per locus. |
| `locus.tree_choice` | When recombination splits a locus into several trees, write the one at its `middle` or its `first` base. |
| `ancestry_model` | `hudson`, `smc`, `smc_prime`, `dtwf`, `{beta: {alpha}}` or `{dirac: {psi, c}}`. `dtwf` steps through every generation, so its cost grows with the tree's depth in generations. |
| `clock.model` | `ucln` (one lognormal rate per gene-tree branch) or `strict`. |
| `clock.mean_rate` | Mean substitutions per site per generation. |
| `clock.sigma` | Standard deviation of the log-rate under `ucln`. |
| `clock.among_loci_gamma_shape` | When set, each gene tree is also scaled by a Gamma draw with this shape and mean 1. |
| `output.branch_length_units` | `normalized`, `substitutions` (scaled by `scale`), or the time units `time_units`, `generations`, `coalescent`, which ignore the clock. |
| `output.scale` | Multiplies every length written. |
| `output.decimals` | Digits after the point in every Newick length. |
| `output.write_ghostparser_config`, `output.plot` | Whether to write `ghostparser_config.yaml`; the figure's format (`svg`, `png`, `pdf`) or `false`. |

Unknown keys are rejected. A number written as `1e-8` is read as a number.

## Species tree models

| Model | How the tree is made |
|---|---|
| `newick` | The given tree. Every species not listed as extinct must reach the present. An extinct species that also reaches the present is truncated by `extinction_fraction`. |
| `yule` | Pure birth over the ingroup, scaled so the ingroup root sits at `root_age`. Nothing dies unless species are listed as extinct, which are then truncated by `extinction_fraction`. |
| `birth_death` | A birth-death tree over the ingroup with the extinct lineages pruned, scaled to `root_age`, with listed extinct species truncated as under `yule`. With `extinct_from_process: true` the extinct species are instead lineages the process drove extinct, so their placement and timing follow the model. |

Random trees are drawn under the process conditioned on the number of tips,
with the general sampling approach
([Hartmann et al. 2010](https://doi.org/10.1093/sysbio/syq026)); tip names are
assigned at random.

## Scenario figure

Time runs up the vertical axis from the present, with the tips along the bottom.
Branch width follows population size, so a size change shows as a step.
Sampled lineages are solid, the outgroup grey, and lineages with no sampled
descendant dashed grey, with † where an extinct lineage ends. Pulses are arrows
from donor to recipient labeled with their proportion (hollow heads for random
ones), gene flow is a shaded band with an arrow and its rate, and bottlenecks
are diamonds. A plotting failure is reported as a warning and leaves the
simulated data in place.
