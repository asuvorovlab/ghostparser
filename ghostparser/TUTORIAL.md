# Tutorial

Every command here runs from the top of a clone of the repository, which
holds the example datasets.

## 1. Install

You need Python 3.10 or later and git. Clone the repository, which also
brings the example datasets, and install the package into a fresh
environment:

```bash
git clone https://github.com/asuvorovlab/ghostparser.git
cd ghostparser
python -m venv .venv
source .venv/bin/activate
pip install .
```

`python -m ghostparser.orchestrator --help` should now list the options.

## 2. What goes in

GhostParser needs three things: a species tree, a file of gene trees, and the
name of the outgroup. The first example is *Thuja*, in
[data/empirical_data/Thuja/](https://github.com/asuvorovlab/ghostparser/tree/main/data/empirical_data/Thuja):
four species, of which `Tdol` is the outgroup.

### The species tree

One tree in Newick format, alone in its file, holding every species you want
analyzed plus the outgroup(s). `species.tree` for *Thuja*:

```
((Tpli:0.0039302,(Tsut:0.00397085,Tst:0.00260092)1:0.00137507):0.00955777,Tdol:0.00955777);
```

- **It must contain the outgroup(s)**, whether it is rooted or not.
  GhostParser always does the rooting itself, on the outgroups you name, and
  then removes them. A tree you rooted yourself and stripped of its outgroups
  cannot be used: put the outgroups back in.
- With several outgroups, they must branch off together: one branch has to
  separate all of them from all the other species. If another species sits
  between them, the run stops and names the groups the species fall into.
- Branch lengths are optional. They only break ties when a gene tree's
  outgroup is chosen; the topology alone defines the triplets.
- Every label must be unique.

### The gene trees

Many trees in Newick format, one per line, in a single file. The first two
lines of `combined.trees`:

```
((Tst:0.0020071977,Tpli:0.0028124659)88:0.0007922821,Tdol:0.0121390378,Tsut:0.0024197772);
((Tst:0.0000026203,Tsut:0.0016602326)90:0.0033244607,Tpli:0.0000026203,Tdol:0.0167773355);
```

- **Labels must match the species tree exactly** (case included), and each
  label can appear at most once per tree. If you sampled several individuals
  of a species, keep one per species under the species tree's label.
- **Every branch needs a length.** The tests compare how tall the gene trees
  are, so the lengths are the data. Any unit works (substitutions per site
  from IQ-TREE or RAxML, for example) as long as every gene tree uses the
  same one. A missing length is read as 0 and the run reports which trees
  lacked one.
- **Each gene tree must carry an outgroup**, rooted or not (the lines above
  are unrooted). As with the species tree, GhostParser roots every gene tree
  itself on its outgroup and then removes it; a gene tree with no outgroup
  cannot be rooted and is dropped.
- **Trees may lack some species.** Each triplet uses only the gene trees that
  carry all three of its species.
- **Support values may stay** (the `88` and `90` above). A tree whose average
  support falls below `min_support_value` (default `0.5`) is dropped. On a
  0-100 bootstrap scale this practically never happens; on a 0-1 scale, such
  as posterior probabilities, it removes the poorly supported trees, so check
  which scale your trees use.
- When a gene tree leaves a triplet's three species in an unresolved
  polytomy, that tree is skipped for that triplet.

### The outgroup

One or more labels from the species tree, given with `-og`; separate several
with commas (`-og Out1,Out2`). It is required even when your trees are
already rooted, since it is what GhostParser roots them on. Choose species
clearly outside the group you are studying. With several outgroups, list the
one you trust most first: how each tree is rooted, and when that order
matters, is under [Rooting](orchestrator/README.md#rooting).

## 3. Check the inputs

Before analyzing anything, let GhostParser check that the trees fit together:

```bash
python -m ghostparser.orchestrator \
  -st data/empirical_data/Thuja/species.tree \
  -gt data/empirical_data/Thuja/combined.trees \
  -og Tdol \
  --output-folder results/thuja_check \
  --preflight-data-check
```

It reads and roots every tree exactly as a run would, writes
`results/thuja_check/preflight_data_check.txt`, and stops. For *Thuja* the
report ends with:

```
missing_lengths:    0  (gene trees lacking some branch length; ...)
usable_pairs:       2968 usable, 0 unresolved, 0 with a taxon absent

No blocking data issues detected for the checked trees/triplets.
```

What to look at in the report:

| Line | What a problem looks like |
| --- | --- |
| `missing_lengths` | Above 0: those trees, listed by line number, lack branch lengths. |
| `tangled` | An outgroup often found nested inside the ingroup in gene trees, which makes it a poor outgroup. |
| `usable_pairs` | Many `with a taxon absent`: labels that do not match the species tree, or sparse sampling. Many `unresolved`: polytomies. |

The command exits with status 0 when the check passes and 3 when it finds
problems, which the report lists.

## 4. Run the analysis

The same command without `--preflight-data-check`:

```bash
python -m ghostparser.orchestrator \
  -st data/empirical_data/Thuja/species.tree \
  -gt data/empirical_data/Thuja/combined.trees \
  -og Tdol \
  --output-folder results/thuja \
  --seed 42
```

`--seed` makes the run reproducible: the same seed gives the same numbers
whatever the number of processes. **The output folder is emptied before each run**, so never point it at the
folder holding your input trees.

`results/thuja/` now holds:

| File | What it is |
| --- | --- |
| `orchestrator_triplet_results.tsv` | The answer: one row per triplet. Open it in a spreadsheet. |
| `consolidation/` | The introgression maps (figures), with their numbers under `consolidation_data/`. |
| `metrics.txt` | The settings used and what each stage did, including the seed. |
| `processed_*` | The trees as analyzed: rooted, outgroups removed. |

## 5. Read the result

GhostParser splits the species into every set of three (a *triplet*) and
analyzes each one separately. *Thuja* has three species besides the outgroup,
so there is one triplet. Its row, abbreviated:

```
triplet           Tst,Tsut,Tpli
species_tree      (Tpli,(Tsut,Tst));
n_con             1346
n_dis1            930
n_dis2            692
dis1_topology     AC
decision_gate     PERM
classification    ambiguous
bootstrap_value   0.53
all_bootstrap     ambiguous=0.53,ghost_introgression=0.29,inflow_introgression=0.18,...
```

The triplet is written `(A, B, C)`, where A and B are sisters in the species
tree: here A = `Tst`, B = `Tsut`, C = `Tpli`. Each gene tree shows these
three species in one of three shapes. `n_con` counts the gene trees that
agree with the species tree, and `n_dis1` and `n_dis2` count the two ways to
disagree; `dis1_topology` names the more common one (`AC`: `Tst` paired with
`Tpli`).

Incomplete lineage sorting alone makes the two disagreeing shapes equally
common. GhostParser asks up to three questions in turn and stops at the first
answer that is no:

1. **Is one disagreeing shape more common than the other?** No: `no_introgression`.
2. **Do the disagreeing gene trees differ in height from the agreeing ones?**
   No: `inflow_introgression`.
3. **Taller or shorter?** Shorter: `outflow_or_bidirectional_introgression`.
   Taller: `ghost_introgression`. Neither shown: `ambiguous`.

| Classification | Meaning |
| --- | --- |
| `no_introgression` | Incomplete lineage sorting explains the disagreement. |
| `inflow_introgression` | Gene flow from C into one of the sisters. |
| `outflow_or_bidirectional_introgression` | Gene flow from one of the sisters into C, possibly in both directions. |
| `ghost_introgression` | Gene flow from an unsampled ("ghost") lineage older than the triplet. |
| `ambiguous` | Gene flow is likely, but its direction could not be told. |

The `inference` column says the same in words, naming the species, and
`decision_gate` names the question that settled it (`DCT`, `THT`, `PERM` for
questions 1, 2, 3).

For *Thuja*, 930 against 692 is a clear excess, and the heights differ, but
the two groups' average heights were shown to be close (`perm_decision` is
`equivalent`), so no direction is called: `ambiguous`.

`bootstrap_value` is how often the same answer came back when the gene trees
were resampled 100 times: 0.53 here, with `all_bootstrap` showing the rest
split between ghost (0.29) and inflow (0.18). Values near 1 are robust; a
value near 0.5 means the call could easily have gone another way.

Every column is described in the
[orchestrator guide](orchestrator/README.md#results-columns), and the tests
behind the three questions are explained there too.

## 6. A larger dataset, with a config file

*Jaltomata*, in
[data/empirical_data/Jaltomata/](https://github.com/asuvorovlab/ghostparser/tree/main/data/empirical_data/Jaltomata),
has 14 species plus the outgroup `Solyc`, so 364 triplets. Instead of a long
command line, put the settings in a YAML file, `jaltomata.yaml`:

```yaml
species_tree_path: data/empirical_data/Jaltomata/species.tree
gene_trees_path: data/empirical_data/Jaltomata/combined.trees
outgroup: Solyc
output_folder: results/jaltomata
species_filter: jaltomata_species.txt
seed: 42
```

Paths are read from the directory you run the command in, not from where the
YAML file sits.

`species_filter` limits the run to every triplet among the species listed in
`jaltomata_species.txt`, which keeps this first try small (10 triplets):

```text
JA0726, JA0432, JA0702
JA0456, JA0701
```

Run it:

```bash
python -m ghostparser.orchestrator -c jaltomata.yaml
```

Nine triplets come back `no_introgression`, and one,
`JA0432,JA0726,JA0702`, comes back
`outflow_or_bidirectional_introgression`: "introgression from JA0432 to
JA0702 with probable JA0702 to JA0432", with a bootstrap value of 0.55.

The map in `results/jaltomata/consolidation/introgression_bootstrap.png`
combines every triplet's call. Columns are the donor species and rows the
recipient, under the species tree, and each coloured cell is an introgression
event, shaded by its average bootstrap value. A dashed box marks the probable
reverse of an outflow call. Bars to the right mark species receiving ghost
introgression. `introgression_counts.png` draws the same map with the number
of supporting triplets instead.

To analyze all 364 triplets, delete the `species_filter` line. A run's time
grows with the number of triplets, and by default it uses every CPU available
to it.

A flag given beside `-c` overrides the file, so one YAML can serve several
datasets:

```bash
python -m ghostparser.orchestrator -c jaltomata.yaml --seed 7 --output-folder results/jaltomata_seed7
```

## 7. Your own data

1. Write the species tree and the gene trees as described in
   [step 2](#2-what-goes-in), with matching labels.
2. Copy `jaltomata.yaml`, change the paths and the outgroup, and drop
   `species_filter`.
3. Run the check (`--preflight-data-check`) and fix what it reports.
4. Run the analysis with a fixed `seed`, and record it with your results.

Useful settings, all described in the [configuration guide](README.md):

| Setting | Use it to |
| --- | --- |
| `species_rename_map` | Show readable names (`Homo sapiens`) in the results and figures instead of tree labels. |
| `species_filter`, `triplet_filter` | Analyze only some species, or only listed triplets. |
| `processes` | Cap the CPUs a run uses, on a shared machine (default: all available). |
| `bootstrap_threshold` | Keep weakly supported triplets out of the maps. |
| `generate_summary_stats` | Write `summary_statistics.tsv`, the input of the [machine-learning module](ml/README.md). |

## When something goes wrong

A failed run prints one line saying why and exits with a non-zero status. The
messages and statuses are listed under [Errors](../README.md#errors). The
most common at first:

- `Could not root the species tree ...`: the species tree does not contain
  the outgroup (removed when it was rooted beforehand, or misspelled in
  `-og`), or the outgroups do not branch off together. The message lists the groups the
  other species fall into, so you can see which ones to add to `-og`.
- `Missing required config field: outgroup`: `-og` was left out; it is
  needed even for trees that are already rooted.
- `Invalid Newick format in ...`: a tree file has a syntax error, often a
  missing `;` or an unbalanced parenthesis.
- A preflight report full of `with a taxon absent`: gene-tree labels that do
  not match the species tree.

Add `--debug` to see the full error trace when reporting a bug.
