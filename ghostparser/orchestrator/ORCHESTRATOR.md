# ghostparser.orchestrator

`ghostparser.orchestrator` detects introgression from a species tree and a set
of gene trees. It decomposes the ingroup into species triplets, reads every
gene tree's version of each triplet, and runs a three-gate cascade of tests
per triplet that labels it `no_introgression`, `inflow_introgression`,
`outflow_introgression`, `ghost_introgression` or `ambiguous`. This document
explains the method; every flag and key is in [CONFIG.md](../../CONFIG.md).

## Running it

```bash
# minimal run
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup

# several outgroups, a species filter, a fixed seed, all available CPUs
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og Out1,Out2 \
    --species-filter species.txt --seed 42 --processes 0

# a config file; a flag given beside it overrides the file's value
python -m ghostparser.orchestrator -c run_config.yaml --seed 42

# one of the shipped scenario configs (CONFIG.md lists them)
python -m ghostparser.orchestrator -c sample_configs/orchestrator_preflight.yaml

# check the inputs and exit without analysis
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup --preflight-data-check
```

## The pipeline

1. **Rooting.** The species tree is rooted where the outgroups branch off and
   pruned of them; the remaining taxa are the ingroup. Each gene tree is
   rooted from its farthest outgroup, at the common ancestor of the
   outgroups outside the ingroup, and pruned of them.
2. **Triplets.** Every ingroup triplet is enumerated (or the ones a
   `triplet_filter` names, or every triplet among the species a
   `species_filter` names) and written `(A, B, C)` with A and B the
   species-tree sisters.
3. **Observations.** Each gene tree is measured once. For every triplet the
   engine reads, from each gene tree that carries all three taxa, the rooted
   topology of the three and a tree height `H(T)`.
4. **Tests.** Per triplet: the discordant count test, the tree-height test
   and the direction test, described below.
5. **Decision.** Once every triplet is measured, the count and tree-height
   p-values are corrected across the run and the cascade classifies each
   triplet. Bootstrap resampling of the gene trees gives each classification a
   support value.
6. **Outputs.** The results TSV, optionally `summary_statistics.tsv`, and the
   consolidated introgression maps.

## Why triplets

For three species with the rooted species tree `((A,B),C)`, a gene tree
carries one of three rooted topologies: the concordant `((A,B),C)` and the two
discordant ones, `((B,C),A)` and `((A,C),B)`. Under the multispecies
coalescent with no gene flow, discordance arises only from incomplete lineage
sorting (ILS): the A and B lineages fail to coalesce in the branch above their
split and are then equally likely to pair with C in either way. With that
branch `T` coalescent units long,

```
P(concordant) = 1 - (2/3) e^(-T)      P(each discordant) = (1/3) e^(-T)
```

(Hudson 1983, *Evolution* 37(1), 203-217; Pamilo & Nei 1988, *Molecular
Biology and Evolution* 5(5), 568-583; Degnan & Rosenberg 2009, *Trends in
Ecology & Evolution* 24(6), 332-340, https://doi.org/10.1016/j.tree.2009.01.009).
The two discordant topologies are exchangeable under ILS, so their counts are
expected to be equal. Gene flow between one of the sisters and C breaks that
symmetry, the same asymmetry the ABBA-BABA *D*-statistic reads from site
patterns (Green et al. 2010, *Science* 328(5979), 710-722,
https://doi.org/10.1126/science.1188021; Durand et al. 2011, *Molecular
Biology and Evolution* 28(8), 2239-2252,
https://doi.org/10.1093/molbev/msr048), and is the signal the first test
looks for.

Counts alone do not say what kind of gene flow produced the excess. The
*timing* of the discordant gene trees does: ILS-discordant trees coalesce
above the species-tree root of the triplet, a gene tree produced by
introgression between the sampled taxa coalesces at the time of that exchange,
and a gene tree carrying DNA from an unsampled ("ghost") lineage that diverged
before the triplet's root coalesces deeper still (Hibbins & Hahn 2019,
*Genetics* 211(3), 1059-1073, https://doi.org/10.1534/genetics.118.301831;
Tricou, Tannier & de Vienne 2022, *Systematic Biology* 71(5), 1147-1158,
https://doi.org/10.1093/sysbio/syac011). Comparing the heights of the
discordant gene trees with the heights of the concordant ones is what the
second and third tests do.

The more frequent discordant topology is `discordant1` (the `dis1_topology`
column, `BC` or `AC`), and its height sample against the concordant sample is
what the last two tests compare.

## Rooting

Read unrooted, a tree is rooted by its outgroups when one branch (or, inside
a polytomy, one node) parts every outgroup from every other taxon, in
whatever orientation the file was written. The species tree must root this
way; if other taxa sit between the outgroups the run stops and names the
groups those taxa fall into, so the ones that are outgroups can be added to
the list.

The rooted species tree ranks the outgroups by their distance from the
ingroup root, the summed branch lengths along the path, with the listed order
breaking exact ties. A species tree lacking any branch length keeps the
listed order; nothing else in a run reads species-tree lengths. An outgroup
the species tree lacks ranks after the rest.

A gene tree is rooted from its farthest outgroup, taken as the most reliable
because a closer outgroup can carry genes that sit nearer the ingroup through
incomplete lineage sorting or introgression. The farthest is the outgroup
with the longest mean path to the ingroup taxa in that gene tree. The mean
needs no root, and for outgroups outside the ingroup it ranks them as their
distance from the ingroup root does. A missing branch length counts as 0 here
and in every tree height; when two outgroups tie, as all do in a tree without
lengths, the species-tree rank decides.

Rooted on the farthest, any other outgroup that sits among the ingroup taxa
(inside the ingroup's common ancestor, on a branch that also leads to
ingroup taxa) is pruned without being used, because rooting on it would move
the root into the ingroup and change the rooted shape of every triplet
spanning the two. The tree is then rooted at the common ancestor of the
farthest and the other outgroups outside the ingroup, read unrooted, as for
the species tree; for the ingroup this is the same root as the farthest
outgroup's own branch. All outgroups are then pruned: no triplet contains
one, and removing a leaf changes no other taxon's rooted shape or heights.

`metrics.txt` and the preflight report give the ranking and count, per
outgroup, the trees in which it was the farthest, the trees rooted using it
and the trees in which it was pruned unused. A distant outgroup on a long
branch can attach inside the ingroup in a single gene and still be the
farthest, so a closer outgroup that is often pruned unused is the sign to
compare with a run that leaves the distant one out. A gene tree carrying no
outgroup, or nothing but outgroups, is dropped. `metrics.txt` and the
preflight report count the gene trees lacking some branch length.

## Reading a triplet from a gene tree

Each gene tree is cached once as three tables: the parent of every node, the
length of the edge above every node, and the lowest common ancestor (LCA) of
every pair of taxa. The LCA table fills in one post-order walk (a node is the
LCA of exactly the pairs drawn from two of its different children) and after
that no gene tree is touched again. A triplet's observation in a gene tree is
then read from three table lookups and a few short walks up the parent chain,
rather than by copying out a subtree.

**Topology.** For taxa `a`, `b`, `c` let `m(a,b)`, `m(a,c)`, `m(b,c)` be the
pairwise LCAs. In a binary rooted tree two of the three are the same node
(the LCA of all three, the triplet's root `r`) and the third lies strictly
below it: the pair whose LCA is the odd one out is the sister pair, and its
LCA is the sisters' node `s`. All three coinciding means the three taxa hang
off one polytomy, and the triplet is skipped in that gene tree. Because the
topology is read from node identity rather than from distances, a zero-length
internal branch is still a resolved topology, not a polytomy.

**Heights.** With `d(x)` the sum of edge lengths from leaf `x` up to `r`, and
`d(s)` the same sum from `s` up to `r`, the tree-height strategies are

| `tree_height_calculation_strategy` | `H(T)` |
| --- | --- |
| `AVG` (default) | `(d(A) + d(B) + d(C)) / 3` |
| `A`, `B`, `C` | `d(A)`, `d(B)`, `d(C)` |
| `INT` | `d(s)`, the internal branch from the triplet's root to the sisters' node |
| `SIS` | `d(s1) + d(s2) - 2 d(s)`, the patristic distance between the two sisters |

Every quantity is a sum along one path, never the difference of two depths
measured from the tree's root, so nothing cancels in floating point. The
same tables serve the species tree, from which each triplet's `(A, B, C)`
order and its `species_tree` subtree are read.

## The tests

### Gate 1: Discordant count test (DCT)

*Are the two discordant topologies equally frequent?* Under ILS alone they
are; an excess of one is gene flow. The test compares `n_dis1` and `n_dis2`
against a 50/50 split and ignores the concordant count.

- `chi-square` (default): Pearson's goodness-of-fit statistic against equal
  expected counts, `(n_dis1 - n_dis2)^2 / (n_dis1 + n_dis2)`, on one degree
  of freedom (Pearson 1900, *Philosophical Magazine* 50(302), 157-175,
  https://doi.org/10.1080/14786440009463897).
- `z-test`: the two-proportion z-test of `n_dis1/n` against `n_dis2/n` with
  the pooled variance at `p = 1/2`. The two proportions are complements, so
  `z^2 = 2 chi^2` on the same counts and the z-test rejects more readily.

A `0/0` split gives `p = 1`. The p-value is corrected across every triplet in
the run (see *Correction*), and the gate is significant at `corrected p <=
alpha_dct`. A triplet that fails it is `no_introgression`.

### Gate 2: Tree-height test (THT)

*Do the concordant and discordant1 height distributions differ at all?* The
two-sample Kolmogorov-Smirnov statistic is the largest gap between the two
empirical distribution functions, `D = sup_x |F_con(x) - F_dis1(x)|`, with
SciPy's exact or asymptotic p-value by sample size (Massey 1951, *Journal of
the American Statistical Association* 46(253), 68-78,
https://doi.org/10.1080/01621459.1951.10500769). The test is deliberately
sensitive to any difference: location, spread or shape.

If the discordant trees are distributed like the concordant ones, they
coalesce on the same timescale, which is what gene flow between the sampled
taxa produces: the triplet is `inflow_introgression`. The p-value is
corrected across every triplet, and the gate is significant at `corrected p
<= alpha_ks`. A significant KS test says the heights differ but not which way;
that is the third gate's question.

### Gate 3: Direction test

*Is the mean concordant height greater than, less than or indistinguishable
from the mean discordant1 height?* Concordant trees taller than the discordant
ones (`greater`) is `outflow_introgression`; discordant trees taller than the
concordant ones (`less`) is `ghost_introgression`, the signature of a lineage
that diverged before the triplet's root; neither is `ambiguous`.

**Statistic.** For concordant sample `x` (size `nx`) and discordant1 sample
`y` (size `ny`),

```
T = (mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)
```

with unbiased variances. The Welch standard error (Welch 1947, *Biometrika*
34(1-2), 28-35, https://doi.org/10.1093/biomet/34.1-2.28) matters because the
concordant sample is usually much the larger and the two variances differ: a
permutation test of the raw mean difference does not hold its level under
that combination, while the studentized one stays asymptotically valid
(Janssen 1997, *Statistics & Probability Letters* 36(1), 9-21,
https://doi.org/10.1016/S0167-7152(97)00043-6).

**Permutation null.** The pooled heights are randomly reassigned to two
groups of the original sizes and `T` is recomputed, variances included; the
distribution of those values is the null under exchangeability. Each tail
gets the add-one estimator

```
p = (1 + count) / (1 + resamples)
```

which counts the observed arrangement itself and is the correctly sized
Monte Carlo p-value (Phipson & Smyth 2010, *Statistical Applications in
Genetics and Molecular Biology* 9(1), Article 39,
https://doi.org/10.2202/1544-6115.1585). `perm_p_greater` is the upper tail,
`perm_p_less` the lower. The two form a family of size two and are corrected
against each other with the run's `p_value_correction` (under `bfn`,
`min(1, 2p)`); a direction is called when exactly one corrected tail is at or
below `alpha_perm`. Permutation p-values are never corrected across triplets:
a Monte Carlo p-value cannot fall below `1/(resamples + 1)`, and a Bonferroni
threshold of `alpha / (2 n_triplets)` sits below that floor for any sizable
run. The direction test therefore carries no run-wide error control; a
direction is conditional on the triplet having passed the first two gates.

**Adaptive stopping.** The first batch draws `min_resamples`. After each
batch a 95% binomial confidence interval (`ci_method`, default Wilson; Wilson
1927, *Journal of the American Statistical Association* 22(158), 209-212,
https://doi.org/10.1080/01621459.1927.10502953; Brown, Cai & DasGupta 2001,
*Statistical Science* 16(2), 101-133, https://doi.org/10.1214/ss/1009213286)
is placed around each tail's p-value and mapped onto the corrected scale. If
`alpha_perm` lies outside both intervals the decision cannot flip and the test
stops with `perm_converged = True`; otherwise the batch grows by a quarter and
the test continues until the total reaches `max_resamples`, the last batch
drawn whole. The interval is computed on the same `count` and `n` as the
p-value, so it brackets the value reported.

**Guards.** Four conditions skip the test and report `inconclusive` with a
`perm_note`: a group with fewer than two observations
(`insufficient_group_size`), a pooled sample with no spread
(`zero_pooled_variance`), two internally constant groups with different means
(`degenerate_observed_scale`), and fewer distinct group assignments than
`min_resamples` (`insufficient_permutation_support`, checked up to a pooled
size of 40). `both_tails_significant` marks an inconsistency a coherent null
cannot produce; no direction is reported.

**Null skewness.** `perm_null_skew` is the sample skewness of the permuted
statistics. It never drives a decision; a large magnitude means a few extreme
heights in the smaller group dominate the resampling, which is worth a look
at the underlying alignments.

**Equivalence.** A non-significant pair of tails means only that no direction
was established. Whether the two means were *shown* to be close needs its own
test, because absence of evidence is not evidence of absence. TOST (two
one-sided tests; Schuirmann 1987, *Journal of Pharmacokinetics and
Biopharmaceutics* 15(6), 657-680, https://doi.org/10.1007/BF01068419) turns
the null around: with margin `delta`, it tests `mean(x) - mean(y) <= -delta`
and `>= +delta` one-sidedly, and rejecting both confines the difference to
`(-delta, +delta)`. Each null is tested at its boundary by shifting the
concordant sample by `±delta` on the same permutations the direction test
drew; the TOST p-value is the larger of the two, and `equivalent` needs it at
or below `alpha_perm`. No multiplicity correction is needed: rejecting a union
of nulls only when every component rejects is an intersection-union test,
which holds its level whenever its components do (Berger 1982,
*Technometrics* 24(4), 295-300, https://doi.org/10.2307/1267823). See Lakens
2017 (*Social Psychological and Personality Science* 8(4), 355-362,
https://doi.org/10.1177/1948550617697177) on pairing it with a significance
test.

The margin is `delta = 0.5 sqrt((var(x) + var(y)) / 2)`: a Cohen's *d* of 0.5,
the conventional medium effect (Cohen 1988, *Statistical Power Analysis for the
Behavioral Sciences*, 2nd ed., Lawrence Erlbaum). An effect-size margin
shrinks relative to the standard error as gene trees accumulate, so more data
makes equivalence easier to establish; a margin in standard-error units would
sit permanently inside the null's own spread, because the studentized
statistic is a pivot, and equivalence could never be reached. `equivalent` and
`inconclusive` both classify as `ambiguous`; the distinction is what the
column reports.

**Known limitation.** The studentized permutation test is asymptotically
valid, not exact: its type-I error rate inflates when the smaller group falls
below roughly 30 observations *and* carries the larger spread. Treat
directional calls on triplets with very few discordant1 gene trees as
provisional; `n_dis1` is in the results for that reason.

## Decisions

### The cascade

"Significant" means `corrected p <= alpha` at every gate. The first gate that
is not significant settles the call:

| Gate | Not significant | Significant |
| --- | --- | --- |
| 1 · DCT | `no_introgression` | continue |
| 2 · THT | `inflow_introgression` | continue |
| 3 · direction | `equivalent`/`inconclusive` → `ambiguous` | `greater` → `outflow_introgression`; `less` → `ghost_introgression` |

`decision_gate` names the gate that settled the triplet (`DCT`, `THT`,
`PERM`). A test reported on a row below its gate took no part in the call.
Nothing is classified while the triplets are measured: the decision is made
once, after the run-wide correction, and is the single decision point.

Only what the cascade can read is measured. The direction test is skipped
below a settled gate, and under `no`/`bfn` the tree-height test is skipped
below a settled count gate; those rows leave the columns empty, with
`perm_note = direction_test_not_consulted` on a skipped direction test.
`diagnostic: true` measures every test for every triplet and changes no
result. Both skips are safe because every supported correction is monotone
(no corrected p-value is below its raw value, so a gate that failed raw cannot
clear corrected) and because permutation p-values are corrected within the
test only, so an unrun direction test moves no other triplet's numbers.

### Correction across the run

The DCT p-values of all triplets form one family and the KS p-values another;
each is corrected once, by `p_value_correction`, with the family size equal
to the triplet count. For the `j`-th smallest of `n` p-values:

| Method | Adjusted value | Controls |
| --- | --- | --- |
| `bfn` | `min(1, n p)` (Dunn 1961, *JASA* 56(293), 52-64, https://doi.org/10.1080/01621459.1961.10482090) | family-wise error rate |
| `holm` | `(n - j + 1) p_(j)` under a running maximum (Holm 1979, *Scandinavian Journal of Statistics* 6(2), 65-70) | family-wise error rate |
| `fdr_bh` | `(n / j) p_(j)` under a running minimum (Benjamini & Hochberg 1995, *JRSS B* 57(1), 289-300, https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) | false discovery rate, independent or positively dependent tests |
| `fdr_by` | BH times `sum_(i=1..n) 1/i` (Benjamini & Yekutieli 2001, *Annals of Statistics* 29(4), 1165-1188, https://doi.org/10.1214/aos/1013699998) | false discovery rate under any dependence |
| `no` | `p` | nothing |

Every multiplier is at least 1, which is the monotonicity the skips above rely
on. Under `no` and `bfn` a member's corrected value follows from its own raw
value and the count, so a member nothing reads can be left unmeasured; the
rank-based methods read every member, so the tree-height test is measured
for every triplet under them.

### Bootstrap support

`bootstrap_value` is the fraction of `bootstrap_options.iterations` resamples
of the triplet's gene trees (drawn with replacement; Efron 1979, *Annals of
Statistics* 7(1), 1-26, https://doi.org/10.1214/aos/1176344552) whose rerun of
the cascade lands on the reported classification; `all_bootstrap` gives the
fraction for every class. Iterations are judged against the same corrected
thresholds as the reported classification (under `no`/`bfn` as they run,
under the rank-based methods once every triplet's p-value for the same
iteration index is in), so the support measures the decision actually made.
Each iteration's direction test runs at a fifth of the resample budget without
the equivalence step, and a guard vote counts as `ambiguous`.

`bootstrap_perm_stat_ci_low` / `bootstrap_perm_stat_ci_high` bracket the
observed `T` between the `alpha_perm` and `1 - alpha_perm` percentiles of its
bootstrap distribution: the level at which an interval and a one-sided test
agree, so an interval excluding zero corresponds to a directional rejection,
and a direction can be read as an effect size in standard errors rather than
only as a threshold crossing.

`seed` drives every random draw. Each triplet derives its own stream from
`(seed, triplet)`, so a run is reproducible at any worker count.

## Shape diagnostics

`shape_diagnostics: true` adds five descriptive columns per height group
(`con_*`, `dis1_*`, `dis2_*`) to the results TSV; nothing in the cascade reads
them.

| Column | Meaning |
| --- | --- |
| `<group>_n_modes` | Peaks in a Gaussian kernel density estimate at Scott's bandwidth (Scott 1992, *Multivariate Density Estimation*, Wiley). A count at one bandwidth: a strongly skewed unimodal sample can read as several. |
| `<group>_modes_p` | Silverman's critical-bandwidth bootstrap p-value for unimodality (Silverman 1981, *JRSS B* 43(1), 97-99, https://doi.org/10.1111/j.2517-6161.1981.tb01155.x): the smallest bandwidth giving one mode is found, and resamples from that smoothed density are counted for still needing more. At or below `alpha` is evidence of a second mode; it needs about three standard deviations of separation to see one. |
| `<group>_skew` | Sample skewness. |
| `<group>_excess_kurtosis` | Sample excess kurtosis, `0` for a normal. |
| `<group>_tail_xi` | Shape parameter of a generalized Pareto distribution fitted to the exceedances over the group's 90th percentile (peaks over threshold; Pickands 1975, *Annals of Statistics* 3(1), 119-131, https://doi.org/10.1214/aos/1176343003): positive is a power-law tail, `0` exponential decay, negative a bounded tail. |

Groups with fewer than 20 observations, or no spread, leave the columns empty;
`tail_xi` also needs 10 exceedances. The modality bootstrap is the expensive
part, which is why the diagnostics are off by default.

## Preflight data check

`--preflight-data-check` roots and reads the trees exactly as a run would,
collects every structural problem instead of stopping at the first, writes
`preflight_data_check.txt` and exits without analysis. The report groups the
issues by category with counts and examples, attributes them to the species
tree or the gene trees, accounts for every triplet/gene-tree pair (measured,
unresolved, or skipped for a missing taxon), and reports the gene-tree
rooting counts described under *Rooting*. It walks at most
`preflight_triplet_cap` triplets (default 15,000; `0` lifts it) and says so
when the cap binds. Passing means the data can be processed, not that the
result is meaningful. `sample_configs/orchestrator_preflight.yaml` is a
ready-made check.

## Outputs

The output folder is reset before the run under the default `overwrite: true`,
so give it a directory of its own, never one holding the input trees.

- `orchestrator_triplet_results.tsv`: one row per triplet (columns below).
- `summary_statistics.tsv`: with `generate_summary_stats`: per triplet, the
  mean, median, mode, variance, entropy, minimum and maximum of the average
  tree height, the internal branch and the sister distance, for the
  concordant, discordant1 and discordant2 gene trees (63 columns), plus the
  identity columns, the counts, the classification and `bootstrap_value`.
- `processed_<species tree>` / `processed_<gene trees>`: the cleaned, rooted
  and pruned trees, in the input labels.
- `metrics.txt`: the run parameters, then one block per stage (species tree,
  gene trees, inference, introgression maps) giving what it processed, its
  timings and what it found: the rooting counts, the classification and gate
  counts and the permutation-test convergence summary.
- `consolidation/`: the introgression maps and their TSV matrices.

### Results columns

| Column | Meaning |
| --- | --- |
| `triplet`, `abc_mapping`, `species_tree` | The `(A, B, C)` labels with A and B the species-tree sisters, the same as `A=…;B=…;C=…`, and the triplet's species subtree. |
| `dis1_topology` | `BC` or `AC`: the more frequent discordant topology (ties to `BC`). |
| `n_con`, `n_dis1`, `n_dis2`, `analyzed_trees` | Gene trees per topology, and the number carrying all three taxa with a resolved topology. |
| `most_frequent_matches_concordant` | Whether the concordant count is at least both discordant counts. |
| `dct_chi_stats` or `dct_z_score`, `dct_p_value` | Count-test statistic (named after the backend) and raw p-value. |
| `dct_p_val_<method>_corr`, `dct_significant` | The run-wide corrected p-value (omitted under `no`) and the gate-1 flag. |
| `ks_statistic`, `ks_p_value`, `ks_p_val_<method>_corr`, `ks_significant` | The KS test, its correction and the gate-2 flag; empty where the test was not measured. |
| `perm_statistic`, `perm_p_greater`, `perm_p_less` | The observed `T` and the raw one-tailed p-values. |
| `perm_p_greater_<method>_corr`, `perm_p_less_<method>_corr` | The tails corrected against each other; the values compared to `alpha_perm`. |
| `perm_p_tost` | The equivalence p-value, populated only when neither tail was significant. |
| `perm_n_resamples`, `perm_converged`, `perm_null_skew`, `perm_note` | Permutations drawn, whether the stopping rule settled before the budget, the null's skewness, and a guard or skip note. |
| `perm_decision` | `greater`, `less`, `equivalent` or `inconclusive`; consulted only when `decision_gate` is `PERM`. |
| `decision_gate`, `classification`, `inference` | The gate that settled the call, the class, and the direction in words naming the species. |
| `bootstrap_value`, `all_bootstrap` | Support for the classification and the fraction per class. |
| `bootstrap_perm_stat_ci_low`, `bootstrap_perm_stat_ci_high` | The bootstrap percentile interval on `T`. |
| `con_*`, `dis1_*`, `dis2_*` | The shape diagnostics, with `shape_diagnostics`. |
| `bootstrap_*` | Per-iteration records, with `bootstrap_options.diagnostic`. |

The results TSV carries statistics and p-values only; descriptive per-group
statistics are in `summary_statistics.tsv`, whose `*_avg_tree_height_*`
columns always average the three root-to-tip distances whatever the height
strategy.

## Consolidation

The final stage turns the classified triplets into a directed introgression
map over the ingroup taxa. For a triplet `(A, B, C)` whose discordant1
topology pairs `C` with one sister `S`:

| Classification | Event recorded |
| --- | --- |
| `inflow_introgression` | edge `C → S` |
| `outflow_introgression` | edge `S → C` |
| `ghost_introgression` | ghost target: the sister `discordant1` leaves out |
| `no_introgression`, `ambiguous` | nothing |

Each heatmap cell and ghost bar is the mean `bootstrap_value` over the
triplets that produced that edge or target, and nothing else:

```
avg(source → target) = sum(bootstrap_value over triplets producing the edge) / count(those triplets)
```

The `*_raw_sum.tsv` and `*_supporting_count.tsv` matrices hold the numerator
and denominator, and `introgression_matrix_sampled_non_sister.tsv` counts, per
taxon pair, the triplets in which the two are not the species-tree sisters,
so the averages can be reweighted against co-occurrence. Without a bootstrap
every classified triplet weighs 1. The figure `introgression_combined.png`
draws the heatmap (rows targets, columns sources, `cividis`) under the species
tree, with the ghost bars beside it: bar length is the ghost support and bar
colour says whether the taxon is also the target of a sampled edge. Taxa in
`outgroup` never appear.

## Parallelization

Triplets are split into chunks across `processes` workers on one machine;
`0` uses the CPUs the process may run on, which under a scheduler or container
is the allocation rather than the machine's cores. The gene-tree cache is
built once and shared with the workers. A job spanning several machines uses
only the one the run starts on.
