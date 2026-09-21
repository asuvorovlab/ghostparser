# ghostparser.orchestrator

`ghostparser.orchestrator` is GhostParser's introgression engine. It fuses triplet
subtree extraction and per-triplet inference into a single streaming pass, so
the intermediate triplet-gene-trees dataset is never written to disk or reloaded
into memory.

This document explains how the module works. The complete reference for every
flag and config key lives in [CONFIG.md](../../CONFIG.md#orchestrator-primary-module).

## Running it

```bash
# minimal run
python -m ghostparser.orchestrator -st species.tree -gt genes.tree -og OutGroup

# multiple outgroups, custom output folder, all cores
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og Out1,Out2 \
    --output-folder results --processes 0

# a filtered set of triplets, no consolidation plots
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup \
    --triplet-filter triplets.txt --no-consolidation

# every triplet among a filtered set of species
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup \
    --species-filter species.txt

# config-file mode (JSON or YAML); other CLI flags are ignored
python -m ghostparser.orchestrator -c run_config.yaml

# start from a shipped sample: orchestrator_minimal.yaml has just the required
# inputs, orchestrator_full.yaml lists every key at its default
python -m ghostparser.orchestrator -c sample_configs/orchestrator_minimal.yaml

# check the input data and exit, without running any analysis
python -m ghostparser.orchestrator \
    -st species.tree -gt genes.tree -og OutGroup --preflight-data-check
```

## Configuration

Three inputs are required — the species tree (`-st`), the gene trees (`-gt`),
and the outgroup(s) (`-og`) — and everything else has a default.
`-c/--config-file` is the only CLI-only option; when given, the file supplies
every setting and the other CLI flags are ignored with a warning.

A handful of settings are config-file-only (`discordant_test`,
`tree_height_calculation_strategy`, `min_support_value`,
`generate_summary_stats`, `shape_diagnostics`, and the `bootstrap_options` and
`permutation_options` blocks). One default is worth calling out: `p_value_correction` defaults to
`bfn` (Bonferroni), and it is applied both run-wide across triplets and inside
each permutation test across its pair of one-tailed p-values.

See [CONFIG.md](../../CONFIG.md#orchestrator-primary-module) for every key, its
default, and its allowed values.

## Preflight data check

`--preflight-data-check` (or `preflight_data_check: true`) turns the run into a
data validation pass. It short-circuits `run_orchestrator` immediately after
the output directory is prepared, so no analysis runs and the only artifact is
`preflight_data_check.txt`.

`preflight.run_preflight_data_check` replays the same structural logic the
engine uses, but collects failures instead of raising on the first one:

1. Standardize the species tree, root it where the outgroups branch off and
   prune them exactly as the run does, build one `triplet_geometry.TripletGeometry`
   over the ingroup that remains, and normalize each triplet to
   A/B/C from the sister pair `triplet_geometry.triplet_subtree_shape` reports,
   recording triplets whose rooted topology cannot be resolved.
2. Root every gene tree on the first outgroup label it carries, recording the
   ones where none is present, how many rooted on each outgroup, and how many
   carry outgroups that do not all lie on one side of the other taxa.
3. Cache each rooted gene tree's geometry and ask
   `triplet_geometry.triplet_resolution` about every triplet, recording the ones
   it reports as missing a taxon or unresolved.

Both tree sides go through the cached geometry the engine itself uses, so the
skip decisions reported here are the ones a run would make.

Every failure becomes an `Issue` with a dotted category. The report groups them
by category with counts, up to 25 examples each naming the gene-tree index and
triplet plus the offending input line, and an attribution summary separating
species-tree causes from gene-tree causes. A `usable_pairs` line accounts for
every triplet/gene-tree pair the check looked at, splitting them into the ones a
run could measure, the ones it would drop as unresolved, and the ones skipped
because that gene tree does not carry all three taxa. The last group is not a
defect -- a gene tree need not be complete, and the engine skips those pairs
too -- so it is counted rather than reported as an issue. `PreflightResult.passed` is `True`
only when nothing was detected.

`rooted_on` and `tangled` lines count, per outgroup, the gene trees rooted
using it and the gene trees in which it sat among the ingroup taxa and was
pruned without being used for rooting; a `tangled_trees` line counts the
trees with at least one tangled outgroup, and an `order_decided` line counts
those in which no set of outgroups held a majority, so the listed order chose
which to root on. None is a defect: a run keeps such a tree, rooted where the
largest set of its outgroups branches off. They are reported because a
frequently tangled outgroup is usually a distant one on a long branch, which
single genes place unreliably, and because the `order_decided` trees are the
ones whose rooting depends on the listed order, so the outgroup you trust
most belongs first.

Generated triplets are capped at `preflight_triplet_cap` (default 15,000; `0`
lifts it), which the runner passes through as `max_triplets`; a bound cap is
itself reported as `analysis.triplet_cap_applied`, so a partial check never
reads as a complete one. Triplets named by a `triplet_filter` are not capped;
those generated from a `species_filter` are, like the full ingroup's.

Three conditions make the check itself impossible and raise `ValueError`
instead: no outgroups given, a species-tree file that does not hold exactly one
tree, and outgroups that do not root the species tree — none of them present,
every taxon an outgroup, or outgroups that branch off at more than one point,
reported with the same message the run gives (see "Species preprocessing"
below).

## Orchestrator summary

`runner.run_orchestrator(config)` coordinates the run:

1. **Species preprocessing** — `trees.clean_and_save_trees` standardizes the
   species tree and drops trees whose mean internal support is below
   `min_support_value`, then roots it where the outgroups branch off and
   prunes them, leaving the ingroup taxa. Read unrooted, the tree must have one
   branch -- or, inside a polytomy, one node -- that parts every outgroup from
   every other taxon; it is rooted there whatever orientation the file wrote
   it in. Otherwise other taxa sit between the outgroups, and the run stops
   with a message that lists the groups those taxa fall into, largest first,
   so the groups that are outgroups can be added to the outgroup list.
2. **Triplet setup** — `trees.generate_triplets` enumerates every ingroup
   triplet, or `trees.read_triplet_filter_file` plus
   `trees.filter_triplets_by_taxa` restricts the run to the triplets a filter
   names, or `trees.read_species_filter_file` restricts the enumeration to the
   ingroup species a filter names (names that are not ingroup taxa are skipped
   with a warning; fewer than three left is an error). The two filters
   exclude each other at config time. One `triplet_geometry.TripletGeometry`
   is built over the species tree and each triplet is read out of it:
   `triplet_geometry.triplet_subtree_shape` gives the sister pair, which
   normalizes the triplet to `(A, B, C)` with A and B the
   species-tree sisters, and the edges the triplet's induced subtree would
   carry, written out as the species subtree.
3. **Gene-tree preprocessing** — `trees.clean_and_save_gene_trees` cleans each
   gene tree, roots it with `trees.root_gene_tree` and writes it without its
   outgroups; a gene tree carrying no outgroup, or nothing but outgroups, is
   dropped. The rooting is the species tree's: the outgroups present root the
   tree together when they part from the other taxa at one point. When they
   do not, the largest subset that does roots the tree and an outgroup
   tangled among the ingroup taxa is pruned without being used, so the
   outgroups that still sit together outvote one that a single gene placed
   among the ingroup; when no subset holds a majority the one listed earliest
   wins. Pruning the outgroups changes nothing downstream -- no triplet
   contains one, and removing a leaf leaves every other taxon's rooted shape
   and heights as they are. `metrics.txt` reports, per outgroup, how many
   trees were rooted using it and in how many it was tangled and pruned
   unused, and how many trees the listed order settled.
4. **Fused extraction + inference** — `stream.stream_triplet_results` caches
   every gene tree's geometry once (`stream.build_run_geometry`), then walks the
   triplets, reads each one's observation out of that cache
   (`triplet_geometry.geometry_observation`), and immediately runs
   `inference.analyze_triplet_from_observations`. Only the small result object
   is retained.
5. **Run-wide correction** — the multiple-testing correction is applied once
   across all triplets, because a global correction needs every p-value in a
   single pass, and every triplet is classified there.
6. **Display names** — every step so far works in the trees' own labels. With a
   `species_rename_map`, each result is rebuilt under its display names:
   `triplet` through `trees.rename_taxon_labels`, and the `species_tree` Newick
   through `trees.rename_newick_labels`, which quotes any name the format
   cannot carry bare. Nothing read back later — the
   processed tree files, the outgroup, a triplet or species filter — ever
   holds a display name, so the names are free to contain spaces, dots or
   quotes.
7. **Writing** — `inference.write_pipeline_results` emits
   `orchestrator_triplet_results.tsv`; `inference.write_summary_statistics_tsv`
   emits `summary_statistics.tsv` when `generate_summary_stats` is set.
8. **Consolidation** — `consolidation.generate_introgression_maps`
   writes the map artifacts into a `consolidation/` subfolder. It orders the
   taxa by the processed species tree on disk, which is in tree labels, and
   takes the same `rename_map` to put the display names on that tree in
   memory so it lines up with the results.

## Per-triplet inference

For each triplet the engine classifies every gene tree's subtree into one of
three topologies — concordant (matching the species tree) plus two discordant
alternatives — and records a tree height H(T) per the configured strategy. The
more frequent of the two discordant topologies is `discordant1`, and its height
sample against the concordant one is what the second and third gates compare.
Three tests are then read in a fixed order, and the first one that comes back
non-significant settles the call:

1. **Discordant count test (DCT)** — compares the two discordant counts
   (`inference.run_discordant_count_test`, chi-square or z-test). A triplet
   whose corrected p-value exceeds `alpha_dct` is `no_introgression`, and
   nothing below is consulted.
2. **Tree-height test (THT)** — a two-sample KS test between the concordant and
   discordant1 height distributions (`inference.run_two_sample_ks_test`). A
   corrected p-value above `alpha_ks` makes the triplet `inflow_introgression`.
3. **Direction test** — the adaptive studentized permutation test on the same
   two height samples (`permutation.run_studentized_permutation_test`):
   `greater` gives `outflow_introgression`, `less` gives
   `ghost_introgression`, and `equivalent` and `inconclusive` both give
   `ambiguous`.

"Significant" means `corrected p-value <= alpha` at every gate, so a value
landing exactly on the threshold counts as significant and the cascade moves on
to the next test.

The classification and the name of the test that settled it — `DCT`, `THT`
or `PERM` — are produced together as one pair, which populates the
`classification` and `decision_gate` columns, so the two are derived in a
single pass and cannot drift apart.

Nothing is classified while the triplets stream: the per-triplet pass returns
measurements only, with the decided fields left empty. The classification is
made once, in the run-wide correction pass, after every triplet's p-values
are in hand and corrected — see "The two run-wide correction families". That
pass is the single decision point in the pipeline.

How much of a triplet actually gets measured is set by `diagnostic`. Off (the
default), a gate that settles the call stops the work there, so the `perm_*`
block is empty on those rows and carries
`perm_note = direction_test_not_consulted`, and under `no`/`bfn` the
tree-height columns are empty below a settled count gate; on, all three tests
run for every triplet and every column is populated. The classification is
the same either way under every correction method — see "The point estimate's
short-circuit" for why.

Either way, read `decision_gate` before reading `perm_decision`: only `PERM`
means the direction result produced the classification. On a diagnostic run,
a row carrying `decision_gate = THT`, `perm_decision = ambiguous`, and
`classification = inflow_introgression` is consistent — the direction test ran
and was recorded, but the tree-height test had already settled the call.

Bootstrap resampling (on by default) redraws the triplet's gene-tree
observations with replacement, reruns the cascade on each resample, and
aggregates the per-iteration classifications into `bootstrap_value`. An
iteration measures only what its own vote reads, whatever `diagnostic` says,
and runs the direction test at a fifth of the configured resample budget with
the equivalence step off — see "Inside the bootstrap" under gate 3 and
"Correction inside the bootstrap".

## The statistical tests

Each gate answers a different question, and each is computed by a named library
routine rather than by hand. This section states what each test measures, how
its value is obtained, when it runs, and where its assumptions bite.

### Gate 1 — Discordant count test

**Question.** Are the two discordant topologies equally frequent?

Under incomplete lineage sorting alone, the two discordant histories are
exchangeable and should appear about equally often; an excess of one of them is
the signal that something other than ILS — introgression — has acted (Huson et
al. 2005, *RECOMB*, https://doi.org/10.1007/11415770_18). The test therefore
asks only whether `n_dis1` and `n_dis2` depart from a 50/50 split, and ignores
the concordant count entirely.

Two backends, selected by `discordant_test`. The statistic lands in
`dct_chi_stats` or `dct_z_score` and the p-value in `dct_p_value`:

- `chi-square` (default) — `scipy.stats.chisquare([n_dis1, n_dis2])`, a
  goodness-of-fit test against equal expected counts. The statistic is
  `sum((observed - expected)^2 / expected)` with `expected = (n_dis1 + n_dis2) / 2`,
  which simplifies to `(n_dis1 - n_dis2)^2 / (n_dis1 + n_dis2)`, compared
  against a chi-square distribution on one degree of freedom.
- `z-test` — `statsmodels.stats.proportion.proportions_ztest` with
  `count=[n_dis1, n_dis2]`, `nobs=[total, total]`, `alternative="two-sided"`:
  the two-proportion z-test comparing `n_dis1 / total` with `n_dis2 / total`
  under the pooled variance `p (1 - p) (2 / total)` at `p = 1/2`. The two
  proportions are complements of one another, so that variance is half the
  binomial variance of their difference and the statistic works out to
  `z^2 = 2 * chi^2` on the same counts: the z-test rejects more readily than
  the chi-square does.

A zero/zero split short-circuits to `(0.0, 1.0)` rather than dividing by zero.

**What is compared to `alpha_dct`.** The raw p-value is corrected once across
every triplet in the run — the family is the triplet count, the method is
`p_value_correction` — and the gate is significant when the corrected value is
at or below `alpha_dct`. That comparison is made in the run-wide pass. Before
it, the stream judges the gate provisionally to decide what else to measure:
exactly under `no`/`bfn`, whose correction follows from the family size alone,
and on the raw value under the rank-based methods, which is the conservative
side (see "Skipping a settled gate"). The count test itself is measured for
every triplet under every setting — it is the top of the cascade and a member
of a run-wide family.

### Gate 2 — Tree-height test (KS)

**Question.** Do the concordant and discordant1 tree-height distributions differ
at all — in any respect, not just in location?

`scipy.stats.ks_2samp(dis1_heights, con_heights, alternative="two-sided",
method="auto")` computes the two-sample Kolmogorov–Smirnov statistic: the
largest absolute gap between the two empirical cumulative distribution
functions, `D = sup_x |F_con(x) - F_dis1(x)|`. SciPy chooses an exact or
asymptotic p-value automatically based on the sample sizes. An empty sample
yields `(0.0, 1.0)`.

The KS test is deliberately omnidirectional. If the two height distributions are
indistinguishable, the discordant gene trees coalesce on the same timescale as
the concordant ones, which is what introgression between the *sampled* taxa
looks like — hence `inflow_introgression` when this gate is not significant.

Because the KS statistic responds to differences in shape, spread, and tails as
well as location, a significant result does not by itself say which direction
the heights moved. That is gate 3's job.

**When it runs, and what is compared to `alpha_ks`.** Like the count p-value,
the KS p-value is corrected across every triplet in the run, and the gate is
significant at `corrected <= alpha_ks`. The test is measured for every triplet
under `holm`, `fdr_bh` and `fdr_by`, whose corrections read every member's
value, and on a diagnostic run. Under `no` and `bfn` a non-diagnostic run
leaves it unmeasured below a count gate the stream has already found failed:
the family is still the triplet count, so no other member's corrected value
moves, and this member's own value would have gone unread. Such rows carry
empty `ks_*` columns, with `ks_significant` empty rather than false — see "The
two run-wide correction families".

### Gate 3 — Adaptive studentized permutation test

**Question.** Is the *mean* concordant height greater than, less than, or
indistinguishable from the mean discordant1 height?

The answer carries a p-value and a confidence statement, so a direction is
reported only when the separation is larger than sampling noise accounts for.
The test lives in `ghostparser/orchestrator/permutation.py`.

**When it runs.** The cascade reads the direction only when both earlier gates
are significant, and this test is much the most expensive of the three, so a
non-diagnostic run declines it for any triplet whose count or tree-height gate
the stream has already found failed. Those rows carry an empty `perm_*` block
and `perm_note = direction_test_not_consulted`. A diagnostic run measures it
for every triplet. The test also cannot run when either height sample is
empty: `perm_decision` is then `inconclusive`, the rest of the block stays
empty with no note, and `metrics.txt` counts those triplets separately (the
`n_con` and `n_dis1` columns show which). Neither omission moves any other
triplet's result, because the p-values here are corrected inside the test and
never across triplets — see "Why this correction is within-triplet only".

**The statistic.** For concordant sample `x` (size `nx`) and discordant1 sample
`y` (size `ny`), the Welch-studentized mean difference is

```
T = (mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)
```

with `var` the unbiased sample variance (`ddof=1`). The denominator is the Welch
standard error, following the unequal-variance standardization introduced by
Welch (1947, *Biometrika* 34(1–2), 28–35,
https://doi.org/10.1093/biomet/34.1-2.28). Using it rather than a pooled
standard error is essential here: the concordant sample is normally much larger
than the discordant1 sample and the two have different variances. Under that
combination a permutation test of the raw mean difference does *not* hold its
nominal level, whereas the studentized version remains asymptotically valid
(Janssen 1997, *Statistics & Probability Letters* 36(1), 9–21,
https://doi.org/10.1016/S0167-7152(97)00043-6). This is the permutation
analogue of the Behrens–Fisher problem.

**The null.** Pool all `nx + ny` heights and randomly reassign them to two
groups of the original sizes. Recompute `T` — including recomputing both
variances from the permuted groups, which is what preserves the studentization —
and repeat. The resulting distribution is the null distribution of `T` under
the hypothesis that group membership carries no information.

**p-values.** Two counts accumulate over the resamples: how many permuted
statistics are `>= T_obs` and how many are `<= T_obs`. Each becomes a
one-tailed p-value with the add-one estimator

```
p = (1 + count) / (1 + resamples)
```

which counts the observed arrangement itself. `perm_p_greater` is the upper
tail (the concordant mean is the larger) and `perm_p_less` the lower. The naive
`count / resamples` ratio can report exactly zero and understates the true
type-I error rate; the add-one form is the correctly-sized estimator for a
Monte Carlo permutation p-value (Phipson & Smyth 2010, *Statistical
Applications in Genetics and Molecular Biology* 9(1), Article 39,
https://doi.org/10.2202/1544-6115.1585).

**Correction inside the test.** The two one-tailed p-values form a testing
family of size two and are corrected against each other
(`correction.adjust_p_values`) with the configured `p_value_correction` method
before being compared to `alpha_perm`; under `no` the raw values are compared.
With the default `bfn` this compares `min(1, 2p)` to `alpha_perm`, which is the
conventional relationship between a two-sided level and its two one-sided
halves. A direction is called when exactly one corrected tail is at or below
`alpha_perm`. This correction is separate from and additional to the run-wide
correction applied across triplets, which covers only the DCT and KS p-values.

**Adaptive stopping.** A Monte Carlo p-value is an estimate, so the run keeps
resampling until the *decision* is safe rather than until a fixed budget is
spent. The first batch draws `min_resamples`. After each batch, a binomial
confidence interval at 95% (`ci_method`, default `wilson`) is placed around
each one-tailed p-value and rescaled onto the corrected scale. If `alpha_perm`
lies outside both intervals, no further resampling can flip the comparison and
the run stops with `perm_converged = True`. Otherwise the batch size grows by
25% and the run continues while the total is below `max_resamples` (raised to
`min_resamples` if configured smaller), after which it stops with
`perm_converged = False` and a `max_resamples_reached` note. Only the two
directional intervals take part in the stopping decision. `metrics.txt`
reports how many triplets spent the budget; which ones is in the results TSV's
`perm_converged` and `perm_note` columns.

`max_resamples` is the point at which the run stops asking for more, not a hard
cap on the total. The batch that crosses it is drawn at its full grown size
rather than trimmed to the remaining budget: sampling is vectorized, so a batch
costs the same per permutation however large it is, and the extra draws tighten
the interval that decides convergence instead of being spent on a stub that can
barely move it. `perm_n_resamples` can therefore exceed `max_resamples` by up to
one batch — about a quarter of the total when the budget spans several batches,
and more when `max_resamples` sits close to `min_resamples`, where a single
grown batch is comparable to the whole budget.

The confidence level is fixed at 95%. It governs how sure the stopping rule must
be before it commits — an internal precision knob rather than a statistical
choice the analysis turns on — so it is not exposed as a config key. It is also
a different quantity from `alpha_perm`, which is the threshold the p-value is
compared *against*.

#### How the interval is computed, and why it matches the p-value

The randomness in a Monte Carlo permutation test lives entirely in one place:
how many of the `n` drawn permutations landed beyond the observed statistic.
That count is `Binomial(n, p_true)`, where `p_true` is the exact permutation
p-value that full enumeration would give. Everything the stopping rule needs is
a statement about how far the count could be from `n · p_true`.

`statsmodels.stats.proportion.proportion_confint(count, nobs, alpha, method)`
answers exactly that. Its point estimate is `q = count / nobs`, and the default
`wilson` method inverts the **score test**: it returns every `p₀` for which the
score statistic stays inside the normal critical value,

```
|q - p₀| / sqrt(p₀(1 - p₀) / n)  ≤  z
```

Solving that quadratic in `p₀` gives the closed form statsmodels implements:

```
center = (q + z²/2n) / (1 + z²/n)
half   = z · sqrt( q(1-q)/n + z²/4n² ) / (1 + z²/n)
```

Note the interval is *not* centered on `q` — it is pulled toward ½ by the
`z²/2n` term, which is precisely why Wilson keeps close-to-nominal coverage near
0 and 1 where the plain normal ("Wald") interval fails and can even run outside
[0, 1] (Brown, Cai & DasGupta 2001, *Statistical Science* 16(2), 101–133,
https://doi.org/10.1214/ss/1009213286). That matters here because the p-values
this test produces are routinely near the floor. Clopper–Pearson (`beta`) is
also available and is guaranteed-coverage rather than approximate, but it is
conservative, so it would keep resampling past the point where the decision is
already settled.

**On method alignment.** The p-value and its interval are not two independent
estimates that might disagree — they are two summaries of the *same* pair of
numbers, `count` and `n`. GhostParser reports the add-one p-value
`(1 + count) / (1 + n)`, and passes `(count + 1, n + 1)` to
`proportion_confint`, whose own point estimate `count / nobs` then works out to
that identical value. So the interval brackets the quantity actually reported,
not one differing from it by `1/n`. The one approximation is that the interval
treats all `n + 1` arrangements as random when one of them — the observed
arrangement the add-one term accounts for — is fixed; this makes the interval
very slightly conservative, which is the safe direction for a stopping rule.

Because every supported correction is monotone in each p-value, the interval
is mapped onto the corrected scale by the same factor the point estimate
received before being compared against `alpha_perm`.

**The sampling optimization.** A naive implementation shuffles the pooled array
once per permutation in Python, which is far too slow to run inside every
bootstrap iteration. GhostParser instead draws a whole batch at once with
three changes:

1. *Only the smaller group is sampled; the larger one is derived exactly.* A
   permutation partitions the pooled values into two groups, so the larger group
   is precisely the complement of the smaller — nothing about it is unknown or
   estimated. The Welch statistic needs only a mean and a variance from each
   group, and both are functions of two running totals. Writing `S` and `Q` for
   the pooled sum and sum-of-squares (computed once, before the loop) and `s`
   and `q` for the drawn group's:

   ```
   S = sum(pooled)        Q = sum(pooled²)
   s = sum(drawn)         q = sum(drawn²)
   ```

   the drawn group of size `k` and its complement of size `m = n - k` have

   ```
   mean_drawn      = s / k
   var_drawn       = (q - k · mean_drawn²) / (k - 1)

   mean_complement = (S - s) / m
   var_complement  = ((Q - q) - m · mean_complement²) / (m - 1)
   ```

   The variance lines are the `E[X²] - E[X]²` identity in unbiased (`ddof=1`)
   form. This is an exact algebraic rearrangement, not an approximation: the
   recovered `var_complement` equals what `numpy.var(complement, ddof=1)` would
   return from the values themselves, to floating-point rounding. The larger
   group's values are therefore never gathered at all. Since the discordant1
   sample is usually the smaller one, the gathered data shrinks by roughly the
   size ratio.
2. *A partial partition replaces the shuffle.* Taking the `k` smallest of `n`
   uniform random keys yields a uniformly random size-`k` subset, and
   `numpy.argpartition` finds them in one linear pass rather than sorting the
   row or running Fisher–Yates over it.
3. *The pooled array is mean-centered once up front.* Centering leaves the mean
   difference and both variances unchanged, but it is what makes the subtraction
   above safe in floating point. `Q - q` is a difference of two positive sums of
   squares, and `m · mean²` is near zero once the pooled mean is zero, so
   neither step cancels significant digits. Without centering,
   `sum_of_squares - m · mean²` becomes a difference of two nearly equal large
   numbers and loses most of its precision — the same failure mode as the
   degenerate-scale guard, but silent.

The exactness of step 1 is not taken on trust: the test suite enumerates all
`C(9, 4) = 126` group assignments of a small case, evaluates each through the
plain scalar statistic, and requires the vectorized sampler to emit exactly that
set of values and nothing else.

Batches are chunked so the matrix of random keys stays near 16 MB, which matters
because every pool worker runs its own tests concurrently.

**Guards.** Four conditions short-circuit the test to `inconclusive` before any
resampling, each recorded in the `perm_note` column with the statistic and
p-values left empty and `perm_n_resamples` at `0`:

| Note | Condition | Why |
| --- | --- | --- |
| `insufficient_group_size` | Either group has fewer than 2 observations | No unbiased variance exists, so the statistic is undefined. |
| `zero_pooled_variance` | All pooled values are effectively identical | Nothing to detect and no scale to measure it on. |
| `degenerate_observed_scale` | Both groups internally constant, means differ | The statistic divides a real difference by numerical noise and reports the p-value floor regardless of how little data backs it. |
| `insufficient_permutation_support` | `C(n, k) < min_resamples`, checked when the pooled size `n` is at most 40 | The permutation distribution has fewer distinct values than the requested batch, so its resolution is capped well short of `alpha_perm`. |

The scale guards compare against a small fraction of the data's own magnitude
rather than against exact zero, because a sample of nominally identical values
such as `[0.9] * 10` has a floating-point variance around `1e-33`, not `0`.
Above a pooled size of 40 the support check is skipped, since `C(n, k)` is
then far beyond any resample budget.

Two further notes can appear on a test that did resample. `max_resamples_reached`
is the stopping rule running out of budget, described above.
`both_tails_significant` flags a state a coherent permutation distribution
cannot produce: the two one-tailed counts overlap on ties, so the two p-values
sum to more than one and cannot both sit at or below `alpha_perm` at once.
Reaching it means the inputs or the accumulators are inconsistent, so it is
surfaced rather than collapsed into a direction, and the test reports no
direction.

**Null skewness.** The direction comes from the pair of corrected one-tailed
tests, each of which is valid whatever shape the permutation null takes. The
shape is still worth knowing, so the sample skewness of the null is accumulated
across the batches and reported as `perm_null_skew`. It never drives a decision.

The value is a plain skewness of the drawn statistics: `0` for a symmetric null,
positive for a long right tail, negative for a long left tail. It is computed
from running power sums rather than by retaining the draws, so it costs nothing
in memory, and it is reported for **every** test that resampled — including the
`equivalent` and `inconclusive` ones, where the null can be just as asymmetric as
anywhere else.

What makes a null asymmetric here is a handful of extreme heights in the smaller
group. A permutation reassigns them, and the statistic depends sharply on *how
many* land in the small group, so the null separates into clusters rather than
forming one smooth curve. Take a 700-versus-19 split in which four heights sit
far above the rest: roughly nine permutations in ten leave all four in the
large group and land together on one side, about one in ten puts one of them
in the small group and pulls the statistic the other way, and a few in a
thousand put two there and pull it further. A strongly negative skewness is
that cluster structure showing up in the third moment.

Reading it: values near zero mean the null behaved like a symmetric reference
distribution and the p-values can be read at face value. Large magnitudes mean a
few gene trees dominate the smaller group's height distribution, so the group is
heterogeneous — a subset of gene trees far taller or shorter than the rest. That
is worth following up on the underlying alignments, but it does not invalidate
the directional call: the one-tailed permutation p-values are computed against
the true null whatever its shape. Tree heights are bounded below by zero and
routinely right-skewed, so mild asymmetry is the norm; `metrics.txt` reports the
median and maximum magnitude and how many triplets reach `|skew| >= 0.5`.

**Equivalence: separating "no difference found" from "shown to be the same".**
A non-significant directional pair means only that the data did not establish a
direction. It does not mean the two means are alike — that conclusion needs its
own test, because absence of evidence is not evidence of absence. When neither
tail is significant, the run therefore performs TOST — two one-sided tests, the
standard equivalence procedure. TOST turns the null around: instead of asking
whether the difference is zero, it states two nulls, that the difference is at
least a margin to one side of zero and that it is at least a margin to the
other, and tests each one-sidedly; rejecting both confines the difference to
within the margin. The margin is `EQUIVALENCE_DELTA = 0.5` pooled standard
deviations, and the outcome splits:

| `perm_decision` | Meaning |
| --- | --- |
| `equivalent` | Both shifted nulls rejected: the mean heights were *shown* to differ by less than the margin. |
| `inconclusive` | At least one was not: nothing was established in either direction. |

Both classify the triplet as `ambiguous`; the distinction is what the column
reports, not a fourth classification. `perm_p_tost` is populated only when the
step ran.

The two nulls are `mean(x) - mean(y) <= -delta` and `>= +delta`, with `delta`
the margin times the pooled standard deviation `sqrt((var(x) + var(y)) / 2)`.
Each is tested at its boundary by shifting the concordant sample: adding
`delta` to every `x` makes the two samples exchangeable under the lower null,
whose rejection region is the upper tail, so `p_lower` is the add-one share of
permuted statistics at or above the observed `T(x + delta, y)`; subtracting
`delta` mirrors it in the lower tail for `p_upper`. The TOST p-value is
`max(p_lower, p_upper)`, and `equivalent` needs it at or below `alpha_perm`. A
shifted observed statistic that is not finite leaves its own tail at 1, so it
cannot be rejected. This pair needs *no* multiplicity correction: rejecting a
union of nulls only when every component test rejects is an intersection-union
test, which holds its nominal level whenever its components do (Berger 1982,
*Technometrics* 24(4), 295–300, https://doi.org/10.2307/1267823). See
Schuirmann 1987 (*Journal of Pharmacokinetics and Biopharmaceutics* 15(6),
657–680, https://doi.org/10.1007/BF01068419) for the procedure and Lakens 2017
(*Social Psychological and Personality Science* 8(4), 355–362,
https://doi.org/10.1177/1948550617697177) for its use as a routine companion to
a significance test.

**The equivalence step rides on the directional draws.** Shifting the
concordant sample does not call for a fresh permutation pass. Let `v` be the
centered pooled vector with the concordant values in its first `nx` positions
(the x block), and `c` the shift. The shifted vector `w` has `w_i = v_i + c`
inside the x block and `w_i = v_i` elsewhere. For a drawn subset `S` of size
`k`, write `a` for how many of its entries came from the x block and `b` for
the unshifted sum of those entries. Then

```
sum_w(S)   = sum_v(S)   + c · a
sumsq_w(S) = sumsq_v(S) + 2c · b + c² · a
```

and the pooled totals move by constants, `S_w = S_v + c · nx` and
`Q_w = Q_v + 2c · sum(v over the x block) + c² · nx`. `a` and `b` are two extra
reductions over the gather the directional statistic already made — the drawn
indices below `nx` mark the x-block entries — so the studentization above runs
a second and third time on the adjusted sums, and every permutation yields
three statistics, unshifted, `+delta` and `-delta`, at a fraction of the cost of
one more pass. Because the studentized statistic is invariant to adding one
constant to every value, shifting the x block of the centered vector is the
same as permuting `(x + c, y)` outright, and the test suite checks that
identity against a literal shift of the data.

Reusing one permutation set across the three hypotheses is standard practice —
the maxT and minP procedures do it deliberately — and each p-value remains a
valid permutation p-value for its own null. The TOST counts take no part in
the stopping rule; they simply accumulate alongside, so the equivalence
p-values land at the same resample count, and the same `1 / (n + 1)`
resolution, as the directional pair.

**Why the margin is an effect size and not a number of standard errors.** The
margin is `0.5` *pooled standard deviations* — a Cohen's *d* of 0.5, the
conventional "medium" effect. Being dimensionless, it applies unchanged to every
triplet however large its tree heights are, which is what a per-triplet margin
has to do.

Expressing it in standard-error units would be dimensionless too, and it is the
obvious thing to reach for given that the test statistic is already studentized
— but it does not work. `T` is a pivot: its null distribution has a spread near
1 at *every* sample size. A margin of half a standard error would therefore sit
permanently inside the null's own scatter, and no quantity of data would move
it: equivalence would be unreachable by construction. The effect-size margin
shrinks relative to the standard error as gene trees accumulate, so more data
makes equivalence easier to establish, which is the behaviour the test needs.
The test suite pins both ends of that: two small samples drawn from one
distribution come back `inconclusive`, and two large ones `equivalent`.

The equivalence step runs only in the point estimate, and only when no direction
was found. Bootstrap iterations pass `equivalence_test=False` and skip it:
`equivalent` and `inconclusive` classify identically, so an iteration's vote
can never depend on which of the two it is.

**The interval on the studentized difference.** `bootstrap_perm_stat_ci_low` and
`bootstrap_perm_stat_ci_high` bracket `perm_statistic` — the observed `T` — at the
`1 - 2 * alpha_perm` percentile level of its bootstrap distribution. Each
bootstrap iteration recomputes `T` on its own resample of the gene trees
(`permutation.studentized_mean_diff`, whichever gate settled that iteration),
and the interval is the `alpha_perm` and `1 - alpha_perm` percentiles of the
finite values.

This is a different object from the `perm_p_*` intervals used by the stopping
rule, which are binomial intervals around a p-value. This one is an interval on
the *effect*, and it answers the question a p-value cannot: how large the
separation might plausibly be, in standard errors. The `1 - 2 * alpha_perm`
level is the one at which an interval and a one-sided test agree, so an interval
excluding zero corresponds to a directional rejection at `alpha_perm`. Reading
it alongside `perm_decision` distinguishes a direction that is well separated
from one that only just cleared the threshold. It is reported only when the
bootstrap is enabled and at least one iteration produced two observations in
both groups; otherwise both columns are empty.

**Known limitation.** The studentized permutation test is *asymptotically*
valid, not exact. Its type-I error rate sits close to `alpha_perm` across most
sample shapes — the test suite holds it there under unequal sizes and unequal
spreads — but it inflates when
the smaller group falls below roughly 30 observations *and* carries the larger
spread. Treat directional calls on triplets with very few discordant1 gene
trees as provisional; the `n_dis1` column is in the results TSV for exactly
this reason.

**Inside the bootstrap.** Each bootstrap iteration re-runs the cascade on its
resample, so the direction test runs there too — at one fifth of the configured
`min_resamples` and `max_resamples` (`permutation.bootstrap_resample_budget`,
floored at one), without the equivalence step, and only in iterations whose
count and tree-height gates both cleared, judged the way the stream judges the
point estimate's (exactly under `no`/`bfn`, on the raw value under a
rank-based method). Below a
failed gate the iteration's vote is already fixed, so no test is drawn; the
per-iteration `T` that feeds the studentized interval is computed either way.
The bootstrap aggregates many iterations into a single support value, which
absorbs the extra per-iteration Monte Carlo noise the reduced budget
introduces. `bootstrap_options.diagnostic` adds a direction test to every
iteration that would otherwise skip one, drawn from a separate seed child so
the votes and the interval do not move, and records all three tests per
iteration. How the iterations' own DCT and KS p-values are corrected is
described under [Correction inside the bootstrap](#correction-inside-the-bootstrap).

**Why this correction is within-triplet only.** The one-tailed pair is
corrected against itself and never across triplets. That is a hard constraint,
not a preference: a Monte Carlo p-value cannot fall below `1/(n_resamples + 1)`
(the add-one estimator's floor), which is `4.0e-4` at `min_resamples = 2500`.
Correcting across `n` triplets would require the raw p-value to clear
`alpha / (2n)` — already `2.5e-4` at only 100 triplets, below what the test can
express, so every triplet would come back `ambiguous` no matter how strong the
signal. At `C(83, 3) = 91,881` triplets the threshold is `2.7e-7` and would need
roughly 3.7 million resamples per triplet.

The DCT and KS p-values have no such floor because they are analytic rather than
sampled: a chi-square on counts `1000/50` gives `p = 6.2e-189`, and multiplying
by 91,881 still leaves `5.7e-184`. Correcting an exact p-value costs nothing;
correcting a sampled one spends resolution that had to be bought with compute.
That is why those two are corrected run-wide and this one is not.

A consequence worth stating plainly: gate 3 therefore carries no across-triplet
error control. Selecting triplets on gates 1 and 2 and then testing gate 3 at
`alpha` is post-selection inference and does not inherit the earlier gates'
control, so a direction call is conditional on that selection — descriptive
rather than confirmatory.

**What a direction means.** Every triplet reaching gate 3 is decided here, so
an `outflow` or `ghost` call always rests on a corrected one-tailed p-value at
or below `alpha_perm`. A separation the data cannot resolve at that threshold
is reported `equivalent` or `inconclusive` and classified `ambiguous`.

### The two run-wide correction families

The DCT and KS p-values are each corrected once across triplets, as two
separate families that both contain **every triplet**. Each family is simply
that test's column of raw p-values, so a corrected value is a function of its
own column and the method — the count test's outcome plays no part in how the
tree-height p-values are corrected, and vice versa. The only place the two
meet is the cascade, which reads the corrected flags in order and stops at the
first gate that settles the call.

The family size is the triplet count, which is known before the first triplet
is measured; that is what lets the inline methods below vote as they go, and
what lets a non-diagnostic run leave a member unmeasured without touching the
rest. Under `no` and `bfn` each member's corrected value follows from its own
raw value and the count alone, so a triplet whose count gate failed can skip
the tree-height test: it stays a member — the measured members are still
corrected by the full triplet count, never by the number measured — and its
own corrected value would have gone unread. Under `holm`, `fdr_bh` and
`fdr_by` every member's value moves the others' ranks, so the tree-height test
is measured for every triplet whether or not the run is diagnostic. Every
triplet reports `dct_p_value_corrected` and `dct_significant`;
`ks_p_value_corrected` and `ks_significant` are empty only on a triplet that
never measured the test.

### Correction inside the bootstrap

`bootstrap_value` is only meaningful if the iterations answer to the same
decision rule the reported classification does. The point estimate compares
*corrected* p-values to `alpha_dct` and `alpha_ks`, so every bootstrap iteration
must too — judging iterations on raw p-values while reporting a classification
made on corrected ones measures two different rules and produces rows whose
support contradicts their own classification.

The obstacle is that a correction is a property of a *family*, not of a single
p-value, and the family here spans triplets: iteration `i` of triplet A belongs
with iteration `i` of every other triplet — one column of the grid is one
resampled run. A streaming engine that finishes one triplet before starting the
next does not have the rest of the column in hand. Two tiers handle this:

- **Inline** (`no`, `bfn`): the multiplier follows from the family size alone,
  and that is the triplet count, known before the stream starts. Each iteration
  classifies itself on the spot with exactly the corrected value the run-wide
  pass would produce, and only the tally leaves the loop. Nothing per
  iteration is stored, which is what keeps a 300,000-triplet run's memory
  flat.
- **Deferred** (`holm`, `fdr_bh`, `fdr_by`): a rank-based multiplier depends on
  every other member's value. Each iteration parks its raw DCT and KS p-values
  and its direction code in a `DeferredBootstrapRecord` (about 17 bytes per
  iteration per triplet); the deferred pass then corrects each column across
  every triplet at once, classifies the whole grid, and tallies
  each triplet's votes.

The tree-height test is measured wherever the correction will read it, and
nowhere else — `diagnostic` reaches only the point estimate. Under a rank-based
method that is every iteration, because every member's value moves the others'
ranks. Under an inline method a value below a failed count gate is never read —
the family size is fixed — so the iteration leaves it unmeasured and votes
`no_introgression` directly. The bootstrap's own
switch is `bootstrap_options.diagnostic`, which measures all three tests in
every iteration and writes them per iteration. The direction tests it adds
below a failed gate draw from a fifth child of the triplet's seed sequence, so
the vote's own tests draw exactly what they would have without the record:
`bootstrap_value` and the studentized interval do not move when it is on.

### Skipping a settled gate

Once a raw gate has failed, the cascade's answer is already fixed: `p_raw > alpha`
implies `p_adjusted > alpha` for every supported correction, so the corrected gate
fails too and nothing below it can change the classification. The direction
test below such a gate can therefore be skipped outright — and it is much the
most expensive of the three. This licenses both short-circuits: a bootstrap
iteration skipping its own direction test, and a non-diagnostic point estimate
skipping a triplet's. In both places the gate is judged the same way: under an
inline method on the exactly corrected value, under a rank-based method on the
raw one, which is the conservative side.

The monotonicity it rests on holds for every supported method by construction:

- **Bonferroni** multiplies by the family size `n >= 1`, so `p_adj = min(1, np) >= p`.
- **Holm** (Holm 1979, *Scandinavian Journal of Statistics* 6(2), 65–70,
  https://www.jstor.org/stable/4615733) adjusts the `j`-th smallest by
  `(n - j + 1) · p_(j)` under a running maximum; the multiplier is at least 1 for
  every `j <= n`.
- **Benjamini–Hochberg** (Benjamini & Hochberg 1995, *JRSS B* 57(1), 289–300,
  https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) uses `n/j · p_(j)` under a
  running minimum, and `n/j >= 1` for every `j <= n`.
- **Benjamini–Yekutieli** (Benjamini & Yekutieli 2001, *Annals of Statistics*
  29(4), 1165–1188, https://doi.org/10.1214/aos/1013699998) multiplies BH by the
  harmonic factor `sum(1/i) >= 1`, so it is uniformly larger still.

This is a constraint on the supported set, not an accident of it. The property
is easy to lose: a procedure that *estimates* the number of true null hypotheses
`n₀ <= n` and substitutes it for `n` has a multiplier `n₀/j` that can fall below
1, so an adjusted p-value can land beneath its raw one. Adding such a method
would silently break both short-circuits — every skipped test would become a
test that might have changed the answer — so monotonicity is a precondition
for anything entering `P_VALUE_CORRECTION_CHOICES`. The test suite asserts it
over the whole choice list rather than a fixed set of names, so a method that
violates it fails immediately on being added.

The family size is the triplet count regardless of any skipping, so the
correction never depends on the optimization. Whichever tier applies, the
resample stream is untouched: the bootstrap draws its resamples from a generator
independent of the permutation tests', so a fixed `seed` reproduces the
same resamples — and the same `bootstrap_perm_stat_ci_*` interval — under every
correction method.

### The point estimate's short-circuit

With `diagnostic` off (the default) the same argument applies to the point
estimate: the per-triplet pass runs the direction test only when both earlier
gates cleared, and under an inline correction runs the tree-height test only
when the count gate cleared. `diagnostic: true` runs every test for every
triplet. The bootstrap iterations take the lean path either way: the switch
reaches the point estimate alone, and every iteration judges its gates the
same way whether or not the run is diagnostic.

Expect the default to be somewhat faster than a diagnostic run, not
dramatically so. What it declines is the point estimate's direction test and,
under `no`/`bfn`, the point estimate's tree-height test, while most of a run's
time is the bootstrap's permutation tests: every iteration of a triplet that
reaches gate 3 runs one either way. The wall time moves less than the CPU
time, because the triplets that reach gate 3 set the length of the run.

Skipping the direction test never changes a decision, under any correction
method, for two reasons that have to hold together:

- **Monotonicity**, above: a raw-failed gate cannot clear once corrected, so a
  skipped triplet could never have reached gate 3. The cascade reads
  `perm_decision` only when *both* corrected gates are significant, and a
  skipped triplet has at least one gate that already failed on the value the
  correction can only raise. Its own classification, gate and
  bootstrap votes are therefore fixed before the direction test would run.
- **Permutation p-values are corrected within the test only**, across its pair
  of one-tailed p-values, never across triplets (see "Why this correction is
  within-triplet only" under gate 3). Omitting one triplet's direction test
  therefore changes nothing for any *other* triplet: there is no run-wide
  permutation family for it to have been a member of. That is what makes the
  skip safe to take triplet by triplet while the stream is still running,
  before the run-wide correction pass has seen the whole family.

The DCT and KS p-values are different precisely on the second point — their
families span the run — which is why the count test is never skipped and the
tree-height test only under a correction whose family size alone corrects it.
Under `no`/`bfn` the gate is judged on the exactly corrected value
(`p × n_triplets`), so the point estimate skips exactly the tests the pass
would find unread; under a rank-based method it is judged on the raw value, so
it may run a direction test the pass then ignores, but never skips one it
needs.

The run-wide pass replays the stored `perm_decision` after correction rather
than recomputing it. A skipped test leaves that field
`None`, while a guard or an empty group still records a string, so reaching gate
3 without a decision is an invariant violation rather than an ambiguous call;
the pass raises instead of silently classifying such a triplet `ambiguous`. It
likewise raises if a rank-based family arrives with members missing.

Skipped triplets report `perm_note = direction_test_not_consulted` and leave the
rest of the `perm_*` block empty, which distinguishes a deliberate skip from a
test that ran and hit a guard. `metrics.txt` reports the setting, spells out
what it skips under the run's correction, and counts the skipped tree-height and
direction tests alongside the triplets clearing each gate.

## Shape diagnostics

Off by default, enabled with `shape_diagnostics: true`. Fifteen columns
describing the shape of each height group — `con_*`, `dis1_*`, `dis2_*`, the
same three groups the tests are built on. Nothing here feeds a classification;
they exist to say what the height distributions actually look like.

They are written to the results TSV only. `summary_statistics.tsv` is consumed
as a feature matrix, and several of these columns are undefined below their
observation floors (see the thresholds below), so carrying them there would
leave a hole in every row that hit one.

| Column | Meaning |
| --- | --- |
| `<group>_n_modes` | Peaks in a Gaussian KDE at Scott's bandwidth. |
| `<group>_modes_p` | Silverman p-value for "the density is unimodal". |
| `<group>_skew` | Sample skewness. |
| `<group>_excess_kurtosis` | Sample excess kurtosis (`0` for a normal). |
| `<group>_tail_xi` | Generalized-Pareto shape of the upper decile. |

**Counting peaks.** `n_modes` is a raw count at one bandwidth, and a KDE mode
count is a property of the bandwidth as much as of the data: a strongly skewed
but unimodal sample, an exponential or a Pareto, reads as several modes at the
default bandwidth. Read it as a description of the smoothed density, not as a
number of components.

Two unrelated resampling procedures are both called a bootstrap in this
codebase. The **run bootstrap** (`bootstrap_options.iterations`) resamples gene
trees to produce `bootstrap_value`, and never touches these diagnostics. The
**smoothed bootstrap** below is internal to the modality test: it resamples from
a smoothed density to calibrate one p-value, and is what
`MODALITY_BOOTSTRAP_RESAMPLES` counts.

`modes_p` is the inferential column. It is Silverman's critical-bandwidth
bootstrap test (Silverman 1981, *Journal of the Royal Statistical Society B*
43(1), 97–99, https://doi.org/10.1111/j.2517-6161.1981.tb01155.x): find the
smallest bandwidth `h` whose KDE has one mode, then resample from the density
smoothed at exactly that `h` and count how often a resample still needs more
smoothing. A `p` at or below `alpha` says the observed `h` was implausibly large
for a unimodal density — evidence of a second mode. Each replicate costs one
mode count rather than a second bisection, because mode count is monotone in
bandwidth: a resample needs a larger critical bandwidth exactly when it still
has too many modes at `h`.

It holds its level on strongly skewed unimodal shapes — normal, lognormal,
exponential, gamma and Pareto samples are all left unrejected, the case a mode
count or a BIC-selected Gaussian mixture gets wrong — at the price of needing
roughly three standard deviations of separation between two components before
it sees them; a mixture two standard deviations apart is missed. A
non-significant `modes_p` is therefore weak evidence of unimodality and a
significant one is strong evidence against it.

One known failure: a uniform density has no interior mode at all, and its hard
edges cannot be reproduced by a Gaussian smoothed bootstrap, so the test rejects
it. Tree heights are not uniform, but a group whose heights are near-constant
across an interval will read as multimodal.

**Tail weight.** `tail_xi` fits `scipy.stats.genpareto` to the exceedances over
the group's 90th percentile, with the location pinned at the threshold. The
shape parameter is the tail index: positive is a power-law tail whose moments
above order `1 / xi` do not exist, `0` is exponential decay, negative is a tail
with a finite upper endpoint at `threshold - sigma / xi`. Peaks-over-threshold is
the estimator to use here rather than Hill's, which is non-negative by
construction — it cannot represent a bounded tail at all, and it reads clearly
positive for an exponential and a lognormal, both of which have `xi = 0`.

**Guards.** Groups with fewer than `SHAPE_MIN_OBSERVATIONS` (20) observations,
or with no spread, leave all five columns empty. `tail_xi` additionally needs
`SHAPE_MIN_TAIL_EXCEEDANCES` (10) points above the threshold, so it stays empty
below 100 observations even when the moments are reported. Below a few hundred
observations the tail index is noisy enough that its sign is the only part worth
reading.

**Where they appear.** The results TSV carries them as `con_*`/`dis1_*`/`dis2_*`,
matching its other per-group naming, and nowhere else: `summary_statistics.tsv`
is a feature matrix, and these columns are undefined for groups below their
observation floors, so including them there would leave holes in it. The
columns are derived from what the run measured rather than from the setting,
so the file cannot advertise a measurement that did not happen.

**Cost.** The modality bootstrap is the expensive part: it resamples once per
height group, three groups per triplet. That is why the default is off —
across a large taxon set it dominates everything else the run does. The
diagnostics are measured once per triplet, in the point estimate, from the
observed heights. Bootstrap iterations never recompute them: an iteration
resamples the gene trees, so its groups are a different sample, and measuring
each one would multiply the cost by the iteration count.

## Outputs

The output folder is reset before the run under the default `overwrite: true`,
so it must be a directory of its own. Pointing it at a directory that holds
input trees deletes them. Consolidation resets its own `consolidation/`
subfolder the same way, which is why it gets a subfolder rather than writing
beside the results TSV.

Written under the output folder:

- `orchestrator_triplet_results.tsv` — one row per triplet.
- `summary_statistics.tsv` — only when `generate_summary_stats` is set;
  per-triplet topology/metric summary statistics (63 metric columns covering
  mean/median/mode/variance/entropy/min/max over avg-tree-height/internal-branch/
  sister-distance for concordant/discordant1/discordant2). The `discordant1_*`
  columns describe whichever discordant topology is more frequent — the same
  group named by the `dis1_topology` column and used by all three tests — and
  `discordant2_*` the other one. The roles are resolved per triplet from the
  observed counts, not fixed to a topology label, so a triplet where `AC|B`
  outnumbers `BC|A` has its `AC|B` gene trees under `discordant1_*`. The shape
  diagnostics are not repeated here — see below.
- `processed_<species tree>` / `processed_<gene trees>` — cleaned, rooted trees,
  in the input trees' own labels whatever `species_rename_map` says. A label a
  bare Newick token cannot hold (a space, a quote, a bracket) is written
  single-quoted.
- `metrics.txt` — per-stage wall/CPU timing and run parameters.
- `consolidation/` — the combined heatmap/bar-chart plot, with the TSV
  matrices under `consolidation/consolidation_data/`, unless
  `--no-consolidation` is given. Consolidation writes into this dedicated
  subfolder so its own output-directory reset never
  removes the run folder's results TSV, processed trees, or the open
  `metrics.txt`.

### How each results column is produced

| Column | Source | Method |
| --- | --- | --- |
| `triplet` | Triplet setup | The normalized `(A, B, C)` labels, A and B being the species-tree sisters. |
| `abc_mapping` | Triplet setup | The same labels spelled out as `A=<taxon>;B=<taxon>;C=<taxon>`. |
| `species_tree` | Triplet setup | The triplet's species subtree, serialized topology-only (branch lengths omitted). |
| `dis1_topology` | Topology ranking | `BC` or `AC` — whichever discordant topology is more frequent; ties resolve to the first listed. |
| `most_frequent_matches_concordant` | Topology counts | True when the concordant count is at least both discordant counts. |
| `n_con` / `n_dis1` / `n_dis2` | Topology counts | Gene trees observed with each topology. |
| `analyzed_trees` | Extraction | Gene trees from which a subtree for this triplet was extracted. |
| `dct_chi_stats` or `dct_z_score` / `dct_p_value` | DCT | SciPy chi-square or statsmodels z-test over `[n_dis1, n_dis2]`; the statistic column is named after the backend. An all-zero discordant split short-circuits to `(0.0, 1.0)`. |
| `dct_p_val_<method>_corr` | Correction | Run-wide correction over every triplet's DCT p-value. Omitted under `no`. |
| `dct_significant` | Decision gate 1 | Corrected DCT p-value at or below `alpha_dct`. |
| `ks_statistic` / `ks_p_value` | Tree-height test | Two-sample KS between concordant and discordant1 heights; an empty sample yields `(0.0, 1.0)`. Empty under `no`/`bfn` when the count gate settled the triplet, unless the run is diagnostic. |
| `ks_p_val_<method>_corr` | Correction | Run-wide correction over every triplet's KS p-value, by the triplet count. Omitted under `no`; empty where the test was not measured. |
| `ks_significant` | Decision gate 2 | Corrected KS p-value at or below `alpha_ks`. Empty where the test was not measured. |
| `perm_statistic` | Direction test | Observed Welch-studentized mean difference. Empty when a guard fired. |
| `perm_p_greater` / `perm_p_less` | Direction test | Raw one-tailed p-values, add-one estimator. |
| `perm_p_greater_<method>_corr` / `perm_p_less_<method>_corr` | Direction test | The one-tailed p-values corrected against each other, and the values compared to `alpha_perm`. Omitted under `no`. |
| `perm_p_tost` | Equivalence test | TOST (two one-sided tests) p-value: `max(p_lower, p_upper)` over the two shifted-null tests, one per side of the equivalence margin, both read off the direction test's own permutations. At or below `alpha_perm` the two mean heights were shown to differ by less than half a pooled standard deviation, which is what makes `perm_decision` read `equivalent` rather than `inconclusive`. Populated only when neither direction was significant. |
| `bootstrap_perm_stat_ci_low` / `bootstrap_perm_stat_ci_high` | Bootstrap | Percentile interval on the studentized difference at the `1 - 2 * alpha_perm` level. Empty without bootstrap, or when no iteration had two observations in both groups. |
| `perm_n_resamples` | Direction test | Permutations drawn; `0` when a guard fired. Can exceed `max_resamples` by up to one batch, since the final batch is not trimmed. |
| `perm_converged` | Direction test | True when the confidence interval excluded `alpha_perm` before the budget ran out. |
| `perm_null_skew` | Direction test | Sample skewness of the permutation null: the third standardized moment of the `perm_n_resamples` studentized statistics drawn while testing this triplet. It describes the *reference distribution the test built*, not the tree heights themselves. `0` is a symmetric null and the p-values behave like a textbook two-sample test; a large magnitude means a few extreme heights in the smaller group dominate the resampling, so the null breaks into clusters by how many of them land where, and the sign names the long tail. Reported for every test that resampled, including `equivalent` and `inconclusive` ones. Never consulted by any decision — see "Null skewness" above for a worked example. |
| `perm_note` | Direction test | Guard slug, `max_resamples_reached`, or `direction_test_not_consulted` when a non-diagnostic run skipped a test an earlier gate had settled; empty on a clean run that resampled. |
| `perm_decision` | Decision gate 3 | `greater`, `less`, `equivalent`, or `inconclusive` for concordant relative to discordant1. Consulted only when `decision_gate` is `PERM`, and populated only there unless the run is diagnostic, which populates it for every triplet. |
| `decision_gate` | Decision logic | Which test settled the classification: `DCT`, `THT`, or `PERM`. |
| `classification` | Decision logic | `no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, or `ambiguous`. |
| `inference` | Reporting | Human-readable direction naming the actual species. |
| `bootstrap_value` / `all_bootstrap` | Bootstrap | Fraction of iterations agreeing with the final classification, plus the full class-fraction map. Iterations are judged against the same corrected thresholds as the point estimate. Present unless `--no-bootstrap`. |
| `con_*` / `dis1_*` / `dis2_*` shape columns | Shape diagnostics | Mode count, Silverman modality p-value, skewness, excess kurtosis and generalized-Pareto tail index per height group. Present only with `shape_diagnostics`; see above for how to read each. |
| `bootstrap_*` diagnostic columns | Bootstrap diagnostic | Per-iteration DCT and KS statistics and p-values, the direction test's studentized statistic, raw one-tailed p-values and decision, con/dis means, and the gene-tree heights. Present only with `bootstrap_diagnostic`; summaries instead of lists with `bootstrap_summary_only`. |

This file reports the *tests*: the direction is read off `perm_p_greater` and
`perm_p_less` after correction, so it carries no per-group mean or median
columns. Descriptive per-group statistics live in `summary_statistics.tsv`,
which `generate_summary_stats` enables, under their own per-topology headers
(`concordant_avg_tree_height_mean`, `discordant1_avg_tree_height_median`, and
so on).

Note that the two files measure different things and will not agree numerically:
the `*_avg_tree_height_*` summary columns always average the three root-to-tip
distances, whereas the result fields average H(T) as selected by
`tree_height_calculation_strategy`. They coincide only under the default `AVG`.

The results TSV is named `orchestrator_triplet_results.tsv`. Consolidation
writes its own artifacts into a `consolidation/` subfolder, so the two never
collide in the run's output folder.

## Parallelization

Triplets are split into chunks and dispatched across `--processes` workers; each
worker runs the fused extract-then-infer loop for its chunk over the shared
gene-tree list and the shared geometry cache. Triplet chunks are the only unit
worth distributing, because a run's cost is per-triplet inference.

`processes: 0` resolves to the CPUs the process may actually run on
(`available_cpu_count`), not the machine's core count. Under a container
limit or a job scheduler the two can differ widely, and a worker per machine
core would then time-slice many processes over few CPUs while each carries
its own per-chunk state. `metrics.txt` reports the resolved count as
`Worker processes`.

Within a chunk, a worker extracts observations for a fixed batch of triplets
at a time and analyzes them before extracting the next batch. One
observation is a tuple and a float per gene tree, so a triplet's list grows
with the gene-tree count and a whole chunk's worth would be far larger than a
batch. The results, which are small and slotted, still accumulate for the
chunk and return to the parent in one payload.

With one worker (or one triplet) the engine runs the fused loop serially in the
parent process.

Everything above happens on one machine. The workers are forked from the
parent, share its geometry cache copy-on-write and return their results over
pipes, none of which crosses a machine boundary. On a cluster, a job spanning
several nodes therefore runs the orchestrator on the node the script started
on and leaves the others idle: ask for the CPUs and memory on a single node.
No multi-node execution is implemented. Memory is not divided among workers
either -- a job's memory limit is a ceiling on the sum over the parent and
every worker, and the peak is the parent's geometry cache plus
the accumulated results (held twice during the decision pass) plus each
worker's small private state.

## Internal design

### File layout

```
ghostparser/orchestrator/
  __init__.py    exports run_orchestrator
  __main__.py    python -m ghostparser.orchestrator entry point: main() wires parsing -> run_orchestrator
  config.py      orchestrator defaults/choices, validation, CLI parser, and CLI/config resolution
  correction.py  multiple-testing correction shared by inference.py and permutation.py
  permutation.py adaptive studentized permutation test (decision gate 3)
  shape.py       optional per-group modality/skew/tail diagnostics
  trees.py       tree/triplet preprocessing
  triplet_geometry.py  per-tree LCA/parent tables and the triplet read-out over them
  inference.py   per-triplet inference + summary stats + result type + TSV writers
  stream.py      fused extract+infer streaming engine
  preflight.py   structural data check reached via --preflight-data-check
  consolidation.py  introgression maps and TSV matrices (the final stage)
  runner.py      run_orchestrator coordinator (cleaning -> streaming -> correction -> writing -> consolidation)
  ORCHESTRATOR.md    this document
```

### Shared dependencies

The orchestrator owns its tree preprocessing, inference, consolidation, and
configuration. It imports only three things from the rest of the package:

- `ghostparser.config` — the shared configuration trunk (`ConfigError`, path
  resolution, raw config-file loading, required-path validation, overwrite
  resolution, `prepare_output_directory`). Orchestrator-specific defaults, choices,
  and validators live in `orchestrator/config.py`, which is why the orchestrator can set
  its own defaults without affecting the ML subpackage.
- `ghostparser.cli_config` — the generic `resolve_cli_or_config_args` resolver
  implementing config-file-wins precedence.
- `ghostparser.triplet_utils` — pure topology helpers.

### Fused streaming engine

The per-triplet unit is `inference.analyze_triplet_from_observations(triplet,
observations, species_subtree, ...)`, which takes precomputed `(topology,
tree-height, metrics)` observations. The third element carries the per-tree
summary metrics and is `None` unless `generate_summary_stats` is enabled.

Every gene tree is parsed and cached once per run, before any triplet is
touched: `stream.build_run_geometry` walks each tree twice — pre-order for a
parent/edge-length array, post-order for a pairwise LCA table — and keeps only
those arrays. The processing unit is then a chunk of triplets read out of that
cache. For each gene tree the engine looks up the triplet's three pairwise LCAs
and sums a few short paths up the parent chain
(`triplet_geometry.geometry_observation`) — no subtree is copied, and there is
**no** serialize-to-Newick and reparse round trip. It then runs
`analyze_triplet_from_observations` for each triplet in the chunk and drops the
observations. The cache is built once in the parent and shared read-only to
workers via a fork/forkserver initializer, so a triplet chunk costs table
lookups rather than tree copies.

`trees.extract_triplet_subtree` and `inference.observation_from_subtree` still
exist as the DendroPy reference implementation of that read-out. Nothing in a
run calls them; the parity tests hold the cached path to agreement with them.

Bootstrap resampling is vectorized with NumPy: per-triplet resample indices are
drawn with a seeded `numpy.random.Generator`, and topology counts and
per-topology height groups are computed with array operations. Bootstrap values
are deterministic under a fixed `seed` — the per-triplet seed is
derived from the run seed and the triplet, so any worker count agrees exactly.

### Memory rationale

- No intermediate triplet-gene-trees file is written or reloaded.
- Peak memory is roughly the geometry cache (one `array` LCA table per gene
  tree, shared copy-on-write with the workers), plus one extraction batch of
  observations per worker, plus the accumulating list of result objects.
- Result objects must be accumulated because global p-value correction needs all
  p-values in a single pass. `TripletPipelineResult` is a slotted dataclass for
  that reason: a per-instance dict would be a large share of each result, and
  the decision pass holds the measured and the decided lists at once.

### Consolidation interface contract

`generate_introgression_maps` reads only four fields from each result
(`triplet`, `classification`, `dis1_topology`, `bootstrap_value`), read
duck-typed. `TripletPipelineResult` keeps those four
fields; treat field parity on them as a maintenance constraint. The taxon
names in `results`, `plot_taxa` and `outgroups` must agree with each other;
`rename_map` maps the species tree file's labels onto them when they differ,
which is the case for a run with a `species_rename_map`.
