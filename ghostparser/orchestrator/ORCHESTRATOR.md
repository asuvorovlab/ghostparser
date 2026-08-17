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
`generate_summary_stats`, and the `bootstrap_options` and `permutation_options`
blocks). One default is worth calling out: `p_value_correction` defaults to
`bfn` (Bonferroni), and it is applied both run-wide across triplets and inside
each permutation test across its pair of one-tailed p-values.

See [CONFIG.md](../../CONFIG.md#orchestrator-primary-module) for every key, its
default, and its allowed values.

## Preflight data check

`--preflight-data-check` (or `preflight_data_check: true`) turns the run into a
data validation pass. `runner._run_preflight_only` short-circuits
`run_orchestrator` immediately after the output directory is prepared, so no
analysis runs and the only artifact is `preflight_data_check.txt`.

`preflight.run_preflight_data_check` replays the same structural logic the
engine uses, but collects failures instead of raising on the first one:

1. Root the species tree on the outgroup and normalize each triplet to A/B/C
   via `find_sister_pair` and `normalize_abc_from_sister_pair`, recording
   triplets whose rooted topology cannot be resolved.
2. Root every gene tree with `trees._root_tree_on_any_outgroup`, recording the
   ones where no outgroup label is present.
3. For each triplet contained in a gene tree, extract the subtree and replay
   `triplet_taxa_labels` → `find_sister_pair` → `topology_from_sister_pair` →
   sister-pair MRCA lookup, recording whichever step fails.

Every failure becomes an `Issue` with a dotted category. The report groups them
by category with counts, up to 25 examples each naming the gene-tree index and
triplet plus the offending input line, and an attribution summary separating
species-tree causes from gene-tree causes. `PreflightResult.passed` is `True`
only when nothing was detected.

Three conditions make the check itself impossible and raise `ValueError`
instead: no outgroups given, a species-tree file that does not hold exactly one
tree, and a species tree containing none of the outgroups.

## Orchestrator summary

`runner.run_orchestrator(config)` coordinates the run:

1. **Species preprocessing** — `trees.clean_and_save_trees` standardizes the
   species tree and drops trees whose mean internal support is below
   `min_support_value`. `trees._root_tree_on_outgroup` roots on the outgroup
   MRCA and prunes the outgroup, returning the ingroup taxa.
2. **Triplet setup** — `trees.generate_triplets` enumerates every ingroup
   triplet (or `trees.read_triplet_filter_file` plus
   `trees.filter_triplets_by_taxa` restricts them).
   `trees._build_species_triplet_metadata` normalizes each triplet to
   `(A, B, C)` with A and B the species-tree sisters, and builds the triplet's
   species subtree.
3. **Gene-tree preprocessing** — `trees.clean_and_save_gene_trees` cleans each
   gene tree and roots it on the outgroup.
4. **Fused extraction + inference** — `stream.stream_triplet_results` walks the
   triplets, extracts each one's subtree per gene tree, converts it directly to
   an observation (`inference.observation_from_subtree`), and immediately runs
   `inference.analyze_triplet_from_observations`. Only the small result object
   is retained; the subtrees are discarded.
5. **Run-wide correction** —
   `inference._apply_triplet_result_p_value_correction` applies the
   multiple-testing correction once across all triplets, because a global
   correction needs every p-value in a single pass.
6. **Writing** — `inference.write_pipeline_results` emits
   `orchestrator_triplet_results.tsv`; `inference.write_summary_statistics_tsv`
   emits `summary_statistics.tsv` when `generate_summary_stats` is set.
7. **Consolidation** — `introgression_mapper.generate_introgression_maps`
   writes the map artifacts into a `consolidation/` subfolder.

## Per-triplet inference

For each triplet the engine classifies every gene tree's subtree into one of
three topologies — concordant (matching the species tree) plus two discordant
alternatives — and records a tree height H(T) per the configured strategy. It
then applies a three-gate decision:

1. **Discordant count test (DCT)** — compares the two discordant counts
   (`inference.run_discordant_count_test`, chi-square or z-test). If the
   corrected p-value is not below `alpha_dct`, the triplet is
   `no_introgression` and the remaining gates are skipped.
2. **Tree-height test (THT)** — a two-sample KS test between the concordant and
   discordant1 height distributions (`inference.run_two_sample_ks_test`). If it
   is *not* significant, the triplet is `inflow_introgression`.
3. **Direction test** — otherwise the concordant and discordant1 heights are
   compared directionally by the studentized permutation test
   (`permutation.run_studentized_permutation_test`): `greater` gives
   `outflow_introgression`, `less` gives `ghost_introgression`, and no
   resolvable direction gives `ambiguous`.

`inference._classify_introgression` implements this decision table directly,
returning the classification and the name of the test that settled it as one
pair. That pair populates the `classification` and `decision_gate` columns, so
the two are derived in a single pass and cannot drift apart.

The cascade *consults* the tests in order, but the point estimate *computes* all
three for every triplet, so the `perm_*` columns are populated even where the
cascade never reached them. Read `decision_gate` before reading `perm_decision`:
only `Permutation` means the direction result produced the classification. A row
carrying `decision_gate = THT`, `perm_decision = ambiguous`, and
`classification = inflow_introgression` is consistent — the direction test ran
and was recorded, but the tree-height test had already settled the call.

Bootstrap resampling (on by default) repeats the analysis over resampled
observations and aggregates the per-iteration classifications into
`bootstrap_value`.

## The statistical tests

Each gate answers a different question, and each is computed by a named library
routine rather than by hand. This section states what each test measures, how
its value is obtained, and where its assumptions bite.

### Gate 1 — Discordant count test

**Question.** Are the two discordant topologies equally frequent?

Under incomplete lineage sorting alone, the two discordant histories are
exchangeable and should appear about equally often; an excess of one of them is
the signal that something other than ILS — introgression — has acted (Huson et
al. 2005, *RECOMB*, https://doi.org/10.1007/11415770_18). The test therefore
asks only whether `n_dis1` and `n_dis2` depart from a 50/50 split, and ignores
the concordant count entirely.

Two backends, selected by `discordant_test`:

- `chi-square` (default) — `scipy.stats.chisquare([n_dis1, n_dis2])`, a
  goodness-of-fit test against equal expected counts. The statistic is
  `sum((observed - expected)^2 / expected)` with `expected = (n_dis1 + n_dis2) / 2`,
  compared against a chi-square distribution on one degree of freedom.
- `z-test` — `statsmodels.stats.proportion.proportions_ztest` with
  `count=[n_dis1, n_dis2]`, `nobs=[total, total]`, `alternative="two-sided"`,
  a two-proportion z-test on the same counts.

A zero/zero split short-circuits to `(0.0, 1.0)` rather than dividing by zero.
Both backends test the same null and agree closely; the chi-square statistic is
approximately the square of the z-score.

### Gate 2 — Tree-height test (KS)

**Question.** Do the concordant and discordant1 tree-height distributions differ
at all — in any respect, not just in location?

`scipy.stats.ks_2samp(dis1_heights, con_heights, alternative="two-sided",
method="auto")` computes the two-sample Kolmogorov–Smirnov statistic: the
largest absolute gap between the two empirical cumulative distribution
functions, `D = sup_x |F_con(x) - F_dis1(x)|`. SciPy chooses an exact or
asymptotic p-value automatically based on the sample sizes.

The KS test is deliberately omnidirectional. If the two height distributions are
indistinguishable, the discordant gene trees coalesce on the same timescale as
the concordant ones, which is what introgression between the *sampled* taxa
looks like — hence `inflow_introgression` when this gate is not significant. An
empty sample yields `(0.0, 1.0)`.

Because the KS statistic responds to differences in shape, spread, and tails as
well as location, a significant result does not by itself say which direction
the heights moved. That is gate 3's job.

### Gate 3 — Adaptive studentized permutation test

**Question.** Is the *mean* concordant height greater than, less than, or
indistinguishable from the mean discordant1 height?

The answer carries a p-value and a confidence statement, so a direction is
reported only when the separation is larger than sampling noise accounts for.
The test lives in `ghostparser/orchestrator/permutation.py`.

**The statistic.** For concordant sample `x` (size `nx`) and discordant1 sample
`y` (size `ny`), the Welch-studentized mean difference is

```
T = (mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)
```

with `var` the unbiased sample variance (`ddof=1`). The denominator is the Welch
standard error, and using it rather than a pooled one is essential here: the
concordant sample is normally much larger than the discordant1 sample and the
two have different variances. Under that combination a permutation test of the
raw mean difference does *not* hold its nominal level, while the studentized
version remains asymptotically valid (Janssen 1997, *Statistics & Probability
Letters* 36(1), 9–21, https://doi.org/10.1016/S0167-7152(97)00043-6). This is
the permutation analogue of the Behrens–Fisher problem.

**The null.** Pool all `nx + ny` heights and randomly reassign them to two
groups of the original sizes. Recompute `T` — including recomputing both
variances from the permuted groups, which is what preserves the studentization —
and repeat. The resulting distribution is the null distribution of `T` under the
hypothesis that group membership carries no information.

**p-values.** Three counts accumulate over the resamples: how many permuted
statistics are `>= T_obs`, how many are `<= T_obs`, and how many exceed
`|T_obs|` in absolute value. Each becomes a p-value with the add-one estimator

```
p = (1 + count) / (1 + resamples)
```

which counts the observed arrangement itself. The naive `count / resamples`
ratio can report exactly zero and understates the true type-I error rate;
the add-one form is the correctly-sized estimator for a Monte Carlo permutation
p-value (Phipson & Smyth 2010, *Statistical Applications in Genetics and
Molecular Biology* 9(1), Article 39, https://doi.org/10.2202/1544-6115.1585).

**Correction inside the test.** The two one-tailed p-values form a testing
family of size two and are corrected against each other with the configured
`p_value_correction` method before being compared to `alpha_perm`. With the
default `bfn` this compares `2p` to `alpha_perm`, which is the conventional
relationship between a two-sided level and its two one-sided halves. This
correction is separate from and additional to the run-wide correction applied
across triplets, which covers only the DCT and KS p-values.

**Adaptive stopping.** A Monte Carlo p-value is an estimate, so the run keeps
resampling until the *decision* is safe rather than until a fixed budget is
spent. The first batch draws `min_resamples`. After each batch, a binomial
confidence interval at 95% is placed around each one-tailed p-value and
rescaled onto the corrected scale. If `alpha_perm` lies outside both intervals,
no further resampling can flip the comparison and the run stops with
`perm_converged = True`. Otherwise the batch size grows by 25% and the run
continues until the total reaches `max_resamples`, after which it stops with
`perm_converged = False` and a `max_resamples_reached` note, listed in
`metrics.txt`.

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

Because both corrections are monotone in each p-value, the interval is mapped
onto the corrected scale by the same factor the point estimate received before
being compared against `alpha_perm`.

**The sampling optimization.** A naive implementation shuffles the pooled array
once per permutation in Python, which is far too slow to run inside every
bootstrap iteration. `permutation._permutation_statistics` instead draws a whole
batch at once with three changes:

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
   return from the values themselves. The larger group's values are therefore
   never gathered at all. Since the discordant1 sample is usually the smaller
   one, the gathered data shrinks by roughly the size ratio.
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

The exactness of step 1 is not taken on trust:
`test_permutation_statistics_match_exhaustive_enumeration` enumerates all
`C(9, 4) = 126` group assignments of a small case, evaluates each through the
plain scalar statistic, and requires the vectorized sampler to emit exactly that
set of values and nothing else. Measured directly on a 40-element pooled sample,
the reconstructed complement variance differs from `numpy.var(complement,
ddof=1)` by about `7e-18` — floating-point rounding, not method error — and
centering roughly halves even that.

Batches are chunked so the matrix of random keys stays near 16 MB, which matters
because every pool worker runs its own tests concurrently.

**Guards.** Four conditions short-circuit the test to `ambiguous` before any
resampling, each recorded in the `perm_note` column:

| Note | Condition | Why |
| --- | --- | --- |
| `insufficient_group_size` | Either group has fewer than 2 observations | No unbiased variance exists, so the statistic is undefined. |
| `zero_pooled_variance` | All pooled values are effectively identical | Nothing to detect and no scale to measure it on. |
| `degenerate_observed_scale` | Both groups internally constant, means differ | The statistic divides a real difference by numerical noise and reports the p-value floor regardless of how little data backs it. |
| `insufficient_permutation_support` | `C(n, k) < min_resamples` | The permutation distribution has fewer distinct values than the requested batch, so its resolution is capped well short of `alpha_perm`. |

The scale guards compare against a small fraction of the data's own magnitude
rather than against exact zero, because a sample of nominally identical values
such as `[0.9] * 10` has a floating-point variance around `1e-33`, not `0`.

**The skew flag.** The direction comes from the pair of corrected one-tailed
tests. The absolute-value two-tailed p-value is reported as `perm_p_two_sided`
for diagnosis only and never drives a decision, because it is valid only when
the permutation null is symmetric.
When it fails to resolve a comparison that the directional tails do resolve,
`perm_null_skewed` is set. That happens when a small discordant1 group carries a
few extreme heights: reassigning them produces a strongly right-skewed null
whose fat right tail swamps a left-tail observation. Tree heights are bounded
below by zero and routinely right-skewed, so this flag is common rather than
exceptional, and it never overturns the directional call — `metrics.txt` reports
only a count.

**Known limitation.** The studentized permutation test is *asymptotically*
valid, not exact. Its type-I error rate sits close to `alpha_perm` across most
sample shapes, but it inflates when the smaller group falls below roughly 30
observations *and* carries the larger spread. At `n_dis1 = 20` against
`n_con = 200` with a 3× spread ratio the empirical rate reaches about 10% at a
nominal 5%. Treat directional calls on triplets with very few discordant1 gene
trees as provisional; the `n_dis1` column is in the results TSV for exactly this
reason.

**Inside the bootstrap.** Each bootstrap iteration re-runs the whole decision, so
the direction test runs there at one fifth of the configured `min_resamples` and
`max_resamples` (`permutation.bootstrap_resample_budget`), and only for
iterations that reach gate 3 at all. The bootstrap aggregates many iterations
into a single support value, which absorbs the extra per-iteration Monte Carlo
noise the reduced budget introduces. The point estimate runs the full budget and
runs for every triplet, including ones the earlier gates already settled, so the
permutation columns are populated throughout the results TSV.

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
an `outflow` or `ghost` call always rests on a corrected one-tailed p-value
below `alpha_perm`. A separation the data cannot resolve at that threshold is
reported `ambiguous`.

## Outputs

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
  outnumbers `BC|A` has its `AC|B` gene trees under `discordant1_*`.
- `processed_<species tree>` / `processed_<gene trees>` — cleaned, rooted trees.
- `metrics.txt` — per-stage wall/CPU timing and run parameters.
- `consolidation/` — the combined heatmap/bar-chart plot and TSV matrices from
  `introgression_mapper`, unless `--no-consolidation` is given. Consolidation
  writes into this dedicated subfolder so its own output-directory reset never
  removes the run folder's results TSV, processed trees, or the open
  `metrics.txt`.

### How each results column is produced

| Column | Source | Method |
| --- | --- | --- |
| `triplet` | Triplet setup | The normalized `(A, B, C)` labels, A and B being the species-tree sisters. |
| `species_tree` | Triplet setup | The triplet's species subtree, serialized topology-only (branch lengths omitted). |
| `dis1_topology` | Topology ranking | `BC` or `AC` — whichever discordant topology is more frequent; ties resolve to the first listed. |
| `most_frequent_matches_concordant` | Topology counts | True when the concordant count is at least both discordant counts. |
| `n_con` / `n_dis1` / `n_dis2` | Topology counts | Gene trees observed with each topology. |
| `analyzed_trees` | Extraction | Gene trees from which a subtree for this triplet was extracted. |
| `dct_statistic` / `dct_p_value` | DCT | SciPy chi-square or statsmodels z-test over `[n_dis1, n_dis2]`. An all-zero discordant split short-circuits to `(0.0, 1.0)`. |
| `dct_p_value_<method>_corr` | Correction | Run-wide correction over every triplet's DCT p-value. |
| `dct_significant` | Decision gate 1 | Corrected DCT p-value below `alpha_dct`. |
| `ks_statistic` / `ks_p_value` | Tree-height test | Two-sample KS between concordant and discordant1 heights; an empty sample yields `(0.0, 1.0)`. |
| `ks_p_value_<method>_corr` | Correction | Run-wide correction over every triplet's KS p-value. |
| `ks_significant` | Decision gate 2 | Corrected KS p-value below `alpha_ks`. |
| `perm_statistic` | Direction test | Observed Welch-studentized mean difference. Empty when a guard fired. |
| `perm_p_greater` / `perm_p_less` | Direction test | Raw one-tailed p-values, add-one estimator. |
| `perm_p_two_sided` | Direction test | Raw two-tailed p-value from the absolute statistic; diagnostic only, see `perm_null_skewed`. |
| `perm_p_greater_<method>_corr` / `perm_p_less_<method>_corr` | Direction test | The one-tailed p-values corrected against each other, and the values compared to `alpha_perm`. |
| `perm_n_resamples` | Direction test | Permutations drawn; `0` when a guard fired. Can exceed `max_resamples` by up to one batch, since the final batch is not trimmed. |
| `perm_converged` | Direction test | True when the confidence interval excluded `alpha_perm` before the budget ran out. |
| `perm_null_skewed` | Direction test | True when the absolute-value two-tailed p-value fails to resolve a comparison the directional tails resolve. Common and harmless; see above. |
| `perm_note` | Direction test | Guard slug, or `max_resamples_reached`; empty on a clean run. |
| `perm_decision` | Decision gate 3 | `greater`, `less`, or `ambiguous` for concordant relative to discordant1. Always populated; consulted only when `decision_gate` is `Permutation`. |
| `decision_gate` | Decision logic | Which test settled the classification: `DCT`, `THT`, or `Permutation`. |
| `classification` | Decision logic | `no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, or `ambiguous`. |
| `inference_description` | Reporting | Human-readable direction naming the actual species. |
| `bootstrap_value` / `all_bootstrap` | Bootstrap | Fraction of iterations agreeing with the final classification, plus the full class-fraction map. Present unless `--no-bootstrap`. |
| `bootstrap_*` debug columns | Bootstrap debug | Per-iteration DCT/KS statistics, con/dis means, and gene-tree heights. Present only with `bootstrap_debug_mode`. |

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

The results TSV is named `orchestrator_triplet_results.tsv`, distinct from the
inputs a standalone `introgression_mapper` run consumes, so both can share an
output folder without colliding.

## Parallelization modes

- `taxon` — triplets are split into chunks and dispatched across workers; each
  worker runs the fused extract-then-infer loop for its chunk over the shared
  gene-tree list.
- `gene` — triplets are processed serially in the parent; within a single
  triplet, per-gene-tree subtree extraction is parallelized across cores.
- `auto` — selects `gene` when the ingroup taxa count is below
  `AUTO_TAXA_SMALL_THRESHOLD` (15) or the gene-tree count exceeds
  `AUTO_GENE_TREES_THRESHOLD` (3500), otherwise `taxon`.

With one worker (or one triplet) the engine runs the fused loop serially in the
parent process regardless of mode.

## Internal design

### File layout

```
ghostparser/orchestrator/
  __init__.py    exports run_orchestrator
  __main__.py    python -m ghostparser.orchestrator entry point: main() wires parsing -> run_orchestrator
  config.py      orchestrator defaults/choices, validation, CLI parser, and CLI/config resolution
  correction.py  multiple-testing correction shared by inference.py and permutation.py
  permutation.py adaptive studentized permutation test (decision gate 3)
  trees.py       tree/triplet preprocessing
  inference.py   per-triplet inference + summary stats + result type + TSV writers
  stream.py      fused extract+infer streaming engine
  preflight.py   structural data check reached via --preflight-data-check
  runner.py      run_orchestrator coordinator (cleaning -> streaming -> correction -> writing -> consolidation)
  ORCHESTRATOR.md    this document
```

### Shared dependencies

The orchestrator owns its tree preprocessing, inference, and configuration. It
imports only four things from the rest of the package:

- `ghostparser.config` — the shared configuration trunk (`ConfigError`, path
  resolution, raw config-file loading, required-path validation, overwrite
  resolution, `prepare_output_directory`). Orchestrator-specific defaults, choices,
  and validators live in `orchestrator/config.py`, which is why the orchestrator can set
  its own `bfn`/`mean` defaults without affecting the ML subpackage.
- `ghostparser.cli_config` — the generic `resolve_cli_or_config_args` resolver
  implementing config-file-wins precedence.
- `ghostparser.triplet_utils` — pure topology helpers.
- `ghostparser.introgression_mapper` — `generate_introgression_maps` for the
  consolidation stage.

### Fused streaming engine

The per-triplet unit is `inference.analyze_triplet_from_observations(triplet,
observations, species_subtree, ...)`, which takes precomputed `(topology,
tree-height, metrics)` observations. The third element carries the per-tree
summary metrics and is `None` unless `generate_summary_stats` is enabled.

The processing unit is a chunk of triplets that share one parse pass over the
gene trees. For each parsed gene tree the engine extracts every in-chunk
triplet's subtree and computes its observation directly from the subtree object
(`inference.observation_from_subtree`) — there is **no** serialize-to-Newick and
reparse round trip. It then runs `analyze_triplet_from_observations` for each
triplet in the chunk and drops the observations. This bounds live memory to one
chunk while amortizing the DendroPy parse cost across the chunk's triplets. Gene
trees are loaded once in the parent and shared read-only to workers via a
fork/forkserver initializer.

Bootstrap resampling is vectorized with NumPy: per-triplet resample indices are
drawn with a seeded `numpy.random.Generator`, and topology counts and
per-topology height groups are computed with array operations. Bootstrap values
are deterministic under a fixed `bootstrap_seed` — the per-triplet seed is
derived from the run seed and the triplet, so every parallelization mode agrees
exactly.

### Memory rationale

- No intermediate triplet-gene-trees file is written or reloaded.
- Peak memory is roughly the shared gene-tree Newick list, plus one chunk of
  transient subtrees, plus the accumulating list of small result objects.
- Result objects must be accumulated because global p-value correction needs all
  p-values in a single pass.

### Consolidation interface contract

`generate_introgression_maps` reads only four fields from each result
(`triplet`, `classification`, `dis1_topology`, `bootstrap_value`) via its
duck-typed `_extract_result_fields`. `TripletPipelineResult` keeps those four
fields; treat field parity on them as a maintenance constraint.
