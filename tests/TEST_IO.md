# Test Input/Output Derivations

This document is the granular companion to [TESTS.md](TESTS.md). For each test
it records **exactly what the inputs are** (literal values, not just fixture
names) and **how the expected output is derived**: the arithmetic, the
definition, or the reference call that produces it.

The suite's guiding rule is that expectations are computed from first
principles, never copied from a previous run. This file is where that
derivation is written down.

## Shared Inputs

### Orchestrator fixture: species tree

`tests/orchestrator/conftest.py` writes one species tree:

```
(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);
```

Rooting on `OUT` and pruning it dissolves the `(D, OUT)` node. `D` keeps its own
`0.1` edge, and the `((A,B),C)` clade absorbs both its own `0.1` and the
dissolved node's `0.2`, giving `0.3`. `OUT`'s `0.5` becomes the root edge:

```
(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;
```

Ingroup taxa are therefore `A, B, C, D`, and the species topology is
`(((A,B),C),D)`.

### Orchestrator fixture: gene trees

The same conftest writes 12 gene trees, all of the shape
`((((X:l,Y:l):l,Z:l):l,D:l):l,OUT:l);`. Rooting on `OUT` folds `OUT`'s edge into
the ingroup clade's edge, and pruning `OUT` leaves that clade as the root.
For tree 0 the ingroup edge `0.10` plus `OUT`'s `0.50` gives `0.6`:

```
(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6;
```

Reading the innermost sister pair off each tree gives, in input order:

| Tree | Innermost pair | Tree | Innermost pair |
| --- | --- | --- | --- |
| 0 | (A,B) | 6 | (A,B) |
| 1 | (A,B) | 7 | (A,C) |
| 2 | (A,C) | 8 | (A,B) |
| 3 | (B,C) | 9 | (B,C) |
| 4 | (A,B) | 10 | (A,B) |
| 5 | (B,C) | 11 | (A,B) |

So for triplet `(A,B,C)`: `(A,B)` appears in trees 0, 1, 4, 6, 8, 10, 11 = **7**;
`(B,C)` in trees 3, 5, 9 = **3**; `(A,C)` in trees 2, 7 = **2**.

For any triplet containing `D`, the two non-`D` taxa always sit together inside
the `((X,Y),Z)` clade with `D` outside, so **every** gene tree yields the
concordant topology: counts are **12 / 0 / 0**.

### Inference fixture: 10 gene subtrees

`tests/orchestrator/test_orchestrator_inference.py` uses 10 three-taxon subtrees with
species subtree `((A,B),C);`, the topology-only Newick a run hands each
triplet (concordant topology `((A,B),C)`).
Each subtree's geometry is tabulated in `_LEAF_GEOMETRY` as
`(topology, dist_A, dist_B, dist_C, internal_branch)`, where a leaf's
root-to-tip distance is its own edge plus the internal branch if it is in the
sister pair, or just its own edge otherwise.

Worked example, index 0, `((A:0.10,B:0.10):0.10,C:0.30);`:
`A = 0.10 + 0.10 = 0.20`, `B = 0.10 + 0.10 = 0.20`, `C = 0.30`, internal `0.10`.

| Idx | Topology | A | B | C | internal |
| --- | --- | --- | --- | --- | --- |
| 0 | ((A,B),C) | 0.20 | 0.20 | 0.30 | 0.10 |
| 1 | ((A,B),C) | 0.21 | 0.20 | 0.32 | 0.09 |
| 2 | ((A,B),C) | 0.30 | 0.25 | 0.40 | 0.10 |
| 3 | ((B,C),A) | 0.50 | 0.30 | 0.30 | 0.10 |
| 4 | ((B,C),A) | 0.55 | 0.35 | 0.32 | 0.10 |
| 5 | ((A,C),B) | 0.40 | 0.70 | 0.40 | 0.10 |
| 6 | ((A,B),C) | 0.27 | 0.27 | 0.35 | 0.12 |
| 7 | ((B,C),A) | 0.65 | 0.45 | 0.45 | 0.10 |
| 8 | ((A,C),B) | 0.28 | 0.60 | 0.26 | 0.10 |
| 9 | ((A,B),C) | 0.21 | 0.23 | 0.31 | 0.10 |

Counts: concordant `((A,B),C)` = **5** (0, 1, 2, 6, 9), `((B,C),A)` = **3**
(3, 4, 7), `((A,C),B)` = **2** (5, 8). Ranking the discordants puts `BC` first
(3 >= 2), so `dis1_topology == "BC"`.

### Tree-height strategies

Each strategy is a function of the row above:

- `AVG` = `(dist_A + dist_B + dist_C) / 3`
- `A` / `B` / `C` = that taxon's root-to-tip distance
- `SIS` = `dist(left sister) + dist(right sister) - 2 x internal`
- `INT` = `internal`

Worked example, index 5, `((A:0.30,C:0.30):0.10,B:0.70);` (sisters A, C):
`AVG = (0.40 + 0.70 + 0.40)/3 = 0.50`; `A = 0.40`; `B = 0.70`; `C = 0.40`;
`SIS = 0.40 + 0.40 - 2(0.10) = 0.60`; `INT = 0.10`.

### Summary statistics

- `mean`: arithmetic mean.
- `median`: middle value (mean of the two middle values for an even count).
- `mode`: used only by the 63-column `summary_statistics.tsv` output, not by
  the decision logic. Values are **binned to 3 decimals**, the most frequent bin
  wins, and ties resolve to the **largest** value. With all-distinct values every
  bin has count 1, so the mode is the maximum. For the AVG concordant sample
  `[0.233, 0.243, 0.317, 0.297, 0.250]` the mode is therefore `0.317`.

## tests/orchestrator/test_orchestrator_inference.py

### `test_inference_matches_derived_expectation`

**Inputs:** the 10 gene subtrees above serialized with the `AVG` strategy,
`alpha_dct = alpha_ks = 0.05`, one of the 2 discordant tests, and `diagnostic`
off or on. Bootstrap runs with 40 iterations at seed `20240724`. The
result is measured by `analyze_triplet_from_observations` and decided by
`_apply_triplet_result_p_value_correction` under `no` as a family of one.

**Expected-output derivation**, performed in `_expected_result`:

1. Group the hand-derived heights by topology using `_LEAF_GEOMETRY`.
2. Rank the discordants by count: `n_dis1 = 3` (BC), `n_dis2 = 2` (AC),
   `n_con = 5`.
3. `most_frequent_matches_concordant` = `5 >= 3 and 5 >= 2` = `True`.
4. **DCT**, chi-square: `scipy.stats.chisquare([3, 2])` with implied expected
   frequencies `[2.5, 2.5]`, giving
   `(3-2.5)^2/2.5 + (2-2.5)^2/2.5 = 0.1 + 0.1 = 0.2`, and with df 1,
   `p ~ 0.6547`. z-test: `proportions_ztest(count=[3,2], nobs=[5,5])`, giving
   `z ~ 0.6325`, `p ~ 0.5271`. Both p-values exceed 0.05, so
   `dct_significant is False`; gate 1 settles the call.
5. **KS**, with `diagnostic` on, measured regardless:
   `scipy.stats.ks_2samp(concordant_heights, dis1_heights)` on
   `[0.2333, 0.2433, 0.3167, 0.2967, 0.2500]` and `[0.3667, 0.4067, 0.5167]`;
   they are completely separated, so `D = 1.0` and `p ~ 0.0357`, giving
   `ks_significant is True`. The cascade never reads it, which is what the row
   demonstrates. With `diagnostic` off the correction is `no`, an inline
   method whose family of one is fixed, so the test below the failed count
   gate is never measured: `ks_statistic`, `ks_p_value` and `ks_significant`
   are all `None`.
6. **Direction**, with `diagnostic` on, measured regardless: the pooled sample is
   `n_con + n_dis1 = 5 + 3 = 8` observations, admitting only
   `C(8, 3) = 56` distinct group assignments. That is far below the 2500
   `min_resamples` floor, so the `insufficient_permutation_support` guard fires
   and `perm_decision` is `inconclusive` with zero resamples. With
   `diagnostic` off the failed count gate settles the cascade first, so the
   test is skipped and `perm_decision` is `None`.
7. **Classification**, gate 1 fails, so `no_introgression` in every case.

Also asserted: `triplet == ("A","B","C")`, `species_tree == "((A,B),C);"` (the
topology-only species subtree, stored as given), and that the bootstrap class
fractions sum to 1.

### `test_observation_heights_match_derived_geometry`

**Inputs:** the same 10 subtrees, serialized under each of the six strategies
in turn.

**Derivation:** for each observation, the topology must equal the tabulated
topology and H(T) must equal `_expected_height(entry, strategy)` computed from
the four geometry primitives, i.e. the strategy formulas above applied to the
table, independent of the implementation.

### `test_empty_observations_decide_no_introgression`

**Inputs:** an empty observation list, measured and decided as above.

**Derivation:** with no observations there are no topology counts, so
`analyzed_trees = 0` and `n_con = 0`; the DCT short-circuits on a zero total to
`(0.0, 1.0)`, which is not significant, so gate 1 returns `no_introgression`.

## tests/orchestrator/test_orchestrator_decision.py

Observations are constructed directly as `(topology, height, None)` tuples, which
lets each test place the triplet on a chosen branch. Bootstrap is disabled
(`iterations: 0`) so results are deterministic.

### `test_decision_cascade_lands_on_each_classification`

**Inputs:** one crafted observation set per row, measured with no bootstrap
with `diagnostic` off and on and decided under `no` correction as a family of one,
so every corrected p-value equals its raw one.

| id | con | dis1 | dis2 |
| --- | --- | --- | --- |
| `no_introgression` | 20 at 0.1 | 10 at 0.9 | 10 at 0.9 |
| `inflow` | 10 at 0.5 | 30 at 0.5 | 2 at 0.5 |
| `outflow` | `_HIGH` = `[0.85, 0.86, ..., 0.94]` | `_LOW` = `[0.05, 0.06, ..., 0.34]` | 2 at 0.1 |
| `ghost` | `_LOW[:10]` = `[0.05, ..., 0.14]` | `_HIGH * 3` | 2 at 0.9 |
| `ambiguous` | 30 alternating `0.5 +- 0.30` | 30 alternating `0.5 +- 0.02` | 2 at 0.5 |

**Derivation:**

- **`no_introgression`.** `chisquare([10, 10])` has expected `[10, 10]`, so the
  statistic is exactly `0.0` and `p = 1.0 > 0.05`: gate 1 fails. The heights are
  fully separated (20 at 0.1 against 10 at 0.9, `D = 1.0`), so KS *is*
  significant and with `diagnostic` on the row asserts `ks_significant is
  True`, which is why this row also shows the DCT gate stopping the cascade
  before a later gate can be consulted; the direction test runs as well and
  reports a decision the classification never reads. With `diagnostic` off
  neither is measured: `no` is an inline correction, so the failed count gate is final
  and the row asserts `ks_p_value`, `ks_significant` and `perm_decision` all
  `None` with `perm_note == "direction_test_not_consulted"`.
- **`inflow`.** `chisquare([30, 2])` has expected `[16, 16]`, so the statistic is
  `(30-16)^2/16 + (2-16)^2/16 = 24.5`, and with df 1 `p ~ 7.4e-07 < 0.05`: gate 1
  passes. Every concordant and dis1 height is 0.5, so the two empirical CDFs
  coincide: `D = 0.0`, `p = 1.0`, not significant → gate 2 returns
  `inflow_introgression`.
- **`outflow`.** The DCT is the same significant 30-vs-2 split.
  `max(_LOW) = 0.34 < min(_HIGH) = 0.85`, so the CDFs separate completely:
  `D = 1.0` and gate 2 passes. The spread matters for gate 3: constant samples
  would trip the `degenerate_observed_scale` guard, but these carry real
  within-group variance, and the pooled 40 observations give
  `C(40, 10) = 847,660,528` assignments, far above the 2500 floor. No permutation
  reproduces the observed statistic, so `p_greater` sits at the add-one floor
  `1/(2500+1) = 4.0e-4`; Bonferroni over the tail pair doubles it to
  `8.0e-4 < 0.05` → `greater` → `outflow_introgression`.
- **`ghost`.** The mirror image: identical DCT and KS reasoning, with gate 3
  finding a negative studentized difference and `p_less` at the floor → `less` →
  `ghost_introgression`.
- **`ambiguous`.** The DCT sees a 30-vs-2 split → significant. The two samples
  share the mean 0.5 but their CDFs differ sharply in spread, so KS is large and
  gate 2 passes. Gate 3 finds no difference in means: the observed statistic is
  near 0 and both one-tailed p-values are far above `alpha_perm`, so the outcome
  is non-directional and the classification is `ambiguous`. The distributions
  differ in shape, not location.

The `decision_gate` assertion is what makes each row specific: a case that
reached the same classification by a different route would fail. With
`diagnostic` off, `perm_decision` is asserted non-`None` only on the three
`PERM` rows and `None` on the `inflow` row too, whose failed tree-height gate
settles the cascade before the direction test.

### `test_permutation_guards_surface_on_the_triplet_result`

**Inputs and derivation**, one row per guard, measured with
`diagnostic=True`: the first row's 2-vs-2 discordant split gives
`chisquare([2, 2])` a statistic of `0.0` and `p = 1.0`, so its count gate
fails and a non-diagnostic run would skip the direction test before the guard
could fire:

| con | dis1 | Guard | Why |
| --- | --- | --- | --- |
| `[0.9, 0.8, 0.7, 0.6]` | `[0.1, 0.2]` | `insufficient_permutation_support` | `C(6, 2) = 15 < 2500`, so the permutation distribution cannot resolve `alpha_perm`. |
| `[0.9] * 10` | `[0.1] * 30` | `degenerate_observed_scale` | Both groups are internally constant, so the standard error is floating-point noise near `1e-17`. |

Each returns `perm_n_resamples == 0` and `perm_decision == "inconclusive"`: a
guard means nothing was established, which is distinct from having shown the
means equivalent.

### `test_summary_statistics_discordant_roles_follow_the_counts`

**Inputs:** 10 concordant subtrees at height 0.30, plus 3 and 9 discordant
subtrees. Row `AC_more_frequent` puts `[0.10, 0.12, 0.14]` on `BC|A` and
`[0.50, 0.51, ..., 0.58]` on `AC|B`; row `BC_more_frequent` swaps them.

**Derivation:** each subtree is `((X:h,Y:h):0.10, Z:h+0.2)`, so the two sisters
sit at `h + 0.10` and the outlier at `h + 0.20`, giving an average tree height of
`(3h + 0.4) / 3 = h + 0.4/3`. Every group mean is therefore its input mean
shifted by exactly `0.4/3`.

`dis1_topology` names whichever discordant topology is more frequent (`AC` on
the first row, `BC` on the second) and `discordant1_*` must describe that same
group of gene trees, with `discordant2_*` describing the other. Counts are
`(n_dis1, n_dis2) == (9, 3)` either way, so the roles follow the counts rather
than the topology label.

### `test_classify_introgression_truth_table`

**Inputs:** `_classify_introgression(dct_significant, ks_significant,
direction)` called directly with 9 explicit rows; it returns the
`(classification, decision_gate)` pair asserted below. The same row goes to
the array form `_classification_codes` as one-element arrays, with a `None`
tree-height flag passed as `False` and the direction mapped through
`_DIRECTION_CODES` (anything but `greater`/`less` is the skipped code).

**Derivation:** straight from the decision definition. The gate is the name of
the test whose branch returned, so it is fixed by `dct_sig` and `ks_sig` alone:
`direction` only ever selects among the three `PERM` classifications. The
array form writes integer codes indexing `_BOOTSTRAP_CLASSES`, so
`_BOOTSTRAP_CLASSES[code]` must equal the scalar label on every row; a `None`
tree-height flag takes the not-significant branch in the scalar form, which
is why `False` stands in for it.

| dct_sig | ks_sig | direction | Expected | Gate | Reason |
| --- | --- | --- | --- | --- | --- |
| False | True | `greater` | `no_introgression` | `DCT` | DCT fails first |
| False | False | `less` | `no_introgression` | `DCT` | DCT fails first |
| True | False | `greater` | `inflow_introgression` | `THT` | tree-height test not significant |
| True | `None` | `greater` | `inflow_introgression` | `THT` | no THT ran; `not None` takes the same branch |
| True | True | `greater` | `outflow_introgression` | `PERM` | con > dis |
| True | True | `less` | `ghost_introgression` | `PERM` | con < dis |
| True | True | `equivalent` | `ambiguous` | `PERM` | means shown close, no direction |
| True | True | `inconclusive` | `ambiguous` | `PERM` | no direction resolved |
| True | True | `None` | `ambiguous` | `PERM` | no direction available |

### `test_adjust_p_values_matches_statsmodels_and_never_lowers_a_value`

**Inputs:** `[0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6]`
and then the family `[0.001] * 8 + [0.4, 0.9]`, under each method in
`P_VALUE_CORRECTION_CHOICES`, `alpha = 0.05`.

**Derivation:** on the first list `no` must return the input unchanged and
every other method must equal
`statsmodels.stats.multitest.multipletests(p_values, alpha=0.05, method=m)[1]`
where `m` maps `bfn → bonferroni` and the rest to their own names; Bonferroni
is `min(1, n x p)`, so the reference also pins that arithmetic.

On the second list each of `no`, `bfn`, `holm`, `fdr_bh`, and `fdr_by`
applies a multiplier of at least 1 to the `j`-th smallest of `n` (`n` for
Bonferroni, `n - j + 1` for Holm, `n/j` for BH, and `(n/j) x sum(1/i)` for BY
), so no adjusted value can fall below its raw one. The family is deliberately
the hostile case for that property: a method that estimates the number of
true nulls `n0` and substitutes it for `n` would reject 8 at the first stage,
giving `n0 = 2` and a multiplier of `2/j`, below 1 for every `j > 2`, so the
later strong p-values would land beneath their raw ones. Parametrizing over
the choice list rather than a fixed set of names means such a method fails
here the moment it is added. The `1e-12` slack absorbs floating-point
rounding in the running max/min sweeps.

### `test_inline_correction_matches_the_family_pass_or_refuses`

**Inputs:** `p = 0.004` under each method in `P_VALUE_CORRECTION_CHOICES`;
for the inline methods, padded with `0.5` entries to families of 1, 7 and 250.

**Derivation:** `no` and `bfn` depend on the family only through its size, so
the inline form must equal the full pass exactly: `0.004` at every size under
`no`, and `min(1, n x p)` under `bfn`: `0.004`, `0.028`, and `1.0`
(`250 x 0.004 = 1.0`). Holm's, BH's and BY's multipliers depend on a p-value's
rank within its family, which a single value does not determine, so
`_adjust_p_value_inline` must raise for them rather than silently pick a wrong
multiplier.

### `test_unknown_methods_are_rejected`

**Inputs:** `_adjust_p_values([0.1], method="bogus")` and
`run_discordant_count_test(3, 2, method="bogus")`.

**Derivation:** an unsupported name is outside the choice tuple, so a
`ValueError` naming the valid options is raised before any computation.

### `test_degenerate_samples_are_non_significant`

**Inputs:** `run_discordant_count_test(0, 0)` under `chi-square` and `z-test`;
`run_two_sample_ks_test([], [1.0])` and `run_two_sample_ks_test([], [])`.

**Derivation:** `n_dis1 + n_dis2 == 0` means there is nothing to test, so the
count test returns the neutral `(0.0, 1.0)` before touching SciPy or
statsmodels; an empty sample on either side makes the KS statistic undefined,
so that function returns the same neutral pair.

### `test_inline_and_deferred_correction_agree_on_a_single_triplet`

**Inputs:** one triplet with 25 concordant heights at 0.9, 22 discordant1 at 0.35
plus 3 at 0.2, and 4 discordant2 at 0.3; `family_size=1`, 40 bootstrap
iterations, `triplet_seed=3`; each of `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`,
measured first and decided afterwards.

**Derivation:** with one test in the family every correction is the identity:
Bonferroni multiplies by 1, and the rank-based methods adjust the single value
`p_(1)` by `(n - 1 + 1)/1 = 1` (Holm), `n/1 = 1` (BH), or `1 x sum(1/i) = 1` for
`n = 1` (BY). So all five must reach the same per-iteration verdicts. `no` and
`bfn` reach them inline during the stream, so the measured result already
carries `all_bootstrap` and no deferred record; the other three cannot rank a
p-value without the rest of its family, so the measured result has
`all_bootstrap is None` and a `bootstrap_deferred` record with one entry per
iteration, 40, whatever the family size.
`_apply_triplet_result_p_value_correction` then corrects iteration `i` across
the (here single-member) family, tallies the votes and clears the record, so
after the pass `bootstrap_deferred is None` for every method and the tally
equals the `no` run's. Equality is a statement about the two code paths rather
than about arithmetic. The resample stream is seeded identically and is
independent of the permutation generator, so the iterations themselves are
the same draws in every case.

### `test_bootstrap_votes_answer_to_the_corrected_threshold`

**Inputs:** 40 concordant heights at 0.9, 18 discordant1 and 6 discordant2 at
0.55; run once with `no` at family size 1 and once with `bfn` at family size
5000, 40 iterations, `triplet_seed=3`.

**Derivation:** the discordant split 18 vs 6 gives a chi-square of
`(18-12)^2/12 + (6-12)^2/12 = 6.0` on 1 df, `p = 0.0143 <= 0.05`, so the raw DCT
gate passes. Bonferroni over 5000 tests gives `0.0143 x 5000 = 71.5 → 1.0`, far
above `alpha_dct`, so the corrected gate fails and the point estimate is
`no_introgression`. Every bootstrap iteration is corrected by the same factor;
no resample of a p-value near 0.014 survives a 5000x multiplier, so all 40
iterations vote `no_introgression` and the fraction is exactly 1.0. Under `no`
the same resamples are judged raw, and enough of them clear 0.05 that the
fraction falls below 1.0, which is what makes the corrected agreement a real
check rather than a tautology.

### `test_studentized_interval_brackets_the_observed_statistic`

**Inputs:** concordant heights `0.50, 0.51, ..., 1.09` (60 values) against
discordant1 `0.20, 0.21, ..., 0.59` (40 values), 5 discordant2 at 0.3, 200
iterations, `no` correction.

**Derivation:** the interval is the empirical `[100 x alpha, 100 x (1 - alpha)]`
percentile pair over the per-iteration studentized differences, i.e. the 5th and
95th percentiles at `alpha_perm = 0.05`. Each iteration resamples the observed
gene trees with replacement, so the bootstrap distribution is centred on the
statistic computed from the observed data; a 90% percentile range of a
distribution centred on that value contains it. Ordering is immediate from the
percentile definition.

Two concordant heights at 0.9 against a single discordant1 height at 0.2 give the
degenerate half of the same test: every observation within a group carries the
same height, so both group variances are zero in any resample that manages two
draws from each group, and resamples that do not have fewer than 2 observations
somewhere. Either way the statistic is `nan`, the percentile is taken over an
empty set, and both bounds are reported `None` rather than fabricated.

### `test_results_tsv_carries_corrected_columns_only_when_correcting`

**Inputs:** one triplet of 25 concordant heights at 0.9, 20 discordant1 at 0.35,
and 4 discordant2 at 0.3, analyzed and corrected under `no`, `bfn`, and
`fdr_bh`, then written with `write_pipeline_results`.

**Derivation:** correction is applied to the p-value, so under `no` the
corrected value is the raw value by definition and a `dct_p_val_no_corr` column
would repeat `dct_p_value` exactly. The writer therefore emits
`dct_p_val_<method>_corr`, `ks_p_val_<method>_corr`,
`perm_p_greater_<method>_corr`, and `perm_p_less_<method>_corr` only when the
method is not `no`: 36 columns with a correction against 32 without. The
significance flags are unconditional because the cascade reads them whatever the
method is. The field-count assertion catches the failure mode this change could
introduce: the header and the row build their conditional sections separately,
so a mismatch would silently shift every later column by one.

## tests/orchestrator/test_orchestrator_permutation.py

Samples here are drawn from seeded generators rather than written out, so the
derivations below name the distribution and the property being checked instead
of a literal array.

### `_random_samples`

The shared input generator. Sizes are drawn uniform on `[8, 120)`, the location
shift uniform on `[-1, 1]`, and both scales uniform on `[0.2, 2.0]`. One of
three families is picked with equal probability: two normals, two lognormals, or
two exponentials with the second shifted. This deliberately spans symmetric,
right-skewed, and heavy-tailed data at unequal sizes and variances, which is the
regime the studentization is there to handle.

### `test_statistic_p_values_and_verdict_match_scipy`

**Inputs:** five seeded sample pairs from `_random_samples`,
`min_resamples = max_resamples = 4000` (pinning the adaptive stopping off),
`correction="bfn"`, `alpha = 0.05`.

**Derivation:** the observed statistic is a deterministic function of the
inputs, not of the resampling, so it must equal SciPy's `result.statistic` for
either alternative and the in-test reference
`(mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)` exactly (to floating-point
tolerance). Any disagreement would mean the two are not testing the same
quantity.

The raw one-tailed p-values (`p_greater`, `p_less` are reported uncorrected
whatever the correction) and SciPy's estimate the same quantity from
independent resampling streams, so exact equality is not expected. Each
estimate has binomial standard error `sqrt(p(1-p)/n)`, and the difference of
two independent estimates has `sqrt(2)` times that. The tolerance is 5 such
standard errors, which at `p = 0.5` and `n = 4000` is about `0.056` and at
`p = 0.01` about `0.011`. A systematic error in the sampler or the counting
would exceed this; ordinary Monte Carlo scatter will not.

With Bonferroni over a family of two the corrected p-value is `2p`, so
`2p <= alpha` is the same condition as `p <= alpha/2`. The expected verdict is
therefore computed from SciPy's raw one-tailed p-values at `alpha/2 = 0.025`:
`greater` if the greater-tail p-value clears it, else `less` if the other
does, else non-directional. It is asserted only when both of SciPy's p-values
sit further than the Monte Carlo tolerance from `0.025`, since inside that band
the two streams can legitimately land on opposite sides of the threshold.

### `test_permutation_statistics_match_exhaustive_enumeration`

**Inputs:** `x = [0.11, 0.24, 0.37, 0.52]` (nx=4),
`y = [0.63, 0.71, 0.88, 0.95, 1.10]` (ny=5), 5000 sampler draws at seed 7.

**Derivation:** with 9 pooled values split 4/5 there are exactly
`C(9, 4) = 126` distinct group assignments. Enumerating all of them with
`itertools.combinations` and evaluating each through the scalar reference
`_studentized_mean_diff` gives the exact support of the permutation
distribution, rounded to 9 decimals. The vectorized sampler (which never
materializes the 5-element group, recovering its sum and sum-of-squares by
subtracting from the pooled totals) must emit only values from that set
(soundness). Over 5000 draws from 126 equally likely assignments the chance of
missing any one is about `126 x (125/126)^5000 ~ 1e-15`, so all 126 must also
appear (completeness). Together these pin the optimization to ground truth.

### `test_random_inputs_preserve_test_invariants`

**Inputs:** twelve seeded pairs from `_random_samples`, `min_resamples=600`,
`max_resamples=3000`, each run twice with the same resampling seed; then the
same 8 values `[0.10, 0.22, 0.31, 0.44, 0.55, 0.61, 0.78, 0.83]` as both
samples at a fixed 600 resamples.

**Derivation:** four properties follow from the construction regardless of data.
p-values use the add-one estimator so they are strictly positive and at most 1.
`count_greater` counts `T_perm >= T_obs` and `count_less` counts
`T_perm <= T_obs`; every resample satisfies at least one and ties satisfy both,
so the counts sum to at least `n_done` and
`p_greater + p_less = (2 + count_greater + count_less)/(n_done+1) > 1`. The
resample total lies within the configured bounds by the loop's construction.
And a directional decision requires the corresponding tail to be small, which
requires the observed statistic to sit on that side of the null, so `greater`
implies a positive statistic and `less` a negative one.

Every random draw in the test comes from the passed generator, so two runs
seeded alike must agree on every field, not approximately, but exactly,
including the resample count the adaptive rule stopped at. This is the
property the whole per-triplet seeding scheme rests on: the orchestrator hands
each triplet a stream derived from `(seed, triplet)`, so a triplet's result
cannot depend on how many workers ran or in what order they finished.

Identical samples have identical means, so the numerator is exactly 0 and the
statistic is 0. Half the permutation distribution lies on either side of 0, so
both one-tailed p-values are near 0.5 and neither clears `alpha`. The
equivalence step then decides between `equivalent` and `inconclusive`; at 8
observations per group it has too little power to rule out a medium effect, so
the assertion accepts either.

### `test_equivalence_step_decides_from_the_directional_resamples`

**Inputs:** `x` and `y` both drawn from `N(1.0, 0.2^2)` at n=8 and at n=400 per
group, 1000 resamples, fixed seeds; and the n=400 pair with
`equivalence_test=False`.

**Derivation:** both samples come from one distribution, so neither directional
tail can be significant and the TOST step decides. TOST is two one-sided
tests: each null puts the mean difference at least a margin to one side of
zero, and rejecting both confines the difference to within the margin, which
is what `equivalent` asserts. The margin is `0.5` pooled
standard deviations, so the shift applied to each null is `0.5 x SD` in raw
units, and the studentized size of that shift is `0.5 x SD / SE`, which grows
like `sqrt(n)`: with equal groups `SE = SD x sqrt(2/n)`, so the shifted null is
displaced by `0.5 / sqrt(2/n) = sqrt(n/8)` null standard deviations: one at
n=8, where the shifted nulls still cover the observed statistic and TOST
reports `inconclusive`, and about 7 at n=400, where neither shifted null
reaches it and the two one-sided tests reject, giving `equivalent` with
`p_tost <= 0.05`. A margin expressed in standard-error units would displace
the null by the same amount at every n, because the studentized statistic is a
pivot whose null spread stays near 1, which is why the margin is an effect
size rather than a number of standard errors.

The equivalence p-value is an add-one estimator over the same draws the
directional test made, so it is `(1 + k)/(n_resamples + 1)` for some count
`k`: at least `1/(n_resamples + 1)` and an exact multiple of it, which the test
checks by dividing it out and comparing to the nearest integer. No further
resampling happens for it.

Switched off, there is nothing left to distinguish `equivalent` from
`inconclusive` when no direction is significant, so `p_tost` must be `None`
and the decision must fall through to `inconclusive`, never to
`equivalent`, which would be an equivalence claim no test supported.

### `test_type_one_error_rate_tracks_alpha_under_unequal_variance`

**Inputs:** 300 replicates, `x ~ N(1.0, 1.0^2)` at n=60 against
`y ~ N(1.0, 0.3^2)` at n=180, `alpha = 0.05`, 1000 resamples, `correction="bfn"`.

**Derivation:** both samples share the mean 1.0, so every *directional* decision
is a type-I error; `equivalent` and `inconclusive` are not rejections of the
directional null and are not counted. A correctly sized test rejects with probability `alpha`, giving
`300 x 0.05 = 15` expected rejections with standard deviation
`sqrt(300 x 0.05 x 0.95) = 3.8`. The assertion band `[3, 30]` spans 1% to 10%,
roughly `+-4` standard deviations, so fixed seeds make it stable while a test
that had stopped controlling its error rate would fall outside. Unequal sizes
paired with unequal variances is the configuration where a permutation test of
the *raw* mean difference loses its level, so this measures what the
studentization buys. The level degrades once the smaller group falls below
roughly 30 observations while carrying the larger spread; that limitation is
recorded in ORCHESTRATOR.md rather than asserted here.

### `test_guards_short_circuit_without_resampling`

**Inputs and derivation**, one row per guard:

| Input | Guard | Why |
| --- | --- | --- |
| `x=[0.4]`, `y=[0.1, 0.2, 0.3]` | `insufficient_group_size` | `nx = 1`, and `np.var(ddof=1)` needs at least 2 observations. |
| `x=y=[0.3, 0.3, 0.3]` | `zero_pooled_variance` | All six pooled values equal, so both the spread and the mean difference are zero. |
| `x=[0.9]*10`, `y=[0.1]*30` | `degenerate_observed_scale` | Both groups internally constant with different means; the standard error is floating-point noise near `1e-17`. |
| `x=[0.9,0.8,0.7,0.6]`, `y=[0.1,0.2]` | `insufficient_permutation_support` | `C(6, 2) = 15 < 2500`. |

Each returns `statistic is None`, `n_resamples == 0`, `converged is False`, and
decision `inconclusive`: a guard means nothing was established, which is
distinct from having shown the means to be equivalent.

### `test_null_skewness_matches_scipy_and_the_null_shape`

**Inputs:** two fixtures. *Asymmetric:* 700 concordant heights from
`N(0.42, 0.10)`, 19 discordant1 of which 15 are from `N(0.5, 0.1)` and 4 from
`N(25, 5)`; 2500 resamples, `bfn`. *Symmetric:* two 150-observation samples
from `N(1.0, 1.0)`, 4000 resamples.

**Derivation:** in both cases the expected skewness is not derived
independently: it is `scipy.stats.skew` over the statistics drawn from the
same seed. That is the point: `null_skew` is accumulated from running power
sums so batches can be discarded, and the test pins that accumulation against
the reference computed from the retained draws, with `n_resamples_skew` equal
to the draw count and no guard tripped.

For the asymmetric fixture the permutation statistic depends sharply on how
many of the four extreme heights land in the 19-slot group, so the null
separates into clusters rather than one smooth curve: drawing 19 of 719
without replacement, the chance that none of the four extremes is among them
is `C(715, 19) / C(719, 19) ~ 0.90`, that exactly one is about `0.10`, and that
two are about `0.005`; the three clusters sit at `T` near `+1.9`, `-0.9` and
`-1.4`. A third moment over that mixture is strongly negative, so
`|null_skew| > 1`; the four extremes pull the discordant1 mean far above the
concordant one, so the statistic is negative and the decision `less`.

For the symmetric fixture equal group sizes drawn from one symmetric family
give a permutation null that is symmetric about zero, so its population
skewness is 0. The `0.15` band is Monte Carlo slack: the standard error of a
sample skewness is about `sqrt(6/n) = 0.039` at n = 4000, so the bound is
roughly 4 standard errors.

### `test_adaptive_run_converges_or_exhausts_its_budget`

**Inputs:** *converges:* `x ~ N(2.0, 0.2^2)` and `y ~ N(0.5, 0.2^2)`, both
n=80, `min_resamples=1000`, `max_resamples=20000`. *Exhausts the budget:* 30
draws from `N(1.4, 1.0^2)` against 30 from `N(1.0, 1.0^2)` (sample seed 3),
`min_resamples=100`, `max_resamples=1000`, resampling seed 103.

**Derivation:** the first pair is separated by more than seven pooled standard
deviations, so no permutation approaches the observed statistic and
`p_greater` lands at the floor `1/1001 ~ 1.0e-3`. The Wilson interval around a
count of 1 in 1001 is roughly `[0.0003, 0.0056]`; doubled by Bonferroni it stays
far below `alpha = 0.05`, so `alpha` is outside it after the very first batch.
The run therefore stops with `batches == 1`, `n_resamples == 1000`, and decision
`greater`: an easy case must not spend the ceiling.

For the second pair the 0.4 shift is marginal at these sample sizes, so the
corrected p-value stays near `alpha` and the Wilson interval never excludes it
: the run draws every batch it is allowed. It therefore ends `converged=False`
with `note = "max_resamples_reached"`, and takes 6 batches rather than 1.
Two properties are asserted, neither pinning the growth factor. **Growth:** a
schedule that repeated the opening batch would total exactly
`batches x 100 = 600`, so a larger total shows the batches grew. **Budget:** the
first five batches are 100, 125, 156, 195, 243, cumulating to 819, short of the
1000 budget, so a sixth batch of 303 is drawn *whole*, ending at 1122. The
assertion is the bracket `1000 <= n_resamples < 2000`: the run does not stop
before the budget is met, and overshoots it by at most one batch. The exact
1122 is left unasserted because it encodes the growth factor.

### `test_bootstrap_resample_budget_shrinks_but_stays_usable_and_monotone`

**Inputs:** the configured pairs `(2500, 25000)`, `(2, 3)`, `(1, 1)`,
`(100, 100)`, `(7, 1000)`, `(999, 1001)`; then configured minima of 1, 2, 5,
10, 100, 2500 and 25000, each paired with ten times itself as the maximum.

**Derivation:** the expected values are not computed from the divisor, because
the divisor is a performance knob; the assertions are the three properties an
iteration budget has to satisfy whatever it is set to. *At least 1*: integer
division sends a small configured budget to 0, so the floor is what stops an
iteration drawing no resamples at all; `(2, 3)` and `(1, 1)` are the cases that
reach it. *Strictly below the configured minimum once that is at least 2*: any
divisor above 1 reduces such a value, so this is what would catch a divisor of 1
silently restoring the full per-iteration cost. *Bounds in order*: `(100, 100)`
and `(999, 1001)` are the cases where the scaled maximum would otherwise land
below the scaled minimum, leaving an adaptive run no range to grow through, and
the implementation holds the maximum at the minimum instead. A divisor of 1
fails five of the six pairs.

Floor division is monotonic and the floor is a constant, so both ends of the
scaled budget must be non-decreasing across the rising sequence. This is the
property a reader actually depends on, that configuring a larger budget
cannot give the bootstrap a smaller one, and it holds for any divisor.

## tests/orchestrator/test_orchestrator_trees.py

### `test_clean_and_save_trees_keeps_well_supported_trees_and_drops_the_rest`

**Inputs:** the species-tree fixture, then `low_support_tree_file`, each with
`min_avg_support=0.5`; the second contains

```
(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup);
(((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);
```

**Derivation:** the species tree carries no internal support labels, so there
is nothing for the filter to reject and nothing to strip. The output must
therefore be the input string plus a trailing newline:
`(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);\n`.

In the second file the mean internal support for tree 0 is
`(0.95 + 0.99 + 0.98)/3 = 2.92/3 = 0.9733 >= 0.5` → kept. For tree 1 it is
`(0.30 + 0.20 + 0.40)/3 = 0.90/3 = 0.3000 < 0.5` → dropped. Exactly one tree
survives, with leaf set `{TaxaC, TaxaD, TaxaF, TaxaG, OutGroup}`, and because
cleaning strips support labels the substring `0.95` must not appear in the
output.

### `test_clean_and_save_trees_quotes_labels_the_format_needs`

**Input:** `(('Homo sapiens':0.1,'Pan sp.':0.2):0.3,'O''Brien':0.4,Mus_musculus:0.5);`,
`min_avg_support=0.5`.

**Derivation:** Newick reserves whitespace and `()[]{}':;,` for structure, and
Bio.Phylo and DendroPy between them also trip on `"`, `\` and `=`, so a label
holding any of those must be single-quoted, with an inner quote written twice.
The three quoted input labels each hold such a character (a space, a space and
a dot, a quote) and so must come out quoted again (`'O''Brien'` with its
doubled quote intact) while `Mus_musculus` holds none and stays bare; with no
support labels to strip and the lengths already short, the file is the input
plus a newline, byte for byte. Read back, Bio.Phylo must give
`Homo sapiens`, `Pan sp.`, `O'Brien`, `Mus_musculus`, and DendroPy under
`preserve_underscores=True` the same four, or the run would be measuring
different taxa from the ones it wrote. Unquoted, the same file reads as
`sapiens`, `sp.`, `Brien` in Bio.Phylo and fails to parse in DendroPy.

### `test_root_species_tree_roots_where_the_outgroups_branch_off`

**Inputs (parametrized):** one of

- `(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);` with outgroups `["OUT"]`
- `(OUT1:0.3,(OUT2:0.2,(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4):0.5);` with `["OUT2", "OUT1"]`
- `(OUT1:0.3,OUT2:0.2,(A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6);` with `["OUT2", "OUT1"]`
- `(OUT1,(OUT2,(OUT3,((A,B),(C,D)))));` with `["OUT3", "OUT2", "OUT1", "OUTX"]`

**Derivation:** the rooting counts the outgroup leaves under every clade and
finds where outgroup-free subtrees hang off the paths joining the outgroups.
With a single outgroup the only such subtree is everything else, hanging off
`OUT` itself, so the tree is rerooted on `OUT`'s edge (see "Orchestrator
fixture: species tree" above): `OUT` is excluded, the ingroup is
`[A, B, C, D]`, and the pruned Newick is
`(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;` because the `((A,B),C)` clade
absorbs `0.1 + 0.2` when the `(D,OUT)` node dissolves.

In the second tree the file's root has `OUT1` on one side and `(OUT2, ingroup)`
on the other, so the outgroups are not a clade as written and their common
ancestor is the whole tree. The only outgroup-free subtree is the ingroup
clade `(((A,B),C),D):0.4`, hanging off the node that joins `OUT2` to it, so
that node is the single host: the tree is rerooted there (its children become
`OUT1`, `OUT2` and the ingroup clade), the two outgroups are pruned, and the
node collapses onto its one remaining child, which keeps its own edge:
`(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4;`. In the third tree the root is a
polytomy joining both outgroups and both ingroup clades; both clades hang off
that same node, so it is the host, rerooting is a no-op, and pruning the
outgroups leaves `((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6);` with the two clades
still under the root. The fourth tree roots the same way at the node joining
`OUT3` to `((A,B),(C,D))` and, having no branch lengths, prunes to
`((A,B),(C,D));`. In every case the ingroup is `[A, B, C, D]`.

Distances run from the ingroup root, the node the pruned tree is rooted at,
along the path to each outgroup. In the first tree that root is the dissolved
`(D,OUT)` node, `0.5` from `OUT`. In the second it is the ingroup clade's own
node, whose `0.4` edge starts both paths: `OUT2` is `0.4 + 0.2 = 0.6` away and
`OUT1`, past the file's root, `0.4 + 0.5 + 0.3 = 1.2`, so `OUT1` ranks first
although listed second. In the third the root polytomy is the ingroup root and
each outgroup's distance is its own edge, `0.3` and `0.2`. The fourth tree has
no lengths, so it has no distances (`None`) and the outgroups keep their listed
order, with `OUTX`, absent from the tree, last: `(OUT3, OUT2, OUT1, OUTX)`.

### `test_root_species_tree_rejects_outgroups_that_branch_off_twice`

**Inputs (parametrized):** outgroups `["OUT1", "OUT2"]` and one of

- `(((A:1,B:1):1,C:1):1,(D:1,(OUT1:1,(OUT2:1,X:1):1):1):1);`
- `((OUT1:1,(A:1,B:1):1):1,(OUT2:1,(C:1,D:1):1):1);`

**Derivation:** in the first tree the outgroups' lowest common clade is
`(OUT1,(OUT2,X))`; the taxa outside it, `A,B,C,D`, hang off that clade, while
`X` hangs off the lower node `(OUT2,X)`, whose branch lies between `OUT1` and
`OUT2`. Two hosts, so the rooting raises, and `separated_groups` lists the
groups largest first: `(("A","B","C","D"), ("X",))`. In the second tree the
common clade is the root, `(A,B)` hangs off the node joining it to `OUT1` and
`(C,D)` off the node joining it to `OUT2` (again two hosts) giving
`(("A","B"), ("C","D"))`, the tie broken by name.

### `test_generate_triplets_and_species_subtrees`

**Inputs (parametrized):** the orchestrator fixture species tree, once as
written and once with every branch length removed; either way it prunes to
`(((A,B),C),D)`.

**Derivation:** 4 ingroup taxa give `C(4,3) = 4` triplets, enumerated over the
sorted taxa: `(A,B,C)`, `(A,B,D)`, `(A,C,D)`, `(B,C,D)`. In each of these the
first two listed taxa are already the species-tree sister pair, so ABC
normalization is the identity and `triplets == raw_triplets`.

One cached geometry over the species tree answers all four triplets:
`triplet_subtree_shape` reports each one's sister pair and the order its induced
subtree lists the children in, and `_format_triplet_subtree_newick` writes the
topology: `((A,B),C);`, `((A,B),D);`, `((A,C),D);` and `((B,C),D);`. In the
pruned tree `(((A,B),C),D)` the sister clade comes first under every triplet's
LCA, and within it the leaves keep the tree's order. Only the topology is
written, so the tree with its branch lengths and the tree without them give the
same four strings.

### `test_read_tree_file_rejects_a_repeated_leaf_label`

**Inputs (parametrized):** `((A:1,B:1):1,(A:1,C:1):1);`, one tree with two
leaves labelled `A`; and a two-tree file whose second tree,
`((A:1,B:1):1,(B:1,C:1):1);`, has two leaves labelled `B`.

**Derivation:** the reader counts each tree's leaf labels and refuses any used
more than once, numbering trees from 1 in file order. The first row's message
names `A` (`uses the leaf label(s) A more than once`); in the second, tree 1 is
clean, so the message names `Tree 2`.

### `test_read_species_filter_file_collects_names_in_order`

**Input:** a file whose lines are ` A, B `, an empty line, `C`, `B,,D`, a line
of spaces, and `A`.

**Derivation:** each line is split on commas and every piece stripped, so the
first line yields `A` and `B` (surrounding spaces gone), the empty and
all-space lines yield nothing, `C` yields itself, `B,,D` yields `B`, an empty
piece that is dropped, and `D`, and the last line yields `A`. Repeats keep
their first position: `B` and `A` are already present, so the result is
`["A", "B", "C", "D"]` in first-seen order.

### `test_clean_and_save_gene_trees_roots_each_tree_from_its_farthest_outgroup`

**Inputs (parametrized over `processes` 1 and 2):** `orchestrator_gene_trees`
(12 trees, outgroup `OUT`); then a gene-tree file of nine trees, checked
against `["OUT1", "OUT2", "OUT3"]`, the species-tree rank:

1. `((((A:1,B:1):1,C:1):1,D:1):1,OUT1:1);`: `OUT1` alone.
2. `(((A:1,B:1):1,(C:1,OUT2:1):1):1,OUT1:1);`: `OUT1` at the root, `OUT2`
   sister to `C`.
3. `((((A:1,B:1):1,(C:1,OUT1:1):1):1,D:1):1,(OUT2:1,OUT3:1):1);`: `OUT1`
   sister to `C`, `OUT2` and `OUT3` sisters.
4. The same as tree 3 with `OUT1:9`.
5. `((((A:1,B:1):1,C:1):1,D:1):1,OUT2:1);`: `OUT2` alone.
6. `(((A,B),(C,OUT2)),OUT1);`: no branch lengths.
7. `(((A:1,B:1),(C:1,OUT2:1):1):1,OUT1:1);`: the `(A,B)` edge lacks a length.
8. `(((A:1,B:1):1,C:1):1,D:1);`: no outgroup.
9. `(OUT1:1,(OUT2:1,OUT3:1):1);`: nothing but outgroups.

**Derivation:** none of the 12 fixture trees carry support labels, so all 12
survive; each carries `OUT` alone, which is the farthest and roots it
(`farthest == rooted_on == {"OUT": 12}`) with nothing to tangle
(`tangled == {"OUT": 0}`, `tangled_trees == 0`). Rooting on `OUT` moves its edge into the ingroup
clade's edge, and pruning `OUT` leaves that clade as the root, so no line
carries `OUT`. Tree 0: `0.10 + 0.50 = 0.6` →
`(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6;`. Tree 3: `0.10 + 0.55 = 0.65` →
`(((B:0.4,C:0.4):0.1,A:0.6):0.1,D:0.35):0.65;`.

For the nine trees, the farthest outgroup has the longest mean path to the
ingroup taxa; an outgroup under a child of the ingroup's common ancestor that
also holds ingroup taxa, once the tree is rooted on the farthest, is tangled;
and the tree is rooted where the remaining outgroups part from the rest.

- Tree 1 carries `OUT1` alone, so it is the farthest and roots the tree;
  its edge 1 folds into the ingroup edge 1: `(((A:1,B:1):1,C:1):1,D:1):2;`.
- Tree 2: `OUT1` is 4 from each of `A`, `B` and `C` (mean 4); `OUT2` is 4
  from `A` and `B` and 2 from `C` (mean 10/3). Rooted on `OUT1`, the
  ingroup's common ancestor is `((A,B),(C,OUT2))`, and `OUT2` sits under its
  child `(C,OUT2)` with `C`, so it is tangled. Pruning it leaves `C` on
  `1 + 1 = 2`: `((A:1,B:1):1,C:2):2;`.
- Tree 3: `OUT1` is 4 from `A`, `B` and `D` and 2 from `C` (mean 3.5);
  `OUT2` is 6 from `A`, `B` and `C` and 4 from `D` (mean 5.5), and so is
  `OUT3`. The tie keeps the species-tree rank, so `OUT2` is the farthest
  although the species tree ranks `OUT1` first. Rooted on `OUT2`, `OUT3` is
  its sister, outside the ingroup, and `OUT1` sits with `C`, so it is
  tangled. The tree is rooted where the pair joins the ingroup, whose edge
  takes the pair's edge (`1 + 1 = 2`), and pruning `OUT1` lengthens `C`'s
  edge to 2: `(((A:1,B:1):1,C:2):1,D:1):2;`.
- Tree 4: `OUT1`'s branch of 9 puts it 12 from `A`, `B` and `D` and 10 from
  `C` (mean 11.5), farther than the pair's 5.5. Rooted on `OUT1`, the
  ingroup's common ancestor is the node `OUT1` shared with `C`, and the pair
  sits under its other child with `A`, `B` and `D`, so both are tangled.
  Rooting on a leaf gives the leaf's branch to the other side, so the root
  keeps 9; pruning the pair dissolves `D`'s node, and `D` takes `1 + 1 = 2`:
  `((D:2,(A:1,B:1):1):1,C:1):9;`.
- Tree 5 carries `OUT2` alone: `(((A:1,B:1):1,C:1):1,D:1):2;`.
- Tree 6 has no branch lengths, so every mean path is 0 and the tie keeps
  the species-tree rank: `OUT1` is the farthest, and `OUT2`, sitting with
  `C`, is tangled. Rerooting gives the missing root edge the value 0, and
  the rest stay unwritten: `((A,B),C):0;`.
- Tree 7 lacks only the `(A,B)` edge, read as 0. `OUT1` is `1 + 1 + 0 + 1 = 3`
  from `A` and `B` and 4 from `C` (mean 10/3); `OUT2` is 3 from `A` and `B`
  and 2 from `C` (mean 8/3). So `OUT1` roots it and `OUT2` is tangled. The
  root edge is the ingroup's 1 plus `OUT1`'s 1, and `C` takes `OUT2`'s parent
  edge: `((A:1,B:1),C:2):2;`.
- Trees 6 and 7 are kept and counted: `missing_length_indices == [6, 7]`.
- Trees 8 and 9 are dropped: `unrootable_indices == [8, 9]`.

Hence `rooted_count == 7`, `farthest == {"OUT1": 5, "OUT2": 2, "OUT3": 0}`
(trees 1, 2, 4, 6 and 7; trees 3 and 5), `rooted_on == {"OUT1": 5,
"OUT2": 2, "OUT3": 1}` (tree 3 counts for `OUT2` and `OUT3`), `tangled ==
{"OUT1": 1, "OUT2": 4, "OUT3": 1}` (`OUT1` in tree 3; `OUT2` in trees 2, 4, 6
and 7; `OUT3` in tree 4) and `tangled_trees == 5`. No tree carries support values, so none is
dropped for support.

## tests/orchestrator/test_orchestrator.py

All runs use the shared species tree and 12 gene trees, seed `20240724`, 40
bootstrap iterations, and the orchestrator defaults (`p_value_correction = bfn`,
the permutation test enabled).

### `test_run_orchestrator_matches_derived_expectation`

**Expected-output derivation**, from `_EXPECTED_COUNTS`:

| Triplet | n_con | n_dis1 | n_dis2 | species_tree |
| --- | --- | --- | --- | --- |
| (A,B,C) | 7 | 3 | 2 | `((A,B),C);` |
| (A,B,D) | 12 | 0 | 0 | `((A,B),D);` |
| (A,C,D) | 12 | 0 | 0 | `((A,C),D);` |
| (B,C,D) | 12 | 0 | 0 | `((B,C),D);` |

(see "Orchestrator fixture: gene trees" for how the counts are read off).

- **DCT for (A,B,C):** `chisquare([3, 2])` → statistic `0.2`, `p ~ 0.6547`.
- **DCT for the D triplets:** `n_dis1 + n_dis2 == 0`, so the short circuit gives
  `(0.0, 1.0)`.
- **Bonferroni correction:** 4 triplets, so `p_corrected = min(1.0, p x 4)`.
  `0.6547 x 4 = 2.62 → 1.0`; `1.0 x 4 → 1.0`. Every corrected DCT p-value is
  `1.0`, far above `alpha_dct = 0.05`, so `dct_significant is False` and every
  triplet classifies as `no_introgression`.
- **KS correction, `diagnostic` on only:** asserted as a relation rather than a
  literal: `ks_p_value_corrected == min(1.0, ks_p_value x 4)`. (For (A,B,C)
  the raw KS p-value is `~0.01667`, giving `~0.06667`; for the D triplets dis1
  is empty so the KS test short-circuits to `1.0` and stays `1.0`.)
  `perm_decision` is populated on every row because the direction test ran.
- **`diagnostic` off:** the run's default correction is `bfn`, an inline
  method, so the stream judges each count gate on the exactly corrected
  `p × 4 = 1.0`, finds it failed, and measures nothing below it: `ks_p_value`,
  `ks_p_value_corrected`, `ks_significant` and `perm_decision` are `None` and
  `perm_note` is `direction_test_not_consulted`. Everything asserted above the
  gate is identical between the two runs.
- `most_frequent_matches_concordant` is `True` everywhere, since `7 >= 3, 2` and
  `12 >= 0, 0`.
- `analyzed_trees == 12` for every triplet, because all 12 gene trees contain
  all five taxa.

### Output-shape tests

- `test_run_orchestrator_matches_derived_expectation` also reads the results
  TSV the run wrote: the file must exist, its header must begin with `triplet`
  and contain `classification`, `bootstrap_value`, `perm_p_greater`,
  `perm_p_less`, `decision_gate`, `perm_p_tost` and the two
  `bootstrap_perm_stat_ci_*` columns, and the data-row count must equal
  `len(results)` (4).
- `test_no_bootstrap_skips_the_bootstrap_and_its_columns`: with
  `bootstrap=False` no iteration runs, so every result's `bootstrap_value`,
  `all_bootstrap` and both `bootstrap_perm_stat_ci_*` bounds are `None` while
  the point estimate still classifies `no_introgression` (every triplet stops
  at the count gate); the writer skips the bootstrap columns, so
  `bootstrap_value` and `all_bootstrap` must be absent from the header while
  `classification` remains and all 4 triplets are produced. The studentized
  interval is produced only by the bootstrap loop, so its absence is what
  shows the loop did not run; `diagnostic` on measures the tests the cascade
  cannot consult and must not reinstate the bootstrap, so both settings give
  the same picture.
- `test_run_outputs_follow_the_settings`: one run with consolidation, summary
  statistics and the bootstrap diagnostic on. Consolidation writes into
  `consolidation/`, so the run's own files (results TSV, `metrics.txt`, both
  processed trees) must all still exist afterwards and the subfolder must be
  non-empty. The summary column count is derived from the contract: 7
  statistics (mean, median, mode, variance, entropy, min, max) x 3 metrics
  (avg_tree_height, internal_branch, sister_distance) x 3 topology classes
  (concordant, discordant1, discordant2) = **63** metric columns. The nine
  named per-iteration bootstrap columns must appear in the results header and
  at least one result must have a populated `bootstrap_dct_stats` and
  `bootstrap_perm_decisions`.
- `test_parallel_runs_match_serial`: bootstrap seeding is per-triplet and
  derived from the run seed, so the worker count (2 or 4) cannot change any
  value; every compared field must be equal to the serial run's, bootstrap
  included.

## tests/orchestrator/test_orchestrator_preflight.py

### Shared inputs

The species tree is `(((A,B),C),(D,OUT));`, without branch lengths, which the
check needs no more than a run does. Removing the
outgroup `OUT` leaves ingroup `{A, B, C, D}`, so the check enumerates
`C(4,3) = 4` triplets: `A,B,C`, `A,B,D`, `A,C,D`, `B,C,D`.

The clean gene-tree file holds two trees, both containing all five taxa:

```
((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);
((((A:1,C:1):1,B:1):1,D:1):1,OUT:1);
```

The defective file holds four, each planted with exactly one problem class:

```
((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);   # well formed
(((A:1,B:1,C:1):1,D:1):1,OUT:1);       # polytomy over A,B,C
(((A:1,B:1):1,C:1):1,MISSING:1);       # no OUT label
((((A,B),C),D),OUT);                   # no branch lengths
```

### `test_clean_inputs_pass_and_the_report_lands_where_documented`

**Inputs:** the clean pair above, `outgroups=["OUT"]`, checked once with
`output_dir` set to a created directory and once with `output_dir=None`.

**Derivation:** both gene trees root on `OUT`, so
`counters["gene_tree.rooted"] == 2` and `gene_tree.total_checked == 2`. Every
one of the 4 triplets resolves a sister pair in the species tree, so
`triplets_checked == 4`. Each gene tree is fully resolved, so no triplet check
fails and `issues == []`, which makes `passed` `True`.

The writer joins `output_dir` with the module constant
`PREFLIGHT_REPORT_FILENAME` (`preflight_data_check.txt`) and writes
`report_text` verbatim, so the file content and `report_text` must be equal and
`report_path` must equal that joined path. The write branch is guarded on
`output_dir is not None`, so the second call leaves `report_path` `None` while
still building the same `report_text`: the check itself does not depend on
where its output goes.

### `test_species_tree_is_rooted_where_the_outgroups_branch_off`

**Inputs:** the species tree `(OUT1:1,(OUT2:1,(((A:1,B:1):1,C:1):1,D:1):1):1);`,
`outgroups=["OUT2", "OUT1"]`, and three gene trees:
`((((A:1,B:1):1,C:1):1,D:1):1,OUT1:1);`,
`(((A:1,B:1):1,(C:1,OUT2:1):1):1,OUT1:1);` and
`((((A:1,C:1):1,B:1):1,D:1):1,OUT2:1);`.

**Derivation:** `OUT1` and `OUT2` sit on either side of the file's root, so they
are not a clade as written; the only outgroup-free subtree is the ingroup
clade `(((A,B),C),D)`, hanging off the node that joins `OUT2` to it, so the
tree is rooted there, both outgroups are pruned and the ingroup is
`A, B, C, D`, giving `C(4,3) = 4` triplets. From the ingroup root `OUT2` is
`1 + 1 = 2` away and `OUT1`, past the file's root, `1 + 1 + 1 = 3`, so `OUT1`
ranks first although listed second. Gene tree 1 carries `OUT1` alone and
tree 3 `OUT2` alone, so each is the farthest in its tree and roots it. In
tree 2 `OUT1`'s mean path to `A`, `B` and `C` is 4 against `OUT2`'s 10/3, so
`OUT1` is the farthest; rooted on it, `OUT2`, `C`'s sister, sits among the
ingroup and is set aside: `gene_tree.rooted == 3`, `farthest.OUT1 == 2`,
`farthest.OUT2 == 1`, `rooted_on.OUT1 == 2`, `rooted_on.OUT2 == 1`,
`tangled.OUT2 == 1` and `tangled_trees == 1`. Tree 2 lacks
`D`, so its three `D` triplets are skipped as absent
(`triplet.taxa_absent_from_gene_tree == 3`); the other `4 + 1 + 4 = 9` pairs
resolve (tree 2's `A,B,C` has `A,B` as sisters with `C` outside), so
`triplet.resolved == 9`, no issue is raised and `passed` is `True`.

### `test_detects_polytomy_and_missing_outgroup`

**Inputs (parametrized over `processes` 1 and 2):** the defective four,
`outgroups=["OUT"]`.

**Derivation:** gene tree 3 contains no `OUT`, so it has no outgroup to root
on → one `gene_tree.rooting_failed`, and that tree is skipped before any
triplet check. Gene tree 4 has no branch lengths, which a run reads as 0: it
roots and is counted (`gene_tree.missing_branch_lengths == 1`) without an
issue, and it is the 4th tree in the file, so `missing_length_indices == [4]`
and the report lists it as
`Gene trees lacking some branch length (1; 1-based, in input file order):`
followed by `  4`. Each tree's checks depend on that tree alone, so two
workers give the same counters, issues and listing as one. Trees 1, 2 and 4 root, so `gene_tree.rooted == 3` while
`gene_tree.total_checked == 4`. Tree 2 collapses
A, B and C into a single polytomous clade, so all three pairwise LCAs of triplet
`A,B,C` are the same node and `triplet_resolution` reports it unresolved → one
`triplet.unresolved_rooted_sister_pair`. The other three triplets each contain
`D`, which sits outside the polytomy, so they still resolve, hence a count of
exactly 1, not 4. The message is formatted with the enumeration index (`Gene
tree #2`) and the comma-joined triplet (`A,B,C`).

The pair accounting follows from the same reading. Three trees root, and each
carries all four ingroup taxa, so the check looks at `4 x 3 = 12` triplet/
gene-tree pairs: trees 1 and 4 resolve all 4 each, tree 2 resolves the three
containing `D` and fails on `A,B,C`. That gives 11 usable, 1 unresolved, and 0
skipped for an absent taxon, and the three must sum to the 12 pairs seen, which is the
property worth pinning: a pair that is neither measured nor reported would
otherwise vanish silently between the counters.

### `test_filter_entries_are_validated`

**Inputs (parametrized):** the clean pair plus either a triplet filter file
containing `A,B,C`, `A,B,NOPE`, `A,B,OUT`, or a species filter file whose
lines are `A,B`, `C`, `NOPE`, `OUT`.

**Derivation:** in the triplet filter `NOPE` is not among the species-tree
labels, so set difference against them is non-empty → one
`triplet_filter.taxa_missing_in_species_tree`. `OUT` is in the outgroup set →
one `triplet_filter.includes_outgroup`. Both lines are dropped rather than
checked, leaving only `A,B,C`, which normalizes successfully, so
`triplets_checked == 1`.

The species filter yields `A`, `B`, `C`, `NOPE`, `OUT`. `NOPE` is not among
the species-tree labels → one `species_filter.taxa_missing_in_species_tree`;
`OUT` is a label but in the outgroup set → one
`species_filter.includes_outgroup`. The three usable species `A`, `B`, `C`
form `C(3, 3) = 1` triplet, against the `C(4, 3) = 4` the unfiltered tree
would give, and it normalizes successfully, so `triplets_checked == 1`.

### `test_impossible_checks_raise`

**Inputs (parametrized):** `outgroups=[]`, `["NOT_PRESENT"]`,
`["A", "B", "C", "D", "OUT"]` and `["OUT", "C"]`, each against the clean
species tree `(((A,B),C),(D,OUT))`; and `["OUT"]` against that species-tree
string written twice.

**Derivation:** the empty list fails the explicit guard at the top of
`run_preflight_data_check`. The next three fail inside the run's own species
rooting, whose `OutgroupRootingError` is a `ValueError` and propagates: with
`NOT_PRESENT` no outgroup is in the tree; with all five taxa named no ingroup
remains; with `OUT` and `C`, `(A,B)` hangs off the node joining `C` to the
tree and `D` off the node joining `OUT`, so the outgroups branch off at two
points and the message lists the groups `[2 taxa: A, B]; [1 taxon: D]`, which
the test matches on `1 taxon: D`. Without a rooted species tree there is no
triplet normalization to perform, so there is nothing to report on. For the
doubled file `read_tree_file` returns 2 trees, and `_load_single_species_tree`
requires exactly 1.

### `test_runner_preflight_mode_skips_analysis`

**Inputs (parametrized):** a config dict with `preflight_data_check: True`
and `preflight_triplet_cap: 3` over the defective four, with
`outgroup="OUT"` and then `outgroup="NOT_PRESENT"`.

**Derivation:** `run_orchestrator` prepares the output directory and then
returns `_run_preflight_only(...)` before any tree cleaning, so the only write
into that directory is the report. Listing the directory must therefore yield
exactly `["preflight_data_check.txt"]`: no `metrics.txt`, no processed trees,
no results TSV. `passed` is `False` because the defective four yield issues.

The cap: the species tree's ingroup is `A, B, C, D`, so `combinations(·, 3)`
yields `C(4, 3) = 4` triplets. `_load_target_triplets` keeps the first
`max_triplets` when `0 < max_triplets < 4`, so a cap of `3` gives
`triplets_checked == 3` and appends one `analysis.triplet_cap_applied` issue.
Had the runner ignored the config and used the module default (15,000), all 4
would be checked and the cap issue would be absent, so both assertions would
fail, which is what makes `3` rather than the default the right value here.

With `NOT_PRESENT` the species tree cannot be rooted, so
`run_preflight_data_check` raises `OutgroupRootingError`, an `InputError`,
before any report is written; the runner lets it propagate, and the prepared
output directory stays empty.

### `test_index_ranges_collapse_consecutive_runs`

**Inputs (parametrized):** `[4]`, `[3, 17, 18, 19, 250]`, `[1, 2, 4, 5, 6]`.

**Derivation:** a run of consecutive positions is written `first-last` and a
position with no neighbour on its own, joined by `, `: `4`; `3`, then
`17, 18, 19` as `17-19`, then `250`; `1, 2` as `1-2` and `4, 5, 6` as `4-6`.

## tests/orchestrator/test_orchestrator_config.py

No test here pins an individual default. Two invariants stand in for all of
them, and the remaining tests cover precedence, parsing, validation, and the
shipped samples.

### `test_cli_and_config_file_share_one_set_of_defaults`

**Inputs:** `build_argument_parser().parse_args(["-st", "species.tree", "-gt",
"genes.tree", "-og", "OUT"])`, and a YAML file carrying only those three keys.

**Derivation:** `resolve_config` on a CLI namespace and `load_orchestrator_config`
on a file both end in `normalize_orchestrator_payload`, which fills every absent
key from the constants in `ghostparser/orchestrator/config.py`. With the same
three inputs and nothing else, the two resolved dicts are therefore equal in
every key, including the resolved paths, since both resolve the same relative
strings against the same working directory. The comparison is whole-dict
equality, so a key that resolved differently on the two paths, or a default
that one path filled and the other did not, fails without the test naming
either. It was checked before being written that the two dicts are in fact
equal today.

### `test_config_only_keys_and_nested_blocks_flatten_from_a_file`

**Input:** a YAML file with `discordant_test: z-test`,
`tree_height_calculation_strategy: SIS`, `min_support_value: 0.9`,
`generate_summary_stats: true`, `alpha_dct: 0.02`, `seed: 7`, and
`bootstrap_options: {iterations: 25, diagnostic: true, summary_only: true}`.

**Derivation:** the first four keys have no CLI flag, so the file is the only
way to set them, and each is read back as written. The nested block is
flattened onto `bootstrap_iterations = 25`, `bootstrap_diagnostic = True`,
`bootstrap_summary_only = True`. The nine values are compared as one dict
against the nine written, so any one of them resolving to something else fails
naming the key.

### `test_cli_flags_override_the_config_file`

**Input:** a config file setting file-specific tree paths, `alpha_dct: 0.03`,
`p_value_correction: "no"` and `overwrite: true`, plus the flags
`alpha_dct=0.5`, `alpha_perm=0.5`, `no_overwrite=True` and
`consolidation=False` with the path flags left unset.

**Derivation:** the given flags are laid over the file's payload before it
is normalized, so `alpha_dct` is the flag's `0.5` and `alpha_perm`, which the
file omits, is the flag's `0.5` rather than the default. `--no-overwrite` is
translated to `overwrite: false` before the merge, so it replaces the file's
`overwrite: true` instead of losing to it in the overwrite validator, which
prefers the canonical key. `--no-consolidation` gives `consolidation: false`
the same way. `p_value_correction` and the paths have no flag given, so they
keep the file's values.

### `test_outgroup_key_is_normalized_or_rejected`

**Inputs and derivation:** every accepted shape flattens the same way: each
entry is coerced to a string, split on commas, stripped, and empty pieces
dropped; a value that leaves no label raises `ConfigError` naming `outgroup`
rather than silently producing an unrooted run.

| Input | Result | Why |
| --- | --- | --- |
| `"OUT"` | `["OUT"]` | Single label, no comma to split on. |
| `"Out1,Out2"` | `["Out1", "Out2"]` | Comma-separated string. |
| `" Out1 , Out2 ,"` | `["Out1", "Out2"]` | Padding stripped, trailing empty piece dropped. |
| `["Out1,Out2", "Out3"]` | `["Out1", "Out2", "Out3"]` | List entries are themselves split, so the two forms compose. |
| `None` | `ConfigError` | No value at all. |
| `"  "` | `ConfigError` | Every piece is blank after stripping. |
| `[""]` | `ConfigError` | The one entry is blank. |
| `42` | `ConfigError` | Neither a string nor a list/tuple/set, so it contributes no entries. |

### `test_invalid_values_are_rejected_by_field_name`

**Inputs:** the required-keys file with one key overridden per row.

| Override | Validator reached | Why it fails | Message must contain |
| --- | --- | --- | --- |
| `species_tree_path: null` | `_validate_required_path` | `None` is treated as absent. | `species_tree_path` |
| `p_value_correction: true` | `_validate_choice` | YAML 1.1 turns a bare `yes` into `True`; the boolean-to-choice mapping covers `False` → `no` but `True` matches no choice. Written as `true` here since the payload is dumped with `yaml.safe_dump`. | `must be one of` |
| `diagnostic: "yes"` | `_validate_optional_bool` | Accepts `None` or a `bool`; a string is rejected even when it spells a boolean, since the YAML forms already resolve to one. | `diagnostic` |
| `seed: "abc"` | `_validate_optional_int` | Accepts `None` or an `int`; also rejects `bool`, since `isinstance(True, int)` holds and `seed: true` is a mistake rather than a seed of 1. | `seed` |
| `preflight_triplet_cap: -1` | `_validate_non_negative_int` | Accepts `int >= 0`; `0` is the documented "no cap", so `-1` is the nearest value with no meaning. | `preflight_triplet_cap` |
| `species_rename_map: absent.tsv` | `_validate_species_rename_map` | The path resolves, but the file is read on the spot (`load_species_rename_map`) and does not exist; its `FileNotFoundError` is re-raised as a config error. | `species_rename_map` |
| `triplet_filter: triplets.txt` beside `species_filter: species.txt` | `_validate_triplet_selection` | Each path resolves on its own (neither file is read here), but the two select triplets in ways that cannot be reconciled (named triplets against every triplet among named species), so both set is an error before any file is opened. | `cannot both be set` |

The `p_value_correction` row is the guard that keeps the boolean mapping from
laundering an invalid value into a valid one.

### `test_shipped_sample_configs_resolve`

**Inputs:** the six orchestrator sample configs under `sample_configs/`.

**Derivation:** these go through the same `load_orchestrator_config` a user
invokes with `-c`, so anything the validator would reject surfaces here. Each
asserted value is what its sample states literally: the outgroup lists
`["OutGroup"]` (minimal) and `["Out1", "Out2"]` (full, written as a YAML list
to exercise that form), `preflight_data_check: true` (preflight),
`generate_summary_stats: true` (species filter), `diagnostic: true` (triplet
filter) and `"bootstrap": false` (screen, JSON). The placeholder tree paths
need not exist, since `_validate_required_path` only checks that the field is
a non-empty string before resolving it, and neither need the filter files,
which are read at run time; a rename map is read at load time, which is why
the species-filter sample carries it commented out.

### `test_full_sample_config_names_every_runtime_key_at_its_default`

**Inputs:** `orchestrator_full.yaml` parsed twice: once as raw YAML for the set
of documented keys, once through `load_orchestrator_config` for the set of
runtime keys, plus a required-keys-only file loaded the same way.

**Derivation:** two checks. First, coverage: the normalizer renames three inputs
(`species_tree_path` → `species_tree`, `gene_trees_path` → `gene_trees`,
`output_folder` → `output`) and flattens the two nested blocks with a prefix
(`permutation_options.ci_method` → `permutation_ci_method`). After undoing both
transformations the runtime key set must be a subset of the documented one, so
the difference is empty. This fails the moment a config key is added without
the sample gaining it. Second, values: outside the input keys (the three
renamed paths and `outgroup`), the sample must resolve to exactly what the
required-keys-only file resolves to. Coverage alone would pass if the sample
omitted a key (it would take the default and compare equal), and equality alone
would pass if the sample listed a key at a non-default value it then also
omitted elsewhere; together they pin that every key is both present and at its
default. It was checked that the sample satisfies this today.

### `test_p_value_correction_accepts_yaml_bare_word_no`

**Inputs:** `p_value_correction` written as the bare word `no`, as quoted
`"no"`, and as `NO`.

**Derivation:** YAML 1.1 resolves the bare words `no` and `NO` to boolean
`False`, so the normalizer maps `False` back to the string `"no"` before the
choice check; quoted `"no"` arrives as a string and passes straight through.
All three therefore resolve to `"no"`. Writing `no` unquoted is the natural
spelling for "no correction", which is why the mapping exists.

### `test_parser_flags_resolve_into_their_config_values`

Parsing `-st s -gt g -og OUT --alpha-dct 0.01 --alpha-ks 0.2
--p-value-correction fdr_bh --diagnostic --species-filter species.txt
--alpha-perm 0.02 --no-overwrite --preflight-data-check
--preflight-triplet-cap 0` must yield those exact values with
`config_file is None`; `0` is chosen for the cap because it is the one value
`_validate_non_negative_int` accepts that differs from the default and is also
the documented "no cap" spelling. `--species-filter` goes through
`_resolve_path`, so only the tail of the resolved path is pinned, and
`triplet_filter` must stay `None` because the flag was not given.

## tests/test_cli.py

### `test_run_cli_maps_each_outcome_to_its_exit_status`

**Inputs (parametrized):** a bare `argparse.ArgumentParser` and a `run`
callable that returns `None` or `EXIT_CHECK_FAILED`, or raises one of
`ConfigError`, `InputError`, `OutgroupRootingError`, `TypeError` or
`KeyboardInterrupt`; parsed from `[]` and, for the raising rows, again from
`["--debug"]`.

**Derivation:** the wrapper adds `--debug` to the parser, so both argument
lists parse. A `None` return is success, `0`; a returned status passes
through, so `EXIT_CHECK_FAILED` stays `3`. The handlers are tried from the
most specific: `ConfigError` → `2`; any other `GhostParserError` → `1`, which
`OutgroupRootingError` reaches through `InputError`; `KeyboardInterrupt`,
which is not an `Exception`, has its own handler → `130`; and `TypeError`,
outside the package's errors, is a bug → `70`. Each failure is printed as one
line, so stderr holds no `Traceback`. With `--debug` every handler re-raises,
so the original exception type leaves the wrapper.

## tests/test_config_trunk.py

### `test_resolve_path_handles_absolute_relative_and_home`

**Inputs:** an absolute `tmp_path/species.nwk`; the relative `genes.nwk` after
`chdir` into `tmp_path`; and `~/data.nwk`.

**Derivation:** absolute paths resolve to themselves; relative paths resolve
against the current working directory, giving `tmp_path/genes.nwk`; `~` expands
to `Path.home()`, so the result is `Path.home()/data.nwk` and contains no `~`.

### `test_load_raw_config_reads_json_and_yaml_and_rejects_the_rest`

**Inputs:** `{"input_path": "data.tsv", "cv_folds": 5}` written as JSON, and
the equivalent two-line YAML under `.yaml` and `.yml`; then a nonexistent
path, a `.txt` file, and a JSON file whose root is `["a","b"]`.

**Derivation:** all three formats must parse to the identical Python mapping:
the loader's only job is format dispatch. For the rejections the existence
check runs first → `ConfigError("Config file not found: ...")`; the suffix
check runs next →
`ConfigError("Config file must be .json, .yaml, or .yml")`; the root-type
check runs last → `ConfigError("Config root must be a key/value object")`.

### `test_validate_required_path_resolves_or_raises`

**Inputs:** `{"species_tree_path": "s.nwk"}` under a chdir'd cwd, then `{}`,
`{"species_tree_path": ""}`, and `{"species_tree_path": "   "}`.

**Derivation:** a non-empty string resolves to an absolute path. Absent, empty,
and whitespace-only values all fail the `not value.strip()` guard → three
`ConfigError`s.

### `test_validate_overwrite_flag_precedence_and_validation`

**Inputs / derivation**, straight from the precedence rule (canonical
`overwrite` first, then negated `no_overwrite`, then the default):

| Payload | Expected | Reason |
| --- | --- | --- |
| `{}` | `DEFAULT_OVERWRITE` (True) | neither key present |
| `{"overwrite": False}` | `False` | canonical key used directly |
| `{"no_overwrite": True}` | `False` | negated |
| `{"no_overwrite": False}` | `True` | negated |
| `{"overwrite": True, "no_overwrite": True}` | `True` | canonical wins |
| `{}` with `default=False` | `False` | explicit default |
| `{"overwrite": "yes"}` | `ConfigError` | not a boolean |
| `{"no_overwrite": "yes"}` | `ConfigError` | not a boolean |

### `test_prepare_output_directory_overwrites_or_suffixes`

**Inputs:** `tmp_path/results` containing `stale.txt`; then the same directory
containing `fresh.txt`, with siblings `results_1` and `results_3` present, called
with `overwrite=False`.

**Derivation:** with `overwrite=True` (the default) the directory is removed and
recreated, so `stale.txt` must be gone. With `overwrite=False` the allocator
scans the parent, finds suffixes `{1, 3}` in use, and returns the smallest
missing positive suffix, `2`, creating `results_2` and leaving `fresh.txt`
untouched in the original directory.

## tests/test_ml_labels_and_metrics.py

### `test_bit_labels_and_their_titles_cover_every_bit`

**Inputs:** the module constants `BIT_LABELS` and `BIT_COUNT`; the six
`BIT_LABELS` entries are `ghost_into_A`, `ghost_into_B`,
`inflow_into_A_from_C`, `inflow_into_B_from_C`, `outflow_from_A_to_C`,
`outflow_from_B_to_C`.

**Derivation:** `BIT_COUNT` is defined as `len(BIT_LABELS)`, and the six
names are distinct, so a bit index reads back to exactly one name. The title
transform is `replace("_", " ")` followed by upper-casing character 0 only, so
`inflow_into_A_from_C` → `inflow into A from C` → `Inflow into A from C`. The
expected titles are written out per label rather than computed, because the
point is the one spelling the obvious implementation gets wrong:
`str.capitalize()` upper-cases the first character *and lower-cases the
rest*, which would yield `Ghost into a` and rename taxon `A`. Every label in
the set carries at least one trailing capital, so any label would catch it;
they are all listed so the failure names which one broke, and comparing the
whole mapping at once also shows six distinct titles, none keeping an
underscore, so no label falls through to a figure as a raw slug.

### `test_64_class_matrix_orders_classes_by_set_bits`

**Inputs:** `y_true` rows `[1,1,0,0,0,0]`, `[0,0,0,0,0,1]`, `[1,1,1,1,1,1]` and
`y_pred` rows `[0,0,0,0,1,1]`, `[0,0,0,0,0,1]`, `[1,1,1,1,1,1]`: classes
`110000 → 000011`, `000001 → 000001`, `111111 → 111111`.

**Derivation:** the builder sorts the 64 six-bit strings on
`(count("1"), label)`, so the list opens `000000`, runs the six single-bit
classes `000001 … 100000`, then the fifteen two-bit classes, and so on to
`111111`: set-bit counts `[1, 6, 15, 20, 15, 6, 1]` (the binomial row for
`n = 6`) and never decreasing. Within a count the tie-break is the string
itself, so each group is lexically sorted. The three rows were chosen so that
each of the three (true, predicted) pairs sits somewhere binary order and
set-bit order disagree: in binary order `110000` is index 48 and `000011` index
3, but after the sort they are at positions 21 and 7. The matrix is indexed by
`labels.index(...)` on both axes, so the assertion passes only if rows *and*
columns were permuted with the labels; a sort that reordered the label list but
left the counts at their binary indices would report the `110000 → 000011`
count at position `(48, 3)`, which the reordered labels name `(101011, 000100)`. The
total of 3 confirms nothing was dropped or duplicated by the reindexing.

### `test_bitstring_labels_are_validated_and_parsed`

**Inputs (parametrized):** `000000`, `111111`, `100001`, ` 010010 `, `10000`,
`1000010`, `100002`, `10000a` and `""`.

**Derivation:** validity is `len(value) == 6 and set(value) <= {"0","1"}`.
Hence `000000`, `111111`, `100001` pass, as does ` 010010 ` once stripped;
`10000` (5 chars), `1000010` (7 chars), `100002` and `10000a` (non-binary),
and `""` all fail. A valid label is stripped and expanded into six integers,
giving a `1 x 6` matrix whose one row re-joins to the trimmed label and whose
sum is the label's count of `1`s: `100001 → 2`, `000000 → 0`, `111111 → 6`,
`010010 → 2`. A malformed label raises `ValueError("Invalid classes label
...")` even when it follows the valid `100001` in the same list, so one bad
row cannot hide behind good ones.

### `test_prediction_metrics_and_rows_match_hand_computation`

**Inputs:**

```
y_true = [[1,0,0,0,0,1],
          [0,1,0,0,0,0]]
y_pred = [[1,0,0,0,0,1],
          [0,0,0,0,0,0]]
```

**Derivation:** row 0 matches exactly; row 1 differs at bit index 1 only.

- `exact_match_accuracy` = 1 of 2 rows = `0.5`.
- Total bit positions = `2 x 6 = 12`, mismatches = 1, so
  `bitwise_accuracy = 11/12` and `hamming_loss = 1/12`.
- `ghost_into_A` (bit 0) is predicted correctly in both rows → accuracy `1.0`.
- `ghost_into_B` (bit 1) has one true positive that was predicted 0 → recall
  `0/1 = 0.0`, support `1`.

For the prediction rows, `matched_label_count` counts positionwise agreement:
row 0 → 6, row 1 → 5. `exact_match` is 1 only for row 0. Label strings are the
joined bits: row 0 `100001`/`100001`; row 1 `010000` true vs `000000`
predicted. The per-bit columns expose that difference as
`true_ghost_into_B = 1`, `pred_ghost_into_B = 0`.

### `test_distributions_count_labels_and_positives_per_bit`

**Inputs:** `["100001", "000000", "100001", "111111"]` (4 labels, one
repeated), and the target matrix `[[1,0,0,0,0,1], [1,1,0,0,0,0]]`.

**Derivation:** label counts are `100001 → 2`, `000000 → 1`, `111111 → 1`;
fractions divide by 4, giving `0.5`, `0.25`, `0.25`. Output keys are sorted,
so the order is `000000`, `100001`, `111111`. The bit distribution is the
column sums over 2 rows: bit 0 is set twice → count 2, fraction `2/2 = 1.0`;
bit 1 once → count 1, fraction `0.5`; bit 2 never → count 0.

### `test_build_feature_importance_rows_sorts_descending`

**Input:** features `("feature_1","feature_2","feature_3")` with scores
`[0.2, 0.5, 0.3]`.

**Derivation:** pair each name with its score and sort by score descending →
`feature_2` (0.5), `feature_3` (0.3), `feature_1` (0.2).

### `test_correlation_feature_groups_collects_the_redundant_columns`

**Input:** the matrix `[[1,1,3],[2,4,1],[3,9,4],[4,16,2],[5,25,5]]`, names
`("height_mean","height_square","unrelated")`, threshold `0.7`.

**Derivation:** column 1 is column 0 squared over positive values, so the two
rank the rows identically and their Spearman correlation is exactly `1`, for a
distance of `1 - 1 = 0`. Column 2's ranks are `3,1,4,2,5` against column 0's
`1,2,3,4,5`: the Spearman correlation is `0.5`, a distance of `0.5`. The cut
sits at `1 - 0.7 = 0.3`, above the first distance and below the second, so
columns 0 and 1 merge and column 2 stays alone, giving `((0, 1), (2,))`.

### `test_grouped_permutation_importance_scores_the_group_that_carries_the_label`

**Input:** 40 rows; column 0 alternates `0.0, 1.0`, column 1 is twice column
0, column 2 is `0..39`. The model predicts `column 0 > 0.5` on all six bits,
and the targets are that same prediction, so the unpermuted micro-F1 is `1.0`.

**Derivation:** group `(2,)` holds no column the model reads, so every
permuted matrix yields the identical prediction and a micro-F1 of `1.0` in
each of the ten repeats: the mean drop is `1.0 - 1.0 = 0` and the standard
deviation over the repeats is `0`, both exactly. Group `(0, 1)` holds the
column the model reads, and a random permutation of a balanced 0/1 column
agrees with the original on about half the rows, so about a quarter of the
rows become false positives and a quarter false negatives, putting micro-F1
near `0.5` and the drop near `0.5`; the assertion asks only for a drop above
`0.2`, which the permutation distribution clears with the bound
comfortably slack.

### `test_feature_importance_rows_refuses_impurity_without_a_tree`

**Input:** a model whose `estimators_` holds one plain object, asked for
method `mdi`.

**Derivation:** the impurity path needs `feature_importances_` on every
sub-estimator; the object has none, so the call raises `ConfigError` before
reading any data, and the message names the tree-based model the measure
requires.

### `test_auto_cv_folds_follows_the_rare_class_policy`

`auto_cv_folds(labels, requested_folds, policy)` caps folds at the smallest class
count, because a stratified split cannot produce more folds than the rarest class
has members.

| Labels | requested | policy | Expected | Derivation |
| --- | --- | --- | --- | --- |
| 5x`a`, 5x`b` | 5 | `warn_reduce_cv` | `5`, no warnings | smallest class 5 >= 5 |
| 10x`a`, 3x`b` | 5 | `warn_reduce_cv` | `3` + "Reduced CV folds from 5 to 3" | capped at 3 |
| 10x`a`, 3x`b` | 5 | `warn_skip_cv` | `None` + "Skipped cross-validation" | reduction needed, so skip instead |
| 10x`a`, 1x`b` | 5 | `warn_reduce_cv` | `None` + "fewer than 2 samples" | a singleton cannot be stratified |
| 10x`a`, 1x`b` | 5 | `error` | `ValueError` | same condition, escalated |
| `[]` | 5 | `error` | `None` + "No labels available for cross-validation" | empty input handled before the policy |

### `test_read_tsv_rows_reads_records_and_rejects_a_header_only_file`

**Inputs:** a TSV with the header `class\tfeature_1\tdis1_topology` and two
data rows, `100001\t1.5\tBC` and `000000\t2.5\tAC`; and a TSV containing
only `class\tfeature_1`.

**Derivation:** the two-row file yields one dict per row keyed by the header
fields, values kept as the strings read. Feature selection on that header
with target `class` removes the target and preserves the remaining order →
`("feature_1", "dis1_topology")`. The reader requires at least one data row,
so the header-only file raises `ValueError("...no data rows")`.

## tests/orchestrator/test_orchestrator_consolidation.py

### Mapping rules

Consolidation reads four fields off each result: `triplet` normalized as
`(A, B, C)` with `A`, `B` the species-tree sisters, `classification`,
`dis1_topology` and `bootstrap_value`, and `_map_event` turns each row into
at most one event:

| Classification | `dis1_topology` | Produces |
| --- | --- | --- |
| `inflow_introgression` | `BC` | sampled edge `C → B` |
| `inflow_introgression` | `AC` | sampled edge `C → A` |
| `outflow_introgression` | `BC` | sampled edge `B → C` |
| `outflow_introgression` | `AC` | sampled edge `A → C` |
| `ghost_introgression` | `BC` | ghost target `A` |
| `ghost_introgression` | `AC` | ghost target `B` |
| anything else | any | nothing |

Each event carries the row's `bootstrap_value` as its weight. A sampled edge's
raw sum is the total weight over the rows that produced it, its supporting
count is the number of those rows, and its average is the quotient; ghost
targets are treated the same way. Non-sister counts come from every row,
whatever its classification: for a normalized triplet the pairs `{A, C}` and
`{B, C}` each gain one. A taxon's `has_sampled_introgression` flag is `1`
when it is the target of a sampled edge with non-zero average weight.

### `test_generate_introgression_maps_writes_every_sheet_from_the_results`

**Inputs:** the species tree `(((A:1,B:1):1,C:1):1,D:1);` (taxa order
`A, B, C, D`) and nine rows:

| # | Triplet | Classification | `dis1_topology` | Weight | Event |
| --- | --- | --- | --- | --- | --- |
| 1 | `(A,B,C)` | inflow | `BC` | 0.6 | edge `C → B` |
| 2 | `(A,B,C)` | inflow | `BC` | 0.2 | edge `C → B` |
| 3 | `(A,B,C)` | no_introgression | `BC` | 0.0 | nothing |
| 4 | `(A,B,C)` | inflow | `AC` | 0.5 | edge `C → A` |
| 5 | `(A,B,C)` | outflow | `BC` | 0.3 | edge `B → C` |
| 6 | `(A,B,D)` | ghost | `BC` | 0.8 | ghost target `A` |
| 7 | `(A,B,D)` | ghost | `BC` | 0.4 | ghost target `A` |
| 8 | `(D,B,C)` | ghost | `BC` | 0.4 | ghost target `D` |
| 9 | `(A,C,D)` | ghost | `AC` | 0.5 | ghost target `C` |

Then the same rows written again into an existing directory `introgression`
holding `stale.txt`, with `overwrite=False`.

**Derivation:** three distinct sampled edges (`C → B`, `C → A`, `B → C`) and
three ghost targets (`A`, `D`, `C`) give `non_ghost_edge_count == 3` and
`ghost_target_count == 3`; the four taxa of the tree give `taxa_count == 4`.
The stage writes `introgression_combined.png` beside a `consolidation_data/`
folder holding the eight TSVs.

*Sampled edges* (matrix rows are targets, columns sources): `C → B` gathers
`0.6 + 0.2 = 0.8` over 2 supporting rows → average `0.4`; row 3 contains the
same taxa but produces no edge, so it raises neither sum nor count. `C → A`
is `0.5` over 1 row → `0.5`, and `B → C` is `0.3` over 1 → `0.3`. The raw sums
total `0.8 + 0.5 + 0.3 = 1.6`, and the largest average is `0.5`, so every
value stays at or below 1 without any rescaling.

*Ghost targets*: `A` gathers `0.8 + 0.4 = 1.2` over 2 rows → `0.6`; `D` is
`0.4` over 1 → `0.4`; `C` is `0.5` over 1 → `0.5`; `B` is never a ghost target
→ raw sum 0, count 0, strength 0, but it is still listed because the sheet is
keyed by the taxa order. Flags: `A` is the target of `C → A`, `B` of `C → B`,
`C` of `B → C`, all with non-zero averages → `1`; no edge points at `D` →
`0`.

*Non-sister pairs*: rows 1-5 are `(A,B,C)` five times → `{A,C}` +5, `{B,C}`
+5. Rows 6-7 are `(A,B,D)` twice → `{A,D}` +2, `{B,D}` +2. Row 8 `(D,B,C)`
has sisters `D`, `B`, so `{C,D}` +1 and `{B,C}` +1 → 6. Row 9 `(A,C,D)` has
sisters `A`, `C`, so `{A,D}` +1 → 3 and `{C,D}` +1 → 2. Hence `{A,C}: 5`,
`{B,C}: 6`, `{A,D}: 3`, `{B,D}: 2`, `{C,D}: 2`, `{A,B}: 0`, written
symmetrically with a zero diagonal.

*Overwrite*: the output directory goes through `prepare_output_directory`,
which with `overwrite=False` leaves `introgression` and its `stale.txt` alone
and picks the smallest free suffix, `introgression_1`, so the plot path
returned points there.

### `test_generate_introgression_maps_selects_the_requested_taxa`

**Inputs (parametrized):** the balanced tree `((A:1,B:1):1,(C:1,D:1):1);` with
one ghost row on `(A,B,C)`, once unfiltered and once with
`plot_taxa=["A","B","C"]`; and `(((A:1,B:1):1,C:1):1,OG:1);` with that row plus
a ghost row on `(A,B,OG)`, with `outgroups=["OG"]`.

**Derivation:** the taxa order is the species tree's leaf order restricted to
`plot_taxa` when given, with outgroups removed and any taxa seen only in the
results appended; the matrix header and the ghost sheet are both written from
that one order. Unfiltered, all four taxa are listed. With `plot_taxa`, `D` is
dropped from the order and so from both sheets, leaving 3. With `OG` an
outgroup, it is removed even though a result names it, leaving `A, B, C`. In
each case `taxa_count` equals the length of that order and the two sheets list
the same names.

### `test_count_helpers_follow_the_classification`

**Inputs:** five rows:

| Triplet | Classification | `dis1_topology` | Event |
| --- | --- | --- | --- |
| `(A,B,C)` | inflow | `BC` | edge `C → B` |
| `(A,B,C)` | no_introgression | `BC` | nothing |
| `(A,B,D)` | inflow | `BC` | edge `D → B` |
| `(A,B,C)` | ghost | `BC` | ghost target `A` |
| `(A,C,D)` | ghost | `AC` | ghost target `C` |

For the presence flag, `taxa_order = ["A","B","C","D"]` with the averages
`{("C","A"): 0.6, ("A","B"): 0.0, ("D","C"): 0.3}`.

**Derivation:** supporting counts follow the events: `(C, B)` once, `(D, B)`
once, nothing for the `no_introgression` row; ghost counts `A` once and `C`
once. Non-sister counts come from every row: three `(A,B,C)` rows give
`{A,C}: 3` and `{B,C}: 3`; `(A,B,D)` gives `{A,D}` and `{B,D}` one each;
`(A,C,D)` gives `{A,D}` (now 2) and `{C,D}: 1`. For the flag, edges are keyed
`(source, target)`, so it looks at the second element: `A` is the target of
`(C,A)` with weight `0.6` → `1`; `B` is the target of `(A,B)` but the weight
is `0.0`, which is falsy, so an edge that carries no support does not count →
`0`; `C` is the target of `(D,C)` with `0.3` → `1`; `D` is never a target →
`0`.

## tests/orchestrator/test_orchestrator_shape.py

### `test_shape_moments_match_scipy`

- **Input**: 500 draws from `numpy.random.default_rng(4).lognormal(0, 0.7)`.
- **Expected `skew`**: `scipy.stats.skew(values)`, the biased sample third
  standardized moment `m3 / m2^1.5`. Exact equality, since the module calls the
  same function.
- **Expected `excess_kurtosis`**: `scipy.stats.kurtosis(values, fisher=True)`,
  i.e. `m4 / m2^2 - 3`, so a normal reads `0`.
- **Signs**: a lognormal with `sigma = 0.7` has population skewness
  `(e^{sigma^2} + 2) sqrt(e^{sigma^2} - 1) ~ 2.5` and positive excess kurtosis,
  so both must come out above zero.

### `test_modality_test_rejects_only_a_well_separated_mixture`

- **Inputs**: 400 draws each of `normal(0, 1)`, `lognormal(0, 0.6)`,
  `exponential(1)`, `gamma(2, 1)`, and a concatenation of 200 `normal(0, 1)`
  plus 200 `normal(4, 1)` draws, all from one `default_rng(17)` stream. The
  modality bootstrap is driven by a separate `default_rng(2)`.
- **Derivation of the expectation**: the first four are unimodal by
  construction (a gamma with shape `> 1` has its mode at `(k-1) * theta`; a
  lognormal at `e^{mu - sigma^2}`), so a test holding its nominal level must not
  reject them. The mixture's components are 4 pooled SD apart, which puts a
  genuine valley between them; measured `modes_p` values are 0.53/0.38/0.16/0.35
  for the unimodal families and 0.005 (the `1 / (n_resamples + 1)` floor at 200
  replicates) for the mixture. The assertion is the side of `alpha = 0.05` each
  falls on, not the value.
- **Why not a mode count**: the same four unimodal samples give KDE peak counts
  of 1, 4, 3 and 2 at Scott's bandwidth, so the count alone would call three of
  them multimodal.

### `test_tail_index_recovers_known_tail_shapes`

- **Inputs**: 4000 draws each: `pareto(3) + 1`, `exponential(1)`,
  `uniform(0, 1)`, from `default_rng(23)`.
- **Expected `xi`**: a Pareto with index `a` has survival `x^{-a}`, whose
  generalized-Pareto tail index is `1 / a = 1/3`. An exponential tail is the
  `xi -> 0` limit of the GPD. A uniform on `[0, 1]` has a finite upper endpoint
  reached linearly, which is the GPD with `xi = -1`.
- **Tolerance**: `abs=0.25`. The fit sees only the upper decile, so 400 points
  back each estimate; the GPD shape estimator's asymptotic standard error is
  roughly `(1 + xi) / sqrt(k)`, which is about `0.07` here, and the tolerance
  leaves room for the threshold-choice bias on top of that.

### `test_shape_is_not_described_for_small_or_flat_groups`

- **`[1.0] * 5`**: 5 observations, below `SHAPE_MIN_OBSERVATIONS = 20`, so
  every field is `None`.
- **`[2.0] * 50`**: above the size floor but with zero variance, so the KDE has
  no scale and the moments are undefined; every field is `None`.
- **25 normal draws**: clears the size floor, so the moments are reported. The
  upper decile holds `25 - ceil(0.9 * 25) = 3` points, below
  `SHAPE_MIN_TAIL_EXCEEDANCES = 10`, so `tail_xi` alone is `None`.

### `test_shape_is_measured_once_and_not_per_bootstrap_iteration`

- **Input**: one 60/40/30 lognormal triplet analyzed twice, at 5 and 60
  bootstrap iterations, with the same `triplet_seed`.
- **Derivation**: the diagnostics are computed in
  `_run_triplet_pipeline_from_observations`, which runs once for the point
  estimate; `_run_bootstrap_iterations` reaches the data through
  `_iteration_outcome`, which never calls it. The seed for the modality
  bootstrap is the fourth child of the per-triplet `SeedSequence`, and spawning
  a fourth child does not disturb the first three, so both runs draw the same
  stream. The two results must therefore agree exactly, not approximately.

### `test_shape_columns_reach_the_results_tsv_only_and_only_when_enabled`

- **Input**: one triplet with 60 concordant, 40 discordant1 and 30 discordant2
  lognormal heights, analyzed with `shape_diagnostics` off and on and written
  through both `write_pipeline_results` and `write_summary_statistics_tsv`.
- **`shape_statistics`**: populated exactly when enabled, which is what the
  "on" arm of the summary check rests on: the result carries the dict and the
  file still ignores it.
- **Results TSV, expected column set**: the product of `SHAPE_GROUP_LABELS`
  (3) and `SHAPE_FIELD_NAMES` (5), so 15 columns, present exactly when enabled.
- **`con_skew > 0`**: the heights are lognormal, which is right-skewed.
- **`con_modes_p` in `(0, 1]`**: the add-one estimator is bounded below by
  `1 / (n_resamples + 1)` and above by `1`.
- **`dis2_tail_xi` empty**: 30 observations leave `30 - 27 = 3` above the 90th
  percentile, under the 10-exceedance floor.
- **Summary TSV, expected column names**: none. No header entry ends in any
  of `SHAPE_FIELD_NAMES`, under either setting, while
  `concordant_avg_tree_height_mean` and `classification` are still present:
  the 63 descriptive columns and the label are untouched.
- **Row width equals header width, in both files and both arms**: the
  columns are appended conditionally in two places (header and row) that must
  stay in step, and dropping a column block is exactly the kind of change that
  leaves a row misaligned.
- **Why the summary file carries none**: the diagnostics are undefined below
  their observation floors: a group under `SHAPE_MIN_OBSERVATIONS` (20) has no
  modality p-value or moments, and a group whose upper decile holds under
  `SHAPE_MIN_TAIL_EXCEEDANCES` (10) points has no tail index.
  `ml_utils.rows_to_matrix` raises on any empty feature cell, so those gaps
  made whole runs unusable as training data.

## Remaining suites

`tests/test_ml_config.py`, `tests/test_ml_utils.py`,
`tests/test_ml_random_forest.py`, `tests/test_ml_multi_knn.py`, and
`tests/test_ml_hyper_tune.py` assert structural outcomes (files written,
columns present, errors raised, a shipped sample loading), and their inputs
are the fixtures described in [TESTS.md](TESTS.md). No test pins a config
default; the samples and CONFIG.md state those. The derivations worth stating
explicitly:

- **KNN neighbor capping**: `n_neighbors` cannot exceed the number of training
  samples, so the builder clamps it and the metrics report states the effective
  value.
- **`max_features` / `class_weight` forms**: the accepted set is read off
  scikit-learn's own parameter constraint: `max_features` takes `'sqrt'`,
  `'log2'`, an `int >= 1`, a `float` in `(0.0, 1.0]`, or `None`; `class_weight`
  takes `'balanced'`, `'balanced_subsample'`, a dict, a list of dicts, or
  `None`. `'auto'` is expected to fail because scikit-learn removed it in 1.3.
  `true` is expected to fail because `bool` is a subclass of `int` in Python
  and would otherwise satisfy the `int >= 1` branch as the value `1`. The
  rejected numerics `0` and `1.5` sit just outside the `int >= 1` and
  `(0.0, 1.0]` bounds respectively.
- **Confusion-matrix row normalization**: the input counts are
  `[[3, 1, 0], [0, 0, 0], [1, 1, 2]]`. Row totals are `4`, `0` and `4`. Rows 0
  and 2 divide through by `4`, giving `[0.75, 0.25, 0]` and
  `[0.25, 0.25, 0.5]`; row 1 has a zero total, so the `where=row_totals > 0`
  guard leaves it at the `[0, 0, 0]` the output buffer was initialized with
  rather than producing `nan`. That is also why the row sums come out `1, 0, 1`:
  a class with no test samples contributes nothing, and the plot masks the whole
  row.
- **W&B artifact routing**: the expected split is read straight off the
  `write_bulk_artifacts = not use_wandb` switch: the model pickle, the plaintext
  report and the search plot are written unconditionally, while the results
  JSON, ranked-candidate TSV, marginals TSV and predictions TSV sit behind the
  switch and are logged as `tables/*` payloads plus the `results_json` run
  summary in the other arm. `artifact_paths` is assembled from the same switch,
  so the three unconditional entries are exactly what it holds under
  `use_wandb: true`.
- **Hyperparameter marginals**: four candidates cross a single parameter
  `n_estimators`, with scores `(5, 0.7)`, `(5, 0.4)`, `(10, 0.9)`, `(10, 0.6)`;
  ranks 1-4 are assigned by sorting on score in the objective's direction. Under
  `max` the ordering is `0.9, 0.7, 0.6, 0.4`, so `n_estimators=10` holds ranks 1
  and 3: its `best_score` is `max(0.9, 0.6) = 0.9`, its `best_rank` is `1`, and
  it takes `value_rank` 1 ahead of `n_estimators=5` (`best_score 0.7`, `best_rank
  2`). Under `min` the ordering reverses to `0.4, 0.6, 0.7, 0.9`, so
  `n_estimators=5` wins with `best_score = min(0.7, 0.4) = 0.4` at `best_rank`
  1. Both values are numeric, so the searched range is `[5, 10]` and the winning
  value is flagged `upper_bound` under `max` (10 is the largest tried) and
  `lower_bound` under `min` (5 is the smallest); the flag is set only on the
  `value_rank` 1 row, so the runner-up's is empty.
- **Parameter influence**: the same four scores `0.9, 0.7, 0.6, 0.4` crossed
  over `n_estimators` (10, 5, 10, 5) and `max_depth` (5, 5, 3, 3). Per-value best
  scores are `n_estimators`: 10 → 0.9, 5 → 0.7, a spread of `0.2`; `max_depth`:
  5 → 0.9, 3 → 0.6, a spread of `0.3`. Mean scores give the same ordering
  (`n_estimators`: 0.75 vs 0.55 = 0.2; `max_depth`: 0.8 vs 0.5 = 0.3), so
  `max_depth` sorts first on both. Sub-`1e-12` differences are flattened to an
  exact `0.0` so a dimension whose values all score alike reports zero spread
  rather than floating-point noise, which is what lets the guidance block call it
  out as having no effect.

### `test_diagnostic_changes_what_is_measured_and_nothing_concluded`

**Inputs:** the five crafted observation sets from the cascade table (one per
outcome), measured with `family_size=5`, a 30-iteration bootstrap and seed
11, once with `diagnostic=False` and once with `diagnostic=True`, each passed
through `_apply_triplet_result_p_value_correction` under each of `bfn`,
`holm` and `fdr_bh`; then the five diagnostic results with every raw DCT
p-value replaced by 1.0, passed through the same method.

**Derivation:** a non-diagnostic run skips a test only below a gate that has
already failed on the value the correction can only raise (under `bfn` the
exactly corrected `p × 5`, under `holm`/`fdr_bh` the raw value), and the
permutation p-values are corrected inside each test rather than across
triplets, so the skipped work could neither reach its own cascade nor move
another triplet's numbers. Every field the cascade reads must therefore be
equal between the two runs, and the bootstrap votes with them: `diagnostic`
does not reach the bootstrap at all, so both runs draw and judge the same
iterations. The five cases land on `DCT`, `THT` and `PERM` gates, which is
asserted so that the comparison exercises the permutation gate rather than
only rows where neither run resamples.

What may differ is confined to what the non-diagnostic run declined to
measure. Under `holm` and `fdr_bh` the tree-height test is a rank-based family
member and is measured on every row, so the raw and corrected KS columns must
agree everywhere. Under `bfn` the `DCT`-gate row (the 10/10 split, `p = 1.0`)
is left unmeasured and carries `None` in all three KS fields; the four
measured survivors must still be corrected as `min(1, p × 5)`: the family is
the triplet count, and a family shrunk to the four measured would give
`min(1, p × 4)` and lower every survivor's corrected value. On `PERM` rows the
permutation fields agree; on the `DCT` and `THT` rows the non-diagnostic
result carries `None` for `perm_decision` and `perm_statistic` with
`perm_note == "direction_test_not_consulted"`, while the diagnostic result
still records a decision nothing read.

The diagnostic family, where every raw value exists, also shows what a
corrected column is. The pass corrects the DCT column over all five results
and the KS column over all five results, so each corrected column must equal
`_adjust_p_values` applied directly to the raw column, asserted rather than
any particular number, so the check holds for a rank-based method whose
values depend on the whole family as much as for Bonferroni. Once every count
p-value is 1.0 the tree-height column is corrected identically, since the KS
family reads nothing of the DCT column, and the cascade then calls every row
`no_introgression` at the first gate.

### `test_bootstrap_measures_the_tree_height_test_its_correction_reads`

**Inputs:** the 10/10-split observation set (`[0.1] * 20` concordant,
`[0.9] * 10` in each discordant group), `family_size=1`, a 20-iteration
bootstrap at seed 11, once per `(method, diagnostic)` pair over `bfn`/`holm`
and off/on, with `inference.run_two_sample_ks_test` replaced by a wrapper
that counts its calls and delegates; the result is then passed through
`_apply_triplet_result_p_value_correction` under the same method so that
`all_bootstrap` is available under `holm` too, whose votes are deferred.

**Derivation:** with 10 vs 10 discordant trees the chi-square p-value is 1.0,
so the point estimate's count gate fails; as a family of one, `bfn` judges it
on `min(1, 1.0 × 1)` and `holm` on the raw value, the same number. A cleared
count gate is the only way an iteration votes anything but
`no_introgression`, so the number of iterations whose resampled gate cleared
is `20 × (1 − all_bootstrap["no_introgression"])`, and it must be below 20 so
that the two counts below cannot coincide by accident.

The point estimate measures KS once unless the method is inline and the run
is not diagnostic: `bfn` off contributes 0 calls and `ks_p_value is None`;
`bfn` on and `holm` either way contribute 1 and carry a value. The bootstrap
is the same whatever `diagnostic` says: under `bfn` (`_CorrectionPolicy.inline`)
an iteration measures KS only when its own count gate cleared, so the total
is the point estimate's contribution plus the cleared count; under `holm`
every iteration's value is a rank-based family member and is measured, so the
total is the point estimate's contribution plus 20.

### `test_diagnostic_bootstrap_records_every_test_without_moving_a_vote`

**Inputs:** two observation sets: `(_HIGH, _LOW[:10], [0.5] * 10)`, ten
concordant heights spread over 0.85-0.94 against 10 and 10 discordant trees,
and the outflow set `(_HIGH, _LOW, [0.1] * 2)`, each measured with
`family_size=1` under `bfn`, a 30-iteration bootstrap at seed 11, with
`bootstrap_options.diagnostic` off, on, and on with `summary_only`, then
decided as a family of one.

**Derivation:** the first set's 10/10 split gives `chisquare([10, 10])` a
p-value of 1.0, and a resample clears the count gate only when its split
reaches 15/5 or wider (`p ≈ 0.025`), so most iterations fail at the first
gate and their direction tests are ones only the record asked for. The second
set clears every gate in most resamples, so most direction tests are the
vote's own. In both, every group keeps its spread under resampling (the
heights are distinct), so no guard fires and every iteration records a
statistic and a p-value pair.

The vote's direction tests draw from the bootstrap permutation stream; the
record's extra tests draw from a fifth child of the triplet's seed sequence,
so the vote stream is consumed identically with the record on or off and
`all_bootstrap` and the studentized interval must be equal. As a family of
one under `bfn`, `_gate_p_value` is `min(1, p × 1) = p`, so replaying
`_classify_introgression(dct_p <= 0.05, ks_p <= 0.05, decision)` over the
recorded iterations is the cascade the vote applied, and the tally it rebuilds
must equal `all_bootstrap` exactly, which also proves the recorded decision
in a gate-cleared iteration is the one the vote read. With `summary_only`,
`_decision_summary` counts the labels, so `greater + less + inconclusive` is
30 and `greater` matches the list's count; `_numeric_summary` on the
`p_greater` list reports `count` 30 and `non_null_count` 30.

### `test_shifted_statistics_match_an_explicit_shift`

**Inputs:** `x ~ Gamma(2.0, 1.0)` and `y ~ Gamma(2.5, 1.2)` at five `(nx, ny)`
splits (`(10, 30)`, `(30, 10)`, `(7, 7)`, `(4, 25)`, `(2, 60)`), pooled,
mean-centered, and drawn 2000 times under seed 31 with `shifts=(0.75, -0.75)`.

**Derivation:** adding a constant `c` to the first `nx` entries of the pooled
vector is what makes the samples exchangeable under
`H0: mean(x) - mean(y) == -c`. For a permutation gathering subset `S`, with
`a = |S ∩ X|` the count of sampled entries from the x block and `b` the sum of
their values:

    sum_w(S)   = sum_v(S)   + c * a
    sumsq_w(S) = sumsq_v(S) + 2c * b + c^2 * a

and the pooled totals shift by `c * nx` and `2c * sum(v over X) + c^2 * nx`. The
studentized statistic is a function of those power sums alone, so the fused row
is algebraically identical to re-running the draw on shifted data, not an
approximation. Requesting the same count under the same seed produces the same
`indices`, since the shifts consume no randomness, so the comparison is
element-wise. Both group orderings and lopsided splits are covered because the
kernel samples whichever group is smaller and recovers the other by subtraction,
which is where a sign or role error would surface. The algebra is exact, so
the `1e-9` tolerance is slack for floating-point reassociation only.

## tests/orchestrator/test_orchestrator_triplet_geometry.py

All derivations use one reference tree:

```
(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);
```

Nodes are numbered in pre-order, which is the numbering
`build_triplet_geometry` assigns: `root=0, X=1, W=2, P=3, Q=4, R=5, Y=6, S=7,
T=8`. `X` is the clade `((P,Q),R)`, `W` is `(P,Q)`, and `Y` is `(S,T)`. The
cache therefore holds `parent = [-1,0,1,2,2,1,0,6,6]` and
`edge_len = [0.0,1.0,2.0,1.0,1.0,3.0,2.0,2.0,2.0]`, and the pairwise LCA table
records `PQ->2`, `PR=QR->1`, `ST->6`, and every P/Q/R-to-S/T pair `->0`.

### `test_geometry_matches_each_reference`

**Inputs (parametrized over eleven trees, each read under both references
and all six strategies):** the reference tree with a nested sister pair (`P,Q,R`), a pair spanning the root (`P,R,S`), a
triplet drawn from both sides (`P,S,T`) and one whose odd taxon is listed first
(`S,P,Q`); a five-taxon ladder read at two depths; a tree carrying three taxa
that get pruned away; `(((A,B),C),(D,E));` with no branch lengths at all; one
with a length missing from a single edge; one whose internal branch is exactly
`0.0`; and one at the edge of double precision
(`1e-12` tips under a `1e-13` internal branch against a `1.0000000000001`
sister).

**Derivation:** there is no closed form to compare against here: the expected
value *is* what an independent implementation produces, which is the point,
and there are two of them, both in `tests/orchestrator/tree_references.py`.
`dendropy_observation` copies the triplet's subtree out with
`extract_triplet_subtree` and measures it with `observation_from_subtree`.
`biopython_observation` never prunes: on the
`Bio.Phylo` tree it takes the common ancestor of all three leaves and of each
pair, calls the one pair whose ancestor is a different node the sisters (all
three coinciding is a polytomy, so `None`), and reads every distance as a
`Bio.Phylo` path sum from the three-way ancestor (the depth of each leaf, the
depth of the sisters' ancestor as the internal branch, and the leaf-to-leaf
distance of the sisters), from which `AVG`/`A`/`B`/`C`/`SIS`/`INT` follow by
definition. `_geometry_observation` builds the cache and reads the same triplet
out of it.
The topology must be equal exactly, because all three derive it from discrete
structure rather than arithmetic: extraction from the copied subtree's sister
clade, BioPython from which pair's common ancestor is not the three-way one,
the cache from which two of the three pairwise LCAs coincide. Heights and
summary metrics compare at `rel=1e-12`, comfortably wider than the only
disagreement three correct implementations can have: floating-point rounding
from summing the same edges in a different order. The cases are chosen for
what they break rather than for coverage:
the zero-length internal branch would be read as a polytomy by any
depth-comparing rule, the missing lengths must count as `0.0` rather than
propagate `None`, the pruned taxa must not enter any path sum, and the
`1e-12`/`1.0000000000001` case puts the two paths' summation orders as far apart
as the fixture set can.

### `test_geometry_matches_each_reference_across_a_nine_taxon_tree`

**Inputs:** the nine-taxon tree

```
((((T1:0.11,T2:0.19):0.23,(T3:0.07,T4:0.31):0.0):0.17,((T5:0.29,T6:0.13):0.41,T7:0.53):0.09):0.37,(T8:0.61,OUT:0.71):0.43);
```

and all `C(9,3) = 84` triplets, read out of one cache built over the whole
tree, under each of the six tree-height strategies against each reference.

**Derivation:** the expected values are whatever the reference produces (
`extract_triplet_subtree` + `observation_from_subtree`, or
`biopython_observation`), so the test is a differential one: it asserts the
implementations agree rather than restating the arithmetic. The tree is shaped
so the sweep covers the cases that distinguish them: `(T3,T4)` sit above a
zero-length internal branch, `T1..T4` and `T5..T7` sit in sibling clades so
many triplets have their sister pair on one side and the odd taxon on the
other, `T7` hangs off a ladder at a different depth from its clade-mates, and
`T8`/`OUT` sit across the root so triplets drawn from them resolve at the seed
node. Every triplet of a nine-taxon rooted binary tree is resolved, so
`triplet_resolution` must call all 84 resolved and every reference and the
cache must yield an observation for each, asserted per triplet, so none is
silently skipped.

One cache answers every triplet, which is the property the design rests on:
the cache is built per tree, not per triplet, costing `O(k^2)` in the pairwise
LCA table while the tree supplies `C(k,3)` triplets. A cache that answered
only the triplets it was built from, or went stale after the first read,
fails here.

The tolerance is `rel=1e-12`, matching the parity tolerance used elsewhere in
the suite. Every path sums the same edge lengths, in an order that can differ
between implementations, so floating-point rounding is the only disagreement
possible and the tolerance sits well above it.

### `test_geometry_matches_hand_derived_values`

**Inputs (parametrized):** the reference tree with triplets `(P,Q,R)` and
`(P,R,S)` under `AVG`; `((A:1.0,B:1.0):0.0,C:1.0);` with `(A,B,C)` under
`INT`; and `(((A,B),C),(D,E));` with `(A,B,C)` under `AVG`; summary metrics
collected, and each read once more without them.

**Derivation, `(P,Q,R)`:** `LCA(P,Q)=2` while `LCA(P,R)=LCA(Q,R)=1`. The two
that agree name the triplet root, so `r=1` and the sister LCA is `s=2`, giving
topology `((A,B),C)`. Walking edges up the parent chain, `internal_branch =
edge_len[2] = 2.0`; `P` and `Q` are each `edge_len[3] = 1.0` below `W`, so each
sits `1.0 + 2.0 = 3.0` below the subtree root, and `R` is `edge_len[5] = 3.0`
below `r` directly. Hence `avg = (3.0+3.0+3.0)/3 = 3.0` and `sister_distance =
3.0 + 3.0 - 2(2.0) = 2.0`, which is the real P-to-Q path `1.0 + 1.0`. Extracting
the subtree gives `((P:1,Q:1):2,R:3);` and the same four numbers.

**Derivation, `(P,R,S)`:** `LCA(P,R)=1` while `LCA(P,S)=LCA(R,S)=0`, so `r=0`
(the root) and `s=1`, again `((A,B),C)` but with `P` and `R` as the sisters.
`internal_branch = edge_len[1] = 1.0`. `P` is `edge_len[3] + edge_len[2] = 3.0`
below `X`, so `3.0 + 1.0 = 4.0` below the root; `R` is `3.0 + 1.0 = 4.0`; `S` is
`edge_len[7] + edge_len[6] = 4.0`. So `avg = 4.0` and `sister_distance =
4.0 + 4.0 - 2(1.0) = 6.0`, the real P-to-R path `1.0 + 2.0 + 3.0`.

**Derivation, zero-length internal branch:** the internal node sits at the
same depth as the root, so any rule that picked the sister pair by comparing
LCA *depths* would see a three-way tie and drop the observation. The cached
path compares LCA node indices and resolves `((A,B),C)`. `INT` is then the
zero-length edge itself, `0.0`; `A` and `B` each sit `1.0 + 0.0 = 1.0` below
the root and `C` `1.0`, so `avg_tree_height = 1.0` and
`sister_distance = 1.0 + 1.0 - 2(0.0) = 2.0`.

**Derivation, no lengths:** `build_triplet_geometry` stores `0.0` wherever the
Newick omits a length, matching `_distance_to_root`, which skips a `None`
edge. Every walk therefore sums zeros and all four derived values are `0.0`.

**Derivation, metrics slot:** the third element of an observation carries the
per-tree summary metrics, and a run only needs them under
`generate_summary_stats`. Read without them the slot must be `None` rather
than an empty dict, since the per-triplet measurement tests that slot for
`None` to decide whether to aggregate.

### `test_geometry_skips_exactly_what_each_reference_skips`

**Inputs (parametrized):** `(A:1.0,B:1.0,C:1.0);`,
`((A:1.0,B:1.0):1.0,D:2.0);` and `((A:1.0,B:1.0):0.0,C:1.0);`, each with
triplet `(A,B,C)`.

**Derivation:** in the polytomy all three pairwise LCAs are the root, so
`triplet_resolution` reports it unresolved, the cached path sees three equal
ids and returns `None`; the extraction path reaches `find_sister_pair`, finds
no pair whose MRCA differs from the root, and raises `ValueError`, which
`observation_from_subtree` converts to `None`; the BioPython path finds no
pair whose common ancestor differs from the three-way one and returns `None`
itself. In the second tree `C` is absent, so `leaf_node[C] = -1` in the cache
(reported as a missing taxon), `set(triplet).issubset(tree_taxa)` fails in
extraction, and the label is missing from BioPython's terminals. The third
tree is resolved on every path: the zero-length branch is not a polytomy,
because all three compare node identity rather than depth, so the diagnostic
reports it resolved and every path yields an observation. `triplet_resolution`
restates `geometry_observation`'s guards without sharing code with it, so
agreeing on all three rows is what stops the preflight's reasons drifting
from the hot path's decisions.

A duplicated taxon label is deliberately not covered: reading a tree file
refuses one (see `test_read_tree_file_rejects_a_repeated_leaf_label`), so no
cached geometry can be built from such a tree.

## tests/orchestrator/test_orchestrator_rename_map.py

### `test_rename_map_reads_a_tsv_or_a_yaml_mapping`

**Inputs (parametrized):** the same two pairs (`T1 -> Homo sapiens`,
`T2 -> Pan troglodytes`) written once as a TSV carrying a `#` comment line
and a blank line, and once as a YAML mapping.

**Derivation:** the loader picks its parser from the file extension, so both
files must yield the identical dict. The comment and blank lines are dropped
before parsing, which is why the TSV's four lines produce two entries.

### `test_rename_map_rejects_malformed_or_missing_files`

**Inputs (parametrized):** `T1\tA\textra` (three columns), `T1` (one
column), `T1\tA` twice with different values, `T1\tA` and `T2\tA`, the same
shared name as YAML, a YAML list, two display names holding a delimiter (
`"Homo\tsapiens"` in YAML (whose double-quoted scalar turns `\t` into a tab)
and `Homo, sapiens` in TSV), and a path that does not exist.

**Derivation:** a rename map is a bijection from tree label to display name.
Three columns and one column both fail the two-column requirement. A repeated
label is ambiguous about which name wins. Two labels sharing a name is the case
worth singling out, and is checked in both file formats: it would rename two
distinct taxa to the same string, so every output would merge them: one row
for both in the consolidation matrices, an `A=X;B=X` triplet in the results,
with nothing pointing back at the map. A YAML list carries no keys, so it
cannot be a mapping. A tab or line break ends a TSV cell and a comma separates
the `triplet` column's names (as `;` and `=` structure `abc_mapping`, refused
by the same lookup), so a name holding one would corrupt the column it lands
in; the tab has to come through YAML because a tab inside a TSV value is read
as a third column. Each of these raises `ValueError` naming the problem. The
missing file raises `FileNotFoundError` instead, because a mistyped path is a
different mistake from a malformed map, and the message names the path so it
can be corrected without opening anything.

### `test_renaming_labels_maps_leaves_and_quotes_as_needed`

**Inputs (parametrized):** the map `{"T1": "Alpha", "T2": "Beta sp."}` over
four Newick strings as the run writes them, and the first string under an
empty map:

| Newick | Map | Expected |
| --- | --- | --- |
| `((T1:0.1,T2:0.2):0.3,T3:0.4);` | full | `((Alpha:0.1,'Beta sp.':0.2):0.3,T3:0.4);` |
| `(T3:0.4,(T1:0.1,T2:0.2):0.3):0.5;` | full | `(T3:0.4,(Alpha:0.1,'Beta sp.':0.2):0.3):0.5;` |
| `((T1:0.1,T10:0.2):0.3,XT1:0.4);` | full | `((Alpha:0.1,T10:0.2):0.3,XT1:0.4);` |
| `(('O''Brien':0.1,T2:0.2):0.3,T3:0.4);` | full | `(('O''Brien':0.1,'Beta sp.':0.2):0.3,T3:0.4);` |
| `((T1:0.1,T2:0.2):0.3,T3:0.4);` | `{}` | `((T1:0.1,T2:0.2):0.3,T3:0.4);` |

**Derivation:** the `species_tree` column is renamed as text, so the renamer
must find exactly the leaf labels: the token after `(` or `,` and before `:`,
`,` or `)`. Branch lengths follow `:` and so are never candidates (the `0.1`
and `0.3` survive every row), and a token is matched whole, so `T10` and `XT1`,
which merely contain the key `T1`, are left alone (row 3). `T2 -> Beta sp.`
puts a space and a dot into the label, which a bare token cannot hold, so it
comes out `'Beta sp.'`; `T1 -> Alpha` needs no quotes and gets none. Row 2 puts
the odd taxon first and adds a root edge, covering both child orders the shape
writer produces. Row 4 starts from an already-quoted label: it is unquoted
before the lookup (`O'Brien` is not in the map) and written quoted again, so a
label the writer had to quote round-trips. Row 5 is the no-rename case every
run without a `species_rename_map` goes through: the string must come back
unchanged. Each output is then parsed with DendroPy and must give, leaf by
leaf, the mapped name and the input's edge length (the check that the string
is still valid Newick and names the right taxa, independent of the exact
spelling asserted above), and the plain label helper, which renames the
results' `triplet` tuples and the taxon lists consolidation is handed, must
map the input's leaf list to those same names, a label absent from the map
passing through unchanged.

### `test_species_filter_runs_every_triplet_among_the_named_species`

**Inputs:** the shared 4-triplet orchestrator fixture (ingroup `A`, `B`, `C`,
`D` after pruning the outgroup `OUT`), serial, with a species filter file per
row: one whose lines are `A, B`, `D`, an empty line, `OUT`, `NOPE`, `B`; and
one reading `A,B,NOPE`.

**Derivation:** the first file yields `A`, `B`, `D`, `OUT`, `NOPE` (the second
`B` is a repeat). `OUT` is an outgroup taxon and `NOPE` is not in the tree, so
neither is an ingroup taxon and both are skipped; `A`, `B`, `D` remain, and
`C(3, 3) = 1` triplet, `(A, B, D)`, is generated in place of the fixture's
four. Its counts are `12/0/0` with species subtree `((A,B),D);`, exactly as in
`test_run_orchestrator_matches_derived_expectation`: the filter changes which
triplets are set up, not the extraction or inference that follows. The second
file leaves only `A` and `B`, which cannot form a triplet, so the runner
raises `InputError` ("at least 3 are needed") before any gene tree is
measured.

### `test_species_rename_map_reaches_every_output`

**Inputs:** the shared 4-triplet orchestrator fixture (taxa `A`, `B`, `C`, `D`,
outgroup `OUT`) with a TSV mapping `A -> Homo sapiens` and `B -> Pan sp.`,
consolidation enabled.

**Derivation:** the run works in the trees' own labels and
`runner._rename_result_taxa` rebuilds the results under the display names after
the decision pass, before anything is written. The display names hold a space
and a dot on purpose: a bare Newick label cannot, so had they been written into
the processed trees and reread, Bio.Phylo would read `sapiens` and `sp.` and
DendroPy would refuse the file: the failure this design avoids. Downstream of
the rename every output sees only display names, which is why the assertions
can span outputs written by unrelated code paths: the triplet tuples
(`{Homo sapiens, Pan sp., C, D}`; `C` and `D` are absent from the map and so
unchanged); the `(A,B,C)` species subtree, topology-only in this fixture, which
becomes `(('Homo sapiens','Pan sp.'),C);` with the two names quoted and `C`
bare, and every `species_tree` column parsing back to exactly its triplet's
names; the `triplet` and `abc_mapping` columns of the results TSV; and the
consolidation artifacts. The two `processed_*.tree` files are the outputs the
run reads back, so their leaves must be drawn from `{A, B, C, D, OUT}` and
nothing else. The consolidation check reads `introgression_taxa_order.tsv`
because that file holds the taxon ordering used to label the heatmap axes and
the bar chart, so it standing in display names is the evidence the plot labels
do too; consolidation gets there by reading the processed species tree in tree
labels and mapping them in memory through the same map.
