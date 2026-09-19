# Test Input/Output Derivations

This document is the granular companion to [TESTS.md](TESTS.md). For each test
it records **exactly what the inputs are** (literal values, not just fixture
names) and **how the expected output is derived** — the arithmetic, the
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
the ingroup clade's edge and leaves `OUT` at length 0. For tree 0 the ingroup
edge `0.10` plus `OUT`'s `0.50` gives `0.6`:

```
((((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6,OUT:0);
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
species subtree `((A:1.0,B:1.0):1.0,C:2.0);` (concordant topology `((A,B),C)`).
Each subtree's geometry is tabulated in `_LEAF_GEOMETRY` as
`(topology, dist_A, dist_B, dist_C, internal_branch)`, where a leaf's
root-to-tip distance is its own edge plus the internal branch if it is in the
sister pair, or just its own edge otherwise.

Worked example, index 0 — `((A:0.10,B:0.10):0.10,C:0.30);`:
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

Worked example, index 5 — `((A:0.30,C:0.30):0.10,B:0.70);` (sisters A, C):
`AVG = (0.40 + 0.70 + 0.40)/3 = 0.50`; `A = 0.40`; `B = 0.70`; `C = 0.40`;
`SIS = 0.40 + 0.40 - 2(0.10) = 0.60`; `INT = 0.10`.

### Summary statistics

- `mean` — arithmetic mean.
- `median` — middle value (mean of the two middle values for an even count).
- `mode` — used only by the 63-column `summary_statistics.tsv` output, not by
  the decision logic. Values are **binned to 3 decimals**, the most frequent bin
  wins, and ties resolve to the **largest** value. With all-distinct values every
  bin has count 1, so the mode is the maximum. For the AVG concordant sample
  `[0.233, 0.243, 0.317, 0.297, 0.250]` the mode is therefore `0.317`.

## tests/orchestrator/test_orchestrator_inference.py

### `test_inference_matches_derived_expectation`

**Inputs:** the 10 gene subtrees above serialized with the `AVG` strategy,
`alpha_dct = alpha_ks = 0.05`, one of the 2 discordant tests, and one of the 2
pipeline modes. Bootstrap runs with 40 iterations at seed `20240724`. The
result is measured by `analyze_triplet_from_observations` and decided by
`_apply_triplet_result_p_value_correction` under `no` as a family of one.

**Expected-output derivation**, performed in `_expected_result`:

1. Group the hand-derived heights by topology using `_LEAF_GEOMETRY`.
2. Rank the discordants by count: `n_dis1 = 3` (BC), `n_dis2 = 2` (AC),
   `n_con = 5`.
3. `most_frequent_matches_concordant` = `5 >= 3 and 5 >= 2` = `True`.
4. **DCT** — chi-square: `scipy.stats.chisquare([3, 2])` with implied expected
   frequencies `[2.5, 2.5]`, giving
   `(3-2.5)^2/2.5 + (2-2.5)^2/2.5 = 0.1 + 0.1 = 0.2`, and with df 1,
   `p ~ 0.6547`. z-test: `proportions_ztest(count=[3,2], nobs=[5,5])`, giving
   `z ~ 0.6325`, `p ~ 0.5271`. Both p-values exceed 0.05, so
   `dct_significant is False` — gate 1 settles the call.
5. **KS** — under `detailed`, measured regardless:
   `scipy.stats.ks_2samp(concordant_heights, dis1_heights)` on
   `[0.2333, 0.2433, 0.3167, 0.2967, 0.2500]` and `[0.3667, 0.4067, 0.5167]`;
   they are completely separated, so `D = 1.0` and `p ~ 0.0357`, giving
   `ks_significant is True`. The cascade never reads it, which is what the row
   demonstrates. Under `efficient` the correction is `no`, an inline method
   whose family of one is fixed, so the test below the failed count gate is
   never measured: `ks_statistic`, `ks_p_value` and `ks_significant` are all
   `None`.
6. **Direction** — under `detailed`, measured regardless: the pooled sample is
   `n_con + n_dis1 = 5 + 3 = 8` observations, admitting only
   `C(8, 3) = 56` distinct group assignments. That is far below the 2500
   `min_resamples` floor, so the `insufficient_permutation_support` guard fires
   and `perm_decision` is `inconclusive` with zero resamples. Under `efficient`
   the failed count gate settles the cascade first, so the test is skipped and
   `perm_decision` is `None`.
7. **Classification** — gate 1 fails, so `no_introgression` in every case.

Also asserted: `triplet == ("A","B","C")`, `species_tree == "((A,B),C);"` (the
species subtree serialized topology-only), and that the bootstrap class
fractions sum to 1.

### `test_observation_heights_match_derived_geometry`

**Inputs:** the same 10 subtrees, one strategy per parametrization.

**Derivation:** for each observation, the topology must equal the tabulated
topology and H(T) must equal `_expected_height(entry, strategy)` computed from
the four geometry primitives — i.e. the strategy formulas above applied to the
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
under each pipeline mode and decided under `no` correction as a family of one,
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
  significant and under `detailed` the row asserts `ks_significant is True` —
  which is why this row also shows the DCT gate stopping the cascade before a
  later gate can be consulted; the direction test runs as well and reports a
  decision the classification never reads. Under `efficient` neither is
  measured: `no` is an inline correction, so the failed count gate is final
  and the row asserts `ks_p_value`, `ks_significant` and `perm_decision` all
  `None` with `perm_note == "direction_test_not_consulted"`.
- **`inflow`.** `chisquare([30, 2])` has expected `[16, 16]`, so the statistic is
  `(30-16)^2/16 + (2-16)^2/16 = 24.5`, and with df 1 `p ~ 7.4e-07 < 0.05`: gate 1
  passes. Every concordant and dis1 height is 0.5, so the two empirical CDFs
  coincide: `D = 0.0`, `p = 1.0`, not significant → gate 2 returns
  `inflow_introgression`.
- **`outflow`.** The DCT is the same significant 30-vs-2 split.
  `max(_LOW) = 0.34 < min(_HIGH) = 0.85`, so the CDFs separate completely:
  `D = 1.0` and gate 2 passes. The spread matters for gate 3 — constant samples
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
reached the same classification by a different route would fail. Under
`efficient`, `perm_decision` is asserted non-`None` only on the three `PERM`
rows and `None` on the `inflow` row too, whose failed tree-height gate settles
the cascade before the direction test.

### `test_permutation_guards_surface_on_the_triplet_result`

**Inputs and derivation**, one row per guard, measured under
`pipeline_mode="detailed"` — the first row's 2-vs-2 discordant split gives
`chisquare([2, 2])` a statistic of `0.0` and `p = 1.0`, so its count gate
fails and the efficient mode would skip the direction test before the guard
could fire:

| con | dis1 | Guard | Why |
| --- | --- | --- | --- |
| `[0.9, 0.8, 0.7, 0.6]` | `[0.1, 0.2]` | `insufficient_permutation_support` | `C(6, 2) = 15 < 2500`, so the permutation distribution cannot resolve `alpha_perm`. |
| `[0.9] * 10` | `[0.1] * 30` | `degenerate_observed_scale` | Both groups are internally constant, so the standard error is floating-point noise near `1e-17`. |

Each returns `perm_n_resamples == 0` and `perm_decision == "inconclusive"` — a
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

`dis1_topology` names whichever discordant topology is more frequent — `AC` on
the first row, `BC` on the second — and `discordant1_*` must describe that same
group of gene trees, with `discordant2_*` describing the other. Counts are
`(n_dis1, n_dis2) == (9, 3)` either way, so the roles follow the counts rather
than the topology label.

### `test_classify_introgression_truth_table`

**Inputs:** `_classify_introgression(dct_significant, ks_significant,
direction)` called directly with 9 explicit rows; it returns the
`(classification, decision_gate)` pair asserted below.

**Derivation:** straight from the decision definition. The gate is the name of
the test whose branch returned, so it is fixed by `dct_sig` and `ks_sig` alone —
`direction` only ever selects among the three `PERM` classifications.

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

### `test_inline_and_deferred_correction_agree_on_a_single_triplet`

**Inputs:** one triplet with 25 concordant heights at 0.9, 22 discordant1 at 0.35
plus 3 at 0.2, and 4 discordant2 at 0.3; `family_size=1`, 40 bootstrap
iterations, `triplet_seed=3`; each of `no`, `bfn`, `holm`, `fdr_bh`, `fdr_by`.

**Derivation:** with one test in the family every correction is the identity —
Bonferroni multiplies by 1, and the rank-based methods adjust the single value
`p_(1)` by `(n - 1 + 1)/1 = 1` (Holm), `n/1 = 1` (BH), or `1 x sum(1/i) = 1` for
`n = 1` (BY). So all five must reach the same per-iteration verdicts. `no` and
`bfn` reach them inline during the stream while the other three park raw
p-values and are corrected afterwards, so equality is a statement about the two
code paths rather than about arithmetic. The resample stream is seeded
identically and is independent of the permutation generator, so the iterations
themselves are the same draws in every case.

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
fraction falls below 1.0 — which is what makes the corrected agreement a real
check rather than a tautology.

### `test_deferred_bootstrap_record_is_cleared_after_correction`

**Inputs:** the same observation set as the agreement test, `holm`, 20
iterations.

**Derivation:** `holm` is rank-based, so `analyze_triplet_from_observations`
cannot classify the iterations and must park them: `all_bootstrap is None` and
`bootstrap_deferred.dct_p_values` has one entry per iteration, 20.
`_apply_triplet_result_p_value_correction` then corrects iteration `i` across
the (here single-member) family, tallies the votes, and clears the record, so
`bootstrap_deferred is None` and the 20 votes normalize to fractions summing to
1.

### `test_vectorized_bootstrap_codes_match_classify_introgression`

**Inputs:** all 2 x 2 x 4 = 16 combinations of `dct_significant`,
`ks_significant`, and direction in `{greater, less, equivalent, inconclusive}`.

**Derivation:** `_classification_codes` writes integer codes over whole arrays
while `_classify_introgression` returns a label for one row; the code indexes
`_BOOTSTRAP_CLASSES`, so `_BOOTSTRAP_CLASSES[code]` must equal the label for
every combination. The expected value is whatever the scalar function returns —
the point is agreement between the two implementations, not a third derivation.

### `test_every_supported_correction_is_monotone`

**Inputs:** the family `[0.001] * 8 + [0.4, 0.9]` under each method in
`P_VALUE_CORRECTION_CHOICES`.

**Derivation:** each of `no`, `bfn`, `holm`, `fdr_bh`, and `fdr_by` applies a
multiplier of at least 1 to the `j`-th smallest of `n` — `n` for Bonferroni,
`n - j + 1` for Holm, `n/j` for BH, and `(n/j) x sum(1/i)` for BY — so no
adjusted value can fall below its raw one.

The family is deliberately the hostile case for that property: a method that
estimates the number of true nulls `n0` and substitutes it for `n` would reject
8 at the first stage, giving `n0 = 2` and a multiplier of `2/j` — below 1 for
every `j > 2`, so the later strong p-values would land beneath their raw ones.
Parametrizing over the choice list rather than a fixed set of names means such a
method fails here the moment it is added.

The `1e-12` slack absorbs floating-point rounding in the running max/min sweeps.

### `test_inline_bonferroni_matches_the_family_correction`

**Inputs:** `p = 0.004` with family sizes 1, 7, and 250, padded to that length
with `0.5` entries.

**Derivation:** Bonferroni is `min(1, n x p)` and depends on the family only
through its size, so the inline form must equal the full pass exactly:
`0.004`, `0.028`, and `1.0` respectively (`250 x 0.004 = 1.0`).

### `test_inline_correction_rejects_a_rank_based_method`

**Inputs:** `holm` passed to `_adjust_p_value_inline`.

**Derivation:** Holm's multiplier depends on a p-value's rank within its family,
which a single value does not determine, so the call must raise rather than
silently pick a wrong multiplier.

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
method is not `no` — 36 columns with a correction against 32 without. The
significance flags are unconditional because the cascade reads them whatever the
method is. The field-count assertion catches the failure mode this change could
introduce: the header and the row build their conditional sections separately,
so a mismatch would silently shift every later column by one.

### `test_adjust_p_values_matches_statsmodels`

**Inputs:** `[0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6]`
and one of the six methods, `alpha = 0.05`.

**Derivation:** `no` must return the input list unchanged. Every other method
must equal `statsmodels.stats.multitest.multipletests(p_values, alpha=0.05,
method=m)[1]` where `m` maps `bfn → bonferroni`, `holm → holm`,
`fdr_bh → fdr_bh`, `fdr_by → fdr_by`.

### `test_adjust_p_values_bonferroni_by_definition`

**Inputs:** `[0.01, 0.2, 0.5]`, method `bfn`.

**Derivation:** Bonferroni multiplies each p-value by the number of tests (3) and
clamps at 1: `0.01x3 = 0.03`, `0.2x3 = 0.6`, `0.5x3 = 1.5 → 1.0`.

### Degenerate-input tests

- `test_discordant_count_test_with_no_discordant_observations` — `n_dis1 +
  n_dis2 == 0` means there is nothing to test, so the function returns the
  neutral `(0.0, 1.0)` before touching SciPy.
- `test_ks_test_with_an_empty_sample` — an empty sample on either side makes the
  KS statistic undefined, so the function returns the neutral `(0.0, 1.0)`.
- `test_adjust_p_values_rejects_unknown_method` and
  `test_discordant_count_test_rejects_unknown_method` — an unsupported name is
  outside the choice tuple, so a `ValueError` naming the valid options is
  raised.

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

### `test_observed_statistic_matches_scipy`

**Inputs:** five seeded sample pairs, 500 resamples.

**Derivation:** the observed statistic is a deterministic function of the inputs,
not of the resampling, so it must equal SciPy's `result.statistic` and the
in-test reference `(mean(x) - mean(y)) / sqrt(var(x)/nx + var(y)/ny)` exactly
(to floating-point tolerance). Any disagreement would mean the two are not
testing the same quantity.

### `test_p_values_match_scipy_within_monte_carlo_error`

**Inputs:** eight seeded sample pairs, `min_resamples = max_resamples = 4000`
(pinning the adaptive stopping off), `correction="no"` so the comparison is
against SciPy's uncorrected p-values.

**Derivation:** both implementations estimate the same quantity from independent
resampling streams, so exact equality is not expected. Each estimate has
binomial standard error `sqrt(p(1-p)/n)`, and the difference of two independent
estimates has `sqrt(2)` times that. The tolerance is 5 such standard errors,
which at `p = 0.5` and `n = 4000` is about `0.056` and at `p = 0.01` about
`0.011`. A systematic error in the sampler or the counting would exceed this;
ordinary Monte Carlo scatter will not.

### `test_decision_matches_scipy_directional_verdict`

**Inputs:** three seeded pairs, `x ~ N(1.0, 0.4^2)` at n=60 against
`y ~ N(0.4, 0.6^2)` at n=45, 4000 resamples, `correction="bfn"`.

**Derivation:** with Bonferroni over a family of two, the corrected p-value is
`2p`, so `2p <= alpha` is the same condition as `p <= alpha/2`. The expected
verdict is therefore computed from SciPy's raw one-tailed p-values at
`alpha/2 = 0.025`, and must match `result.decision`.

### `test_permutation_statistics_match_exhaustive_enumeration`

**Inputs:** `x = [0.11, 0.24, 0.37, 0.52]` (nx=4),
`y = [0.63, 0.71, 0.88, 0.95, 1.10]` (ny=5), 5000 sampler draws at seed 7.

**Derivation:** with 9 pooled values split 4/5 there are exactly
`C(9, 4) = 126` distinct group assignments. Enumerating all of them with
`itertools.combinations` and evaluating each through the scalar reference
`_studentized_mean_diff` gives the exact support of the permutation
distribution, rounded to 9 decimals. The vectorized sampler — which never
materializes the 5-element group, recovering its sum and sum-of-squares by
subtracting from the pooled totals — must emit only values from that set
(soundness). Over 5000 draws from 126 equally likely assignments the chance of
missing any one is about `126 x (125/126)^5000 ~ 1e-15`, so all 126 must also
appear (completeness). Together these pin the optimization to ground truth.

### `test_random_inputs_preserve_test_invariants`

**Inputs:** twelve seeded pairs, `min_resamples=600`, `max_resamples=3000`.

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

### `test_equal_samples_give_a_zero_statistic_and_no_direction`

**Inputs:** the same 8 values `[0.10, 0.22, 0.31, 0.44, 0.55, 0.61, 0.78, 0.83]`
as both samples.

**Derivation:** identical samples have identical means, so the numerator is
exactly 0 and the statistic is 0. Half the permutation distribution lies on
either side of 0, so both one-tailed p-values are near 0.5 and neither clears
`alpha`. The equivalence step then decides between `equivalent` and
`inconclusive`; at 8 observations per group it has too little power to rule out
a medium effect, so the assertion accepts either.

### `test_equivalence_needs_enough_data_to_conclude`

**Inputs:** `x` and `y` both drawn from `N(1.0, 0.2^2)` at n=8 and at n=400 per
group, 1000 resamples, fixed seeds.

**Derivation:** both samples come from one distribution, so neither directional
tail can be significant and the TOST step decides. The margin is `0.5` pooled
standard deviations, so the shift applied to each null is `0.5 x SD` in raw
units, and the studentized size of that shift is `0.5 x SD / SE`, which grows
like `sqrt(n)`. Measured on samples from one distribution at 4000 resamples:

| n per group | `p_tost`, margin in SE units | `p_tost`, margin in pooled SD |
| --- | --- | --- |
| 8 | 0.997 | 0.993 |
| 30 | 0.464 | 0.071 |
| 100 | 0.623 | 0.004 |
| 400 | 0.709 | 0.0002 |
| 2000 | 0.918 | 0.0002 |

The pooled-SD column crosses `alpha = 0.05` between n=30 and n=100, so n=8 gives
`inconclusive` and n=400 gives `equivalent`. The SE column is flat in n because
the studentized statistic is a pivot whose null spread stays near 1 at every
sample size, which is why the margin is an effect size rather than a number of
standard errors.

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
studentization buys. Measured behaviour across shapes: about 4% here, about 3%
when the smaller group has the smaller variance, and rising to about 10% at
n=20 against n=200 when the smaller group carries a 3x larger spread — the
small-sample limitation recorded in ORCHESTRATOR.md.

### `test_guards_short_circuit_without_resampling`

**Inputs and derivation**, one row per guard:

| Input | Guard | Why |
| --- | --- | --- |
| `x=[0.4]`, `y=[0.1, 0.2, 0.3]` | `insufficient_group_size` | `nx = 1`, and `np.var(ddof=1)` needs at least 2 observations. |
| `x=y=[0.3, 0.3, 0.3]` | `zero_pooled_variance` | All six pooled values equal, so both the spread and the mean difference are zero. |
| `x=[0.9]*10`, `y=[0.1]*30` | `degenerate_observed_scale` | Both groups internally constant with different means; the standard error is floating-point noise near `1e-17`. |
| `x=[0.9,0.8,0.7,0.6]`, `y=[0.1,0.2]` | `insufficient_permutation_support` | `C(6, 2) = 15 < 2500`. |

Each returns `statistic is None`, `n_resamples == 0`, `converged is False`, and
decision `inconclusive` — a guard means nothing was established, which is
distinct from having shown the means to be equivalent.

### `test_null_skewness_is_measured_and_matches_scipy`

**Inputs:** 700 concordant heights from `N(0.42, 0.10)`, 19 discordant1 of which
15 are from `N(0.5, 0.1)` and 4 from `N(25, 5)`; 2500 resamples, `bfn`.

**Derivation:** the permutation statistic depends sharply on how many of the four
extreme heights land in the 19-slot group, so the null separates into clusters
rather than one smooth curve: about 90% of permutations put all four in the large
group (`T` near `+1.9`), 10% put one in the small group (`T` near `-0.9`), and
0.5% put two (`T` near `-1.4`). A third moment over that mixture is strongly
negative, so `|null_skew| > 1`.

The expected value is not derived independently — it is `scipy.stats.skew` over
the statistics drawn from the same seed. That is the point: `null_skew` is
accumulated from running power sums so batches can be discarded, and the test
pins that accumulation against the reference computed from the retained draws.

### `test_null_skewness_is_near_zero_for_a_symmetric_null`

**Inputs:** two 150-observation samples from `N(1.0, 1.0)`, 4000 resamples.

**Derivation:** equal group sizes drawn from one symmetric family give a
permutation null that is symmetric about zero, so its population skewness is 0.
The `0.15` band is Monte Carlo slack: the standard error of a sample skewness is
about `sqrt(6/n) = 0.039` at n = 4000, so the bound is roughly 4 standard errors.

### `test_adaptive_run_grows_batches_until_it_converges`

**Inputs:** `x ~ N(2.0, 0.2^2)` and `y ~ N(0.5, 0.2^2)`, both n=80,
`min_resamples=1000`, `max_resamples=20000`.

**Derivation:** the samples are separated by more than seven pooled standard
deviations, so no permutation approaches the observed statistic and
`p_greater` lands at the floor `1/1001 ~ 1.0e-3`. The Wilson interval around a
count of 1 in 1001 is roughly `[0.0003, 0.0056]`; doubled by Bonferroni it stays
far below `alpha = 0.05`, so `alpha` is outside it after the very first batch.
The run therefore stops with `batches == 1`, `n_resamples == 1000`, and decision
`greater` — an easy case must not spend the ceiling.

### `test_undecided_runs_grow_their_batches_until_the_budget_is_reached`

**Inputs:** 30 draws from `N(1.4, 1.0^2)` against 30 from `N(1.0, 1.0^2)`
(sample seed 3), `min_resamples=100`, `max_resamples=1000`, resampling seed 103.

**Derivation:** the 0.4 shift is marginal at these sample sizes, so the
corrected p-value stays near `alpha` and the Wilson interval never excludes it
— the run draws every batch it is allowed. It therefore ends `converged=False`
with `note = "max_resamples_reached"`, and takes 6 batches rather than 1.

Two properties are asserted, neither pinning the growth factor. **Growth:** a
schedule that repeated the opening batch would total exactly
`batches x 100 = 600`, so a larger total shows the batches grew. **Budget:** the
first five batches are 100, 125, 156, 195, 243, cumulating to 819 — short of the
1000 budget, so a sixth batch of 303 is drawn *whole*, ending at 1122. The
assertion is the bracket `1000 <= n_resamples < 2000`: the run does not stop
before the budget is met, and overshoots it by at most one batch. The exact
1122 is left unasserted because it encodes the growth factor.

### `test_seeded_runs_are_reproducible`

**Inputs:** one random sample pair, run twice through
`run_studentized_permutation_test` with `min_resamples=800`,
`max_resamples=2000`, and a fresh `default_rng(42)` each time.

**Derivation:** every random draw in the test comes from the passed generator,
so two runs seeded alike must agree on every field -- not approximately, but
exactly, including the resample count the adaptive rule stopped at. This is the
property the whole per-triplet seeding scheme rests on: the orchestrator hands
each triplet a stream derived from `(seed, triplet)`, so a triplet's result
cannot depend on how many workers ran or in what order they finished.

### `test_equivalence_test_disabled_reports_no_tost`

**Inputs:** two 40-point standard-normal samples (so neither direction is
significant), with `equivalence_test=False` and a fixed 2000 resamples.

**Derivation:** the equivalence step is what distinguishes `equivalent` from
`inconclusive` when no direction is significant. Switched off, there is nothing
left to make that distinction, so `p_tost` must be `None` and the decision must
fall through to `inconclusive` -- never to `equivalent`, which would be an
equivalence claim no test supported.

### `test_bootstrap_resample_budget_is_reduced_but_always_usable`

**Inputs:** the configured pairs `(2500, 25000)`, `(2, 3)`, `(1, 1)`,
`(100, 100)`, `(7, 1000)`, `(999, 1001)`.

**Derivation:** the expected values are not computed from the divisor, because
the divisor is a performance knob — the assertions are the three properties an
iteration budget has to satisfy whatever it is set to. *At least 1*: integer
division sends a small configured budget to 0, so the floor is what stops an
iteration drawing no resamples at all — `(2, 3)` and `(1, 1)` are the cases that
reach it. *Strictly below the configured minimum once that is at least 2*: any
divisor above 1 reduces such a value, so this is what would catch a divisor of 1
silently restoring the full per-iteration cost. *Bounds in order*: `(100, 100)`
and `(999, 1001)` are the cases where the scaled maximum would otherwise land
below the scaled minimum, leaving an adaptive run no range to grow through, and
the implementation holds the maximum at the minimum instead.

Confirmed to have teeth by mutation: setting the divisor to 1 fails five of the
six parametrized cases.

### `test_bootstrap_resample_budget_never_shrinks_as_the_budget_grows`

**Inputs:** configured minima of 1, 2, 5, 10, 100, 2500 and 25000, each paired
with ten times itself as the maximum.

**Derivation:** floor division is monotonic and the floor is a constant, so both
ends of the scaled budget must be non-decreasing across that sequence. This is
the property a reader actually depends on — that configuring a larger budget
cannot give the bootstrap a smaller one — and it holds for any divisor.

## tests/orchestrator/test_orchestrator_trees.py

### `test_clean_and_save_trees_preserves_a_well_supported_tree`

**Input:** the species-tree fixture, `min_avg_support=0.5`.

**Derivation:** the tree carries no internal support labels, so there is nothing
for the filter to reject and nothing to strip. The output must therefore be the
input string plus a trailing newline:
`(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);\n`.

### `test_clean_and_save_trees_drops_low_average_support`

**Input:** `low_support_tree_file`, `min_avg_support=0.5`, containing

```
(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup);
(((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);
```

**Derivation:** mean internal support for tree 0 is
`(0.95 + 0.99 + 0.98)/3 = 2.92/3 = 0.9733 >= 0.5` → kept. For tree 1 it is
`(0.30 + 0.20 + 0.40)/3 = 0.90/3 = 0.3000 < 0.5` → dropped. Exactly one tree
survives, with leaf set `{TaxaC, TaxaD, TaxaF, TaxaG, OutGroup}`, and because
cleaning strips support labels the substring `0.95` must not appear in the
output.

### `test_root_tree_on_outgroup_prunes_and_reports_ingroup`

**Input:** the cleaned species tree, outgroup `["OUT"]`.

**Derivation:** see "Orchestrator fixture: species tree" above — `OUT` is excluded,
no requested outgroup is missing, the ingroup is `[A, B, C, D]`, and the pruned
Newick is `(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;` because the `((A,B),C)`
clade absorbs `0.1 + 0.2` when the `(D,OUT)` node dissolves.

### `test_generate_triplets_and_species_subtrees`

**Input:** the pruned species tree `(((A,B),C),D)`.

**Derivation:** 4 ingroup taxa give `C(4,3) = 4` triplets, enumerated over the
sorted taxa: `(A,B,C)`, `(A,B,D)`, `(A,C,D)`, `(B,C,D)`. In each of these the
first two listed taxa are already the species-tree sister pair, so ABC
normalization is the identity and `triplets == raw_triplets`.

One cached geometry over the species tree answers all four triplets:
`triplet_subtree_shape` reports each one's sister pair and the edges its induced
subtree carries, and `_format_triplet_subtree_newick` writes that out. Each
subtree keeps the triplet's taxa and sums the edges along the paths that
collapse when the taxa between them are dropped, which is what the `Why` column
below states. The subtree's root edge is its own case: suppressing the
unifurcations above the triplet's LCA collapses the whole path from the tree
root into one edge, so `(A,B,C)` carries `0.3 + 0.5 = 0.8` -- the LCA's depth
plus the tree root's own edge -- while the three triplets containing `D` are
rooted at the tree root itself and keep its `0.5` unchanged:

| Triplet | Species subtree | Why |
| --- | --- | --- |
| (A,B,C) | `((A:0.1,B:0.1):0.1,C:0.2):0.8;` | root edge `0.3 + 0.5 = 0.8` |
| (A,B,D) | `((A:0.1,B:0.1):0.4,D:0.1):0.5;` | A,B clade edge `0.1 + 0.3 = 0.4` |
| (A,C,D) | `((A:0.2,C:0.2):0.3,D:0.1):0.5;` | A absorbs `0.1 + 0.1 = 0.2` |
| (B,C,D) | `((B:0.2,C:0.2):0.3,D:0.1):0.5;` | B absorbs `0.1 + 0.1 = 0.2` |

### `test_clean_and_save_gene_trees_roots_every_tree_on_the_outgroup`

**Input:** the 12 gene trees, outgroup `["OUT"]`, `min_avg_support=0.5`.

**Derivation:** none of the trees carry support labels, so all 12 survive.
Rerooting at the outgroup attachment point moves `OUT`'s edge into the ingroup
clade's edge and zeroes `OUT`, so every line must end `OUT:0);`. Tree 0:
`0.10 + 0.50 = 0.6` → `((((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6,OUT:0);`.
Tree 3: `0.10 + 0.55 = 0.65` →
`((((B:0.4,C:0.4):0.1,A:0.6):0.1,D:0.35):0.65,OUT:0);`.

### `test_extract_triplet_subtree_selects_the_triplet_taxa`

**Input:** the first cleaned gene tree, triplet `("A","B","C")`.

**Derivation:** extraction retains only the requested taxa, so the resulting
leaf-label set must be exactly `{A, B, C}` — `D` and `OUT` are dropped. No run
calls this function; it is the reference the geometry parity tests measure the
cached path against, so it needs a correctness check of its own rather than only
being compared to.

## tests/orchestrator/test_orchestrator_tree_parity.py

**Inputs:** the three `triplet_comparison_cases`, each a
`(newick, triplet)` pair:

1. `((A:1.0,B:1.0):2.0,C:3.0,D:4.0);` with `(A,B,C)`
2. `(((A:0.1,X:0.1):0.2,(B:0.1,Y:0.1):0.2):0.3,(C:0.1,Z:0.1):0.4);` with `(A,B,C)`
3. `((A:0.5,B:0.5):0.5,(C:0.2,D:0.2):0.8);` with `(A,B,C)`

**Derivation:** these are the suite's only tests whose expectation comes from a
second implementation rather than a closed-form definition — that is the point.
Patristic distance is well defined (`d(a,b) = d(a,root) + d(b,root) -
2 x d(mrca,root)`), so DendroPy's `extract_triplet_subtree` and BioPython's
`Clade.prune` must agree on every pairwise distance. `test_triplet_branch_lengths_match`
checks agreement after a DendroPy → Newick → BioPython round trip (catching
serialization precision loss); `test_triplet_collapse_consistency_dendropy_vs_biopython`
compares the two collapse implementations directly to `abs=1e-12`.

Worked check for case 1, pair `(A,B)`: both taxa hang off the `(A,B)` node at
depth 2.0, so `d = 1.0 + 1.0 = 2.0`. Pair `(A,C)`: `d = (1.0 + 2.0) + 3.0 = 6.0`.

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
- **KS correction, `detailed` only:** asserted as a relation rather than a
  literal — `ks_p_value_corrected == min(1.0, ks_p_value x 4)`. (For (A,B,C)
  the raw KS p-value is `~0.01667`, giving `~0.06667`; for the D triplets dis1
  is empty so the KS test short-circuits to `1.0` and stays `1.0`.)
  `perm_decision` is populated on every row because the direction test ran.
- **`efficient`:** the run's default correction is `bfn`, an inline method, so
  the stream judges each count gate on the exactly corrected `p × 4 = 1.0`,
  finds it failed, and measures nothing below it: `ks_p_value`,
  `ks_p_value_corrected`, `ks_significant` and `perm_decision` are `None` and
  `perm_note` is `direction_test_not_consulted`. Everything asserted above the
  gate is identical between the modes.
- `most_frequent_matches_concordant` is `True` everywhere, since `7 >= 3, 2` and
  `12 >= 0, 0`.
- `analyzed_trees == 12` for every triplet, because all 12 gene trees contain
  all five taxa.

### Output-shape tests

- `test_run_orchestrator_writes_results_tsv` — the file must exist, its header must
  begin with `triplet` and contain `classification` and `bootstrap_value`, and
  the data-row count must equal `len(results)` (4).
- `test_no_bootstrap_omits_the_bootstrap_columns` — with `bootstrap=False` the
  writer skips the bootstrap columns, so `bootstrap_value` and `all_bootstrap`
  must be absent while `classification` remains and all 4 triplets are produced.
- `test_generate_summary_stats_writes_tsv` — the column count is derived from
  the contract: 7 statistics (mean, median, mode, variance, entropy, min, max)
  x 3 metrics (avg_tree_height, internal_branch, sister_distance) x 3 topology
  classes (concordant, discordant1, discordant2) = **63** metric columns.
- `test_bootstrap_debug_mode_writes_debug_columns` — the five named debug
  columns must appear and at least one result must have a populated
  `bootstrap_dct_stats`.
- `test_consolidation_preserves_run_outputs` — consolidation writes into
  `consolidation/`, so the run's own files (results TSV, `metrics.txt`, both
  processed trees) must all still exist afterwards.
- `test_parallel_runs_match_serial` — bootstrap seeding is per-triplet and
  derived from the run seed, so mode and worker count cannot change any value;
  every compared field must be equal, bootstrap included.

## tests/orchestrator/test_orchestrator_preflight.py

### Shared inputs

The species tree is `(((A:1,B:1):1,C:1):1,(D:1,OUT:1):1);`. Removing the
outgroup `OUT` leaves ingroup `{A, B, C, D}`, so the check enumerates
`C(4,3) = 4` triplets: `A,B,C`, `A,B,D`, `A,C,D`, `B,C,D`.

The clean gene-tree file holds two trees, both containing all five taxa:

```
((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);
((((A:1,C:1):1,B:1):1,D:1):1,OUT:1);
```

The defective file holds three, each planted with exactly one problem class:

```
((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);   # well formed
(((A:1,B:1,C:1):1,D:1):1,OUT:1);       # polytomy over A,B,C
(((A:1,B:1):1,C:1):1,MISSING:1);       # no OUT label
```

### `test_clean_inputs_pass_with_no_issues`

**Inputs:** the clean pair above, `outgroups=["OUT"]`.

**Derivation:** both gene trees root on `OUT`, so
`counters["gene_tree.rooted"] == 2` and `gene_tree.total_checked == 2`. Every
one of the 4 triplets resolves a sister pair in the species tree, so
`triplets_checked == 4`. Each gene tree is fully resolved, so no triplet check
fails and `issues == []`, which makes `passed` `True`.

### `test_report_is_written_only_when_an_output_dir_is_given`

**Inputs:** the clean pair, checked once with `output_dir` set to a created
directory and once with `output_dir=None`.

**Derivation:** the writer joins `output_dir` with the module constant
`PREFLIGHT_REPORT_FILENAME` (`preflight_data_check.txt`) and writes
`report_text` verbatim, so the file content and `report_text` must be equal and
`report_path` must equal that joined path. The write branch is guarded on
`output_dir is not None`, so the second call leaves `report_path` `None` while
still building the same `report_text` -- the check itself does not depend on
where its output goes.

### `test_detects_polytomy_and_missing_outgroup`

**Inputs:** the defective trio, `outgroups=["OUT"]`.

**Derivation:** gene tree 3 contains no `OUT`, so `_root_tree_on_any_outgroup`
returns no used outgroup → one `gene_tree.rooting_failed`, and that tree is
skipped before any triplet check. Trees 1 and 2 root, so
`gene_tree.rooted == 2` while `gene_tree.total_checked == 3`. Tree 2 collapses
A, B and C into a single polytomous clade, so all three pairwise LCAs of triplet
`A,B,C` are the same node and `triplet_resolution` reports it unresolved → one
`triplet.unresolved_rooted_sister_pair`. The other three triplets each contain
`D`, which sits outside the polytomy, so they still resolve — hence a count of
exactly 1, not 4. The message is formatted with the enumeration index (`Gene
tree #2`) and the comma-joined triplet (`A,B,C`).

The pair accounting follows from the same reading. Two trees root, and each
carries all four ingroup taxa, so the check looks at `4 x 2 = 8` triplet/
gene-tree pairs: tree 1 resolves all 4, tree 2 resolves the three containing
`D` and fails on `A,B,C`. That gives 7 usable, 1 unresolved, and 0 skipped for
an absent taxon — and the three must sum to the 8 pairs seen, which is the
property worth pinning: a pair that is neither measured nor reported would
otherwise vanish silently between the counters.

### `test_triplet_filter_entries_are_validated`

**Inputs:** the clean pair plus a filter file containing `A,B,C`, `A,B,NOPE`,
`A,B,OUT`.

**Derivation:** `NOPE` is not among the species-tree labels, so set difference
against them is non-empty → one
`triplet_filter.taxa_missing_in_species_tree`. `OUT` is in the outgroup set →
one `triplet_filter.includes_outgroup`. Both lines are dropped rather than
checked, leaving only `A,B,C`, which normalizes successfully, so
`triplets_checked == 1`.

### `test_impossible_checks_raise`

**Inputs (parametrized):** `outgroups=[]`, and `outgroups=["NOT_PRESENT"]`.

**Derivation:** the empty list fails the explicit guard at the top of
`run_preflight_data_check`. `NOT_PRESENT` is absent from the species tree, so
`_root_tree_on_any_outgroup` returns `used_outgroup=None` and the function
raises rather than reporting — without a rooted species tree there is no
triplet normalization to perform, so there is nothing to report on.

### `test_multi_tree_species_file_raises`

**Input:** the species-tree string written twice.

**Derivation:** `read_tree_file` returns 2 trees, and
`_load_single_species_tree` requires exactly 1.

### `test_runner_preflight_mode_skips_analysis`

**Input:** a config dict with `preflight_data_check: True` and
`preflight_triplet_cap: 3` over the defective trio.

**Derivation:** `run_orchestrator` prepares the output directory and then
returns `_run_preflight_only(...)` before any tree cleaning, so the only write
into that directory is the report. Listing the directory must therefore yield
exactly `["preflight_data_check.txt"]` — no `metrics.txt`, no processed trees,
no results TSV. `passed` is `False` because the defective trio yields issues.

The cap: the species tree's ingroup is `A, B, C, D`, so `combinations(·, 3)`
yields `C(4, 3) = 4` triplets. `_load_target_triplets` keeps the first
`max_triplets` when `0 < max_triplets < 4`, so a cap of `3` gives
`triplets_checked == 3` and appends one `analysis.triplet_cap_applied` issue.
Had the runner ignored the config and used the module default (15,000), all 4
would be checked and the cap issue would be absent, so both assertions would
fail — which is what makes `3` rather than the default the right value here.

### `test_runner_returns_none_when_preflight_cannot_run`

**Input:** the same config with `outgroup="NOT_PRESENT"`.

**Derivation:** the `ValueError` raised above is caught in
`_run_preflight_only`, which returns `None`, so the runner returns `None`
rather than propagating.

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
every key — including the resolved paths, since both resolve the same relative
strings against the same working directory. The comparison is whole-dict
equality, so a key that resolved differently on the two paths, or a default
that one path filled and the other did not, fails without the test naming
either. It was checked before being written that the two dicts are in fact
equal today.

### `test_config_only_keys_and_nested_blocks_flatten_from_a_file`

**Input:** a YAML file with `discordant_test: z-test`,
`tree_height_calculation_strategy: SIS`, `min_support_value: 0.9`,
`generate_summary_stats: true`, `alpha_dct: 0.02`, `seed: 7`, and
`bootstrap_options: {iterations: 25, debug_mode: true, summary_only: true}`.

**Derivation:** the first four keys have no CLI flag, so the file is the only
way to set them, and each is read back as written. The nested block is
flattened onto `bootstrap_iterations = 25`, `bootstrap_debug_mode = True`,
`bootstrap_summary_only = True`. The nine values are compared as one dict
against the nine written, so any one of them resolving to something else fails
naming the key.

### `test_config_file_wins_over_cli`

**Input:** a config file setting `alpha_dct: 0.03` and file-specific tree paths,
plus conflicting CLI flags `alpha_dct=0.5` and `alpha_perm=0.5`.

**Derivation:** in config-file mode the file supplies everything, so
`alpha_dct` is `0.03` (not `0.5`). The decisive check is `alpha_perm`: the file
omits it, so it takes the default — whatever that is — and must not be the CLI's
`0.5`. The assertion is `!= 0.5` rather than `== 0.05` so the test does not
pin the default; it pins only that the CLI value was ignored rather than
merged. The resolved species path must come from the file.

### `test_outgroup_accepts_single_comma_separated_and_list_forms`

**Inputs and derivation:** every accepted shape flattens the same way — each
entry is coerced to a string, split on commas, stripped, and empty pieces
dropped.

| Input | Result | Why |
| --- | --- | --- |
| `"OUT"` | `["OUT"]` | Single label, no comma to split on. |
| `"Out1,Out2"` | `["Out1", "Out2"]` | Comma-separated string. |
| `" Out1 , Out2 ,"` | `["Out1", "Out2"]` | Padding stripped, trailing empty piece dropped. |
| `["Out1", "Out2"]` | `["Out1", "Out2"]` | List of labels. |
| `("Out1", "Out2")` | `["Out1", "Out2"]` | Tuples accepted alongside lists. |
| `["Out1,Out2", "Out3"]` | `["Out1", "Out2", "Out3"]` | List entries are themselves split, so the two forms compose. |

### `test_outgroup_rejects_empty_and_non_label_values`

**Inputs and derivation:** `None`, `""`, `"  "`, `","`, `[]`, `["", "  "]`, and
`42` all reduce to an empty label list — the first six because every piece is
blank after stripping, and `42` because it is neither a string nor a
list/tuple/set and so contributes no entries. An empty result raises
`ConfigError` naming `outgroup` rather than silently producing an unrooted run.

### `test_invalid_values_are_rejected_by_field_name`

**Inputs:** the required-keys file with one key overridden per row.

| Override | Validator reached | Why it fails | Message must contain |
| --- | --- | --- | --- |
| `species_tree_path: null` | `_validate_required_path` | `None` is treated as absent. | `species_tree_path` |
| `p_value_correction: true` | `_validate_choice` | YAML 1.1 turns a bare `yes` into `True`; the boolean-to-choice mapping covers `False` → `no` but `True` matches no choice. Written as `true` here since the payload is dumped with `yaml.safe_dump`. | `must be one of` |
| `pipeline_mode: "fast"` | `_validate_choice` | Not one of `efficient`/`detailed`. | `pipeline_mode` |
| `seed: "abc"` | `_validate_optional_int` | Accepts `None` or an `int`; also rejects `bool`, since `isinstance(True, int)` holds and `seed: true` is a mistake rather than a seed of 1. | `seed` |
| `preflight_triplet_cap: -1` | `_validate_non_negative_int` | Accepts `int >= 0`; `0` is the documented "no cap", so `-1` is the nearest value with no meaning. | `preflight_triplet_cap` |

The `p_value_correction` row is the guard that keeps the boolean mapping from
laundering an invalid value into a valid one.

### `test_shipped_sample_configs_resolve`

**Inputs:** the two orchestrator sample configs under `sample_configs/`.

**Derivation:** these go through the same `load_orchestrator_config` a user
invokes with `-c`, so anything the validator would reject surfaces here. The
asserted value is what the samples state literally — a non-empty outgroup list
(`["OutGroup"]` for the minimal sample, `["Out1", "Out2"]` for the full one,
written there as a YAML list to exercise that form). The placeholder tree paths
need not exist: `_validate_required_path` only checks that the field is a
non-empty string before resolving it.

### `test_full_sample_config_names_every_runtime_key_at_its_default`

**Inputs:** `orchestrator_full.yaml` parsed twice — once as raw YAML for the set
of documented keys, once through `load_orchestrator_config` for the set of
runtime keys — plus a required-keys-only file loaded the same way.

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
omitted elsewhere — together they pin that every key is both present and at its
default. It was checked that the sample satisfies this today.

### `test_p_value_correction_accepts_yaml_bare_word_no`

YAML 1.1 resolves the bare word `no` (and `No`, `NO`) to boolean `False`, so the
normalizer maps `False` back to the string `"no"` before the choice check;
quoted `"no"` and `bfn` arrive as strings and pass straight through. Writing
`no` unquoted is the natural spelling for "no correction", which is why the
mapping exists.

### `test_parser_flags_resolve_into_their_config_values`

Parsing `-st s -gt g -og OUT --alpha-dct 0.01 --alpha-ks 0.2
--p-value-correction fdr_bh --pipeline-mode detailed --alpha-perm 0.02
--no-overwrite --preflight-data-check --preflight-triplet-cap 0` must yield
those exact values
with `config_file is None`; `0` is chosen for the cap because it is the one
value `_validate_non_negative_int` accepts that differs from the default and is
also the documented "no cap" spelling.

## tests/test_config_trunk.py

### `test_resolve_path_handles_absolute_relative_and_home`

**Inputs:** an absolute `tmp_path/species.nwk`; the relative `genes.nwk` after
`chdir` into `tmp_path`; and `~/data.nwk`.

**Derivation:** absolute paths resolve to themselves; relative paths resolve
against the current working directory, giving `tmp_path/genes.nwk`; `~` expands
to `Path.home()`, so the result is `Path.home()/data.nwk` and contains no `~`.

### `test_load_raw_config_reads_json_and_yaml`

**Inputs:** `{"input_path": "data.tsv", "cv_folds": 5}` written as JSON, and the
equivalent two-line YAML.

**Derivation:** both formats must parse to the identical Python mapping — the
loader's only job is format dispatch.

### `test_load_raw_config_rejects_missing_unsupported_and_non_mapping`

**Inputs:** a nonexistent path; a `.txt` file; a JSON file whose root is
`["a","b"]`.

**Derivation:** the existence check runs first → `FileNotFoundError`. The suffix
check runs next → `ConfigError("Config file must be .json, .yaml, or .yml")`. The
root-type check runs last → `ConfigError("Config root must be a key/value
object")`.

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
missing positive suffix — `2` — creating `results_2` and leaving `fresh.txt`
untouched in the original directory.

## tests/test_ml_labels_and_metrics.py

### `test_parse_classes_round_trips_bitstrings`

**Input:** `["100001", "000000", "111111", " 010010 "]`.

**Derivation:** each label is stripped and expanded into six integers, giving a
4x6 matrix. Re-joining each row must reproduce the trimmed label. Set bits per
label: `100001 → 2`, `000000 → 0`, `111111 → 6`, `010010 → 2`, so the matrix sum
is `2 + 0 + 6 + 2 = 10`.

### `test_bit_label_titles_keep_the_taxon_letters_upper_case`

**Inputs:** the six `BIT_LABELS` entries — `ghost_into_A`, `ghost_into_B`,
`inflow_into_A_from_C`, `inflow_into_B_from_C`, `outflow_from_A_to_C`,
`outflow_from_B_to_C`.

**Derivation:** the transform is `replace("_", " ")` followed by upper-casing
character 0 only, so `inflow_into_A_from_C` → `inflow into A from C` →
`Inflow into A from C`. The expected values are written out per label rather
than computed, because the point is the one spelling the obvious implementation
gets wrong: `str.capitalize()` upper-cases the first character *and lower-cases
the rest*, which would yield `Ghost into a` and rename taxon `A`. Every label in
the set carries at least one trailing capital, so any label would catch it —
they are all listed so the failure names which one broke.

### `test_every_bit_label_has_a_title`

**Inputs:** the whole `BIT_LABELS` tuple.

**Derivation:** six inputs must give six distinct outputs — a collision would
put the same title on two panels of the per-bit figure. No output may keep an
underscore (the transform is total, not a lookup table with gaps) and each must
start upper-case.

### `test_64_class_matrix_orders_classes_by_set_bits`

**Inputs:** `y_true` rows `[1,1,0,0,0,0]`, `[0,0,0,0,0,1]`, `[1,1,1,1,1,1]` and
`y_pred` rows `[0,0,0,0,1,1]`, `[0,0,0,0,0,1]`, `[1,1,1,1,1,1]` — classes
`110000 → 000011`, `000001 → 000001`, `111111 → 111111`.

**Derivation:** the builder sorts the 64 six-bit strings on
`(count("1"), label)`, so the list opens `000000`, runs the six single-bit
classes `000001 … 100000`, then the fifteen two-bit classes, and so on to
`111111`: set-bit counts `[1, 6, 15, 20, 15, 6, 1]` — the binomial row for
`n = 6` — and never decreasing. Within a count the tie-break is the string
itself, so each group is lexically sorted. The three rows were chosen so that
each of the three (true, predicted) pairs sits somewhere binary order and
set-bit order disagree: in binary order `110000` is index 48 and `000011` index
3, but after the sort they are at positions 21 and 7. The matrix is indexed by
`labels.index(...)` on both axes, so the assertion passes only if rows *and*
columns were permuted with the labels; a sort that reordered the label list but
left the counts at their binary indices would report the `110000 → 000011`
count at position `(48, 3)`, which the labels now name `(101011, 000100)`. The
total of 3 confirms nothing was dropped or duplicated by the reindexing.

### `test_is_valid_bitstring`

**Derivation:** validity is `len(value) == 6 and set(value) <= {"0","1"}`. Hence
`000000`, `111111`, `100001` pass; `10000` (5 chars), `1000010` (7 chars),
`100002` and `10000a` (non-binary), and `""` all fail.

### `test_evaluate_predictions_matches_hand_computed_metrics`

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

### `test_build_prediction_rows_reports_matched_label_count`

**Inputs:** the same two arrays.

**Derivation:** `matched_label_count` counts positionwise agreement: row 0 → 6,
row 1 → 5. `exact_match` is 1 only for row 0. Label strings are the joined bits:
row 0 `100001`/`100001`; row 1 `010000` true vs `000000` predicted. The per-bit
columns expose that difference as `true_ghost_into_B = 1`,
`pred_ghost_into_B = 0`.

### `test_summarize_distribution_counts_and_fractions`

**Input:** `["100001", "000000", "100001", "111111"]` (4 labels, one repeated).

**Derivation:** counts are `100001 → 2`, `000000 → 1`, `111111 → 1`; fractions
divide by 4, giving `0.5`, `0.25`, `0.25`. Output keys are sorted, so the order
is `000000`, `100001`, `111111`.

### `test_bit_distribution_counts_positives_per_bit`

**Input:** `[[1,0,0,0,0,1], [1,1,0,0,0,0]]`.

**Derivation:** column sums over 2 rows. Bit 0 is set twice → count 2, fraction
`2/2 = 1.0`. Bit 1 once → count 1, fraction `0.5`. Bit 2 never → count 0.

### `test_build_feature_importance_rows_sorts_descending`

**Input:** features `("feature_1","feature_2","feature_3")` with scores
`[0.2, 0.5, 0.3]`.

**Derivation:** pair each name with its score and sort by score descending →
`feature_2` (0.5), `feature_3` (0.3), `feature_1` (0.2).

### `TestAutoCvFolds`

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

### `test_read_tsv_rows_*`

**Inputs:** a TSV containing only `class\tfeature_1`; and a TSV with that header
plus two data rows.

**Derivation:** the reader requires at least one data row, so the header-only
file raises `ValueError("...no data rows")`. The two-row file yields one dict per
row keyed by the header fields.

### `test_select_feature_names_excludes_the_target_column`

**Input:** `["class", "feature_1", "feature_2", "dis1_topology"]`, target
`"class"`.

**Derivation:** the target is removed and the remaining header order is
preserved → `("feature_1", "feature_2", "dis1_topology")`.

## tests/orchestrator/test_orchestrator_consolidation.py — ghost bar colouring

### Shared scenario

`_ghost_colour_scenario_results()` builds three results against the species tree
`(((A:1,B:1):1,C:1):1,D:1);`. The mapping rules in `_map_event` turn them into:

| Result | Rule | Produces |
| --- | --- | --- |
| `(A,B,C)` ghost, dis1 `BC` | ghost + `BC` → ghost target is the A-taxon | ghost target `A`, weight 0.8 |
| `(A,B,C)` inflow, dis1 `AC` | inflow + `AC` → edge `(c_taxon, a_taxon)` | sampled edge `(C, A)`, weight 0.5 |
| `(D,B,C)` ghost, dis1 `BC` | ghost + `BC` → ghost target is the A-taxon | ghost target `D`, weight 0.4 |

So `A` is both a ghost target and the target of sampled edge `(C, A)` → flag
`1` → `GHOST_WITH_SAMPLED_BAR_COLOR` (cividis low end, dark blue). `D` is a
ghost target with no sampled edge pointing at it → flag `0` →
`GHOST_ONLY_BAR_COLOR` (cividis high end, yellow).

### `test_sampled_introgression_presence_flags_targets_with_sampled_edges`

**Inputs:** `taxa_order = ["A","B","C","D"]` and
`{("C","A"): 0.6, ("A","B"): 0.0, ("D","C"): 0.3}`.

**Derivation:** edges are keyed `(source, target)`, so the flag looks at the
second element. `A` is the target of `(C,A)` with weight `0.6` → `1`. `B` is the
target of `(A,B)` but the weight is `0.0`, which is falsy, so an edge that
carries no support does not count → `0`. `C` is the target of `(D,C)` with
`0.3` → `1`. `D` is never a target → `0`. Expected:
`{"A": 1, "B": 0, "C": 1, "D": 0}`.

### `test_ghost_strength_tsv_records_sampled_introgression_flag`

**Inputs:** the shared scenario above.

**Derivation:** the header gains a third field, so it must be exactly
`["target_taxon", "raw_strength", "has_sampled_introgression"]`. Every taxon in
the plot order is listed (`{A, B, C, D}`), because the sheet is keyed by
`taxa_order` rather than by which taxa have ghost signal. From the table above,
`A → 1` and `D → 0`. Both `A` and `D` have non-zero ghost strength, confirming
the flag is an added column rather than a replacement for the strength value.

## tests/orchestrator/test_orchestrator_shape.py

### `test_shape_moments_match_scipy`

- **Input** — 500 draws from `numpy.random.default_rng(4).lognormal(0, 0.7)`.
- **Expected `skew`** — `scipy.stats.skew(values)`, the biased sample third
  standardized moment `m3 / m2^1.5`. Exact equality, since the module calls the
  same function.
- **Expected `excess_kurtosis`** — `scipy.stats.kurtosis(values, fisher=True)`,
  i.e. `m4 / m2^2 - 3`, so a normal reads `0`.
- **Signs** — a lognormal with `sigma = 0.7` has population skewness
  `(e^{sigma^2} + 2) sqrt(e^{sigma^2} - 1) ~ 2.5` and positive excess kurtosis,
  so both must come out above zero.

### `test_modality_test_rejects_only_a_well_separated_mixture`

- **Inputs** — 400 draws each of `normal(0, 1)`, `lognormal(0, 0.6)`,
  `exponential(1)`, `gamma(2, 1)`, and a concatenation of 200 `normal(0, 1)`
  plus 200 `normal(4, 1)` draws, all from one `default_rng(17)` stream. The
  modality bootstrap is driven by a separate `default_rng(2)`.
- **Derivation of the expectation** — the first four are unimodal by
  construction (a gamma with shape `> 1` has its mode at `(k-1) * theta`; a
  lognormal at `e^{mu - sigma^2}`), so a test holding its nominal level must not
  reject them. The mixture's components are 4 pooled SD apart, which puts a
  genuine valley between them; measured `modes_p` values are 0.53/0.38/0.16/0.35
  for the unimodal families and 0.005 — the `1 / (n_resamples + 1)` floor at 200
  replicates — for the mixture. The assertion is the side of `alpha = 0.05` each
  falls on, not the value.
- **Why not a mode count** — the same four unimodal samples give KDE peak counts
  of 1, 4, 3 and 2 at Scott's bandwidth, so the count alone would call three of
  them multimodal.

### `test_tail_index_recovers_known_tail_shapes`

- **Inputs** — 4000 draws each: `pareto(3) + 1`, `exponential(1)`,
  `uniform(0, 1)`, from `default_rng(23)`.
- **Expected `xi`** — a Pareto with index `a` has survival `x^{-a}`, whose
  generalized-Pareto tail index is `1 / a = 1/3`. An exponential tail is the
  `xi -> 0` limit of the GPD. A uniform on `[0, 1]` has a finite upper endpoint
  reached linearly, which is the GPD with `xi = -1`.
- **Tolerance** — `abs=0.25`. The fit sees only the upper decile, so 400 points
  back each estimate; the GPD shape estimator's asymptotic standard error is
  roughly `(1 + xi) / sqrt(k)`, which is about `0.07` here, and the tolerance
  leaves room for the threshold-choice bias on top of that.

### `test_shape_is_not_described_for_small_or_flat_groups`

- **`[1.0] * 5`** — 5 observations, below `SHAPE_MIN_OBSERVATIONS = 20`, so
  every field is `None`.
- **`[2.0] * 50`** — above the size floor but with zero variance, so the KDE has
  no scale and the moments are undefined; every field is `None`.
- **25 normal draws** — clears the size floor, so the moments are reported. The
  upper decile holds `25 - ceil(0.9 * 25) = 3` points, below
  `SHAPE_MIN_TAIL_EXCEEDANCES = 10`, so `tail_xi` alone is `None`.

### `test_shape_is_measured_once_and_not_per_bootstrap_iteration`

- **Input** — one 60/40/30 lognormal triplet analyzed twice, at 5 and 60
  bootstrap iterations, with the same `triplet_seed`.
- **Derivation** — the diagnostics are computed in
  `_run_triplet_pipeline_from_observations`, which runs once for the point
  estimate; `_run_bootstrap_iterations` reaches the data through
  `_iteration_outcome`, which never calls it. The seed for the modality
  bootstrap is the fourth child of the per-triplet `SeedSequence`, and spawning
  a fourth child does not disturb the first three, so both runs draw the same
  stream. The two results must therefore agree exactly, not approximately.

### `test_summary_statistics_tsv_never_carries_shape_columns`

- **Input** — the same triplet, written through `write_summary_statistics_tsv`
  with the diagnostics off and on.
- **Expected column names** — none. No header entry ends in any of
  `SHAPE_FIELD_NAMES`, under either setting. The "on" case is the one that
  matters: the result carries a populated `shape_statistics` dict and the file
  still ignores it.
- **Row width equals header width** — dropping a column block is exactly the
  kind of change that leaves a row misaligned.
- **`concordant_avg_tree_height_mean` and `classification` still present** — the
  63 descriptive columns and the label are untouched.
- **Why** — the diagnostics are undefined below their observation floors: a
  group under `SHAPE_MIN_OBSERVATIONS` (20) has no modality p-value or moments,
  and a group whose upper decile holds under `SHAPE_MIN_TAIL_EXCEEDANCES` (10)
  points has no tail index. `ml_utils.rows_to_matrix` raises on any empty
  feature cell, so those gaps made whole runs unusable as training data.

### `test_results_tsv_carries_shape_columns_only_when_enabled`

- **Input** — one triplet with 60 concordant, 40 discordant1 and 30 discordant2
  lognormal heights, analyzed with `shape_diagnostics` off and on.
- **Expected column set** — the product of `SHAPE_GROUP_LABELS` (3) and
  `SHAPE_FIELD_NAMES` (5), so 15 columns, present exactly when enabled.
- **`con_skew > 0`** — the heights are lognormal, which is right-skewed.
- **`con_modes_p` in `(0, 1]`** — the add-one estimator is bounded below by
  `1 / (n_resamples + 1)` and above by `1`.
- **`dis2_tail_xi` empty** — 30 observations leave `30 - 27 = 3` above the 90th
  percentile, under the 10-exceedance floor.
- **Row alignment** — asserted in both arms, since the columns are appended
  conditionally in two places (header and row) that must stay in step.

## Remaining suites

`tests/orchestrator/test_orchestrator_consolidation.py`, `tests/test_ml_config.py`,
`tests/test_ml_utils.py`, `tests/test_ml_random_forest.py`,
`tests/test_ml_multi_knn.py`, and `tests/test_ml_hyper_tune.py` assert
structural outcomes -- files written, columns present, errors raised, a shipped
sample loading -- and their inputs are the fixtures described in
[TESTS.md](TESTS.md). No test pins a config default; the samples and CONFIG.md
state those. The derivations worth stating explicitly:

- **Introgression map averaging** — a directed pair's denominator is every
  triplet in which both taxa co-occur, which for `n` ingroup taxa is `n - 2`
  triplets; a ghost target's denominator is every triplet containing that taxon,
  which is `(n-1) choose 2`. The reported value is the summed `bootstrap_value`
  over classified triplets divided by that denominator.
- **KNN neighbor capping** — `n_neighbors` cannot exceed the number of training
  samples, so the builder clamps it and the metrics report states the effective
  value.
- **`max_features` / `class_weight` forms** — the accepted set is read off
  scikit-learn's own parameter constraint: `max_features` takes `'sqrt'`,
  `'log2'`, an `int >= 1`, a `float` in `(0.0, 1.0]`, or `None`; `class_weight`
  takes `'balanced'`, `'balanced_subsample'`, a dict, a list of dicts, or
  `None`. `'auto'` is expected to fail because scikit-learn removed it in 1.3.
  `true` is expected to fail because `bool` is a subclass of `int` in Python
  and would otherwise satisfy the `int >= 1` branch as the value `1`. The
  rejected numerics `0` and `1.5` sit just outside the `int >= 1` and
  `(0.0, 1.0]` bounds respectively.
- **Confusion-matrix row normalization** — the input counts are
  `[[3, 1, 0], [0, 0, 0], [1, 1, 2]]`. Row totals are `4`, `0` and `4`. Rows 0
  and 2 divide through by `4`, giving `[0.75, 0.25, 0]` and
  `[0.25, 0.25, 0.5]`; row 1 has a zero total, so the `where=row_totals > 0`
  guard leaves it at the `[0, 0, 0]` the output buffer was initialized with
  rather than producing `nan`. That is also why the row sums come out `1, 0, 1`:
  a class with no test samples contributes nothing, and the plot masks the whole
  row.
- **W&B artifact routing** — the expected split is read straight off the
  `write_bulk_artifacts = not use_wandb` switch: the model pickle, the plaintext
  report and the search plot are written unconditionally, while the results
  JSON, ranked-candidate TSV, marginals TSV and predictions TSV sit behind the
  switch and are logged as `tables/*` payloads plus the `results_json` run
  summary in the other arm. `artifact_paths` is assembled from the same switch,
  so the three unconditional entries are exactly what it holds under
  `use_wandb: true`.
- **Hyperparameter marginals** — four candidates cross a single parameter
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
- **Parameter influence** — the same four scores `0.9, 0.7, 0.6, 0.4` crossed
  over `n_estimators` (10, 5, 10, 5) and `max_depth` (5, 5, 3, 3). Per-value best
  scores are `n_estimators`: 10 → 0.9, 5 → 0.7, a spread of `0.2`; `max_depth`:
  5 → 0.9, 3 → 0.6, a spread of `0.3`. Mean scores give the same ordering
  (`n_estimators`: 0.75 vs 0.55 = 0.2; `max_depth`: 0.8 vs 0.5 = 0.3), so
  `max_depth` sorts first on both. Sub-`1e-12` differences are flattened to an
  exact `0.0` so a dimension whose values all score alike reports zero spread
  rather than floating-point noise, which is what lets the guidance block call it
  out as having no effect.

### `test_each_test_is_corrected_over_every_triplet`

**Inputs:** the five crafted observation sets from the cascade table (one per
outcome), measured under `pipeline_mode="detailed"` with `family_size=5`, no
bootstrap and seed 11, then passed through
`_apply_triplet_result_p_value_correction` under each of `bfn`, `holm` and
`fdr_bh`; then the same five results with every raw DCT p-value replaced by
1.0, passed through the same method.

**Derivation:** the pass corrects the DCT column over all five results and the
KS column over all five results, so each corrected column must equal
`_adjust_p_values` applied directly to the raw column — which is asserted
rather than any particular number, so the check holds for a rank-based method
whose values depend on the whole family as much as for Bonferroni. Every row
carries a `ks_significant` and a `perm_decision` because every test ran. The
five cases were built to land on `DCT`, `THT` and `PERM` gates, so the family
has members the count gate settled and members it did not; the tree-height
column is corrected identically once every count p-value is 1.0, since the KS
family reads nothing of the DCT column, and the cascade then calls every row
`no_introgression` at the first gate.

### `test_efficient_and_detailed_modes_agree_on_every_classification`

**Inputs:** the same five observation sets, measured with `family_size=5`, a
30-iteration bootstrap and seed 11 under each pipeline mode, then passed through
`_apply_triplet_result_p_value_correction` under each of `bfn`, `holm` and
`fdr_bh`.

**Derivation:** the efficient mode skips a test only below a gate that has
already failed on the value the correction can only raise — under `bfn` the
exactly corrected `p × 5`, under `holm`/`fdr_bh` the raw value — and the
permutation p-values are corrected inside each test rather than across
triplets, so the skipped work could neither reach its own cascade nor move
another triplet's numbers. Every field the cascade reads must therefore be
equal between the modes, and the bootstrap votes with them, since each
iteration is judged against the same corrected thresholds. The five cases land
on `DCT`, `THT` and `PERM` gates, which is asserted so that the comparison
exercises the permutation gate rather than only rows where neither mode
resamples.

What may differ is confined to what the efficient mode declined to measure.
Under `holm` and `fdr_bh` the tree-height test is a rank-based family member
and is measured on every row, so the raw and corrected KS columns must agree
everywhere. Under `bfn` the `DCT`-gate row (the 10/10 split, `p = 1.0`) is
left unmeasured and carries `None` in all three KS fields; the four measured
survivors must still be corrected as `min(1, p × 5)` — the family is the
triplet count, and a family shrunk to the four measured would give
`min(1, p × 4)` and lower every survivor's corrected value. On `PERM` rows the
permutation fields agree; on the `DCT` and `THT` rows the efficient result
carries `None` for `perm_decision` and `perm_statistic` with
`perm_note == "direction_test_not_consulted"`, while the detailed result still
records a decision nothing read.

### `test_shifted_statistics_match_an_explicit_shift`

**Inputs:** `x ~ Gamma(2.0, 1.0)` and `y ~ Gamma(2.5, 1.2)` at five `(nx, ny)`
splits — `(10, 30)`, `(30, 10)`, `(7, 7)`, `(4, 25)`, `(2, 60)` — pooled,
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
is algebraically identical to re-running the draw on shifted data — not an
approximation. Requesting the same count under the same seed produces the same
`indices`, since the shifts consume no randomness, so the comparison is
element-wise. Both group orderings and lopsided splits are covered because the
kernel samples whichever group is smaller and recovers the other by subtraction,
which is where a sign or role error would surface. Measured agreement is at
`1e-14` relative; the `1e-9` tolerance is slack.

### `test_equivalence_reuses_the_directional_resamples`

**Inputs:** two `Normal(0, 1)` samples of 40 under seed 8, tested with
`min_resamples = max_resamples = 2000` under seed 2.

**Derivation:** equal means leave neither tail significant, so the equivalence
step runs. Its p-value is `max` of two add-one estimators `(1 + count) / (n + 1)`
over the directional test's own `n_done` draws, so it is an exact multiple of
`1 / (n_resamples + 1)` and cannot fall below that floor. Both properties would
break if the step drew its own resamples at a different count.

### `test_no_bootstrap_skips_the_bootstrap_itself`

**Inputs:** the shared 4-triplet orchestrator fixture with `bootstrap=False`,
under each pipeline mode.

**Derivation:** `bootstrap_stat_ci_low`/`_high` are a percentile interval over
the per-iteration studentized differences, computed inside the bootstrap loop.
If any iteration ran they would be populated, so `None` on every triplet is the
observable proof the loop was skipped rather than merely having its columns
dropped. The classification comes from the point estimate, which the bootstrap
does not feed, so it stays `no_introgression` as in
`test_run_orchestrator_matches_derived_expectation`. The `detailed` row pins
that the mode reinstates only work the cascade could have read, not work the
user switched off.

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

### `test_geometry_matches_subtree_extraction`

**Inputs (parametrized over the six strategies x eleven trees):** the reference
tree with a nested sister pair (`P,Q,R`), a pair spanning the root (`P,R,S`), a
triplet drawn from both sides (`P,S,T`) and one whose odd taxon is listed first
(`S,P,Q`); a five-taxon ladder read at two depths; a tree carrying three taxa
that get pruned away; `(((A,B),C),(D,E));` with no branch lengths at all; one
with a length missing from a single edge; one whose internal branch is exactly
`0.0`; and one at the edge of double precision
(`1e-12` tips under a `1e-13` internal branch against a `1.0000000000001`
sister).

**Derivation:** there is no closed form to compare against here -- the expected
value *is* what extracting the subtree and measuring it produces, which is the
point. `_dendropy_observation` copies the triplet's subtree out with
`extract_triplet_subtree` and measures it with `observation_from_subtree`;
`_geometry_observation` builds the cache and reads the same triplet out of it.
The topology must be equal exactly, because both derive it from discrete
structure rather than arithmetic: extraction from the copied subtree's sister
clade, the cache from which two of the three pairwise LCAs coincide. Heights and
summary metrics compare at `rel=1e-12`, about three orders of magnitude looser
than the largest disagreement measured on real data (one unit in the last place,
`AVG` only). The cases are chosen for what they break rather than for coverage:
the zero-length internal branch would be read as a polytomy by any
depth-comparing rule, the missing lengths must count as `0.0` rather than
propagate `None`, the pruned taxa must not enter any path sum, and the
`1e-12`/`1.0000000000001` case puts the two paths' summation orders as far apart
as the fixture set can.

### `test_geometry_omits_summary_metrics_when_not_collecting`

**Inputs:** the reference tree's `(P,Q,R)` under `AVG`, with
`collect_summary_statistics` false on both paths.

**Derivation:** the third element of an observation carries the per-tree summary
metrics, and a run only needs them under `generate_summary_stats`. Both paths
must return `None` there rather than an empty dict, since the per-triplet
measurement tests that slot for `None` to decide whether to aggregate.

### `test_one_cache_serves_every_triplet_in_the_tree`

**Inputs:** one cache built over the five-taxon reference tree, then all
`C(5,3) = 10` triplets read out of it under `AVG`.

**Derivation:** this is the property the whole design rests on -- the cache is
built per tree, not per triplet, so one cache has to answer every triplet the
tree can supply. The cache costs `O(k^2)` in the pairwise LCA table while a tree
supplies `C(k,3)` triplets, which is what makes building it worthwhile; a
per-triplet cache would cost more than the extraction it replaces. Each triplet
is compared against its own extracted subtree, and the skip decisions are
compared as well (`(actual is None) == (expected is None)`), so a cache that
answered only the triplets it was built from -- or that went stale after the
first read -- fails here.

### `test_geometry_on_a_hand_derived_tree`

**Inputs:** the reference tree, triplets `(P,Q,R)` and `(P,R,S)`, strategy
`AVG`, summary metrics collected.

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

### `test_zero_length_internal_branch_still_resolves`

**Inputs:** `((A:1.0,B:1.0):0.0,C:1.0);`, triplet `(A,B,C)`, strategy `INT`.

**Derivation:** the internal node sits at the same depth as the root, so any
rule that picked the sister pair by comparing LCA *depths* would see a
three-way tie and drop the observation. `find_sister_pair` compares node
identity (`mrca is not root`) and resolves `((A,B),C)`; the cached path compares
LCA node indices and resolves it the same way. `INT` is then the zero-length
edge itself, and `sister_distance = 1.0 + 1.0 - 2(0.0) = 2.0`.

### `test_missing_edge_lengths_count_as_zero`

**Inputs:** `(((A,B),C),(D,E));`, triplet `(A,B,C)`, strategy `AVG`.

**Derivation:** `build_triplet_geometry` stores `0.0` wherever the Newick omits
a length, matching `_distance_to_root`, which skips a `None` edge. Every walk
therefore sums zeros and all four derived values are `0.0`.

### `test_geometry_skips_exactly_what_extraction_skips`

**Inputs:** `(A:1.0,B:1.0,C:1.0);` and `((A:1.0,B:1.0):1.0,D:2.0);`, triplet
`(A,B,C)`.

**Derivation:** in the polytomy all three pairwise LCAs are the root, so the
cached path sees three equal ids and returns `None`; the extraction path reaches
`find_sister_pair`, finds no pair whose MRCA differs from the root, and raises
`ValueError`, which `observation_from_subtree` converts to `None`. In the second
tree `C` is absent, so `leaf_node[C] = -1` on one side and
`set(triplet).issubset(tree_taxa)` fails on the other.

A duplicated taxon label is deliberately not covered: DendroPy raises
`NewickReaderDuplicateTaxonError` while parsing, so neither path can be reached
with one.

### `test_geometry_matches_dendropy_across_a_nine_taxon_tree`

**Inputs:** the nine-taxon tree

```
((((T1:0.11,T2:0.19):0.23,(T3:0.07,T4:0.31):0.0):0.17,((T5:0.29,T6:0.13):0.41,T7:0.53):0.09):0.37,(T8:0.61,OUT:0.71):0.43);
```

and all `C(9,3) = 84` triplets, under each of the six tree-height strategies.

**Derivation:** the expected values are whatever
`extract_triplet_subtree` + `observation_from_subtree` produce, so the test is a
differential one -- it asserts the two implementations agree rather than
restating the arithmetic. The tree is shaped so the sweep covers the cases that
distinguish them: `(T3,T4)` sit above a zero-length internal branch, `T1..T4`
and `T5..T7` sit in sibling clades so many triplets have their sister pair on
one side and the odd taxon on the other, `T7` hangs off a ladder at a different
depth from its clade-mates, and `T8`/`OUT` sit across the root so triplets
drawn from them resolve at the seed node. Every triplet of a nine-taxon rooted
binary tree is resolved, so all 84 must yield an observation on both paths;
`compared == 84` pins that none were silently skipped.

The tolerance is `rel=1e-12`, matching the parity tolerance used elsewhere in
the suite. Measured on real data the two paths agree exactly for the `A`, `B`,
`C`, `SIS` and `INT` strategies and to within one unit in the last place for
`AVG`, so the tolerance is roughly three orders of magnitude looser than the
observed difference.

## tests/orchestrator/test_orchestrator_rename_map.py

### `test_rename_map_reads_a_two_column_tsv` / `test_rename_map_reads_a_yaml_mapping`

**Inputs:** the same two pairs (`T1 -> Homo sapiens`, `T2 -> Pan troglodytes`)
written once as a TSV carrying a `#` comment line and a blank line, and once as
a YAML mapping.

**Derivation:** the loader picks its parser from the file extension, so both
files must yield the identical dict. The comment and blank lines are dropped
before parsing, which is why the TSV's four lines produce two entries.

### `test_rename_map_rejects_malformed_files`

**Inputs:** `T1\tA\textra` (three columns), `T1` (one column), `T1\tA` twice
with different values, `T1\tA` and `T2\tA`, and a YAML list.

**Derivation:** a rename map is a bijection from tree label to display name.
Three columns and one column both fail the two-column requirement. A repeated
label is ambiguous about which name wins. Two labels sharing a name is the case
worth singling out, and is checked in both file formats: it would rename two
distinct taxa to the same string, and since DendroPy raises
`NewickReaderDuplicateTaxonError` on duplicate labels the run would fail later
in the parser with nothing pointing back at the map. A YAML list carries no
keys, so it cannot be a mapping.

### `test_rename_map_rejects_a_missing_file`

**Inputs:** a path that does not exist.

**Derivation:** `FileNotFoundError` rather than `ValueError`, because a mistyped
path is a different mistake from a malformed map and the message names the path
so it can be corrected without opening anything.

### `test_renaming_a_tree_touches_only_mapped_terminals`

**Inputs:** `((T1:0.1,T2:0.2):0.3,T3:0.4);` with a map covering `T1` and `T2`
only.

**Derivation:** the rename is applied per terminal against the map, so the two
mapped taxa become `Alpha` and `Beta` while `T3` keeps its label -- a partial
map is the normal case, since a study usually renames only the taxa it reports
on. The return value counts terminals actually renamed, so it must be `2`, not
the map's size or the tree's terminal count; and the resulting label set
`["Alpha", "Beta", "T3"]` confirms nothing was dropped or duplicated in the
process.

### `test_renaming_labels_leaves_unmapped_names_alone`

**Inputs:** `("T1", "T3")` under `{"T1": "Alpha"}`, and `("T1",)` under `{}`.

**Derivation:** the label helper handles the outgroup and triplet-filter
entries, which arrive as plain strings rather than tree nodes, and it has to
agree with the tree helper or those keys would stop matching the renamed trees.
An empty map is the no-rename case and must return the labels unchanged rather
than an empty list.

### `test_species_rename_map_reaches_every_output`

**Inputs:** the shared 4-triplet orchestrator fixture (taxa `A`, `B`, `C`, `D`,
outgroup `OUT`) with a TSV mapping `A -> Homo` and `B -> Pan`, consolidation
enabled.

**Derivation:** the rename is applied inside `clean_and_save_trees` and
`clean_and_save_gene_trees`, immediately after the Newick is read and before
anything else runs. Everything downstream therefore sees only display names,
which is why the assertions can span outputs written by unrelated code paths:
the triplet tuples (`{Homo, Pan, C, D}` -- `C` and `D` are absent from the map
and so unchanged), the `triplet` and `abc_mapping` columns of the results TSV,
the `species_tree` column whose Newick is rebuilt from the renamed species tree,
the two `processed_*.tree` files, and the consolidation artifacts. The
consolidation check reads `introgression_taxa_order.tsv` because that file holds
the taxon ordering used to label the heatmap axes and the bar chart, so it
standing in display names is the evidence the plot labels do too.

### `test_triplet_resolution_agrees_with_the_observation_guards`

**Inputs:** four hand-picked cases -- a normally resolved triplet, one whose
sister pair sits above a zero-length internal branch, a root polytomy
`(A:1.0,B:1.0,C:1.0);`, and a tree missing taxon `C` -- then every one of the 84
triplets of the nine-taxon tree.

**Derivation:** `geometry_observation` returns `None` for exactly two reasons: a
triplet taxon with no leaf in the cache, and all three pairwise LCAs coinciding.
`triplet_resolution` restates those two guards so preflight can report *which*
one fired, and deliberately shares no code with the hot path, which is kept free
of reason tracking. The expected statuses follow from the guards directly: the
zero-length internal branch is `resolved` because the sister pair is chosen by
LCA node identity rather than depth; the polytomy is `unresolved` because its
three pair LCAs are all the root; the tree without `C` is `missing_taxon`. The
whole-tree sweep then asserts the weaker but essential property -- that
`resolved` and "yields an observation" coincide on every triplet -- so a future
edit to one function that is not mirrored in the other fails here rather than
silently changing what preflight reports.
