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
- `mode` — values are **binned to 3 decimals**, the most frequent bin wins, and
  ties resolve to the **largest** value. With all-distinct values every bin has
  count 1, so the mode is the maximum. For the AVG concordant sample
  `[0.233, 0.243, 0.317, 0.297, 0.250]` the mode is therefore `0.317`.

## tests/orchestrator/test_orchestrator_inference.py

### `test_analyze_triplet_matches_derived_expectation`

**Inputs:** the 10 gene subtrees above, `alpha_dct = alpha_ks = 0.05`, and one of
36 parameter combinations (6 strategies x 3 statistics x 2 discordant tests).
Bootstrap runs with 40 iterations at seed `20240724`.

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
   `dct_significant is False` and the classification is `no_introgression` for
   **every** one of the 36 cases — the decision stops at gate 1.
5. **KS** — `scipy.stats.ks_2samp(concordant_heights, dis1_heights)`. For `AVG`
   the samples are `[0.2333, 0.2433, 0.3167, 0.2967, 0.2500]` and
   `[0.3667, 0.4067, 0.5167]`; they are completely separated, so `D = 1.0` and
   `p ~ 0.0357`. For strategy `C` they overlap
   (`[0.30, 0.32, 0.40, 0.35, 0.31]` vs `[0.30, 0.32, 0.45]`), giving
   `D = 1/3` and `p ~ 0.9643`. For `INT` the values are near-constant, giving
   `D = 0.2`, `p = 1.0`.
6. **Summaries** — apply the chosen statistic to each group. For `AVG` + `mean`:
   concordant `(0.2333+0.2433+0.3167+0.2967+0.2500)/5 = 0.268`; discordant1
   `(0.3667+0.4067+0.5167)/3 = 0.43`.
7. **Classification** — gate 1 fails, so `no_introgression`.

Also asserted: `triplet == ("A","B","C")`, `species_tree == "((A,B),C);"` (the
species subtree serialized topology-only), and that the bootstrap class
fractions sum to 1.

### `test_observation_heights_match_derived_geometry`

**Inputs:** the same 10 subtrees, one strategy per parametrization.

**Derivation:** for each observation, the topology must equal the tabulated
topology and H(T) must equal `_expected_height(entry, strategy)` computed from
the four geometry primitives — i.e. the strategy formulas above applied to the
table, independent of the implementation.

### `test_analyze_triplet_from_observations_matches_newick_path`

**Inputs:** the same triplet, once from `_serialize_triplet_gene_trees` output
and once from raw Newick strings, both with `AVG`/`mean`/`chi-square`.

**Derivation:** both must equal the `_expected_result("AVG","mean","chi-square")`
expectation above. Because the observations and seed are identical, the NumPy
bootstrap draws the same indices, so `bootstrap_value` and `all_bootstrap` must
match exactly (not just approximately).

### `test_bootstrap_is_deterministic_under_seed`

**Inputs:** two identical calls, seed `20240724`, 40 iterations.

**Derivation:** the per-triplet seed is derived deterministically from the run
seed and the triplet, so a repeat call must reproduce the aggregates bit for
bit.

### `test_analyze_triplet_empty_observations`

**Inputs:** an empty gene-subtree list.

**Derivation:** with no observations there are no topology counts, so
`analyzed_trees = 0` and `n_con = 0`; the DCT short-circuits on a zero total to
`(0.0, 1.0)`, which is not significant, so gate 1 returns `no_introgression`.

## tests/orchestrator/test_orchestrator_decision.py

Observations are constructed directly as `(topology, height, None)` tuples, which
lets each test place the triplet on a chosen branch. Bootstrap is disabled
(`iterations: 0`) so results are deterministic.

### `test_classify_no_introgression_when_dct_not_significant`

**Inputs:** 20 concordant observations at height 0.1, 10 dis1 at 0.9, 10 dis2 at
0.9.

**Derivation:** `chisquare([10, 10])` has expected `[10, 10]`, so the statistic
is exactly `0.0` and `p = 1.0`. `1.0 > 0.05`, so gate 1 fails →
`no_introgression`.

### `test_classify_inflow_when_tree_height_test_not_significant`

**Inputs:** 10 concordant at 0.5, 30 dis1 at 0.5, 2 dis2 at 0.5.

**Derivation:** `chisquare([30, 2])` has expected `[16, 16]`, so the statistic is
`(30-16)^2/16 + (2-16)^2/16 = 12.25 + 12.25 = 24.5`, and with df 1
`p ~ 7.4e-07 < 0.05` → gate 1 passes. Every concordant and dis1 height is
identical (0.5), so the two empirical CDFs coincide: `D = 0.0`, `p = 1.0`, not
significant → gate 2 returns `inflow_introgression`.

### `test_classify_outflow_when_concordant_heights_exceed_discordant`

**Inputs:** 10 concordant at 0.9, 30 dis1 at 0.1, 2 dis2 at 0.1.

**Derivation:** the DCT is the same significant 30-vs-2 split. The two height
samples are disjoint with no overlap, so the CDFs separate completely:
`D = 1.0`, and with these sample sizes `p < 0.05` → gate 2 passes. Gate 3
compares `summary_con = 0.9` against `summary_dis = 0.1`; con > dis →
`outflow_introgression`.

### `test_classify_ghost_when_discordant_heights_exceed_concordant`

**Inputs:** the mirror image — 10 concordant at 0.1, 30 dis1 at 0.9, 2 dis2 at
0.9.

**Derivation:** identical DCT and KS reasoning; gate 3 now sees
`summary_con = 0.1 < summary_dis = 0.9` → `ghost_introgression`.

### `test_classify_introgression_truth_table`

**Inputs:** `_classify_introgression(dct_significant, ks_significant,
summary_con, summary_dis)` called directly with 8 explicit rows.

**Derivation:** straight from the decision definition —

| dct_sig | ks_sig | con | dis | Expected | Reason |
| --- | --- | --- | --- | --- | --- |
| False | True | 1.0 | 2.0 | `no_introgression` | gate 1 fails first |
| False | False | 1.0 | 2.0 | `no_introgression` | gate 1 fails first |
| True | False | 1.0 | 2.0 | `inflow_introgression` | gate 2 not significant |
| True | True | 2.0 | 1.0 | `outflow_introgression` | con > dis |
| True | True | 1.0 | 2.0 | `ghost_introgression` | con < dis |
| True | True | 1.0 | 1.0 | `unresolved` | con == dis |
| True | True | None | 1.0 | `unresolved` | missing summary |
| True | True | 1.0 | None | `unresolved` | missing summary |

### `test_adjust_p_values_matches_statsmodels`

**Inputs:** `[0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6]`
and one of the six methods, `alpha = 0.05`.

**Derivation:** `no` must return the input list unchanged. Every other method
must equal `statsmodels.stats.multitest.multipletests(p_values, alpha=0.05,
method=m)[1]` where `m` maps `bfn → bonferroni`, `holm → holm`,
`fdr_bh → fdr_bh`, `fdr_by → fdr_by`, `fdr_tsbh → fdr_tsbh`.

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
normalization is the identity and `triplets == raw_triplets`. Each species
subtree keeps the triplet's taxa and sums the edges along collapsed paths:

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
leaf-label set must be exactly `{A, B, C}` — `D` and `OUT` are dropped.

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
`summary_statistic = mean`).

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
- **KS correction:** asserted as a relation rather than a literal —
  `ks_p_value_corrected == min(1.0, ks_p_value x 4)`. (For (A,B,C) the raw KS
  p-value is `~0.01667`, giving `~0.06667`; for the D triplets dis1 is empty so
  the KS test short-circuits to `1.0` and stays `1.0`.)
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
- `test_parallel_modes_match_serial` — bootstrap seeding is per-triplet and
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
fails and `issues == []`, which makes `passed` `True` and selects the
"No blocking data issues detected" branch of the report.

### `test_report_is_written_to_output_dir`

**Inputs:** the clean pair, with `output_dir` set to a created directory.

**Derivation:** the writer joins `output_dir` with the module constant
`PREFLIGHT_REPORT_FILENAME` (`preflight_data_check.txt`) and writes
`report_text` verbatim, so the file content and `report_text` must be equal and
`report_path` must equal that joined path.

### `test_no_output_dir_skips_writing`

**Inputs:** the clean pair with `output_dir=None`.

**Derivation:** the write branch is guarded on `output_dir is not None`, so
`report_path` stays `None` while `report_text` is still built.

### `test_detects_polytomy_and_missing_outgroup`

**Inputs:** the defective trio, `outgroups=["OUT"]`.

**Derivation:** gene tree 3 contains no `OUT`, so `_root_tree_on_any_outgroup`
returns no used outgroup → one `gene_tree.rooting_failed`, and that tree is
skipped before any triplet check. Trees 1 and 2 root, so
`gene_tree.rooted == 2` while `gene_tree.total_checked == 3`. Tree 2 collapses
A, B and C into a single polytomous clade, so `find_sister_pair` cannot pick a
rooted pair for triplet `A,B,C` → one
`triplet.unresolved_rooted_sister_pair`. The other three triplets each contain
`D`, which sits outside the polytomy, so they still resolve — hence a count of
exactly 1, not 4. The message is formatted with the enumeration index (`Gene
tree #2`) and the comma-joined triplet (`A,B,C`).

### `test_report_attributes_issues_to_gene_trees`

**Inputs:** the defective trio.

**Derivation:** the two categories above start with `gene_tree.` and `triplet.`,
both in `_GENE_CATEGORY_PREFIXES`, and neither starts with `species_tree.`,
`species_triplet.`, or `triplet_filter.`. Summing gives 0 species-tree issues
and 2 gene-tree issues, which the report renders as `NO (count=0)` and
`YES (count=2)` with the fixed column padding shown in the assertions.

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

**Input:** a config dict with `preflight_data_check: True` over the defective
trio.

**Derivation:** `run_orchestrator` prepares the output directory and then
returns `_run_preflight_only(...)` before any tree cleaning, so the only write
into that directory is the report. Listing the directory must therefore yield
exactly `["preflight_data_check.txt"]` — no `metrics.txt`, no processed trees,
no results TSV. `passed` is `False` because the defective trio yields 2 issues.

### `test_runner_preflight_reports_unrootable_species_tree`

**Input:** the same config with `outgroup="NOT_PRESENT"`.

**Derivation:** the `ValueError` raised above is caught in
`_run_preflight_only`, which prints the message and returns `None`, so the
runner returns `None` rather than propagating.

## tests/orchestrator/test_orchestrator_config.py

### `test_cli_defaults_resolve`

**Input:** an `argparse.Namespace` with the three required paths set and every
optional argument `None`.

**Derivation:** each key falls back to its constant in
`ghostparser/orchestrator/config.py`: `alpha_dct`/`alpha_ks` `0.05`,
`p_value_correction` `"bfn"`, `summary_statistic` `"mean"`, `overwrite` `True`,
`discordant_test` `"chi-square"`, `tree_height_calculation_strategy` `"AVG"`,
`min_support_value` `0.5`, `bootstrap_iterations` `100`, `bootstrap_seed`
`None`, and the boolean feature flags `False` — including
`preflight_data_check`, whose default `DEFAULT_PREFLIGHT_DATA_CHECK` is `False`
so that an ordinary run is never turned into a check-only run by accident.
`stats_backend` must be absent entirely, since the custom backend no longer
exists.

### `test_cli_overrides_for_config_plus_cli_options`

**Input:** `alpha_dct=0.01`, `alpha_ks=0.2`, `summary_statistic="mean"`,
`p_value_correction="fdr_bh"`, `no_overwrite=True`.

**Derivation:** each supplied value replaces its default. `no_overwrite=True`
negates to `overwrite is False`.

### `test_config_only_keys_read_from_config_file`

**Input:** a JSON file with `discordant_test: "z-test"`,
`tree_height_calculation_strategy: "SIS"`, `min_support_value: 0.9`,
`generate_summary_stats: true`, `alpha_dct: 0.02`, and
`bootstrap_options: {iterations: 25, seed: 7, debug_mode: true,
summary_only: true}`.

**Derivation:** these keys have no CLI flag, so the file is the only way to set
them. The nested block is flattened onto `bootstrap_iterations = 25`,
`bootstrap_seed = 7`, `bootstrap_debug_mode = True`,
`bootstrap_summary_only = True`.

### `test_config_file_wins_over_cli`

**Input:** a config file setting `alpha_dct: 0.03` and file-specific tree paths,
plus conflicting CLI flags `alpha_dct=0.5` and `summary_statistic="mode"`.

**Derivation:** in config-file mode the file supplies everything, so
`alpha_dct` is `0.03` (not `0.5`). The decisive check is `summary_statistic`:
the file omits it, so it must fall back to the orchestrator default `"mean"` — if
the CLI were consulted it would be `"mode"`. The resolved species path must come
from the file, and a warning naming the ignored flags must be printed.

### Validation tests

- `test_missing_required_field_raises` — `species_tree_path=None` fails
  required-path validation → `ConfigError`.
- `test_parser_exposes_config_file_and_new_flags` — parsing
  `-st s -gt g -og OUT --alpha-dct 0.01 --alpha-ks 0.2 --p-value-correction
  fdr_bh --summary-statistic mean --no-overwrite` must yield those exact values
  with `config_file is None`.

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

### `test_prepare_output_directory_creates_missing_parents`

**Input:** `tmp_path/a/b/results`, which does not exist.

**Derivation:** the helper creates parents on demand, so the whole chain must
exist afterwards.

## tests/test_ml_labels_and_metrics.py

### `test_parse_classes_round_trips_bitstrings`

**Input:** `["100001", "000000", "111111", " 010010 "]`.

**Derivation:** each label is stripped and expanded into six integers, giving a
4x6 matrix. Re-joining each row must reproduce the trimmed label. Set bits per
label: `100001 → 2`, `000000 → 0`, `111111 → 6`, `010010 → 2`, so the matrix sum
is `2 + 0 + 6 + 2 = 10`.

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

## tests/test_introgression_mapper.py — ghost bar colouring

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

### `test_ghost_bars_use_constant_colours_by_sampled_presence`

**Inputs:** the shared scenario, with `matplotlib.axes.Axes.barh` monkeypatched
to record the `color` list and the bar widths before delegating to the original.

**Derivation:** the bar chart is the only `barh` call in the figure, so the
captured `color` list is exactly the ghost bar colours in `taxa_order`. Colour
is a two-valued function of the flag, so the captured set must be a subset of
`{GHOST_ONLY_BAR_COLOR, GHOST_WITH_SAMPLED_BAR_COLOR}`; a per-magnitude colormap
would instead yield a distinct RGBA per bar. `A` maps to
`GHOST_WITH_SAMPLED_BAR_COLOR` and `D` to `GHOST_ONLY_BAR_COLOR` per the table.
The widths for `A` (0.8) and `D` (0.4) differ, which together with ghost-only
taxa sharing one colour demonstrates that magnitude lives in length alone.

### `test_zero_heatmap_cells_are_masked`

**Inputs:** the shared ghost-colour scenario, with `sns.heatmap` monkeypatched
to record the `data` and `mask` it receives.

**Derivation:** the scenario produces exactly one sampled edge, `(C, A)` with
weight 0.5, so in the 4x4 target-by-source matrix over `taxa_order` only the
cell at row `A`, column `C` is non-zero; the diagonal and all 14 remaining
off-diagonal cells are 0. The plotting code passes `mask = heat_values == 0.0`,
so the recorded mask must equal that comparison elementwise, and the count of
unmasked cells must be exactly 1.

### `test_metrics_txt_leads_with_hyperparameters` (random forest)

**Inputs:** `summary_statistics_tsv` with `n_estimators=25`, `random_state=7`,
`test_size=0.25`, `max_depth=None`, `cv_folds=3`,
`rare_class_policy="warn_reduce_cv"`.

**Derivation:** `format_hyperparameter_section` emits the title
`Hyperparameters:` followed by one `  {name:<width}  {value}` line per entry,
where `width` is the longest key (`cv_folds_requested` /
`cv_folds_effective`, 18 characters). The test splits each line on whitespace
rather than asserting column positions, because that width shifts whenever a
key is added. `max_depth` is `None`, which the formatter renders as the literal
`none` so an unset knob is distinguishable from an empty string. The block is
prepended to `text_lines` before the `Test metrics:` group, so its index in the
file is strictly smaller. `cv_folds_requested` echoes the configured `3` while
`cv_folds_effective` carries whatever `auto_cv_folds` resolved for the fixture's
class distribution, which is why only the requested value is asserted.

## Remaining suites

`tests/test_introgression_mapper.py`, `tests/test_ml_config.py`,
`tests/test_ml_utils.py`, `tests/test_ml_random_forest.py`,
`tests/test_ml_multi_knn.py`, and `tests/test_ml_hyper_tune.py` assert either
structural outcomes (files written, columns present, errors raised) or documented
config defaults, and their inputs are the fixtures described in
[TESTS.md](TESTS.md). The two derivations worth stating explicitly:

- **Introgression map averaging** — a directed pair's denominator is every
  triplet in which both taxa co-occur, which for `n` ingroup taxa is `n - 2`
  triplets; a ghost target's denominator is every triplet containing that taxon,
  which is `(n-1) choose 2`. The reported value is the summed `bootstrap_value`
  over classified triplets divided by that denominator.
- **KNN neighbor capping** — `n_neighbors` cannot exceed the number of training
  samples, so the builder clamps it and the metrics report states the effective
  value.
