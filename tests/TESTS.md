# Test Suite Documentation

This is the test map: every test function, the inputs and fixtures it uses, and
the behavior it asserts. For the concrete input values and a step-by-step
derivation of each expected result, see the companion [TEST_IO.md](TEST_IO.md).

## Running Tests

```bash
# Everything
pytest

# One category (see below); combine with -m "core and not integration" etc.
pytest -m core           # the statistics and decisions
pytest -m config         # config loading, resolution and validation
pytest -m output         # what gets written: files, columns, report fields
pytest -m integration    # entry points driven end to end
pytest -m parity         # our geometry vs an independent BioPython reference

# One file / one test
pytest tests/orchestrator/test_orchestrator_inference.py
pytest tests/orchestrator/test_orchestrator_decision.py::test_classify_introgression_truth_table
```

## Categories

Every test carries at least one of five markers, so a change to one part of the
code can be checked with the tests that cover that part rather than the whole
suite. The criteria are meant to be objective:

| Marker | Criterion |
|---|---|
| `core` | The statistics and decisions: the cascade, corrections, permutation and its TOST equivalence step (two one-sided tests against a margin; rejecting both shows the mean difference to be smaller than the margin), geometry, tree preprocessing, metrics math, tuner marginals. Every test with a hand derivation in TEST_IO.md. |
| `config` | Exercises a loader, normalizer or validator; asserts on the resolved config dict or a `ConfigError`. |
| `output` | Asserts on *what was written* (files present, TSV columns, report fields, artifact routing, a figure saved), not on the numbers in them. |
| `integration` | Drives an entry point end to end: `run_orchestrator`, `train_random_forest`, `train_multi_knn`, `tune_hyperparameters`. |
| `parity` | Cross-implementation agreement: the cached geometry against an independent BioPython implementation. |

`config`, `output` and `integration` are applied in the test files, module-wide
with `pytestmark` where a file is homogeneous, per test otherwise. `core` is
added by `tests/conftest.py` to every test that carries none of those three, so
a new logic test needs no mark and `pytest -m core` never silently drops one.
A test may carry several: the consolidation tests that check a computed average
by reading the TSV it was written to are `core` and `output`; the trainer
smoke tests are `integration` and `output`.

`core` is the slowest category (two `holm` cases in the decision tests run
the full deferred bootstrap), and that is the right shape: touching inference
means paying for inference tests. What the split buys is the other direction:
a changed default runs `config` in seconds, a changed report line runs
`output`, and the end-to-end runs wait until commit time.

Config tests do not pin individual default values. CONFIG.md and the sample
configs state those, and two invariants keep them in step with the code
without a per-key test: CLI mode and config-file mode resolve to one and the
same set of defaults, and each full sample (`orchestrator_full.yaml`,
`ml_trainer.yaml`, `hyper_tune.yaml`) names every key at its default.

## Test Philosophy

Expected values are **derived from definitions**, not captured from a reference
implementation: topology counts are read off the input Newick strings by hand,
tree heights are recomputed from root-to-tip distances, and test statistics are
recomputed inline with SciPy/statsmodels. The one intentional exception is the
`@pytest.mark.parity` tests, which exist specifically to check one implementation
against another: the cached geometry against a BioPython implementation written
from the definitions in `tests/orchestrator/tree_references.py`.

## Fixtures In Use

Shared fixtures live in `tests/conftest.py`; the orchestrator's trees live in
`tests/orchestrator/conftest.py`.

- `summary_statistics_tsv`. Inputs: 16 rows, four `class` bitstrings four times
  each, a `dis1_topology` string column and `feature_1..4`. Expected usage: the
  ML trainer and tuner runs (a 25% hold-out keeps one row per class, and three
  CV folds stay feasible).
- `orchestrator_species_tree`. Inputs: the 5-taxon species tree
  `(((A,B),C),(D,OUT))`. Expected usage: orchestrator preprocessing and end-to-end
  runs.
- `orchestrator_gene_trees`. Inputs: 12 gene trees, each containing
  `A,B,C,D,OUT` so rooting always succeeds. Expected usage: orchestrator
  preprocessing and end-to-end runs.

`tests/orchestrator/tree_references.py` holds the BioPython reference, the
observation builders the inference, decision and shape tests share, and the
0.01 threshold those tests pass to every gate, so their derivations, written
at it, hold whatever the run defaults are.

## Function-Level Coverage

### tests/orchestrator/test_orchestrator.py

End-to-end `run_orchestrator` behavior on the shared 5-taxon / 12-gene-tree fixture.

Marked `integration` throughout (every test drives `run_orchestrator`) and `output` where it asserts on files or columns (`run_orchestrator_matches_derived_expectation`, `every_optional_output_lands_where_documented_under_display_names`).

- `test_run_orchestrator_matches_derived_expectation`. Inputs (parametrized,
  2 rows): `run_orchestrator` (serial, fixed bootstrap seed, consolidation off)
  with the defaults, and with `diagnostic` on and the bootstrap off. Expected outputs: all 4 triplets present with the
  hand-derived topology counts (7/3/2 for `(A,B,C)`; 12/0/0 for each
  D-containing triplet), the SciPy chi-square DCT statistic and p-value, a
  Bonferroni-corrected DCT p-value (the raw value times the triplet count),
  `decision_gate == "DCT"` and `no_introgression` for every triplet. With
  `diagnostic` on the KS p-value is likewise corrected by the triplet count,
  the `ks_significant` flag follows it and `perm_decision` is populated; off,
  the raw and corrected KS columns, `ks_significant` and `perm_decision` are
  all `None` and `perm_note` reads `direction_test_not_consulted`.
  With the bootstrap off, `bootstrap_value`, `all_bootstrap` and both
  `bootstrap_perm_stat_ci_*` values are `None` and the first two columns are
  absent from the header. `orchestrator_triplet_results.tsv` has a header
  starting with `triplet` and including `classification`, `perm_p_greater`,
  `perm_p_less`, `decision_gate`, `perm_p_tost` and the two
  `bootstrap_perm_stat_ci_*` columns, and one row per result. Purpose:
  end-to-end correctness against values derived from the fixture rather than
  another module; `diagnostic` changes only what is measured below a settled
  gate; switching the bootstrap off stops its work (the studentized interval
  only the bootstrap loop produces is absent), even in a diagnostic run.
- `test_species_filter_runs_every_triplet_among_the_named_species`. Inputs
  (parametrized, 2 rows): the shared fixture run with a `species_filter` file
  reading `A, B` / `D` / a blank line / `OUT` / `NOPE` / `B`, and one reading
  `A,B,NOPE`. Expected outputs: the first run's results are exactly the one
  triplet `(A, B, D)`, with the hand-derived counts and species subtree of the
  unfiltered run; the second raises `InputError`. Purpose: the filter's names are
  matched against the pruned ingroup (the outgroup and an unknown name are
  skipped, a repeat counts once); every triplet among the survivors and no
  other is run, unchanged in how it is measured, and fewer than three
  survivors is an input error.
- `test_every_optional_output_lands_where_documented_under_display_names`.
  Inputs: one run with `species_rename_map` mapping `A -> Homo sapiens` and
  `B -> Pan sp.`, consolidation, `generate_summary_stats` and
  `bootstrap_diagnostic` on. Expected outputs: the triplet taxa are
  `{Homo sapiens, Pan sp., C, D}`; the `(A,B,C)` species subtree is
  `(('Homo sapiens','Pan sp.'),C);` and every `species_tree` column parses back
  to its triplet's names; the display names appear in the results TSV, the
  consolidation taxa-order file and the inflow/outflow matrix, and the combined
  plot is written; both processed tree files still carry only `A`-`D`/`OUT`;
  `metrics.txt` survives consolidation; `summary_statistics.tsv` has exactly 63
  metric columns (including `discordant2_sister_distance_max`) and every result
  carries `topology_metric_statistics`; the nine per-iteration `bootstrap_*`
  columns are in the results header and every result has `bootstrap_dct_stats`.
  Purpose: every optional output honours its location and column contract, the
  map reaches every named output (quoted where the output is Newick) and
  nothing the run reads back is renamed, all from one run.
- `test_parallel_runs_match_serial`. Inputs: the same fixture run serially
  and with 3 workers under a fixed seed. Expected outputs: every compared
  field, bootstrap values included, is identical. Purpose: distributing
  triplet chunks (three workers over four triplets split unevenly) must not
  change results.

### tests/orchestrator/test_orchestrator_inference.py

Per-triplet inference on a 10-gene-subtree fixture, with expectations recomputed
from the tabulated tree geometry.

All `core`.

- `test_inference_matches_derived_expectation`. Inputs (parametrized over the
  2 discordant tests and `diagnostic` off and on): the 10 gene subtrees serialized
  with the `AVG` strategy, measured by `analyze_triplet_from_observations` and
  decided by the run-wide pass as a family of one under `no`. Expected outputs:
  the counts, DCT/KS statistics and p-values, both significance flags,
  `perm_decision`, and the classification all equal values derived in-test from
  `_LEAF_GEOMETRY` plus direct SciPy/statsmodels calls; bootstrap fractions sum
  to 1. The fixture's C(8, 3) = 56 possible group assignments fall below
  `min_resamples`, so with `diagnostic` on the support guard fires and the
  direction is `inconclusive`; off, the expectation carries `None` for the KS
  fields below a failed count gate and for `perm_decision` below any failed
  gate. Purpose: the measurement and the decision reproduce their definitions
  with the setting off and on; the heights per strategy are pinned separately
  by the geometry test below, so one strategy suffices here.
- `test_observation_heights_match_derived_geometry`. Inputs:
  `serialize_triplet_gene_trees` (the BioPython reference serializer in
  `tests/orchestrator/tree_references.py`) under each of the 6 strategies
  in turn.
  Expected outputs: each observation's topology and H(T) match the
  hand-derived geometry for that strategy. Purpose: pins each tree-height
  strategy to its definition.
- `test_empty_observations_decide_no_introgression`. Inputs: an empty
  observation list, measured and decided. Expected outputs: `analyzed_trees ==
  0`, `n_con == 0`, classification `no_introgression`. Purpose:
  degenerate-input safety.

### tests/orchestrator/test_orchestrator_decision.py

The decision logic and p-value correction, driven with crafted observation sets
because the shared fixture never produces a significant DCT.

All `core` except `test_results_tsv_carries_corrected_columns_only_when_correcting`, which is `output`.

- `test_summary_statistics_discordant_roles_follow_the_counts`. Inputs
  (parametrized, 2 rows): 10 concordant plus 3 and 9 discordant gene subtrees
  serialized with `collect_summary_statistics=True`, run once with `AC|B` the
  more frequent discordant and once with `BC|A`. Expected outputs: `dis1_topology` naming the
  more frequent group, `(n_dis1, n_dis2) == (9, 3)`, and the
  `discordant1_*`/`discordant2_*` means matching the corresponding groups.
  Purpose: the summary columns name the same gene trees as `dis1_topology` and
  the tests, in either direction of the count.
- `test_classify_introgression_truth_table`. Inputs (parametrized, 9 rows):
  every combination of DCT/KS significance and direction, including `equivalent`,
  `inconclusive`, and `None`, put to `_classify_introgression` and, with an
  unmeasured KS flag read as `False`, to the array-form
  `_classification_codes`. Expected outputs: the documented classification
  and the terminating gate (`DCT`/`THT`/`PERM`) for each row, and the same
  label from the array form. Purpose: exhaustive coverage of the cascade,
  which returns the pair, so the gate column cannot drift out of step with
  the classification it explains, and the array form used to tally deferred
  votes cannot drift from the scalar one.
- `test_corrections_match_statsmodels_never_lower_a_value_and_agree_inline`.
  Inputs (parametrized over `P_VALUE_CORRECTION_CHOICES`): a fixed 10-value
  p-value list; a family of eight p-values at 0.001 plus 0.4 and 0.9; and
  p=0.004 padded out to families of 1, 7 and 250. Expected outputs: `no`
  returns the first list unchanged and every other method equals
  `statsmodels.multipletests` with the mapped method name; on the second list
  no adjusted value falls below its raw one; `_adjust_p_value_inline` equals
  the full pass on each padded family under `no` and `bfn`, and raises
  `ValueError` under `holm`, `fdr_bh` and `fdr_by`. Purpose: correction
  correctness against the reference library; the monotonicity that licenses
  skipping a settled gate, asserted over the whole choice list so a
  non-monotone method fails here; and that the inline shortcut is exact.
- `test_unknown_methods_are_rejected`. Inputs (parametrized, 2 rows): an
  unsupported correction method name, and an unsupported discordant-test
  method name. Expected outputs: `ValueError` naming the unsupported method.
  Purpose: input validation.
- `test_degenerate_samples_are_non_significant`. Inputs (parametrized, 4
  rows): `(0, 0)` counts under both discordant tests, and the KS test with one
  empty sample and with two. Expected outputs: `(0.0, 1.0)`. Purpose: the
  zero-discordant and empty-sample short circuits.
- `test_inline_and_deferred_correction_agree_on_a_single_triplet`. Inputs
  (parametrized over `P_VALUE_CORRECTION_CHOICES`): one triplet, family size
  1, 40 bootstrap iterations at a fixed seed, measured and then decided.
  Expected outputs: before the decision pass an inline method has already
  tallied `all_bootstrap` and carries no deferred record, while a rank-based
  one has `all_bootstrap is None` and a 40-element `bootstrap_deferred`; after
  it the record is cleared and `all_bootstrap` and `bootstrap_value` equal the
  `no` run's. Purpose: a family of one leaves each correction as the
  identity, so the inline path (`no`, `bfn`) and the deferred path (the other
  three) must implement one decision rule, and the hand-off from parked votes
  to a tally completes.
- `test_bootstrap_votes_answer_to_the_corrected_threshold`. Inputs: a triplet
  whose raw DCT p-value clears 0.01, analyzed once with `no` and once with `bfn`
  at family size 5000. Expected outputs: the corrected run classifies
  `no_introgression` with `all_bootstrap["no_introgression"] == 1.0`, while the
  uncorrected run's same iterations put it below 1.0. Purpose: the regression
  test for bootstrap votes being judged on raw p-values while the classification
  used corrected ones.
- `test_diagnostic_bootstrap_records_every_test_without_moving_a_vote`.
  Inputs (parametrized, 2 rows): an observation set whose even 10/10
  discordant split fails the count gate in most resamples (`_HIGH` concordant,
  `_LOW[:10]` and `[0.5] * 10` discordant), and the fully separated outflow
  set that clears every gate; each measured under `bfn` as a family of one
  with a 30-iteration bootstrap at seed 11, once with
  `bootstrap_options.diagnostic` off, once on, and once on with
  `summary_only`. Expected outputs: `all_bootstrap` and the
  `bootstrap_perm_stat_ci_*` bounds are equal with the record on and off, and
  the lean run carries no record; the diagnostic run's six recorded lists
  (DCT and KS p-values, permutation statistic, both one-tailed p-values,
  decisions) each hold one entry per iteration, every statistic and p-value
  is present, and every decision is `greater`/`less`/`inconclusive`;
  replaying `_classify_introgression` over the recorded p-values (at
  `alpha = 0.01`, the raw value being the corrected one as a family of one)
  and decisions rebuilds `all_bootstrap` exactly; in summary form the
  decision counts sum to 30 and match the list, and the p-value summary
  counts 30 of 30. Purpose: a diagnostic bootstrap measures all three tests
  in every iteration (the direction test included where the vote never
  reads it), records the vote's own numbers where the vote did, and changes
  no vote, because the extra direction tests draw from their own stream.
- `test_studentized_interval_brackets_the_observed_statistic`. Inputs: 60
  concordant against 40 discordant1 heights on evenly spaced ramps at 200
  iterations, then a degenerate pair of two concordant heights against one
  discordant1 height. Expected outputs: an ordered interval containing
  `perm_statistic` for the first, both bounds `None` for the second. Purpose: the
  percentile interval describes the statistic it accompanies, and is left
  unreported when no resample yields two observations in both groups.
- `test_results_tsv_carries_corrected_columns_only_when_correcting`. Inputs
  (parametrized over `no`/`bfn`/`fdr_bh`): one analyzed triplet written through
  `write_pipeline_results`. Expected outputs: the four `*_corr` columns are
  present for `bfn`/`fdr_bh` and absent for `no`, no header ends in `_no_corr`,
  the raw p-value and significance columns are present in every case, and the
  data row has exactly as many fields as the header. Purpose: the column
  contract, including that the two conditional lists stay in step so a row never
  shifts against its header.
- `test_diagnostic_changes_what_is_measured_and_nothing_concluded`. Inputs
  (parametrized over `bfn`, `holm`, `fdr_bh`): one crafted observation set per
  cascade outcome (an even 10/10 discordant split; a 30/2 split with identical
  con/dis1 heights; fully separated heights in either direction; groups
  sharing a mean but differing in spread), measured with a 30-iteration
  bootstrap as a family of five with `diagnostic` off and again on, each run
  through the run-wide pass. Expected outputs: each set lands on its
  classification at its gate (`no_introgression`/`DCT`,
  `inflow_introgression`/`THT`, `outflow_introgression`, `ghost_introgression`
  and `ambiguous` at `PERM`, with `perm_decision` `greater`, `less`, and
  `equivalent` or `inconclusive`); the counts, `dis1_topology`, every `dct_*` field, `classification`,
  `decision_gate`, `bootstrap_value` and `all_bootstrap` are equal field for
  field between the two runs; the raw and corrected KS fields agree wherever
  the non-diagnostic run measured the test, and are `None` only on `bfn`'s
  `DCT`-gate rows, where under `bfn` the measured survivors' corrected value
  is `min(1, p × 5)`, the triplet count, not the number measured; the
  `perm_*` fields agree on `PERM` rows and are `None` with
  `perm_note == "direction_test_not_consulted"` on the rest, where the
  diagnostic row still carries a decision. On the diagnostic family the
  corrected DCT column and the corrected KS column each equal
  `_adjust_p_values` applied to the whole raw column, and setting every raw
  DCT p-value to 1.0 and re-running the pass leaves the corrected KS column
  identical and classifies every row `no_introgression`. Purpose: the
  every branch of the cascade is reached by its intended route (the gate is
  asserted, so a case cannot pass by coincidence); the setting changes what
  is computed and never what is concluded, across an inline correction and
  two deferred ones, bootstrap votes included; each
  test's family is all triplets, the two corrections read nothing of each
  other, and the cascade order alone decides which one a classification
  rests on.
- `test_bootstrap_measures_the_tree_height_test_its_correction_reads`. Inputs
  (parametrized over `bfn`/`holm` × `diagnostic` off/on): the 10/10-split
  observation set, a 20-iteration bootstrap at seed 11 as a family of one, with
  `run_two_sample_ks_test` wrapped to count its calls, then decided by the
  run-wide pass to read `all_bootstrap`. Expected outputs: fewer than 20
  iterations cleared the count gate (read off `all_bootstrap` as
  `20 × (1 − no_introgression)`); the point estimate carries a KS value except
  under `bfn` with `diagnostic` off; under `bfn` the call count is that point
  estimate (1 or 0) plus the number of iterations that cleared the gate, under
  `holm` it is that point estimate plus every iteration. Purpose: the bootstrap
  measures the tree-height test exactly where its correction reads it (below
  a failed count gate only under a rank-based method), and `diagnostic`
  reaches the point estimate alone, adding at most that one measurement.

### tests/orchestrator/test_orchestrator_permutation.py

The adaptive studentized permutation test, checked against SciPy, against
exhaustive enumeration, and over randomized inputs.

All `core`.

- `test_statistic_p_values_and_verdict_match_scipy`. Inputs (parametrized
  over 5 seeds): random sample pairs at a pinned 4000 resamples under `bfn`.
  Expected outputs: the observed statistic equals
  `scipy.stats.permutation_test`'s and the in-test reference formula exactly;
  both raw one-tailed p-values match SciPy's corresponding single-alternative
  runs within 5 sigma of the binomial standard error of the difference of two
  independent Monte Carlo estimates; and, wherever both of SciPy's p-values
  sit more than that tolerance away from `alpha / 2`, the directional decision
  equals the verdict SciPy's p-values give at `alpha / 2` (the threshold
  Bonferroni over the one-tailed family produces) or is non-directional when
  neither clears it. Purpose: the statistic, the p-value machinery (allowing
  for the two independent resampling streams) and the decision rule all agree
  with the reference.
- `test_permutation_statistics_match_exhaustive_enumeration`. Inputs: a 4-vs-5
  sample pair and 5000 draws from the vectorized sampler. Expected outputs:
  every sampled statistic lies in the exact set obtained by enumerating all
  `C(9, 4) = 126` group assignments through the scalar reference statistic, and
  all 126 appear. Purpose: validates the sample-the-smaller-group-and-subtract
  optimization against ground truth.
- `test_random_inputs_preserve_test_invariants`. Inputs: 12 seeded random
  pairs of random sizes (8-120), spreads, separations, and normal / lognormal /
  exponential families, each run twice under the same seed; then the same 8
  values as both samples. Expected outputs: for every pair a valid decision
  label, p-values in (0, 1], `p_greater + p_less > 1` (both tails count ties),
  the resample count at least the minimum and under twice the budget (the
  last batch is drawn whole), a directional decision
  agreeing with the sign of the statistic, and a second run identical to the
  first; for the identical samples a statistic of 0 and a non-directional
  decision (`equivalent` or `inconclusive`). Purpose: structural invariants
  on shapes no fixed fixture covers, reproducibility under a seed, and that
  identical inputs cannot produce a direction.
- `test_equivalence_step_decides_from_the_directional_resamples`. Inputs
  (parametrized, 3 rows): two samples drawn from one normal distribution at
  n=8 and n=400 per group at a fixed 1000 resamples, and the n=400 pair with
  `equivalence_test=False`. Expected outputs: `inconclusive` at n=8 and
  `equivalent` at n=400, with `p_tost <= 0.01` exactly on the `equivalent`
  row and `p_tost` at or above `1 / (n_resamples + 1)` and an exact multiple
  of it on both; with the step off, `p_tost is None` and the decision falls
  through to `inconclusive`. Purpose: the equivalence margin is an effect
  size, so it shrinks relative to the standard error as data accrues and TOST
  gains power (a margin in standard-error units would report `inconclusive`
  at every n); TOST is an add-one estimator over the directional test's own
  draws, so it answers at that resolution and draws nothing extra; and the
  bootstrap path, which disables the step, pays nothing for it.
- `test_type_one_error_rate_tracks_alpha_under_unequal_variance`. Inputs: 1500
  null replicates, n=60 at sd 1.0 against n=180 at sd 0.3, alpha 0.01. Expected
  outputs: between 3 and 30 *directional* decisions (0.2%-2%; nominal is 15);
  non-directional outcomes are not rejections of the directional null. Purpose: the
  studentization holds the nominal level under unequal sizes and variances,
  which is the regime where an unstudentized permutation test fails.
- `test_guards_short_circuit_without_resampling`. Inputs (parametrized, 4
  rows): one input per guard condition. Expected outputs: the matching
  `note`, an `inconclusive` decision, zero resamples, no statistic, and
  `converged is False`. Purpose: each guard is reachable and inert.
- `test_null_skewness_matches_scipy_and_the_null_shape`. Inputs
  (parametrized, 2 rows): 700 concordant heights against 19 discordant1 of
  which 4 are extreme, at 2500 resamples; and two balanced 150-observation
  samples from one normal family at 4000. Expected outputs: for both,
  `note is None`, `n_resamples_skew` equal to the resample count and
  `null_skew` equal to `scipy.stats.skew` over the same draws; for the first
  `|null_skew| > 1`, a negative statistic and decision `less`; for the second
  `|null_skew| < 0.15`. Purpose: the running power-sum accumulation reproduces
  the reference skewness; an outlier-driven null is measurably asymmetric
  without that disturbing the directional call, and the measure reads near
  zero when the null is symmetric, so a large value means something.
- `test_adaptive_run_converges_or_exhausts_its_budget`. Inputs (parametrized,
  2 rows): a clearly separated pair with `min_resamples=2500`; and a marginal
  0.6 mean shift between two 30-element samples with `min_resamples=100`,
  `max_resamples=1000`. Expected outputs: the first converges after exactly
  one batch of 2500 with decision `greater`; the second has `batches > 1`, a
  total exceeding `batches x 100`, `1000 <= n_resamples < 2000`,
  `converged is False`, and `note == "max_resamples_reached"`. Purpose: an
  easy case stops at the minimum budget instead of spending the ceiling; an
  undecided run escalates its effort and stops once the budget is met,
  drawing the final batch whole rather than trimming it, so the total meets
  or slightly overshoots `max_resamples`, without pinning the growth factor
  (a performance knob).
- `test_bootstrap_resample_budget_shrinks_but_stays_usable_and_monotone`.
  Inputs: six configured `(min, max)` pairs from `(1, 1)` to `(2500, 25000)`,
  then configured budgets from 1 to 25000. Expected outputs: for each pair the
  scaled minimum is at least 1 and, wherever the configured minimum is at
  least 2, strictly below it, and the scaled maximum stays at or above the
  scaled minimum and no higher than the configured one; over the rising
  budgets neither end of the scaled budget ever decreases. Purpose: a
  bootstrap iteration must cost less than the point estimate while still
  leaving an adaptive run a valid range to grow through, and raising the
  configured budget cannot lower the bootstrap's, asserted as properties, so
  the divisor stays a tunable rather than a pinned constant.
- `test_shifted_statistics_match_an_explicit_shift`. Inputs: five group-size
  splits including both orderings, each a gamma-drawn pooled sample drawn
  once with `shifts=(0.75, -0.75)` and once per shift on explicitly shifted
  data under the same seed. Expected outputs: the fused rows equal the
  explicit re-draws to `rel=1e-9`. Purpose: the equivalence test folds its
  shift into the power sums the directional draw already computed; same seed
  means same permutations, so the two must agree to floating point, not
  merely in distribution.

### tests/orchestrator/test_orchestrator_shape.py

The optional distribution-shape diagnostics: moment parity against SciPy, the
modality test's calibration on samples of known modality, and the results-TSV
column contract.

The statistics are `core`; the TSV column-contract test is `output`.

- `test_shape_moments_match_scipy`. Inputs: 500 lognormal draws (seed 4).
  Expected outputs: `skew` and `excess_kurtosis` equal `scipy.stats.skew` and
  `scipy.stats.kurtosis(fisher=True)` called on the same sample, and both are
  positive. Purpose: the moment columns are the SciPy estimators, not a
  re-derivation.
- `test_modality_test_rejects_only_a_well_separated_mixture`. Inputs
  (parametrized over 5 shapes): 400 draws each from a normal, lognormal,
  exponential and gamma(2), plus a 200+200 mixture of normals 4 SD apart.
  Expected outputs: `modes_p > 0.05` for the four unimodal families and
  `<= 0.05` for the mixture. Purpose: the test holds its level on skewed
  unimodal shapes (where a raw mode count reports spurious peaks) while still
  detecting a real mixture.
- `test_tail_index_recovers_known_tail_shapes`. Inputs (parametrized over 3
  shapes): 4000 draws from a Pareto(3), an exponential, and a uniform. Expected
  outputs: `tail_xi` within 0.25 of `1/3`, `0`, and `-1` respectively. Purpose:
  the peaks-over-threshold fit recovers positive, zero, and bounded tails, the
  last of which a Hill estimator cannot represent at all.
- `test_shape_is_not_described_for_small_or_flat_groups`. Inputs: five
  identical values, fifty identical values, and 25 normal draws. Expected
  outputs: every field `None` for the first two; for the third the moments are
  populated but `tail_xi` is `None`. Purpose: the two guards
  (`SHAPE_MIN_OBSERVATIONS`, `SHAPE_MIN_TAIL_EXCEEDANCES`) are independent.
- `test_shape_is_measured_once_and_not_per_bootstrap_iteration`. Inputs: the
  same triplet analyzed with 5 and with 60 bootstrap iterations. Expected
  outputs: identical `con_skew` and `con_modes_p`. Purpose: the diagnostics come
  from the point estimate only, so iteration count cannot move them and the
  modality bootstrap is not paid per iteration.
- `test_shape_columns_reach_the_results_tsv_only_and_only_when_enabled`.
  Inputs (parametrized over `shape_diagnostics` on/off): a 60/40/30 lognormal
  triplet written through `write_pipeline_results` and through
  `write_summary_statistics_tsv`. Expected outputs: `shape_statistics` is
  populated exactly when enabled; in the results TSV the fifteen
  `<group>_<field>` columns are present exactly when enabled, the row aligned
  with the header either way, and with the diagnostics on a positive
  `con_skew`, a `con_modes_p` in `(0, 1]`, and an empty `dis2_tail_xi` (30
  observations leave only 3 in the upper decile); in the summary TSV no column
  ending in any `SHAPE_FIELD_NAMES` suffix appears under either setting, the
  row stays aligned, and the descriptive per-topology columns and
  `classification` are still present. Purpose: the column set is read off the
  results rather than plumbed separately, so it cannot claim a measurement
  that did not happen, and the summary file (consumed as a feature matrix,
  where the diagnostics are undefined below their observation floors) never
  receives them even when measured.

### tests/orchestrator/test_orchestrator_trees.py

Tree preprocessing, asserted against explicit Newick literals.

All `core`, except `test_clean_and_save_trees_filters_on_support_and_quotes_labels`
(`output`).

- `test_clean_and_save_trees_filters_on_support_and_quotes_labels`. Inputs
  (parametrized, 3 rows), each cleaned with `min_avg_support=0.5`: the
  support-free fixture species tree; two supported trees, one averaging 0.97
  and one 0.30; and a four-leaf tree whose labels hold a space, a dot, a quote
  (`'O''Brien'`) and an underscore. Expected outputs, byte for byte: the
  species tree verbatim; only the well-supported tree, its supports stripped;
  the quoted tree verbatim (three labels quoted, the underscore one bare),
  reading back as `Homo sapiens`, `Pan sp.`, `O'Brien`, `Mus_musculus`.
  Purpose: the mean support filter, that a kept tree is written unchanged
  apart from its supports, and that the writer quotes exactly the labels a
  bare token cannot hold, since the run rereads the processed trees.
- `test_root_species_tree_roots_where_the_outgroups_branch_off`. Inputs
  (parametrized, 4 rows): the fixture species tree
  `(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2)` with outgroup `OUT`;
  a tree written with the outgroups on either side of the file's root,
  `(OUT1:0.3,(OUT2:0.2,(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4):0.5)`, and
  one with a root polytomy joining both outgroups and both ingroup clades,
  `(OUT1:0.3,OUT2:0.2,(A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6)`, both with
  outgroups `["OUT2", "OUT1"]`; and a tree without branch lengths,
  `(OUT1,(OUT2,(OUT3,((A,B),(C,D)))))`, with outgroups
  `["OUT3", "OUT2", "OUT1", "OUTX"]`. Expected outputs: ingroup
  `[A, B, C, D]`; the pruned Newick
  `(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;`,
  `(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4;`,
  `((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6);` and `((A,B),(C,D));`; distances
  from the ingroup root, farthest first, `{OUT: 0.5}`,
  `{OUT1: 1.2, OUT2: 0.6}`, `{OUT1: 0.3, OUT2: 0.2}` and none for the tree
  without lengths; and the order that breaks gene-tree rooting ties, `(OUT,)`,
  `(OUT1, OUT2)`, `(OUT1, OUT2)` and `(OUT3, OUT2, OUT1, OUTX)`. Purpose: rooting folds the removed node's edge
  into its sibling; the outgroups need not be a clade as written, only
  branch off at one point; a polytomy at that point keeps every ingroup
  clade under the root; and the outgroups are ranked by path length from
  the ingroup root whatever the listed order, the listed order breaking ties
  and standing in whole when a length is missing, and an outgroup the tree
  lacks ranked last.
- `test_root_species_tree_rejects_outgroups_that_branch_off_twice`.
  Inputs (parametrized): `(((A,B),C),(D,(OUT1,(OUT2,X))))`, where `X` nests
  among the outgroups, and `((OUT1,(A,B)),(OUT2,(C,D)))`, where each outgroup
  carries its own pocket of ingroup taxa, both with outgroups `OUT1`, `OUT2`.
  Expected outputs: `OutgroupRootingError` whose `separated_groups` is
  `(("A","B","C","D"), ("X",))` and `(("A","B"), ("C","D"))` respectively.
  Purpose: outgroups that branch off at more than one point are refused, and
  the error carries the groups the other taxa fall into, largest first.
- `test_generate_triplets_and_species_subtrees`. Inputs (parametrized): the
  pruned species tree, with its branch lengths and with them removed.
  Expected outputs: the 4 sorted triplets from 4 ingroup taxa, no skipped
  triplets, identity ABC normalization, and the exact per-triplet species
  subtree topology Newick strings, the same either way. Purpose: triplet
  enumeration and species-subtree construction, which need no species-tree
  branch lengths.
- `test_read_tree_file_rejects_a_repeated_leaf_label`. Inputs (parametrized,
  2 rows): one tree using leaf label `A` twice; a two-tree file whose second
  tree uses `B` twice. Expected outputs: `InputError` naming the label, and
  `Tree 2` for the second file. Purpose: every geometry and rooting step maps
  a label to one leaf, so a repeated label is refused when the file is read,
  as an input the user can fix.
- `test_read_species_filter_file_collects_names_in_order`. Inputs: a file
  reading ` A, B ` / a blank line / `C` / `B,,D` / a whitespace line / `A`.
  Expected outputs: `["A", "B", "C", "D"]`. Purpose: the species filter's file
  contract; names may share a line or take one each, surrounding whitespace
  and empty entries are dropped, and a repeat keeps its first position.
- `test_clean_and_save_gene_trees_roots_each_tree_from_its_farthest_outgroup`.
  Inputs: `orchestrator_gene_trees` with outgroup `OUT`; then nine gene
  trees with outgroups `["OUT1", "OUT2", "OUT3"]` in species-tree rank:
  `OUT1` alone; `OUT1` with `OUT2` nested among the ingroup; `OUT1` nested
  beside `C` with `OUT2` and `OUT3` sisters; the same with `OUT1` on a
  branch of 9; `OUT2` alone; a tree without branch lengths and one with a
  single internal branch lacking a length; no outgroup; nothing but
  outgroups.
  Expected outputs: all 12 fixture trees survive with no `OUT` leaf, trees 0
  and 3 match their expected rerooted Newick (the ingroup edge absorbs OUT's
  original edge length), `farthest == {"OUT": 12}`, `rooted_on ==
  {"OUT": 12}`, `tangled == {"OUT": 0}` and `tangled_trees == 0`; for the
  nine trees `rooted_count == 7`, `farthest == {"OUT1": 5, "OUT2": 2,
  "OUT3": 0}`, `rooted_on == {"OUT1": 5, "OUT2": 2, "OUT3": 1}`, `tangled ==
  {"OUT1": 1, "OUT2": 4, "OUT3": 1}`, `tangled_trees == 5`,
  `missing_length_indices == [6, 7]`, `unrootable_indices == [8, 9]`, no
  support drops, the seven written trees as exact Newick, and the returned
  `trees` equal to the written lines. Purpose:
  gene-tree rooting semantics; the farthest outgroup is read from the gene
  tree's own branch lengths, a missing one counting as 0, so it can differ
  from the species-tree rank either way, and the rank breaks ties, as in a
  tree without lengths; a tree lacking some length is kept and counted; an
  outgroup among the ingroup taxa once the tree is rooted
  on the farthest is pruned without being used, even when two others sit
  together; the tree is rooted at the common ancestor of the outgroups
  outside the ingroup; every outgroup is pruned; and the counts report, per
  outgroup, the trees it was the farthest in, rooted and was tangled in.

### tests/orchestrator/test_orchestrator_triplet_geometry.py

Covers `orchestrator/triplet_geometry.py`, which reads a triplet's geometry out
of a cached tree instead of extracting its subtree, the path every run takes,
for gene trees and the species tree alike; the cache is built from a
`Bio.Phylo` tree. The parity tests compare it against
`biopython_observation` in `tests/orchestrator/tree_references.py`, written
with `Bio.Phylo` alone: the sister pair is the one pair whose common ancestor
is not the common ancestor of all three, and every distance is a `Bio.Phylo`
path sum from that three-way ancestor.

All `core`; the three that compare against the reference
(`test_geometry_matches_the_reference`,
`test_geometry_skips_exactly_what_the_reference_skips` and
`test_geometry_matches_the_reference_across_a_nine_taxon_tree`) are `parity`
as well.

- `test_geometry_matches_the_reference` - Inputs (parametrized over eleven
  `(newick, triplet)` cases covering nested pairs, pairs spanning the root,
  pruned extra taxa, a ladder, absent and mixed edge lengths, a zero-length
  internal branch, and near-degenerate lengths): the same tree and triplet
  through the reference and through `build_triplet_geometry` +
  `geometry_observation`, under each of the six tree-height strategies.
  Expected outputs: identical topology, and tree height and all three summary
  metrics equal within `rel=1e-12`, for every strategy. Purpose: the cached
  path agrees, across every strategy and tree shape, with an implementation
  that shares no code with it.
- `test_geometry_matches_the_reference_across_a_nine_taxon_tree` - Inputs:
  all 84 triplets of a nine-taxon tree carrying an outgroup, nested clades, a
  ladder, uneven branch lengths and one zero-length internal branch, read out
  of one cache built over the whole tree, against the reference under each of
  the six strategies. Expected outputs: `triplet_resolution` calls every
  triplet resolved, the reference yields an observation for it,
  and topology, height and metrics agree within `rel=1e-12`. Purpose: a
  whole-tree sweep rather than hand-picked shapes, so sister pairs on either
  side of the root and across the zero-length branch are all covered; one build answers every triplet, which is the reuse the
  cache exists for; and the diagnostic agrees with the hot path on every
  triplet.
- `test_geometry_matches_hand_derived_values` - Inputs (parametrized, 4 rows):
  a five-taxon tree with the nested pair `(P, Q, R)` and the root-spanning
  `(P, R, S)` under `AVG`; `((A:1.0,B:1.0):0.0,C:1.0);` under `INT`; and a
  Newick with no lengths under `AVG`. Expected outputs: the stated topology,
  height and all three summary metrics, worked out by hand: 3.0 / 2.0 / 2.0
  and 4.0 / 1.0 / 6.0 for the first two, a zero internal branch with
  `((A,B),C)` resolved for the third, and every value zero for the fourth,
  and a `None` metrics slot when summary statistics are not requested.
  Purpose: pins the arithmetic to a derivation rather than only to the
  reference; the sister pair is chosen by LCA node identity, so a
  zero-length internal branch is not mistaken for a polytomy; a missing
  length counts as zero; and the observation contract leaves the metrics
  empty unless asked for.
- `test_geometry_skips_exactly_what_the_reference_skips` - Inputs
  (parametrized, 3 rows): a root polytomy, a triplet with an absent taxon,
  and a triplet across a zero-length internal branch. Expected outputs:
  `triplet_resolution` returns `unresolved`, `missing taxon` and `resolved`
  respectively, and `geometry_observation` and the reference return an
  observation exactly on the resolved row and `None` on the other two.
  Purpose: the skip decisions match, so observation counts do, and the
  diagnostic the preflight uses restates the hot path's guards without
  sharing code; this is what stops the two drifting apart.

### tests/orchestrator/test_orchestrator_rename_map.py

Covers loading and validating a species rename map, and the two label helpers
the outputs are renamed with. The end-to-end effect on the outputs is covered by
`test_every_optional_output_lands_where_documented_under_display_names` in
`test_orchestrator.py`.

The two map-loading tests are `config`; the renaming test is `core`.

- `test_rename_map_reads_a_tsv_or_a_yaml_mapping`. Inputs (parametrized, 2
  rows): a TSV with a comment line, a blank line, and two entries; the same
  pairs as YAML. Expected outputs: the two-entry mapping from either. Purpose:
  the TSV form, that blanks and comments are ignored, and that the format is
  chosen by extension.
- `test_rename_map_rejects_malformed_or_missing_files`. Inputs (parametrized,
  9 rows): a TSV row with three columns, one with a single column, a repeated
  label, two labels sharing a display name in each of the TSV and YAML forms,
  a YAML list rather than a mapping, a display name holding a tab (YAML
  `"Homo\tsapiens"`) or a comma, and a path that does not exist. Expected
  outputs: `ValueError` naming the problem, or `FileNotFoundError` naming the
  path. Purpose: malformed maps, names the results TSV could not carry, and a
  mistyped path fail at load rather than silently renaming nothing or
  corrupting a column.
- `test_renaming_labels_maps_leaves_and_quotes_as_needed`. Inputs
  (parametrized, 5 rows): three-leaf Newick strings as the run writes them:
  sisters first, odd taxon first with a root edge, labels that contain a map
  key as a substring (`T10`, `XT1`), and an already-quoted label, with the
  map `T1 -> Alpha`, `T2 -> Beta sp.`, and the first string with an empty
  map. Expected outputs: the exact renamed string, with `'Beta sp.'` quoted
  and the substring labels untouched, and the input unchanged under the empty
  map; parsed with Bio.Phylo, the leaves read as the mapped names in the
  input's order with the input's edge lengths, and the plain label helper
  maps the input's leaf list to the same names. Purpose: the `species_tree`
  column is renamed as text by whole leaf token, never inside a branch
  length, and stays valid Newick; the helper used on the results' `triplet`
  and on consolidation's taxon lists agrees with it; and no map is a no-op.

### tests/orchestrator/test_orchestrator_preflight.py

The structural preflight data check and the runner short-circuit that reaches it.

The checks themselves are `core`; `test_clean_inputs_pass_and_the_report_lands_where_documented` is `output`; the runner test is `integration` and `output`.

- `test_clean_inputs_pass_and_the_report_lands_where_documented`. Inputs: a
  5-taxon species tree `(((A,B),C),(D,OUT))` and two well-formed gene trees,
  checked once with an explicit output directory and once with
  `output_dir=None`. Expected outputs: `passed is True`, an empty `issues`
  list, `triplets_checked == 4`, and both gene trees rooted; with a directory,
  `report_path` points at `preflight_data_check.txt` inside it and the file
  content equals `report_text`; without one, `report_path is None` and the
  same report text comes back. Purpose: a clean dataset produces no false
  positives, the report is persisted where documented, and the check is
  usable without touching disk.
- `test_species_tree_is_rooted_where_the_outgroups_branch_off`. Inputs: the
  species tree `(OUT1,(OUT2,(((A,B),C),D)))` with `outgroups=["OUT2", "OUT1"]`,
  and three gene trees: one carrying `OUT1` only, one carrying `OUT1` with
  `OUT2` nested among the ingroup, one carrying `OUT2` only. Expected outputs:
  `passed is True`, `triplets_checked == 4`, three gene trees rooted with
  `gene_tree.farthest.OUT1 == 2`, `gene_tree.farthest.OUT2 == 1`,
  `gene_tree.rooted_on.OUT1 == 2`, `gene_tree.rooted_on.OUT2 == 1`,
  `gene_tree.tangled.OUT2 == 1`, `gene_tree.tangled_trees == 1`,
  `triplet.resolved == 9` and
  `triplet.taxa_absent_from_gene_tree == 3`. Purpose: the
  check roots the species tree with the run's own rooting, so outgroups that
  are not a clade as written but branch off the ingroup at one node still give
  the ingroup `A,B,C,D` and its four triplets; each gene tree roots as the run
  roots it, and in the tree carrying both outgroups its own branch lengths
  make `OUT1` the farthest, which leaves `OUT2` among the ingroup, counted
  rather than reported as a defect.
- `test_detects_polytomy_and_missing_outgroup`. Inputs (parametrized over
  `processes` 1 and 2): gene tree 1 well formed, gene tree 2 a polytomy over
  A/B/C, gene tree 3 with no outgroup label, gene tree 4 without branch
  lengths. Expected outputs: exactly one `gene_tree.rooting_failed` and one
  `triplet.unresolved_rooted_sister_pair`, a
  `gene_tree.missing_branch_lengths` count of 1 with no issue raised for it,
  `missing_length_indices == [4]` and the report section listing `4`,
  three trees rooted out of four checked, the polytomy message naming
  `Gene tree #2` and `A,B,C`, and the pair counters accounting for all 12 triplet/gene-tree pairs as 11 usable, 1
  unresolved, 0 with an absent taxon, with either worker count. Purpose: each
  defect class is detected once and located precisely, the tree lacking a
  length is named by its input position, every pair the check looked at is
  accounted for, and workers change nothing.
- `test_filter_entries_are_validated`. Inputs (parametrized, 2 rows): a
  triplet filter with one valid line, one naming an unknown taxon and one
  naming the outgroup; and a species filter reading `A,B` / `C` / `NOPE` /
  `OUT`, each against the clean fixture. Expected outputs: one
  `<filter>.taxa_missing_in_species_tree`, one `<filter>.includes_outgroup`,
  and `triplets_checked == 1`, the one valid triplet line, or the one
  triplet the three usable species form out of the tree's four. Purpose:
  filter entries are validated rather than silently dropped, and the rest are
  checked.
- `test_impossible_checks_raise`. Inputs (parametrized, 5 rows): an empty
  outgroup list; an outgroup absent from the species tree; every taxon of the
  tree named as an outgroup; outgroups `OUT`, `C`, which branch off the clean
  species tree at two points; and a species-tree file holding two trees.
  Expected outputs: `ValueError` matching "No outgroup taxa were provided",
  "none of the outgroup taxa", "every taxon", "1 taxon: D" and "exactly one
  tree" respectively. Purpose: conditions that make the check impossible
  fail loudly, outgroups that do not root the tree fail with the groups
  named, and the single-tree precondition is enforced.
- `test_runner_preflight_mode_skips_analysis`. Inputs (parametrized, 2
  rows): a config dict with `preflight_data_check: True` and
  `preflight_triplet_cap: 3` against the defective dataset, whose ingroup
  yields four triplets, with outgroup `OUT` and then with an outgroup absent
  from the species tree. Expected outputs: with `OUT` the returned result has
  `passed is False`, `triplets_checked == 3`, an
  `analysis.triplet_cap_applied` issue, and the output directory contains only
  `preflight_data_check.txt`; with the absent outgroup `run_orchestrator`
  raises `InputError` and the output directory is empty. Purpose: the flag
  runs the check and nothing else, the configured cap is what the check
  receives (set below the triplet count so the module default could not pass
  in its place), and an impossible check raises out of the runner.
- `test_index_ranges_collapse_consecutive_runs`. Inputs (parametrized, 3
  rows): `[4]`, `[3, 17, 18, 19, 250]`, `[1, 2, 4, 5, 6]`. Expected outputs:
  `4`, `3, 17-19, 250`, `1-2, 4-6`. Purpose: the report's list of gene trees
  lacking a branch length writes each run of consecutive positions as one
  range and every other position on its own.

### tests/orchestrator/test_orchestrator_config.py

Orchestrator config resolution and config-file precedence. Marked `config`
throughout. No individual default is pinned; two invariants stand in for all of
them.

- `test_cli_and_config_file_share_one_set_of_defaults`. Inputs: the real
  parser given only `-st`/`-gt`/`-og` with the values of the shipped
  `orchestrator_minimal.yaml`, and that file. Expected outputs: the two
  resolved dicts are equal. Purpose: CLI mode and config-file mode fill every
  key from the same constants, pinned without naming any default so a changed
  default never needs a test edit; the minimal sample stays loadable.
- `test_config_only_keys_and_nested_blocks_flatten_from_a_file`. Inputs: a
  file setting `discordant_test`, `tree_height_calculation_strategy`,
  `min_support_value`, `generate_summary_stats`, `alpha_dct`, `seed`, and a
  nested `bootstrap_options` block. Expected outputs: the nine resolved values,
  compared as one dict, with the block flattened to `bootstrap_*`. Purpose:
  config-file-only keys and nested-block flattening.
- `test_p_value_correction_accepts_yaml_bare_word_no`. Inputs (parametrized,
  3 rows): `p_value_correction` written as the bare word `no`, as quoted
  `"no"`, and as `NO`. Expected outputs: each resolves to the `no` choice.
  Purpose: YAML 1.1 resolves bare `no` and `NO` to boolean `False`, so the
  value never reaches validation as a string; the mapping back to the written
  choice is what makes the unquoted spellings work.
- `test_cli_flags_override_the_config_file`. Inputs: a config file setting
  the species tree path, `alpha_dct`, `alpha_ks`, `p_value_correction: no`,
  `overwrite: true` and `consolidation: true`, parsed with the real flag
  strings `-c`, `--alpha-dct`, `--alpha-perm`, `--p-value-correction`,
  `--diagnostic`, `--species-filter`, `--no-overwrite`, `--no-consolidation`,
  `--no-bootstrap`, `--preflight-data-check`, `--preflight-triplet-cap 0`,
  `--processes 3` and `--seed 5`. Expected outputs: each flag's value in its
  key (`triplet_filter` left `None` beside the species filter); `overwrite`,
  `consolidation` and `bootstrap` all `False`; the file's species tree path
  and `alpha_ks` kept. Purpose: the one place the flag strings are pinned;
  a flag beside the file replaces the file's value or adds a key it lacks, the
  negated switches beat the file's canonical `true`, and untouched keys keep
  the file's value.
- `test_outgroup_key_is_normalized_or_rejected`. Inputs (parametrized, 8
  rows), each written to a config file: `outgroup` given as a single label, a comma-separated string, a
  padded string with a trailing comma, and a list whose entries are
  themselves comma-separated; then as `None`, a whitespace string, a list
  holding only an empty string, and an integer. Expected outputs: the first
  four resolve to the flat label list; the rest raise `ConfigError` naming
  `outgroup`. Purpose: one key covers the single- and multiple-outgroup cases
  in every accepted shape, and a value that yields no labels is an error
  rather than an empty outgroup list.
- `test_invalid_values_are_rejected_by_field_name`. Inputs (parametrized, 7
  rows): `species_tree_path: null`, `p_value_correction: true` (what YAML makes
  of a bare `yes`), `diagnostic: "yes"`, `seed: "abc"`,
  `preflight_triplet_cap: -1`, `species_rename_map: absent.tsv`, and
  `triplet_filter` set beside `species_filter`. Expected outputs: `ConfigError`
  whose message names the field (or, for the boolean-typed choice, says `must
  be one of`; for the pair, `cannot both be set`). Purpose: one case per
  validator shape: required path, choice list, optional bool, optional int,
  non-negative int, a path whose file is read when the config resolves, and
  two paths that exclude each other, each failing by name rather than
  surfacing later.
- `test_full_sample_config_names_every_runtime_key_at_its_default`. Inputs:
  `orchestrator_full.yaml` read both as raw YAML and through the loader, plus a
  required-keys-only file. Expected outputs: every key the normalizer produces
  is named in the sample (after allowing for the three renamed path keys and
  the two prefix-flattened blocks), and the sample resolves to the same values
  as the required-keys-only file outside the input paths. Purpose: a new key
  cannot be added without the sample gaining it, and the value the sample shows
  for each key is the one the code would use anyway.
### tests/test_cli.py

The shared command-line wrapper in `ghostparser.cli_config`.

- `test_run_cli_maps_each_outcome_to_its_exit_status`. Inputs (parametrized,
  7 rows): a run returning `None`, returning `EXIT_CHECK_FAILED`, or raising
  `ConfigError`, `InputError`, `OutgroupRootingError`, `TypeError` or
  `KeyboardInterrupt`. Expected outputs: exit statuses 0, 3, 2, 1, 1, 70 and
  130; no traceback on stderr for a failure; and with `--debug` the same
  exception raised out of the wrapper. Purpose: every entry point shares one
  mapping from how a run ends to its exit status, the package's own errors
  are told from bugs, and the debug switch restores the traceback.

### tests/test_config_trunk.py

The shared configuration trunk in `ghostparser.config`.

Marked `config` throughout.

- `test_resolve_path_handles_absolute_relative_and_home`. Inputs: an absolute
  path, a relative path resolved from a chdir'd cwd, and a `~/` path. Expected
  outputs: each resolves to the correct absolute path with no `~` remaining.
  Purpose: path-resolution rules shared by all modules.
- `test_load_raw_config_reads_json_and_yaml_and_rejects_the_rest`. Inputs:
  equivalent payloads as `.json`, `.yaml` and `.yml`; then a missing path, a
  `.txt` file, and a JSON list. Expected outputs: the same mapping from each
  format; then `ConfigError` three times with the documented messages. Purpose: config-file loading and its validation.
- `test_validate_required_path_resolves_or_raises`. Inputs: a present path, then
  absent/empty/whitespace values. Expected outputs: resolution, then
  `ConfigError` for each invalid case. Purpose: required-path validation.
- `test_validate_overwrite_flag_precedence_and_validation`. Inputs: every
  combination of `overwrite`/`no_overwrite`, an explicit default, and
  non-boolean values. Expected outputs: `overwrite` wins over `no_overwrite`,
  negation is applied correctly, and non-booleans raise `ConfigError`. Purpose:
  the shared overwrite semantics.
- `test_prepare_output_directory_overwrites_or_suffixes`. Inputs: an existing
  directory with a stale file, then the same directory with `results_1` and
  `results_3` already taken and `overwrite=False`. Expected outputs: the
  directory is reset in the first case; the second returns `results_2` (the
  smallest missing suffix) and leaves the original intact; a nested path with
  missing parents is created on demand. Purpose: output directory preparation,
  suffix allocation, and parent creation.

### tests/orchestrator/test_orchestrator_consolidation.py

Consolidation outputs and count aggregation.

The count-helper test is `core`; the two that call `generate_introgression_maps` are `output`, and `core` as well since they check computed values through the written TSVs.

- `test_generate_introgression_maps_writes_every_sheet_from_the_results`.
  Inputs: nine result rows on the species tree `(((A,B),C),D)`, with the edge
  `C → B` produced twice (weights 0.6 and 0.2) beside a `no_introgression`
  row over the same taxa, `C → A` once (0.5), `B → C` once (0.3), the ghost
  target `A` twice (0.8 and 0.4), and ghost targets `D` (0.4) and `C` (0.5)
  once each; then the same results written again into an existing directory
  with `overwrite=False`. Expected outputs: `taxa_count == 4`, three sampled
  edges, three ghost targets, the combined PNG and the eight TSVs under
  `consolidation_data/`; in the inflow/outflow sheets `B ← C` averages 0.4
  over a raw sum of 0.8 and a supporting count of 2, `A ← C` 0.5 / 0.5 / 1
  and `C ← B` 0.3 / 0.3 / 1, the raw sums totalling 1.6 and no average above
  0.5; ghost strengths `A: 0.6, B: 0, C: 0.5, D: 0.4` over raw sums
  `1.2, 0, 0.5, 0.4` and supporting counts `2, 0, 1, 1`, with
  `has_sampled_introgression` `1` for `A`, `B` and `C` and `0` for `D`; the
  non-sister sheet holding `{A,C}: 5`, `{B,C}: 6`, `{A,D}: 3`, `{B,D}: 2`,
  `{C,D}: 2` and zero elsewhere; and the second call writing to
  `introgression_1` while the existing directory keeps its stale file.
  Purpose: artifact generation, undiluted averages (each divided by the
  triplets that produced the edge or target, never by the triplets that merely
  contain the taxa), with the raw-sum and supporting-count sheets carrying
  the two parts, the flag that drives ghost bar colour, the non-sister
  counts, and overwrite behavior.
- `test_generate_introgression_maps_selects_the_requested_taxa`. Inputs
  (parametrized, 3 rows): a balanced four-taxon tree with no taxon filter, the
  same tree with `plot_taxa=["A", "B", "C"]`, and an outgroup tree with
  `outgroups=["OG"]`. Expected outputs: `taxa_count` of 4, 3, and 3; the matrix
  header and the ghost rows listing the same taxa; and the excluded taxon absent
  from both. Purpose: `plot_taxa` and `outgroups` decide which taxa reach every
  output, and no output disagrees with another about the set.
- `test_count_helpers_follow_the_classification`. Inputs: five result rows
  covering a classified edge, a `no_introgression` row over the same taxa, an
  unrelated triplet, and two ghost rows on either discordant topology; and a
  taxa order with a `(source, target)` weight map holding one zero-weight
  edge. Expected outputs: supporting counts `{(C, B): 1, (D, B): 1}` and ghost
  counts `{A: 1, C: 1}`; non-sister counts `{A,C}: 3`, `{B,C}: 3`, `{A,D}: 2`,
  `{B,D}: 1`, `{C,D}: 1`; and presence flags `{"A": 1, "B": 0, "C": 1, "D": 0}`.
  Purpose: supporting counts follow the classification, not mere
  co-occurrence; non-sister counts come from every row, the
  `no_introgression` one included; and only taxa that are the target of a
  non-zero sampled edge are flagged.

### tests/test_ml_config.py

Trainer config resolution. Marked `config` throughout; no individual default is
pinned.

- `test_trainer_sample_config_names_every_key_at_its_default`. Inputs:
  `sample_configs/ml_trainer.yaml` read as raw YAML and through the loader,
  plus a file holding only its two required paths. Expected outputs: every
  resolved key is named in the sample (top level, `model` or `evaluation`),
  and the sample resolves to exactly what the paths-only file resolves to.
  Purpose: the one trainer sample loads, covers both trainers' keys, and shows
  each key at the value the code would use anyway.
- `test_ml_config_accepts_explicit_values`. Inputs (parametrized, 14 rows):
  `overwrite: false`, an explicit `target_column` and `seed`; `n_estimators` and
  `min_samples_leaf` under `model`; `model.max_features` as `null`, `log2`, `3`
  and `0.5`; `model.class_weight` as `balanced`, a mapping and a list of
  mappings; and `feature_importance_method` and its correlation threshold
  under `evaluation`. Expected outputs: each value resolves as written, under
  its top-level name. Purpose: precedence over the defaults, the flattening of
  the `model` and `evaluation` blocks, and that every form scikit-learn
  accepts for the two passthrough keys survives the validator.
- `test_ml_config_rejects_invalid_values`. Inputs (parametrized, 17 rows):
  `max_features` as `auto`, the string `None`, `sqrt2`, `0`, `1.5` and `true`;
  `class_weight` as `nope` and `5`; `cv_folds` of `1`, `0`, `null`, `2.5` and
  `true`; `feature_importance_method: shapley`; a threshold of `1.4`; the
  method key at the top level; and `seed` under `model`. Expected outputs:
  `ConfigError` carrying the listed fragments: the dotted location and the
  received value with the accepted forms for the estimator keys (`auto` told
  scikit-learn removed it, the string `None` told to omit the key or write
  `null`), `cv_folds must be an integer >= 2`, the three methods, a fraction
  between 0 and 1, the misplaced key's name, and `top level`. Purpose: a bad
  value fails in config naming the key the user wrote, not as an estimator
  error at fit time; `true` pins that a bool is not taken for the integer 1;
  stratified k-fold cannot split into one fold; and a key in the wrong
  section is refused rather than ignored.
- `test_flags_beside_a_config_file_override_it`. Inputs (parametrized over the
  trainer and the tuner): a config file with `seed: 7`, `overwrite: true` and
  `test_size: 0.3`, parsed by that entry point's own parser with `-c`,
  `--seed 9` and `--no-overwrite`, and resolved by its own resolver. Expected
  outputs: `seed == 9`, `overwrite is False`, `test_size == 0.3`. Purpose:
  both ML entry points wire their two flags through the shared resolver, and
  the flags replace only their own keys.

### tests/test_ml_utils.py

The ML label contract, feature encoding, metrics, distributions, CV-fold
policy and feature importance. All `core`, except the figure-folder test.

- `test_rows_to_matrix_one_hot_encodes_strings_and_drops_the_target`. Inputs:
  two rows with the target, two low-cardinality string columns and two numeric
  ones. Expected outputs: the feature names `topology_a=AC`, `topology_a=BC`,
  `topology_b=XY`, `topology_b=XZ`, `feature_a`, `feature_x`, the matching
  one-hot and numeric values, and the labels in row order. Purpose: the saved
  model is reused on new data, so the column layout (one-hot blocks per column
  in sorted-category order, then the numerics) is a contract.
- `test_rows_to_matrix_rejects_string_features`. Inputs: a string column
  exceeding the cardinality limit. Expected outputs: an error rather than an
  invented ordering.
- `test_row_normalize_confusion_matrix_turns_counts_into_per_class_fractions`.
  Inputs: a 3x3 count matrix `[[3, 1, 0], [0, 0, 0], [1, 1, 2]]` whose middle
  true class has no samples. Expected outputs: cells become
  `[[0.75, 0.25, 0], [0, 0, 0], [0.25, 0.25, 0.5]]`. Purpose: the 64-class
  confusion-matrix plot is drawn on a fixed 0-1 scale, so the normalization has
  to put populated rows on a common scale and leave an empty class at zero (drawn
  at the low end of the scale) instead of dividing by zero.
- `test_per_class_recall_skips_classes_without_hold_out_rows`. Inputs: hold-out
  rows `000000` twice (predicted `000000` and `000001`) and `111111` once
  (predicted correctly). Expected outputs: set-bit counts `[0, 6]` and recalls
  `[0.5, 1.0]`, nothing for the 62 absent classes. Purpose: the per-class
  accuracy figure must not score a class with no hold-out rows as 0, which
  would drag its group down for want of data.
- `test_classes_are_balanced_bounds_the_count_ratio`. Parametrized over
  training/hold-out label sets: every class equally often; ten training
  classes at exactly 1.5 times the others; ten training classes above 1.5
  times; ten hold-out classes at twice the others; a class missing from the
  hold-out set; a class missing from the training set. Expected outputs:
  `True`, `True`, `False`, `False`, `False`, `False`. Purpose: the 1/64 chance
  line assumes equal class weights, so it is drawn only when every partition
  holds all 64 classes and no class outnumbers another by more than 1.5 times.
- `test_save_evaluation_figures_creates_no_folder_without_figures` (`output`).
  Inputs: neither the per-bit nor the 64-class matrix. Expected outputs: an
  empty path mapping and no `figures/` folder. Purpose: a run that reports no
  confusion matrices leaves no empty folder behind.
- `test_correlation_feature_groups_collects_the_redundant_columns`. Inputs: a
  5-row matrix whose second column is the square of the first and whose third
  ranks the rows in an unrelated order, at a threshold of `0.7`. Expected
  outputs: the groups `((0, 1), (2,))`. Purpose: grouping has to catch
  features that rank the rows alike whatever their scale, since that is what
  divides one signal's credit, and leave an unrelated feature on its own.
- `test_grouped_permutation_importance_scores_the_group_that_carries_the_label`.
  Inputs: a 40-row matrix whose first two columns hold one binary signal and
  whose third counts the rows, labels equal to that signal on every bit, and a
  stand-in model predicting from the first column alone. Expected outputs: the
  first group's mean drop exceeds 0.2, the second group's drop and its
  standard deviation are exactly zero. Purpose: a group is scored by the
  micro-F1 it costs to shuffle it, so a group the model does not read must
  cost nothing at all, not merely little.
- `test_feature_importance_rows_refuses_impurity_without_a_tree`. Inputs: a
  model whose estimators expose no `feature_importances_`, asked for `mdi`.
  Expected outputs: `ConfigError` naming a tree-based model. Purpose: a model
  with no impurity to read has to say so rather than fall back to another
  measure, which would leave the reported ranking unexplained.

- `test_64_class_matrix_orders_classes_by_set_bits`. Inputs: three rows with
  true classes `110000`, `000001`, `111111` and predicted classes `000011`,
  `000001`, `111111`. Expected outputs: 64 distinct `class_labels` whose set-bit
  count never decreases along the list and which are lexically sorted within
  each count, `000000` first and `111111` last; the matrix sums to 3 with a 1 at
  each (true, predicted) position looked up by label. Purpose: both axes of the
  64-class matrix run from no bits set to all six, in the same order, and the
  counts are permuted along with the labels: a label-only sort would put every
  off-diagonal count under the wrong pair of names.
- `test_bitstring_labels_are_validated_and_parsed`. Inputs (parametrized, 9
  rows): `000000`, `111111`, `100001`, a whitespace-padded ` 010010 `, and
  the malformed `10000`, `1000010`, `100002`, `10000a` and `""`. Expected
  outputs: `is_valid_bitstring` is `True` for the first four (trimmed) and
  `False` for the rest; a valid label parses to a `1 × 6` binary row that
  re-joins to the trimmed label with as many set bits as it has `1`s; a
  malformed one makes `parse_classes` raise `ValueError` even beside a valid
  neighbour. Purpose: label validation and the bitstring/matrix round trip.
- `test_prediction_metrics_and_rows_match_hand_computation`. Inputs: a 2x6
  true/predicted pair differing in one bit. Expected outputs: exact-match 0.5,
  bitwise 11/12, hamming 1/12, the expected per-bit recall/support; per-row
  matched counts 6 and 5, correct exact-match flags, label strings, and
  per-bit columns. Purpose: metric definitions and the predictions.tsv
  contract.
- `test_distributions_count_labels_and_positives_per_bit`. Inputs: four
  labels with one repeat, and a 2x6 target matrix. Expected outputs: sorted
  labels with correct counts and fractions; per-bit positive counts and
  fractions. Purpose: class and bit distributions.
- `test_build_feature_importance_rows_sorts_descending`. Inputs: three named
  features with scores. Expected outputs: rows ranked most-important first.
  Purpose: importance reporting.
- `test_auto_cv_folds_follows_the_rare_class_policy`. Inputs (parametrized,
  6 rows): label arrays whose smallest class varies, under each
  `rare_class_policy`. Expected outputs: folds kept, reduced to the smallest
  class count, skipped, or raising for a singleton class, each with its
  warning; empty labels return no folds. Purpose: full CV-fold policy
  coverage.
- `test_read_tsv_rows_reads_records_and_rejects_a_header_only_file`. Inputs:
  a two-row TSV with `class`, `feature_1` and `dis1_topology` columns, and a
  header-only TSV. Expected outputs: one dict per data row, feature selection
  on its header keeping order and dropping the target column, and
  `ValueError` for the empty file. Purpose: input reading and feature
  selection.

### tests/test_ml_trainers.py

The smoke tests are `integration` and `output`; the neighbour-capping test is
`core`.

- `test_trainer_runs_end_to_end_and_writes_every_artifact`. Inputs
  (parametrized over `random_forest` and `multi_knn`): `summary_statistics_tsv`
  with a 25% hold-out and three CV folds. Expected outputs: training completes;
  the returned metrics carry both metric tiers, the dataset summary, the
  cross-validation block and the 64-class matrix; the metrics JSON flags
  exact match as diagnostic, carries the bit-label order and the per-stage
  timings; and the predictions TSV (with `matched_label_count`), the model
  pickle and the four figures under `figures/` are written. Purpose: the one
  end-to-end smoke test per trainer, plus its output-file and metrics-field
  contract.
- `test_multi_knn_caps_neighbors_to_the_training_size`. Inputs:
  `n_neighbors=20` against training sets of 2 and of 50. Expected outputs: the
  effective count and the estimator's own `n_neighbors` are both `2` in the
  first case, and `20` passes through untouched in the second. Purpose: asking
  KNN for more neighbors than training samples raises at fit time, so the clamp
  has to reach the estimator and must not fire when it isn't needed.

### tests/test_ml_hyper_tune.py

The loader tests are marked `config`; the `tune_hyperparameters` runs are
`integration` and `output`; the marginal/influence math is `core`; the
run-logger routing test is `output`.

- `test_load_hyper_tune_config_accepts_hyperparameter_tuning_section`. Inputs:
  a tuning section with `model`, `method: random`, `objective: macro_f1` and a
  search space holding a list (`n_estimators`) and a lone value
  (`max_features: log2`). Expected outputs: `model_name`, `search_method` and
  `objective_metric` carry the written values, and the search space keeps the
  shape each entry was written in. Purpose: the rename a caller reads the
  resolved config by, and that candidates passing the `model`-block rules keep
  their container shape.
- `test_load_hyper_tune_config_rejects_malformed_configs`. Inputs
  (parametrized, 7 rows): no `hyperparameter_tuning` section at all;
  `use_wandb: "yes"`; `wandb_detailed_payloads: true` beside `use_wandb: false`;
  an `evaluation` block at the top level; `search_space.max_features:
  [sqrt, auto]`; and `cv_folds` of `1` and of `null`. Expected outputs:
  `ConfigError` naming, respectively, `hyperparameter_tuning`, `must be a
  boolean`, `wandb_detailed_payloads requires`, `evaluation`,
  `search_space.max_features` together with `'auto'`, and `cv_folds must be an
  integer >= 2` for both fold counts. Purpose: one case per structural rule; the section is required
  because nothing about a search can be inferred; `use_wandb` may be omitted
  but not mistyped, and gates the detailed payloads rather than letting them
  be silently ignored; an `evaluation` block is refused rather than ignored,
  since accepting it would let a user believe their settings applied to every
  candidate; a bad candidate is named against its `search_space` key
  before any fit; and every candidate is scored by stratified k-fold, which
  needs at least two folds.
- `test_tuning_sample_config_names_every_key_at_its_default`. Inputs:
  `sample_configs/hyper_tune.yaml` read as raw YAML and through the loader,
  plus a file holding only its required paths and its search space. Expected
  outputs: the sample's top level is exactly the runtime keys plus
  `hyperparameter_tuning`, that block names all nine tuner keys, and the sample
  resolves to exactly what the minimal file resolves to. Purpose: the tuner
  sample loads and shows every key at the value the code would use anyway.
- `test_tune_hyperparameters_runs_end_to_end`. Inputs (parametrized over grid
  search and random search with `n_iter=1`): `summary_statistics_tsv` over a
  space crossing `n_estimators` and `max_depth` at two values each, with the
  `use_wandb` attribute deleted from the namespace. Expected outputs:
  `search_method` echoes the request and `use_wandb` is `False`; exactly 4
  candidates under grid and 1 under random; the best-model pickle is written
  and no `wandb/` directory is; the marginals TSV and the search-report figure
  exist and the plot path is returned; the ranked candidate TSV carries `rank`,
  `is_best`, `cv_score` and `elapsed_seconds`; the marginals header matches
  the documented column order; the text report has its search-space,
  per-parameter, guidance and artifact sections and records W&B as disabled;
  and the results JSON carries `parameter_marginals`, `artifact_paths` and one
  `parameter_influence` entry per search-space key. Purpose: enumeration
  against sampling; the programmatic entry point takes the same W&B default as
  the config loader; and with W&B off the local artifacts are the only way to
  navigate a search, so their file and column contract is pinned.
- `test_tune_hyperparameters_routes_bulk_artifacts_to_wandb`. Inputs:
  `summary_statistics_tsv`, grid search, `use_wandb=True`, with
  `_import_wandb` patched to a stub module. Expected outputs: the output
  directory holds `hyper_tune_best_model.pkl`, `hyper_tune_results.txt` and
  `figures/hyper_tune_search_report.png` but none of `hyper_tune_results.json`,
  `hyper_tune_results.tsv`, `hyper_tune_parameter_marginals.tsv` or
  `predictions.tsv`; the returned paths for those four are `None`;
  `artifact_paths` names only the three written files; the stub received the
  `tables/ranked_candidates`, `tables/parameter_marginals` and
  `tables/predictions` log payloads; the run summary carries
  `bulk_artifacts_written_locally` false and a `results_json` blob that parses
  back to the full payload; and the text report records W&B as enabled and
  points at the run for the rest. Purpose: turning W&B on has to move the bulk
  outputs rather than duplicate them, and nothing may be silently dropped.
- `test_compute_parameter_marginals_summarizes_each_value`. Inputs
  (parametrized over `max` and `min`): four candidates over
  `n_estimators` in `{5, 10}` with scores `0.7/0.4` and `0.9/0.6`, ranked in the
  objective's direction. Expected outputs: values are ranked against each other
  in that direction, the winner's `best_score` is `0.9` under `max` and `0.4`
  under `min`, each value counts 2 candidates and reaches `best_rank` 1, and the
  winner is flagged `upper_bound` / `lower_bound` respectively while the runner-up
  carries no flag. Purpose: the marginal summary drives the guidance block, so
  the direction handling and the bound flag are checked in both directions.
- `test_parameter_influence_ranks_by_best_score_spread`. Inputs: four candidates
  crossing `n_estimators` and `max_depth` with scores `0.9/0.7/0.6/0.4`. Expected
  outputs: `max_depth` ranks first with a best- and mean-score spread of `0.3`,
  ahead of `n_estimators` at `0.2`. Purpose: the influence ordering tells a user
  which dimension to keep searching.
- `test_compute_parameter_marginals_summarizes_each_value`. Inputs
  (parametrized over `max` and `min`): four candidates over
  `n_estimators` in `{5, 10}` with scores `0.7/0.4` and `0.9/0.6`, ranked in the
  objective's direction. Expected outputs: values are ranked against each other
  in that direction, the winner's `best_score` is `0.9` under `max` and `0.4`
  under `min`, each value counts 2 candidates and reaches `best_rank` 1, and the
  winner is flagged `upper_bound` / `lower_bound` respectively while the runner-up
  carries no flag. Purpose: the marginal summary drives the guidance block, so
  the direction handling and the bound flag are checked in both directions.
- `test_parameter_influence_ranks_by_best_score_spread`. Inputs: four candidates
  crossing `n_estimators` and `max_depth` with scores `0.9/0.7/0.6/0.4`. Expected
  outputs: `max_depth` ranks first with a best- and mean-score spread of `0.3`,
  ahead of `n_estimators` at `0.2`. Purpose: the influence ordering tells a user
  which dimension to keep searching.
- `test_create_run_logger_routes_on_the_use_wandb_choice`. Inputs
  (parametrized): `use_wandb` false and true, with `_import_wandb` patched to a
  stub module. Expected outputs: with `false` the logger reports
  `enabled is False`, forwards nothing, and creates no `wandb/` directory; with
  `true` it forwards the log payload, writes the run summary, finishes the run,
  and creates the directory. Purpose: one config flag is the only thing that
  decides whether the run touches W&B.
- `test_import_wandb_raises_config_error_when_unavailable`. Inputs: `wandb`
  patched to raise `ModuleNotFoundError` on import. Expected outputs:
  `ConfigError` whose message names `use_wandb: false`. Purpose: a missing
  optional dependency surfaces as an actionable config error rather than an
  import traceback.

## Parity Tests (`@pytest.mark.parity`)

Run with `pytest -m parity`. Only the three geometry tests in
`tests/orchestrator/test_orchestrator_triplet_geometry.py` that compare the cached
path against the BioPython reference carry this marker; every
other test derives its expectations from definitions instead of comparing
against a second implementation.
