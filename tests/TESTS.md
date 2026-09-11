# Test Suite Documentation

This is the test map: every test function, the inputs and fixtures it uses, and
the behavior it asserts. For the concrete input values and a step-by-step
derivation of each expected result, see the companion [TEST_IO.md](TEST_IO.md).

## Running Tests

```bash
# Everything
pytest

# One file / one test
pytest tests/orchestrator/test_orchestrator_inference.py
pytest tests/orchestrator/test_orchestrator_decision.py::test_classify_ghost_when_discordant_heights_exceed_concordant

# Only the cross-library parity tests (DendroPy vs BioPython)
pytest -m parity
```

On the SLURM cluster, route the suite onto a compute node rather than running it
on a login node.

## Test Philosophy

Expected values are **derived from definitions**, not captured from a reference
implementation: topology counts are read off the input Newick strings by hand,
tree heights are recomputed from root-to-tip distances, and test statistics are
recomputed inline with SciPy/statsmodels. The one intentional exception is the
`@pytest.mark.parity` pair, which exists specifically to check that DendroPy's
triplet extraction agrees with an independent BioPython implementation.

## Fixtures In Use

Shared fixtures live in `tests/fixtures.py` and are re-exported to the whole
suite by `tests/conftest.py`. Orchestrator-specific fixtures live in
`tests/orchestrator/conftest.py`.

- `simple_newick_file` — Inputs: one rooted Newick tree containing `OutGroup`.
  Expected usage: basic tree reading/standardization.
- `newick_with_support_file` — Inputs: a Newick tree with internal support
  labels. Expected usage: support parsing and stripping.
- `multiple_trees_file` — Inputs: three Newick trees, one per line. Expected
  usage: multi-tree file reading.
- `low_support_tree_file` — Inputs: two trees, one with supports
  (0.95, 0.99, 0.98) and one with (0.30, 0.20, 0.40). Expected usage: mean
  support filtering at a 0.5 threshold.
- `simple_species_tree` / `simple_gene_trees` — Inputs: a minimal species tree
  and three gene trees. Expected usage: small parser integration checks.
- `triplet_comparison_cases` — Inputs: three `(newick, triplet)` pairs with
  known branch lengths. Expected usage: the DendroPy-vs-BioPython parity tests.
- `summary_statistics_tsv` / `summary_statistics_tsv_tuning` — Inputs: feature
  tables with a `class` bitstring column, `dis1_topology`, and `feature_1..4`.
  Expected usage: ML trainer and tuner tests.
- `orchestrator_species_tree` — Inputs: the 5-taxon species tree
  `(((A,B),C),(D,OUT))`. Expected usage: orchestrator preprocessing and end-to-end
  runs.
- `orchestrator_gene_trees` — Inputs: 12 gene trees, each containing
  `A,B,C,D,OUT` so rooting always succeeds. Expected usage: orchestrator
  preprocessing and end-to-end runs.

## Function-Level Coverage

### tests/orchestrator/test_orchestrator.py

End-to-end `run_orchestrator` behavior on the shared 5-taxon / 12-gene-tree fixture.

- `test_run_orchestrator_matches_derived_expectation` — Inputs: `run_orchestrator`
  (serial, `taxon` mode, fixed bootstrap seed, consolidation off). Expected
  outputs: all 4 triplets present with the hand-derived topology counts
  (7/3/2 for `(A,B,C)`; 12/0/0 for each D-containing triplet), the SciPy
  chi-square DCT statistic and p-value, Bonferroni-corrected p-values, and
  `no_introgression` for every triplet. Purpose: end-to-end correctness against
  values derived from the fixture rather than another module.
- `test_run_orchestrator_writes_results_tsv` — Inputs: the same serial run.
  Expected outputs: `orchestrator_triplet_results.tsv` exists, its header starts
  with `triplet` and includes `classification`, `bootstrap_value`,
  `perm_p_greater`, `perm_p_less`, and `decision_gate`, and its row count
  equals the number of results. Purpose: TSV shape and naming.
- `test_no_bootstrap_omits_the_bootstrap_columns` — Inputs: a run with
  `bootstrap=False`. Expected outputs: `bootstrap_value` and `all_bootstrap` are
  absent from the header while `classification` remains, and all 4 triplets are
  still produced. Purpose: the bootstrap toggle only removes bootstrap output.
- `test_no_bootstrap_skips_the_bootstrap_itself` — Inputs (parametrized over both
  pipeline modes): a run with `bootstrap=False`. Expected outputs: every result
  has `bootstrap_value`, `all_bootstrap`, `bootstrap_stat_ci_low` and
  `bootstrap_stat_ci_high` at `None`, and the classification is unchanged.
  Purpose: the toggle stops the work rather than only the columns, and the
  detailed mode does not reinstate it — the studentized interval comes from the
  bootstrap loop, so its absence is the observable proof no iterations ran.
- `test_consolidation_preserves_run_outputs` — Inputs: a run with
  `consolidation=True`. Expected outputs: the results TSV, `metrics.txt`, and
  both processed tree files survive, and a non-empty `consolidation/` subfolder
  exists. Purpose: regression guard against consolidation wiping the run folder.
- `test_generate_summary_stats_writes_tsv` — Inputs: a run with
  `generate_summary_stats=True`. Expected outputs: `summary_statistics.tsv`
  exists with exactly 63 topology/metric columns (including
  `concordant_avg_tree_height_mean` and `discordant2_sister_distance_max`), and
  at least one result carries populated `topology_metric_statistics`. Purpose:
  the summary-statistics feature and its column contract.
- `test_bootstrap_debug_mode_writes_debug_columns` — Inputs: a run with
  `bootstrap_debug_mode=True`. Expected outputs: the debug columns
  (`bootstrap_dct_stats`, `bootstrap_dct_p_value`, `bootstrap_ks_stats`,
  `bootstrap_ks_p_value`, `bootstrap_gene_tree_heights`) appear in the header
  and at least one result has a populated `bootstrap_dct_stats`. Purpose: the
  bootstrap-debug output path.
- `test_species_rename_map_reaches_every_output` — Inputs: the shared fixture
  run with `species_rename_map` mapping `A -> Homo` and `B -> Pan`, with
  consolidation on. Expected outputs: the triplet taxa are `{Homo, Pan, C, D}`,
  the display names appear in the results TSV, the per-triplet species subtree,
  both processed tree files, the consolidation taxa-order file and the inflow
  /outflow matrix, and the combined plot is written. Purpose: renaming at read
  time reaches every named output, and unmapped taxa are left alone.
- `test_parallel_runs_match_serial` — Inputs (parametrized over 2 and 4
  workers): the same fixture run serially and in parallel with a fixed seed.
  Expected outputs: every compared field, bootstrap values included, is
  identical. Purpose: distributing triplet chunks must not change results, at
  either worker count.

Also in `test_orchestrator.py`:


### tests/orchestrator/test_orchestrator_inference.py

Per-triplet inference on a 10-gene-subtree fixture, with expectations recomputed
from the tabulated tree geometry.

- `test_analyze_triplet_matches_derived_expectation` — Inputs (parametrized over
  6 tree-height strategies x permutation on/off x 2 discordant tests = 24
  cases): `analyze_triplet` over the 10 gene subtrees. Expected outputs: the
  counts, DCT/KS statistics, mean/median con/dis values, `perm_decision`, and
  classification all equal values derived in-test from `_LEAF_GEOMETRY` plus
  direct SciPy/statsmodels calls; bootstrap fractions sum to 1. With the
  permutation test on, the fixture's C(8, 3) = 56 possible group assignments
  fall below `min_resamples`, so the support guard fires and the direction is
  `ambiguous`; with it off, the median comparison decides. Purpose: the full
  inference surface across every parameter combination, on both direction paths.
- `test_observation_heights_match_derived_geometry` — Inputs (parametrized over
  the 6 strategies): `_serialize_triplet_gene_trees`. Expected outputs: each
  observation's topology and H(T) match the hand-derived geometry for that
  strategy. Purpose: pins each tree-height strategy to its definition.
- `test_analyze_triplet_from_observations_matches_newick_path` — Inputs: the same
  triplet analyzed from precomputed observations and from Newick strings.
  Expected outputs: both equal the derived expectation and share identical
  bootstrap aggregates. Purpose: the observation fast path is equivalent to the
  Newick path.
- `test_bootstrap_is_deterministic_under_seed` — Inputs: two identical
  `analyze_triplet` calls with the same seed. Expected outputs: identical
  `all_bootstrap` and `bootstrap_value`. Purpose: seeded reproducibility.
- `test_analyze_triplet_empty_observations` — Inputs: an empty gene-subtree list.
  Expected outputs: `analyzed_trees == 0`, `n_con == 0`, classification
  `no_introgression`. Purpose: degenerate-input safety.

### tests/orchestrator/test_orchestrator_decision.py

The decision logic and p-value correction, driven with crafted observation sets
because the shared fixture never produces a significant DCT.

- `test_decision_cascade_lands_on_each_classification` — Inputs (parametrized, 5
  rows): crafted observation sets, one per outcome — an even 10/10 discordant
  split; a significant 30/2 split with identical con/dis1 heights; the same split
  with fully separated spread-out heights in either direction; and a split whose
  groups share a mean but differ sharply in spread. Expected outputs: the DCT and
  KS significance flags, the classification, `decision_gate`, and (where a
  direction is expected) `perm_decision`, with `perm_decision` non-None
  throughout. Purpose: every branch of the cascade is reached by its intended
  route, and the gate column shows which test settled each call rather than
  letting a case pass by coincidence.
- `test_permutation_guards_surface_on_the_triplet_result` — Inputs (parametrized,
  2 rows): 4 concordant against 2 discordant1 heights (C(6, 2) = 15 assignments),
  and internally constant groups with different means (`[0.9] * 10` versus
  `[0.1] * 30`). Expected outputs: the matching `perm_note`, zero resamples, and
  `perm_decision == "inconclusive"`. Purpose: a guarded direction test reports
  its reason on the triplet result instead of a direction.
- `test_summary_statistics_discordant_roles_follow_the_counts` — Inputs
  (parametrized, 2 rows): 10 concordant plus 3 and 9 discordant gene subtrees
  with `collect_summary_statistics=True`, run once with `AC|B` the more frequent
  discordant and once with `BC|A`. Expected outputs: `dis1_topology` naming the
  more frequent group, `(n_dis1, n_dis2) == (9, 3)`, and the
  `discordant1_*`/`discordant2_*` means matching the corresponding groups.
  Purpose: the summary columns name the same gene trees as `dis1_topology` and
  the tests, in either direction of the count.
- `test_classify_introgression_truth_table` — Inputs (parametrized, 9 rows):
  every combination of DCT/KS significance and direction, including `equivalent`,
  `inconclusive`, and `None`. Expected outputs: the documented classification and the terminating gate
  (`DCT`/`THT`/`PERM`) for each row. Purpose: exhaustive coverage of
  `_classify_introgression`, which returns the pair, so the gate column cannot
  drift out of step with the classification it explains.
- `test_inline_and_deferred_correction_agree_on_a_single_triplet` — Inputs
  (parametrized over `no`/`bfn`/`holm`/`fdr_bh`/`fdr_by`): one triplet, family
  size 1, 40 bootstrap iterations at a fixed seed. Expected outputs: identical
  `all_bootstrap` and `bootstrap_value` for every method. Purpose: a family of
  one leaves each correction as the identity, so the inline path (`no`, `bfn`)
  and the deferred path (the other three) must implement one decision rule.
- `test_bootstrap_votes_answer_to_the_corrected_threshold` — Inputs: a triplet
  whose raw DCT p-value clears 0.05, analyzed once with `no` and once with `bfn`
  at family size 5000. Expected outputs: the corrected run classifies
  `no_introgression` with `all_bootstrap["no_introgression"] == 1.0`, while the
  uncorrected run's same iterations put it below 1.0. Purpose: the regression
  test for bootstrap votes being judged on raw p-values while the classification
  used corrected ones.
- `test_deferred_bootstrap_record_is_cleared_after_correction` — Inputs: one
  triplet under `holm` with 20 iterations. Expected outputs: `all_bootstrap is
  None` and a 20-element `bootstrap_deferred` before correction; after
  `_apply_triplet_result_p_value_correction`, `bootstrap_deferred is None` and
  the fractions sum to 1. Purpose: the deferral hand-off completes.
- `test_vectorized_bootstrap_codes_match_classify_introgression` — Inputs
  (parametrized, 16 rows): every DCT/KS/direction combination. Expected outputs:
  `_classification_codes` maps to the same label `_classify_introgression`
  returns. Purpose: the array-form cascade used to tally deferred votes cannot
  drift from the scalar one.
- `test_every_supported_correction_is_monotone` — Inputs (parametrized over
  `P_VALUE_CORRECTION_CHOICES`): a family of eight p-values at 0.001 plus 0.4 and
  0.9. Expected outputs: no adjusted value falls below its raw one, for every
  method. Purpose: asserts the property that licenses both short-circuits — the
  point estimate's and the bootstrap's — over the whole choice list rather than a
  fixed set of names, so adding a non-monotone method fails here immediately.
- `test_efficient_and_detailed_modes_agree_on_every_classification` — Inputs
  (parametrized over `bfn` and `holm`): one observation set per cascade outcome,
  analyzed in both pipeline modes, then run through the
  run-wide correction pass as one family. Expected outputs: the two modes agree
  field for field on counts, DCT/KS values, `classification`, `decision_gate`,
  `bootstrap_value`, and `all_bootstrap`; the efficient run has an empty `perm_*`
  block with `perm_note = direction_test_not_consulted` on exactly the rows whose
  gate is not `PERM`. Purpose: the efficient mode changes what is computed, never
  what is concluded, across an inline correction and a deferred one.
- `test_tree_height_family_holds_only_the_count_gate_survivors` — Inputs
  (parametrized over both pipeline modes): three triplets under `bfn`, two
  settled by the count gate and one not. Expected outputs: the settled two report
  no `ks_p_value_corrected` and no `ks_significant` and classify
  `no_introgression`; the survivor's corrected value equals its raw one, since it
  is the only family member. Purpose: pins that the tree-height correction family
  is the count-gate survivors, and that the set does not depend on the mode.
- `test_inline_bonferroni_matches_the_family_correction` — Inputs (parametrized
  over family sizes 1, 7, 250): p=0.004 padded out to that family. Expected
  outputs: `_adjust_p_value_inline` equals the full `_adjust_p_values` pass on
  the same family. Purpose: the inline shortcut is exact, not an approximation.
- `test_inline_correction_rejects_a_rank_based_method` — Inputs: `holm` passed to
  the inline corrector. Expected outputs: `ValueError`. Purpose: rank-based
  methods cannot be applied without the whole family.
- `test_studentized_interval_brackets_the_observed_statistic` — Inputs: 60
  concordant against 40 discordant1 heights on evenly spaced ramps at 200
  iterations, then a degenerate pair of two concordant heights against one
  discordant1 height. Expected outputs: an ordered interval containing
  `perm_statistic` for the first, both bounds `None` for the second. Purpose: the
  percentile interval describes the statistic it accompanies, and is left
  unreported when no resample yields two observations in both groups.
- `test_results_tsv_carries_corrected_columns_only_when_correcting` — Inputs
  (parametrized over `no`/`bfn`/`fdr_bh`): one analyzed triplet written through
  `write_pipeline_results`. Expected outputs: the four `*_corr` columns are
  present for `bfn`/`fdr_bh` and absent for `no`, no header ends in `_no_corr`,
  the raw p-value and significance columns are present in every case, and the
  data row has exactly as many fields as the header. Purpose: the column
  contract, including that the two conditional lists stay in step so a row never
  shifts against its header.
- `test_adjust_p_values_matches_statsmodels` — Inputs (parametrized over all 6
  correction methods): a fixed 10-value p-value list. Expected outputs: `no`
  returns the input unchanged; every other method equals
  `statsmodels.multipletests` called directly with the mapped method name.
  Purpose: correction correctness against the reference library.
- `test_adjust_p_values_bonferroni_by_definition` — Inputs: `[0.01, 0.2, 0.5]`
  with `bfn`. Expected outputs: `[0.03, 0.6, 1.0]` — multiply by 3, clamp at 1.
  Purpose: pins Bonferroni to its arithmetic definition.
- `test_adjust_p_values_rejects_unknown_method` — Inputs: an unsupported method
  name. Expected outputs: `ValueError`. Purpose: input validation.
- `test_discordant_count_test_with_no_discordant_observations` — Inputs
  (parametrized over both tests): `(0, 0)` counts. Expected outputs:
  `(0.0, 1.0)`. Purpose: the zero-discordant short circuit.
- `test_discordant_count_test_rejects_unknown_method` — Inputs: an unsupported
  method. Expected outputs: `ValueError`. Purpose: input validation.
- `test_ks_test_with_an_empty_sample` — Inputs (parametrized over three
  empty/non-empty combinations). Expected outputs: `(0.0, 1.0)`. Purpose: the
  empty-sample short circuit.

### tests/orchestrator/test_orchestrator_permutation.py

The adaptive studentized permutation test, checked against SciPy, against
exhaustive enumeration, and over randomized inputs.

- `test_observed_statistic_matches_scipy` — Inputs (parametrized over 5 seeds):
  random sample pairs. Expected outputs: the observed statistic equals
  `scipy.stats.permutation_test`'s and the in-test reference formula exactly.
  Purpose: the statistic itself is the Welch-studentized mean difference.
- `test_p_values_match_scipy_within_monte_carlo_error` — Inputs (parametrized
  over 8 seeds): random sample pairs at a pinned 4000 resamples with correction
  off. Expected outputs: both one-tailed p-values match SciPy's corresponding
  single-alternative runs within 5 sigma of the binomial standard error of the
  difference of two independent Monte Carlo estimates. Purpose: parity of the
  p-value machinery, allowing for the two independent resampling streams.
- `test_decision_matches_scipy_directional_verdict` — Inputs (parametrized over
  3 seeds): a separated sample pair. Expected outputs: the directional decision
  equals the verdict from SciPy's two one-tailed p-values at `alpha / 2`, which
  is the threshold Bonferroni over the one-tailed family produces. Purpose: the
  decision rule, not just the p-values, agrees with the reference.
- `test_permutation_statistics_match_exhaustive_enumeration` — Inputs: a 4-vs-5
  sample pair and 5000 draws from the vectorized sampler. Expected outputs:
  every sampled statistic lies in the exact set obtained by enumerating all
  `C(9, 4) = 126` group assignments through the scalar reference statistic, and
  all 126 appear. Purpose: validates the sample-the-smaller-group-and-subtract
  optimization against ground truth.
- `test_random_inputs_preserve_test_invariants` — Inputs (parametrized over 12
  seeds): random sizes (8-120), spreads, separations, and normal / lognormal /
  exponential families. Expected outputs: a valid decision label, p-values in
  (0, 1], `p_greater + p_less > 1` (both tails count ties), the resample count
  inside its configured bounds, and a directional decision agreeing with the
  sign of the statistic. Purpose: structural invariants on shapes no fixed
  fixture covers.
- `test_equal_samples_give_a_zero_statistic_and_no_direction` — Inputs: the same
  8 values as both samples. Expected outputs: statistic 0 and a non-directional
  decision (`equivalent` or `inconclusive`). Purpose: identical inputs cannot
  produce a direction.
- `test_equivalence_needs_enough_data_to_conclude` — Inputs (parametrized, 2
  rows): two samples drawn from one normal distribution at n=8 and n=400 per
  group. Expected outputs: `inconclusive` at n=8 and `equivalent` at n=400, with
  `p_tost <= 0.05` exactly on the `equivalent` row. Purpose: the equivalence
  margin is an effect size, so it shrinks relative to the standard error as data
  accrues and TOST gains power; a margin in standard-error units would report
  `inconclusive` at every n.
- `test_type_one_error_rate_tracks_alpha_under_unequal_variance` — Inputs: 300
  null replicates, n=60 at sd 1.0 against n=180 at sd 0.3, alpha 0.05. Expected
  outputs: between 3 and 30 *directional* decisions (1%-10%; nominal is 15);
  non-directional outcomes are not rejections of the directional null. Purpose: the
  studentization holds the nominal level under unequal sizes and variances,
  which is the regime where an unstudentized permutation test fails.
- `test_seeded_runs_are_reproducible` — Inputs: two identical calls with the
  same seed. Expected outputs: identical results. Purpose: reproducibility.
- `test_guards_short_circuit_without_resampling` — Inputs (parametrized, 4
  rows): one input per guard condition. Expected outputs: the matching
  `note`, an `inconclusive` decision, zero resamples, no statistic, and
  `converged is False`. Purpose: each guard is reachable and inert.
- `test_null_skewness_is_measured_and_matches_scipy` — Inputs: 700 concordant
  heights against 19 discordant1 of which 4 are extreme, 2500 resamples.
  Expected outputs: `null_skew` equal to `scipy.stats.skew` over the same draws,
  `n_resamples_skew == 2500`, `|null_skew| > 1`, a negative statistic, decision
  `less`, and `note is None`. Purpose: the running power-sum accumulation
  reproduces the reference skewness, and an outlier-driven null is measurably
  asymmetric without that disturbing the directional call.
- `test_null_skewness_is_near_zero_for_a_symmetric_null` — Inputs: two balanced
  150-observation samples from one normal family. Expected outputs:
  `|null_skew| < 0.15`. Purpose: the measure reads near zero when the null is
  symmetric, so a large value means something.
- `test_adaptive_run_grows_batches_until_it_converges` — Inputs: a clearly
  separated pair with `min_resamples=1000`. Expected outputs: converged after
  exactly one batch of 1000 with decision `greater`. Purpose: an easy case stops
  at the minimum budget instead of spending the ceiling.
- `test_undecided_runs_grow_their_batches_until_the_budget_is_reached` —
  Inputs: a marginal 0.4 mean shift between two 30-element samples,
  `min_resamples=100`, `max_resamples=1000`. Expected outputs: `batches > 1`, a
  total exceeding `batches x 100`, `1000 <= n_resamples < 2000`,
  `converged is False`, and `note == "max_resamples_reached"`. Purpose: an
  undecided run escalates its effort and stops once the budget is met, drawing
  the final batch whole rather than trimming it — so the total meets or slightly
  overshoots `max_resamples` — without pinning the growth factor (a performance
  knob).
- `test_bootstrap_resample_budget_is_reduced_but_always_usable` — Inputs
  (parametrized over six configured `(min, max)` pairs from `(1, 1)` to
  `(2500, 25000)`). Expected outputs: the scaled minimum is at least 1 and,
  wherever the configured minimum is at least 2, strictly below it; the scaled
  maximum stays at or above the scaled minimum and no higher than the configured
  one. Purpose: a bootstrap iteration must cost less than the point estimate
  while still leaving an adaptive run a valid range to grow through — asserted
  as properties, so the divisor stays a tunable rather than a pinned constant.
- `test_bootstrap_resample_budget_never_shrinks_as_the_budget_grows` — Inputs:
  configured budgets from 1 to 25000. Expected outputs: neither end of the
  scaled budget ever decreases. Purpose: monotonicity, so raising the configured
  budget cannot lower the bootstrap's.

- `test_shifted_statistics_match_an_explicit_shift` — Inputs (parametrized over
  five group-size splits including both orderings): a gamma-drawn pooled sample,
  drawn once with `shifts=(0.75, -0.75)` and once per shift on explicitly shifted
  data under the same seed. Expected outputs: the fused rows equal the explicit
  re-draws to `rel=1e-9`. Purpose: the equivalence test folds its shift into the
  power sums the directional draw already computed; same seed means same
  permutations, so the two must agree to floating point, not merely in
  distribution.
- `test_equivalence_reuses_the_directional_resamples` — Inputs: two same-mean
  normal samples of 40 at a fixed 2000 resamples. Expected outputs: a `p_tost` at
  or above `1 / (n_resamples + 1)` and an exact multiple of it. Purpose: TOST is
  an add-one estimator over the directional test's own draws, so it answers at
  that resolution and draws nothing extra.
- `test_equivalence_test_disabled_reports_no_tost` — Inputs: the same samples with
  `equivalence_test=False`. Expected outputs: `p_tost is None` and the decision
  falls through to `inconclusive`. Purpose: the bootstrap path, which disables the
  equivalence step, pays nothing for it.

### tests/orchestrator/test_orchestrator_shape.py

The optional distribution-shape diagnostics: moment parity against SciPy, the
modality test's calibration on samples of known modality, and the results-TSV
column contract.

- `test_shape_moments_match_scipy` — Inputs: 500 lognormal draws (seed 4).
  Expected outputs: `skew` and `excess_kurtosis` equal `scipy.stats.skew` and
  `scipy.stats.kurtosis(fisher=True)` called on the same sample, and both are
  positive. Purpose: the moment columns are the SciPy estimators, not a
  re-derivation.
- `test_modality_test_rejects_only_a_well_separated_mixture` — Inputs
  (parametrized over 5 shapes): 400 draws each from a normal, lognormal,
  exponential and gamma(2), plus a 200+200 mixture of normals 4 SD apart.
  Expected outputs: `modes_p > 0.05` for the four unimodal families and
  `<= 0.05` for the mixture. Purpose: the test holds its level on skewed
  unimodal shapes — where a raw mode count reports spurious peaks — while still
  detecting a real mixture.
- `test_tail_index_recovers_known_tail_shapes` — Inputs (parametrized over 3
  shapes): 4000 draws from a Pareto(3), an exponential, and a uniform. Expected
  outputs: `tail_xi` within 0.25 of `1/3`, `0`, and `-1` respectively. Purpose:
  the peaks-over-threshold fit recovers positive, zero, and bounded tails, the
  last of which a Hill estimator cannot represent at all.
- `test_shape_is_not_described_for_small_or_flat_groups` — Inputs: five
  identical values, fifty identical values, and 25 normal draws. Expected
  outputs: every field `None` for the first two; for the third the moments are
  populated but `tail_xi` is `None`. Purpose: the two guards
  (`SHAPE_MIN_OBSERVATIONS`, `SHAPE_MIN_TAIL_EXCEEDANCES`) are independent.
- `test_shape_is_measured_once_and_not_per_bootstrap_iteration` — Inputs: the
  same triplet analyzed with 5 and with 60 bootstrap iterations. Expected
  outputs: identical `con_skew` and `con_modes_p`. Purpose: the diagnostics come
  from the point estimate only, so iteration count cannot move them and the
  modality bootstrap is not paid per iteration.
- `test_summary_statistics_tsv_never_carries_shape_columns` — Inputs
  (parametrized over `shape_diagnostics` on/off): the same triplet written
  through `write_summary_statistics_tsv`. Expected outputs: no column ending in
  any `SHAPE_FIELD_NAMES` suffix appears in the header under either setting, the
  row stays aligned with it, and the descriptive per-topology columns and
  `classification` are still present. Purpose: that file is consumed as a
  feature matrix and the diagnostics are undefined below their observation
  floors, so they must not reach it even when measured.
- `test_results_tsv_carries_shape_columns_only_when_enabled` — Inputs
  (parametrized over `shape_diagnostics` on/off): a 60/40/30 lognormal triplet
  written through `write_pipeline_results`. Expected outputs: the fifteen
  `<group>_<field>` columns present exactly when enabled, the row aligned with
  the header either way, and with the diagnostics on a positive `con_skew`, a
  `con_modes_p` in `(0, 1]`, and an empty `dis2_tail_xi` (30 observations leave
  only 3 in the upper decile). Purpose: the column set is read off the results
  rather than plumbed separately, so it cannot claim a measurement that did not
  happen.

### tests/orchestrator/test_orchestrator_trees.py

Tree preprocessing, asserted against explicit Newick literals.

- `test_clean_and_save_trees_preserves_a_well_supported_tree` — Inputs:
  `orchestrator_species_tree` with `min_avg_support=0.5`. Expected outputs: the
  cleaned file round-trips the input Newick verbatim and reads back as one tree.
  Purpose: support-free trees pass the filter unchanged.
- `test_clean_and_save_trees_drops_low_average_support` — Inputs:
  `low_support_tree_file` with `min_avg_support=0.5`. Expected outputs: exactly
  one tree survives (mean 0.973 kept, mean 0.300 dropped) carrying the expected
  taxa, and support values are stripped from the output. Purpose: the mean
  support filter.
- `test_root_tree_on_outgroup_prunes_and_reports_ingroup` — Inputs: the cleaned
  species tree and outgroup `OUT`. Expected outputs: `excluded == {"OUT"}`, no
  missing taxa, ingroup `[A, B, C, D]`, and the pruned Newick
  `(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;`. Purpose: outgroup rooting folds
  the removed node's edge into its sibling.
- `test_generate_triplets_and_species_subtrees` — Inputs: the pruned species
  tree. Expected outputs: the 4 sorted triplets from 4 ingroup taxa, no skipped
  triplets, identity ABC normalization, and the exact per-triplet species
  subtree Newick strings. Purpose: triplet enumeration and species-subtree
  construction.
- `test_clean_and_save_gene_trees_roots_every_tree_on_the_outgroup` — Inputs:
  `orchestrator_gene_trees` with outgroup `OUT`. Expected outputs: all 12 trees
  survive, each ends in `OUT:0);`, and trees 0 and 3 match their expected
  rerooted Newick (the ingroup edge absorbs OUT's original edge length).
  Purpose: gene-tree rooting semantics.
- `test_extract_triplet_subtree_selects_the_triplet_taxa` — Inputs: the first
  cleaned gene tree and triplet `(A, B, C)`. Expected outputs: a subtree whose
  leaf set is exactly `{A, B, C}`. Purpose: the DendroPy reference extraction
  the geometry parity tests compare against — no run calls it, so its own
  correctness has to be pinned here.

### tests/orchestrator/test_orchestrator_tree_parity.py

The suite's only parity tests, both marked `@pytest.mark.parity`.

- `test_triplet_branch_lengths_match` — Inputs (parametrized over the pairs
  `(A,B)`, `(A,C)`, `(B,C)`): each `triplet_comparison_cases` tree, extracted
  with DendroPy then re-read through BioPython. Expected outputs: the pairwise
  patristic distance is the same in both libraries. Purpose: the Newick round
  trip preserves branch lengths across libraries.
- `test_triplet_collapse_consistency_dendropy_vs_biopython` — Inputs: the same
  cases collapsed by `extract_triplet_subtree` and by a BioPython pruning
  reference. Expected outputs: all three pairwise distances agree to `1e-12`.
  Purpose: DendroPy's triplet collapsing matches standard BioPython pruning.

### tests/orchestrator/test_orchestrator_triplet_geometry.py

Covers `orchestrator/triplet_geometry.py`, which reads a triplet's geometry out
of a cached tree instead of extracting its subtree — the path every run takes,
for gene trees and the species tree alike. Most tests compare against
`extract_triplet_subtree` + `observation_from_subtree`, the DendroPy reference
behaviour.

- `test_geometry_matches_subtree_extraction` - Inputs (parametrized over the six
  tree-height strategies and eleven `(newick, triplet)` cases covering nested
  pairs, pairs spanning the root, pruned extra taxa, a ladder, absent and mixed
  edge lengths, a zero-length internal branch, and near-degenerate lengths): the
  same tree and triplet through `extract_triplet_subtree` +
  `observation_from_subtree` and through `build_triplet_geometry` +
  `geometry_observation`. Expected outputs: identical topology, and tree height
  and all three summary metrics equal within `rel=1e-12`. Purpose: the cached
  path is a drop-in for extraction across every strategy and tree shape.
- `test_geometry_omits_summary_metrics_when_not_collecting` - Inputs: one case
  with `collect_summary_statistics` false. Expected outputs: both paths leave
  the metrics slot `None`. Purpose: the observation contract is unchanged.
- `test_geometry_on_a_hand_derived_tree` - Inputs: a five-taxon tree whose
  distances are worked out by hand. Expected outputs: the stated topology,
  height, internal branch, and sister distance. Purpose: pins the arithmetic to
  a derivation rather than only to the other implementation.
- `test_geometry_skips_exactly_what_extraction_skips` - Inputs (parametrized): a
  root polytomy and a triplet with an absent taxon. Expected outputs: both paths
  return `None`. Purpose: the skip decisions match, so observation counts do.
- `test_zero_length_internal_branch_still_resolves` - Inputs: `((A:1.0,B:1.0):0.0,C:1.0);`.
  Expected outputs: both paths resolve `((A,B),C)` with a zero internal branch.
  Purpose: the sister pair is chosen by LCA node identity, so a zero-length
  internal branch is not mistaken for a polytomy.
- `test_missing_edge_lengths_count_as_zero` - Inputs: a Newick with no lengths.
  Expected outputs: every derived value is zero. Purpose: matches
  `_distance_to_root` treating a missing length as zero.
- `test_one_cache_serves_every_triplet_in_the_tree` - Inputs: one cache queried
  for all ten triplets of a five-taxon tree. Expected outputs: each matches the
  extraction path. Purpose: one build answers every triplet, which is the reuse
  the cache exists for.
- `test_triplet_resolution_agrees_with_the_observation_guards` - Inputs: a
  resolved triplet, one across a zero-length internal branch, a root polytomy,
  an absent taxon, and then all 84 triplets of the nine-taxon tree. Expected
  outputs: `triplet_resolution` returns the named status, and returns
  `resolved` exactly when `geometry_observation` returns an observation.
  Purpose: the diagnostic preflight uses restates the hot path's guards without
  sharing code, so this is what stops the two drifting apart.
- `test_geometry_matches_dendropy_across_a_nine_taxon_tree` - Inputs
  (parametrized over the six strategies): all 84 triplets of a nine-taxon tree
  carrying an outgroup, nested clades, a ladder, uneven branch lengths and one
  zero-length internal branch. Expected outputs: identical topology and skip
  decision on every triplet, with heights and metrics equal within `rel=1e-12`.
  Purpose: a whole-tree sweep rather than hand-picked shapes, so sister pairs on
  either side of the root and across the zero-length branch are all covered.

### tests/orchestrator/test_orchestrator_rename_map.py

Covers loading and validating a species rename map. The end-to-end effect on
the outputs is covered by `test_species_rename_map_reaches_every_output` in
`test_orchestrator.py`.

- `test_rename_map_reads_a_two_column_tsv` — Inputs: a TSV with a comment line,
  a blank line, and two entries. Expected outputs: the two-entry mapping.
  Purpose: the TSV form, and that blanks and comments are ignored.
- `test_rename_map_reads_a_yaml_mapping` — Inputs: the same pairs as YAML.
  Expected outputs: the same mapping. Purpose: format chosen by extension.
- `test_rename_map_rejects_malformed_files` — Inputs (parametrized): a TSV row
  with three columns, one with a single column, a repeated label, two labels
  sharing a display name in each of the TSV and YAML forms, and a YAML list
  rather than a mapping. Expected outputs: `ValueError` naming the problem.
  Purpose: malformed maps fail at load rather than silently renaming nothing.
- `test_rename_map_rejects_a_missing_file` — Inputs: a path that does not exist.
  Expected outputs: `FileNotFoundError` naming the path. Purpose: a mistyped
  path is reported as such.
- `test_renaming_a_tree_touches_only_mapped_terminals` — Inputs: a three-taxon
  tree and a map covering two of them. Expected outputs: two renames, the third
  label unchanged. Purpose: partial maps are supported.
- `test_renaming_labels_leaves_unmapped_names_alone` — Inputs: label tuples with
  a partial map and an empty map. Expected outputs: unmapped labels pass
  through. Purpose: the helper used for the outgroup and triplet-filter entries.

### tests/orchestrator/test_orchestrator_preflight.py

The structural preflight data check and the runner short-circuit that reaches it.

- `test_clean_inputs_pass_with_no_issues` — Inputs: a 5-taxon species tree
  `(((A,B),C),(D,OUT))` and two well-formed gene trees. Expected outputs:
  `passed is True`, an empty `issues` list, `triplets_checked == 4`, and both
  gene trees rooted. Purpose: a clean dataset produces no false positives.
- `test_report_is_written_only_when_an_output_dir_is_given` — Inputs: the clean
  dataset checked twice, once with an explicit output directory and once with
  `output_dir=None`. Expected outputs: with a directory, `report_path` points at
  `preflight_data_check.txt` inside it and the file content equals
  `report_text`; without one, `report_path is None` and the same report text
  comes back. Purpose: the report is persisted where documented, and the check
  is usable without touching disk.
- `test_detects_polytomy_and_missing_outgroup` — Inputs: gene tree 1 well
  formed, gene tree 2 a polytomy over A/B/C, gene tree 3 with no outgroup
  label. Expected outputs: exactly one `gene_tree.rooting_failed` and one
  `triplet.unresolved_rooted_sister_pair`, two trees rooted out of three
  checked, the polytomy message naming `Gene tree #2` and `A,B,C`, and the
  pair counters accounting for all 8 triplet/gene-tree pairs as 7 usable, 1
  unresolved, 0 with an absent taxon. Purpose: each defect class is detected
  once and located precisely, and every pair the check looked at is accounted
  for.
- `test_triplet_filter_entries_are_validated` — Inputs: a filter file with one
  valid line, one naming an unknown taxon, one naming the outgroup. Expected
  outputs: one `triplet_filter.taxa_missing_in_species_tree`, one
  `triplet_filter.includes_outgroup`, and `triplets_checked == 1`. Purpose:
  filter entries are validated rather than silently dropped.
- `test_impossible_checks_raise` — Inputs (parametrized): an empty outgroup
  list, and an outgroup absent from the species tree. Expected outputs:
  `ValueError` matching "No outgroup taxa were provided" and "Could not root
  species tree". Purpose: conditions that make the check impossible fail loudly.
- `test_multi_tree_species_file_raises` — Inputs: a species-tree file holding
  two trees. Expected outputs: `ValueError` matching "exactly one tree".
  Purpose: the single-tree precondition is enforced.
- `test_runner_preflight_mode_skips_analysis` — Inputs: a config dict with
  `preflight_data_check: True` against the defective dataset. Expected outputs:
  the returned result has `passed is False` and the output directory contains
  only `preflight_data_check.txt`. Purpose: the flag runs the check and nothing
  else.
- `test_runner_returns_none_when_preflight_cannot_run` — Inputs: the same
  config with an outgroup absent from the species tree. Expected outputs:
  `run_orchestrator` returns `None`. Purpose: an impossible check is reported,
  not raised out of the runner.

### tests/orchestrator/test_orchestrator_config.py

Orchestrator config resolution and config-file precedence.

- `test_cli_defaults_resolve` — Inputs: a CLI namespace with every optional arg
  `None`. Expected outputs: `alpha_dct`/`alpha_ks` 0.05, the orchestrator-specific
  `p_value_correction == "bfn"`, `alpha_perm == 0.05`, the permutation
  resample/CI defaults (2500, 25000, `wilson`), `overwrite is True`, the config-file-only keys at their defaults,
  `preflight_data_check is False`, `pipeline_mode == "efficient"`, and
  and `pipeline_mode == "efficient"`. Purpose: default resolution in CLI mode.
- `test_pipeline_mode_rejects_an_unknown_value` — Inputs: `pipeline_mode: fast`.
  Expected outputs: `ConfigError` naming the field. Purpose: the choice list is
  enforced.
- `test_seed_rejects_a_non_integer` — Inputs: `seed: "abc"`. Expected outputs:
  `ConfigError` naming the field. Purpose: the run-wide seed must be an integer
  when provided.
- `test_parser_flags_resolve_into_their_config_values` — Inputs: the CLI flag
  strings parsed by `build_argument_parser`, then resolved. Expected outputs:
  each flag's value reaches its config key, `--no-overwrite` gives
  `overwrite is False`, and `--preflight-data-check` gives `True`. Purpose: every
  other config test builds a namespace directly, so this is the only place the
  flag names are pinned; resolving covers the override path in the same pass.
- `test_config_only_keys_read_from_config_file` — Inputs: a JSON config setting
  `discordant_test`, `tree_height_calculation_strategy`, `min_support_value`,
  `generate_summary_stats`, `alpha_dct`, and a nested `bootstrap_options` block.
  Expected outputs: every key, including the flattened bootstrap options, is
  read. Purpose: config-file-only keys and nested bootstrap parsing.
- `test_p_value_correction_accepts_yaml_bare_word_no` — Inputs (parametrized):
  a YAML config writing `p_value_correction` as the bare word `no`, as quoted
  `"no"`, as `No`/`NO`, and as the unaffected `bfn`/`fdr_bh`. Expected outputs:
  each resolves to its own choice. Purpose: YAML 1.1 resolves bare `no` to
  boolean `False`, so
  the value never reaches validation as a string; the mapping back to the
  written choice is what makes the unquoted spelling work.
- `test_p_value_correction_rejects_a_value_with_no_matching_choice` — Inputs: a
  YAML config writing `p_value_correction: yes`. Expected outputs: `ConfigError`
  matching `must be one of`. Purpose: the boolean mapping must not turn an
  invalid value into a valid one — `yes` has no corresponding choice and still
  fails, with the received value named.
- `test_config_file_wins_over_cli` — Inputs: a config file plus conflicting CLI
  flags. Expected outputs: the file's `alpha_dct` wins, `alpha_perm`
  falls back to the orchestrator default (proving the CLI value was ignored),
  and the file's paths are used. Purpose: config-file precedence.
- `test_outgroup_accepts_single_comma_separated_and_list_forms` — Inputs
  (parametrized, 6 rows): `outgroup` given as a single label, a comma-separated
  string, a padded string with a trailing comma, a list, a tuple, and a list
  whose entries are themselves comma-separated. Expected outputs: each resolves
  to the same flat label list. Purpose: one key covers the single- and
  multiple-outgroup cases in every accepted shape.
- `test_outgroup_rejects_empty_and_non_label_values` — Inputs (parametrized, 7
  rows): `None`, empty and whitespace strings, a lone comma, empty and
  blank-only lists, and an integer. Expected outputs: `ConfigError` naming
  `outgroup`. Purpose: a value that yields no labels is an error rather than an
  empty outgroup list.
- `test_shipped_sample_configs_resolve` — Inputs (parametrized):
  `sample_configs/orchestrator_minimal.yaml` and `orchestrator_full.yaml`.
  Expected outputs: each loads without error, yields a non-empty list of
  string outgroup labels, and carries the current `alpha_perm` and
  `permutation_ci_method` defaults. Purpose: the shipped samples cannot drift
  out of step with the validator and leave users copying a rejected config.
- `test_full_sample_config_covers_every_runtime_key` — Inputs:
  `orchestrator_full.yaml` read both as raw YAML and through the loader.
  Expected outputs: every key the normalizer produces is documented in the
  sample, after allowing for the three renamed path keys and the two
  prefix-flattened nested blocks. Purpose: a new config key cannot be added
  without the sample gaining it too.
- `test_missing_required_field_raises` — Inputs: a namespace missing the species
  tree. Expected outputs: `ConfigError`. Purpose: required-field validation.
- `test_preflight_data_check_resolves_from_config_file` — Inputs
  (parametrized): a JSON config setting `preflight_data_check` to `true`,
  `false`, or omitting it. Expected outputs: `True`, `False`, and `False`
  respectively. Purpose: the flag is settable from a config file and defaults
  to off.
### tests/test_config_trunk.py

The shared configuration trunk in `ghostparser.config`.

- `test_resolve_path_handles_absolute_relative_and_home` — Inputs: an absolute
  path, a relative path resolved from a chdir'd cwd, and a `~/` path. Expected
  outputs: each resolves to the correct absolute path with no `~` remaining.
  Purpose: path-resolution rules shared by all modules.
- `test_load_raw_config_reads_json_and_yaml` — Inputs (parametrized over
  `.json`, `.yaml`, `.yml`): equivalent payloads. Expected outputs: the same
  mapping from each format. Purpose: config-file loading.
- `test_load_raw_config_rejects_missing_unsupported_and_non_mapping` — Inputs: a
  missing path, a `.txt` file, and a JSON list. Expected outputs:
  `FileNotFoundError`, then `ConfigError` twice with the documented messages.
  Purpose: loader validation.
- `test_validate_required_path_resolves_or_raises` — Inputs: a present path, then
  absent/empty/whitespace values. Expected outputs: resolution, then
  `ConfigError` for each invalid case. Purpose: required-path validation.
- `test_validate_overwrite_flag_precedence_and_validation` — Inputs: every
  combination of `overwrite`/`no_overwrite`, an explicit default, and
  non-boolean values. Expected outputs: `overwrite` wins over `no_overwrite`,
  negation is applied correctly, and non-booleans raise `ConfigError`. Purpose:
  the shared overwrite semantics.
- `test_prepare_output_directory_overwrites_or_suffixes` — Inputs: an existing
  directory with a stale file, then the same directory with `results_1` and
  `results_3` already taken and `overwrite=False`. Expected outputs: the
  directory is reset in the first case; the second returns `results_2` (the
  smallest missing suffix) and leaves the original intact; a nested path with
  missing parents is created on demand. Purpose: output directory preparation,
  suffix allocation, and parent creation.
### tests/orchestrator/test_orchestrator_consolidation.py

Consolidation outputs, count aggregation, and plot rendering.

- `test_generate_introgression_maps_creates_expected_outputs` — Inputs: synthetic
  triplet results and a species tree. Expected outputs: the combined PNG and all
  three TSVs are written. Purpose: artifact generation.
- `test_collect_counts_correct_avg_in_generate_introgression_maps` — Inputs:
  classified triplet results. Expected outputs: averages match the documented
  co-occurrence denominators. Purpose: bootstrap averaging.
- `test_collect_non_sister_counts_counts_non_sister_pairs` — Inputs: results with
  non-sister pairs. Expected outputs: only non-sister pairs are counted.
  Purpose: pair selection.
- `test_generate_introgression_maps_selects_the_requested_taxa` — Inputs
  (parametrized, 3 rows): a balanced four-taxon tree with no taxon filter, the
  same tree with `plot_taxa=["A", "B", "C"]`, and an outgroup tree with
  `outgroups=["OG"]`. Expected outputs: `taxa_count` of 4, 3, and 3; the matrix
  header and the ghost rows listing the same taxa; and the excluded taxon absent
  from both. Purpose: `plot_taxa` and `outgroups` decide which taxa reach every
  output, and no output disagrees with another about the set.
- `test_collect_counts_counts_only_the_rows_that_produced_an_edge` — Inputs: five
  result rows covering a classified edge, a `no_introgression` row over the same
  taxa, an unrelated triplet, and two ghost rows on either discordant topology.
  Expected outputs: the supporting count for `(C, B)` is 1 and for `(B, D)` is 0;
  ghost counts are 1 for `A` and `C`, 0 for `B` and `D`. Purpose: supporting
  counts follow the classification, not mere co-occurrence.
- `test_generate_introgression_maps_uses_raw_values_with_separate_scales` —
  Inputs: results spanning a value range. Expected outputs: raw values with
  per-plot scales. Purpose: colour scaling.
- `test_generate_introgression_maps_appends_suffix_when_overwrite_disabled` —
  Inputs: an existing output directory with `overwrite=False`. Expected outputs:
  a suffixed sibling directory. Purpose: overwrite behavior.
- `test_sampled_introgression_presence_flags_targets_with_sampled_edges` —
  Inputs: a taxa order and a `(source, target)` weight map with one zero-weight
  edge. Expected outputs: `{"A": 1, "B": 0, "C": 1, "D": 0}` — only taxa that
  are the target of a non-zero sampled edge are flagged. Purpose: the flag that
  drives ghost bar colour.
- `test_ghost_strength_tsv_records_sampled_introgression_flag` — Inputs: results
  where taxon A has both ghost and sampled introgression and taxon D has ghost
  only. Expected outputs: the ghost TSV header is
  `target_taxon / raw_strength / has_sampled_introgression`, with `A → 1` and
  `D → 0`, and the strengths are unchanged by the flag. Purpose: the new column
  and its cross-referencing against the sampled sheet.

### tests/test_ml_labels_and_metrics.py

The ML label contract, evaluation metrics, distributions, and CV-fold policy.

- `test_bit_labels_stay_in_step_with_the_bit_count` — Inputs: the module
  constants. Expected outputs: `len(BIT_LABELS) == BIT_COUNT`, all distinct.
  Purpose: a label can be read back positionally only while the two constants
  agree; the bitstring width itself is pinned behaviorally by
  `test_is_valid_bitstring`.
- `test_bit_label_titles_keep_the_taxon_letters_upper_case` — Inputs
  (parametrized over all six bit labels): each `BIT_LABELS` entry. Expected
  outputs: underscores become spaces and only the first character is
  upper-cased, so `ghost_into_A` renders `Ghost into A`. Purpose: the plot-title
  form, and specifically that `str.capitalize` is not used — it would lower-case
  the taxon letters and rename the taxon.
- `test_every_bit_label_has_a_title` — Inputs: the whole `BIT_LABELS` tuple.
  Expected outputs: six distinct titles, none retaining an underscore and each
  starting upper-case. Purpose: no label falls through the formatter and reaches
  a figure as a raw slug.
- `test_64_class_matrix_orders_classes_by_set_bits` — Inputs: three rows with
  true classes `110000`, `000001`, `111111` and predicted classes `000011`,
  `000001`, `111111`. Expected outputs: 64 distinct `class_labels` whose set-bit
  count never decreases along the list and which are lexically sorted within
  each count, `000000` first and `111111` last; the matrix sums to 3 with a 1 at
  each (true, predicted) position looked up by label. Purpose: both axes of the
  64-class matrix run from no bits set to all six, in the same order, and the
  counts are permuted along with the labels — a label-only sort would put every
  off-diagonal count under the wrong pair of names.
- `test_is_valid_bitstring` — Inputs (parametrized, 8 cases): valid and invalid
  strings. Expected outputs: only six-character 0/1 strings validate. Purpose:
  label validation.
- `test_parse_classes_round_trips_bitstrings` — Inputs: four labels including a
  whitespace-padded one. Expected outputs: a 4x6 binary matrix whose rows
  re-join to the trimmed labels, with the expected total set-bit count. Purpose:
  the bitstring/matrix round trip.
- `test_parse_classes_rejects_malformed_labels` — Inputs (parametrized): wrong
  length or non-binary labels. Expected outputs: `ValueError`. Purpose: label
  validation.
- `test_select_feature_names_excludes_the_target_column` — Inputs: a header list.
  Expected outputs: header order preserved, target column dropped. Purpose:
  feature selection.
- `test_evaluate_predictions_matches_hand_computed_metrics` — Inputs: a 2x6
  true/predicted pair differing in one bit. Expected outputs: exact-match 0.5,
  bitwise 11/12, hamming 1/12, and the expected per-bit recall/support. Purpose:
  metric definitions.
- `test_build_prediction_rows_reports_matched_label_count` — Inputs: the same
  pair. Expected outputs: per-row matched counts 6 and 5, correct exact-match
  flags, label strings, and per-bit columns. Purpose: the predictions.tsv
  contract.
- `test_summarize_distribution_counts_and_fractions` — Inputs: four labels with
  one repeat. Expected outputs: sorted labels with correct counts and fractions.
  Purpose: class distribution.
- `test_bit_distribution_counts_positives_per_bit` — Inputs: a 2x6 target matrix.
  Expected outputs: per-bit positive counts and fractions. Purpose: bit
  distribution.
- `test_build_feature_importance_rows_sorts_descending` — Inputs: three named
  features with scores. Expected outputs: rows ranked most-important first.
  Purpose: importance reporting.
- `TestAutoCvFolds` (6 tests) — Inputs: label arrays whose smallest class varies,
  under each `rare_class_policy`. Expected outputs: folds kept, reduced to the
  smallest class count, skipped, or raising for a singleton class; empty labels
  return no folds. Purpose: full CV-fold policy coverage.
- `test_read_tsv_rows_rejects_a_header_only_file` /
  `test_read_tsv_rows_reads_records` — Inputs: a header-only TSV and a
  two-row TSV. Expected outputs: `ValueError`, then one dict per data row.
  Purpose: input reading.

### tests/test_ml_config.py

- `test_ml_config_fills_defaults_for_omitted_keys` — Inputs (parametrized over
  4 keys): a config carrying only `input_path` and `output_dir`. Expected
  outputs: `target_column` `class`, `overwrite` `True`, `min_samples_split` `2`,
  `min_samples_leaf` `1`. Purpose: the documented defaults are what an otherwise
  empty config resolves to.
- `test_ml_config_explicit_values_win_over_defaults` — Inputs (parametrized over
  4 payloads): `overwrite: false`, an explicit `target_column`, and two
  hyperparameters given under the nested `model` block. Expected outputs: each
  value replaces its default. Purpose: precedence, plus the `model` block
  flattening — a key given there must surface at the top level of the resolved
  config, including one that also has a default.
- `test_ml_config_normalizes_estimator_passthrough_values` — Inputs
  (parametrized over 12 cases): `model.max_features` as `null`, the strings
  `None`/`none`/`null`, `log2`, `3` and `0.5`; `model.class_weight` as `null`,
  `none`, `balanced`, a mapping, and a list of mappings. Expected outputs: every
  null spelling resolves to `None` and every other form passes through
  unchanged. Purpose: YAML reads the bare words `None`/`none` as strings, so
  without this they reach scikit-learn as literal text and are rejected there
  instead of in config.
- `test_ml_config_rejects_invalid_estimator_values` — Inputs (parametrized over
  7 cases): `max_features` as `auto`, `sqrt2`, `0`, `1.5` and `true`;
  `class_weight` as `nope` and `5`. Expected outputs: `ConfigError` whose
  message names the dotted config location, the received value, and the
  accepted forms — `auto` gets its own message saying scikit-learn removed it.
  Purpose: the failure has to name the config key the user wrote, not surface
  as a bare estimator error at fit time. `true` also pins that a bool is
  rejected rather than passing the `int` check.

### tests/test_ml_utils.py

- `test_rows_to_matrix_uses_numeric_features_and_excludes_target_column` —
  Inputs: rows with numeric features plus the target column. Expected outputs: a
  numeric matrix excluding the target.
- `test_rows_to_matrix_encodes_multiple_string_columns` — Inputs: rows with
  low-cardinality string columns. Expected outputs: one-hot encoded features.
- `test_rows_to_matrix_rejects_string_features` — Inputs: a string column
  exceeding the cardinality limit. Expected outputs: an error rather than an
  invented ordering.
- `test_row_normalize_confusion_matrix_turns_counts_into_per_class_fractions` —
  Inputs: a 3x3 count matrix `[[3, 1, 0], [0, 0, 0], [1, 1, 2]]` whose middle
  true class has no samples. Expected outputs: cells become
  `[[0.75, 0.25, 0], [0, 0, 0], [0.25, 0.25, 0.5]]`, every value lies in
  `[0, 1]`, and the row sums are `1, 0, 1`. Purpose: the 64-class
  confusion-matrix plot is drawn on a fixed 0-1 scale, so the normalization has
  to put populated rows on a common scale and leave an empty class at zero — the
  value the plot masks — instead of dividing by zero.

### tests/test_ml_random_forest.py

- `test_train_random_forest_smoke` — Inputs: `summary_statistics_tsv`. Expected
  outputs: training completes; the returned metrics name the objective and carry
  both metric tiers, the dataset summary and the 64-class confusion matrix; the
  written metrics JSON carries the bit-label order and one timing per stage; and
  the predictions TSV, model pickle and both confusion-matrix figures are
  written. Purpose: the one end-to-end smoke test for this entry point, plus its
  output-file and metrics-field contract.

### tests/test_ml_multi_knn.py

- `test_multi_knn_train_smoke` — Inputs: `summary_statistics_tsv`. Expected
  outputs: training completes; the returned metrics carry both metric tiers and
  the cross-validation block, and the predictions TSV, model pickle and both
  confusion-matrix figures are written. Purpose: the one end-to-end smoke test
  for this entry point, plus its output-file contract.
- `test_multi_knn_build_model_caps_neighbors_to_training_size` — Inputs:
  `n_neighbors=20` against training sets of 2 and of 50. Expected outputs: the
  effective count and the estimator's own `n_neighbors` are both `2` in the
  first case, and `20` passes through untouched in the second. Purpose: asking
  KNN for more neighbors than training samples raises at fit time, so the clamp
  has to reach the estimator and must not fire when it isn't needed.

### tests/test_ml_hyper_tune.py

- `test_load_hyper_tune_config_accepts_hyperparameter_tuning_section` — Inputs: a
  tuning config with `use_wandb: false`. Expected outputs: the section loads and
  `use_wandb` resolves to `False`.
- `test_load_hyper_tune_config_fills_model_defaults` — Inputs (parametrized):
  a tuning block omitting model parameters, over `use_wandb` absent/false/true
  and `wandb_detailed_payloads` absent or `true`. Expected outputs: trainer
  defaults are filled in, an absent `use_wandb` resolves to `False`, and both
  flags otherwise follow what was written. Purpose: defaults and explicit values
  both resolve out of the nested `hyperparameter_tuning` block, and W&B logging
  is opt-in rather than reached for by a config that never mentions it.
- `test_load_hyper_tune_config_rejects_invalid_wandb_choice` — Inputs
  (parametrized): a tuning block with `use_wandb` set to the string `"yes"`, and
  one with `wandb_detailed_payloads: true` alongside `use_wandb: false`. Expected
  outputs: `ConfigError` in each case. Purpose: `use_wandb` may be omitted but
  not mistyped, and it gates the detailed-payload flag rather than letting it be
  silently ignored.
- `test_shipped_tuning_sample_configs_resolve` — Inputs (parametrized):
  `sample_configs/hyperparameter_tuning_random_forest.yaml` and
  `sample_configs/hyperparameter_tuning_multi_knn.json`. Expected outputs: both
  load, resolve to their stated model, report `use_wandb` and
  `wandb_detailed_payloads` as `False`, and use only search-space keys their
  model supports. Purpose: the shipped samples have not drifted out of step with
  the validator in either config format.
- `test_load_hyper_tune_config_rejects_evaluation_section` — Inputs: a tuning
  config containing `evaluation`. Expected outputs: `ConfigError`.
- `test_load_hyper_tune_config_requires_hyperparameter_tuning_section` — Inputs:
  a config without the section. Expected outputs: `ConfigError`.
- `test_tune_hyperparameters_defaults_to_wandb_off` — Inputs: a tuning namespace
  with the `use_wandb` attribute deleted. Expected outputs: the search runs,
  `results.use_wandb` is `False`, no `wandb/` directory is written, and
  `hyper_tune_results.json` is. Purpose: the programmatic entry point takes the
  same default as the config loader, so a caller building its own namespace need
  not know the key exists.
- `test_tune_hyperparameters_grid_search_smoke` /
  `test_tune_hyperparameters_random_search_smoke` — Inputs:
  `summary_statistics_tsv_tuning` with each search method and `use_wandb=False`.
  Expected outputs: the search completes, `model_name`/`search_method` echo the
  request, `use_wandb` echoes to the results payload, the candidate count matches
  the search space (2 for the grid, `n_iter=1` for the random search), the
  best-model pickle and results JSON are written, and no `wandb/` directory is
  created.
- `test_tune_hyperparameters_writes_local_navigation_artifacts` — Inputs:
  `summary_statistics_tsv_tuning`, grid search, `use_wandb=False`, `top_k=3`.
  Expected outputs: `hyper_tune_parameter_marginals.tsv` and
  `hyper_tune_search_report.png` exist and the plot path is returned; the ranked
  candidate TSV carries `rank`, `is_best`, `cv_score` and `elapsed_seconds`; the
  marginals TSV header matches the documented column order; the text report
  contains the search-space, top-candidate, per-parameter, guidance and artifact
  sections and records W&B as disabled; and the results JSON carries
  `parameter_marginals`, `artifact_paths`, and one `parameter_influence` entry
  per search-space key. Purpose: with W&B off, the local artifacts are the only
  way to navigate a search, so their file and column contract is pinned.
- `test_tune_hyperparameters_routes_bulk_artifacts_to_wandb` — Inputs:
  `summary_statistics_tsv_tuning`, grid search, `use_wandb=True`, with
  `_import_wandb` patched to a stub module. Expected outputs: the output
  directory holds `hyper_tune_best_model.pkl`, `hyper_tune_results.txt` and
  `hyper_tune_search_report.png` but none of `hyper_tune_results.json`,
  `hyper_tune_results.tsv`, `hyper_tune_parameter_marginals.tsv` or
  `predictions.tsv`; the returned paths for those four are `None`;
  `artifact_paths` names only the three written files; the stub received the
  `tables/ranked_candidates`, `tables/parameter_marginals` and
  `tables/predictions` log payloads; the run summary carries
  `bulk_artifacts_written_locally` false and a `results_json` blob that parses
  back to the full payload; and the text report records W&B as enabled and
  points at the run for the rest. Purpose: turning W&B on has to move the bulk
  outputs rather than duplicate them, and nothing may be silently dropped.
- `test_compute_parameter_marginals_summarizes_each_value` — Inputs
  (parametrized over `max` and `min`): four candidates over
  `n_estimators` in `{5, 10}` with scores `0.7/0.4` and `0.9/0.6`, ranked in the
  objective's direction. Expected outputs: values are ranked against each other
  in that direction, the winner's `best_score` is `0.9` under `max` and `0.4`
  under `min`, each value counts 2 candidates and reaches `best_rank` 1, and the
  winner is flagged `upper_bound` / `lower_bound` respectively while the runner-up
  carries no flag. Purpose: the marginal summary drives the guidance block, so
  the direction handling and the bound flag are checked in both directions.
- `test_parameter_influence_ranks_by_best_score_spread` — Inputs: four candidates
  crossing `n_estimators` and `max_depth` with scores `0.9/0.7/0.6/0.4`. Expected
  outputs: `max_depth` ranks first with a best- and mean-score spread of `0.3`,
  ahead of `n_estimators` at `0.2`. Purpose: the influence ordering tells a user
  which dimension to keep searching.
- `test_search_space_guidance_flags_a_dimension_with_no_effect` — Inputs: two
  candidates differing only in `max_features`, both scoring `0.8`. Expected
  outputs: the guidance reports a score spread of zero and names `max_features`
  as having no effect on the objective. Purpose: an inert search dimension is
  the finding a user most needs called out.
- `test_create_run_logger_routes_on_the_use_wandb_choice` — Inputs
  (parametrized): `use_wandb` false and true, with `_import_wandb` patched to a
  stub module. Expected outputs: with `false` the logger reports
  `enabled is False`, forwards nothing, and creates no `wandb/` directory; with
  `true` it forwards the log payload, writes the run summary, finishes the run,
  and creates the directory. Purpose: one config flag is the only thing that
  decides whether the run touches W&B.
- `test_import_wandb_raises_config_error_when_unavailable` — Inputs: `wandb`
  patched to raise `ModuleNotFoundError` on import. Expected outputs:
  `ConfigError` whose message names `use_wandb: false`. Purpose: a missing
  optional dependency surfaces as an actionable config error rather than an
  import traceback.
- `test_load_hyper_tune_config_normalizes_search_space_values` — Inputs
  (parametrized): `search_space.max_features` as the list
  `["sqrt", "none", null, 4]` and as the lone value `"None"`. Expected outputs:
  `["sqrt", None, None, 4]` and `None` respectively. Purpose: search-space
  candidates reach the estimator one at a time, so they need the same per-value
  rules as the `model` block, and the container shape a user wrote (list or
  scalar) has to survive normalization.
- `test_load_hyper_tune_config_rejects_invalid_search_space_value` — Inputs:
  `search_space.max_features: ["sqrt", "auto"]`. Expected outputs: `ConfigError`
  naming `hyperparameter_tuning.search_space.max_features` and `'auto'`.
  Purpose: one bad candidate among valid ones is caught before the search
  starts, rather than partway through as the grid reaches it.

## Parity Tests (`@pytest.mark.parity`)

Run with `pytest -m parity`. Only the two DendroPy-vs-BioPython tests in
`tests/orchestrator/test_orchestrator_tree_parity.py` carry this marker; every other
test derives its expectations from definitions instead of comparing against a
second implementation.
