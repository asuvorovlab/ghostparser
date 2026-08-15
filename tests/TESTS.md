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
  with `triplet` and includes `classification` and `bootstrap_value`, and its
  row count equals the number of results. Purpose: TSV shape and naming.
- `test_no_bootstrap_omits_the_bootstrap_columns` — Inputs: a run with
  `bootstrap=False`. Expected outputs: `bootstrap_value` and `all_bootstrap` are
  absent from the header while `classification` remains, and all 4 triplets are
  still produced. Purpose: the bootstrap toggle only removes bootstrap output.
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
- `test_parallel_modes_match_serial` — Inputs (parametrized over
  `("taxon", 2)` and `("gene", 2)`): the same fixture run serially and in
  parallel with a fixed seed. Expected outputs: every compared field, bootstrap
  values included, is identical. Purpose: parallelization must not change
  results.

Also in `test_orchestrator.py`:

- `test_results_tsv_carries_no_group_summary_columns` — Inputs: a run with
  `generate_summary_stats=True`. Expected outputs:
  `mean_con`/`mean_dis`/`median_con`/`median_dis` are absent from the results
  TSV header even under that flag, while `perm_p_greater`, `perm_p_less`, and
  `permutation_consistency_flag` are present. Purpose: the direction comes from
  the permutation p-values, so per-group summaries belong to
  `summary_statistics.tsv` alone and never appear in the results TSV.

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

- `test_classify_no_introgression_when_dct_not_significant` — Inputs: an even
  10/10 discordant split. Expected outputs: DCT statistic 0, p-value 1.0, not
  significant, `no_introgression`. Purpose: gate 1 short-circuits.
- `test_classify_inflow_when_tree_height_test_not_significant` — Inputs: a
  significant 30/2 discordant split with identical concordant and discordant
  heights. Expected outputs: DCT significant, KS statistic 0 and not
  significant, `inflow_introgression`. Purpose: gate 2 maps to inflow.
- `test_classify_outflow_when_concordant_heights_exceed_discordant` — Inputs: a
  significant 30/2 split with fully separated, spread-out heights (10 concordant
  in [0.85, 0.94], 30 discordant1 in [0.05, 0.34]). Expected outputs: KS
  statistic 1.0 and significant, `mean_con > mean_dis`, no permutation guard,
  `perm_decision == "greater"`, `outflow_introgression`. Purpose: gate 3,
  con > dis, with the permutation test actually resampling.
- `test_classify_ghost_when_discordant_heights_exceed_concordant` — Inputs: the
  mirror case, discordant above concordant. Expected outputs:
  `mean_con < mean_dis`, `perm_decision == "less"`, `ghost_introgression`.
  Purpose: gate 3, con < dis.
- `test_classify_ambiguous_when_direction_is_undetectable` — Inputs: a
  significant split where concordant and discordant1 share a mean but differ
  sharply in spread (+-0.30 versus +-0.02 around 0.5). Expected outputs: DCT and
  KS both significant, `perm_decision == "ambiguous"`, classification
  `ambiguous`. Purpose: the KS test separates distributions that the direction
  test cannot order, which is the case that has no directional answer.
- `test_permutation_guard_reports_insufficient_support` — Inputs: 4 concordant
  and 2 discordant1 heights, giving C(6, 2) = 15 assignments. Expected outputs:
  `perm_note == "insufficient_permutation_support"`, zero resamples, ambiguous.
  Purpose: too-small samples are refused rather than decided.
- `test_permutation_guard_reports_degenerate_scale` — Inputs: internally
  constant groups with different means (`[0.9] * 10` versus `[0.1] * 30`).
  Expected outputs: `perm_note == "degenerate_observed_scale"`, ambiguous
  classification. Purpose: a zero standard error is caught relative to the
  data's magnitude rather than against exact zero.
- `test_median_fallback_decides_direction_when_permutation_disabled` — Inputs:
  the same degenerate heights with `permutation_test=False`. Expected outputs:
  `perm_decision == "greater"`, `perm_statistic is None`,
  `outflow_introgression`. Purpose: the fallback path still resolves direction.
- `test_summary_statistics_discordant1_follows_frequency_not_topology_name` —
  Inputs: 10 concordant, 3 `BC|A`, and 9 `AC|B` gene subtrees with
  `collect_summary_statistics=True`, so `AC|B` is the more frequent discordant.
  Expected outputs: `dis1_topology == "AC"`, `(n_dis1, n_dis2) == (9, 3)`, and
  `discordant1_avg_tree_height_mean` equal to the AC group's derived mean (with
  `discordant2_*` the BC group's), matching `result.mean_dis`. Purpose: the
  summary columns name the same gene trees as `dis1_topology` and the tests,
  rather than a fixed topology label.
- `test_summary_statistics_discordant_roles_swap_with_the_counts` — Inputs: the
  same fixture with the two discordant groups exchanged. Expected outputs:
  `dis1_topology == "BC"` and the two summary column families swap accordingly.
  Purpose: the role assignment tracks the counts in both directions.
- `test_classify_introgression_truth_table` — Inputs (parametrized, 7 rows):
  every combination of DCT/KS significance and direction, including `None`.
  Expected outputs: the documented classification for each row. Purpose:
  exhaustive coverage of `_classify_introgression`.
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
- `test_decision_rules_agree_on_random_inputs` — Inputs (parametrized over 8
  seeds): random sample pairs with Bonferroni correction. Expected outputs:
  `consistent is True` and the two decision rules produce the same label.
  Purpose: the one-tailed rule and the two-tailed-gate-then-sign rule coincide.
- `test_equal_samples_give_a_zero_statistic_and_no_direction` — Inputs: the same
  8 values as both samples. Expected outputs: statistic 0, ambiguous. Purpose:
  identical inputs cannot produce a direction.
- `test_type_one_error_rate_tracks_alpha_under_unequal_variance` — Inputs: 300
  null replicates, n=60 at sd 1.0 against n=180 at sd 0.3, alpha 0.05. Expected
  outputs: between 3 and 30 rejections (1%-10%; nominal is 15). Purpose: the
  studentization holds the nominal level under unequal sizes and variances,
  which is the regime where an unstudentized permutation test fails.
- `test_seeded_runs_are_reproducible` — Inputs: two identical calls with the
  same seed. Expected outputs: identical results. Purpose: reproducibility.
- `test_guards_short_circuit_without_resampling` — Inputs (parametrized, 4
  rows): one input per guard condition. Expected outputs: the matching
  `note`, an ambiguous decision, zero resamples, no statistic, and
  `converged is False`. Purpose: each guard is reachable and inert.
- `test_skewed_null_keeps_the_directional_call_and_flags_it` — Inputs: 700
  concordant heights against 19 discordant1 heights of which 4 are extreme —
  the shape observed on real data. Expected outputs: a negative statistic,
  `decision == "less"`, `p_two_sided > alpha`, `null_skewed is True`,
  `note is None`, and `consistent is True`. Purpose: an asymmetric null is
  recorded without overturning the directional call, and the cross-check's
  doubled-smaller-tail gate stays valid under that asymmetry.
- `test_max_resamples_reached_is_reported` — Inputs: two near-identical samples
  at a 200-resample ceiling. Expected outputs: the budget is respected and, if
  unconverged, `note == "max_resamples_reached"`. Purpose: budget exhaustion is
  surfaced rather than silently treated as convergence.
- `test_adaptive_run_grows_batches_until_it_converges` — Inputs: a clearly
  separated pair with `min_resamples=1000`. Expected outputs: converged after
  exactly one batch of 1000 with decision `greater`. Purpose: an easy case stops
  at the minimum budget instead of spending the ceiling.
- `test_batches_grow_by_one_quarter_until_the_budget_is_spent` — Inputs: two
  30-element samples from the same distribution, `min_resamples=100`,
  `max_resamples=1000`. Expected outputs: the resample total lands on the
  cumulative schedule 100, 225, 381, 576, 819, 1000 produced by `int(previous x
  1.25)` with the last batch clipped, at the index matching `batches`. Purpose:
  pins the 1.25 growth factor and the final-batch clipping.
- `test_median_sign_decision_matches_definition` — Inputs: four median
  comparisons including a tie and an empty sample. Expected outputs:
  `greater` / `less` / `ambiguous` / `ambiguous`. Purpose: the fallback path.
- `test_bootstrap_resample_budget_scales_by_one_fifth` — Inputs: `(2500, 25000)`
  and `(2, 3)`. Expected outputs: `(500, 5000)` and `(1, 1)`. Purpose: the
  bootstrap budget divisor and its floor.

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
  leaf set is exactly `{A, B, C}`. Purpose: subtree extraction.

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

### tests/orchestrator/test_orchestrator_preflight.py

The structural preflight data check and the runner short-circuit that reaches it.

- `test_clean_inputs_pass_with_no_issues` — Inputs: a 5-taxon species tree
  `(((A,B),C),(D,OUT))` and two well-formed gene trees. Expected outputs:
  `passed is True`, an empty `issues` list, `triplets_checked == 4`, both gene
  trees rooted, and the "No blocking data issues detected" line in the report.
  Purpose: a clean dataset produces no false positives.
- `test_report_is_written_to_output_dir` — Inputs: the clean dataset with an
  explicit output directory. Expected outputs: `report_path` points at
  `preflight_data_check.txt` inside it and the file content equals
  `report_text`. Purpose: the report is persisted where documented.
- `test_no_output_dir_skips_writing` — Inputs: the clean dataset with
  `output_dir=None`. Expected outputs: `report_path is None` and non-empty
  `report_text`. Purpose: the check is usable without touching disk.
- `test_detects_polytomy_and_missing_outgroup` — Inputs: gene tree 1 well
  formed, gene tree 2 a polytomy over A/B/C, gene tree 3 with no outgroup
  label. Expected outputs: exactly one `gene_tree.rooting_failed` and one
  `triplet.unresolved_rooted_sister_pair`, two trees rooted out of three
  checked, and the polytomy message naming `Gene tree #2` and `A,B,C`.
  Purpose: each defect class is detected once and located precisely.
- `test_report_attributes_issues_to_gene_trees` — Inputs: the same defective
  dataset. Expected outputs: the report attributes 0 issues to the species tree
  and 2 to the gene trees. Purpose: the attribution summary is correct.
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
- `test_runner_preflight_reports_unrootable_species_tree` — Inputs: the same
  config with an outgroup absent from the species tree. Expected outputs:
  `run_orchestrator` returns `None` and prints "Preflight data check could not
  run". Purpose: an impossible check is reported, not raised out of the runner.

### tests/orchestrator/test_orchestrator_config.py

Orchestrator config resolution and config-file precedence.

- `test_cli_defaults_resolve` — Inputs: a CLI namespace with every optional arg
  `None`. Expected outputs: `alpha_dct`/`alpha_ks` 0.05, the orchestrator-specific
  `p_value_correction == "bfn"`, `alpha_perm == 0.05`,
  `permutation_test is True` with its resample/CI defaults (2500, 25000,
  `wilson`), `overwrite is True`, the config-file-only keys at their defaults,
  `preflight_data_check is False`, and no `stats_backend` key. Purpose: default
  resolution in CLI mode.
- `test_cli_overrides_for_config_plus_cli_options` — Inputs: CLI values for
  alpha-dct/alpha-ks/alpha-perm/p-value-correction/no-overwrite. Expected
  outputs: each override is honored and `overwrite` becomes `False`. Purpose:
  the config+CLI options are wired.
- `test_config_only_keys_read_from_config_file` — Inputs: a JSON config setting
  `discordant_test`, `tree_height_calculation_strategy`, `min_support_value`,
  `generate_summary_stats`, `alpha_dct`, and a nested `bootstrap_options` block.
  Expected outputs: every key, including the flattened bootstrap options, is
  read. Purpose: config-file-only keys and nested bootstrap parsing.
- `test_config_file_wins_over_cli` — Inputs: a config file plus conflicting CLI
  flags. Expected outputs: the file's `alpha_dct` wins, `alpha_perm`
  falls back to the orchestrator default (proving the CLI value was ignored), the
  file's paths are used, and a warning is printed. Purpose: config-file
  precedence.
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
- `test_parser_exposes_config_file_and_new_flags` — Inputs: an argv list using
  the config+CLI flags including `--preflight-data-check`. Expected outputs:
  each parses to its expected value, `preflight_data_check is True`, and
  `config_file` defaults to `None`. Purpose: parser surface.

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
  smallest missing suffix) and leaves the original intact. Purpose: output
  directory preparation and suffix allocation.
- `test_prepare_output_directory_creates_missing_parents` — Inputs: a nested
  path. Expected outputs: the full directory chain is created. Purpose: parent
  creation.

### tests/test_introgression_mapper.py

Consolidation outputs, count aggregation, and plot rendering.

- `test_generate_introgression_maps_creates_expected_outputs` — Inputs: synthetic
  triplet results and a species tree. Expected outputs: the combined PNG and all
  three TSVs are written. Purpose: artifact generation.
- `test_collect_counts_correct_avg_in_generate_introgression_maps` — Inputs:
  classified triplet results. Expected outputs: averages match the documented
  co-occurrence denominators. Purpose: bootstrap averaging.
- `test_collect_counts_non_ghost_denominator_is_all_co_occurring_triplets` —
  Inputs: directed-pair results. Expected outputs: the denominator counts every
  triplet containing both taxa. Purpose: population-level normalization.
- `test_collect_counts_ghost_denominator_is_all_triplets_containing_taxon` —
  Inputs: ghost-classified results. Expected outputs: the denominator counts
  every triplet containing the target taxon. Purpose: ghost normalization.
- `test_collect_non_sister_counts_counts_non_sister_pairs` — Inputs: results with
  non-sister pairs. Expected outputs: only non-sister pairs are counted.
  Purpose: pair selection.
- `test_generate_introgression_maps_excludes_outgroups` — Inputs: an outgroup
  list. Expected outputs: outgroup taxa are absent from plots and TSVs. Purpose:
  outgroup exclusion.
- `test_generate_introgression_maps_prunes_requested_plot_taxa` /
  `test_generate_introgression_maps_uses_full_species_tree_by_default` — Inputs:
  with and without `plot_taxa`. Expected outputs: the plotted tree is pruned or
  left full. Purpose: plot taxa selection.
- `test_generate_introgression_maps_uses_raw_values_with_separate_scales` —
  Inputs: results spanning a value range. Expected outputs: raw values with
  per-plot scales. Purpose: colour scaling.
- `test_generate_introgression_maps_appends_suffix_when_overwrite_disabled` —
  Inputs: an existing output directory with `overwrite=False`. Expected outputs:
  a suffixed sibling directory. Purpose: overwrite behavior.
- `test_generate_introgression_maps_preserves_run_dir_when_reset_disabled` —
  Inputs: `reset_output_dir=False`. Expected outputs: pre-existing run files
  survive. Purpose: the orchestrator's consolidation contract.
- `test_draw_species_tree_strip_shows_leaf_labels_by_default` /
  `test_draw_species_tree_strip_suppresses_leaf_labels` — Inputs: the tree strip
  renderer with and without label suppression. Expected outputs: labels present
  or absent. Purpose: plot layout.
- `test_scaled_consolidation_text_sizes_grow_with_taxa_count` — Inputs: taxa
  counts across a range. Expected outputs: text sizes scale and stay capped.
  Purpose: readability on large figures.
- `test_sampled_introgression_presence_flags_targets_with_sampled_edges` —
  Inputs: a taxa order and a `(source, target)` weight map with one zero-weight
  edge. Expected outputs: `{"A": 1, "B": 0, "C": 1, "D": 0}` — only taxa that
  are the target of a non-zero sampled edge are flagged. Purpose: the flag that
  drives ghost bar colour.
- `test_zero_heatmap_cells_are_masked` — Inputs: the ghost-colour scenario, with
  `sns.heatmap` monkeypatched to capture its `mask` argument. Expected outputs:
  the mask equals `data == 0` elementwise and exactly one cell is unmasked (the
  single sampled edge `(C, A)`). Purpose: empty cells are left unpainted rather
  than drawn at the colormap's low end.
- `test_ghost_strength_tsv_records_sampled_introgression_flag` — Inputs: results
  where taxon A has both ghost and sampled introgression and taxon D has ghost
  only. Expected outputs: the ghost TSV header is
  `target_taxon / raw_strength / has_sampled_introgression`, with `A → 1` and
  `D → 0`, and the strengths are unchanged by the flag. Purpose: the new column
  and its cross-referencing against the sampled sheet.
- `test_ghost_bars_use_constant_colours_by_sampled_presence` — Inputs: the same
  results, with `Axes.barh` monkeypatched to capture the colours actually
  passed. Expected outputs: only the two constants are used, `A` is
  `GHOST_WITH_SAMPLED_BAR_COLOR` (cividis low end) and `D` is
  `GHOST_ONLY_BAR_COLOR` (cividis high end), while their bar widths differ.
  Purpose: colour encodes co-occurrence, not magnitude.

### tests/test_ml_labels_and_metrics.py

The ML label contract, evaluation metrics, distributions, and CV-fold policy.

- `test_bit_labels_define_a_six_bit_contract` — Inputs: the module constants.
  Expected outputs: `BIT_COUNT == 6` with six unique labels. Purpose: pins the
  label contract.
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

- `test_load_ml_config_defaults_target_column_to_class` — Inputs: a config
  omitting `target_column`. Expected outputs: it defaults to `class`.
- `test_load_ml_config_accepts_explicit_class_target_column` — Inputs: an
  explicit `target_column`. Expected outputs: it is honored.
- `test_load_ml_config_defaults_min_samples_parameters` — Inputs: a config
  omitting the min-samples keys. Expected outputs: the documented defaults.
- `test_load_ml_config_honors_overwrite_flag` — Inputs: `overwrite: false`.
  Expected outputs: the flag is carried into the resolved config.

### tests/test_ml_utils.py

- `test_rows_to_matrix_uses_numeric_features_and_excludes_target_column` —
  Inputs: rows with numeric features plus the target column. Expected outputs: a
  numeric matrix excluding the target.
- `test_rows_to_matrix_encodes_multiple_string_columns` — Inputs: rows with
  low-cardinality string columns. Expected outputs: one-hot encoded features.
- `test_rows_to_matrix_rejects_string_features` — Inputs: a string column
  exceeding the cardinality limit. Expected outputs: an error rather than an
  invented ordering.

### tests/test_ml_random_forest.py

- `test_parse_classes_returns_binary_matrix` — Inputs: label strings. Expected
  outputs: the binary target matrix.
- `test_train_random_forest_smoke` — Inputs: `summary_statistics_tsv`. Expected
  outputs: training completes and writes its artifacts.
- `test_train_random_forest_creates_bitwise_metrics_report` — Inputs: the same
  fixture. Expected outputs: the metrics report contains the bitwise section.
- `test_metrics_txt_leads_with_hyperparameters` — Inputs: the same fixture with
  `n_estimators=25`, `random_state=7`, `test_size=0.25`, `max_depth=None`,
  `cv_folds=3`. Expected outputs: `Hyperparameters:` appears before
  `Test metrics:`; the parsed block reports those configured values with
  `max_depth` as `none`; and `metrics["hyperparameters"]` carries the same
  values with `cv_folds_requested == 3`. Purpose: the run's hyperparameters are
  recorded in both the text and JSON reports.

### tests/test_ml_multi_knn.py

- `test_multi_knn_train_smoke` — Inputs: `summary_statistics_tsv`. Expected
  outputs: training completes and writes its artifacts.
- `test_multi_knn_build_model_caps_neighbors_to_training_size` — Inputs: a
  configured `n_neighbors` larger than the training set. Expected outputs: the
  effective neighbor count is capped.
- `test_multi_knn_metrics_report_mentions_effective_neighbors` — Inputs: the same
  capped run. Expected outputs: the report carries a `Hyperparameters:` block
  naming both `n_neighbors_requested` and `n_neighbors_effective`, and
  `metrics["hyperparameters"]` reports 20 requested with at most 20 effective.
  Purpose: the cap is visible in both reports.

### tests/test_ml_hyper_tune.py

- `test_load_hyper_tune_config_accepts_hyperparameter_tuning_section` — Inputs: a
  tuning config. Expected outputs: the section loads.
- `test_load_hyper_tune_config_fills_model_defaults` — Inputs: a config omitting
  model parameters. Expected outputs: trainer defaults are filled in.
- `test_load_hyper_tune_config_accepts_wandb_detailed_payloads` — Inputs:
  `wandb_detailed_payloads: true`. Expected outputs: the flag is honored.
- `test_load_hyper_tune_config_rejects_evaluation_section` — Inputs: a tuning
  config containing `evaluation`. Expected outputs: `ConfigError`.
- `test_load_hyper_tune_config_requires_hyperparameter_tuning_section` — Inputs:
  a config without the section. Expected outputs: `ConfigError`.
- `test_tune_hyperparameters_grid_search_smoke` /
  `test_tune_hyperparameters_random_search_smoke` — Inputs:
  `summary_statistics_tsv_tuning` with each search method. Expected outputs: the
  search completes and reports ranked candidates.

## Parity Tests (`@pytest.mark.parity`)

Run with `pytest -m parity`. Only the two DendroPy-vs-BioPython tests in
`tests/orchestrator/test_orchestrator_tree_parity.py` carry this marker; every other
test derives its expectations from definitions instead of comparing against a
second implementation.
