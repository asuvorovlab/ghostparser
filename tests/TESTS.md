# Test Suite Documentation

This is the test map: every test function, the inputs and fixtures it uses, and
the behavior it asserts. For the concrete input values and a step-by-step
derivation of each expected result, see the companion [TEST_IO.md](TEST_IO.md).

## Running Tests

```bash
# Everything
pytest

# One file / one test
pytest tests/pipeline/test_pipeline_inference.py
pytest tests/pipeline/test_pipeline_decision.py::test_classify_ghost_when_discordant_heights_exceed_concordant

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
suite by `tests/conftest.py`. Pipeline-specific fixtures live in
`tests/pipeline/conftest.py`.

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
- `pipeline_species_tree` (pipeline) — Inputs: the 5-taxon species tree
  `(((A,B),C),(D,OUT))`. Expected usage: pipeline preprocessing and end-to-end
  runs.
- `pipeline_gene_trees` (pipeline) — Inputs: 12 gene trees, each containing
  `A,B,C,D,OUT` so rooting always succeeds. Expected usage: pipeline
  preprocessing and end-to-end runs.

## Function-Level Coverage

### tests/pipeline/test_pipeline.py

End-to-end `run_pipeline` behavior on the shared 5-taxon / 12-gene-tree fixture.

- `test_run_pipeline_matches_derived_expectation` — Inputs: `run_pipeline`
  (serial, `taxon` mode, fixed bootstrap seed, consolidation off). Expected
  outputs: all 4 triplets present with the hand-derived topology counts
  (7/3/2 for `(A,B,C)`; 12/0/0 for each D-containing triplet), the SciPy
  chi-square DCT statistic and p-value, Bonferroni-corrected p-values, and
  `no_introgression` for every triplet. Purpose: end-to-end correctness against
  values derived from the fixture rather than another module.
- `test_run_pipeline_writes_results_tsv` — Inputs: the same serial run.
  Expected outputs: `pipeline_triplet_results.tsv` exists, its header starts
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

### tests/pipeline/test_pipeline_inference.py

Per-triplet inference on a 10-gene-subtree fixture, with expectations recomputed
from the tabulated tree geometry.

- `test_analyze_triplet_matches_derived_expectation` — Inputs (parametrized over
  6 tree-height strategies x 3 summary statistics x 2 discordant tests = 36
  cases): `analyze_triplet` over the 10 gene subtrees. Expected outputs: the
  counts, DCT/KS statistics, summary values, and classification all equal values
  derived in-test from `_LEAF_GEOMETRY` plus direct SciPy/statsmodels calls;
  bootstrap fractions sum to 1. Purpose: the full inference surface across every
  parameter combination.
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

### tests/pipeline/test_pipeline_decision.py

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
  significant split with fully separated heights, concordant above discordant.
  Expected outputs: KS statistic 1.0 and significant, `summary_con >
  summary_dis`, `outflow_introgression`. Purpose: gate 3, con > dis.
- `test_classify_ghost_when_discordant_heights_exceed_concordant` — Inputs: the
  mirror case, discordant above concordant. Expected outputs:
  `summary_con < summary_dis`, `ghost_introgression`. Purpose: gate 3, con < dis.
- `test_classify_introgression_truth_table` — Inputs (parametrized, 8 rows):
  every combination of DCT/KS significance and summary ordering, including
  `None` summaries. Expected outputs: the documented classification for each
  row. Purpose: exhaustive coverage of `_classify_introgression`.
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

### tests/pipeline/test_pipeline_trees.py

Tree preprocessing, asserted against explicit Newick literals.

- `test_clean_and_save_trees_preserves_a_well_supported_tree` — Inputs:
  `pipeline_species_tree` with `min_avg_support=0.5`. Expected outputs: the
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
  `pipeline_gene_trees` with outgroup `OUT`. Expected outputs: all 12 trees
  survive, each ends in `OUT:0);`, and trees 0 and 3 match their expected
  rerooted Newick (the ingroup edge absorbs OUT's original edge length).
  Purpose: gene-tree rooting semantics.
- `test_extract_triplet_subtree_selects_the_triplet_taxa` — Inputs: the first
  cleaned gene tree and triplet `(A, B, C)`. Expected outputs: a subtree whose
  leaf set is exactly `{A, B, C}`. Purpose: subtree extraction.

### tests/pipeline/test_pipeline_tree_parity.py

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

### tests/pipeline/test_pipeline_config.py

Pipeline config resolution and config-file precedence.

- `test_cli_defaults_resolve` — Inputs: a CLI namespace with every optional arg
  `None`. Expected outputs: `alpha_dct`/`alpha_ks` 0.05, the pipeline-specific
  `p_value_correction == "bfn"` and `summary_statistic == "mean"`,
  `overwrite is True`, the config-file-only keys at their defaults, and no
  `stats_backend` key. Purpose: default resolution in CLI mode.
- `test_cli_overrides_for_config_plus_cli_options` — Inputs: CLI values for
  alpha-dct/alpha-ks/summary-statistic/p-value-correction/no-overwrite. Expected
  outputs: each override is honored and `overwrite` becomes `False`. Purpose:
  the config+CLI options are wired.
- `test_config_only_keys_read_from_config_file` — Inputs: a JSON config setting
  `discordant_test`, `tree_height_calculation_strategy`, `min_support_value`,
  `generate_summary_stats`, `alpha_dct`, and a nested `bootstrap_options` block.
  Expected outputs: every key, including the flattened bootstrap options, is
  read. Purpose: config-file-only keys and nested bootstrap parsing.
- `test_config_file_wins_over_cli` — Inputs: a config file plus conflicting CLI
  flags. Expected outputs: the file's `alpha_dct` wins, `summary_statistic`
  falls back to the pipeline default (proving the CLI value was ignored), the
  file's paths are used, and a warning is printed. Purpose: config-file
  precedence.
- `test_missing_required_field_raises` — Inputs: a namespace missing the species
  tree. Expected outputs: `ConfigError`. Purpose: required-field validation.
- `test_parser_exposes_config_file_and_new_flags` — Inputs: an argv list using
  the config+CLI flags. Expected outputs: each parses to its expected value and
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
  survive. Purpose: the pipeline's consolidation contract.
- `test_draw_species_tree_strip_shows_leaf_labels_by_default` /
  `test_draw_species_tree_strip_suppresses_leaf_labels` — Inputs: the tree strip
  renderer with and without label suppression. Expected outputs: labels present
  or absent. Purpose: plot layout.
- `test_scaled_consolidation_text_sizes_grow_with_taxa_count` — Inputs: taxa
  counts across a range. Expected outputs: text sizes scale and stay capped.
  Purpose: readability on large figures.

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

### tests/test_ml_multi_knn.py

- `test_multi_knn_train_smoke` — Inputs: `summary_statistics_tsv`. Expected
  outputs: training completes and writes its artifacts.
- `test_multi_knn_build_model_caps_neighbors_to_training_size` — Inputs: a
  configured `n_neighbors` larger than the training set. Expected outputs: the
  effective neighbor count is capped.
- `test_multi_knn_metrics_report_mentions_effective_neighbors` — Inputs: the same
  capped run. Expected outputs: the report states the effective neighbor count.

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
`tests/pipeline/test_pipeline_tree_parity.py` carry this marker; every other
test derives its expectations from definitions instead of comparing against a
second implementation.
