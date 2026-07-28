# Test Suite Documentation

This document lists the current tests, fixtures, and marker-based slices, with function-level input/output expectations.
Each test is documented either as an individual entry or inside a grouped entry that provides inputs, expected behavior, and test purpose for that set of tests.

## Running Tests

Run all tests:

```bash
pytest
```

Run one file:

```bash
pytest tests/test_triplet_processor.py
```

Run one function:

```bash
pytest tests/test_triplet_processor.py::test_adjust_p_values_standard_matches_statsmodels_for_supported_methods
```

Run backend parity tests only:

```bash
pytest -m backend_parity
```

Run only config tests:

```bash
pytest tests/test_config.py
```

## Fixtures In Use

Shared fixtures are exported via `tests/conftest.py` from `tests/fixtures.py`.

- `simple_newick_file`
Inputs: one rooted Newick with branch lengths and one outgroup.
Expected usage: tree read/format/clean path, no support labels present.

- `newick_with_support_file`
Inputs: Newick containing internal support labels.
Expected usage: support extraction/removal and support-aware cleaning checks.

- `multiple_trees_file`
Inputs: file containing 3 trees.
Expected usage: multi-tree read/write and count-preservation checks.

- `low_support_tree_file`
Inputs: mixed high-support and low-support trees.
Expected usage: support-threshold filtering checks.

- `simple_species_tree`, `simple_gene_trees`
Inputs: minimal species/gene trees for integration-style parser workflows.
Expected usage: end-to-end triplet extraction paths.

- `triplet_comparison_cases`
Inputs: `(newick, triplet)` cases for DendroPy vs BioPython consistency checks.
Expected usage: numeric branch-length/partristic-distance parity checks.

Local fixture in `tests/test_tree_parser.py`:

- `gene_trees_missing_outgroup_file`
Inputs: one valid gene tree plus one missing required outgroup.
Expected output: missing-outgroup tree is discarded during cleaning.

## Function-Level Coverage

### `tests/test_config.py`

- `test_load_orchestrator_config_json`
Inputs: full orchestrator JSON config (paths, methods, thresholds, bootstrap settings).
Expected outputs: normalized absolute paths, parsed outgroups list, preserved explicit values, and default parquet-orchestrator settings.

- `test_load_orchestrator_config_yaml`
Inputs: minimal YAML config with comma-separated outgroup string.
Expected outputs: `outgroup` normalized to list.

- `test_load_orchestrator_config_single_outgroup_string_is_single_taxon`
Inputs: single outgroup string.
Expected outputs: one-element outgroup list.

- `test_load_orchestrator_config_missing_required`
Inputs: config missing required fields.
Expected outputs: `ConfigError` with missing-field message.

- `test_load_orchestrator_config_invalid_choice_fields` (parametrized)
Inputs: invalid values for `discordant_test`, `summary_statistic`, `stats_backend`.
Expected outputs: `ConfigError` referencing offending field.

- `test_load_orchestrator_config_tree_height_strategy_validation` (parametrized)
Inputs: supported values (`AVG`, `A`, `B`, `C`, `SIS`, `INT`) and an invalid value (`D`).
Expected outputs: supported values load successfully; invalid value raises `ConfigError`.

- `test_load_orchestrator_config_defaults_processes_to_zero`
Inputs: no `processes` key.
Expected outputs: `processes == 0`.

- `test_load_orchestrator_config_honors_overwrite_flag`
Inputs: orchestrator JSON config with `overwrite: false`.
Expected outputs: normalized config preserves `overwrite == False`.

- `test_prepare_output_directory_overwrites_or_suffixes`
Inputs: an existing results directory with stale files, plus a second existing target path.
Expected outputs: the first directory is cleared when overwrite is enabled and the second call returns a suffixed directory path.

- `test_load_orchestrator_config_allows_disabling_consolidation`
Inputs: config with `consolidation: false`.
Expected outputs: normalized orchestrator config preserves `consolidation == False`.

- Only orchestrator config loading is covered in `tests/test_config.py`.

- `test_path_resolution_for_absolute_relative_and_home_paths` (parametrized)
Inputs: absolute paths, relative paths, and `~` paths.
Expected outputs: all normalized to resolved absolute paths.

### `tests/test_orchestrator.py`

- `test_resolve_processes_zero_uses_all_cores`
Inputs: `processes` in `{0, None, 4}` with monkeypatched CPU count.
Expected outputs: `0 -> cpu_count`, `None -> None`, explicit value preserved.

- `test_resolve_parallel_mode`
Inputs: `{0, 1, 4}` with monkeypatched CPU count.
Expected outputs: `(processes, use_multiprocessing)` toggles correctly (`1` disables multiprocessing).

- `test_resolve_runtime_args_cli_defaults_and_overrides` (parametrized)
Inputs: CLI defaults and a CLI override scenario (`processes`, bootstrap options).
Expected outputs: default statistical settings/path resolution, default parquet output settings, default `consolidation=True`, and preserved CLI overrides.

- `test_resolve_runtime_args_config_with_cli_warns_and_ignores`
Inputs: config file + conflicting CLI args.
Expected outputs: warning emitted; config values/defaults win over CLI extras, including default parquet orchestrator settings.

- `test_resolve_runtime_args_config_processes_behavior` (parametrized)
Inputs: config with explicit `processes` and config without `processes`.
Expected outputs: preserves configured value or defaults to `0`.

- `test_write_triplet_gene_trees_parquet_multiprocess`
Inputs: one triplet, a small set of gene-tree Newick strings, and a species-triplet map.
Expected outputs: parquet dataset directory is created with `triplets/` and `observations/` subdirectories, parquet part files are written, and the returned counts reflect extracted subtrees.

- `test_parse_and_analyze_triplet_gene_trees_parquet`
Inputs: parquet triplet dataset created from one triplet and three gene trees.
Expected outputs: parquet parser returns a populated triplet entry with cached observation rows, and analysis succeeds when `input_format="parquet"` is selected.

### `tests/test_tree_parser.py`

Runtime-arg resolution:

- `test_resolve_runtime_args_tree_parser_cli_processes` (parametrized)
Inputs: CLI defaults and explicit `processes=6`.
Expected outputs: path resolution, outgroup parsing, default `min_support_value`, and expected process behavior.

- Tree-parser runtime arg tests cover CLI defaults and overrides.

Tree IO and cleaning:

- `test_read_tree_file_single_tree`, `test_read_tree_file_multiple_trees`
Inputs: one-tree and multi-tree files.
Expected outputs: correct tree counts and tree object types.

- `test_read_tree_file_not_found`
Inputs: nonexistent file path.
Expected outputs: `FileNotFoundError`.

- `test_read_tree_file_invalid_inputs_raise_value_error` (parametrized)
Inputs: malformed Newick, random invalid text, empty file.
Expected outputs: `ValueError("Invalid Newick format")`.

- `test_calculate_average_support_with_values` / `test_calculate_average_support_no_values`
Inputs: tree with support labels vs tree without labels.
Expected outputs: numeric average support vs `None`.

- `test_remove_support_values`, `test_standardize_tree_removes_support`
Inputs: supported tree.
Expected outputs: supports removed (`None` average support afterwards).

- `test_standardize_tree_preserves_branch_lengths`
Inputs: branch-length tree.
Expected outputs: branch lengths unchanged within tolerance.

Triplet generation/writer/integration tests:

- Covers `generate_triplets`, filter parsing, taxa filtering, and triplet writing.
Inputs: taxa sets, outgroup forms (single/list/comma-separated), optional filter files.
Expected outputs: correct triplet counts/content, outgroup exclusion, valid sectioned output formats.

- Covers single-process, streaming, and multiprocess triplet writers.
Inputs: triplet collections with varying sizes/empties.
Expected outputs: stable output shape, correct separators/headers, no crashes on edge cases.

- Cross-library consistency tests (`test_triplet_branch_lengths_match`, `test_triplet_collapse_consistency_dendropy_vs_biopython`).
Inputs: fixture case set from `triplet_comparison_cases`.
Expected outputs: pairwise distances agree across implementations.

### `tests/test_triplet_processor.py`

Core statistic helpers:

- `test_compute_tree_height_statistic_matches_definition`
Inputs: known tree with explicit branch lengths.
Expected outputs: exact formula match for AVG strategy.

- `test_compute_tree_height_statistic_supports_extended_strategies` (parametrized)
Inputs: rooted triplets with expected values for `A`, `B`, `C`, `SIS`, and `INT`.
Expected outputs: each strategy returns its expected branch-length-based statistic.

- `test_compute_tree_height_statistic_rejects_unknown_strategy`
Inputs: invalid strategy.
Expected outputs: `ValueError`.

- `test_compute_tree_height_statistic_requires_species_triplet_for_taxon_specific_strategies`
Inputs: taxon-specific strategy without `species_triplet`.
Expected outputs: `ValueError`.

Topology classification and pipeline behavior:

- Tests cover concordant/discordant label mapping, tie behavior, deterministic discordant role assignment, and the 5 output classes (`no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, `unresolved` where applicable).
Inputs: controlled synthetic topology distributions and tree-height profiles.
Expected outputs: deterministic role counts, significance states, and final classification strings.

- `test_run_triplet_pipeline_bootstrap_unresolved_when_metrics_missing`
Inputs: concordant-only synthetic trees with bootstrap enabled.
Expected outputs: bootstrap unresolved fraction is 1.0 and `bootstrap_value` is 1.0.

- `test_run_bootstrap_iterations_joins_tied_classes`
Inputs: monkeypatched per-iteration classifications split evenly across two classes.
Expected outputs: bootstrap support equals the tied fraction (`0.5`).

Parser/writer behavior:

- Roundtrip parsing/writing and dynamic column tests.
Inputs: sectioned `unique_triplets_gene_trees.txt` test content and generated pipeline results.
Expected outputs: header validation, dynamic summary columns, dynamic corrected columns (`dct_p_val_<method>_corr`, `ks_p_val_<method>_corr`), required error paths for unsupported settings.

- `test_write_pipeline_results_adds_bootstrap_columns_when_enabled`
Inputs: bootstrap-enabled triplet result written to TSV.
Expected outputs: base bootstrap columns (`bootstrap_value`, `all_bootstrap`) are present, debug bootstrap columns are present when debug mode is enabled. The `all_bootstrap` cell is formatted as comma-separated `classification=value` pairs (e.g. `no_introgression=0.5,ghost_introgression=0.5`). Debug array columns (e.g. `bootstrap_gene_tree_heights`) remain JSON-serialized.

- `test_write_pipeline_results_adds_bootstrap_gene_tree_heights_when_summary_only_false`
Inputs: bootstrap-enabled result written with debug mode enabled and `summary_only=false`.
Expected outputs: TSV includes `bootstrap_gene_tree_heights` and stores the raw per-triplet tree-height list.

### `tests/pipeline/`

Shared fixtures (`tests/pipeline/conftest.py`): `pipeline_species_tree` and `pipeline_gene_trees` write a small 5-taxon species tree and a set of gene trees that all carry the `OUT` outgroup so rooting produces real observations.

`tests/pipeline/test_pipeline_inference.py`:

- `test_analyze_triplet_matches_analyze_triplet_entry`
Inputs: one triplet with fixed gene subtrees and species subtree, parametrized across `discordant_test` (chi-square/z-test), `summary_statistic` (mean/median/mode), `stats_backend` (custom/standard), and tree-height `strategy` (AVG/A/B/C/SIS/INT), with a fixed bootstrap seed.
Expected outputs: `pipeline.inference.analyze_triplet` equals `triplet_processor.analyze_triplet_entry` on every non-bootstrap field (counts, DCT/KS statistics and p-values, classification, summaries), and the pipeline result's bootstrap aggregates are well formed. Bootstrap values are excluded from the equality because the pipeline resamples with NumPy. Purpose: prove the ported inference math stays bit-for-bit equal.

- `test_analyze_triplet_from_observations_matches_newick_path`
Inputs: observations serialized from fixed gene subtrees, and the same subtrees as Newick, with a fixed seed.
Expected outputs: `analyze_triplet_from_observations` equals `analyze_triplet` on all non-bootstrap fields, and (same observations + same seed) identical bootstrap aggregates.

- `test_bootstrap_is_deterministic_under_seed`
Inputs: two `analyze_triplet` calls with the same fixed seed.
Expected outputs: identical `all_bootstrap` and `bootstrap_value`. Purpose: the NumPy bootstrap is reproducible under a seed.

- `test_analyze_triplet_empty_observations`
Inputs: a triplet with zero gene subtrees.
Expected outputs: `analyzed_trees == 0`, `n_con == 0`, classification `no_introgression`.

`tests/pipeline/test_pipeline_trees.py`:

- `test_species_preprocessing_matches_tree_parser`
Inputs: the `pipeline_species_tree` fixture.
Expected outputs: pipeline cleaning, outgroup rooting/pruning, ingroup set, normalized A/B/C triplets, and species-subtree map all equal `tree_parser`'s.

- `test_gene_tree_cleaning_matches_tree_parser`
Inputs: the `pipeline_gene_trees` fixture.
Expected outputs: pipeline gene-tree cleaning/rooting output file and parsed Newick list equal `tree_parser`'s.

`tests/pipeline/test_pipeline.py`:

- `test_run_pipeline_matches_reference`
Inputs: `run_pipeline` (serial) on the fixtures with a fixed bootstrap seed and consolidation disabled.
Expected outputs: per-triplet results equal a reference built independently from `tree_parser` extraction plus `triplet_processor.analyze_triplet_entry` and run-wide correction. Purpose: end-to-end parity with the orchestrator path.

- `test_run_pipeline_writes_results_tsv`
Inputs: a `run_pipeline` run.
Expected outputs: `pipeline_triplet_results.tsv` exists with the expected header (`triplet`, `classification`, `bootstrap_value`) and one row per result.

- `test_consolidation_preserves_run_outputs`
Inputs: a `run_pipeline` run with consolidation enabled.
Expected outputs: the run folder still contains `pipeline_triplet_results.tsv`, `metrics.txt`, and both processed trees, and the consolidation artifacts land in a non-empty `consolidation/` subfolder. Purpose: regression guard that consolidation's output-directory reset does not delete the run's primary outputs.

- `test_parallel_modes_match_serial`
Inputs: `run_pipeline` in `taxon` and `gene` modes with two workers, versus a serial run.
Expected outputs: identical per-triplet results across modes (bootstrap is deterministic per triplet under a fixed seed).

### `tests/test_ml_config.py`

- `test_load_ml_config_defaults_target_column_to_class`
Inputs: minimal ML config containing only `input_path` and `output_dir`.
Expected outputs: `target_column` defaults to `class` and `overwrite` defaults to `True`.

- `test_load_ml_config_accepts_explicit_class_target_column`
Inputs: minimal ML config with `target_column: class` and one model parameter.
Expected outputs: explicit `class` target is preserved and model parameters load correctly.

- `test_load_ml_config_defaults_min_samples_parameters`
Inputs: minimal ML config without `min_samples_split` or `min_samples_leaf`.
Expected outputs: the loader fills in the default RF values for both keys.

- `test_load_ml_config_honors_overwrite_flag`
Inputs: ML config with `overwrite: false`.
Expected outputs: normalized config preserves `overwrite == False`.

### `tests/test_ml_utils.py`

- `test_rows_to_matrix_uses_numeric_features_and_excludes_target_column`
Inputs: TSV-like rows containing only numeric feature columns and the `class` target bitstring.
Expected outputs: the target column is excluded from features, numeric columns are used directly, and the label bitstrings are preserved for training.

- `test_rows_to_matrix_rejects_string_features`
Inputs: TSV-like rows containing one string-valued feature column and the `class` target bitstring.
Expected outputs: `ValueError` indicating that the non-target feature must be numeric.

P-value correction behavior:

- `test_adjust_p_values_custom_fdr_matches_known_bh_example`
Inputs: fixed BH example p-values.
Expected outputs: known corrected values.

- `test_adjust_p_values_standard_matches_statsmodels_for_supported_methods`
Inputs: same p-values for `bfn`, `holm`, `fdr_bh`, `fdr_by`, `fdr_tsbh`.
Expected outputs: exact match to `statsmodels.multipletests` outputs.

- `test_adjust_p_values_custom_matches_standard_randomized` (parametrized)
Inputs: randomized p-values across correction methods and optional alpha.
Expected outputs: custom backend equals standard backend within tight tolerance.

Runtime-arg resolution:

- Covers CLI defaults and `processes` default/preservation behavior.
Inputs: CLI args.
Expected outputs: resolved defaults and expected process semantics.

### Explicit Grouped Test Names

This addendum lists tests that are intentionally grouped in the narrative sections above and named explicitly.

#### tests/test_orchestrator.py

- Tests: `test_main_uses_file_backed_pipeline`
Inputs: orchestrator runtime with normalized species/gene trees and triplet metadata.
Expected outputs/behavior: orchestrator always runs file-backed triplet extraction to `unique_triplets_gene_trees.txt`, then runs inference from that file, and runs consolidation stage by default.
Purpose: verify the orchestrator executes the canonical two-stage file-backed pipeline.

- Tests: `test_main_skips_consolidation_stage_when_disabled`
Inputs: orchestrator runtime with `consolidation=False`.
Expected outputs/behavior: extraction and inference still run, while consolidation/map generation is skipped.
Purpose: verify configuration-controlled enable/disable behavior for consolidation artifacts.

#### tests/test_introgression_mapper.py

- Tests: `test_generate_introgression_maps_creates_expected_outputs`, `test_generate_introgression_maps_appends_suffix_when_overwrite_disabled`, `test_generate_introgression_maps_uses_full_species_tree_by_default`, `test_generate_introgression_maps_prunes_requested_plot_taxa`, `test_generate_introgression_maps_uses_raw_values_with_separate_scales`
  Inputs: synthetic triplet results with inflow/outflow/ghost classifications and bootstrap weights.
  Expected outputs/behavior: mapper writes expected plot/TSV artifacts (including `introgression_matrix_sampled_non_sister.tsv`), uses the full processed species tree by default, optionally prunes to requested plot taxa when supplied, average bootstrap values use population-level denominators, source taxon labels appear on top of the heatmap (between the tree strip and the heatmap cells), and the species tree strip is drawn above that.
  Purpose: validate consolidation artifact generation, plot layout semantics, and denominator correctness.

- Tests: `test_generate_introgression_maps_excludes_outgroups`
Inputs: results containing a triplet with a taxon designated as outgroup via the `outgroups` parameter.
Expected outputs/behavior: outgroup taxon is absent from the matrix TSV column headers, ghost strength TSV rows, and the reported `taxa_count`.
Purpose: verify that the `outgroups` parameter correctly filters taxa from all consolidation outputs.

- Tests: `test_generate_introgression_maps_preserves_run_dir_when_reset_disabled`
Inputs: an output directory containing a pre-existing `orchestrator_triplet_results.tsv`, with `reset_output_dir=False`.
Expected outputs/behavior: the pre-existing file survives and the combined plot is written into the same directory (no reset, no suffixed sibling).
Purpose: regression guard that consolidation invoked by the orchestrator/pipeline does not delete the caller's already-written run outputs.

- Tests: `test_collect_counts_non_ghost_denominator_is_all_co_occurring_triplets`
Inputs: three synthetic results — one classified inflow, one no_introgression, one unrelated triplet (ABD).
Expected outputs/behavior: `non_ghost_counts[(C, B)]` equals 2 (both ABC rows, regardless of classification); `non_ghost_counts[(B, D)]` equals 1 (only ABD).
Purpose: verify that the non-ghost denominator counts all triplets where both taxa co-occur, not just classified ones.

- Tests: `test_collect_counts_ghost_denominator_is_all_triplets_containing_taxon`
Inputs: three synthetic results across triplets (A,B,C) ×2 and (A,C,D) ×1.
Expected outputs/behavior: `ghost_counts[A]` = 3, `ghost_counts[C]` = 3, `ghost_counts[D]` = 1, `ghost_counts[B]` = 2.
Purpose: verify that the ghost denominator counts all triplets where a taxon appears in any position.

- Tests: `test_collect_counts_correct_avg_in_generate_introgression_maps`
Inputs: triplets (A,B,C) with one inflow (weight 0.6) and one no_introgression; triplet (A,B,D) with one ghost (weight 0.8).
Expected outputs/behavior: matrix TSV cell `B←C` = 0.3 (0.6/2); ghost TSV cell `A` ≈ 0.2667 (0.8/3).
Purpose: end-to-end verification that population-level denominators flow through to TSV output values.

- Tests: `test_draw_species_tree_strip_suppresses_leaf_labels`, `test_draw_species_tree_strip_shows_leaf_labels_by_default`
  Inputs: three-taxon species tree; `show_leaf_labels=False` vs default (`True`).
  Expected outputs/behavior: with `False`, no Text artists with taxon names appear on the axis; with default `True`, one Text artist per leaf taxon is present.
  Purpose: verify the `show_leaf_labels` parameter controls leaf annotation rendering on the tree strip axis.

- Tests: `test_collect_non_sister_counts_counts_non_sister_pairs`
  Inputs: two results for triplet (A,B,C) (inflow and no_introgression) and one result for triplet (A,C,D) (ghost).
  Expected outputs/behavior: `counts[(A,C)]` = 2; `counts[(B,C)]` = 2; `counts[(A,B)]` = 0 (sister pair, never incremented); `counts[(A,D)]` = 1; `counts[(C,D)]` = 1.
  Purpose: verify that `_collect_non_sister_counts` increments only non-sister pairs and accumulates counts across multiple triplet rows.

#### tests/test_tree_parser.py

- Tests: `test_write_clean_trees_outputs_expected_tree_count`, `test_clean_and_save_trees_filters_low_support`, `test_clean_and_save_trees_no_filters`, `test_clean_and_save_trees_creates_output_file`, `test_clean_and_save_gene_trees_discards_missing_outgroup`
Inputs: single/multiple trees, low-support trees, gene trees missing outgroup.
Expected outputs/behavior: clean outputs are written, support filtering behaves correctly, invalid/missing-outgroup trees are excluded where required.
Purpose: validate cleaned-tree persistence and support/outgroup filtering behavior.

- Tests: `test_get_taxa_from_tree_correct_names`, `test_generate_triplets_count`, `test_generate_triplets_excludes_outgroup`, `test_generate_triplets_content`, `test_generate_triplets_large_set`, `test_generate_triplets_multiple_outgroups`, `test_generate_triplets_outgroup_comma_separated_with_spaces`
Inputs: rooted species trees with varying taxa sets and outgroup forms.
Expected outputs/behavior: taxa extraction is correct; triplets are generated with correct count/content and outgroup exclusions.
Purpose: verify triplet generation semantics across small and larger taxa sets.

- Tests: `test_write_triplets_to_file`, `test_write_triplets_to_file_empty`, `test_read_triplet_filter_file_parses_valid_and_skips_invalid`, `test_filter_triplets_by_taxa_skips_missing_taxa`
Inputs: generated triplet collections, empty collections, valid/invalid triplet-filter file lines, taxa-subset filters.
Expected outputs/behavior: triplet files are written in expected format; empty handling is stable; filter parsing and taxa-based filtering are correct.
Purpose: validate triplet-file IO and filter utility behavior.

- Tests: `test_format_newick_with_precision_trailing_zeros`, `test_format_newick_with_precision_default_places`, `test_format_newick_with_custom_precision`, `test_format_newick_with_precision_triplet_parser`
Inputs: branch-length Newick trees with precision options.
Expected outputs/behavior: formatted Newick strings preserve intended precision and representation.
Purpose: ensure deterministic and configurable Newick formatting.

- Tests: `test_extract_triplet_subtree_all_taxa_present`, `test_extract_triplet_subtree_missing_taxa`, `test_extract_triplet_subtree_preserves_branch_lengths`, `test_process_gene_trees_for_triplets`, `test_process_gene_trees_for_triplets_empty`, `test_build_species_triplet_metadata_normalizes_abc`
Inputs: gene-tree triplet extraction requests with complete/missing taxa and species-triplet metadata setup.
Expected outputs/behavior: extraction succeeds only when all taxa are present, preserves branch lengths, handles empty cases, and normalizes species metadata to A/B/C conventions.
Purpose: validate extraction core and metadata normalization used by downstream inference.

- Tests: `test_write_triplet_gene_trees`, `test_write_triplet_gene_trees_includes_species_tree_header`, `test_write_triplet_gene_trees_empty_triplet`, `test_triplet_gene_trees_separator_format`, `test_write_triplet_gene_trees_streaming`
Inputs: triplet-to-gene-tree mappings in normal, empty, and streaming write modes.
Expected outputs/behavior: mapping file sections, species-tree header, and separators are correctly serialized.
Purpose: verify canonical serialization format for triplet gene-tree mapping output.

- Tests: `test_write_triplet_gene_trees_multiprocess_with_workers`, `test_write_triplet_gene_trees_multiprocess_includes_species_header`, `test_multiprocessing_triplet_writer_handles_empty_triplets`, `test_write_triplet_gene_trees_multiprocess_triplets_single_worker`, `test_write_triplet_gene_trees_multiprocess_accepts_list`
Inputs: multiprocess writer invocations across worker-count and input-shape variants.
Expected outputs/behavior: output format remains valid; species header persists; empty and list-based inputs are handled safely; optional worker CPU telemetry (`return_worker_cpu=True`) returns a non-negative CPU-seconds value.
Purpose: validate robust multiprocess mapping-file writer behavior.

- Tests: `test_get_clean_filename_variants`, `test_metrics_logger_context_manager`, `test_metrics_logger_file_not_opened_before_enter`
Inputs: filename variants and metrics-logger lifecycle usage.
Expected outputs/behavior: cleaned output filenames are formed correctly; logger opens/writes only in expected context-manager lifecycle.
Purpose: verify utility helpers that support parser CLI workflows.

- Tests: `test_integration_full_workflow`, `test_integration_triplets_workflow`, `test_integration_full_triplet_extraction_workflow`
Inputs: integration-style species/gene tree fixtures and output destinations.
Expected outputs/behavior: end-to-end parsing, triplet generation/extraction, and file outputs complete successfully.
Purpose: ensure combined parser workflow remains functional.

#### tests/test_triplet_processor.py

- Tests: `test_classify_triplet_topology_string_for_all_three_topologies`, `test_classify_triplet_topology_labels_concordant_and_discordants`, `test_balanced_discordant_count_tests_are_not_significant`
Inputs: representative topology strings and balanced discordant count scenarios.
Expected outputs/behavior: topology labels map correctly and balanced discordant tests remain non-significant.
Purpose: validate baseline topology classification and count-test behavior.

- Tests: `test_run_triplet_pipeline_uses_species_concordant_and_frequency_ranked_discordants`, `test_run_triplet_pipeline_supports_z_test_for_discordant_counts`, `test_run_triplet_pipeline_supports_standard_stats_backend`, `test_run_triplet_pipeline_supports_median_summary_statistic`, `test_run_triplet_pipeline_supports_mode_summary_statistic`, `test_run_triplet_pipeline_supports_taxon_specific_tree_height_strategy`, `test_run_triplet_pipeline_breaks_discordant_ties_by_first_topology`, `test_run_triplet_pipeline_selects_ac_as_discordant1_when_ac_is_more_frequent`, `test_run_triplet_pipeline_no_introgression_when_dct_not_significant`, `test_run_triplet_pipeline_inflow_when_ks_not_significant`, `test_run_triplet_pipeline_outflow_when_con_summary_higher`, `test_run_triplet_pipeline_ghost_when_dis_summary_higher`
Inputs: synthetic per-triplet topology/tree-height distributions, configurable test/stat backends, and strategy variants.
Expected outputs/behavior: discordant role assignment, statistical backend selection, summary-stat selection, and final classification outcomes match expected logic.
Purpose: validate triplet inference decision logic across major branches.

- Tests: `test_analyze_triplet_gene_tree_file_with_multiprocessing`, `test_parse_analyze_and_write_pipeline_roundtrip_with_species_header`, `test_collect_triplet_statistics_returns_dict_list`
Inputs: mapping files and pipeline run settings, including multiprocessing.
Expected outputs/behavior: analyze/parse/write pipeline roundtrips successfully, multiprocessing analysis can return optional worker CPU telemetry (`return_worker_cpu=True`) with a non-negative value, and statistics collection returns expected dictionary-list structures.
Purpose: validate end-to-end processing API behavior.

- Tests: `test_analyze_triplet_gene_tree_file_rejects_unsupported_runtime_options`, `test_parse_triplet_gene_trees_file_rejects_malformed_sections`
Inputs: invalid configuration values and malformed mapping-file headers/content.
Expected outputs/behavior: parser/analyzer rejects invalid inputs with explicit error paths.
Purpose: verify input validation and defensive error handling.

- Tests: `test_analyze_triplet_gene_tree_file_applies_selected_correction`, `test_two_sample_ks_test_hybrid_uses_scipy_near_threshold`, `test_two_sample_ks_test_hybrid_keeps_custom_when_not_borderline`, `test_two_sample_ks_test_hybrid_rejects_negative_margin`
Inputs: p-value correction selections and KS hybrid-mode threshold conditions.
Expected outputs/behavior: selected correction is applied; KS hybrid dispatches to expected backend and validates margin constraints.
Purpose: validate statistical-dispatch control flow.

- Tests: `test_write_pipeline_results_includes_dis1_topology_and_omits_removed_topology_columns`, `test_write_pipeline_results_uses_dynamic_summary_column_names`, `test_write_pipeline_results_uses_dct_chi_stats_column_for_chi_square`, `test_write_pipeline_results_uses_dct_z_score_column_for_z_test`, `test_write_pipeline_results_uses_mode_summary_columns_for_mode`, `test_write_pipeline_results_uses_dynamic_corrected_p_value_column_names`, `test_write_pipeline_results_includes_abc_mapping_column`, `test_write_pipeline_results_ghost_inference_uses_dis1_outgroup_recipient`, `test_write_pipeline_results_rejects_mixed_discordant_test_outputs`, `test_write_pipeline_results_rejects_unsupported_p_value_correction`, `test_write_pipeline_results_rejects_unsupported_summary_statistic`, `test_write_pipeline_statistics_json`, `test_serialize_bootstrap_value_rejects_non_json_value`
Inputs: synthetic pipeline result rows across discordant-test/summary-stat/correction settings and serialization targets.
Expected outputs/behavior: TSV/JSON outputs contain expected dynamic columns (including `inference`), enforce strict bootstrap JSON serialization, and reject unsupported or mixed output states.
Purpose: validate output-schema stability and writer safeguards.

- Tests: `test_run_triplet_pipeline_skips_summary_metric_collection_when_disabled`, `test_write_summary_statistics_tsv_includes_expected_columns_and_counts`, `test_write_summary_statistics_tsv_includes_bootstrap_value_when_enabled`
Inputs: pipeline runs with summary-stat metric collection disabled/enabled and summary-statistics payloads with bootstrap disabled/enabled.
Expected outputs/behavior: topology summary metrics are skipped when disabled; when enabled, summary TSV includes required 63-stat topology metrics plus identity/count/classification fields and bootstrap_value when enabled.
Purpose: verify summary-statistics gating and summary-statistics file schema/conditional bootstrap column behavior.

- Tests: `test_resolve_runtime_args_triplet_processor_cli_defaults_and_overrides`
Inputs: CLI-mode argument combinations, including processes handling.
Expected outputs/behavior: runtime args resolve defaults/overrides correctly.
Purpose: validate triplet-processor runtime argument resolution behavior.

## Machine Learning tests (`tests/test_ml_random_forest.py`)

- `test_parse_classes_returns_binary_matrix`
  - Inputs: two example 6-bit bitstrings (`"101001"`, `"010010"`).
  - Expected outputs: a (2,6) binary numpy matrix and the original string labels preserved.
  - Purpose: verify `classes` parsing enforces a 6-character 0/1 bitstring and converts to binary targets.

- `test_train_random_forest_smoke`
  - Inputs: `summary_statistics_tsv` fixture (small TSV with `class`, a low-cardinality string feature, and numeric feature columns), runtime config (small forest for speed, `cv_folds=3`, `random_state=7`).
  - Expected outputs: training completes, artifacts exist (`random_forest_model.pkl`, `random_forest_overall_metrics.json`, `random_forest_metrics.txt`, `random_forest_confusion_matrices.png`, `random_forest_confusion_matrix_64_classes.png`, `predictions.tsv`), metrics contain `primary_metrics`, `diagnostic_metrics`, `dataset_summary`, and `timings_seconds`, and the prediction TSV includes `matched_label_count`.
  - Purpose: smoke-test end-to-end training flow, evaluation, and artifact writing.

- `test_train_random_forest_creates_bitwise_metrics_report`
  - Inputs: same TSV fixture, compact training config (`n_estimators=15`, `cv_folds=2`, `random_state=11`).
  - Expected outputs: human-readable metrics file includes per-bit metrics, diagnostic statements, dataset summary, and timings; primary metrics (hamming loss, micro/macro/weighted f1) are present and finite.
  - Purpose: ensure textual and JSON metric artifacts include per-bit breakdowns, consolidated dataset metadata, and primary/diagnostic distinctions.

## Machine Learning config tests (`tests/test_ml_config.py`)

- `test_load_ml_config_defaults_target_column_to_class`
  - Inputs: minimal JSON config with only `input_path` and `output_dir`.
  - Expected outputs: normalized ML config defaults `target_column` to `class`.
  - Purpose: verify the default ML target column matches the dataframe shape used for random forest training.

- `test_load_ml_config_accepts_explicit_class_target_column`
  - Inputs: JSON config with `input_path`, `output_dir`, explicit `target_column: class`, and a simple `model.n_estimators` override.
  - Expected outputs: normalized ML config preserves `target_column: class` and carries through model hyperparameters.
  - Purpose: verify explicit target column handling and basic nested config parsing.

## Machine Learning tests (`tests/test_ml_multi_knn.py`)

- `test_multi_knn_train_smoke`
  - Inputs: `summary_statistics_tsv` fixture with a low-cardinality string column and numeric features, KNN runtime config (`n_neighbors=5`, `cv_folds=3`, `random_state=7`, `weights=uniform`).
  - Expected outputs: training completes, artifacts exist (`multi_knn_model.pkl`, `multi_knn_overall_metrics.json`, `multi_knn_metrics.txt`, `multi_knn_confusion_matrices.png`, `multi_knn_confusion_matrix_64_classes.png`, `predictions.tsv`), metrics include `classifier: multi_knn`, the `knn` details block is present, and the prediction TSV includes `matched_label_count`.
  - Purpose: smoke-test the multi-label KNN baseline end to end.

- `test_multi_knn_build_model_caps_neighbors_to_training_size`
  - Inputs: direct model build request with `n_neighbors=20` and `train_size=2`.
  - Expected outputs: effective neighbor count is capped to 2.
  - Purpose: verify the adaptive neighbor sizing used to avoid KNN failures on small training folds.

- `test_multi_knn_metrics_report_mentions_effective_neighbors`
  - Inputs: same TSV fixture with a larger requested neighbor count (`n_neighbors=20`) and `weights=distance`.
  - Expected outputs: text metrics report includes configured/effective neighbor details and timings, the JSON metrics include the `knn` block and dataset summary, and `feature_importances.tsv` is written.
  - Purpose: validate the KNN-specific reporting, permutation-importance artifact, and timing summary.

## Machine Learning tests (`tests/test_ml_hyper_tune.py`)

- `test_load_hyper_tune_config_accepts_hyperparameter_tuning_section`
  - Inputs: JSON config with a top-level `hyperparameter_tuning` section containing a random-forest grid search space.
  - Expected outputs: tuner config normalizes the model, method, objective, and search space correctly.
  - Purpose: verify the new tuning config header and nested search settings.

- `test_load_hyper_tune_config_fills_model_defaults`
  - Inputs: JSON config with a `hyperparameter_tuning` section for `random_forest` that omits optional model parameters from `search_space`.
  - Expected outputs: normalized tuner config fills in the selected model's default parameters such as `class_weight`, `max_features`, and `min_samples_split`.
  - Purpose: verify that omitted model parameters fall back to trainer defaults during tuning.

- `test_load_hyper_tune_config_rejects_evaluation_section`
  - Inputs: JSON config that includes a top-level `evaluation` section alongside `hyperparameter_tuning`.
  - Expected outputs: `ConfigError` rejecting the unsupported `evaluation` section.
  - Purpose: ensure the tuner config stays isolated from trainer-only settings.

- `test_load_hyper_tune_config_requires_hyperparameter_tuning_section`
  - Inputs: JSON config missing the `hyperparameter_tuning` section.
  - Expected outputs: `ConfigError` requiring the new section.
  - Purpose: ensure the tuner config is explicit and self-contained.

- `test_tune_hyperparameters_grid_search_smoke`
  - Inputs: larger shared ML tuning fixture plus a small random-forest grid search over two candidate values.
  - Expected outputs: tuning writes best-model and results artifacts and returns two ranked candidates.
  - Purpose: exercise the full grid-search tuning path.

- `test_tune_hyperparameters_random_search_smoke`
  - Inputs: larger shared ML tuning fixture plus a small random-forest random search space with `n_iter=1`.
  - Expected outputs: tuning samples one candidate, writes best-model and results artifacts, and returns the best candidate rank.
  - Purpose: exercise the random-search tuning path.

## Backend Parity Tests (`@pytest.mark.backend_parity`)

These tests can be run as a dedicated slice with:

```bash
pytest -m backend_parity
```

Currently marked tests:

- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_z_test_matches_statsmodels_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_standard_z_test_matches_statsmodels_reference_randomized`
- `test_adjust_p_values_custom_matches_standard_randomized`

Expected behavior for this slice:

- Numeric agreement between custom and reference/standard implementations.
- Stable tolerance-bounded parity across randomized samples.
