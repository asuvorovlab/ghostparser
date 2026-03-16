# Test Suite Documentation

This suite provides canonical (non-redundant) coverage grouped by functional area.

## Running Tests

### Run all tests

```bash
pytest
```

### Run one test file

```bash
pytest tests/test_tree_parser.py
```

### Run one test

```bash
pytest tests/test_tree_parser.py::test_extract_triplet_subtree_all_taxa_present
```

### Run reference tests (external parity checks)

```bash
pytest -m reference
```

Run non-reference tests:

```bash
pytest -m "not reference"
```

Reference-marked tests in `tests/test_triplet_processor.py`:
- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_z_test_matches_statsmodels_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_standard_z_test_matches_statsmodels_reference_randomized`
- `test_two_sample_ks_test_hybrid_uses_scipy_near_threshold`
- `test_adjust_p_values_custom_matches_standard_bonferroni_randomized`
- `test_adjust_p_values_custom_matches_standard_fdr_randomized`

## Pipeline-First Coverage Map

### Orchestrator (Primary Pipeline)

- `tests/test_orchestrator.py` covers process resolution and unified CLI/config runtime argument behavior for the end-to-end pipeline entry point.
- End-to-end orchestrator behavior is also exercised via integrated tree parsing + triplet inference flows covered in module tests.

### Tree Parser (Submodule)

- `tests/test_tree_parser.py` covers tree normalization, support filtering, outgroup handling, triplet extraction, and streaming/multiprocessing writers.

### Triplet Processor (Submodule)

- `tests/test_triplet_processor.py` covers DCT/KS/statistics logic, classification outputs, TSV/JSON writing, backend parity checks, and runtime argument resolution.
- Includes correction-method coverage for `none`, `bonferroni`, and `fdr_bh`, with randomized parity checks between custom and standard correction implementations.

### Config Loading (Cross-Cutting)

- `tests/test_config.py` covers config parsing, required fields, defaults, and validation for orchestrator and both submodules.

## Shared Fixtures (`tests/fixtures.py`)

- `simple_newick_file`
  - Input fixture content:
    - `(TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,OutGroup:0.558);`
  - Used in:
    - `test_read_tree_file_single_tree`
    - `test_calculate_average_support_no_values`
    - `test_standardize_tree_preserves_branch_lengths`
    - `test_format_newick_with_precision_trailing_zeros`
    - `test_format_newick_with_precision_default_places`
    - `test_format_newick_with_custom_precision`
    - `test_write_clean_trees`
    - `test_clean_and_save_trees_no_filters`
    - `test_clean_and_save_trees_creates_output_file`
    - `test_get_taxa_from_tree_correct_names`
    - `test_integration_full_workflow`
    - `test_integration_triplets_workflow`
  - Expectations:
    - parsing returns one valid tree
    - average support is `None` when no support labels are present
    - branch lengths are preserved through standardization
    - formatted Newick remains valid and precision behavior is respected
    - cleaned output files are created and non-empty
    - extracted taxa match the expected sorted taxa set
    - integration workflows produce expected counts and output rows

- `newick_with_support_file`
  - Input fixture content:
    - `(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup)0.85;`
  - Used in:
    - `test_calculate_average_support_with_values`
    - `test_remove_support_values`
    - `test_standardize_tree_removes_support`
  - Expectations:
    - average support is computed in the expected range
    - support labels are removed by cleaning/standardization
    - tree remains parseable after support removal

- `multiple_trees_file`
  - Input fixture content:
    - 3 Newick trees (one per line) with shared taxa set and varying topology
  - Used in:
    - `test_read_tree_file_multiple_trees`
    - `test_write_clean_trees_multiple`
  - Expectations:
    - parser returns all trees in the file
    - cleaned write/read roundtrip preserves tree count

- `low_support_tree_file`
  - Input fixture content:
    - one high-support tree and one low-support tree
  - Used in:
    - `test_clean_and_save_trees_filters_low_support`
  - Expectations:
    - low-support tree is dropped
    - dropped index metadata is reported correctly

- `triplet_comparison_cases`
  - Input fixture content:
    - list of `(newick_str, triplet)` pairs for cross-library triplet checks
  - Used in:
    - `test_triplet_branch_lengths_match`
    - `test_triplet_collapse_consistency_dendropy_vs_biopython`
  - Expectations:
    - pairwise patristic distances match between DendroPy and BioPython-derived triplet subtrees
    - DendroPy triplet extraction remains numerically consistent with BioPython prune/collapse behavior

## Triplet Utilities (`tests/test_tree_parser.py`)

- `test_generate_triplets_multiple_outgroups`
  - Confirms list-form outgroups are excluded.
- `test_generate_triplets_outgroup_comma_separated_with_spaces`
  - Confirms comma-separated outgroups with whitespace are excluded.
- `test_read_triplet_filter_file_parses_valid_and_skips_invalid`
  - Confirms valid triplets are parsed and malformed lines are reported.
- `test_filter_triplets_by_taxa_skips_missing_taxa`
  - Confirms triplets with taxa outside available set are skipped.
- `test_write_triplets_to_file`
  - Confirms ordered comma-separated output formatting.
- `test_write_triplet_gene_trees_streaming`
  - Confirms streaming writer emits valid triplet sections.
- `test_write_triplet_gene_trees_multiprocess_triplets_single_worker`
  - Confirms triplet-parallel path with one worker works.
- `test_write_triplet_gene_trees_multiprocess_accepts_list`
  - Confirms list input support and chunk cleanup behavior.
- `test_format_newick_with_precision_triplet_parser`
  - Confirms formatted Newick terminates with `;`.

## Triplet Equivalence (`tests/test_tree_parser.py`)

- `test_triplet_branch_lengths_match` (parametrized)
  - Confirms pairwise patristic distances match between DendroPy-extracted subtree and BioPython distance evaluation.
- `test_triplet_collapse_consistency_dendropy_vs_biopython`
  - Confirms DendroPy triplet extraction is numerically consistent with an independent BioPython prune/collapse path for all three pairwise distances in each triplet.

## Triplet Processing Pipeline (`tests/test_triplet_processor.py`)

- `test_compute_tree_height_statistic_matches_definition`
- `test_classify_triplet_topology_string_for_all_three_topologies`
- `test_pearson_discordant_chi_square_balanced_counts_not_significant`
- `test_two_proportion_discordant_z_test_balanced_counts_not_significant`
- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_z_test_matches_statsmodels_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_run_triplet_pipeline_no_introgression_when_dct_not_significant`
- `test_run_triplet_pipeline_inflow_when_ks_not_significant`
- `test_run_triplet_pipeline_outflow_when_con_summary_higher`
- `test_run_triplet_pipeline_ghost_when_dis_summary_higher`
- `test_parse_analyze_and_write_pipeline_roundtrip_with_species_header`

These tests cover topology classification, discordant-count statistics, KS behavior, and final introgression classification outputs.
They also validate preservation of original p-values alongside corrected p-values in outputs, and dynamic corrected-column naming by correction method.

## Runtime Argument Resolution and Process Defaults

### Orchestrator (`tests/test_orchestrator.py`)

- `test_resolve_runtime_args_cli_custom_processes_preserved`
- `test_resolve_runtime_args_config_processes_preserved_when_set`
- `test_resolve_runtime_args_config_with_cli_warns_and_ignores`

### Tree Parser (`tests/test_tree_parser.py`)

- `test_resolve_runtime_args_tree_parser_cli_custom_processes_preserved`
- `test_resolve_runtime_args_tree_parser_config_warns_and_ignores`
- `test_resolve_runtime_args_tree_parser_config_defaults_processes_to_zero`

### Triplet Processor (`tests/test_triplet_processor.py`)

- `test_resolve_runtime_args_triplet_processor_cli_custom_processes_preserved`
- `test_resolve_runtime_args_triplet_processor_config_processes_preserved_when_set`
- `test_resolve_runtime_args_triplet_processor_config_defaults_processes_to_zero`

### Config loaders (`tests/test_config.py`)

- `test_load_orchestrator_config_defaults_processes_to_zero`
- `test_load_tree_parser_config_defaults_processes_to_zero`
- `test_load_triplet_processor_config_defaults_processes_to_zero`
- `test_path_resolution_absolute_paths`
- `test_path_resolution_relative_paths`
- `test_path_resolution_home_directory`

These tests confirm omitted `processes` defaults to `0`, explicit values are preserved, and config mode precedence is enforced.
They also validate centralized default behavior resolved through normalization (including `discordant_test=chi-square`, `summary_statistic=median`, `stats_backend=standard`, and `p_value_correction=bonferroni`).
Additionally, they validate OS-style path resolution in config values: absolute paths (`/`), relative paths (from current working directory), and home paths (`~`).

### Path-resolution behavior in runtime-arg tests

- Runtime argument resolution tests in `test_orchestrator.py`, `test_tree_parser.py`, and `test_triplet_processor.py` now assert resolved absolute paths for path fields.
- This aligns CLI/config normalization with OS semantics and centralized path handling in `ghostparser.config`.

## Complete Test Function Index

### `tests/test_config.py`

- `test_load_orchestrator_config_json`
- `test_load_orchestrator_config_yaml`
- `test_load_orchestrator_config_missing_required`
- `test_load_orchestrator_config_invalid_discordant_test`
- `test_load_orchestrator_config_invalid_summary_statistic`
- `test_load_orchestrator_config_invalid_stats_backend`
- `test_load_orchestrator_config_defaults_processes_to_zero`
- `test_load_tree_parser_config_json`
- `test_load_tree_parser_config_invalid_no_multiprocessing`
- `test_load_tree_parser_config_defaults_processes_to_zero`
- `test_load_triplet_processor_config_json`
- `test_load_triplet_processor_config_invalid_stats_backend`
- `test_load_triplet_processor_config_missing_input`
- `test_load_triplet_processor_config_defaults_processes_to_zero`
- `test_path_resolution_absolute_paths`
- `test_path_resolution_relative_paths`
- `test_path_resolution_home_directory`

### `tests/test_orchestrator.py`

- `test_resolve_processes_zero_uses_all_cores`
- `test_resolve_parallel_mode`
- `test_resolve_runtime_args_cli_defaults`
- `test_resolve_runtime_args_cli_custom_processes_preserved`
- `test_resolve_runtime_args_config_with_cli_warns_and_ignores`
- `test_resolve_runtime_args_config_processes_preserved_when_set`

### `tests/test_tree_parser.py`

- `test_resolve_runtime_args_tree_parser_cli_defaults`
- `test_resolve_runtime_args_tree_parser_cli_custom_processes_preserved`
- `test_resolve_runtime_args_tree_parser_config_warns_and_ignores`
- `test_resolve_runtime_args_tree_parser_config_defaults_processes_to_zero`
- `test_read_tree_file_single_tree`
- `test_read_tree_file_multiple_trees`
- `test_read_tree_file_not_found`
- `test_read_tree_file_invalid_newick`
- `test_read_tree_file_random_text`
- `test_read_tree_file_empty_file`
- `test_calculate_average_support_with_values`
- `test_calculate_average_support_no_values`
- `test_remove_support_values`
- `test_standardize_tree_removes_support`
- `test_standardize_tree_preserves_branch_lengths`
- `test_format_newick_with_precision_trailing_zeros`
- `test_format_newick_with_precision_default_places`
- `test_format_newick_with_custom_precision`
- `test_write_clean_trees`
- `test_write_clean_trees_multiple`
- `test_clean_and_save_trees_filters_low_support`
- `test_clean_and_save_trees_no_filters`
- `test_clean_and_save_trees_creates_output_file`
- `test_clean_and_save_gene_trees_discards_missing_outgroup`
- `test_get_taxa_from_tree_correct_names`
- `test_generate_triplets_count`
- `test_generate_triplets_excludes_outgroup`
- `test_generate_triplets_content`
- `test_generate_triplets_large_set`
- `test_write_triplets_to_file`
- `test_write_triplets_to_file_empty`
- `test_get_clean_filename_simple`
- `test_get_clean_filename_different_extension`
- `test_get_clean_filename_no_extension`
- `test_integration_full_workflow`
- `test_integration_triplets_workflow`
- `test_extract_triplet_subtree_all_taxa_present`
- `test_extract_triplet_subtree_missing_taxa`
- `test_extract_triplet_subtree_preserves_branch_lengths`
- `test_process_gene_trees_for_triplets`
- `test_process_gene_trees_for_triplets_empty`
- `test_write_triplet_gene_trees`
- `test_write_triplet_gene_trees_includes_species_tree_header`
- `test_write_triplet_gene_trees_empty_triplet`
- `test_integration_full_triplet_extraction_workflow`
- `test_write_triplet_gene_trees_multiprocess_with_workers`
- `test_write_triplet_gene_trees_multiprocess_includes_species_header`
- `test_metrics_logger_context_manager`
- `test_metrics_logger_file_not_opened_before_enter`
- `test_triplet_gene_trees_separator_format`
- `test_multiprocessing_triplet_writer_handles_empty_triplets`
- `test_generate_triplets_multiple_outgroups`
- `test_generate_triplets_outgroup_comma_separated_with_spaces`
- `test_read_triplet_filter_file_parses_valid_and_skips_invalid`
- `test_filter_triplets_by_taxa_skips_missing_taxa`
- `test_write_triplet_gene_trees_streaming`
- `test_write_triplet_gene_trees_multiprocess_triplets_single_worker`
- `test_write_triplet_gene_trees_multiprocess_accepts_list`
- `test_build_species_triplet_metadata_normalizes_abc`
- `test_format_newick_with_precision_triplet_parser`
- `test_triplet_branch_lengths_match`
- `test_triplet_collapse_consistency_dendropy_vs_biopython`

### `tests/test_triplet_processor.py`

- `test_compute_tree_height_statistic_matches_definition`
- `test_classify_triplet_topology_string_for_all_three_topologies`
- `test_classify_triplet_topology_labels_concordant_and_discordants`
- `test_pearson_discordant_chi_square_balanced_counts_not_significant`
- `test_two_proportion_discordant_z_test_balanced_counts_not_significant`
- `test_custom_chi_square_matches_scipy_reference_randomized`
- `test_custom_ks_matches_scipy_asymptotic_reference_randomized`
- `test_run_triplet_pipeline_uses_species_concordant_and_frequency_ranked_discordants`
- `test_run_triplet_pipeline_supports_z_test_for_discordant_counts`
- `test_run_triplet_pipeline_supports_standard_stats_backend`
- `test_standard_z_test_matches_statsmodels_reference_randomized`
- `test_run_triplet_pipeline_supports_median_summary_statistic`
- `test_run_triplet_pipeline_breaks_discordant_ties_by_first_topology`
- `test_run_triplet_pipeline_no_introgression_when_dct_not_significant`
- `test_run_triplet_pipeline_inflow_when_ks_not_significant`
- `test_run_triplet_pipeline_outflow_when_con_summary_higher`
- `test_run_triplet_pipeline_ghost_when_dis_summary_higher`
- `test_parse_analyze_and_write_pipeline_roundtrip_with_species_header`
- `test_analyze_triplet_gene_tree_file_with_multiprocessing`
- `test_analyze_triplet_gene_tree_file_rejects_unknown_discordant_test`
- `test_analyze_triplet_gene_tree_file_rejects_unknown_summary_statistic`
- `test_run_triplet_pipeline_supports_mode_summary_statistic`
- `test_analyze_triplet_gene_tree_file_rejects_unknown_stats_backend`
- `test_two_sample_ks_test_hybrid_uses_scipy_near_threshold`
- `test_two_sample_ks_test_hybrid_keeps_custom_when_not_borderline`
- `test_two_sample_ks_test_hybrid_rejects_negative_margin`
- `test_parse_triplet_gene_trees_file_requires_species_tree_column`
- `test_parse_triplet_gene_trees_file_rejects_empty_species_tree`
- `test_write_pipeline_results_omits_removed_topology_columns`
- `test_write_pipeline_results_uses_dynamic_summary_column_names`
- `test_collect_triplet_statistics_returns_dict_list`
- `test_write_pipeline_results_uses_dct_chi_statistic_column_for_chi_square`
- `test_write_pipeline_results_uses_dct_z_score_column_for_z_test`
- `test_write_pipeline_results_uses_mode_summary_columns_for_mode`
- `test_write_pipeline_statistics_json`
- `test_write_pipeline_results_rejects_mixed_discordant_test_outputs`
- `test_write_pipeline_results_rejects_unsupported_summary_statistic`
- `test_resolve_runtime_args_triplet_processor_cli_defaults`
- `test_resolve_runtime_args_triplet_processor_cli_custom_processes_preserved`
- `test_resolve_runtime_args_triplet_processor_config_warns_and_ignores`
- `test_resolve_runtime_args_triplet_processor_config_processes_preserved_when_set`
- `test_resolve_runtime_args_triplet_processor_config_defaults_processes_to_zero`
