# Test Suite Documentation

This document lists the current tests, fixtures, and marker-based slices, with function-level input/output expectations.

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
Inputs: full orchestrator JSON config (paths, methods, thresholds).
Expected outputs: normalized absolute paths, parsed outgroups list, preserved explicit values.

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

- `test_load_tree_parser_config_json`
Inputs: full tree parser JSON config including `no_multiprocessing`.
Expected outputs: normalized paths and preserved explicit values.

- `test_load_tree_parser_config_invalid_no_multiprocessing`
Inputs: non-boolean `no_multiprocessing`.
Expected outputs: `ConfigError`.

- `test_load_tree_parser_config_defaults_processes_to_zero`
Inputs: no `processes`.
Expected outputs: `processes == 0`.

- `test_load_triplet_processor_config_json`
Inputs: full triplet-processor JSON config.
Expected outputs: normalized paths and preserved explicit methods/thresholds.

- `test_load_triplet_processor_config_invalid_choice_fields` (parametrized)
Inputs: invalid `p_value_correction`, `stats_backend`.
Expected outputs: `ConfigError` referencing offending field.

- `test_load_triplet_processor_config_tree_height_strategy_validation` (parametrized)
Inputs: supported values (`AVG`, `A`, `B`, `C`, `SIS`, `INT`) and an invalid value (`D`).
Expected outputs: supported values load successfully; invalid value raises `ConfigError`.

- `test_load_triplet_processor_config_missing_input`
Inputs: config without `input_path`.
Expected outputs: `ConfigError` on required input.

- `test_load_triplet_processor_config_defaults_processes_to_zero`
Inputs: minimal config.
Expected outputs: `processes == 0`, `tree_height_calculation_strategy == "AVG"`, `p_value_correction == "no"`.

- `test_load_triplet_processor_config_invalid_bootstrap_fields` (parametrized)
Inputs: invalid bootstrap payload variants (`bootstrap` non-bool, invalid `iterations`, invalid `seed`, invalid `debug_mode`, invalid `summary_only`).
Expected outputs: `ConfigError` references the offending bootstrap field.

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
Expected outputs: default statistical settings/path resolution and preserved CLI overrides.

- `test_resolve_runtime_args_config_with_cli_warns_and_ignores`
Inputs: config file + conflicting CLI args.
Expected outputs: warning emitted; config values/defaults win over CLI extras.

- `test_resolve_runtime_args_config_processes_behavior` (parametrized)
Inputs: config with explicit `processes` and config without `processes`.
Expected outputs: preserves configured value or defaults to `0`.

### `tests/test_tree_parser.py`

Runtime-arg resolution:

- `test_resolve_runtime_args_tree_parser_cli_processes` (parametrized)
Inputs: CLI defaults and explicit `processes=6`.
Expected outputs: path resolution, outgroup parsing, default `min_support_value`, and expected process behavior.

- `test_resolve_runtime_args_tree_parser_config_warns_and_ignores`
Inputs: config mode + extra CLI args.
Expected outputs: warning + config precedence.

- `test_resolve_runtime_args_tree_parser_config_processes_behavior` (parametrized)
Inputs: config with explicit `processes` and config without `processes`.
Expected outputs: preserves configured value or defaults to `0`.

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

- Tests cover concordant/discordant label mapping, tie behavior, relabeling under canonicalization, and the 5 output classes (`no_introgression`, `inflow_introgression`, `outflow_introgression`, `ghost_introgression`, `unresolved` where applicable).
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
Expected outputs: base bootstrap columns (`bootstrap_value`, `all_bootstrap`) are present, debug bootstrap columns are present when debug mode is enabled, and summary-mode payload cells use JSON-style object strings.

- `test_write_pipeline_results_adds_bootstrap_gene_tree_heights_when_summary_only_false`
Inputs: bootstrap-enabled result written with debug mode enabled and `summary_only=false`.
Expected outputs: TSV includes `bootstrap_gene_tree_heights` and stores the raw per-triplet tree-height list.

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

- Covers CLI defaults, config precedence, and `processes` default/preservation behavior.
Inputs: CLI-only args and config-file mode args.
Expected outputs: resolved defaults and correct precedence semantics.

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
