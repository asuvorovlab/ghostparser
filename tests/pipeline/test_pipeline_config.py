"""Config resolution and config-file mode for the streaming pipeline CLI."""

import argparse
import json

import pytest

from ghostparser.pipeline.config import (
    ConfigError,
    build_argument_parser,
    load_pipeline_config,
    resolve_config,
)


def _base_cli_args(**overrides):
    """Build a parsed-CLI namespace with all pipeline args defaulted to None.

    Args:
        **overrides: Attribute values to set on the namespace.

    Returns:
        An ``argparse.Namespace`` mirroring what the pipeline parser produces.
    """
    args = argparse.Namespace(
        config_file=None,
        species_tree_path="species.tree",
        gene_trees_path="genes.tree",
        outgroups="OUT",
        output_folder=None,
        triplet_filter=None,
        no_overwrite=None,
        processes=None,
        parallelization_mode=None,
        alpha_dct=None,
        alpha_ks=None,
        p_value_correction=None,
        summary_statistic=None,
        consolidation=None,
        bootstrap=None,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_cli_defaults_resolve():
    """CLI-only mode fills config+CLI and config-only keys with shared defaults."""
    config = resolve_config(_base_cli_args())
    assert config["alpha_dct"] == 0.05
    assert config["alpha_ks"] == 0.05
    # Pipeline-specific defaults (differ from the shared orchestrator defaults).
    assert config["p_value_correction"] == "bfn"
    assert config["summary_statistic"] == "mean"
    assert config["overwrite"] is True
    # Config-file-only keys take their defaults in CLI mode.
    assert config["discordant_test"] == "chi-square"
    assert config["tree_height_calculation_strategy"] == "AVG"
    assert config["min_support_value"] == 0.5
    assert config["bootstrap_iterations"] == 100
    assert config["bootstrap_seed"] is None
    assert config["generate_summary_stats"] is False
    assert config["bootstrap_debug_mode"] is False
    assert config["bootstrap_summary_only"] is False
    # The custom stats backend is gone: no stats_backend key at all.
    assert "stats_backend" not in config


def test_cli_overrides_for_config_plus_cli_options():
    """The config+CLI flags override their defaults when supplied on the CLI."""
    args = _base_cli_args(
        alpha_dct=0.01,
        alpha_ks=0.2,
        summary_statistic="mean",
        p_value_correction="fdr_bh",
        no_overwrite=True,
    )
    config = resolve_config(args)
    assert config["alpha_dct"] == 0.01
    assert config["alpha_ks"] == 0.2
    assert config["summary_statistic"] == "mean"
    assert config["p_value_correction"] == "fdr_bh"
    assert config["overwrite"] is False


def test_config_only_keys_read_from_config_file(tmp_path):
    """Config-file-only keys (and nested bootstrap_options) load from a file."""
    payload = {
        "species_tree_path": "species.tree",
        "gene_trees_path": "genes.tree",
        "outgroups": "OUT",
        "discordant_test": "z-test",
        "tree_height_calculation_strategy": "SIS",
        "min_support_value": 0.9,
        "generate_summary_stats": True,
        "alpha_dct": 0.02,
        "bootstrap_options": {
            "iterations": 25,
            "seed": 7,
            "debug_mode": True,
            "summary_only": True,
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    config = load_pipeline_config(str(config_path))
    assert config["discordant_test"] == "z-test"
    assert config["tree_height_calculation_strategy"] == "SIS"
    assert config["min_support_value"] == 0.9
    assert config["generate_summary_stats"] is True
    assert config["alpha_dct"] == 0.02
    assert config["bootstrap_iterations"] == 25
    assert config["bootstrap_seed"] == 7
    assert config["bootstrap_debug_mode"] is True
    assert config["bootstrap_summary_only"] is True


def test_config_file_wins_over_cli(tmp_path, capsys):
    """In config-file mode the file wins and ignored CLI flags trigger a warning."""
    payload = {
        "species_tree_path": "file_species.tree",
        "gene_trees_path": "file_genes.tree",
        "outgroups": "OUT",
        "alpha_dct": 0.03,
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    args = _base_cli_args(
        config_file=str(config_path), alpha_dct=0.5, summary_statistic="mode"
    )
    config = resolve_config(args)

    assert config["alpha_dct"] == 0.03
    # The file omits summary_statistic, so it takes the pipeline default (not the
    # CLI value "mode"), proving the CLI flag was ignored.
    assert config["summary_statistic"] == "mean"
    assert config["species_tree"].endswith("file_species.tree")
    assert "--config-file provided" in capsys.readouterr().out


def test_missing_required_field_raises():
    """A missing required input raises ConfigError in CLI mode."""
    with pytest.raises(ConfigError):
        resolve_config(_base_cli_args(species_tree_path=None))


def test_parser_exposes_config_file_and_new_flags():
    """The parser exposes -c/--config-file and the new config+CLI flags."""
    parser = build_argument_parser()
    opts = parser.parse_args(
        [
            "-st", "s", "-gt", "g", "-og", "OUT",
            "--alpha-dct", "0.01",
            "--alpha-ks", "0.2",
            "--p-value-correction", "fdr_bh",
            "--summary-statistic", "mean",
            "--no-overwrite",
        ]
    )
    assert opts.config_file is None
    assert opts.alpha_dct == 0.01
    assert opts.alpha_ks == 0.2
    assert opts.p_value_correction == "fdr_bh"
    assert opts.summary_statistic == "mean"
    assert opts.no_overwrite is True
