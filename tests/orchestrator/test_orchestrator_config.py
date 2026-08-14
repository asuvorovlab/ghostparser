"""Config resolution and config-file mode for the orchestrator CLI."""

import argparse
import json

import pytest

from ghostparser.orchestrator.config import (
    ConfigError,
    build_argument_parser,
    load_orchestrator_config,
    resolve_config,
)


def _base_cli_args(**overrides):
    """Build a parsed-CLI namespace with all orchestrator args defaulted to None.

    Args:
        **overrides: Attribute values to set on the namespace.

    Returns:
        An ``argparse.Namespace`` mirroring what the orchestrator parser produces.
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
        alpha_perm=None,
        p_value_correction=None,
        consolidation=None,
        bootstrap=None,
        permutation_test=None,
        preflight_data_check=None,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_cli_defaults_resolve():
    """CLI-only mode fills config+CLI and config-only keys with shared defaults."""
    config = resolve_config(_base_cli_args())
    assert config["alpha_dct"] == 0.05
    assert config["alpha_ks"] == 0.05
    # Defaults specific to this package.
    assert config["p_value_correction"] == "bfn"
    assert config["alpha_perm"] == 0.05
    assert config["permutation_test"] is True
    assert config["permutation_min_resamples"] == 2500
    assert config["permutation_max_resamples"] == 25000
    assert config["permutation_ci_method"] == "wilson"
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
    assert config["preflight_data_check"] is False
    # The custom stats backend is gone: no stats_backend key at all.
    assert "stats_backend" not in config


def test_cli_overrides_for_config_plus_cli_options():
    """The config+CLI flags override their defaults when supplied on the CLI."""
    args = _base_cli_args(
        alpha_dct=0.01,
        alpha_ks=0.2,
        alpha_perm=0.02,
        p_value_correction="fdr_bh",
        no_overwrite=True,
    )
    config = resolve_config(args)
    assert config["alpha_dct"] == 0.01
    assert config["alpha_ks"] == 0.2
    assert config["alpha_perm"] == 0.02
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

    config = load_orchestrator_config(str(config_path))
    assert config["discordant_test"] == "z-test"
    assert config["tree_height_calculation_strategy"] == "SIS"
    assert config["min_support_value"] == 0.9
    assert config["generate_summary_stats"] is True
    assert config["alpha_dct"] == 0.02
    assert config["bootstrap_iterations"] == 25
    assert config["bootstrap_seed"] == 7
    assert config["bootstrap_debug_mode"] is True
    assert config["bootstrap_summary_only"] is True


@pytest.mark.parametrize(
    "payload_value, expected",
    [(True, True), (False, False), (None, False)],
)
def test_preflight_data_check_resolves_from_config_file(
    tmp_path, payload_value, expected
):
    """preflight_data_check reads from a config file and defaults to False."""
    payload = {
        "species_tree_path": "species.tree",
        "gene_trees_path": "genes.tree",
        "outgroups": "OUT",
    }
    if payload_value is not None:
        payload["preflight_data_check"] = payload_value
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    config = load_orchestrator_config(str(config_path))
    assert config["preflight_data_check"] is expected


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
        config_file=str(config_path), alpha_dct=0.5, alpha_perm=0.5
    )
    config = resolve_config(args)

    assert config["alpha_dct"] == 0.03
    # The file omits alpha_perm, so it takes the orchestrator default (not the
    # CLI value 0.5), proving the CLI flag was ignored.
    assert config["alpha_perm"] == 0.05
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
            "--alpha-perm", "0.02",
            "--no-permutation-test",
            "--no-overwrite",
            "--preflight-data-check",
        ]
    )
    assert opts.config_file is None
    assert opts.preflight_data_check is True
    assert opts.alpha_dct == 0.01
    assert opts.alpha_ks == 0.2
    assert opts.p_value_correction == "fdr_bh"
    assert opts.alpha_perm == 0.02
    assert opts.permutation_test is False
    assert opts.no_overwrite is True
