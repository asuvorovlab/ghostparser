"""Config resolution and config-file mode for the orchestrator CLI."""

import argparse
import json
from pathlib import Path

import pytest
import yaml

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
        outgroup="OUT",
        output_folder=None,
        triplet_filter=None,
        no_overwrite=None,
        processes=None,
        alpha_dct=None,
        alpha_ks=None,
        alpha_perm=None,
        p_value_correction=None,
        pipeline_mode=None,
        consolidation=None,
        bootstrap=None,
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
    assert config["permutation_min_resamples"] == 2500
    assert config["permutation_max_resamples"] == 25000
    assert config["permutation_ci_method"] == "wilson"
    assert config["overwrite"] is True
    # Config-file-only keys take their defaults in CLI mode.
    assert config["discordant_test"] == "chi-square"
    assert config["tree_height_calculation_strategy"] == "AVG"
    assert config["min_support_value"] == 0.5
    assert config["bootstrap_iterations"] == 100
    assert config["seed"] is None
    assert config["generate_summary_stats"] is False
    assert config["bootstrap_debug_mode"] is False
    assert config["bootstrap_summary_only"] is False
    assert config["preflight_data_check"] is False
    assert config["pipeline_mode"] == "efficient"


def test_config_only_keys_read_from_config_file(tmp_path):
    """Config-file-only keys (and nested bootstrap_options) load from a file."""
    payload = {
        "species_tree_path": "species.tree",
        "gene_trees_path": "genes.tree",
        "outgroup": "OUT",
        "discordant_test": "z-test",
        "tree_height_calculation_strategy": "SIS",
        "min_support_value": 0.9,
        "generate_summary_stats": True,
        "alpha_dct": 0.02,
        "seed": 7,
        "bootstrap_options": {
            "iterations": 25,
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
    assert config["seed"] == 7
    assert config["bootstrap_debug_mode"] is True
    assert config["bootstrap_summary_only"] is True


@pytest.mark.parametrize(
    "written, expected",
    [
        ("no", "no"),
        ('"no"', "no"),
        ("No", "no"),
        ("NO", "no"),
        ("bfn", "bfn"),
        ("fdr_bh", "fdr_bh"),
    ],
)
def test_p_value_correction_accepts_yaml_bare_word_no(tmp_path, written, expected):
    """`p_value_correction: no` resolves to the "no" choice, quoted or not.

    YAML 1.1 resolves the bare word ``no`` to boolean ``False``, so an unquoted
    value never reaches validation as a string. Every spelling of the choice
    must still select it.
    """
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "species_tree_path: species.tree\n"
        "gene_trees_path: genes.tree\n"
        "outgroup: OUT\n"
        f"p_value_correction: {written}\n"
    )

    config = load_orchestrator_config(str(config_path))
    assert config["p_value_correction"] == expected


def test_p_value_correction_rejects_a_value_with_no_matching_choice(tmp_path):
    """A boolean with no equivalent choice is still an error, and names the value."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "species_tree_path: species.tree\n"
        "gene_trees_path: genes.tree\n"
        "outgroup: OUT\n"
        "p_value_correction: yes\n"
    )

    with pytest.raises(ConfigError, match="must be one of"):
        load_orchestrator_config(str(config_path))


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
        "outgroup": "OUT",
    }
    if payload_value is not None:
        payload["preflight_data_check"] = payload_value
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    config = load_orchestrator_config(str(config_path))
    assert config["preflight_data_check"] is expected


def test_config_file_wins_over_cli(tmp_path):
    """In config-file mode the file wins and every CLI flag is ignored."""
    payload = {
        "species_tree_path": "file_species.tree",
        "gene_trees_path": "file_genes.tree",
        "outgroup": "OUT",
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


@pytest.mark.parametrize(
    "value,expected",
    [
        ("OUT", ["OUT"]),
        ("Out1,Out2", ["Out1", "Out2"]),
        (" Out1 , Out2 ,", ["Out1", "Out2"]),
        (["Out1", "Out2"], ["Out1", "Out2"]),
        (("Out1", "Out2"), ["Out1", "Out2"]),
        (["Out1,Out2", "Out3"], ["Out1", "Out2", "Out3"]),
    ],
)
def test_outgroup_accepts_single_comma_separated_and_list_forms(value, expected):
    """One `outgroup` key covers a single label, a comma-separated string, or a list."""
    config = resolve_config(_base_cli_args(outgroup=value))
    assert config["outgroup"] == expected


@pytest.mark.parametrize("value", [None, "", "  ", ",", [], ["", "  "], 42])
def test_outgroup_rejects_empty_and_non_label_values(value):
    """A missing or label-free `outgroup` is a config error, not an empty list."""
    with pytest.raises(ConfigError, match="outgroup"):
        resolve_config(_base_cli_args(outgroup=value))


@pytest.mark.parametrize(
    "sample_name", ["orchestrator_minimal.yaml", "orchestrator_full.yaml"]
)
def test_shipped_sample_configs_resolve(sample_name):
    """The sample configs load cleanly and use the current key names.

    Guards against the samples drifting out of step with the validator, which
    would leave users copying a config the loader rejects.
    """
    sample_dir = Path(__file__).resolve().parents[2] / "sample_configs"
    config = load_orchestrator_config(str(sample_dir / sample_name))

    assert config["outgroup"]
    assert all(isinstance(label, str) and label for label in config["outgroup"])
    assert config["alpha_perm"] == 0.05
    assert config["permutation_ci_method"] == "wilson"


def test_full_sample_config_covers_every_runtime_key():
    """`orchestrator_full.yaml` documents every key the normalizer produces.

    Resolved-only keys (`species_tree`/`gene_trees`/`output`) are the renamed
    forms of the `*_path`/`output_folder` inputs, and the nested blocks are
    flattened with a prefix, so both are normalized before comparing.
    """
    sample_dir = Path(__file__).resolve().parents[2] / "sample_configs"
    sample_path = sample_dir / "orchestrator_full.yaml"
    raw = yaml.safe_load(sample_path.read_text())

    documented = set(raw)
    for block in ("permutation_options", "bootstrap_options"):
        documented.update(raw.get(block) or {})

    resolved = load_orchestrator_config(str(sample_path))
    renamed = {"species_tree", "gene_trees", "output"}
    undocumented = [
        key
        for key in resolved
        if key not in renamed
        and key not in documented
        and key.removeprefix("permutation_") not in documented
        and key.removeprefix("bootstrap_") not in documented
    ]
    assert undocumented == []


def test_missing_required_field_raises():
    """A missing required input raises ConfigError in CLI mode."""
    with pytest.raises(ConfigError):
        resolve_config(_base_cli_args(species_tree_path=None))


def test_parser_flags_resolve_into_their_config_values():
    """The CLI flag names are wired through the parser to the resolved config.

    Every other config test builds an ``argparse.Namespace`` directly, so this is
    the only place the actual flag strings are pinned; resolving the parsed args
    covers the override path in the same pass.
    """
    parser = build_argument_parser()
    opts = parser.parse_args(
        [
            "-st", "s", "-gt", "g", "-og", "OUT",
            "--alpha-dct", "0.01",
            "--alpha-ks", "0.2",
            "--p-value-correction", "fdr_bh",
            "--pipeline-mode", "detailed",
            "--alpha-perm", "0.02",
            "--no-overwrite",
            "--preflight-data-check",
        ]
    )
    assert opts.config_file is None

    config = resolve_config(opts)
    assert config["alpha_dct"] == 0.01
    assert config["alpha_ks"] == 0.2
    assert config["alpha_perm"] == 0.02
    assert config["pipeline_mode"] == "detailed"
    assert config["p_value_correction"] == "fdr_bh"
    assert config["overwrite"] is False
    assert config["preflight_data_check"] is True


def _payload(tmp_path, **extra):
    """Write a minimal orchestrator config file with extra keys merged in.

    Args:
        tmp_path: The pytest temporary directory.
        **extra: Additional config keys.

    Returns:
        The path to the written YAML file.
    """
    payload = {
        "species_tree_path": "species.tree",
        "gene_trees_path": "genes.tree",
        "outgroup": "OUT",
    }
    payload.update(extra)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(payload))
    return path


def test_pipeline_mode_rejects_an_unknown_value(tmp_path):
    """An unsupported pipeline mode is refused by name."""
    with pytest.raises(ConfigError, match="pipeline_mode"):
        load_orchestrator_config(str(_payload(tmp_path, pipeline_mode="fast")))


def test_seed_rejects_a_non_integer(tmp_path):
    """The run-wide seed must be an integer when provided."""
    with pytest.raises(ConfigError, match="seed"):
        load_orchestrator_config(str(_payload(tmp_path, seed="abc")))

