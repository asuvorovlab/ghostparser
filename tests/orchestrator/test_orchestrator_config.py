"""Config resolution and config-file mode for the orchestrator CLI.

Every test here is marked ``config``: they cover how a config is loaded,
resolved and validated, not what the pipeline computes. Run the pipeline logic
alone with ``pytest -m "not config"``.

Individual default values are deliberately not pinned. CONFIG.md and
``sample_configs/orchestrator_full.yaml`` state them, and two invariants keep
those in step with the code: CLI mode and config-file mode resolve to one and
the same set of defaults, and the full sample names every runtime key at its
default.
"""

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

pytestmark = pytest.mark.config

_SAMPLE_DIR = Path(__file__).resolve().parents[2] / "sample_configs"
_REQUIRED = {
    "species_tree_path": "species.tree",
    "gene_trees_path": "genes.tree",
    "outgroup": "OUT",
}


def _base_cli_args(**overrides):
    """Build a CLI namespace with the required paths set and every option ``None``.

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
        preflight_triplet_cap=None,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _payload(tmp_path, **extra):
    """Write a minimal orchestrator config file with extra keys merged in.

    A key given as ``None`` is written as null, which the loader treats as
    absent -- the way to express a missing required field.

    Args:
        tmp_path: The pytest temporary directory.
        **extra: Additional config keys.

    Returns:
        The path to the written YAML file.
    """
    payload = {**_REQUIRED, **extra}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(payload))
    return path


def test_cli_and_config_file_share_one_set_of_defaults(tmp_path):
    """A bare CLI invocation and a required-keys-only file resolve identically.

    The two entry paths fill config+CLI keys and config-file-only keys from the
    same constants; comparing whole resolved dicts pins that without naming any
    default, so a changed default never needs a test edit.
    """
    opts = build_argument_parser().parse_args(
        ["-st", "species.tree", "-gt", "genes.tree", "-og", "OUT"]
    )

    assert resolve_config(opts) == load_orchestrator_config(str(_payload(tmp_path)))


def test_config_only_keys_and_nested_blocks_flatten_from_a_file(tmp_path):
    """Config-file-only keys load, and ``bootstrap_options`` flattens to a prefix."""
    expected = {
        "discordant_test": "z-test",
        "tree_height_calculation_strategy": "SIS",
        "min_support_value": 0.9,
        "generate_summary_stats": True,
        "alpha_dct": 0.02,
        "seed": 7,
        "bootstrap_iterations": 25,
        "bootstrap_debug_mode": True,
        "bootstrap_summary_only": True,
    }
    config_path = _payload(
        tmp_path,
        discordant_test="z-test",
        tree_height_calculation_strategy="SIS",
        min_support_value=0.9,
        generate_summary_stats=True,
        alpha_dct=0.02,
        seed=7,
        bootstrap_options={"iterations": 25, "debug_mode": True, "summary_only": True},
    )

    config = load_orchestrator_config(str(config_path))

    assert {key: config[key] for key in expected} == expected


@pytest.mark.parametrize(
    "written, expected",
    [("no", "no"), ('"no"', "no"), ("No", "no"), ("NO", "no"), ("bfn", "bfn")],
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


def test_config_file_wins_over_cli(tmp_path):
    """In config-file mode the file wins and every CLI flag is ignored."""
    payload = {**_REQUIRED, "species_tree_path": "file_species.tree", "alpha_dct": 0.03}
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    args = _base_cli_args(
        config_file=str(config_path), alpha_dct=0.5, alpha_perm=0.5
    )
    config = resolve_config(args)

    assert config["alpha_dct"] == 0.03
    # The file omits alpha_perm, so it takes the default rather than the CLI's
    # 0.5, proving the CLI flag was ignored rather than merged.
    assert config["alpha_perm"] != 0.5
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
    "extra, match",
    [
        ({"species_tree_path": None}, "species_tree_path"),
        ({"p_value_correction": True}, "must be one of"),
        ({"pipeline_mode": "fast"}, "pipeline_mode"),
        ({"seed": "abc"}, "seed"),
        ({"preflight_triplet_cap": -1}, "preflight_triplet_cap"),
        ({"species_rename_map": "absent.tsv"}, "species_rename_map"),
    ],
)
def test_invalid_values_are_rejected_by_field_name(tmp_path, extra, match):
    """A missing required key or an out-of-domain value fails naming the field.

    One case per validator shape: a required path, a choice list (fed a
    boolean with no matching choice, which is what YAML makes of a bare
    ``yes``), an optional int, a non-negative int whose ``0`` already means
    "no cap", and a path whose file is read when the config resolves.
    """
    with pytest.raises(ConfigError, match=match):
        load_orchestrator_config(str(_payload(tmp_path, **extra)))


@pytest.mark.parametrize(
    "sample_name", ["orchestrator_minimal.yaml", "orchestrator_full.yaml"]
)
def test_shipped_sample_configs_resolve(sample_name):
    """The sample configs load cleanly and use the current key names.

    Guards against the samples drifting out of step with the validator, which
    would leave users copying a config the loader rejects.
    """
    config = load_orchestrator_config(str(_SAMPLE_DIR / sample_name))

    assert config["outgroup"]
    assert all(isinstance(label, str) and label for label in config["outgroup"])


def test_full_sample_config_names_every_runtime_key_at_its_default(tmp_path):
    """`orchestrator_full.yaml` lists every key the normalizer produces, at its default.

    Two checks. The raw file must name every resolved key (resolved-only keys
    `species_tree`/`gene_trees`/`output` are the renamed `*_path`/`output_folder`
    inputs, and the nested blocks flatten with a prefix, so both are normalized
    before comparing). And loading it must equal loading a required-keys-only
    file, so the value the sample shows for each key is the value the code
    would have used anyway.
    """
    sample_path = _SAMPLE_DIR / "orchestrator_full.yaml"
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

    defaults = load_orchestrator_config(str(_payload(tmp_path)))
    inputs = renamed | {"outgroup"}
    assert {k: v for k, v in resolved.items() if k not in inputs} == {
        k: v for k, v in defaults.items() if k not in inputs
    }


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
            "--preflight-triplet-cap", "0",
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
    assert config["preflight_triplet_cap"] == 0
