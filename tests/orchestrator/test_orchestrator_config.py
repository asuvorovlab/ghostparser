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


def _payload(tmp_path, **extra):
    """Write a minimal orchestrator config file with extra keys merged in.

    A key given as ``None`` is written as null, which the loader treats as
    absent, the way to express a missing required field.

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


def test_cli_and_config_file_share_one_set_of_defaults():
    """A bare CLI invocation and the shipped minimal sample resolve identically.

    The two entry paths fill config+CLI keys and config-file-only keys from the
    same constants; comparing whole resolved dicts pins that without naming any
    default, so a changed default never needs a test edit. Loading the minimal
    sample also proves it stays loadable.
    """
    sample = yaml.safe_load((_SAMPLE_DIR / "orchestrator_minimal.yaml").read_text())
    opts = build_argument_parser().parse_args(
        [
            "-st", sample["species_tree_path"],
            "-gt", sample["gene_trees_path"],
            "-og", sample["outgroup"],
        ]
    )

    assert resolve_config(opts) == load_orchestrator_config(
        str(_SAMPLE_DIR / "orchestrator_minimal.yaml")
    )


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
        "bootstrap_diagnostic": True,
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
        bootstrap_options={"iterations": 25, "diagnostic": True, "summary_only": True},
    )

    config = load_orchestrator_config(str(config_path))

    assert {key: config[key] for key in expected} == expected


@pytest.mark.parametrize("written", ["no", '"no"', "NO"])
def test_p_value_correction_accepts_yaml_bare_word_no(tmp_path, written):
    """`p_value_correction: no` resolves to the "no" choice, quoted or not.

    YAML 1.1 resolves the bare words ``no`` and ``NO`` to boolean ``False``,
    so an unquoted value never reaches validation as a string. Every spelling
    of the choice must still select it.
    """
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "species_tree_path: species.tree\n"
        "gene_trees_path: genes.tree\n"
        "outgroup: OUT\n"
        f"p_value_correction: {written}\n"
    )

    config = load_orchestrator_config(str(config_path))
    assert config["p_value_correction"] == "no"


def test_cli_flags_override_the_config_file(tmp_path):
    """Each flag string reaches its key and, beside a config file, replaces the file's value.

    Keys the command line leaves alone keep the file's value, and the negated
    switches (`--no-overwrite`, `--no-consolidation`) beat the file's
    canonical `true` rather than losing to it inside the normalizer. This is
    the one place the flag strings themselves are pinned.
    """
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                **_REQUIRED,
                "species_tree_path": "file_species.tree",
                "alpha_dct": 0.03,
                "alpha_ks": 0.04,
                "p_value_correction": "no",
                "overwrite": True,
                "consolidation": True,
            }
        )
    )
    opts = build_argument_parser().parse_args(
        [
            "-c", str(config_path),
            "--alpha-dct", "0.01",
            "--alpha-perm", "0.02",
            "--p-value-correction", "fdr_bh",
            "--diagnostic",
            "--species-filter", "species.txt",
            "--no-overwrite",
            "--no-consolidation",
            "--no-bootstrap",
            "--preflight-data-check",
            "--preflight-triplet-cap", "0",
            "--processes", "3",
            "--seed", "5",
        ]
    )

    config = resolve_config(opts)
    assert config["species_tree"].endswith("file_species.tree")
    assert config["alpha_ks"] == 0.04
    assert config["alpha_dct"] == 0.01
    assert config["alpha_perm"] == 0.02
    assert config["p_value_correction"] == "fdr_bh"
    assert config["diagnostic"] is True
    assert config["species_filter"].endswith("species.txt")
    assert config["triplet_filter"] is None
    assert config["overwrite"] is False
    assert config["consolidation"] is False
    assert config["bootstrap"] is False
    assert config["preflight_data_check"] is True
    assert config["preflight_triplet_cap"] == 0
    assert (config["processes"], config["seed"]) == (3, 5)


@pytest.mark.parametrize(
    "value, expected",
    [
        ("OUT", ["OUT"]),
        ("Out1,Out2", ["Out1", "Out2"]),
        (" Out1 , Out2 ,", ["Out1", "Out2"]),
        (["Out1,Out2", "Out3"], ["Out1", "Out2", "Out3"]),
        (None, ConfigError),
        ("  ", ConfigError),
        ([""], ConfigError),
        (42, ConfigError),
    ],
)
def test_outgroup_key_is_normalized_or_rejected(value, expected, tmp_path):
    """One `outgroup` key covers a label, a comma-separated string or a list; anything label-free is an error.

    Whitespace and empty entries are dropped from what is accepted; a missing
    value, blanks alone or a non-string are a config error naming the key,
    never an empty list.
    """
    config_path = _payload(tmp_path, outgroup=value)
    if expected is ConfigError:
        with pytest.raises(ConfigError, match="outgroup"):
            load_orchestrator_config(str(config_path))
        return
    assert load_orchestrator_config(str(config_path))["outgroup"] == expected


@pytest.mark.parametrize(
    "extra, match",
    [
        ({"species_tree_path": None}, "species_tree_path"),
        ({"p_value_correction": True}, "must be one of"),
        ({"diagnostic": "yes"}, "diagnostic"),
        ({"seed": "abc"}, "seed"),
        ({"preflight_triplet_cap": -1}, "preflight_triplet_cap"),
        ({"species_rename_map": "absent.tsv"}, "species_rename_map"),
        (
            {"triplet_filter": "triplets.txt", "species_filter": "species.txt"},
            "cannot both be set",
        ),
    ],
)
def test_invalid_values_are_rejected_by_field_name(tmp_path, extra, match):
    """A missing required key or an out-of-domain value fails naming the field.

    One case per validator shape: a required path, a choice list (fed a
    boolean with no matching choice, which is what YAML makes of a bare
    ``yes``), an optional bool, an optional int, a non-negative int whose
    ``0`` already means "no cap", a path whose file is read when the config
    resolves, and the two triplet-selection paths that exclude each other.
    """
    with pytest.raises(ConfigError, match=match):
        load_orchestrator_config(str(_payload(tmp_path, **extra)))


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
