"""Tests for configuration parsing."""

import json
from pathlib import Path

import pytest

from ghostparser.config import (
    ConfigError,
    load_orchestrator_config,
)


def test_load_orchestrator_config_json(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroups": ["OutA", "OutB"],
                "output_folder": "out",
                "processes": 4,
                "generate_summary_stats": True,
                "triplet_filter": "triplets.txt",
                "min_support_value": 0.7,
                "discordant_test": "z-test",
                "summary_statistic": "median",
                "stats_backend": "standard",
                "tree_height_calculation_strategy": "B",
                "p_value_correction": "bfn",
                "alpha_dct": 0.02,
                "alpha_ks": 0.1,
            }
        )
    )

    config = load_orchestrator_config(str(config_path))

    # Paths are resolved to absolute paths
    assert config["species_tree"] == str(Path("species.nwk").resolve())
    assert config["gene_trees"] == str(Path("genes.nwk").resolve())
    assert config["outgroup"] == ["OutA", "OutB"]
    assert config["output"] == str(Path("out").resolve())
    assert config["processes"] == 4
    assert config["generate_summary_stats"] is True
    assert config["consolidation"] is True
    assert config["triplet_filter"] == str(Path("triplets.txt").resolve())
    assert config["min_support_value"] == 0.7
    assert config["discordant_test"] == "z-test"
    assert config["summary_statistic"] == "median"
    assert config["stats_backend"] == "standard"
    assert config["tree_height_calculation_strategy"] == "B"
    assert config["p_value_correction"] == "bfn"
    assert config["alpha_dct"] == 0.02
    assert config["alpha_ks"] == 0.1


def test_load_orchestrator_config_yaml(tmp_path):
    yaml = pytest.importorskip("yaml")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA,OutB",
            }
        )
    )

    config = load_orchestrator_config(str(config_path))
    assert config["outgroup"] == ["OutA", "OutB"]


def test_load_orchestrator_config_single_outgroup_string_is_single_taxon(tmp_path):
    config_path = tmp_path / "single_outgroup.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
            }
        )
    )

    config = load_orchestrator_config(str(config_path))
    assert config["outgroup"] == ["OutA"]


def test_load_orchestrator_config_missing_required(tmp_path):
    config_path = tmp_path / "bad.json"
    config_path.write_text(json.dumps({"species_tree_path": "species.nwk"}))

    with pytest.raises(ConfigError, match="Missing required config field"):
        load_orchestrator_config(str(config_path))


@pytest.mark.parametrize(
    "field,bad_value",
    [
        ("discordant_test", "invalid"),
        ("summary_statistic", "invalid-summary"),
        ("stats_backend", "numpy"),
    ],
)
def test_load_orchestrator_config_invalid_choice_fields(tmp_path, field, bad_value):
    config_path = tmp_path / f"bad_{field}.json"
    payload = {
        "species_tree_path": "species.nwk",
        "gene_trees_path": "genes.nwk",
        "outgroup": "OutA",
        field: bad_value,
    }
    config_path.write_text(json.dumps(payload))

    with pytest.raises(ConfigError, match=field):
        load_orchestrator_config(str(config_path))


@pytest.mark.parametrize(
    "strategy,should_raise",
    [
        ("AVG", False),
        ("A", False),
        ("B", False),
        ("C", False),
        ("SIS", False),
        ("INT", False),
        ("D", True),
    ],
)
def test_load_orchestrator_config_tree_height_strategy_validation(tmp_path, strategy, should_raise):
    config_path = tmp_path / f"orchestrator_tree_height_{strategy}.json"
    payload = {
        "species_tree_path": "species.nwk",
        "gene_trees_path": "genes.nwk",
        "outgroup": "OutA",
        "tree_height_calculation_strategy": strategy,
    }
    config_path.write_text(json.dumps(payload))

    if should_raise:
        with pytest.raises(ConfigError, match="tree_height_calculation_strategy"):
            load_orchestrator_config(str(config_path))
    else:
        config = load_orchestrator_config(str(config_path))
        assert config["tree_height_calculation_strategy"] == strategy


def test_load_orchestrator_config_defaults_processes_to_zero(tmp_path):
    config_path = tmp_path / "config_default_processes.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
            }
        )
    )

    config = load_orchestrator_config(str(config_path))
    assert config["processes"] == 0
    assert config["generate_summary_stats"] is False
    assert config["bootstrap"] is True
    assert config["consolidation"] is True
    assert config["bootstrap_options"] == {
        "iterations": 100,
        "seed": None,
        "debug_mode": False,
        "summary_only": False,
    }


def test_load_orchestrator_config_allows_disabling_consolidation(tmp_path):
    config_path = tmp_path / "config_consolidation_false.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "consolidation": False,
            }
        )
    )

    config = load_orchestrator_config(str(config_path))
    assert config["consolidation"] is False


@pytest.mark.parametrize(
    "species_path,genes_path,output_path",
    [
        ("/absolute/path/to/species.nwk", "/absolute/path/to/genes.nwk", "/absolute/path/to/output"),
        ("data/species.nwk", "./genes.nwk", "results"),
        ("~/data/species.nwk", "~/data/genes.nwk", None),
    ],
)
def test_path_resolution_for_absolute_relative_and_home_paths(tmp_path, species_path, genes_path, output_path):
    config_path = tmp_path / "path_resolution.json"
    payload = {
        "species_tree_path": species_path,
        "gene_trees_path": genes_path,
        "outgroup": "OutA",
    }
    if output_path is not None:
        payload["output_folder"] = output_path

    config_path.write_text(json.dumps(payload))
    config = load_orchestrator_config(str(config_path))

    assert config["species_tree"] == str(Path(species_path).expanduser().resolve())
    assert config["gene_trees"] == str(Path(genes_path).expanduser().resolve())
    if output_path is not None:
        assert config["output"] == str(Path(output_path).expanduser().resolve())

