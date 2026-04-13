"""Tests for configuration parsing."""

import json
from pathlib import Path

import pytest

from ghostparser.config import (
    ConfigError,
    load_orchestrator_config,
    load_tree_parser_config,
    load_triplet_processor_config,
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
        ("tree_height_calculation_strategy", "D"),
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
    assert config["bootstrap"] is False
    assert config["bootstrap_options"] == {
        "iterations": 100,
        "seed": None,
        "debug_mode": False,
        "summary_only": False,
    }


def test_load_tree_parser_config_json(tmp_path):
    config_path = tmp_path / "tree_parser_config.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroups": ["OutA", "OutB"],
                "output_folder": "out",
                "processes": 2,
                "triplet_filter": "triplets.txt",
                "min_support_value": 0.6,
                "no_multiprocessing": True,
            }
        )
    )

    config = load_tree_parser_config(str(config_path))

    # Paths are resolved to absolute paths
    assert config["species_tree"] == str(Path("species.nwk").resolve())
    assert config["gene_trees"] == str(Path("genes.nwk").resolve())
    assert config["outgroup"] == ["OutA", "OutB"]
    assert config["output"] == str(Path("out").resolve())
    assert config["processes"] == 2
    assert config["triplet_filter"] == str(Path("triplets.txt").resolve())
    assert config["min_support_value"] == 0.6
    assert config["no_multiprocessing"] is True


def test_load_tree_parser_config_invalid_no_multiprocessing(tmp_path):
    config_path = tmp_path / "tree_parser_bad.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "no_multiprocessing": "yes",
            }
        )
    )

    with pytest.raises(ConfigError, match="no_multiprocessing"):
        load_tree_parser_config(str(config_path))


def test_load_tree_parser_config_defaults_processes_to_zero(tmp_path):
    config_path = tmp_path / "tree_parser_default_processes.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
            }
        )
    )

    config = load_tree_parser_config(str(config_path))
    assert config["processes"] == 0


def test_load_triplet_processor_config_json(tmp_path):
    config_path = tmp_path / "triplet_processor_config.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "unique_triplets_gene_trees.txt",
                "output_path": "results.tsv",
                "stats_output": "stats.json",
                "alpha_dct": 0.02,
                "alpha_ks": 0.1,
                "discordant_test": "z-test",
                "summary_statistic": "median",
                "stats_backend": "standard",
                "tree_height_calculation_strategy": "C",
                "p_value_correction": "no",
                "processes": 3,
                "no_multiprocessing": False,
                "bootstrap": True,
                "bootstrap_options": {
                    "iterations": 15,
                    "seed": 42,
                    "debug_mode": True,
                    "summary_only": False,
                },
            }
        )
    )

    config = load_triplet_processor_config(str(config_path))

    # Paths are resolved to absolute paths
    assert config["input"] == str(Path("unique_triplets_gene_trees.txt").resolve())
    assert config["output"] == str(Path("results.tsv").resolve())
    assert config["stats_output"] == str(Path("stats.json").resolve())
    assert config["alpha_dct"] == 0.02
    assert config["alpha_ks"] == 0.1
    assert config["discordant_test"] == "z-test"
    assert config["summary_statistic"] == "median"
    assert config["stats_backend"] == "standard"
    assert config["tree_height_calculation_strategy"] == "C"
    assert config["p_value_correction"] == "no"
    assert config["processes"] == 3
    assert config["no_multiprocessing"] is False
    assert config["bootstrap"] is True
    assert config["bootstrap_options"] == {
        "iterations": 15,
        "seed": 42,
        "debug_mode": True,
        "summary_only": False,
    }


@pytest.mark.parametrize(
    "payload,match",
    [
        ({"bootstrap": "yes"}, "bootstrap"),
        ({"bootstrap_options": {"iterations": 0}}, "iterations"),
        ({"bootstrap_options": {"seed": "abc"}}, "seed"),
        ({"bootstrap_options": {"debug_mode": "no"}}, "debug_mode"),
        ({"bootstrap_options": {"summary_only": "no"}}, "summary_only"),
    ],
)
def test_load_triplet_processor_config_invalid_bootstrap_fields(tmp_path, payload, match):
    config_path = tmp_path / "triplet_processor_bad_bootstrap.json"
    base = {"input_path": "unique_triplets_gene_trees.txt"}
    base.update(payload)
    config_path.write_text(json.dumps(base))

    with pytest.raises(ConfigError, match=match):
        load_triplet_processor_config(str(config_path))


@pytest.mark.parametrize(
    "field,bad_value",
    [
        ("p_value_correction", "sidak"),
        ("stats_backend", "numpy"),
        ("tree_height_calculation_strategy", "D"),
    ],
)
def test_load_triplet_processor_config_invalid_choice_fields(tmp_path, field, bad_value):
    config_path = tmp_path / f"triplet_processor_bad_{field}.json"
    payload = {
        "input_path": "unique_triplets_gene_trees.txt",
        field: bad_value,
    }
    config_path.write_text(json.dumps(payload))

    with pytest.raises(ConfigError, match=field):
        load_triplet_processor_config(str(config_path))


def test_load_triplet_processor_config_missing_input(tmp_path):
    config_path = tmp_path / "triplet_processor_bad.json"
    config_path.write_text(json.dumps({"alpha_dct": 0.02}))

    with pytest.raises(ConfigError, match="input_path"):
        load_triplet_processor_config(str(config_path))


def test_load_triplet_processor_config_defaults_processes_to_zero(tmp_path):
    config_path = tmp_path / "triplet_processor_default_processes.json"
    config_path.write_text(json.dumps({"input_path": "input.tsv"}))

    config = load_triplet_processor_config(str(config_path))
    assert config["processes"] == 0
    assert config["tree_height_calculation_strategy"] == "AVG"
    assert config["p_value_correction"] == "no"


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

