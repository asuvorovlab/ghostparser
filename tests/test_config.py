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
                "p_value_correction": "bonferroni",
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
    assert config["p_value_correction"] == "bonferroni"
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


def test_load_orchestrator_config_invalid_discordant_test(tmp_path):
    config_path = tmp_path / "bad_test_method.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "discordant_test": "invalid",
            }
        )
    )

    with pytest.raises(ConfigError, match="discordant_test"):
        load_orchestrator_config(str(config_path))


def test_load_orchestrator_config_invalid_summary_statistic(tmp_path):
    config_path = tmp_path / "bad_summary.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "summary_statistic": "invalid-summary",
            }
        )
    )

    with pytest.raises(ConfigError, match="summary_statistic"):
        load_orchestrator_config(str(config_path))


def test_load_orchestrator_config_invalid_stats_backend(tmp_path):
    config_path = tmp_path / "bad_stats_backend.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "stats_backend": "numpy",
            }
        )
    )

    with pytest.raises(ConfigError, match="stats_backend"):
        load_orchestrator_config(str(config_path))


def test_load_orchestrator_config_invalid_tree_height_calculation_strategy(tmp_path):
    config_path = tmp_path / "bad_tree_height_strategy.json"
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "species.nwk",
                "gene_trees_path": "genes.nwk",
                "outgroup": "OutA",
                "tree_height_calculation_strategy": "D",
            }
        )
    )

    with pytest.raises(ConfigError, match="tree_height_calculation_strategy"):
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
                "p_value_correction": "none",
                "processes": 3,
                "no_multiprocessing": False,
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
    assert config["p_value_correction"] == "none"
    assert config["processes"] == 3
    assert config["no_multiprocessing"] is False


def test_load_triplet_processor_config_invalid_p_value_correction(tmp_path):
    config_path = tmp_path / "triplet_processor_bad_p_value_correction.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "unique_triplets_gene_trees.txt",
                "p_value_correction": "holm",
            }
        )
    )

    with pytest.raises(ConfigError, match="p_value_correction"):
        load_triplet_processor_config(str(config_path))


def test_load_triplet_processor_config_invalid_stats_backend(tmp_path):
    config_path = tmp_path / "triplet_processor_bad_stats_backend.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "unique_triplets_gene_trees.txt",
                "stats_backend": "numpy",
            }
        )
    )

    with pytest.raises(ConfigError, match="stats_backend"):
        load_triplet_processor_config(str(config_path))


def test_load_triplet_processor_config_invalid_tree_height_calculation_strategy(tmp_path):
    config_path = tmp_path / "triplet_processor_bad_tree_height_strategy.json"
    config_path.write_text(
        json.dumps(
            {
                "input_path": "unique_triplets_gene_trees.txt",
                "tree_height_calculation_strategy": "D",
            }
        )
    )

    with pytest.raises(ConfigError, match="tree_height_calculation_strategy"):
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
    assert config["p_value_correction"] == "bonferroni"


def test_path_resolution_absolute_paths(tmp_path):
    """Test that absolute paths are preserved as-is."""
    config_path = tmp_path / "absolute_paths.json"
    abs_species = "/absolute/path/to/species.nwk"
    abs_genes = "/absolute/path/to/genes.nwk"
    abs_output = "/absolute/path/to/output"
    
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": abs_species,
                "gene_trees_path": abs_genes,
                "outgroup": "OutA",
                "output_folder": abs_output,
            }
        )
    )
    
    config = load_orchestrator_config(str(config_path))
    assert config["species_tree"] == abs_species
    assert config["gene_trees"] == abs_genes
    assert config["output"] == abs_output


def test_path_resolution_relative_paths(tmp_path):
    """Test that relative paths are resolved from current working directory."""
    config_path = tmp_path / "relative_paths.json"
    
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "data/species.nwk",
                "gene_trees_path": "./genes.nwk",
                "outgroup": "OutA",
                "output_folder": "results",
            }
        )
    )
    
    config = load_orchestrator_config(str(config_path))
    # Relative paths should be resolved from current working directory
    assert config["species_tree"] == str(Path("data/species.nwk").resolve())
    assert config["gene_trees"] == str(Path("./genes.nwk").resolve())
    assert config["output"] == str(Path("results").resolve())


def test_path_resolution_home_directory(tmp_path):
    """Test that ~ is expanded to user home directory."""
    config_path = tmp_path / "home_paths.json"
    
    config_path.write_text(
        json.dumps(
            {
                "species_tree_path": "~/data/species.nwk",
                "gene_trees_path": "~/data/genes.nwk",
                "outgroup": "OutA",
            }
        )
    )
    
    config = load_orchestrator_config(str(config_path))
    # ~ should be expanded to home directory
    assert config["species_tree"] == str(Path("~/data/species.nwk").expanduser().resolve())
    assert config["gene_trees"] == str(Path("~/data/genes.nwk").expanduser().resolve())

