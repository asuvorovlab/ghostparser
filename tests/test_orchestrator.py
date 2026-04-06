"""Tests for orchestrator module."""

import argparse
from pathlib import Path

from ghostparser.orchestrator import (
    _resolve_parallel_mode,
    _resolve_runtime_args,
    _resolve_processes,
)


def test_resolve_processes_zero_uses_all_cores(monkeypatch):
    monkeypatch.setattr("ghostparser.orchestrator.cpu_count", lambda: 16)
    assert _resolve_processes(0) == 16
    assert _resolve_processes(None) is None
    assert _resolve_processes(4) == 4


def test_resolve_parallel_mode(monkeypatch):
    monkeypatch.setattr("ghostparser.orchestrator.cpu_count", lambda: 8)

    processes, use_multiprocessing = _resolve_parallel_mode(0)
    assert processes == 8
    assert use_multiprocessing is True

    processes, use_multiprocessing = _resolve_parallel_mode(1)
    assert processes == 1
    assert use_multiprocessing is False

    processes, use_multiprocessing = _resolve_parallel_mode(4)
    assert processes == 4
    assert use_multiprocessing is True


def test_resolve_runtime_args_cli_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(
        config_file=None,
        species_tree_path="species.nwk",
        gene_trees_path="genes.nwk",
        outgroups="Out1,Out2",
        triplet_filter=None,
        output_folder=str(tmp_path / "results"),
        processes=0,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        tree_height_calculation_strategy=None,
        p_value_correction=None,
        alpha_dct=None,
        alpha_ks=None,
    )

    resolved = _resolve_runtime_args(args)

    # Paths are resolved to absolute paths
    assert resolved.species_tree == str((tmp_path / "species.nwk").resolve())
    assert resolved.gene_trees == str((tmp_path / "genes.nwk").resolve())
    assert resolved.outgroup == ["Out1", "Out2"]
    assert resolved.output == str(tmp_path / "results")
    assert resolved.processes == 0
    assert resolved.min_support_value == 0.5
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.bootstrap is False
    assert resolved.bootstrap_options == {"iterations": 100, "seed": None, "summary_only": True}


def test_resolve_runtime_args_cli_custom_processes_preserved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(
        config_file=None,
        species_tree_path="species.nwk",
        gene_trees_path="genes.nwk",
        outgroups="Out1,Out2",
        triplet_filter=None,
        output_folder=str(tmp_path / "results"),
        processes=5,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        tree_height_calculation_strategy=None,
        p_value_correction=None,
        alpha_dct=None,
        alpha_ks=None,
        bootstrap=True,
        bootstrap_iterations=25,
        bootstrap_seed=99,
        bootstrap_summary_only=False,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 5
    assert resolved.bootstrap is True
    assert resolved.bootstrap_options == {"iterations": 25, "seed": 99, "summary_only": False}


def test_resolve_runtime_args_config_with_cli_warns_and_ignores(tmp_path, capsys):
    config_path = tmp_path / "run_config.json"
    config_path.write_text(
        '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA"}'
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        species_tree_path="species.nwk",
        gene_trees_path=None,
        outgroups=None,
        triplet_filter=None,
        output_folder=str(tmp_path / "results"),
        processes=7,
        min_support_value=0.9,
        discordant_test="chi-square",
        summary_statistic="mean",
        stats_backend="standard",
        tree_height_calculation_strategy="C",
            p_value_correction="no",
        alpha_dct=0.2,
        alpha_ks=0.3,
        bootstrap=True,
        bootstrap_iterations=20,
        bootstrap_seed=1,
        bootstrap_summary_only=False,
    )

    resolved = _resolve_runtime_args(args)
    captured = capsys.readouterr()

    assert "Warning: --config-file provided; CLI arguments not in config will be ignored" in captured.out
    # Paths are resolved to absolute paths
    assert resolved.species_tree == str(Path("s.nwk").resolve())
    assert resolved.gene_trees == str(Path("g.nwk").resolve())
    assert resolved.outgroup == ["OutA"]
    assert resolved.processes == 0
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.bootstrap is False
    assert resolved.bootstrap_options == {"iterations": 100, "seed": None, "summary_only": True}


def test_resolve_runtime_args_config_processes_preserved_when_set(tmp_path):
    config_path = tmp_path / "run_config_with_processes.json"
    config_path.write_text(
        '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA","processes":4}'
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        species_tree_path=None,
        gene_trees_path=None,
        outgroups=None,
        triplet_filter=None,
        output_folder=str(tmp_path / "results"),
        processes=None,
        min_support_value=None,
        discordant_test=None,
        summary_statistic=None,
        stats_backend=None,
        tree_height_calculation_strategy=None,
        p_value_correction=None,
        alpha_dct=None,
        alpha_ks=None,
        bootstrap=None,
        bootstrap_iterations=None,
        bootstrap_seed=None,
        bootstrap_summary_only=None,
    )

    resolved = _resolve_runtime_args(args)
    assert resolved.processes == 4
