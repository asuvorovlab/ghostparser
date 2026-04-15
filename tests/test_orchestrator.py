"""Tests for orchestrator module."""

import argparse
from pathlib import Path
from types import SimpleNamespace

import pytest

import ghostparser.orchestrator as orchestrator_module

from ghostparser.orchestrator import (
    _resolve_parallel_mode,
    _resolve_runtime_args,
    _resolve_processes,
)


def _orchestrator_args(tmp_path, **overrides):
    base = {
        "config_file": None,
        "species_tree_path": "species.nwk",
        "gene_trees_path": "genes.nwk",
        "outgroups": "Out1,Out2",
        "triplet_filter": None,
        "output_folder": str(tmp_path / "results"),
        "processes": 0,
        "generate_summary_stats": None,
        "min_support_value": None,
        "discordant_test": None,
        "summary_statistic": None,
        "stats_backend": None,
        "tree_height_calculation_strategy": None,
        "p_value_correction": None,
        "alpha_dct": None,
        "alpha_ks": None,
        "bootstrap": False,
        "bootstrap_iterations": None,
        "bootstrap_seed": None,
        "bootstrap_debug_mode": None,
        "bootstrap_summary_only": None,
    }
    base.update(overrides)
    return argparse.Namespace(**base)


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


@pytest.mark.parametrize(
    "cli_overrides,expected_processes,expected_bootstrap,expected_bootstrap_options",
    [
        (
            {},
            0,
            False,
            {
                "iterations": 100,
                "seed": None,
                "debug_mode": False,
                "summary_only": False,
            },
        ),
        (
            {
                "processes": 5,
                "bootstrap": True,
                "bootstrap_iterations": 25,
                "bootstrap_seed": 99,
                "bootstrap_debug_mode": True,
                "bootstrap_summary_only": False,
            },
            5,
            True,
            {
                "iterations": 25,
                "seed": 99,
                "debug_mode": True,
                "summary_only": False,
            },
        ),
    ],
)
def test_resolve_runtime_args_cli_defaults_and_overrides(
    tmp_path,
    monkeypatch,
    cli_overrides,
    expected_processes,
    expected_bootstrap,
    expected_bootstrap_options,
):
    monkeypatch.chdir(tmp_path)
    resolved = _resolve_runtime_args(_orchestrator_args(tmp_path, **cli_overrides))

    # Paths are resolved to absolute paths
    assert resolved.species_tree == str((tmp_path / "species.nwk").resolve())
    assert resolved.gene_trees == str((tmp_path / "genes.nwk").resolve())
    assert resolved.outgroup == ["Out1", "Out2"]
    assert resolved.output == str(tmp_path / "results")
    assert resolved.processes == expected_processes
    assert resolved.generate_summary_stats is False
    assert resolved.min_support_value == 0.5
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.bootstrap is expected_bootstrap
    assert resolved.bootstrap_options == expected_bootstrap_options


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
        generate_summary_stats=True,
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
        bootstrap_debug_mode=True,
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
    assert resolved.generate_summary_stats is False
    assert resolved.discordant_test == "chi-square"
    assert resolved.summary_statistic == "median"
    assert resolved.stats_backend == "standard"
    assert resolved.tree_height_calculation_strategy == "AVG"
    assert resolved.p_value_correction == "no"
    assert resolved.alpha_dct == 0.01
    assert resolved.alpha_ks == 0.05
    assert resolved.bootstrap is False
    assert resolved.bootstrap_options == {
        "iterations": 100,
        "seed": None,
        "debug_mode": False,
        "summary_only": False,
    }


@pytest.mark.parametrize(
    "config_payload,expected_processes",
    [
        (
            '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA","processes":4}',
            4,
        ),
        (
            '{"species_tree_path":"s.nwk","gene_trees_path":"g.nwk","outgroup":"OutA"}',
            0,
        ),
    ],
)
def test_resolve_runtime_args_config_processes_behavior(tmp_path, config_payload, expected_processes):
    config_path = tmp_path / "run_config_processes.json"
    config_path.write_text(config_payload)

    resolved = _resolve_runtime_args(
        _orchestrator_args(
            tmp_path,
            config_file=str(config_path),
            species_tree_path=None,
            gene_trees_path=None,
            outgroups=None,
            processes=None,
        )
    )
    assert resolved.processes == expected_processes


def _runtime_args(tmp_path):
    species_path = tmp_path / "species.nwk"
    genes_path = tmp_path / "genes.nwk"
    species_path.write_text("(A:1,B:1);\n")
    genes_path.write_text("(A:1,B:1);\n")

    return SimpleNamespace(
        species_tree=str(species_path),
        gene_trees=str(genes_path),
        outgroup=["Out1"],
        triplet_filter=None,
        output=str(tmp_path / "results"),
        processes=1,
        generate_summary_stats=False,
        min_support_value=0.5,
        discordant_test="chi-square",
        summary_statistic="median",
        stats_backend="standard",
        tree_height_calculation_strategy="AVG",
        p_value_correction="no",
        alpha_dct=0.01,
        alpha_ks=0.05,
        bootstrap=False,
        bootstrap_options={
            "iterations": 100,
            "seed": None,
            "debug_mode": False,
            "summary_only": False,
        },
    )


def _patch_orchestrator_runtime_dependencies(monkeypatch):
    class _DummyParser:
        def parse_args(self):
            return SimpleNamespace()

    monkeypatch.setattr(orchestrator_module, "_build_argument_parser", lambda: _DummyParser())

    monkeypatch.setattr(orchestrator_module, "_parse_outgroup_arg", lambda outgroup: ["Out1"])
    monkeypatch.setattr(orchestrator_module, "clean_and_save_trees", lambda *args, **kwargs: (["species_tree"], {}))
    monkeypatch.setattr(orchestrator_module, "read_tree_file", lambda *_args, **_kwargs: ["species_tree"])
    monkeypatch.setattr(orchestrator_module, "get_taxa_from_tree", lambda *_args, **_kwargs: ["A", "B", "C", "Out1"])
    monkeypatch.setattr(
        orchestrator_module,
        "_root_tree_on_outgroup",
        lambda *_args, **_kwargs: ("pruned_tree", set(), set(), {"A", "B", "C"}),
    )
    monkeypatch.setattr(orchestrator_module, "write_clean_trees", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(orchestrator_module, "generate_triplets", lambda *_args, **_kwargs: [("A", "B", "C")])
    monkeypatch.setattr(orchestrator_module, "format_newick_with_precision", lambda *_args, **_kwargs: "((A:1,B:1):1,C:1);")
    monkeypatch.setattr(orchestrator_module.dendropy.Tree, "get", lambda *_args, **_kwargs: "species_dendro_tree")
    monkeypatch.setattr(
        orchestrator_module,
        "_build_species_triplet_metadata",
        lambda *_args, **_kwargs: ([("A", "B", "C")], {("A", "B", "C"): "((A:1,B:1):1,C:1);"}, []),
    )
    monkeypatch.setattr(
        orchestrator_module,
        "clean_and_save_gene_trees",
        lambda *_args, **_kwargs: (["gene_tree"], {}, 1, []),
    )
    monkeypatch.setattr(orchestrator_module, "write_pipeline_results", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(orchestrator_module, "write_summary_statistics_tsv", lambda *_args, **_kwargs: None)


def test_main_uses_file_backed_pipeline(tmp_path, monkeypatch):
    args = _runtime_args(tmp_path)
    _patch_orchestrator_runtime_dependencies(monkeypatch)
    monkeypatch.setattr(orchestrator_module, "_resolve_runtime_args", lambda _parsed: args)

    calls = {"file_extract": 0, "file_infer": 0}

    def _file_extract_stub(*_args, **_kwargs):
        calls["file_extract"] += 1
        return (3, 1, 1)

    def _file_infer_stub(*_args, **_kwargs):
        calls["file_infer"] += 1
        return [SimpleNamespace()]

    monkeypatch.setattr(orchestrator_module, "write_triplet_gene_trees_multiprocess", _file_extract_stub)
    monkeypatch.setattr(orchestrator_module, "analyze_triplet_gene_tree_file", _file_infer_stub)

    orchestrator_module.main()

    assert calls["file_extract"] == 1
    assert calls["file_infer"] == 1
