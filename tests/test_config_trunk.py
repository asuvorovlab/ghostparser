"""Tests for the shared configuration trunk in ``ghostparser.config``.

The trunk holds only the helpers whose behaviour is identical for every
consumer (the orchestrator, the ML subpackage, and the introgression mapper):
``ConfigError``, path resolution, raw config-file loading, required-path
validation, overwrite-flag resolution, and output-directory preparation.
Module-specific defaults and validators live in each module's own config.
"""

import json
from pathlib import Path

import pytest

from ghostparser.config import (
    DEFAULT_OVERWRITE,
    ConfigError,
    _load_raw_config,
    _resolve_path,
    _validate_overwrite_flag,
    _validate_required_path,
    prepare_output_directory,
)

pytestmark = pytest.mark.config


def test_resolve_path_handles_absolute_relative_and_home(tmp_path, monkeypatch):
    """Absolute paths pass through; relative resolve from cwd; ``~`` expands."""
    absolute = tmp_path / "species.nwk"
    assert _resolve_path(str(absolute)) == str(absolute.resolve())

    monkeypatch.chdir(tmp_path)
    assert _resolve_path("genes.nwk") == str((tmp_path / "genes.nwk").resolve())

    home_resolved = _resolve_path("~/data.nwk")
    assert home_resolved == str(Path.home() / "data.nwk")
    assert "~" not in home_resolved


def test_load_raw_config_reads_json_and_yaml_and_rejects_the_rest(tmp_path):
    """JSON and YAML files load into one mapping; missing files, unknown suffixes and non-mapping roots are rejected."""
    payload = {"input_path": "data.tsv", "cv_folds": 5}
    for suffix in (".json", ".yaml", ".yml"):
        config_path = tmp_path / f"config{suffix}"
        if suffix == ".json":
            config_path.write_text(json.dumps(payload))
        else:
            config_path.write_text("input_path: data.tsv\ncv_folds: 5\n")
        assert _load_raw_config(str(config_path)) == payload, suffix

    with pytest.raises(FileNotFoundError):
        _load_raw_config(str(tmp_path / "absent.json"))

    unsupported = tmp_path / "config.txt"
    unsupported.write_text("input_path: data.tsv")
    with pytest.raises(ConfigError, match="must be .json, .yaml, or .yml"):
        _load_raw_config(str(unsupported))

    non_mapping = tmp_path / "list.json"
    non_mapping.write_text(json.dumps(["a", "b"]))
    with pytest.raises(ConfigError, match="key/value object"):
        _load_raw_config(str(non_mapping))


def test_validate_required_path_resolves_or_raises(tmp_path, monkeypatch):
    """A present path is resolved; a missing or blank one raises ConfigError."""
    monkeypatch.chdir(tmp_path)
    assert _validate_required_path({"species_tree_path": "s.nwk"}, "species_tree_path") == str(
        (tmp_path / "s.nwk").resolve()
    )

    for payload in ({}, {"species_tree_path": ""}, {"species_tree_path": "   "}):
        with pytest.raises(ConfigError, match="Missing required config field"):
            _validate_required_path(payload, "species_tree_path")


def test_validate_overwrite_flag_precedence_and_validation():
    """``overwrite`` wins over ``no_overwrite``; both must be booleans."""
    assert _validate_overwrite_flag({}) is DEFAULT_OVERWRITE
    assert _validate_overwrite_flag({"overwrite": False}) is False
    assert _validate_overwrite_flag({"no_overwrite": True}) is False
    assert _validate_overwrite_flag({"no_overwrite": False}) is True
    # The canonical key takes precedence over the CLI-style negated key.
    assert _validate_overwrite_flag({"overwrite": True, "no_overwrite": True}) is True
    # An explicit default is honoured when neither key is present.
    assert _validate_overwrite_flag({}, default=False) is False

    with pytest.raises(ConfigError, match="overwrite must be a boolean"):
        _validate_overwrite_flag({"overwrite": "yes"})
    with pytest.raises(ConfigError, match="no_overwrite must be a boolean"):
        _validate_overwrite_flag({"no_overwrite": "yes"})


def test_prepare_output_directory_overwrites_or_suffixes(tmp_path):
    """Overwrite resets the directory; disabling it picks the smallest free suffix."""
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    stale_file = results_dir / "stale.txt"
    stale_file.write_text("old contents")

    prepared = prepare_output_directory(results_dir)
    assert prepared == str(results_dir.resolve())
    assert results_dir.exists()
    assert not stale_file.exists()

    (results_dir / "fresh.txt").write_text("fresh contents")
    (tmp_path / "results_1").mkdir()
    (tmp_path / "results_3").mkdir()

    # results_1 and results_3 are taken, so the smallest missing suffix is 2.
    suffixed = prepare_output_directory(results_dir, overwrite=False)
    assert suffixed == str((tmp_path / "results_2").resolve())
    assert (tmp_path / "results_2").exists()
    assert (results_dir / "fresh.txt").exists()

    # A path whose parents do not exist yet is created on demand.
    nested = tmp_path / "a" / "b" / "results"
    assert prepare_output_directory(nested) == str(nested.resolve())
    assert nested.is_dir()
