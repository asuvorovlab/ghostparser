"""Config resolution for the ML trainers.

Marked ``config`` throughout. Individual defaults are not pinned: CONFIG.md
and the shipped sample configs state them; what is pinned is that the loader
resolves, that explicit values win, and that the two estimator passthrough keys
accept every form scikit-learn does and nothing else.
"""

import json
from pathlib import Path

import pytest

from ghostparser.config import ConfigError
from ghostparser.ml.config import (
    build_trainer_argument_parser,
    load_ml_config,
    resolve_trainer_runtime_args,
)

pytestmark = pytest.mark.config

_SAMPLE_DIR = Path(__file__).resolve().parents[1] / "sample_configs"

_REQUIRED = {
    "input_path": "./results/summary_statistics.tsv",
    "output_dir": "./results/ml_out",
}


def _load(tmp_path, **extra):
    """Write a config holding the required keys plus ``extra`` and load it.

    Args:
        tmp_path: pytest temporary directory.
        **extra: Additional payload keys.

    Returns:
        The resolved config dict.
    """
    config_path = tmp_path / "ml_config.json"
    config_path.write_text(json.dumps({**_REQUIRED, **extra}))
    return load_ml_config(str(config_path))


@pytest.mark.parametrize(
    "sample_name", ["random_forest_minimal.yaml", "multi_knn_minimal.yaml"]
)
def test_shipped_trainer_sample_configs_resolve(sample_name):
    """The trainer samples load cleanly and use the current key names.

    Guards against a sample drifting out of step with the validator, which
    would leave users copying a config the loader rejects.
    """
    config = load_ml_config(str(_SAMPLE_DIR / sample_name))

    assert config["input_path"].endswith("summary_statistics.tsv")
    assert config["target_column"]


@pytest.mark.parametrize(
    "payload, key, expected",
    [
        ({"overwrite": False}, "overwrite", False),
        ({"target_column": "label"}, "target_column", "label"),
        ({"model": {"n_estimators": 25}}, "n_estimators", 25),
        ({"model": {"min_samples_leaf": 4}}, "min_samples_leaf", 4),
    ],
)
def test_ml_config_explicit_values_win_over_defaults(payload, key, expected, tmp_path):
    """Explicit values override the defaults, from the top level or ``model``.

    The last two also pin the nested-block flattening: hyperparameters given
    under ``model`` surface at the top level of the resolved config.
    """
    assert _load(tmp_path, **payload)[key] == expected


@pytest.mark.parametrize(
    "key, value, expected",
    [
        ("max_features", None, None),
        ("max_features", "log2", "log2"),
        ("max_features", 3, 3),
        ("max_features", 0.5, 0.5),
        ("class_weight", None, None),
        ("class_weight", "balanced", "balanced"),
        ("class_weight", {"0": 1.0}, {"0": 1.0}),
        ("class_weight", [{"0": 1.0}], [{"0": 1.0}]),
    ],
)
def test_ml_config_accepts_every_estimator_value_form(
    key, value, expected, tmp_path
):
    """Every form scikit-learn accepts survives the normalizer unchanged."""
    assert _load(tmp_path, model={key: value})[key] == expected


@pytest.mark.parametrize(
    "key, value, expected_message",
    [
        ("max_features", "auto", "removed"),
        ("max_features", "None", "omit the key or write null"),
        ("max_features", "sqrt2", "'sqrt', 'log2'"),
        ("max_features", 0, "integer >= 1"),
        ("max_features", 1.5, "(0.0, 1.0]"),
        ("max_features", True, "Valid values are"),
        ("class_weight", "nope", "'balanced', 'balanced_subsample'"),
        ("class_weight", 5, "mapping of class label to weight"),
    ],
)
def test_ml_config_rejects_invalid_estimator_values(
    key, value, expected_message, tmp_path
):
    """An invalid value is caught in config, naming the value and what is valid.

    The `"None"` case is the YAML trap: the bare word parses as a string,
    so its rejection has to say how null is actually written.
    """
    with pytest.raises(ConfigError) as excinfo:
        _load(tmp_path, model={key: value})

    message = str(excinfo.value)
    assert f"model.{key}" in message
    assert repr(value) in message
    assert expected_message in message


def test_seed_is_a_top_level_key_with_a_cli_flag(tmp_path):
    """`seed` resolves from the config file or `--seed`, and only from the top level."""
    assert _load(tmp_path, seed=7)["seed"] == 7

    parser = build_trainer_argument_parser("test")
    args = parser.parse_args(
        ["-i", _REQUIRED["input_path"], "-o", _REQUIRED["output_dir"], "--seed", "7"]
    )
    assert resolve_trainer_runtime_args(args).seed == 7

    # Beside a config file the flag replaces the file's seed and nothing else.
    config_path = tmp_path / "ml_config.json"
    config_path.write_text(json.dumps({**_REQUIRED, "seed": 7, "test_size": 0.3}))
    args = parser.parse_args(["-c", str(config_path), "--seed", "9"])
    resolved = resolve_trainer_runtime_args(args)
    assert resolved.seed == 9
    assert resolved.test_size == 0.3

    with pytest.raises(ConfigError, match="top level"):
        _load(tmp_path, model={"seed": 7})


def test_feature_importance_method_resolves_under_evaluation(tmp_path):
    """The importance estimator is an ``evaluation`` key with three choices.

    Null stays null in the resolved config: each trainer reads it as its own
    estimator, impurity for the forest and permutation for the neighbours
    classifier, which measures none.
    """
    resolved = _load(
        tmp_path,
        evaluation={
            "feature_importance_method": "grouped_permutation",
            "feature_importance_correlation_threshold": 0.9,
        },
    )
    assert resolved["feature_importance_method"] == "grouped_permutation"
    assert resolved["feature_importance_correlation_threshold"] == 0.9
    assert _load(tmp_path)["feature_importance_method"] is None

    with pytest.raises(ConfigError, match="mdi, permutation, grouped_permutation"):
        _load(tmp_path, evaluation={"feature_importance_method": "shapley"})
    with pytest.raises(ConfigError, match="fraction between 0 and 1"):
        _load(tmp_path, evaluation={"feature_importance_correlation_threshold": 1.4})
    with pytest.raises(ConfigError, match="feature_importance_method"):
        _load(tmp_path, feature_importance_method="mdi")
