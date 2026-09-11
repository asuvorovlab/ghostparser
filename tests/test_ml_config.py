import json

import pytest

from ghostparser.config import ConfigError
from ghostparser.ml.config import load_ml_config

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
    "key, expected",
    [
        ("target_column", "class"),
        ("overwrite", True),
        ("min_samples_split", 2),
        ("min_samples_leaf", 1),
        ("max_features", None),
    ],
)
def test_ml_config_fills_defaults_for_omitted_keys(key, expected, tmp_path):
    """A config carrying only the required paths takes every other default."""
    assert _load(tmp_path)[key] == expected


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

    The `"None"` case is the YAML trap -- the bare word parses as a string --
    so its rejection has to say how null is actually written.
    """
    with pytest.raises(ConfigError) as excinfo:
        _load(tmp_path, model={key: value})

    message = str(excinfo.value)
    assert f"model.{key}" in message
    assert repr(value) in message
    assert expected_message in message
