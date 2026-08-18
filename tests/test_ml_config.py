import json

import pytest

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
