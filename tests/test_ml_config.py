"""Config resolution for the ML trainers.

Marked ``config`` throughout. Individual defaults are not pinned: CONFIG.md
and the shipped sample configs state them; what is pinned is that the loader
resolves, that explicit values win, that the two estimator passthrough keys
accept every form scikit-learn does and nothing else, and that a flag beside
a config file overrides it.
"""

import json
from pathlib import Path

import pytest
import yaml

from ghostparser.config import ConfigError
from ghostparser.ml import hyper_tune
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


def test_trainer_sample_config_names_every_key_at_its_default(tmp_path):
    """`ml_trainer.yaml` names every trainer key, each at its default.

    The raw file must name every resolved key (the `evaluation.metrics` key
    resolves as `evaluation_metrics`), and loading it must equal loading a
    required-keys-only file with the same paths, so the value the sample shows
    for each key is the value the code would have used anyway.
    """
    sample_path = _SAMPLE_DIR / "ml_trainer.yaml"
    raw = yaml.safe_load(sample_path.read_text())
    documented = set(raw) | set(raw["model"]) | set(raw["evaluation"])
    documented.add("evaluation_metrics")

    resolved = load_ml_config(str(sample_path))
    assert [key for key in resolved if key not in documented] == []

    paths = {key: raw[key] for key in ("input_path", "output_dir")}
    config_path = tmp_path / "required.json"
    config_path.write_text(json.dumps(paths))
    assert resolved == load_ml_config(str(config_path))


@pytest.mark.parametrize(
    "payload, key, expected",
    [
        ({"overwrite": False}, "overwrite", False),
        ({"target_column": "label"}, "target_column", "label"),
        ({"seed": 7}, "seed", 7),
        # Hyperparameters under `model` surface at the top level.
        ({"model": {"n_estimators": 25}}, "n_estimators", 25),
        ({"model": {"min_samples_leaf": 4}}, "min_samples_leaf", 4),
        # Every form scikit-learn accepts for the two passthrough keys.
        ({"model": {"max_features": None}}, "max_features", None),
        ({"model": {"max_features": "log2"}}, "max_features", "log2"),
        ({"model": {"max_features": 3}}, "max_features", 3),
        ({"model": {"max_features": 0.5}}, "max_features", 0.5),
        ({"model": {"class_weight": "balanced"}}, "class_weight", "balanced"),
        ({"model": {"class_weight": {"0": 1.0}}}, "class_weight", {"0": 1.0}),
        ({"model": {"class_weight": [{"0": 1.0}]}}, "class_weight", [{"0": 1.0}]),
        # Reporting controls under `evaluation`.
        (
            {"evaluation": {"feature_importance_method": "grouped_permutation"}},
            "feature_importance_method",
            "grouped_permutation",
        ),
        (
            {"evaluation": {"feature_importance_correlation_threshold": 0.9}},
            "feature_importance_correlation_threshold",
            0.9,
        ),
    ],
)
def test_ml_config_accepts_explicit_values(payload, key, expected, tmp_path):
    """An explicit value survives the normalizer, from the top level, `model` or `evaluation`."""
    assert _load(tmp_path, **payload)[key] == expected


@pytest.mark.parametrize(
    "payload, fragments",
    [
        ({"model": {"max_features": "auto"}}, ("model.max_features", "'auto'", "removed")),
        # The YAML trap: the bare word parses as a string, so the message says
        # how null is actually written.
        ({"model": {"max_features": "None"}}, ("model.max_features", "omit the key or write null")),
        ({"model": {"max_features": "sqrt2"}}, ("model.max_features", "'sqrt', 'log2'")),
        ({"model": {"max_features": 0}}, ("model.max_features", "integer >= 1")),
        ({"model": {"max_features": 1.5}}, ("model.max_features", "(0.0, 1.0]")),
        ({"model": {"max_features": True}}, ("model.max_features", "Valid values are")),
        ({"model": {"class_weight": "nope"}}, ("model.class_weight", "'balanced', 'balanced_subsample'")),
        ({"model": {"class_weight": 5}}, ("model.class_weight", "mapping of class label to weight")),
        # Stratified k-fold cannot split into one fold, so 1 would otherwise
        # pass the loader and fail inside scikit-learn as an internal error.
        *(({"cv_folds": value}, ("cv_folds must be an integer >= 2",)) for value in (1, 0, None, 2.5, True)),
        ({"evaluation": {"feature_importance_method": "shapley"}}, ("mdi, permutation, grouped_permutation",)),
        ({"evaluation": {"feature_importance_correlation_threshold": 1.4}}, ("fraction between 0 and 1",)),
        # A key in the wrong section is refused, not ignored.
        ({"feature_importance_method": "mdi"}, ("feature_importance_method",)),
        ({"model": {"seed": 7}}, ("top level",)),
    ],
)
def test_ml_config_rejects_invalid_values(payload, fragments, tmp_path):
    """An invalid value or a misplaced key fails in config, naming the key and what is valid."""
    with pytest.raises(ConfigError) as excinfo:
        _load(tmp_path, **payload)
    for fragment in fragments:
        assert fragment in str(excinfo.value)


@pytest.mark.parametrize("entry_point", ["trainer", "tuner"])
def test_flags_beside_a_config_file_override_it(entry_point, tmp_path):
    """`--seed` and `--no-overwrite` beside a config file replace its values, and only those."""
    payload = {**_REQUIRED, "seed": 7, "overwrite": True, "test_size": 0.3}
    if entry_point == "trainer":
        parser = build_trainer_argument_parser("test")
        resolve = resolve_trainer_runtime_args
    else:
        payload["hyperparameter_tuning"] = {"search_space": {"n_estimators": [5]}}
        parser = hyper_tune._build_argument_parser()
        resolve = hyper_tune.resolve_tuner_args
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(payload))

    resolved = resolve(parser.parse_args(["-c", str(config_path), "--seed", "9", "--no-overwrite"]))

    assert resolved.seed == 9
    assert resolved.overwrite is False
    assert resolved.test_size == 0.3
