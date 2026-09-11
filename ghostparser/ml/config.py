"""Configuration for Ghostparser ML.

Owns the ML defaults, choices, and validation rules; shared helpers come from
the :mod:`ghostparser.config` trunk.
"""

import argparse

from ..cli_config import resolve_cli_or_config_args
from ..config import (
    DEFAULT_OVERWRITE,
    ConfigError,
    _load_raw_config,
    _resolve_path,
    _validate_overwrite_flag,
    _validate_required_path,
)

DEFAULT_ML_OUTPUT_DIR = "ml_results"
DEFAULT_TARGET_COLUMN = "class"
DEFAULT_TEST_SIZE = 0.2
DEFAULT_CV_FOLDS = 5
DEFAULT_RANDOM_STATE = None
DEFAULT_RARE_CLASS_POLICY = "warn_reduce_cv"
DEFAULT_EVALUATION_METRICS = "all"
DEFAULT_REPORT_CLASS_DISTRIBUTION = True
DEFAULT_REPORT_CONFUSION_MATRIX = True
DEFAULT_REPORT_FEATURE_IMPORTANCE = True
DEFAULT_SAVE_LABEL_MAP = True
DEFAULT_SAVE_PREDICTIONS = True
DEFAULT_N_ESTIMATORS = 200
DEFAULT_MAX_DEPTH = None
DEFAULT_MIN_SAMPLES_SPLIT = 2
DEFAULT_MIN_SAMPLES_LEAF = 1
DEFAULT_MAX_FEATURES = None
DEFAULT_CLASS_WEIGHT = None
DEFAULT_N_JOBS = -1
DEFAULT_N_NEIGHBORS = 5
DEFAULT_KNN_WEIGHTS = "uniform"
DEFAULT_KNN_ALGORITHM = "auto"
DEFAULT_KNN_LEAF_SIZE = 30
DEFAULT_KNN_METRIC = "minkowski"
DEFAULT_KNN_P = 2

RARE_CLASS_POLICY_CHOICES = ("warn_reduce_cv", "warn_skip_cv", "error")
KNN_WEIGHT_CHOICES = ("uniform", "distance")
KNN_ALGORITHM_CHOICES = ("auto", "ball_tree", "kd_tree", "brute")
EVALUATION_METRICS_CHOICES = ("all", "primary", "diagnostic", "per_bit")
MAX_FEATURES_STRING_CHOICES = ("sqrt", "log2")
CLASS_WEIGHT_STRING_CHOICES = ("balanced", "balanced_subsample")

# `auto` was scikit-learn's max_features default until 1.1 and was removed in
# 1.3 for being ambiguous between classifiers and regressors. It is worth its
# own message because it is the value most users reach for first.
_REMOVED_MAX_FEATURES_SPELLINGS = frozenset({"auto"})

# Both keys default to null and accept it, but YAML reads the bare words `None`
# and `none` as strings, so a rejection has to say how null is actually written.
_HOW_TO_WRITE_NULL = (
    "To use null, omit the key or write null (YAML also accepts ~); the bare "
    "words None and none are read as strings and rejected."
)

def _validate_optional_path(payload: dict, key: str, default: str) -> str:
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(
            f"Config field {key} must be a non-empty string when provided"
        )
    return _resolve_path(value.strip())


def _validate_optional_string(
    payload: dict, key: str, default: str | None
) -> str | None:
    value = payload.get(key, default)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(
            f"Config field {key} must be a non-empty string when provided"
        )
    return value.strip()


def _validate_optional_float(payload: dict, key: str, default: float) -> float:
    value = payload.get(key, default)
    if value is None:
        value = default
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"Config field {key} must be numeric") from exc
    if not 0 < value < 1:
        raise ConfigError(f"Config field {key} must be a fraction between 0 and 1")
    return value


def _validate_positive_int(payload: dict, key: str, default: int | None) -> int:
    """Require an integer >= 1, taking ``default`` when the key is absent."""
    value = payload.get(key, default)
    if not isinstance(value, int) or value < 1:
        raise ConfigError(f"Config field {key} must be an integer >= 1")
    return value


def _validate_optional_positive_int(
    payload: dict, key: str, default: int | None
) -> int | None:
    """Allow ``None``, else require an integer >= 1."""
    value = payload.get(key, default)
    if value is None:
        return None
    if not isinstance(value, int) or value < 1:
        raise ConfigError(f"Config field {key} must be an integer >= 1")
    return value


def _validate_optional_int(payload: dict, key: str, default: int | None) -> int | None:
    value = payload.get(key, default)
    if value is None:
        return None
    if not isinstance(value, int):
        raise ConfigError(f"Config field {key} must be an integer when provided")
    return value


def _validate_optional_choice(
    payload: dict, key: str, default: str, choices: tuple[str, ...]
) -> str:
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, str) or value not in choices:
        raise ConfigError(f"Config field {key} must be one of: {', '.join(choices)}")
    return value


def _validate_optional_bool(payload: dict, key: str, default: bool) -> bool:
    value = payload.get(key, default)
    if not isinstance(value, bool):
        raise ConfigError(f"Config field {key} must be a boolean")
    return value


def normalize_max_features(value: object, *, where: str) -> object:
    """Validate one ``max_features`` value against what scikit-learn accepts.

    Args:
        value: The configured value, straight from the config file.
        where: Dotted config location used in the error message.

    Returns:
        ``None``, one of :data:`MAX_FEATURES_STRING_CHOICES`, an ``int >= 1``,
        or a ``float`` in ``(0.0, 1.0]``.

    Raises:
        ConfigError: If the value is anything else, naming what was received,
            every accepted form, and how to write null -- the strings ``"None"``
            and ``"none"`` that YAML produces from the bare words land here.
    """
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped in MAX_FEATURES_STRING_CHOICES:
            return stripped
        if stripped.lower() in _REMOVED_MAX_FEATURES_SPELLINGS:
            raise ConfigError(
                f"Config field {where} got {value!r}, which scikit-learn removed "
                "in 1.3. Use 'sqrt' for the same behaviour on a classifier, or "
                "null to use every feature at each split."
            )
        raise ConfigError(_max_features_message(where, value))
    # bool is a subclass of int, so it has to be rejected before the int check.
    if isinstance(value, bool):
        raise ConfigError(_max_features_message(where, value))
    if isinstance(value, int):
        if value < 1:
            raise ConfigError(_max_features_message(where, value))
        return value
    if isinstance(value, float):
        if not 0.0 < value <= 1.0:
            raise ConfigError(_max_features_message(where, value))
        return value
    raise ConfigError(_max_features_message(where, value))


def _max_features_message(where: str, value: object) -> str:
    choices = ", ".join(repr(choice) for choice in MAX_FEATURES_STRING_CHOICES)
    return (
        f"Config field {where} got {value!r}. Valid values are {choices}, an "
        "integer >= 1 (that many features per split), a float in (0.0, 1.0] "
        "(that fraction of the features), or null to use every feature at each "
        f"split. {_HOW_TO_WRITE_NULL}"
    )


def normalize_class_weight(value: object, *, where: str) -> object:
    """Validate one ``class_weight`` value against what scikit-learn accepts.

    Args:
        value: The configured value, straight from the config file.
        where: Dotted config location used in the error message.

    Returns:
        ``None``, one of :data:`CLASS_WEIGHT_STRING_CHOICES`, a mapping of
        class label to weight, or a list of such mappings (one per label).

    Raises:
        ConfigError: If the value is anything else, naming what was received,
            every accepted form, and how to write null -- the strings ``"None"``
            and ``"none"`` that YAML produces from the bare words land here.
    """
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped in CLASS_WEIGHT_STRING_CHOICES:
            return stripped
        raise ConfigError(_class_weight_message(where, value))
    if isinstance(value, dict):
        return value
    if isinstance(value, list) and all(isinstance(item, dict) for item in value):
        return value
    raise ConfigError(_class_weight_message(where, value))


def _class_weight_message(where: str, value: object) -> str:
    choices = ", ".join(repr(choice) for choice in CLASS_WEIGHT_STRING_CHOICES)
    return (
        f"Config field {where} got {value!r}. Valid values are {choices}, a "
        "mapping of class label to weight, a list of such mappings (one per "
        f"label), or null for no class weighting. {_HOW_TO_WRITE_NULL}"
    )


def normalize_ml_payload(payload: dict) -> dict:
    """Normalize ML config payload using a strict top-level + nested layout.

    Required top-level keys:
      - `input_path`, `output_dir`

    Top-level runtime keys:
      - `target_column`, `test_size`, `cv_folds`, `rare_class_policy`,
        `random_state`, `n_jobs`

    Model hyperparameters must be provided under `model`.
    Evaluation reporting controls must be provided under `evaluation`.
    """
    input_path = _validate_required_path(payload, "input_path")
    output_dir = _validate_required_path(payload, "output_dir")

    target_column = _validate_optional_string(
        payload, "target_column", DEFAULT_TARGET_COLUMN
    )
    if target_column is None:
        raise ConfigError("Config field target_column must be a non-empty string")

    model_section = payload.get("model", {})
    if model_section is not None and not isinstance(model_section, dict):
        raise ConfigError("Config field 'model' must be a mapping/object when provided")
    evaluation_section = payload.get("evaluation", {})
    if evaluation_section is not None and not isinstance(evaluation_section, dict):
        raise ConfigError(
            "Config field 'evaluation' must be a mapping/object when provided"
        )

    forbidden_top_level_keys = {
        "n_estimators",
        "max_depth",
        "min_samples_split",
        "min_samples_leaf",
        "max_features",
        "class_weight",
        "n_neighbors",
        "weights",
        "algorithm",
        "leaf_size",
        "metric",
        "p",
        "metrics",
        "report_class_distribution",
        "report_confusion_matrix",
        "report_feature_importance",
        "save_label_map",
        "save_predictions",
    }
    present_forbidden = forbidden_top_level_keys & set(payload.keys())
    if present_forbidden:
        raise ConfigError(
            "Do not place model hyperparameters or evaluation reporting controls at the top level. "
            "Group model parameters under 'model' and reporting controls under 'evaluation'. Offending keys: "
            f"{', '.join(sorted(present_forbidden))}"
        )

    forbidden_model_runtime_keys = {"random_state", "n_jobs"}
    present_forbidden_model_runtime = forbidden_model_runtime_keys & set(
        model_section.keys()
    )
    if present_forbidden_model_runtime:
        raise ConfigError(
            "Place 'random_state' and 'n_jobs' at the top level, not under 'model'. Offending keys: "
            f"{', '.join(sorted(present_forbidden_model_runtime))}"
        )

    return {
        "input_path": input_path,
        "output_dir": output_dir,
        "overwrite": _validate_overwrite_flag(payload, DEFAULT_OVERWRITE),
        "target_column": target_column,
        "test_size": _validate_optional_float(payload, "test_size", DEFAULT_TEST_SIZE),
        "cv_folds": _validate_positive_int(
            payload, "cv_folds", DEFAULT_CV_FOLDS
        ),
        "rare_class_policy": _validate_optional_choice(
            payload,
            "rare_class_policy",
            DEFAULT_RARE_CLASS_POLICY,
            RARE_CLASS_POLICY_CHOICES,
        ),
        "random_state": _validate_optional_int(
            payload, "random_state", DEFAULT_RANDOM_STATE
        ),
        "n_jobs": _validate_optional_int(payload, "n_jobs", DEFAULT_N_JOBS),
        "n_estimators": _validate_positive_int(
            model_section, "n_estimators", DEFAULT_N_ESTIMATORS
        ),
        "max_depth": _validate_optional_int(
            model_section, "max_depth", DEFAULT_MAX_DEPTH
        ),
        "min_samples_split": _validate_positive_int(
            model_section, "min_samples_split", DEFAULT_MIN_SAMPLES_SPLIT
        ),
        "min_samples_leaf": _validate_positive_int(
            model_section, "min_samples_leaf", DEFAULT_MIN_SAMPLES_LEAF
        ),
        "max_features": normalize_max_features(
            model_section.get("max_features", DEFAULT_MAX_FEATURES),
            where="model.max_features",
        ),
        "class_weight": normalize_class_weight(
            model_section.get("class_weight", DEFAULT_CLASS_WEIGHT),
            where="model.class_weight",
        ),
        "n_neighbors": _validate_positive_int(
            model_section, "n_neighbors", DEFAULT_N_NEIGHBORS
        ),
        "weights": _validate_optional_choice(
            model_section, "weights", DEFAULT_KNN_WEIGHTS, KNN_WEIGHT_CHOICES
        ),
        "algorithm": _validate_optional_choice(
            model_section, "algorithm", DEFAULT_KNN_ALGORITHM, KNN_ALGORITHM_CHOICES
        ),
        "leaf_size": _validate_positive_int(
            model_section, "leaf_size", DEFAULT_KNN_LEAF_SIZE
        ),
        "metric": _validate_optional_string(
            model_section, "metric", DEFAULT_KNN_METRIC
        ),
        "p": _validate_positive_int(model_section, "p", DEFAULT_KNN_P),
        "evaluation_metrics": _validate_optional_choice(
            evaluation_section,
            "metrics",
            DEFAULT_EVALUATION_METRICS,
            EVALUATION_METRICS_CHOICES,
        ),
        "report_class_distribution": _validate_optional_bool(
            evaluation_section,
            "report_class_distribution",
            DEFAULT_REPORT_CLASS_DISTRIBUTION,
        ),
        "report_confusion_matrix": _validate_optional_bool(
            evaluation_section,
            "report_confusion_matrix",
            DEFAULT_REPORT_CONFUSION_MATRIX,
        ),
        "report_feature_importance": _validate_optional_bool(
            evaluation_section,
            "report_feature_importance",
            DEFAULT_REPORT_FEATURE_IMPORTANCE,
        ),
        "save_label_map": _validate_optional_bool(
            evaluation_section,
            "save_label_map",
            DEFAULT_SAVE_LABEL_MAP,
        ),
        "save_predictions": _validate_optional_bool(
            evaluation_section,
            "save_predictions",
            DEFAULT_SAVE_PREDICTIONS,
        ),
    }


def load_ml_config(config_file: str) -> dict:
    """Load and normalize an ML config file."""
    return normalize_ml_payload(_load_raw_config(config_file))


def build_trainer_argument_parser(description: str) -> argparse.ArgumentParser:
    """Build the argument parser both trainers expose.

    Args:
        description: The trainer's own ``--help`` description.

    Returns:
        A parser carrying the shared config/input/output/overwrite flags.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-c", "--config-file", type=str, default=None,
        help="Path to a JSON or YAML config file",
    )
    parser.add_argument(
        "-i", "--input-path", type=str, default=None,
        help="Path to summary_statistics.tsv",
    )
    parser.add_argument(
        "-o", "--output-dir", type=str, default=None,
        help="Directory for ML outputs",
    )
    parser.add_argument(
        "--no-overwrite", dest="no_overwrite", action="store_true", default=None,
        help="Append a numeric suffix when the output directory already exists",
    )
    return parser


def resolve_trainer_runtime_args(args: argparse.Namespace) -> argparse.Namespace:
    """Resolve a trainer's CLI/config arguments under config-file-wins precedence.

    Args:
        args: The parsed CLI namespace.

    Returns:
        The resolved namespace.
    """
    return resolve_cli_or_config_args(
        args,
        load_config=load_ml_config,
        normalize_payload=normalize_ml_payload,
        payload_arg_names=["input_path", "output_dir", "no_overwrite"],
    )
