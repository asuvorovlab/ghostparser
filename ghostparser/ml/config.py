"""Configuration helpers for Ghostparser ML workflows."""

from __future__ import annotations

import json
from pathlib import Path

from ..config import DEFAULT_OVERWRITE, ConfigError, _validate_overwrite_flag

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
DEFAULT_MAX_FEATURES = "sqrt"
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


def _resolve_path(path_str: str) -> str:
    return str(Path(path_str).expanduser().resolve())


def _load_raw_config(config_file: str) -> dict:
    path = Path(config_file)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    elif suffix in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ConfigError("YAML support requires PyYAML to be installed") from exc

        with open(path, "r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle)
    else:
        raise ConfigError("Config file must be .json, .yaml, or .yml")

    if not isinstance(payload, dict):
        raise ConfigError("Config root must be a key/value object")
    return payload


def _validate_required_path(payload: dict, key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return _resolve_path(value.strip())


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


def _validate_optional_positive_int(
    payload: dict, key: str, default: int | None
) -> int | None:
    value = payload.get(key, default)
    if value is None:
        raise ConfigError(f"Config field {key} must be an integer >= 1")
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
        "cv_folds": _validate_optional_positive_int(
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
        "n_estimators": _validate_optional_positive_int(
            model_section, "n_estimators", DEFAULT_N_ESTIMATORS
        ),
        "max_depth": _validate_optional_int(
            model_section, "max_depth", DEFAULT_MAX_DEPTH
        ),
        "min_samples_split": _validate_optional_positive_int(
            model_section, "min_samples_split", DEFAULT_MIN_SAMPLES_SPLIT
        ),
        "min_samples_leaf": _validate_optional_positive_int(
            model_section, "min_samples_leaf", DEFAULT_MIN_SAMPLES_LEAF
        ),
        "max_features": model_section.get("max_features", DEFAULT_MAX_FEATURES),
        "class_weight": model_section.get("class_weight", DEFAULT_CLASS_WEIGHT),
        "n_neighbors": _validate_optional_positive_int(
            model_section, "n_neighbors", DEFAULT_N_NEIGHBORS
        ),
        "weights": _validate_optional_choice(
            model_section, "weights", DEFAULT_KNN_WEIGHTS, KNN_WEIGHT_CHOICES
        ),
        "algorithm": _validate_optional_choice(
            model_section, "algorithm", DEFAULT_KNN_ALGORITHM, KNN_ALGORITHM_CHOICES
        ),
        "leaf_size": _validate_optional_positive_int(
            model_section, "leaf_size", DEFAULT_KNN_LEAF_SIZE
        ),
        "metric": _validate_optional_string(
            model_section, "metric", DEFAULT_KNN_METRIC
        ),
        "p": _validate_optional_positive_int(model_section, "p", DEFAULT_KNN_P),
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
