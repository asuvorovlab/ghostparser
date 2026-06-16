"""Hyperparameter tuning for Ghostparser ML baselines.

This module runs a config-driven grid or random search over the supported ML
trainers. Put the tuning settings under the top-level
`hyperparameter_tuning` section and the module writes the best fitted model
plus candidate rankings.
"""

from __future__ import annotations

import argparse
import itertools
import pickle
import random
import time
from pathlib import Path

import numpy as np

from ..config import ConfigError, prepare_output_directory
from . import ml_utils as shared
from . import multi_knn as knn_module
from . import random_forest as rf_module
from .config import (
    DEFAULT_CLASS_WEIGHT,
    DEFAULT_KNN_ALGORITHM,
    DEFAULT_KNN_LEAF_SIZE,
    DEFAULT_KNN_METRIC,
    DEFAULT_KNN_P,
    DEFAULT_KNN_WEIGHTS,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FEATURES,
    DEFAULT_MIN_SAMPLES_LEAF,
    DEFAULT_MIN_SAMPLES_SPLIT,
    DEFAULT_N_ESTIMATORS,
    DEFAULT_N_JOBS,
    DEFAULT_N_NEIGHBORS,
    DEFAULT_RANDOM_STATE,
    DEFAULT_TARGET_COLUMN,
    _load_raw_config,
)

DEFAULT_TUNER_MODEL = "random_forest"
DEFAULT_SEARCH_METHOD = "grid"
DEFAULT_OBJECTIVE_METRIC = "exact_match_accuracy"
DEFAULT_TOP_K = 10
DEFAULT_RANDOM_ITERATIONS = 20
DEFAULT_MAX_CANDIDATES = 5000

SUPPORTED_MODELS = ("multi_knn", "random_forest")
SUPPORTED_SEARCH_METHODS = ("grid", "random")
SUPPORTED_OBJECTIVES = (
    "exact_match_accuracy",
    "hamming_loss",
    "bitwise_accuracy",
    "micro_f1",
    "macro_f1",
    "weighted_f1",
)
MINIMIZE_OBJECTIVES = {"hamming_loss"}

MODEL_SEARCH_KEYS = {
    "random_forest": {
        "n_estimators",
        "max_depth",
        "min_samples_split",
        "min_samples_leaf",
        "max_features",
        "class_weight",
    },
    "multi_knn": {
        "n_neighbors",
        "weights",
        "algorithm",
        "leaf_size",
        "metric",
        "p",
    },
}

MODEL_DEFAULTS = {
    "random_forest": {
        "n_estimators": DEFAULT_N_ESTIMATORS,
        "max_depth": DEFAULT_MAX_DEPTH,
        "min_samples_split": DEFAULT_MIN_SAMPLES_SPLIT,
        "min_samples_leaf": DEFAULT_MIN_SAMPLES_LEAF,
        "max_features": DEFAULT_MAX_FEATURES,
        "class_weight": DEFAULT_CLASS_WEIGHT,
    },
    "multi_knn": {
        "n_neighbors": DEFAULT_N_NEIGHBORS,
        "weights": DEFAULT_KNN_WEIGHTS,
        "algorithm": DEFAULT_KNN_ALGORITHM,
        "leaf_size": DEFAULT_KNN_LEAF_SIZE,
        "metric": DEFAULT_KNN_METRIC,
        "p": DEFAULT_KNN_P,
    },
}


def _format_seconds(seconds: float) -> str:
    return f"{seconds:.2f}s"


def _pluralize(word: str, count: int) -> str:
    return word if count == 1 else f"{word}s"


def _log_progress(message: str) -> None:
    print(f"[hyper_tune] {message}", flush=True)


RUNTIME_KEYS = {
    "input_path",
    "output_dir",
    "target_column",
    "test_size",
    "cv_folds",
    "rare_class_policy",
    "random_state",
    "n_jobs",
}


def _validate_required_path(payload: dict, key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return str(Path(value.strip()).expanduser().resolve())


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


def _normalize_search_values(value: object) -> list[object]:
    if isinstance(value, (list, tuple)):
        candidates = list(value)
        if not candidates:
            raise ConfigError("Search space lists must contain at least one value")
        return candidates
    return [value]


def _build_candidate_grid(
    search_space: dict[str, object],
) -> tuple[list[dict[str, object]], list[str]]:
    if not search_space:
        raise ConfigError(
            "Config field search_space must define at least one parameter"
        )

    keys = list(search_space.keys())
    values = [_normalize_search_values(search_space[key]) for key in keys]
    candidates = [
        dict(zip(keys, combo, strict=True)) for combo in itertools.product(*values)
    ]
    return candidates, keys


def _objective_direction(metric_name: str) -> str:
    return "min" if metric_name in MINIMIZE_OBJECTIVES else "max"


def _candidate_score(cv_results: dict, objective_metric: str) -> float:
    aggregate = cv_results.get("aggregate", {})
    key = f"{objective_metric}_mean"
    if key not in aggregate:
        raise ConfigError(
            f"Cross-validation results do not include objective metric: {objective_metric}"
        )
    return float(aggregate[key])


def _build_training_namespace(
    base_config: dict[str, object],
    candidate_params: dict[str, object],
    model_name: str,
) -> argparse.Namespace:
    payload = dict(base_config)
    payload.update(MODEL_DEFAULTS[model_name])
    payload.update(candidate_params)
    return argparse.Namespace(**payload)


def load_hyper_tune_config(config_file: str) -> dict[str, object]:
    payload = _load_raw_config(config_file)
    return normalize_hyper_tune_payload(payload)


def normalize_hyper_tune_payload(payload: dict) -> dict[str, object]:
    input_path = _validate_required_path(payload, "input_path")
    output_dir = _validate_required_path(payload, "output_dir")

    allowed_top_level_keys = RUNTIME_KEYS | {"hyperparameter_tuning"}
    unexpected_top_level_keys = set(payload) - allowed_top_level_keys
    if unexpected_top_level_keys:
        raise ConfigError(
            "Unexpected top-level keys in hyperparameter_tuning config: "
            f"{', '.join(sorted(unexpected_top_level_keys))}"
        )

    target_column = _validate_optional_string(
        payload, "target_column", DEFAULT_TARGET_COLUMN
    )
    if target_column is None:
        raise ConfigError("Config field target_column must be a non-empty string")

    if "hyperparameter_tuning" not in payload:
        raise ConfigError("Missing required config field: hyperparameter_tuning")

    tuning_section = payload.get("hyperparameter_tuning")
    if tuning_section is not None and not isinstance(tuning_section, dict):
        raise ConfigError(
            "Config field 'hyperparameter_tuning' must be a mapping/object when provided"
        )
    tuning_section = tuning_section or {}
    forbidden_top_level_keys = {
        "evaluation",
        "model",
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
    }
    present_forbidden_top_level_keys = forbidden_top_level_keys & set(payload)
    if present_forbidden_top_level_keys:
        raise ConfigError(
            "Do not place model sections, evaluation controls, or model hyperparameters at the top level in hyperparameter_tuning configs. "
            "Use only the runtime keys plus 'hyperparameter_tuning'. Offending keys: "
            f"{', '.join(sorted(present_forbidden_top_level_keys))}"
        )

    search_space = tuning_section.get("search_space", {})
    if search_space is not None and not isinstance(search_space, dict):
        raise ConfigError(
            "Config field 'hyperparameter_tuning.search_space' must be a mapping/object when provided"
        )

    model_name = _validate_optional_choice(
        tuning_section,
        "model",
        DEFAULT_TUNER_MODEL,
        SUPPORTED_MODELS,
    )
    search_method = _validate_optional_choice(
        tuning_section,
        "method",
        DEFAULT_SEARCH_METHOD,
        SUPPORTED_SEARCH_METHODS,
    )
    objective_metric = _validate_optional_choice(
        tuning_section,
        "objective",
        DEFAULT_OBJECTIVE_METRIC,
        SUPPORTED_OBJECTIVES,
    )
    top_k = _validate_optional_positive_int(tuning_section, "top_k", DEFAULT_TOP_K)
    n_iter = _validate_optional_positive_int(
        tuning_section, "n_iter", DEFAULT_RANDOM_ITERATIONS
    )
    max_candidates = _validate_optional_positive_int(
        tuning_section, "max_candidates", DEFAULT_MAX_CANDIDATES
    )

    forbidden_search_keys = RUNTIME_KEYS | {"hyperparameter_tuning"}
    present_forbidden_search_keys = forbidden_search_keys & set(search_space)
    if present_forbidden_search_keys:
        raise ConfigError(
            "Do not place runtime fields inside hyperparameter_tuning.search_space. Offending keys: "
            f"{', '.join(sorted(present_forbidden_search_keys))}"
        )

    allowed_search_keys = MODEL_SEARCH_KEYS[model_name]
    unexpected_search_keys = set(search_space) - allowed_search_keys
    if unexpected_search_keys:
        raise ConfigError(
            f"Unsupported search_space keys for {model_name}: "
            f"{', '.join(sorted(unexpected_search_keys))}"
        )

    normalized_config = {
        "input_path": input_path,
        "output_dir": output_dir,
        "target_column": target_column,
        "test_size": _validate_optional_float(payload, "test_size", 0.2),
        "cv_folds": _validate_optional_positive_int(payload, "cv_folds", 5),
        "rare_class_policy": _validate_optional_choice(
            payload,
            "rare_class_policy",
            "warn_reduce_cv",
            ("warn_reduce_cv", "warn_skip_cv", "error"),
        ),
        "random_state": _validate_optional_int(
            payload, "random_state", DEFAULT_RANDOM_STATE
        ),
        "n_jobs": _validate_optional_int(payload, "n_jobs", DEFAULT_N_JOBS),
        "model_name": model_name,
        "search_method": search_method,
        "objective_metric": objective_metric,
        "objective_direction": _objective_direction(objective_metric),
        "top_k": top_k,
        "n_iter": n_iter,
        "max_candidates": max_candidates,
        "search_space": search_space,
    }

    normalized_config.update(MODEL_DEFAULTS[model_name])
    return normalized_config


def _evaluate_candidate(
    base_config: dict[str, object],
    candidate_params: dict[str, object],
    model_name: str,
    x_train: np.ndarray,
    y_train: np.ndarray,
    labels_train: np.ndarray,
    folds: int,
) -> tuple[argparse.Namespace, dict, float]:
    candidate_config = _build_training_namespace(
        base_config, candidate_params, model_name
    )
    if model_name == "random_forest":

        def model_factory() -> object:
            return rf_module._build_model(candidate_config)

        cv_results = rf_module._cross_validate(
            model_factory,
            x_train,
            y_train,
            labels_train,
            folds,
            candidate_config.random_state,
        )
    elif model_name == "multi_knn":
        cv_results = knn_module._cross_validate(
            x_train,
            y_train,
            labels_train,
            folds,
            candidate_config.random_state,
            candidate_config,
        )
    else:
        raise ConfigError(f"Unsupported tuning model: {model_name}")
    score = _candidate_score(cv_results, candidate_config.objective_metric)
    return candidate_config, cv_results, score


def tune_hyperparameters(config: argparse.Namespace) -> dict[str, object]:
    run_start = time.perf_counter()
    _log_progress(
        f"Starting hyperparameter tuning for {config.model_name} with {config.search_method} search"
    )
    load_start = time.perf_counter()
    base_config = dict(vars(config))
    base_config.update(MODEL_DEFAULTS[config.model_name])
    rows = shared.read_tsv_rows(config.input_path)
    matrix = shared.rows_to_matrix(rows, config.target_column)
    labels = shared.combination_labels(matrix.train_labels)
    load_seconds = time.perf_counter() - load_start

    split_start = time.perf_counter()
    x_train, x_test, y_train, y_test, labels_train, labels_test, split_notes = (
        shared.split_dataset(
            matrix.train_features,
            matrix.train_targets,
            labels,
            config.test_size,
            config.random_state,
        )
    )
    split_seconds = time.perf_counter() - split_start

    cv_start = time.perf_counter()
    cv_folds, cv_warnings = shared.auto_cv_folds(
        labels_train, config.cv_folds, config.rare_class_policy
    )
    if cv_folds is None:
        raise ConfigError(
            "Hyperparameter tuning requires a feasible cross-validation split. "
            "Reduce cv_folds or adjust rare_class_policy so every class has at least two samples."
        )
    cv_seconds = time.perf_counter() - cv_start

    candidate_grid, _ = _build_candidate_grid(config.search_space)
    if config.search_method == "grid":
        if len(candidate_grid) > config.max_candidates:
            raise ConfigError(
                "Search space expands to "
                f"{len(candidate_grid)} candidates, which exceeds "
                f"max_candidates={config.max_candidates}. "
                "Switch to hyperparameter_tuning.method: random or shrink the search_space."
            )
        candidates = candidate_grid
    else:
        rng = random.Random(config.random_state)
        sample_size = min(config.n_iter, len(candidate_grid))
        candidates = rng.sample(candidate_grid, sample_size)

    total_candidates = len(candidates)
    estimated_model_fits = total_candidates * cv_folds
    _log_progress(
        "Will evaluate "
        f"{total_candidates} {_pluralize('candidate case', total_candidates)} across {cv_folds} CV folds "
        f"(~{estimated_model_fits} model fits)"
    )
    _log_progress(
        f"Setup completed in {_format_seconds(time.perf_counter() - run_start)}; starting candidate search"
    )

    search_start = time.perf_counter()
    evaluated_candidates: list[dict[str, object]] = []
    best_candidate_index: int | None = None
    best_candidate_params: dict[str, object] | None = None
    best_candidate_config: argparse.Namespace | None = None
    best_cv_results: dict[str, object] | None = None
    best_score: float | None = None
    objective_direction = config.objective_direction

    for candidate_index, candidate_params in enumerate(candidates, start=1):
        candidate_start = time.perf_counter()
        candidate_label = ", ".join(
            f"{key}={value!r}" for key, value in candidate_params.items()
        )
        _log_progress(
            f"[{candidate_index}/{total_candidates}] evaluating {candidate_label}"
        )
        candidate_config, cv_results, score = _evaluate_candidate(
            base_config,
            candidate_params,
            config.model_name,
            x_train,
            y_train,
            labels_train,
            cv_folds,
        )
        objective_key = f"{config.objective_metric}_mean"
        aggregate = cv_results["aggregate"]
        candidate_row = {
            "candidate_index": candidate_index,
            "objective_metric": config.objective_metric,
            "objective_direction": objective_direction,
            "cv_score": score,
            f"{config.objective_metric}_mean": aggregate.get(objective_key),
            f"{config.objective_metric}_std": aggregate.get(
                f"{config.objective_metric}_std"
            ),
            **candidate_params,
        }
        evaluated_candidates.append(candidate_row)

        if best_score is None:
            is_better = True
        elif objective_direction == "max":
            is_better = score > best_score
        else:
            is_better = score < best_score
        if is_better:
            best_score = score
            best_candidate_index = candidate_index
            best_candidate_params = candidate_params
            best_candidate_config = candidate_config
            best_cv_results = cv_results
        _log_progress(
            f"[{candidate_index}/{total_candidates}] done in {_format_seconds(time.perf_counter() - candidate_start)}; cv_score={score:.6f}"
        )
    search_seconds = time.perf_counter() - search_start

    _log_progress(
        f"Candidate search finished in {_format_seconds(search_seconds)}; best candidate so far is #{best_candidate_index}"
    )

    if any(
        value is None
        for value in (
            best_candidate_index,
            best_candidate_params,
            best_candidate_config,
            best_cv_results,
        )
    ):
        raise ConfigError("No candidates were evaluated")

    assert best_candidate_index is not None
    assert best_candidate_params is not None
    assert best_candidate_config is not None
    assert best_cv_results is not None

    if config.model_name == "multi_knn":
        best_model, effective_n_neighbors = knn_module._build_model(
            best_candidate_config, len(x_train)
        )
    else:
        best_model = rf_module._build_model(best_candidate_config)
        effective_n_neighbors = None

    fit_start = time.perf_counter()
    _log_progress("Refitting the best candidate on the training split")
    best_model.fit(x_train, y_train)
    test_predictions = best_model.predict(x_test)
    test_metrics = shared.evaluate_predictions(y_test, test_predictions)
    fit_predict_seconds = time.perf_counter() - fit_start
    _log_progress(
        f"Best candidate fit and prediction completed in {_format_seconds(fit_predict_seconds)}"
    )

    output_dir = Path(
        prepare_output_directory(config.output_dir, overwrite=config.overwrite)
    )

    dataset_summary = shared.build_dataset_summary(
        matrix.train_labels,
        list(labels_train),
        list(labels_test),
        y_train,
        y_test,
        split_notes,
        config.test_size,
        config.random_state,
    )
    prediction_rows = shared.build_prediction_rows(y_test, test_predictions)

    best_model_path = output_dir / "hyper_tune_best_model.pkl"
    results_json_path = output_dir / "hyper_tune_results.json"
    results_txt_path = output_dir / "hyper_tune_results.txt"
    results_tsv_path = output_dir / "hyper_tune_results.tsv"
    predictions_path = output_dir / "predictions.tsv"

    best_candidate_record = {
        "candidate_index": best_candidate_index,
        "candidate_params": best_candidate_params,
        "cv_score": best_score,
        "cv_results": best_cv_results,
    }
    if effective_n_neighbors is not None:
        best_candidate_record["effective_n_neighbors"] = int(effective_n_neighbors)

    results_payload = {
        "objective": "hyperparameter tuning",
        "model_name": config.model_name,
        "search_method": config.search_method,
        "objective_metric": config.objective_metric,
        "objective_direction": objective_direction,
        "candidate_count": len(evaluated_candidates),
        "split_notes": split_notes,
        "cv_notes": cv_warnings,
        "dataset_summary": dataset_summary,
        "best_candidate": best_candidate_record,
        "test_metrics": {
            "primary_metrics": {
                "hamming_loss": test_metrics["hamming_loss"],
                "micro_f1": test_metrics["micro_f1"],
                "macro_f1": test_metrics["macro_f1"],
                "weighted_f1": test_metrics["weighted_f1"],
                "bitwise_accuracy": test_metrics["bitwise_accuracy"],
            },
            "diagnostic_metrics": {
                "exact_match_accuracy": test_metrics["exact_match_accuracy"],
            },
            "per_bit": test_metrics["per_bit"],
            "classification_report": test_metrics["classification_report"],
        },
        "timings_seconds": {
            "load": load_seconds,
            "split": split_seconds,
            "cv_feasibility": cv_seconds,
            "candidate_search": search_seconds,
            "fit_and_predict": fit_predict_seconds,
        },
    }

    if config.model_name == "multi_knn":
        results_payload["best_candidate"]["effective_n_neighbors"] = (
            effective_n_neighbors
        )

    with open(best_model_path, "wb") as handle:
        pickle.dump(best_model, handle)

    artifact_start = time.perf_counter()
    _log_progress("Writing tuning artifacts")
    shared.write_json(results_json_path, results_payload)
    shared.write_tsv(predictions_path, prediction_rows)

    sorted_candidates = sorted(
        evaluated_candidates,
        key=lambda row: row["cv_score"],
        reverse=(objective_direction == "max"),
    )
    for rank, row in enumerate(sorted_candidates, start=1):
        row["rank"] = rank
        row["is_best"] = rank == 1

    shared.write_tsv(results_tsv_path, sorted_candidates)
    artifact_seconds = time.perf_counter() - artifact_start
    results_payload["timings_seconds"]["artifact_write"] = artifact_seconds
    results_payload["timings_seconds"]["total"] = time.perf_counter() - run_start
    shared.write_json(results_json_path, results_payload)
    _log_progress(
        f"Finished in {_format_seconds(results_payload['timings_seconds']['total'])}; artifacts written to {output_dir}"
    )

    text_lines = [
        f"Ghostparser ML hyperparameter tuning ({config.model_name})",
        f"Search method: {config.search_method}",
        f"Objective metric: {config.objective_metric} ({objective_direction})",
        f"Candidates evaluated: {len(evaluated_candidates)}",
        "",
        "Best candidate:",
        f"  Candidate index: {best_candidate_record['candidate_index']}",
    ]
    if effective_n_neighbors is not None:
        text_lines.append(f"  Effective n_neighbors: {effective_n_neighbors}")
    for key, value in best_candidate_record["candidate_params"].items():
        text_lines.append(f"  {key}: {value}")
    text_lines.extend(
        [
            "",
            "Test metrics:",
            f"  Hamming loss: {test_metrics['hamming_loss']:.6f}",
            f"  Bitwise accuracy: {test_metrics['bitwise_accuracy']:.6f}",
            f"  Exact-match accuracy: {test_metrics['exact_match_accuracy']:.6f}",
            f"  Micro F1: {test_metrics['micro_f1']:.6f}",
            f"  Macro F1: {test_metrics['macro_f1']:.6f}",
            f"  Weighted F1: {test_metrics['weighted_f1']:.6f}",
            "",
            "Dataset summary:",
            "  Label map:",
        ]
    )
    for key, value in dataset_summary["label_map"].items():
        text_lines.append(f"    {key}: {value}")
    text_lines.append("  Split:")
    for key, value in dataset_summary["split"].items():
        text_lines.append(f"    {key}: {value}")
    text_lines.extend(
        [
            "",
            "Top candidates:",
        ]
    )
    for row in sorted_candidates[: config.top_k]:
        candidate_pairs = ", ".join(
            f"{key}={value}" for key, value in row.items() if key in config.search_space
        )
        text_lines.append(
            f"  rank={row['rank']} score={row['cv_score']:.6f} {candidate_pairs}"
        )
    if cv_warnings:
        text_lines.extend(["", "CV notes:"] + [f"  {note}" for note in cv_warnings])
    if split_notes:
        text_lines.extend(["", "Split notes:"] + [f"  {note}" for note in split_notes])
    text_lines.extend(
        [
            "",
            "Timings (seconds):",
            f"  load: {load_seconds:.6f}",
            f"  split: {split_seconds:.6f}",
            f"  cv_feasibility: {cv_seconds:.6f}",
            f"  candidate_search: {search_seconds:.6f}",
            f"  fit_and_predict: {fit_predict_seconds:.6f}",
            f"  artifact_write: {artifact_seconds:.6f}",
            f"  total: {results_payload['timings_seconds']['total']:.6f}",
        ]
    )
    shared.write_text(results_txt_path, "\n".join(text_lines) + "\n")

    return {
        "best_model_path": str(best_model_path),
        "results_json_path": str(results_json_path),
        "results_txt_path": str(results_txt_path),
        "results_tsv_path": str(results_tsv_path),
        "predictions_path": str(predictions_path),
        "results": results_payload,
        "candidates": sorted_candidates,
    }


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ghostparser ML hyperparameter tuner for supported model modules."
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        required=True,
        help="Path to a JSON or YAML hyperparameter tuning config file",
    )
    parser.add_argument(
        "--no-overwrite",
        dest="no_overwrite",
        action="store_true",
        default=False,
        help="Append a numeric suffix when the output directory already exists",
    )
    return parser


def main() -> None:
    args = _build_argument_parser().parse_args()
    config = load_hyper_tune_config(args.config_file)
    if args.no_overwrite:
        config["overwrite"] = False
    result = tune_hyperparameters(argparse.Namespace(**config))
    print(result["results_txt_path"])


if __name__ == "__main__":
    main()
__all__ = [
    "load_hyper_tune_config",
    "normalize_hyper_tune_payload",
    "tune_hyperparameters",
    "main",
]
