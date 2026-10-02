"""Config-driven grid or random hyperparameter search over the ML trainers.

Settings live under the ``hyperparameter_tuning`` config section, documented
under "Configuration" in the ML guide.
"""

import argparse
import itertools
import json
import os
import pickle
import random
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from ..config import (
    DEFAULT_OVERWRITE,
    ConfigError,
    _validate_overwrite_flag,
    _validate_required_path,
    prepare_output_directory,
)
from ..cli_config import resolve_cli_or_config_args, run_cli
from . import ml_utils as shared
from . import multi_knn as knn_module
from . import random_forest as rf_module
from . import tuning_report
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
    DEFAULT_SEED,
    DEFAULT_TARGET_COLUMN,
    _load_raw_config,
    _validate_optional_bool,
    _validate_optional_choice,
    _validate_optional_float,
    _validate_optional_int,
    _validate_optional_positive_int,
    _validate_optional_string,
    normalize_class_weight,
    normalize_max_features,
)

# Search-space values reach the estimator one candidate at a time, so the same
# per-value rules the `model` block enforces have to run across each list here.
SEARCH_VALUE_NORMALIZERS = {
    "max_features": normalize_max_features,
    "class_weight": normalize_class_weight,
}

DEFAULT_TUNER_MODEL = "random_forest"
DEFAULT_SEARCH_METHOD = "grid"
DEFAULT_OBJECTIVE_METRIC = "exact_match_accuracy"
DEFAULT_TOP_K = 10
DEFAULT_RANDOM_ITERATIONS = 20
DEFAULT_MAX_CANDIDATES = 5000
DEFAULT_WANDB_PROJECT = "ghostparser-hyper-tune"
DEFAULT_USE_WANDB = False

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
    """Format a duration in seconds for the progress log."""
    return f"{seconds:.2f}s"


def _pluralize(word: str, count: int) -> str:
    """Pluralize ``word`` unless ``count`` is one."""
    return word if count == 1 else f"{word}s"


def _log_progress(message: str) -> None:
    """Print a progress line prefixed with the tuner's name."""
    print(f"[hyper_tune] {message}", flush=True)


def _import_wandb():
    """Import Weights & Biases only when a run actually opted into it."""
    try:
        import wandb
    except ImportError as exc:
        raise ConfigError(
            "hyperparameter_tuning.use_wandb is true but the 'wandb' package could "
            "not be imported. Install it with 'pip install .[wandb]', or set "
            "hyperparameter_tuning.use_wandb: false to run the tuner with local "
            "plaintext, TSV, and plot reporting only."
        ) from exc
    return wandb


def _get_wandb_project() -> str:
    """Read the W&B project from ``WANDB_PROJECT``, falling back to the default."""
    value = os.getenv("WANDB_PROJECT")
    if value is None:
        return DEFAULT_WANDB_PROJECT
    stripped = value.strip()
    return stripped or DEFAULT_WANDB_PROJECT


def _get_wandb_entity() -> str | None:
    """Read the W&B entity from ``WANDB_ENTITY``, or ``None``."""
    value = os.getenv("WANDB_ENTITY")
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def _build_wandb_run_name(config: argparse.Namespace) -> str:
    """Name the W&B run after the model, the search method and the time."""
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return (
        f"ghostparser-{config.model_name}-{config.search_method}-"
        f"{config.objective_metric}-{timestamp}"
    )


class _NullRunLogger:
    """Run logger used when ``use_wandb`` is false; every call is a no-op."""

    enabled = False
    detailed_payloads = False

    def log(self, payload: dict) -> None:
        """Discard the payload."""
        return None

    def log_table(self, name: str, rows: list[dict[str, object]]) -> None:
        """Discard the table."""
        return None

    def set_summary(self, key: str, value: object) -> None:
        """Discard the summary value."""
        return None

    def finish(self) -> None:
        """Do nothing."""
        return None


class _WandbRunLogger:
    """Run logger that forwards tuning metrics to a Weights & Biases run."""

    enabled = True

    def __init__(self, wandb_module, run, detailed_payloads: bool) -> None:
        """Wrap a live W&B run."""
        self._wandb = wandb_module
        self._run = run
        self.detailed_payloads = detailed_payloads

    def log(self, payload: dict) -> None:
        """Log a payload to the run."""
        self._wandb.log(payload)

    def log_table(self, name: str, rows: list[dict[str, object]]) -> None:
        """Send a row-oriented table to the run as a sortable W&B table.

        Args:
            name: Key the table appears under in the run.
            rows: Row dicts; the union of their keys becomes the columns, in
                first-seen order. An empty list logs nothing.
        """
        if not rows:
            return
        columns: list[str] = []
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)
        table = self._wandb.Table(
            columns=columns,
            data=[[row.get(column) for column in columns] for row in rows],
        )
        self._wandb.log({name: table})

    def set_summary(self, key: str, value: object) -> None:
        """Set a run-summary value."""
        self._run.summary[key] = value

    def finish(self) -> None:
        """Finish the run."""
        self._run.finish()


def _create_run_logger(
    config: argparse.Namespace,
    output_dir: Path,
    total_candidates: int,
    cv_folds: int,
    use_wandb: bool,
) -> _NullRunLogger | _WandbRunLogger:
    """Build the run logger the config asked for.

    Args:
        config: Resolved tuning configuration.
        output_dir: Run output directory; hosts the ``wandb/`` scratch folder.
        total_candidates: Number of candidates the search will evaluate.
        cv_folds: Effective cross-validation fold count.
        use_wandb: Resolved ``hyperparameter_tuning.use_wandb`` choice.

    Returns:
        A no-op logger when W&B is off, otherwise a live W&B-backed logger.
    """
    if not use_wandb:
        _log_progress(
            "Weights & Biases logging is disabled; writing local plaintext, TSV, "
            "and plot reports only"
        )
        return _NullRunLogger()

    wandb = _import_wandb()
    wandb_dir = output_dir / "wandb"
    wandb_dir.mkdir(parents=True, exist_ok=True)

    run_config = {
        "model_name": config.model_name,
        "search_method": config.search_method,
        "objective_metric": config.objective_metric,
        "objective_direction": config.objective_direction,
        "candidate_count": int(total_candidates),
        "cv_folds": int(cv_folds),
        "test_size": float(config.test_size),
        "rare_class_policy": config.rare_class_policy,
        "seed": config.seed,
        "n_jobs": config.n_jobs,
    }

    try:
        run = wandb.init(
            project=_get_wandb_project(),
            entity=_get_wandb_entity(),
            name=_build_wandb_run_name(config),
            dir=str(wandb_dir),
            config=run_config,
            tags=["ghostparser", "hyper_tune", config.model_name, config.search_method],
            mode=os.getenv("WANDB_MODE", "online"),
        )
    except Exception as exc:  # noqa: BLE001
        raise ConfigError(
            "Failed to initialize Weights & Biases for hyperparameter tuning. "
            "Run 'wandb login' first, set WANDB_MODE=offline if network is "
            "unavailable, or set hyperparameter_tuning.use_wandb: false."
        ) from exc

    wandb.define_metric("candidate_index")
    wandb.define_metric("candidate/*", step_metric="candidate_index")
    return _WandbRunLogger(
        wandb,
        run,
        bool(getattr(config, "wandb_detailed_payloads", False)),
    )


RUNTIME_KEYS = {
    "input_path",
    "output_dir",
    "overwrite",
    "target_column",
    "test_size",
    "cv_folds",
    "rare_class_policy",
    "seed",
    "n_jobs",
}


def _normalize_search_space_values(search_space: dict[str, object]) -> dict:
    """Apply each key's per-value rules across its candidate list.

    Args:
        search_space: Raw search space, each key mapping to a candidate list or
            a lone value.

    Returns:
        The same mapping with every value of a validated key normalized, and
        the container shape (list or scalar) left as it was written.

    Raises:
        ConfigError: If any candidate value is invalid for its key.
    """
    normalized: dict[str, object] = {}
    for key, value in search_space.items():
        normalizer = SEARCH_VALUE_NORMALIZERS.get(key)
        if normalizer is None:
            normalized[key] = value
            continue
        where = f"hyperparameter_tuning.search_space.{key}"
        if isinstance(value, (list, tuple)):
            normalized[key] = [normalizer(item, where=where) for item in value]
        else:
            normalized[key] = normalizer(value, where=where)
    return normalized


def _normalize_search_values(value: object) -> list[object]:
    """Turn a search-space entry into a non-empty candidate list, accepting a lone value.
    """
    if isinstance(value, (list, tuple)):
        candidates = list(value)
        if not candidates:
            raise ConfigError("Search space lists must contain at least one value")
        return candidates
    return [value]


def _build_candidate_grid(
    search_space: dict[str, object],
) -> tuple[list[dict[str, object]], list[str]]:
    """Expand the search space into every candidate combination.

    Args:
        search_space: Parameter name to candidate list.

    Returns:
        ``(candidates, parameter_names)``: one dict per combination, and the
        parameter names in order.

    Raises:
        ConfigError: If the search space is empty or a list has no values.
    """
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
    """Return ``min`` for the objectives that are losses, else ``max``."""
    return "min" if metric_name in MINIMIZE_OBJECTIVES else "max"


def _candidate_score(cv_results: dict, objective_metric: str) -> float:
    """Read the candidate's cross-validated objective mean.

    Raises:
        ConfigError: If the objective is not among the CV metrics.
    """
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
    """Build a trainer config from the runtime keys, the model defaults and one candidate's parameters.
    """
    payload = dict(base_config)
    payload.update(MODEL_DEFAULTS[model_name])
    payload.update(candidate_params)
    return argparse.Namespace(**payload)


def load_hyper_tune_config(config_file: str) -> dict[str, object]:
    """Load and normalize a tuner config file."""
    payload = _load_raw_config(config_file)
    return normalize_hyper_tune_payload(payload)


def normalize_hyper_tune_payload(payload: dict) -> dict[str, object]:
    """Validate a tuner config: the runtime keys plus the ``hyperparameter_tuning`` block.

    Args:
        payload: The raw config mapping.

    Returns:
        The resolved tuner config dict.

    Raises:
        ConfigError: If a key is missing, misplaced or invalid.
    """
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
            "Do not place model sections, evaluation controls, or model "
            "hyperparameters at the top level in hyperparameter_tuning configs. "
            "Use only the runtime keys plus 'hyperparameter_tuning'. "
            "Offending keys: "
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
    use_wandb = _validate_optional_bool(
        tuning_section,
        "use_wandb",
        DEFAULT_USE_WANDB,
    )
    wandb_detailed_payloads = _validate_optional_bool(
        tuning_section,
        "wandb_detailed_payloads",
        False,
    )
    if wandb_detailed_payloads and not use_wandb:
        raise ConfigError(
            "Config field hyperparameter_tuning.wandb_detailed_payloads requires "
            "hyperparameter_tuning.use_wandb: true"
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
    search_space = _normalize_search_space_values(search_space)

    normalized_config = {
        "input_path": input_path,
        "output_dir": output_dir,
        "overwrite": _validate_overwrite_flag(payload, DEFAULT_OVERWRITE),
        "target_column": target_column,
        "test_size": _validate_optional_float(payload, "test_size", 0.2),
        "cv_folds": _validate_optional_positive_int(payload, "cv_folds", 5),
        "rare_class_policy": _validate_optional_choice(
            payload,
            "rare_class_policy",
            "warn_reduce_cv",
            ("warn_reduce_cv", "warn_skip_cv", "error"),
        ),
        "seed": _validate_optional_int(payload, "seed", DEFAULT_SEED),
        "n_jobs": _validate_optional_int(payload, "n_jobs", DEFAULT_N_JOBS),
        "model_name": model_name,
        "search_method": search_method,
        "objective_metric": objective_metric,
        "objective_direction": _objective_direction(objective_metric),
        "top_k": top_k,
        "n_iter": n_iter,
        "max_candidates": max_candidates,
        "use_wandb": use_wandb,
        "wandb_detailed_payloads": wandb_detailed_payloads,
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
    """Cross-validate one candidate and return its score and results.

    Args:
        base_config: The runtime keys.
        candidate_params: The candidate's hyperparameters.
        model_name: ``random_forest`` or ``multi_knn``.
        x_train: The training features.
        y_train: The training binary targets.
        labels_train: The training combination labels.
        folds: Number of folds.

    Returns:
        ``(candidate_config, cv_results, score)``.
    """
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
            candidate_config.seed,
        )
    elif model_name == "multi_knn":
        cv_results = knn_module._cross_validate(
            x_train,
            y_train,
            labels_train,
            folds,
            candidate_config.seed,
            candidate_config,
        )
    else:
        raise ConfigError(f"Unsupported tuning model: {model_name}")
    score = _candidate_score(cv_results, candidate_config.objective_metric)
    return candidate_config, cv_results, score


def tune_hyperparameters(config: argparse.Namespace) -> dict[str, object]:
    """Search the hyperparameter space, refit the best candidate and write the reports.

    Args:
        config: The resolved tuner config.

    Returns:
        A dict with the ranked candidates, the best parameters, the test
        metrics, the parameter marginals and the artifact paths.
    """
    run_start = time.perf_counter()
    use_wandb = bool(getattr(config, "use_wandb", DEFAULT_USE_WANDB))
    # The ranked-candidate table, the parameter marginals, the per-row
    # predictions and the full results payload are the bulky outputs. With W&B
    # on they go to the run instead of the output directory, which then keeps
    # only the model pickle, the plaintext report and the search plot.
    write_bulk_artifacts = not use_wandb
    _log_progress(
        f"Starting hyperparameter tuning for {config.model_name} with {config.search_method} search"
    )
    load_start = time.perf_counter()
    base_config = dict(vars(config))
    base_config.update(MODEL_DEFAULTS[config.model_name])
    output_dir = Path(
        prepare_output_directory(config.output_dir, overwrite=config.overwrite)
    )
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
            config.seed,
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

    candidate_grid, search_parameter_keys = _build_candidate_grid(config.search_space)
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
        rng = random.Random(config.seed)
        sample_size = min(config.n_iter, len(candidate_grid))
        candidates = rng.sample(candidate_grid, sample_size)

    total_candidates = len(candidates)
    estimated_model_fits = total_candidates * cv_folds
    _log_progress(
        "Will evaluate "
        f"{total_candidates} {_pluralize('candidate case', total_candidates)} across {cv_folds} CV folds "
        f"(~{estimated_model_fits} model fits)"
    )

    run_logger = _create_run_logger(
        config, output_dir, total_candidates, cv_folds, use_wandb
    )
    detailed_wandb_logging = run_logger.detailed_payloads
    run_logger.log(
        {
            "candidate_index": 0,
            "candidate/total": total_candidates,
            "candidate/estimated_model_fits": estimated_model_fits,
            "timing/setup_seconds": time.perf_counter() - run_start,
        }
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
        candidate_seconds = time.perf_counter() - candidate_start
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
            "elapsed_seconds": candidate_seconds,
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

        candidate_log = {
            "candidate_index": candidate_index,
            "candidate/cv_score": float(score),
            "candidate/is_best": int(is_better),
            "candidate/elapsed_seconds": candidate_seconds,
            "candidate/remaining": total_candidates - candidate_index,
        }
        candidate_log["candidate/bitwise_accuracy_mean"] = float(
            aggregate.get("bitwise_accuracy_mean", 0.0)
        )
        candidate_log["candidate/overall_accuracy_mean"] = float(
            aggregate.get("exact_match_accuracy_mean", 0.0)
        )
        candidate_log.update(
            {f"candidate/param/{key}": value for key, value in candidate_params.items()}
        )
        run_logger.log(candidate_log)
        if detailed_wandb_logging:
            run_logger.log(
                {
                    "candidate_index": candidate_index,
                    "candidate/detailed/params_json": json.dumps(
                        candidate_params,
                        sort_keys=True,
                    ),
                    "candidate/detailed/cv_aggregate_json": json.dumps(
                        cv_results.get("aggregate", {}),
                        sort_keys=True,
                    ),
                    "candidate/detailed/cv_folds_json": json.dumps(
                        cv_results.get("folds", []),
                        sort_keys=True,
                    ),
                }
            )

        _log_progress(
            f"[{candidate_index}/{total_candidates}] done in "
            f"{_format_seconds(candidate_seconds)}; "
            f"cv_score={score:.6f}"
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

    dataset_summary = shared.build_dataset_summary(
        matrix.train_labels,
        list(labels_train),
        list(labels_test),
        y_train,
        y_test,
        split_notes,
        config.test_size,
        config.seed,
    )
    prediction_rows = shared.build_prediction_rows(y_test, test_predictions)

    best_model_path = output_dir / "hyper_tune_best_model.pkl"
    results_json_path = output_dir / "hyper_tune_results.json"
    results_txt_path = output_dir / "hyper_tune_results.txt"
    results_tsv_path = output_dir / "hyper_tune_results.tsv"
    marginals_tsv_path = output_dir / "hyper_tune_parameter_marginals.tsv"
    plot_path = output_dir / "figures" / "hyper_tune_search_report.png"
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
        "use_wandb": use_wandb,
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
    if write_bulk_artifacts:
        _log_progress("Writing tuning artifacts")
    else:
        _log_progress(
            "Weights & Biases logging is on; sending the ranked candidates, "
            "parameter marginals, predictions and full results payload to the "
            "run instead of writing them to the output directory"
        )

    sorted_candidates = sorted(
        evaluated_candidates,
        key=lambda row: row["cv_score"],
        reverse=(objective_direction == "max"),
    )
    for rank, row in enumerate(sorted_candidates, start=1):
        row["rank"] = rank
        row["is_best"] = rank == 1

    marginals = tuning_report.compute_parameter_marginals(
        sorted_candidates, search_parameter_keys, objective_direction
    )
    marginal_rows = tuning_report.marginal_tsv_rows(marginals)

    if write_bulk_artifacts:
        shared.write_json(results_json_path, results_payload)
        shared.write_tsv(predictions_path, prediction_rows)
        shared.write_tsv(results_tsv_path, sorted_candidates)
        shared.write_tsv(marginals_tsv_path, marginal_rows)
    else:
        run_logger.log_table("tables/ranked_candidates", sorted_candidates)
        run_logger.log_table("tables/parameter_marginals", marginal_rows)
        run_logger.log_table("tables/predictions", prediction_rows)

    saved_plot_path = tuning_report.save_tuning_plots(
        sorted_candidates,
        search_parameter_keys,
        config.objective_metric,
        objective_direction,
        plot_path,
        top_k=config.top_k,
    )

    artifact_seconds = time.perf_counter() - artifact_start
    results_payload["timings_seconds"]["artifact_write"] = artifact_seconds
    results_payload["timings_seconds"]["total"] = time.perf_counter() - run_start
    results_payload["parameter_marginals"] = marginal_rows
    results_payload["parameter_influence"] = tuning_report.parameter_influence(
        marginals
    )

    artifact_paths = {
        "best_model": str(best_model_path),
        "results_txt": str(results_txt_path),
    }
    if write_bulk_artifacts:
        artifact_paths["results_json"] = str(results_json_path)
        artifact_paths["ranked_candidates_tsv"] = str(results_tsv_path)
        artifact_paths["parameter_marginals_tsv"] = str(marginals_tsv_path)
        artifact_paths["predictions_tsv"] = str(predictions_path)
    if saved_plot_path is not None:
        artifact_paths["search_report_plot"] = saved_plot_path
    results_payload["artifact_paths"] = artifact_paths

    shared.write_text(
        results_txt_path,
        tuning_report.build_results_text(
            config=config,
            ranked_candidates=sorted_candidates,
            marginals=marginals,
            parameter_keys=search_parameter_keys,
            best_candidate_record=best_candidate_record,
            effective_n_neighbors=effective_n_neighbors,
            test_metrics=test_metrics,
            dataset_summary=dataset_summary,
            cv_warnings=cv_warnings,
            split_notes=split_notes,
            timings={
                "load": load_seconds,
                "split": split_seconds,
                "cv_feasibility": cv_seconds,
                "candidate_search": search_seconds,
                "fit_and_predict": fit_predict_seconds,
                "artifact_write": artifact_seconds,
            },
            total_seconds=results_payload["timings_seconds"]["total"],
            cv_folds=cv_folds,
            wandb_enabled=use_wandb,
            artifact_paths=artifact_paths,
        ),
    )
    if write_bulk_artifacts:
        shared.write_json(results_json_path, results_payload)
    else:
        run_logger.set_summary(
            "results_json", json.dumps(results_payload, sort_keys=True, default=str)
        )

    run_logger.log(
        {
            "timing/load_seconds": load_seconds,
            "timing/split_seconds": split_seconds,
            "timing/cv_feasibility_seconds": cv_seconds,
            "timing/candidate_search_seconds": search_seconds,
            "timing/fit_and_predict_seconds": fit_predict_seconds,
            "timing/artifact_write_seconds": artifact_seconds,
            "timing/total_seconds": results_payload["timings_seconds"]["total"],
            "final/exact_match_accuracy": test_metrics["exact_match_accuracy"],
            "final/overall_accuracy": test_metrics["exact_match_accuracy"],
            "final/hamming_loss": test_metrics["hamming_loss"],
            "final/bitwise_accuracy": test_metrics["bitwise_accuracy"],
            "final/micro_f1": test_metrics["micro_f1"],
            "final/macro_f1": test_metrics["macro_f1"],
            "final/weighted_f1": test_metrics["weighted_f1"],
            "final/best_candidate_index": best_candidate_index,
        }
    )
    run_logger.set_summary("artifact_dir", str(output_dir))
    run_logger.set_summary("best_candidate_params", dict(best_candidate_params))
    # Only advertise paths that were actually written, so a W&B run never points
    # at a file the run deliberately skipped.
    run_logger.set_summary(
        "local_artifact_paths", json.dumps(artifact_paths, sort_keys=True)
    )
    run_logger.set_summary("bulk_artifacts_written_locally", write_bulk_artifacts)
    run_logger.set_summary("wandb_detailed_payloads", detailed_wandb_logging)
    if detailed_wandb_logging:
        run_logger.set_summary(
            "dataset_summary_json", json.dumps(dataset_summary, sort_keys=True)
        )
        run_logger.set_summary(
            "best_candidate_cv_results_json",
            json.dumps(best_cv_results, sort_keys=True),
        )
    run_logger.finish()
    _log_progress(
        f"Finished in {_format_seconds(results_payload['timings_seconds']['total'])}; artifacts written to {output_dir}"
    )
    _log_progress(f"Ranked candidate report: {results_txt_path}")
    if write_bulk_artifacts:
        _log_progress(f"Per-parameter value summary: {marginals_tsv_path}")
    if saved_plot_path is not None:
        _log_progress(f"Search report plot: {saved_plot_path}")

    return {
        "best_model_path": str(best_model_path),
        "results_json_path": str(results_json_path) if write_bulk_artifacts else None,
        "results_txt_path": str(results_txt_path),
        "results_tsv_path": str(results_tsv_path) if write_bulk_artifacts else None,
        "marginals_tsv_path": (
            str(marginals_tsv_path) if write_bulk_artifacts else None
        ),
        "plot_path": saved_plot_path,
        "predictions_path": str(predictions_path) if write_bulk_artifacts else None,
        "results": results_payload,
        "candidates": sorted_candidates,
        "parameter_marginals": marginals,
    }


def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the tuner's argument parser."""
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
        "--seed",
        type=int,
        default=None,
        help="RNG seed for the split, the folds, the models and a random "
        "search's draws; overrides the config file's seed",
    )
    parser.add_argument(
        "--no-overwrite",
        dest="no_overwrite",
        action="store_true",
        default=None,
        help="Append a numeric suffix when the output directory already exists",
    )
    return parser


def _run(parsed_args) -> None:
    """Resolve the config, run the search, and print where the results went.

    Args:
        parsed_args: The parsed command line.
    """
    config = resolve_cli_or_config_args(
        parsed_args,
        normalize_payload=normalize_hyper_tune_payload,
        payload_arg_names=["seed", "no_overwrite"],
    )
    result = tune_hyperparameters(config)
    print(result["results_txt_path"])


def main() -> None:
    """Run the tuner from the command line and exit with its status."""
    sys.exit(run_cli(_build_argument_parser(), _run))


if __name__ == "__main__":
    main()
__all__ = [
    "load_hyper_tune_config",
    "normalize_hyper_tune_payload",
    "tune_hyperparameters",
    "main",
]
