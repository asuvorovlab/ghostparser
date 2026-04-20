"""Configuration loading and normalization helpers for GhostParser CLIs."""

from __future__ import annotations

import json
from pathlib import Path


class ConfigError(ValueError):
    """Raised when a config payload is invalid or missing required fields."""


def _resolve_path(path_str: str) -> str:
    """Resolve a path string to an absolute path, handling ~, relative, and absolute paths.
    
    Args:
        path_str: Path string that can be:
                 - Absolute (starts with /): /path/to/file
                 - Relative (no leading /): path/to/file (resolved from cwd)
                 - User home (starts with ~): ~/path/to/file
    
    Returns:
        Absolute path as a string
    """
    path = Path(path_str)
    
    # Expand user home directory (~)
    path = path.expanduser()
    
    # Resolve to absolute path
    # If already absolute, this keeps it as is
    # If relative, resolves from current working directory
    path = path.resolve()
    
    return str(path)


DEFAULT_OUTPUT_FOLDER = "results"
DEFAULT_PROCESSES = 0
DEFAULT_MIN_SUPPORT_VALUE = 0.5
DEFAULT_DISCORDANT_TEST = "chi-square"
DEFAULT_SUMMARY_STATISTIC = "median"
DEFAULT_STATS_BACKEND = "standard"
DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY = "AVG"
DEFAULT_P_VALUE_CORRECTION = "no"
DEFAULT_ALPHA_DCT = 0.05
DEFAULT_ALPHA_KS = 0.05
DEFAULT_BOOTSTRAP = True
DEFAULT_BOOTSTRAP_ITERATIONS = 100
DEFAULT_BOOTSTRAP_DEBUG_MODE = False
DEFAULT_BOOTSTRAP_SUMMARY_ONLY = False
DEFAULT_GENERATE_SUMMARY_STATS = False
DEFAULT_TRIPLET_OUTPUT_FORMAT = "parquet"
DEFAULT_INPUT_FORMAT = "parquet"
DEFAULT_PARQUET_PARTITIONS = 128
DEFAULT_PARQUET_COMPRESSION = "zstd"

TRIPLET_IO_FORMAT_CHOICES = ("txt", "parquet")
INPUT_FORMAT_CHOICES = ("auto", "txt", "parquet")
PARQUET_COMPRESSION_CHOICES = ("zstd", "snappy", "gzip", "brotli", "none")

DISCORDANT_TEST_CHOICES = ("chi-square", "z-test")
SUMMARY_STATISTIC_CHOICES = ("mean", "median", "mode")
STATS_BACKEND_CHOICES = ("custom", "standard")
TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES = ("AVG", "A", "B", "C", "SIS", "INT")
P_VALUE_CORRECTION_CHOICES = ("no", "bfn", "holm", "fdr_bh", "fdr_by", "fdr_tsbh")


def _load_raw_config(config_file: str) -> dict:
    """Load a raw config dictionary from JSON or YAML."""
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


def _validate_required_string(payload: dict, key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return value.strip()


def _validate_required_path(payload: dict, key: str) -> str:
    """Validate and resolve a required path field."""
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return _resolve_path(value.strip())


def _validate_optional_string(payload: dict, key: str) -> str | None:
    value = payload.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Config field {key} must be a non-empty string when provided")
    return value.strip()


def _validate_optional_path(payload: dict, key: str) -> str | None:
    """Validate and resolve an optional path field."""
    value = payload.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Config field {key} must be a non-empty string when provided")
    return _resolve_path(value.strip())


def _validate_non_negative_int(payload: dict, key: str, default: int) -> int:
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, int) or value < 0:
        raise ConfigError(f"Config field {key} must be an integer >= 0")
    return value


def _validate_optional_bool(payload: dict, key: str, default: bool) -> bool:
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, bool):
        raise ConfigError(f"Config field {key} must be a boolean when provided")
    return value


def _validate_optional_float(payload: dict, key: str, default: float) -> float:
    value = payload.get(key, default)
    if value is None:
        value = default
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"Config field {key} must be a numeric value") from exc


def _validate_choice(payload: dict, key: str, default: str, choices: tuple[str, ...]) -> str:
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, str) or value not in choices:
        raise ConfigError(f"Config field {key} must be one of: {', '.join(choices)}")
    return value


def _validate_bootstrap_options(payload: dict) -> tuple[bool, dict]:
    """Validate bootstrap toggle and nested options.

    The canonical shape is:
        bootstrap: bool
        bootstrap_options:
          iterations: int >= 1
          seed: int | null
                    debug_mode: bool
          summary_only: bool

    CLI payloads may provide flat keys (`bootstrap_iterations`, `bootstrap_seed`,
    `bootstrap_debug_mode`, `bootstrap_summary_only`), which are merged into `bootstrap_options`.
    """
    bootstrap = _validate_optional_bool(payload, "bootstrap", DEFAULT_BOOTSTRAP)

    raw_options = payload.get("bootstrap_options")
    if raw_options is None:
        raw_options = {}
    if not isinstance(raw_options, dict):
        raise ConfigError("Config field bootstrap_options must be a key/value object when provided")

    iterations = payload.get("bootstrap_iterations", raw_options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    if iterations is None:
        iterations = DEFAULT_BOOTSTRAP_ITERATIONS
    if not isinstance(iterations, int) or iterations < 1:
        raise ConfigError("Config field bootstrap_options.iterations must be an integer >= 1")

    seed = payload.get("bootstrap_seed", raw_options.get("seed"))
    if seed is not None and not isinstance(seed, int):
        raise ConfigError("Config field bootstrap_options.seed must be an integer when provided")

    debug_mode = payload.get(
        "bootstrap_debug_mode",
        raw_options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE),
    )
    if debug_mode is None:
        debug_mode = DEFAULT_BOOTSTRAP_DEBUG_MODE
    if not isinstance(debug_mode, bool):
        raise ConfigError("Config field bootstrap_options.debug_mode must be a boolean when provided")

    summary_only = payload.get(
        "bootstrap_summary_only",
        raw_options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY),
    )
    if summary_only is None:
        summary_only = DEFAULT_BOOTSTRAP_SUMMARY_ONLY
    if not isinstance(summary_only, bool):
        raise ConfigError("Config field bootstrap_options.summary_only must be a boolean when provided")

    return bootstrap, {
        "iterations": iterations,
        "seed": seed,
        "debug_mode": debug_mode,
        "summary_only": summary_only,
    }


def _parse_outgroups(value) -> list[str]:
    """Parse outgroup(s) value from payload.

    Accepts either:
    - a comma-separated string
    - a list/tuple/set of strings
    """
    if isinstance(value, str):
        outgroups = [part.strip() for part in value.split(",") if part.strip()]
    elif isinstance(value, (list, tuple, set)):
        outgroups = [str(part).strip() for part in value if str(part).strip()]
    else:
        outgroups = []

    if not outgroups:
        raise ConfigError("Missing required config field: outgroup(s)")
    return outgroups


def normalize_orchestrator_payload(payload: dict) -> dict:
    """Normalize orchestrator config/CLI payload to internal runtime keys with defaults."""
    species_tree = _validate_required_path(payload, "species_tree_path")
    gene_trees = _validate_required_path(payload, "gene_trees_path")

    outgroups_source = payload.get("outgroups")
    if outgroups_source is None:
        outgroups_source = payload.get("outgroup")
    outgroups = _parse_outgroups(outgroups_source)

    output = _validate_optional_path(payload, "output_folder")
    if output is None:
        output = _resolve_path(DEFAULT_OUTPUT_FOLDER)

    bootstrap, bootstrap_options = _validate_bootstrap_options(payload)

    return {
        "species_tree": species_tree,
        "gene_trees": gene_trees,
        "outgroup": outgroups,
        "triplet_filter": _validate_optional_path(payload, "triplet_filter"),
        "output": output,
        "triplet_output_format": _validate_choice(
            payload,
            "triplet_output_format",
            DEFAULT_TRIPLET_OUTPUT_FORMAT,
            TRIPLET_IO_FORMAT_CHOICES,
        ),
        "parquet_partitions": _validate_non_negative_int(
            payload,
            "parquet_partitions",
            DEFAULT_PARQUET_PARTITIONS,
        ),
        "parquet_compression": _validate_choice(
            payload,
            "parquet_compression",
            DEFAULT_PARQUET_COMPRESSION,
            PARQUET_COMPRESSION_CHOICES,
        ),
        "processes": _validate_non_negative_int(payload, "processes", DEFAULT_PROCESSES),
        "generate_summary_stats": _validate_optional_bool(
            payload,
            "generate_summary_stats",
            DEFAULT_GENERATE_SUMMARY_STATS,
        ),
        "min_support_value": _validate_optional_float(payload, "min_support_value", DEFAULT_MIN_SUPPORT_VALUE),
        "discordant_test": _validate_choice(
            payload,
            "discordant_test",
            DEFAULT_DISCORDANT_TEST,
            DISCORDANT_TEST_CHOICES,
        ),
        "summary_statistic": _validate_choice(
            payload,
            "summary_statistic",
            DEFAULT_SUMMARY_STATISTIC,
            SUMMARY_STATISTIC_CHOICES,
        ),
        "stats_backend": _validate_choice(
            payload,
            "stats_backend",
            DEFAULT_STATS_BACKEND,
            STATS_BACKEND_CHOICES,
        ),
        "tree_height_calculation_strategy": _validate_choice(
            payload,
            "tree_height_calculation_strategy",
            DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
            TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
        ),
        "p_value_correction": _validate_choice(
            payload,
            "p_value_correction",
            DEFAULT_P_VALUE_CORRECTION,
            P_VALUE_CORRECTION_CHOICES,
        ),
        "alpha_dct": _validate_optional_float(payload, "alpha_dct", DEFAULT_ALPHA_DCT),
        "alpha_ks": _validate_optional_float(payload, "alpha_ks", DEFAULT_ALPHA_KS),
        "bootstrap": bootstrap,
        "bootstrap_options": bootstrap_options,
    }


def load_orchestrator_config(config_file: str) -> dict:
    """Load and normalize orchestrator config values from JSON/YAML file."""
    payload = _load_raw_config(config_file)
    return normalize_orchestrator_payload(payload)
