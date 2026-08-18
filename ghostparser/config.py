"""Shared configuration trunk for GhostParser.

Holds only the helpers whose behaviour is identical for every consumer.
Module-specific defaults, choices, and validators live with their module.
"""

import json
import shutil
from pathlib import Path


class ConfigError(ValueError):
    """Raised when a config payload is invalid or missing required fields."""


DEFAULT_OVERWRITE = True


def _resolve_path(path_str: str) -> str:
    """Resolve a path string to an absolute path.

    Handles ``~`` expansion, relative paths (resolved from the current working
    directory), and absolute paths (kept as-is).

    Args:
        path_str: The path string to resolve.

    Returns:
        The absolute path as a string.
    """
    return str(Path(path_str).expanduser().resolve())


def _load_raw_config(config_file: str) -> dict:
    """Load a raw config mapping from a JSON or YAML file.

    Args:
        config_file: Path to a ``.json``/``.yaml``/``.yml`` file.

    Returns:
        The parsed config as a dict.

    Raises:
        FileNotFoundError: If the file does not exist.
        ConfigError: If the suffix is unsupported, YAML support is unavailable,
            or the root is not a mapping.
    """
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
    """Validate and resolve a required path field.

    Args:
        payload: The config/CLI payload.
        key: The field name.

    Returns:
        The resolved absolute path.

    Raises:
        ConfigError: If the field is missing or empty.
    """
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Missing required config field: {key}")
    return _resolve_path(value.strip())


def _validate_overwrite_flag(payload: dict, default: bool = DEFAULT_OVERWRITE) -> bool:
    """Resolve the overwrite flag from either ``overwrite`` or ``no_overwrite``.

    The canonical ``overwrite`` key takes precedence over the CLI-style negated
    ``no_overwrite`` key.

    Args:
        payload: The config/CLI payload.
        default: The value used when neither key is present.

    Returns:
        The resolved overwrite boolean.

    Raises:
        ConfigError: If either key is present but not a boolean.
    """
    overwrite = payload.get("overwrite")
    if overwrite is not None:
        if not isinstance(overwrite, bool):
            raise ConfigError("Config field overwrite must be a boolean when provided")
        return overwrite

    no_overwrite = payload.get("no_overwrite")
    if no_overwrite is None:
        return default
    if not isinstance(no_overwrite, bool):
        raise ConfigError("Config field no_overwrite must be a boolean when provided")
    return not no_overwrite


def _next_available_suffixed_path(base_path: Path) -> Path:
    """Return the smallest suffixed path ``<name>_<n>`` that does not exist.

    Uses a single parent-directory scan and computes the smallest missing
    positive suffix in memory.

    Args:
        base_path: The base output path.

    Returns:
        The first non-existing suffixed sibling path.
    """
    parent = base_path.parent
    base_name = base_path.name
    prefix = f"{base_name}_"

    used_suffixes: set[int] = set()
    for entry in parent.iterdir():
        name = entry.name
        if not name.startswith(prefix):
            continue
        raw_suffix = name[len(prefix):]
        if raw_suffix.isdigit():
            used_suffixes.add(int(raw_suffix))

    suffix = 1
    while suffix in used_suffixes:
        suffix += 1

    return base_path.with_name(f"{base_name}_{suffix}")


def prepare_output_directory(
    output_dir: str | Path, *, overwrite: bool = DEFAULT_OVERWRITE
) -> str:
    """Resolve an output directory and either reset it or pick a unique suffix.

    Args:
        output_dir: The requested output directory.
        overwrite: When ``True``, reset an existing directory; when ``False``,
            write to an auto-suffixed sibling if the directory already exists.

    Returns:
        The prepared output directory path as a string.
    """
    path = Path(output_dir).expanduser().resolve()

    if path.exists() and overwrite:
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

    if path.exists() and not overwrite:
        path = _next_available_suffixed_path(path)

    path.mkdir(parents=True, exist_ok=True)
    return str(path)
