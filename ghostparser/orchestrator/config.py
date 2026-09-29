"""The orchestrator's defaults, choices, validation, CLI parser and config
resolution.
"""

import argparse

from ..cli_config import resolve_cli_or_config_args
from ..config import (
    ConfigError,
    _load_raw_config,
    _resolve_path,
    _validate_overwrite_flag,
    _validate_required_path,
    prepare_output_directory,
)
from .trees import load_species_rename_map

__all__ = ["ConfigError", "prepare_output_directory", "build_argument_parser",
           "resolve_config", "load_orchestrator_config", "normalize_orchestrator_payload"]

# Runtime defaults specific to the orchestrator.
DEFAULT_OUTPUT_FOLDER = "results"
DEFAULT_PROCESSES = 0
DEFAULT_MIN_SUPPORT_VALUE = 0.5
DEFAULT_DISCORDANT_TEST = "chi-square"
DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY = "AVG"
DEFAULT_P_VALUE_CORRECTION = "bfn"
DEFAULT_ALPHA_DCT = 0.05
DEFAULT_ALPHA_KS = 0.05
DEFAULT_ALPHA_PERM = 0.05
DEFAULT_PERMUTATION_MIN_RESAMPLES = 2500
DEFAULT_PERMUTATION_MAX_RESAMPLES = 25000
DEFAULT_PERMUTATION_CI_METHOD = "wilson"
DEFAULT_BOOTSTRAP = True
DEFAULT_BOOTSTRAP_ITERATIONS = 100
DEFAULT_BOOTSTRAP_DIAGNOSTIC = False
DEFAULT_BOOTSTRAP_SUMMARY_ONLY = False
DEFAULT_GENERATE_SUMMARY_STATS = False
DEFAULT_SHAPE_DIAGNOSTICS = False
DEFAULT_DIAGNOSTIC = False
DEFAULT_CONSOLIDATION = True
DEFAULT_PREFLIGHT_DATA_CHECK = False
# Caps the triplets the preflight check walks so it stays quick on large
# inputs; 0 lifts the cap. Analysis-only: a real run always processes every
# triplet.
DEFAULT_PREFLIGHT_TRIPLET_CAP = 15000

DISCORDANT_TEST_CHOICES = ("chi-square", "z-test")
TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES = ("AVG", "A", "B", "C", "SIS", "INT")
P_VALUE_CORRECTION_CHOICES = ("no", "bfn", "holm", "fdr_bh", "fdr_by")

# Binomial proportion interval methods accepted by statsmodels'
# ``proportion_confint``. ``wilson`` is the default: it stays inside [0, 1] and
# keeps close-to-nominal coverage for the very small proportions this test
# produces, where the normal approximation degrades badly.
PERMUTATION_CI_METHOD_CHOICES = (
    "wilson",
    "beta",
    "agresti_coull",
    "jeffreys",
    "binom_test",
    "normal",
)

# Parallelization knobs specific to this package.

# CLI argument dest names that also map to config-file payload keys. These are
# the config+CLI options; config-file-only keys (discordant_test,
# tree_height_calculation_strategy, min_support_value, bootstrap_iterations,
# generate_summary_stats, shape_diagnostics, bootstrap_diagnostic,
# bootstrap_summary_only, permutation_options) are intentionally absent so they
# are read only from a config file and otherwise take their defaults.
_ORCHESTRATOR_PAYLOAD_ARG_NAMES = [
    "species_tree_path",
    "gene_trees_path",
    "outgroup",
    "output_folder",
    "triplet_filter",
    "species_filter",
    "species_rename_map",
    "seed",
    "no_overwrite",
    "processes",
    "alpha_dct",
    "alpha_ks",
    "alpha_perm",
    "p_value_correction",
    "diagnostic",
    "consolidation",
    "bootstrap",
    "preflight_data_check",
    "preflight_triplet_cap",
]


def _validate_optional_path(payload: dict, key: str) -> str | None:
    """Validate and resolve an optional path field.

    Args:
        payload: The config/CLI payload.
        key: The field name.

    Returns:
        The resolved absolute path, or ``None`` when the field is absent.

    Raises:
        ConfigError: If the field is present but not a non-empty string.
    """
    value = payload.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(
            f"Config field {key} must be a non-empty string when provided"
        )
    return _resolve_path(value.strip())


def _validate_triplet_selection(payload: dict) -> tuple[str | None, str | None]:
    """Resolve the triplet filter and the species filter, at most one of which may be set.

    Args:
        payload: The config/CLI payload.

    Returns:
        The resolved ``(triplet_filter, species_filter)`` paths, each ``None``
        when absent.

    Raises:
        ConfigError: If either path is not a non-empty string, or both are set.
    """
    triplet_filter = _validate_optional_path(payload, "triplet_filter")
    species_filter = _validate_optional_path(payload, "species_filter")
    if triplet_filter is not None and species_filter is not None:
        raise ConfigError(
            "Config fields triplet_filter and species_filter cannot both be set; "
            "name the triplets or the species, not both"
        )
    return triplet_filter, species_filter


def _validate_species_rename_map(payload: dict) -> str | None:
    """Resolve the optional rename-map path and read the file it names.

    The file is read here, not at run time, so a map that is missing,
    malformed, or gives a taxon a name the outputs cannot carry fails before
    any tree is cleaned rather than after.

    Args:
        payload: The config/CLI payload.

    Returns:
        The resolved absolute path, or ``None`` when the field is absent.

    Raises:
        ConfigError: If the field is present but not a non-empty string, or
            the file cannot be loaded as a rename map.
    """
    path = _validate_optional_path(payload, "species_rename_map")
    if path is None:
        return None
    try:
        load_species_rename_map(path)
    except (FileNotFoundError, ValueError) as exc:
        raise ConfigError(f"Config field species_rename_map: {exc}") from exc
    return path


def _validate_non_negative_int(payload: dict, key: str, default: int) -> int:
    """Validate a non-negative integer field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated integer.

    Raises:
        ConfigError: If the value is not an integer >= 0.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, int) or value < 0:
        raise ConfigError(f"Config field {key} must be an integer >= 0")
    return value


def _validate_optional_int(payload: dict, key: str) -> int | None:
    """Validate an optional integer field.

    Args:
        payload: The config/CLI payload.
        key: The field name.

    Returns:
        The integer, or ``None`` when the field is absent.

    Raises:
        ConfigError: If the field is present but not an integer.
    """
    value = payload.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"Config field {key} must be an integer when provided")
    return value


def _validate_optional_bool(payload: dict, key: str, default: bool) -> bool:
    """Validate an optional boolean field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated boolean.

    Raises:
        ConfigError: If the value is present but not a boolean.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    if not isinstance(value, bool):
        raise ConfigError(f"Config field {key} must be a boolean when provided")
    return value


def _validate_optional_float(payload: dict, key: str, default: float) -> float:
    """Validate an optional numeric field.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.

    Returns:
        The validated float.

    Raises:
        ConfigError: If the value cannot be coerced to a float.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"Config field {key} must be a numeric value") from exc


# YAML 1.1 resolves these bare words to booleans, so a config writing a choice
# such as ``p_value_correction: no`` reaches validation as ``False`` rather than
# ``"no"``. Each bool is mapped back to whichever spelling the field actually
# offers, so the unquoted form works as written.
_YAML_BOOL_WORD_CHOICES = {
    False: ("no", "off", "n", "false"),
    True: ("yes", "on", "y", "true"),
}


def _coerce_yaml_bool_choice(value, choices: tuple[str, ...]):
    """Map a YAML-coerced boolean back to the string choice it was written as.

    Args:
        value: The raw config value.
        choices: The allowed values for the field.

    Returns:
        The matching choice string, or ``value`` unchanged when it is not a
        boolean or the field offers no corresponding spelling.
    """
    if not isinstance(value, bool):
        return value
    for word in _YAML_BOOL_WORD_CHOICES[value]:
        if word in choices:
            return word
    return value


def _validate_choice(
    payload: dict, key: str, default: str, choices: tuple[str, ...]
) -> str:
    """Validate an enumerated string field.

    A boolean value is mapped back to the equivalent string choice first, so
    that YAML's bare-word booleans (``no``, ``yes``, ``on``, ``off``) select the
    choice the user wrote rather than failing validation.

    Args:
        payload: The config/CLI payload.
        key: The field name.
        default: The value used when the field is absent or ``None``.
        choices: The allowed values.

    Returns:
        The validated choice.

    Raises:
        ConfigError: If the value is not one of ``choices``.
    """
    value = payload.get(key, default)
    if value is None:
        value = default
    value = _coerce_yaml_bool_choice(value, choices)
    if not isinstance(value, str) or value not in choices:
        raise ConfigError(
            f"Config field {key} must be one of: {', '.join(choices)} "
            f"(got {value!r})"
        )
    return value


def _validate_bootstrap_options(payload: dict) -> tuple[bool, dict]:
    """Validate the ``bootstrap`` toggle and its options, given nested under
    ``bootstrap_options`` or flat as ``bootstrap_<key>``.

    Args:
        payload: The config/CLI payload.

    Returns:
        A tuple ``(bootstrap, options)`` where ``options`` has
        ``iterations``/``seed``/``diagnostic``/``summary_only``.

    Raises:
        ConfigError: If any value is malformed.
    """
    bootstrap = _validate_optional_bool(payload, "bootstrap", DEFAULT_BOOTSTRAP)

    raw_options = payload.get("bootstrap_options")
    if raw_options is None:
        raw_options = {}
    if not isinstance(raw_options, dict):
        raise ConfigError(
            "Config field bootstrap_options must be a key/value object when provided"
        )

    iterations = payload.get(
        "bootstrap_iterations",
        raw_options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS),
    )
    if iterations is None:
        iterations = DEFAULT_BOOTSTRAP_ITERATIONS
    if not isinstance(iterations, int) or iterations < 1:
        raise ConfigError(
            "Config field bootstrap_options.iterations must be an integer >= 1"
        )

    diagnostic = payload.get(
        "bootstrap_diagnostic",
        raw_options.get("diagnostic", DEFAULT_BOOTSTRAP_DIAGNOSTIC),
    )
    if diagnostic is None:
        diagnostic = DEFAULT_BOOTSTRAP_DIAGNOSTIC
    if not isinstance(diagnostic, bool):
        raise ConfigError(
            "Config field bootstrap_options.diagnostic must be a boolean when provided"
        )

    summary_only = payload.get(
        "bootstrap_summary_only",
        raw_options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY),
    )
    if summary_only is None:
        summary_only = DEFAULT_BOOTSTRAP_SUMMARY_ONLY
    if not isinstance(summary_only, bool):
        raise ConfigError(
            "Config field bootstrap_options.summary_only must be a boolean when provided"
        )

    return bootstrap, {
        "iterations": iterations,
        "diagnostic": diagnostic,
        "summary_only": summary_only,
    }


def _validate_permutation_options(payload: dict) -> dict:
    """Validate the config-file-only ``permutation_options`` block.

    Args:
        payload: The config/CLI payload.

    Returns:
        A dict with ``min_resamples``/``max_resamples``/``ci_method``.

    Raises:
        ConfigError: If any value is malformed or out of range.
    """
    raw_options = payload.get("permutation_options")
    if raw_options is None:
        raw_options = {}
    if not isinstance(raw_options, dict):
        raise ConfigError(
            "Config field permutation_options must be a key/value object when provided"
        )

    min_resamples = raw_options.get("min_resamples", DEFAULT_PERMUTATION_MIN_RESAMPLES)
    if min_resamples is None:
        min_resamples = DEFAULT_PERMUTATION_MIN_RESAMPLES
    if not isinstance(min_resamples, int) or min_resamples < 1:
        raise ConfigError(
            "Config field permutation_options.min_resamples must be an integer >= 1"
        )

    max_resamples = raw_options.get("max_resamples", DEFAULT_PERMUTATION_MAX_RESAMPLES)
    if max_resamples is None:
        max_resamples = DEFAULT_PERMUTATION_MAX_RESAMPLES
    if not isinstance(max_resamples, int) or max_resamples < 1:
        raise ConfigError(
            "Config field permutation_options.max_resamples must be an integer >= 1"
        )
    if max_resamples < min_resamples:
        raise ConfigError(
            "Config field permutation_options.max_resamples must be >= min_resamples"
        )

    ci_method = raw_options.get("ci_method", DEFAULT_PERMUTATION_CI_METHOD)
    if ci_method is None:
        ci_method = DEFAULT_PERMUTATION_CI_METHOD
    if not isinstance(ci_method, str) or ci_method not in PERMUTATION_CI_METHOD_CHOICES:
        raise ConfigError(
            "Config field permutation_options.ci_method must be one of: "
            f"{', '.join(PERMUTATION_CI_METHOD_CHOICES)}"
        )

    return {
        "min_resamples": min_resamples,
        "max_resamples": max_resamples,
        "ci_method": ci_method,
    }


def _parse_outgroup(value) -> list[str]:
    """Parse the ``outgroup`` value into a list of taxon labels.

    One key covers both the single- and multiple-outgroup cases: the value may
    be a single label, a comma-separated string, or a list/tuple/set of labels,
    and always resolves to a list.

    Args:
        value: The raw ``outgroup`` value.

    Returns:
        A non-empty list of outgroup labels.

    Raises:
        ConfigError: If no outgroup labels can be parsed.
    """
    if isinstance(value, str):
        entries = [value]
    elif isinstance(value, (list, tuple, set)):
        entries = [str(entry) for entry in value]
    else:
        entries = []

    # Split every entry on commas, so a list of labels, a comma-separated
    # string, and a list containing comma-separated strings all flatten the
    # same way.
    outgroup = [
        label
        for entry in entries
        for label in (part.strip() for part in entry.split(","))
        if label
    ]

    if not outgroup:
        raise ConfigError("Missing required config field: outgroup")
    return outgroup


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the orchestrator CLI parser.

    Returns:
        A configured ``argparse.ArgumentParser`` for the orchestrator entry point.
    """
    parser = argparse.ArgumentParser(
        prog="python -m ghostparser.orchestrator",
        description=(
            "GhostParser orchestrator: fuse triplet extraction and inference "
            "in a single pass without materializing the intermediate triplet file."
        ),
    )
    parser.add_argument(
        "-c",
        "--config-file",
        default=None,
        help="Path to a JSON or YAML config file; any other flag given alongside it overrides the file's value",
    )
    parser.add_argument(
        "-st",
        "--species-tree-path",
        default=None,
        help="Path to the species tree file in Newick format",
    )
    parser.add_argument(
        "-gt",
        "--gene-trees-path",
        default=None,
        help="Path to the gene trees file in Newick format",
    )
    parser.add_argument(
        "-og",
        "--outgroup",
        default=None,
        help="Outgroup species identifier(s), comma-separated for more than one",
    )
    parser.add_argument(
        "--output-folder",
        default=None,
        help=f"Output folder (default: ./{DEFAULT_OUTPUT_FOLDER})",
    )
    parser.add_argument(
        "--triplet-filter",
        default=None,
        help="Path to triplet filter file (comma-separated taxa per line)",
    )
    parser.add_argument(
        "--species-filter",
        default=None,
        help="Path to a species filter file (comma-separated taxa, any number "
        "per line); every triplet among the listed species is analyzed. "
        "Cannot be combined with --triplet-filter",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base RNG seed for the whole run; every random draw derives from "
        "it. Omit for a fresh seed each run (the value used is reported in "
        "metrics.txt either way)",
    )
    parser.add_argument(
        "--species-rename-map",
        default=None,
        help="Path to a species rename map (two-column TSV, or YAML mapping) "
        "giving the name each taxon should appear under in the outputs",
    )
    parser.add_argument(
        "--no-overwrite",
        dest="no_overwrite",
        action="store_true",
        default=None,
        help="Append a numeric suffix when the output folder already exists",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes (0 = all cores)",
    )
    parser.add_argument(
        "--alpha-dct",
        type=float,
        default=None,
        help=f"DCT significance threshold (default: {DEFAULT_ALPHA_DCT})",
    )
    parser.add_argument(
        "--alpha-ks",
        type=float,
        default=None,
        help=f"KS significance threshold (default: {DEFAULT_ALPHA_KS})",
    )
    parser.add_argument(
        "--alpha-perm",
        type=float,
        default=None,
        help=f"Permutation test significance threshold (default: {DEFAULT_ALPHA_PERM})",
    )
    parser.add_argument(
        "--p-value-correction",
        choices=P_VALUE_CORRECTION_CHOICES,
        default=None,
        help=f"Multiple-testing correction for triplet p-values (default: {DEFAULT_P_VALUE_CORRECTION})",
    )
    parser.add_argument(
        "--diagnostic",
        dest="diagnostic",
        action="store_true",
        default=None,
        help=(
            "Measure every test for every triplet so every ks_* and perm_* "
            "column is filled; by default a test the decision cascade cannot "
            "consult is skipped (default: disabled)"
        ),
    )
    parser.add_argument(
        "--no-consolidation",
        dest="consolidation",
        action="store_false",
        default=None,
        help="Disable introgression map consolidation outputs (default: enabled)",
    )
    parser.add_argument(
        "--no-bootstrap",
        dest="bootstrap",
        action="store_false",
        default=None,
        help="Disable bootstrap sampling-with-replacement during triplet inference",
    )
    parser.add_argument(
        "--preflight-data-check",
        dest="preflight_data_check",
        action="store_true",
        default=None,
        help=(
            "Run only the structural preflight check on the input trees, write "
            "preflight_data_check.txt into the output folder, and exit without "
            "running any analysis (default: disabled)"
        ),
    )
    parser.add_argument(
        "--preflight-triplet-cap",
        dest="preflight_triplet_cap",
        type=int,
        default=None,
        help=(
            "Maximum number of ingroup triplets the preflight data check walks; "
            "0 checks every triplet (default: "
            f"{DEFAULT_PREFLIGHT_TRIPLET_CAP})"
        ),
    )
    return parser


def normalize_orchestrator_payload(payload: dict) -> dict:
    """Validate a config-file or CLI payload and fill every runtime key,
    applying the defaults for anything omitted.

    Args:
        payload: Raw config/CLI key-value mapping.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.orchestrator.runner.run_orchestrator`.

    Raises:
        ConfigError: If a required field is missing or a value is invalid.
    """
    species_tree = _validate_required_path(payload, "species_tree_path")
    gene_trees = _validate_required_path(payload, "gene_trees_path")

    outgroup = _parse_outgroup(payload.get("outgroup"))

    output = _validate_optional_path(payload, "output_folder")
    if output is None:
        output = _resolve_path(DEFAULT_OUTPUT_FOLDER)

    bootstrap, bootstrap_options = _validate_bootstrap_options(payload)
    permutation_options = _validate_permutation_options(payload)
    triplet_filter, species_filter = _validate_triplet_selection(payload)
    return {
        "species_tree": species_tree,
        "gene_trees": gene_trees,
        "outgroup": outgroup,
        "triplet_filter": triplet_filter,
        "species_filter": species_filter,
        "species_rename_map": _validate_species_rename_map(payload),
        "output": output,
        "overwrite": _validate_overwrite_flag(payload),
        "processes": _validate_non_negative_int(
            payload, "processes", DEFAULT_PROCESSES
        ),
        "seed": _validate_optional_int(payload, "seed"),
        "consolidation": _validate_optional_bool(
            payload, "consolidation", DEFAULT_CONSOLIDATION
        ),
        "preflight_data_check": _validate_optional_bool(
            payload, "preflight_data_check", DEFAULT_PREFLIGHT_DATA_CHECK
        ),
        "preflight_triplet_cap": _validate_non_negative_int(
            payload, "preflight_triplet_cap", DEFAULT_PREFLIGHT_TRIPLET_CAP
        ),
        "shape_diagnostics": _validate_optional_bool(
            payload, "shape_diagnostics", DEFAULT_SHAPE_DIAGNOSTICS
        ),
        "diagnostic": _validate_optional_bool(
            payload, "diagnostic", DEFAULT_DIAGNOSTIC
        ),
        "generate_summary_stats": _validate_optional_bool(
            payload, "generate_summary_stats", DEFAULT_GENERATE_SUMMARY_STATS
        ),
        "min_support_value": _validate_optional_float(
            payload, "min_support_value", DEFAULT_MIN_SUPPORT_VALUE
        ),
        "discordant_test": _validate_choice(
            payload, "discordant_test", DEFAULT_DISCORDANT_TEST, DISCORDANT_TEST_CHOICES
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
        "alpha_perm": _validate_optional_float(
            payload, "alpha_perm", DEFAULT_ALPHA_PERM
        ),
        "permutation_min_resamples": permutation_options["min_resamples"],
        "permutation_max_resamples": permutation_options["max_resamples"],
        "permutation_ci_method": permutation_options["ci_method"],
        "bootstrap": bootstrap,
        "bootstrap_iterations": bootstrap_options["iterations"],
        "bootstrap_diagnostic": bootstrap_options["diagnostic"],
        "bootstrap_summary_only": bootstrap_options["summary_only"],
    }


def load_orchestrator_config(config_file: str) -> dict:
    """Load and normalize an orchestrator config from a JSON/YAML file.

    Args:
        config_file: Path to the JSON or YAML config file.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.orchestrator.runner.run_orchestrator`.

    Raises:
        ConfigError: If the file is malformed or a value is invalid.
    """
    payload = _load_raw_config(config_file)
    return normalize_orchestrator_payload(payload)


def resolve_config(args: argparse.Namespace) -> dict:
    """Resolve a parsed CLI namespace into the runtime config, laying any given
    flags over the config file when one is named.

    Args:
        args: Parsed ``argparse.Namespace`` from :func:`build_argument_parser`.

    Returns:
        The resolved config dict consumed by
        :func:`ghostparser.orchestrator.runner.run_orchestrator`.

    Raises:
        ConfigError: If a required argument is missing or a value is invalid.
    """
    resolved = resolve_cli_or_config_args(
        args,
        normalize_payload=normalize_orchestrator_payload,
        payload_arg_names=_ORCHESTRATOR_PAYLOAD_ARG_NAMES,
    )
    return vars(resolved)
