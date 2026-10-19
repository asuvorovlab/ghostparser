"""Load a simulation config (YAML or JSON), fill defaults and validate it."""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import yaml


class SimError(Exception):
    """A problem the user can fix; reported as one line with an exit status."""

    exit_status = 1


class SimConfigError(SimError):
    """The config is malformed or inconsistent."""

    exit_status = 2


class SimInputError(SimError):
    """An input or the environment cannot be used (a tree file, a missing package)."""

    exit_status = 1


TIME_UNITS = ("generations", "coalescent", "years")
TREE_MODELS = ("newick", "yule", "birth_death")
TREE_MODEL_ALIASES = {"pure_birth": "yule"}
TREE_CHOICES = ("middle", "first")
ANCESTRY_MODELS = ("hudson", "smc", "smc_prime", "dtwf")
CLOCK_MODELS = ("ucln", "strict")
OUTPUT_UNITS = ("normalized", "substitutions", "time_units", "generations", "coalescent")
# Units whose lengths are substitutions, so the clock applies to them.
CLOCK_UNITS = ("normalized", "substitutions")
PLOT_FORMATS = ("svg", "png", "pdf")
BRANCH_KINDS = ("tips", "ancestral", "any")

_UNSAFE_LABEL = re.compile(r"[\s()\[\]:;,'\"{}]")
# YAML 1.1 reads 1e-8 (no decimal point) as a string; numbers like it are coerced.
_SCIENTIFIC = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)[eE][-+]?\d+$")

DEFAULTS: dict[str, Any] = {
    "seed": None,
    "output_dir": "./sim_out",
    "overwrite": False,
    "replicates": 1,
    "time_units": "generations",
    "generation_time": None,
    "species": None,
    "outgroup": None,
    "ghosts": [],
    "extinct": [],
    "species_tree": {
        "model": None,
        "newick": None,
        "newick_file": None,
        "root_age": None,
        "birth_rate": 1.0,
        "death_rate": None,
        "extinction_fraction": [0.2, 0.8],
        "extinct_from_process": False,
        "outgroup_age": 1.5,
        "max_attempts": 1000,
    },
    "population": {"size": 10000, "sizes": [], "ploidy": 2},
    "introgression": [],
    "random_introgression": [],
    "gene_flow": [],
    "size_changes": [],
    "bottlenecks": [],
    "gene_tree_count": 1000,
    "locus": {"sequence_length": 1, "recombination_rate": 0.0, "tree_choice": "middle"},
    "ancestry_model": "hudson",
    "clock": {"model": "ucln", "mean_rate": 1e-8, "sigma": 0.3, "among_loci_gamma_shape": None},
    "output": {
        "branch_length_units": "normalized",
        "scale": 1.0,
        "decimals": 6,
        "write_ghostparser_config": True,
        "plot": "svg",
    },
}

_EVENT_KEYS = {
    "introgression": {"donor", "recipient", "proportion", "time", "time_fraction"},
    "random_introgression": {
        "count",
        "donors",
        "recipients",
        "proportion",
        "time_fraction",
        "allow_sisters",
        "donor_branches",
        "recipient_branches",
    },
    "gene_flow": {"donor", "recipient", "rate", "start", "end", "symmetric"},
    "size_changes": {"branch", "time", "time_fraction", "size", "growth_rate"},
    "bottlenecks": {"branch", "time", "time_fraction", "strength", "proportion"},
}

_BD_ONLY_KEYS = ("death_rate", "extinct_from_process", "max_attempts")
_RANDOM_TREE_KEYS = ("root_age", "birth_rate", "outgroup_age") + _BD_ONLY_KEYS
_NEWICK_KEYS = ("newick", "newick_file")


@dataclass(frozen=True)
class SpeciesSet:
    """Every species name with its role, plus the prefix groups that name them.

    Attributes:
        names: All species, in declaration order.
        groups: Prefix to the names it expanded to (empty for a plain list).
        outgroup: Sampled outgroup species.
        ghosts: Species alive at present but never sampled.
        extinct: Species that end before the present.
    """

    names: tuple[str, ...]
    groups: dict[str, tuple[str, ...]]
    outgroup: tuple[str, ...]
    ghosts: tuple[str, ...]
    extinct: tuple[str, ...]

    @property
    def ingroup(self) -> tuple[str, ...]:
        return tuple(n for n in self.names if n not in self.outgroup)

    @property
    def sampled(self) -> tuple[str, ...]:
        hidden = set(self.ghosts) | set(self.extinct)
        return tuple(n for n in self.names if n not in hidden)

    def role(self, name: str) -> str:
        if name in self.outgroup:
            return "outgroup"
        if name in self.ghosts:
            return "ghost"
        if name in self.extinct:
            return "extinct"
        return "sampled"

    def resolve(self, tokens: Any, key: str) -> tuple[str, ...]:
        """Expand names and prefixes into species names, in declaration order.

        Args:
            tokens: One name or prefix, or a list of them.
            key: The config key being resolved, for the error message.

        Returns:
            The species the tokens name, without duplicates.

        Raises:
            SimConfigError: A token is neither a species name nor a prefix.
        """
        if tokens is None:
            return ()
        if isinstance(tokens, str):
            tokens = [tokens]
        if not isinstance(tokens, list) or not all(isinstance(t, str) for t in tokens):
            raise SimConfigError(f"{key} must be a name, a prefix, or a list of them")
        chosen: set[str] = set()
        for token in tokens:
            if token in self.names:
                chosen.add(token)
            elif token in self.groups:
                chosen.update(self.groups[token])
            else:
                known = ", ".join(list(self.groups) + list(self.names))
                raise SimConfigError(f"{key}: '{token}' is neither a species nor a prefix (known: {known})")
        return tuple(n for n in self.names if n in chosen)


@dataclass(frozen=True)
class SimConfig:
    """A validated config. ``raw`` is the defaults-filled mapping, kept for the record.

    Times in ``raw`` are in ``time_units``; ``generations_per_unit`` converts them.
    """

    raw: dict[str, Any]
    species: SpeciesSet
    newick: str | None
    output_dir: Path

    def __getitem__(self, key: str) -> Any:
        return self.raw[key]

    @property
    def tree(self) -> dict[str, Any]:
        return self.raw["species_tree"]

    @property
    def generations_per_unit(self) -> float:
        units = self.raw["time_units"]
        if units == "generations":
            return 1.0
        if units == "coalescent":
            return float(self.raw["population"]["ploidy"] * self.raw["population"]["size"])
        return 1.0 / self.raw["generation_time"]


def load_config(path: str | Path, overrides: dict[str, Any] | None = None) -> SimConfig:
    """Read a config file, lay any overrides over it, fill defaults and validate.

    Args:
        path: YAML or JSON config file.
        overrides: Top-level keys that replace the file's values (from the command line).

    Returns:
        The validated config.

    Raises:
        SimConfigError: The config is malformed or inconsistent.
        SimInputError: The config or the Newick file it names cannot be read.
    """
    path = Path(path).expanduser()
    try:
        payload = yaml.safe_load(path.read_text())
    except OSError as exc:
        raise SimInputError(f"Cannot read config {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise SimConfigError(f"Cannot parse config {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SimConfigError(f"Config {path} must be a mapping of keys to values")
    payload.update({k: v for k, v in (overrides or {}).items() if v is not None})
    raw = _merge_defaults(_coerce_scientific(payload))
    _validate_scalars(raw)
    newick = _read_newick(raw)
    species = _build_species_set(raw, newick)
    _validate_events(raw)
    raw["species"] = list(species.names)
    raw["outgroup"] = list(species.outgroup)
    raw["ghosts"] = list(species.ghosts)
    raw["extinct"] = list(species.extinct)
    if newick is not None:
        raw["species_tree"]["newick"] = newick
        raw["species_tree"]["newick_file"] = None
    output_dir = Path(raw["output_dir"]).expanduser().resolve()
    return SimConfig(raw=raw, species=species, newick=newick, output_dir=output_dir)


def _coerce_scientific(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _coerce_scientific(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_coerce_scientific(v) for v in value]
    if isinstance(value, str) and _SCIENTIFIC.match(value):
        return float(value)
    return value


def _merge_defaults(payload: dict[str, Any]) -> dict[str, Any]:
    raw = copy.deepcopy(DEFAULTS)
    _reject_unknown(payload, DEFAULTS, "config")
    for key, value in payload.items():
        default = DEFAULTS[key]
        if isinstance(default, dict):
            if not isinstance(value, dict):
                raise SimConfigError(f"{key} must be a mapping")
            _reject_unknown(value, default, key)
            raw[key].update(value)
        else:
            raw[key] = value
    tree = raw["species_tree"]
    model = tree["model"]
    if model is None:
        raise SimConfigError("species_tree.model is required (newick, yule or birth_death)")
    tree["model"] = TREE_MODEL_ALIASES.get(model, model)
    if tree["model"] not in TREE_MODELS:
        raise SimConfigError(f"species_tree.model must be one of {', '.join(TREE_MODELS)}, got '{model}'")
    # Keys left at their default (as in a resolved_config.yaml) count as not given.
    given = {k for k, v in payload.get("species_tree", {}).items() if v != DEFAULTS["species_tree"][k]}
    if tree["model"] == "newick":
        stray = given.intersection(_RANDOM_TREE_KEYS)
    elif tree["model"] == "yule":
        stray = given.intersection(_NEWICK_KEYS + _BD_ONLY_KEYS)
    else:
        stray = given.intersection(_NEWICK_KEYS)
        if tree["death_rate"] is None:
            tree["death_rate"] = 0.5
    if stray:
        raise SimConfigError(
            f"species_tree.{sorted(stray)[0]} does not apply to model '{tree['model']}'"
        )
    if tree["model"] != "birth_death":
        tree["death_rate"] = None
        tree["extinct_from_process"] = False
    return raw


def _reject_unknown(given: dict[str, Any], allowed: dict[str, Any] | set[str], where: str) -> None:
    unknown = sorted(set(given) - set(allowed))
    if unknown:
        raise SimConfigError(f"Unknown key '{unknown[0]}' in {where}")


def _number(value: Any, key: str, *, low: float | None = None, high: float | None = None,
            low_open: bool = False, integer: bool = False) -> float:
    kind = int if integer else (int, float)
    if isinstance(value, bool) or not isinstance(value, kind):
        raise SimConfigError(f"{key} must be {'an integer' if integer else 'a number'}, got {value!r}")
    if low is not None and (value < low or (low_open and value == low)):
        raise SimConfigError(f"{key} must be {'>' if low_open else '>='} {low}, got {value}")
    if high is not None and value > high:
        raise SimConfigError(f"{key} must be <= {high}, got {value}")
    return value


def _range(value: Any, key: str, *, low: float, high: float, low_open: bool = True) -> tuple[float, float]:
    """Read a scalar or a ``[lo, hi]`` pair inside ``(low, high]`` as a range."""
    pair = value if isinstance(value, list) else [value, value]
    if len(pair) != 2:
        raise SimConfigError(f"{key} must be a number or a [low, high] pair")
    lo, hi = (_number(v, key, low=low, high=high, low_open=low_open) for v in pair)
    if lo > hi:
        raise SimConfigError(f"{key}: low {lo} exceeds high {hi}")
    return lo, hi


def _bool(value: Any, key: str) -> bool:
    if not isinstance(value, bool):
        raise SimConfigError(f"{key} must be true or false, got {value!r}")
    return value


def _validate_scalars(raw: dict[str, Any]) -> None:
    if raw["seed"] is not None:
        _number(raw["seed"], "seed", low=0, integer=True)
    _bool(raw["overwrite"], "overwrite")
    _number(raw["replicates"], "replicates", low=1, integer=True)
    if raw["time_units"] not in TIME_UNITS:
        raise SimConfigError(f"time_units must be one of {', '.join(TIME_UNITS)}")
    if raw["time_units"] == "years":
        if raw["generation_time"] is None:
            raise SimConfigError("time_units: years needs generation_time (years per generation)")
        _number(raw["generation_time"], "generation_time", low=0, low_open=True)

    tree = raw["species_tree"]
    if tree["model"] != "newick":
        if tree["root_age"] is None:
            raise SimConfigError(f"species_tree.root_age is required for model '{tree['model']}'")
        _number(tree["root_age"], "species_tree.root_age", low=0, low_open=True)
        _number(tree["birth_rate"], "species_tree.birth_rate", low=0, low_open=True)
        _number(tree["outgroup_age"], "species_tree.outgroup_age", low=1, low_open=True)
    if tree["model"] == "birth_death":
        _number(tree["death_rate"], "species_tree.death_rate", low=0)
        if tree["death_rate"] >= tree["birth_rate"]:
            raise SimConfigError("species_tree.death_rate must be below birth_rate")
        _bool(tree["extinct_from_process"], "species_tree.extinct_from_process")
        _number(tree["max_attempts"], "species_tree.max_attempts", low=1, integer=True)
        if tree["extinct_from_process"] and tree["death_rate"] == 0:
            raise SimConfigError("species_tree.extinct_from_process needs death_rate > 0")
    tree["extinction_fraction"] = list(
        _range(tree["extinction_fraction"], "species_tree.extinction_fraction", low=0, high=1)
    )
    if tree["extinction_fraction"][1] >= 1:
        raise SimConfigError("species_tree.extinction_fraction must stay below 1")

    pop = raw["population"]
    _number(pop["size"], "population.size", low=0, low_open=True)
    _number(pop["ploidy"], "population.ploidy", low=1, integer=True)
    if not isinstance(pop["sizes"], list):
        raise SimConfigError("population.sizes must be a list of {branch, size}")
    for i, entry in enumerate(pop["sizes"]):
        where = f"population.sizes[{i}]"
        if not isinstance(entry, dict) or set(entry) != {"branch", "size"}:
            raise SimConfigError(f"{where} must have exactly the keys branch and size")
        _number(entry["size"], f"{where}.size", low=0, low_open=True)

    _number(raw["gene_tree_count"], "gene_tree_count", low=1, integer=True)
    locus = raw["locus"]
    _number(locus["sequence_length"], "locus.sequence_length", low=1)
    _number(locus["recombination_rate"], "locus.recombination_rate", low=0)
    if locus["tree_choice"] not in TREE_CHOICES:
        raise SimConfigError(f"locus.tree_choice must be one of {', '.join(TREE_CHOICES)}")
    _validate_ancestry_model(raw["ancestry_model"])

    clock = raw["clock"]
    if clock["model"] not in CLOCK_MODELS:
        raise SimConfigError(f"clock.model must be one of {', '.join(CLOCK_MODELS)}")
    _number(clock["mean_rate"], "clock.mean_rate", low=0, low_open=True)
    _number(clock["sigma"], "clock.sigma", low=0)
    if clock["among_loci_gamma_shape"] is not None:
        _number(clock["among_loci_gamma_shape"], "clock.among_loci_gamma_shape", low=0, low_open=True)

    out = raw["output"]
    if out["branch_length_units"] not in OUTPUT_UNITS:
        raise SimConfigError(f"output.branch_length_units must be one of {', '.join(OUTPUT_UNITS)}")
    _number(out["scale"], "output.scale", low=0, low_open=True)
    _number(out["decimals"], "output.decimals", low=0, integer=True)
    _bool(out["write_ghostparser_config"], "output.write_ghostparser_config")
    if out["plot"] is not False and out["plot"] not in PLOT_FORMATS:
        raise SimConfigError(f"output.plot must be one of {', '.join(PLOT_FORMATS)} or false")


def _validate_ancestry_model(model: Any) -> None:
    if isinstance(model, str):
        if model not in ANCESTRY_MODELS:
            raise SimConfigError(
                f"ancestry_model must be one of {', '.join(ANCESTRY_MODELS)}, beta or dirac, got '{model}'"
            )
        return
    if not isinstance(model, dict) or len(model) != 1:
        raise SimConfigError("ancestry_model must be a name, {beta: {alpha}} or {dirac: {psi, c}}")
    (name, params), = model.items()
    params = params or {}
    if name == "beta":
        _reject_unknown(params, {"alpha", "truncation_point"}, "ancestry_model.beta")
        if "alpha" not in params:
            raise SimConfigError("ancestry_model.beta needs alpha")
        _number(params["alpha"], "ancestry_model.beta.alpha", low=1, high=2, low_open=True)
        if "truncation_point" in params:
            _number(params["truncation_point"], "ancestry_model.beta.truncation_point", low=0, low_open=True)
    elif name == "dirac":
        _reject_unknown(params, {"psi", "c"}, "ancestry_model.dirac")
        if not {"psi", "c"} <= set(params):
            raise SimConfigError("ancestry_model.dirac needs psi and c")
        _number(params["psi"], "ancestry_model.dirac.psi", low=0, high=1, low_open=True)
        _number(params["c"], "ancestry_model.dirac.c", low=0)
    else:
        raise SimConfigError(f"ancestry_model: unknown model '{name}'")


def _read_newick(raw: dict[str, Any]) -> str | None:
    tree = raw["species_tree"]
    if tree["model"] != "newick":
        return None
    if (tree["newick"] is None) == (tree["newick_file"] is None):
        raise SimConfigError("model 'newick' needs exactly one of species_tree.newick or newick_file")
    if tree["newick"] is not None:
        return str(tree["newick"]).strip()
    path = Path(tree["newick_file"]).expanduser()
    try:
        return path.read_text().strip()
    except OSError as exc:
        raise SimInputError(f"Cannot read species_tree.newick_file {path}: {exc}") from exc


def newick_tip_names(newick: str) -> list[str]:
    """Return the tip labels of a Newick string, in file order.

    Raises:
        SimInputError: The string is not one parseable Newick tree.
    """
    from Bio import Phylo

    try:
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception as exc:  # Bio.Phylo raises several unrelated types
        raise SimInputError(f"Invalid Newick species tree: {exc}") from exc
    return [str(tip.name) if tip.name is not None else "" for tip in tree.get_terminals()]


def _expand_species(spec: Any) -> tuple[list[str], dict[str, tuple[str, ...]], dict[str, set[str]]]:
    """Expand ``species`` into names, prefix groups and inline roles."""
    roles: dict[str, set[str]] = {"outgroup": set(), "ghost": set(), "extinct": set()}
    if isinstance(spec, list):
        if not all(isinstance(n, str) for n in spec):
            raise SimConfigError("species list entries must be names")
        return list(spec), {}, roles
    if not isinstance(spec, dict):
        raise SimConfigError("species must be a list of names or a {prefix: count} mapping")
    names: list[str] = []
    groups: dict[str, tuple[str, ...]] = {}
    for prefix, value in spec.items():
        prefix = str(prefix)
        options = value if isinstance(value, dict) else {"count": value}
        _reject_unknown(options, {"count", "outgroup", "ghost", "extinct"}, f"species.{prefix}")
        if "count" not in options:
            raise SimConfigError(f"species.{prefix} needs a count")
        count = _number(options["count"], f"species.{prefix}.count", low=1, integer=True)
        group = tuple(f"{prefix}{i}" for i in range(count))
        groups[prefix] = group
        names.extend(group)
        flags = [r for r in ("outgroup", "ghost", "extinct") if _bool(options.get(r, False), f"species.{prefix}.{r}")]
        if len(flags) > 1:
            raise SimConfigError(f"species.{prefix} can take only one of outgroup, ghost, extinct")
        for role in flags:
            roles[role].update(group)
    return names, groups, roles


def _build_species_set(raw: dict[str, Any], newick: str | None) -> SpeciesSet:
    if raw["species"] is None:
        if newick is None:
            raise SimConfigError("species is required unless species_tree.model is newick")
        names, groups, roles = newick_tip_names(newick), {}, {"outgroup": set(), "ghost": set(), "extinct": set()}
    else:
        names, groups, roles = _expand_species(raw["species"])
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise SimConfigError(f"Duplicate species name '{duplicates[0]}'")
    for name in names:
        if not name or _UNSAFE_LABEL.search(name):
            raise SimConfigError(f"Species name '{name}' is empty or holds a character Newick cannot carry")
    if newick is not None and raw["species"] is not None:
        tips = newick_tip_names(newick)
        if sorted(tips) != sorted(names):
            missing = sorted(set(names) ^ set(tips))
            raise SimConfigError(f"species and the Newick tips differ, first mismatch '{missing[0]}'")

    probe = SpeciesSet(tuple(names), groups, (), (), ())
    outgroup = set(probe.resolve(raw["outgroup"], "outgroup")) | roles["outgroup"]
    ghosts = set(probe.resolve(raw["ghosts"], "ghosts")) | roles["ghost"]
    extinct = set(probe.resolve(raw["extinct"], "extinct")) | roles["extinct"]
    if not outgroup:
        raise SimConfigError("outgroup is required (a name or prefix, or outgroup: true on a species group)")
    for a, b, label in ((outgroup, ghosts, "ghost"), (outgroup, extinct, "extinct"), (ghosts, extinct, "extinct")):
        both = sorted(a & b)
        if both:
            first = "outgroup" if a is outgroup else "ghost"
            raise SimConfigError(f"'{both[0]}' is both {first} and {label}")
    ordered = lambda chosen: tuple(n for n in names if n in chosen)  # noqa: E731
    species = SpeciesSet(tuple(names), groups, ordered(outgroup), ordered(ghosts), ordered(extinct))
    sampled_ingroup = [n for n in species.sampled if n not in outgroup]
    if len(sampled_ingroup) < 3:
        raise SimConfigError(f"At least 3 sampled ingroup species are needed, got {len(sampled_ingroup)}")
    if len(species.ingroup) < 2:
        raise SimConfigError("The ingroup needs at least 2 species")
    return species


def _check_branch(value: Any, where: str) -> None:
    ok = isinstance(value, str) or (
        isinstance(value, list) and value and all(isinstance(v, str) for v in value)
    )
    if not ok:
        raise SimConfigError(f"{where} must be a species name, a prefix, or a list of them")


def _validate_events(raw: dict[str, Any]) -> None:
    for section, allowed in _EVENT_KEYS.items():
        entries = raw[section]
        if not isinstance(entries, list):
            raise SimConfigError(f"{section} must be a list")
        for i, entry in enumerate(entries):
            where = f"{section}[{i}]"
            if not isinstance(entry, dict):
                raise SimConfigError(f"{where} must be a mapping")
            _reject_unknown(entry, allowed, where)
            _VALIDATORS[section](entry, where)


def _need(entry: dict[str, Any], keys: tuple[str, ...], where: str) -> None:
    for key in keys:
        if key not in entry:
            raise SimConfigError(f"{where} needs {key}")


def _one_of(entry: dict[str, Any], keys: tuple[str, str], where: str) -> str:
    present = [k for k in keys if k in entry]
    if len(present) != 1:
        raise SimConfigError(f"{where} needs exactly one of {keys[0]} or {keys[1]}")
    return present[0]


def _check_introgression(entry: dict[str, Any], where: str) -> None:
    _need(entry, ("donor", "recipient", "proportion"), where)
    _check_branch(entry["donor"], f"{where}.donor")
    _check_branch(entry["recipient"], f"{where}.recipient")
    _number(entry["proportion"], f"{where}.proportion", low=0, high=1, low_open=True)
    _check_time(entry, where)


def _check_time(entry: dict[str, Any], where: str) -> None:
    """Exactly one of ``time`` (> 0) or ``time_fraction`` (strictly between 0 and 1)."""
    if _one_of(entry, ("time", "time_fraction"), where) == "time":
        _number(entry["time"], f"{where}.time", low=0, low_open=True)
    else:
        _number(entry["time_fraction"], f"{where}.time_fraction", low=0, high=1, low_open=True)
        if entry["time_fraction"] >= 1:
            raise SimConfigError(f"{where}.time_fraction must be below 1")


def _check_random_introgression(entry: dict[str, Any], where: str) -> None:
    _need(entry, ("count", "donors", "recipients", "proportion"), where)
    _number(entry["count"], f"{where}.count", low=1, integer=True)
    _check_branch(entry["donors"], f"{where}.donors")
    _check_branch(entry["recipients"], f"{where}.recipients")
    entry["proportion"] = list(_range(entry["proportion"], f"{where}.proportion", low=0, high=1))
    fraction = list(_range(entry.get("time_fraction", [0.1, 0.9]), f"{where}.time_fraction", low=0, high=1))
    if fraction[1] >= 1:
        raise SimConfigError(f"{where}.time_fraction must stay below 1")
    entry["time_fraction"] = fraction
    entry["allow_sisters"] = _bool(entry.get("allow_sisters", False), f"{where}.allow_sisters")
    for key in ("donor_branches", "recipient_branches"):
        entry[key] = entry.get(key, "tips")
        if entry[key] not in BRANCH_KINDS:
            raise SimConfigError(f"{where}.{key} must be one of {', '.join(BRANCH_KINDS)}")


def _check_gene_flow(entry: dict[str, Any], where: str) -> None:
    _need(entry, ("donor", "recipient", "rate"), where)
    _check_branch(entry["donor"], f"{where}.donor")
    _check_branch(entry["recipient"], f"{where}.recipient")
    _number(entry["rate"], f"{where}.rate", low=0, high=1, low_open=True)
    for key in ("start", "end"):
        if entry.get(key) is not None:
            _number(entry[key], f"{where}.{key}", low=0)
    if entry.get("start") is not None and entry.get("end") is not None and entry["start"] >= entry["end"]:
        raise SimConfigError(f"{where}.start must be below end (times run backwards from the present)")
    entry["symmetric"] = _bool(entry.get("symmetric", False), f"{where}.symmetric")


def _check_size_change(entry: dict[str, Any], where: str) -> None:
    _need(entry, ("branch",), where)
    _check_branch(entry["branch"], f"{where}.branch")
    _check_time(entry, where)
    if "size" not in entry and "growth_rate" not in entry:
        raise SimConfigError(f"{where} needs size, growth_rate or both")
    if "size" in entry:
        _number(entry["size"], f"{where}.size", low=0, low_open=True)
    if "growth_rate" in entry:
        _number(entry["growth_rate"], f"{where}.growth_rate")


def _check_bottleneck(entry: dict[str, Any], where: str) -> None:
    _need(entry, ("branch",), where)
    _check_branch(entry["branch"], f"{where}.branch")
    _check_time(entry, where)
    if _one_of(entry, ("strength", "proportion"), where) == "strength":
        _number(entry["strength"], f"{where}.strength", low=0, low_open=True)
    else:
        _number(entry["proportion"], f"{where}.proportion", low=0, high=1, low_open=True)


_VALIDATORS = {
    "introgression": _check_introgression,
    "random_introgression": _check_random_introgression,
    "gene_flow": _check_gene_flow,
    "size_changes": _check_size_change,
    "bottlenecks": _check_bottleneck,
}
