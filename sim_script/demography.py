"""Turn a species tree and the configured events into an msprime demography.

Events are written forward in time (genes move from ``donor`` into
``recipient``); msprime runs backwards, so a pulse moves the recipient's
lineages into the donor. Every time here is in generations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sim_config import SimConfig, SimConfigError
from species_tree import Node, find_tip, mrca

EVENT_COLUMNS = (
    "kind", "origin", "donor", "recipient", "branch", "time", "end_time",
    "time_generations", "end_time_generations", "proportion", "rate", "size",
    "growth_rate", "strength",
)


@dataclass
class Event:
    """One realized event. ``time``/``end`` are generations before the present."""

    kind: str  # pulse, gene_flow, size_change, bottleneck, extinction
    origin: str  # config, random or tree
    time: float
    end: float | None = None
    donor: Node | None = None
    recipient: Node | None = None
    branch: Node | None = None
    proportion: float | None = None
    rate: float | None = None
    size: float | None = None
    growth_rate: float | None = None
    strength: float | None = None

    def row(self, per_unit: float) -> dict[str, str]:
        """The event as an ``events.tsv`` row, times also given in ``time_units``."""
        def fmt(value: float | None) -> str:
            return "" if value is None else f"{value:.6g}"

        def lab(node: Node | None) -> str:
            return "" if node is None else node.label()

        return {
            "kind": self.kind, "origin": self.origin, "donor": lab(self.donor),
            "recipient": lab(self.recipient), "branch": lab(self.branch),
            "time": fmt(self.time / per_unit),
            "end_time": fmt(None if self.end is None else self.end / per_unit),
            "time_generations": fmt(self.time), "end_time_generations": fmt(self.end),
            "proportion": fmt(self.proportion), "rate": fmt(self.rate), "size": fmt(self.size),
            "growth_rate": fmt(self.growth_rate), "strength": fmt(self.strength),
        }


def resolve_branch(root: Node, cfg: SimConfig, token: str | list[str], where: str) -> Node:
    """The branch a reference names: a tip, or the stem above the MRCA of the named tips.

    Raises:
        SimConfigError: The reference names no species.
    """
    names = cfg.species.resolve(token, where)
    if len(names) == 1:
        return find_tip(root, names[0])
    return mrca(root, names)


def window(a: Node, b: Node) -> tuple[float, float]:
    """The span ``(low, high)`` in which both branches exist; empty when ``low >= high``."""
    return max(a.age, b.age), min(a.parent_age, b.parent_age)


def _span(low: float, high: float, per_unit: float) -> str:
    top = "the root" if math.isinf(high) else f"{high / per_unit:g}"
    return f"({low / per_unit:g}, {top})"


def _check_inside(time: float, low: float, high: float, where: str, per_unit: float) -> None:
    if not low < time < high:
        raise SimConfigError(f"{where}: time {time / per_unit:g} lies outside the branch window {_span(low, high, per_unit)}")


def _place(entry: dict, low: float, high: float, where: str, per_unit: float) -> float:
    """An event's time in generations, from ``time`` (checked against the window) or ``time_fraction``."""
    if "time" in entry:
        time = entry["time"] * per_unit
        _check_inside(time, low, high, where, per_unit)
        return time
    if math.isinf(high):
        raise SimConfigError(f"{where}: time_fraction needs a finite window; give time instead")
    return low + entry["time_fraction"] * (high - low)


def _pair(root: Node, cfg: SimConfig, entry: dict, where: str) -> tuple[Node, Node, float, float]:
    donor = resolve_branch(root, cfg, entry["donor"], f"{where}.donor")
    recipient = resolve_branch(root, cfg, entry["recipient"], f"{where}.recipient")
    if donor is recipient:
        raise SimConfigError(f"{where}: donor and recipient are the same branch {donor.label()}")
    low, high = window(donor, recipient)
    if low >= high:
        raise SimConfigError(f"{where}: {donor.label()} and {recipient.label()} never coexist")
    return donor, recipient, low, high


def resolve_events(cfg: SimConfig, root: Node, rng: np.random.Generator) -> list[Event]:
    """Place every configured and random event on the tree and check it fits.

    Args:
        cfg: The validated config.
        root: The replicate's full species tree.
        rng: The replicate's event stream.

    Returns:
        The events, extinctions included, sorted by time.

    Raises:
        SimConfigError: An event names an unknown branch, falls outside the
            window its branches share, or too few random pairs exist.
    """
    per_unit = cfg.generations_per_unit
    events: list[Event] = []

    for i, entry in enumerate(cfg["introgression"]):
        where = f"introgression[{i}]"
        donor, recipient, low, high = _pair(root, cfg, entry, where)
        time = _place(entry, low, high, where, per_unit)
        events.append(Event("pulse", "config", time, donor=donor, recipient=recipient,
                            proportion=entry["proportion"]))

    for i, entry in enumerate(cfg["random_introgression"]):
        events.extend(_random_pulses(cfg, root, rng, entry, f"random_introgression[{i}]"))

    for i, entry in enumerate(cfg["gene_flow"]):
        where = f"gene_flow[{i}]"
        donor, recipient, low, high = _pair(root, cfg, entry, where)
        start = low if entry.get("start") is None else entry["start"] * per_unit
        end = high if entry.get("end") is None else entry["end"] * per_unit
        if math.isinf(end):
            raise SimConfigError(f"{where}: give an end time, the window runs to the root")
        if not (low <= start < end <= high):
            raise SimConfigError(f"{where}: span ({start / per_unit:g}, {end / per_unit:g}) lies outside the "
                                 f"branch window {_span(low, high, per_unit)}")
        pairs = [(donor, recipient)] + ([(recipient, donor)] if entry["symmetric"] else [])
        for d, r in pairs:
            events.append(Event("gene_flow", "config", start, end, donor=d, recipient=r, rate=entry["rate"]))

    for section, kind in (("size_changes", "size_change"), ("bottlenecks", "bottleneck")):
        for i, entry in enumerate(cfg[section]):
            where = f"{section}[{i}]"
            branch = resolve_branch(root, cfg, entry["branch"], f"{where}.branch")
            time = _place(entry, branch.age, branch.parent_age, where, per_unit)
            events.append(Event(
                kind, "config", time, branch=branch, size=entry.get("size"),
                growth_rate=entry.get("growth_rate"), proportion=entry.get("proportion"),
                strength=None if entry.get("strength") is None else entry["strength"] * per_unit,
            ))

    for tip in root.tips():
        if tip.age > 0:
            events.append(Event("extinction", "tree", tip.age, branch=tip))

    events.sort(key=lambda e: e.time)
    return events


def _pool(cfg: SimConfig, root: Node, tokens, kind: str, where: str) -> list[Node]:
    """The branches a random-pulse pool admits.

    A tip qualifies when it is in the pool. An ancestral branch qualifies when
    its clade holds at least one pool species and nothing outside the pool but
    unsampled lineages, so a ghost or extinct lineage nested in a clade does
    not disqualify that clade's ancestor. ``kind`` keeps tips, ancestral
    branches or both; the root, which never coexists with another branch, is
    left out.
    """
    chosen = set(cfg.species.resolve(tokens, where))
    unsampled = set(cfg.species.ghosts) | set(cfg.species.extinct)
    branches = []
    for node in root.preorder():
        if node.parent is None or (kind == "tips" and not node.is_tip) or (kind == "ancestral" and node.is_tip):
            continue
        tips = set(node.tip_names())
        if tips & chosen and tips <= chosen | unsampled:
            branches.append(node)
    return branches


def _random_pulses(cfg: SimConfig, root: Node, rng: np.random.Generator, entry: dict, where: str) -> list[Event]:
    donors = _pool(cfg, root, entry["donors"], entry["donor_branches"], f"{where}.donors")
    recipients = _pool(cfg, root, entry["recipients"], entry["recipient_branches"], f"{where}.recipients")
    pairs = []
    for d in donors:
        for r in recipients:
            if d is r or (not entry["allow_sisters"] and d.parent is r.parent):
                continue
            low, high = window(d, r)
            if low < high:
                pairs.append((d, r, low, high))
    if not pairs:
        hint = "" if "tips" not in (entry["donor_branches"], entry["recipient_branches"]) else \
            "; donor_branches or recipient_branches: any admits ancestral branches"
        raise SimConfigError(f"{where}: no donor/recipient pair coexists{hint}")
    events = []
    # Pairs are drawn with replacement: one pair may take several pulses at different times.
    for slot in rng.choice(len(pairs), size=entry["count"], replace=True):
        d, r, low, high = pairs[int(slot)]
        fraction = rng.uniform(*entry["time_fraction"])
        events.append(Event("pulse", "random", low + fraction * (high - low), donor=d, recipient=r,
                            proportion=float(rng.uniform(*entry["proportion"]))))
    return events


def branch_sizes(cfg: SimConfig, root: Node) -> dict[int, float]:
    """Each branch's starting Ne (node index to size), overrides applied in config order."""
    sizes = {node.index: float(cfg["population"]["size"]) for node in root.preorder()}
    for i, entry in enumerate(cfg["population"]["sizes"]):
        node = resolve_branch(root, cfg, entry["branch"], f"population.sizes[{i}].branch")
        sizes[node.index] = float(entry["size"])
    return sizes


def population_name(node: Node) -> str:
    return f"p{node.index}"


def build_demography(cfg: SimConfig, root: Node, events: list[Event]):
    """Build and validate the msprime demography: one population per branch.

    Args:
        cfg: The validated config.
        root: The replicate's full species tree, nodes indexed.
        events: The resolved events.

    Returns:
        The ``msprime.Demography``.
    """
    import msprime

    demography = msprime.Demography()
    sizes = branch_sizes(cfg, root)
    for node in root.preorder():
        demography.add_population(
            name=population_name(node), initial_size=sizes[node.index], description=node.label(),
        )
    for node in root.preorder():
        if not node.is_tip:
            demography.add_population_split(
                time=node.age,
                derived=[population_name(c) for c in node.children],
                ancestral=population_name(node),
            )
    for event in events:
        if event.kind == "pulse":
            demography.add_mass_migration(
                time=event.time, source=population_name(event.recipient),
                dest=population_name(event.donor), proportion=event.proportion,
            )
        elif event.kind == "gene_flow":
            source, dest = population_name(event.recipient), population_name(event.donor)
            if event.time == 0:
                demography.set_migration_rate(source=source, dest=dest, rate=event.rate)
            else:
                demography.add_migration_rate_change(time=event.time, rate=event.rate, source=source, dest=dest)
            # A split zeroes every rate into or out of its derived populations,
            # so a span that ends at one needs no reset.
            if event.end < min(event.donor.parent_age, event.recipient.parent_age):
                demography.add_migration_rate_change(time=event.end, rate=0, source=source, dest=dest)
        elif event.kind == "size_change":
            demography.add_population_parameters_change(
                time=event.time, initial_size=event.size, growth_rate=event.growth_rate,
                population=population_name(event.branch),
            )
        elif event.kind == "bottleneck":
            if event.strength is not None:
                demography.add_instantaneous_bottleneck(
                    time=event.time, population=population_name(event.branch), strength=event.strength,
                )
            else:
                demography.add_simple_bottleneck(
                    time=event.time, population=population_name(event.branch), proportion=event.proportion,
                )
    demography.sort_events()
    demography.validate()
    return demography
