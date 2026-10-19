"""Species trees as ages: Newick in and out, birth-death simulation, extinction and outgroups.

Every node carries the age (time before the present) at which its branch ends:
0 for an extant tip, the extinction time for an extinct tip, the split time for
an internal node. A branch therefore spans ``(node.age, node.parent.age)``.
"""

from __future__ import annotations

import math
from io import StringIO
from typing import Callable, Iterator

import numpy as np

from sim_config import SimConfig, SimInputError

_ULTRAMETRIC_TOLERANCE = 1e-6
# GSA stops a run once the lineage count reaches this many above the target and
# the chance of ever falling back to the target is below this probability.
_GSA_RETURN_PROBABILITY = 1e-6
_GSA_MAX_OVERSHOOT = 5000


class Node:
    """One branch of a species tree, named by the node at its lower end."""

    __slots__ = ("name", "age", "parent", "children", "index")

    def __init__(self, name: str | None = None, age: float = 0.0, children: list[Node] | None = None):
        self.name = name
        self.age = age
        self.parent: Node | None = None
        self.children: list[Node] = []
        self.index = -1
        for child in children or []:
            self.add_child(child)

    def add_child(self, child: Node) -> None:
        child.parent = self
        self.children.append(child)

    @property
    def is_tip(self) -> bool:
        return not self.children

    @property
    def parent_age(self) -> float:
        return self.parent.age if self.parent is not None else math.inf

    def preorder(self) -> Iterator[Node]:
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def postorder(self) -> Iterator[Node]:
        for child in self.children:
            yield from child.postorder()
        yield self

    def tips(self) -> list[Node]:
        return [n for n in self.preorder() if n.is_tip]

    def tip_names(self) -> list[str]:
        return [n.name for n in self.tips()]

    def label(self) -> str:
        """The tip name, or ``{A,B,...}`` (its tips) for an internal branch."""
        return self.name if self.is_tip else "{" + ",".join(self.tip_names()) + "}"


def index_nodes(root: Node) -> None:
    """Number the nodes in preorder, so each has a stable population index."""
    for i, node in enumerate(root.preorder()):
        node.index = i


def find_tip(root: Node, name: str) -> Node:
    for tip in root.tips():
        if tip.name == name:
            return tip
    raise KeyError(name)


def mrca(root: Node, names: list[str] | tuple[str, ...]) -> Node:
    """The most recent common ancestor of the named tips (the tip itself for one name)."""
    wanted = set(names)
    best = root
    for node in root.postorder():
        if wanted <= set(node.tip_names()):
            best = node
            break
    return best


def prune(root: Node, keep: set[str]) -> Node:
    """Copy the tree keeping only the named tips, collapsing nodes left with one child.

    Ages are absolute, so collapsing a node needs no branch-length bookkeeping.
    """

    def copy(node: Node) -> Node | None:
        if node.is_tip:
            return Node(node.name, node.age) if node.name in keep else None
        kids = [c for c in (copy(child) for child in node.children) if c is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        return Node(node.name, node.age, kids)

    pruned = copy(root)
    if pruned is None:
        raise ValueError("pruning removed every tip")
    return pruned


def to_newick(root: Node, length: Callable[[Node], float], decimals: int) -> str:
    """Write a Newick string; ``length(node)`` gives each branch's length."""

    def write(node: Node) -> str:
        body = node.name if node.is_tip else "(" + ",".join(write(c) for c in node.children) + ")"
        if node.parent is None:
            return body
        return f"{body}:{length(node):.{decimals}f}"

    return write(root) + ";"


def parse_newick(newick: str, extinct: set[str], scale: float) -> Node:
    """Read a species tree, turning branch lengths into ages.

    Args:
        newick: The tree; every branch below the root needs a length.
        extinct: Tips allowed to end above the present.
        scale: Generations per Newick length unit.

    Returns:
        The root, with ages in generations.

    Raises:
        SimInputError: A length is missing or negative, an internal branch has
            no length, or a tip that is not extinct ends above the present.
    """
    from Bio import Phylo

    tree = Phylo.read(StringIO(newick), "newick")
    depth: dict[int, float] = {}

    def build(clade, parent_depth: float, is_root: bool) -> Node:
        if is_root:
            here = 0.0
        else:
            if clade.branch_length is None:
                where = clade.name or "an internal branch"
                raise SimInputError(f"Species tree branch {where} has no length")
            if clade.branch_length < 0:
                raise SimInputError(f"Species tree branch {clade.name or 'internal'} has a negative length")
            here = parent_depth + clade.branch_length * scale
        node = Node(str(clade.name) if clade.is_terminal() else None)
        depth[id(node)] = here
        for child in clade.clades:
            node.add_child(build(child, here, False))
        if len(node.children) == 1:  # collapse a unary node into its child
            only = node.children[0]
            only.parent = None
            return only
        return node

    root = build(tree.root, 0.0, True)
    height = max(depth[id(t)] for t in root.tips())
    for node in root.preorder():
        node.age = height - depth[id(node)]
    for tip in root.tips():
        if tip.name not in extinct:
            if tip.age > _ULTRAMETRIC_TOLERANCE * max(height, 1.0):
                raise SimInputError(
                    f"Tip {tip.name} ends {tip.age / scale:g} time units above the present; "
                    "only species listed as extinct may"
                )
            tip.age = 0.0
    for node in root.preorder():
        if not node.is_tip and any(node.age <= c.age for c in node.children):
            raise SimInputError(f"Internal branch above {node.label()} has zero length")
    return root


def _gsa_complete_tree(rng: np.random.Generator, n: int, birth: float, death: float,
                       max_attempts: int) -> Node:
    """Simulate a complete birth-death tree with exactly ``n`` lineages alive at the end.

    General sampling approach (Hartmann, Wong & Stadler 2010): run the process
    forward from one lineage until it holds far more than ``n`` lineages, pick a
    time at which exactly ``n`` were alive with probability proportional to how
    long each such spell lasted, and read the tree off at that time. Extinct
    lineages stay in as tips ending above the present.

    Raises:
        SimInputError: Every attempt died out before reaching ``n`` lineages.
    """
    if death == 0:
        ceiling = n + 1
    else:
        overshoot = math.ceil(math.log(_GSA_RETURN_PROBABILITY) / math.log(death / birth))
        ceiling = n + min(max(overshoot, 1), _GSA_MAX_OVERSHOOT)
    p_birth = birth / (birth + death)
    for _ in range(max_attempts):
        start = [0.0]        # lineage -> time it began
        end: list[float | None] = [None]   # lineage -> time it split or died
        kids: list[tuple[int, int] | None] = [None]
        alive = [0]
        t = 0.0
        spells: list[tuple[float, float]] = []
        while 0 < len(alive) < ceiling:
            k = len(alive)
            dt = rng.exponential(1.0 / (k * (birth + death)))
            if k == n:
                spells.append((t, t + dt))
            t += dt
            slot = int(rng.integers(k))
            lineage = alive[slot]
            end[lineage] = t
            if rng.random() < p_birth:
                left, right = len(start), len(start) + 1
                start += [t, t]
                end += [None, None]
                kids += [None, None]
                kids[lineage] = (left, right)
                alive[slot] = left
                alive.append(right)
            else:
                alive[slot] = alive[-1]
                alive.pop()
        if not spells:
            continue
        lengths = np.array([b - a for a, b in spells])
        a, b = spells[int(rng.choice(len(spells), p=lengths / lengths.sum()))]
        when = rng.uniform(a, b)

        def build(lineage: int) -> Node:
            stop = end[lineage]
            if stop is None or stop > when:
                return Node(age=0.0)
            if kids[lineage] is None:
                return Node(age=when - stop)
            return Node(age=when - stop, children=[build(c) for c in kids[lineage]])

        return build(0)
    raise SimInputError(f"Birth-death process died out in all {max_attempts} attempts")


def _scale(root: Node, target_age: float) -> None:
    factor = target_age / root.age
    for node in root.preorder():
        node.age *= factor


def _reconstructed_tree(rng: np.random.Generator, names: list[str], birth: float, death: float,
                        max_attempts: int, root_age: float) -> Node:
    """A tree over ``names`` with extinct lineages pruned, scaled to ``root_age``, tips labeled at random."""
    if len(names) == 1:
        return Node(names[0], 0.0)
    complete = _gsa_complete_tree(rng, len(names), birth, death, max_attempts)
    extant = [t for t in complete.tips() if t.age == 0.0]
    for tip in complete.tips():
        tip.name = None
    for tip, name in zip(extant, rng.permutation(names)):
        tip.name = str(name)
    root = prune(complete, set(names))
    _scale(root, root_age)
    return root


def _process_tree(rng: np.random.Generator, extant: list[str], extinct: list[str], birth: float,
                  death: float, max_attempts: int, root_age: float) -> Node:
    """A complete birth-death tree whose own extinct lineages carry the extinct names."""
    for _ in range(max_attempts):
        complete = _gsa_complete_tree(rng, len(extant), birth, death, max_attempts)
        dead = [t for t in complete.tips() if t.age > 0.0]
        if len(dead) < len(extinct):
            continue
        living = [t for t in complete.tips() if t.age == 0.0]
        for tip, name in zip(living, rng.permutation(extant)):
            tip.name = str(name)
        for slot, name in zip(rng.choice(len(dead), size=len(extinct), replace=False), rng.permutation(extinct)):
            dead[int(slot)].name = str(name)
        root = prune(complete, set(extant) | set(extinct))
        _scale(root, root_age)
        return root
    raise SimInputError(
        f"No birth-death tree with {len(extinct)} extinct lineage(s) in {max_attempts} attempts; "
        "raise death_rate or max_attempts, or set extinct_from_process: false"
    )


def _truncate(rng: np.random.Generator, root: Node, names: set[str], fraction: list[float]) -> None:
    """End each named tip that reaches the present partway along its terminal branch."""
    lo, hi = fraction
    for tip in root.tips():
        if tip.name in names and tip.age == 0.0:
            survived = rng.uniform(lo, hi)
            tip.age = tip.parent_age * (1.0 - survived)


def build_species_tree(cfg: SimConfig, rng: np.random.Generator) -> Node:
    """Build one replicate's full species tree (ghosts and extinct tips included), ages in generations.

    Args:
        cfg: The validated config.
        rng: The replicate's tree stream.

    Returns:
        The root, with nodes indexed in preorder.

    Raises:
        SimInputError: The Newick tree is unusable, or the birth-death process
            could not produce the requested tree.
    """
    species = cfg.species
    tree = cfg.tree
    extinct = set(species.extinct)
    per_unit = cfg.generations_per_unit
    if tree["model"] == "newick":
        root = parse_newick(cfg.newick, extinct, per_unit)
        _truncate(rng, root, extinct, tree["extinction_fraction"])
        index_nodes(root)
        return root

    birth = tree["birth_rate"]
    death = tree["death_rate"] or 0.0
    attempts = tree["max_attempts"]
    root_age = tree["root_age"] * per_unit
    ingroup = list(species.ingroup)
    if tree["extinct_from_process"]:
        extant = [n for n in ingroup if n not in extinct]
        ingroup_root = _process_tree(rng, extant, [n for n in ingroup if n in extinct],
                                     birth, death, attempts, root_age)
    else:
        ingroup_root = _reconstructed_tree(rng, ingroup, birth, death, attempts, root_age)
        _truncate(rng, ingroup_root, extinct, tree["extinction_fraction"])

    outgroup_root = _reconstructed_tree(rng, list(species.outgroup), birth, death, attempts, root_age)
    root = Node(age=root_age * tree["outgroup_age"], children=[ingroup_root, outgroup_root])
    index_nodes(root)
    return root
