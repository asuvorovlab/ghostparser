"""Triplet observations read from a per-tree cache of parent, edge-length and
pairwise-LCA tables instead of an extracted subtree.
"""

from array import array
from itertools import combinations
from typing import NamedTuple

from ghostparser.triplet_utils import TOPOLOGY_AB, TOPOLOGY_AC, TOPOLOGY_BC

__all__ = [
    "TRIPLET_MISSING_TAXON",
    "TRIPLET_RESOLVED",
    "TRIPLET_UNRESOLVED",
    "TripletGeometry",
    "TripletSubtreeShape",
    "build_taxon_index",
    "build_triplet_geometry",
    "geometry_observation",
    "triplet_resolution",
    "triplet_subtree_shape",
]

# Why a triplet does or does not yield an observation.
TRIPLET_RESOLVED = "resolved"
TRIPLET_MISSING_TAXON = "missing_taxon"
TRIPLET_UNRESOLVED = "unresolved"


class TripletGeometry(NamedTuple):
    """One rooted gene tree reduced to what triplet queries need.

    Attributes:
        parent: Per-node index of the parent node, ``-1`` at the root. Node
            indices are pre-order ranks, which :func:`geometry_observation`
            relies on when it reproduces the subtree's leaf ordering.
        edge_len: Per-node length of the edge up to its parent; ``0.0`` where
            the Newick omits a length and at the root.
        leaf_node: Per-taxon node index of that taxon's leaf, ``-1`` when the
            taxon is absent from this tree.
        mrca: Flattened ``n_taxa x n_taxa`` table of pairwise LCA node indices,
            ``-1`` where either taxon is unusable.
        n_taxa: Row stride of ``mrca``.
    """

    parent: list
    edge_len: list
    leaf_node: list
    mrca: array
    n_taxa: int


def build_taxon_index(triplets):
    """Map every taxon appearing in the triplets to a dense position.

    Args:
        triplets: Iterable of triplet tuples.

    Returns:
        A dict mapping taxon label to its position, ordered by sorted label.
    """
    labels = sorted({label for triplet in triplets for label in triplet})
    return {label: position for position, label in enumerate(labels)}


def build_triplet_geometry(tree, taxon_index):
    """Cache one gene tree as parent, edge-length, leaf and pairwise-LCA
    arrays, from one pre-order walk and one pass back over it.

    Args:
        tree: A rooted ``Bio.Phylo`` tree.
        taxon_index: Mapping of taxon label to position, from
            :func:`build_taxon_index`.

    Returns:
        A :class:`TripletGeometry` for this tree.
    """
    parent = []
    edge_len = []
    children = []
    n_taxa = len(taxon_index)
    # Labels are unique: reading a tree file rejects a repeated leaf label,
    # so a tree that reached here has at most one leaf per label.
    leaf_node = [-1] * n_taxa
    position_of = []

    stack = [(tree.root, -1)]
    while stack:
        clade, parent_index = stack.pop()
        index = len(parent)
        parent.append(parent_index)
        length = clade.branch_length
        edge_len.append(
            0.0 if parent_index < 0 or length is None else float(length)
        )
        children.append([])
        if parent_index >= 0:
            children[parent_index].append(index)
        position = None
        if not clade.clades:
            position = taxon_index.get(clade.name)
            if position is not None:
                leaf_node[position] = index
        position_of.append(position)
        # Pushed in reverse so the first child is popped, and numbered, first.
        stack.extend((child, index) for child in reversed(clade.clades))

    mrca = array("i", [-1]) * (n_taxa * n_taxa)
    members_of = [None] * len(parent)
    # Every descendant has a higher pre-order rank than its ancestor, so
    # walking the ranks backwards reaches a node after all of its children.
    for node_position in range(len(parent) - 1, -1, -1):
        child_indices = children[node_position]
        if not child_indices:
            position = position_of[node_position]
            members_of[node_position] = [] if position is None else [position]
            continue

        groups = [members_of[child] for child in child_indices]
        for child in child_indices:
            members_of[child] = None
        for left_group, right_group in combinations(groups, 2):
            for left in left_group:
                for right in right_group:
                    mrca[left * n_taxa + right] = node_position
                    mrca[right * n_taxa + left] = node_position
        members_of[node_position] = [member for group in groups for member in group]

    return TripletGeometry(parent, edge_len, leaf_node, mrca, n_taxa)


def _path_sum(parent, edge_len, start, stop):
    """Sum edge lengths walking from ``start`` up to, but excluding, ``stop``.

    Args:
        parent: Parent-index array.
        edge_len: Edge-length array.
        start: Node index to walk from.
        stop: Ancestor node index to stop at.

    Returns:
        The accumulated length, or ``None`` if ``stop`` is not an ancestor of
        ``start``.
    """
    total = 0.0
    current = start
    while current != stop:
        if current < 0:
            return None
        total += edge_len[current]
        current = parent[current]
    return total


class _ResolvedTriplet(NamedTuple):
    """Which nodes a resolved triplet's geometry hangs off.

    Attributes:
        topology: The canonical topology string.
        root_node: Node index of the triplet's own LCA, the subtree's root.
        sister_node: Node index of the sister pair's LCA.
        first_leaf: Node index of the sister leaf named first in ``(A, B, C)``.
        second_leaf: Node index of the other sister leaf.
        odd_leaf: Node index of the leaf outside the sister pair.
    """

    topology: str
    root_node: int
    sister_node: int
    first_leaf: int
    second_leaf: int
    odd_leaf: int


def _resolve_triplet(geometry, triplet_positions):
    """Resolve a triplet's sister pair from its three pairwise LCAs: the pair
    whose LCA differs from the other two is the sister pair, and all
    three equal is a polytomy.

    Args:
        geometry: The tree's :class:`TripletGeometry`.
        triplet_positions: The ``(A, B, C)`` taxon positions.

    Returns:
        A :class:`_ResolvedTriplet`, or ``None`` when a taxon is missing from
        this tree or the triplet is unresolved.
    """
    leaf_node = geometry.leaf_node
    mrca = geometry.mrca
    n_taxa = geometry.n_taxa

    a_position, b_position, c_position = triplet_positions
    a_leaf = leaf_node[a_position]
    b_leaf = leaf_node[b_position]
    c_leaf = leaf_node[c_position]
    if a_leaf < 0 or b_leaf < 0 or c_leaf < 0:
        return None

    ab_mrca = mrca[a_position * n_taxa + b_position]
    ac_mrca = mrca[a_position * n_taxa + c_position]
    bc_mrca = mrca[b_position * n_taxa + c_position]
    if ab_mrca < 0 or ac_mrca < 0 or bc_mrca < 0:
        return None

    if ab_mrca == ac_mrca == bc_mrca:
        return None
    if ac_mrca == bc_mrca:
        return _ResolvedTriplet(
            TOPOLOGY_AB, ac_mrca, ab_mrca, a_leaf, b_leaf, c_leaf
        )
    if ab_mrca == bc_mrca:
        return _ResolvedTriplet(
            TOPOLOGY_AC, ab_mrca, ac_mrca, a_leaf, c_leaf, b_leaf
        )
    if ab_mrca == ac_mrca:
        return _ResolvedTriplet(
            TOPOLOGY_BC, ab_mrca, bc_mrca, b_leaf, c_leaf, a_leaf
        )
    return None


def geometry_observation(
    geometry,
    triplet_positions,
    tree_height_calculation_strategy,
    collect_summary_statistics=False,
):
    """Compute one triplet's observation from a cached gene tree.

    Args:
        geometry: The gene tree's :class:`TripletGeometry`.
        triplet_positions: The ``(A, B, C)`` taxon positions.
        tree_height_calculation_strategy: One of ``AVG``/``A``/``B``/``C``/
            ``SIS``/``INT``.
        collect_summary_statistics: When ``True``, also return the per-tree
            ``avg_tree_height``/``internal_branch``/``sister_distance`` metrics.

    Returns:
        A ``(topology, tree_height, summary_metrics)`` tuple, or ``None`` when
        a taxon is missing or the triplet is unresolved: the same cases where
        the subtree path yields no observation.

    Raises:
        ValueError: If the strategy is unsupported.
    """
    resolved = _resolve_triplet(geometry, triplet_positions)
    if resolved is None:
        return None
    topology = resolved.topology
    root_node = resolved.root_node
    sister_node = resolved.sister_node
    first_leaf = resolved.first_leaf
    second_leaf = resolved.second_leaf
    odd_leaf = resolved.odd_leaf

    parent = geometry.parent
    edge_len = geometry.edge_len

    # An extracted subtree puts each sister leaf two edges below its root and
    # the odd leaf one, so the sister heights are summed as two terms to match
    # how a root-distance walk over that subtree accumulates them.
    internal_branch = _path_sum(parent, edge_len, sister_node, root_node)
    first_height = _path_sum(parent, edge_len, first_leaf, sister_node)
    second_height = _path_sum(parent, edge_len, second_leaf, sister_node)
    odd_height = _path_sum(parent, edge_len, odd_leaf, root_node)
    if None in (internal_branch, first_height, second_height, odd_height):
        return None
    first_height += internal_branch
    second_height += internal_branch

    if topology == TOPOLOGY_AB:
        a_height, b_height, c_height = first_height, second_height, odd_height
    elif topology == TOPOLOGY_AC:
        a_height, b_height, c_height = first_height, odd_height, second_height
    else:
        a_height, b_height, c_height = odd_height, first_height, second_height

    # An extracted subtree averages in its own leaf order, which is ascending
    # pre-order index.
    ordered = sorted(
        (
            (first_leaf, first_height),
            (second_leaf, second_height),
            (odd_leaf, odd_height),
        )
    )
    avg_tree_height = (ordered[0][1] + ordered[1][1] + ordered[2][1]) / 3.0
    sister_distance = first_height + second_height - 2.0 * internal_branch

    if tree_height_calculation_strategy == "AVG":
        tree_height = avg_tree_height
    elif tree_height_calculation_strategy == "A":
        tree_height = a_height
    elif tree_height_calculation_strategy == "B":
        tree_height = b_height
    elif tree_height_calculation_strategy == "C":
        tree_height = c_height
    elif tree_height_calculation_strategy == "SIS":
        tree_height = sister_distance
    elif tree_height_calculation_strategy == "INT":
        tree_height = internal_branch
    else:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}"
        )

    summary_metrics = None
    if collect_summary_statistics:
        summary_metrics = {
            "avg_tree_height": avg_tree_height,
            "internal_branch": internal_branch,
            "sister_distance": sister_distance,
        }

    return (topology, tree_height, summary_metrics)


def triplet_resolution(geometry, triplet_positions):
    """Name why a triplet does or does not yield an observation, mirroring the
    guards in :func:`geometry_observation`.

    Args:
        geometry: The gene tree's :class:`TripletGeometry`.
        triplet_positions: The ``(A, B, C)`` taxon positions.

    Returns:
        :data:`TRIPLET_RESOLVED`, :data:`TRIPLET_MISSING_TAXON` when a taxon has
        no leaf in this tree, or :data:`TRIPLET_UNRESOLVED` when all three pair
        LCAs coincide, which is the polytomy the subtree path cannot classify.
    """
    leaf_node = geometry.leaf_node
    mrca = geometry.mrca
    n_taxa = geometry.n_taxa
    a_position, b_position, c_position = triplet_positions

    if (
        leaf_node[a_position] < 0
        or leaf_node[b_position] < 0
        or leaf_node[c_position] < 0
    ):
        return TRIPLET_MISSING_TAXON

    ab_mrca = mrca[a_position * n_taxa + b_position]
    ac_mrca = mrca[a_position * n_taxa + c_position]
    bc_mrca = mrca[b_position * n_taxa + c_position]
    if ab_mrca < 0 or ac_mrca < 0 or bc_mrca < 0:
        return TRIPLET_MISSING_TAXON
    if ab_mrca == ac_mrca == bc_mrca:
        return TRIPLET_UNRESOLVED
    return TRIPLET_RESOLVED


class TripletSubtreeShape(NamedTuple):
    """A triplet's induced subtree, as the ordering needed to write it.

    The induced subtree is always ``((s1,s2),odd)``: two sister leaves under
    their own LCA, and the third leaf attached directly to the triplet's LCA.

    Attributes:
        sister_positions: The two sister taxa's positions, in the order the
            subtree lists them.
        odd_position: The remaining taxon's position.
        sister_clade_first: Whether the sister clade is listed before the odd
            leaf among the triplet LCA's children.
    """

    sister_positions: tuple
    odd_position: int
    sister_clade_first: bool


def triplet_subtree_shape(geometry, triplet_positions):
    """Describe the subtree a triplet induces, without building it.

    Node indices are pre-order ranks and a subtree copy preserves child order,
    so comparing them recovers the order the copied subtree would list its
    children in.

    Args:
        geometry: The tree's :class:`TripletGeometry`.
        triplet_positions: The ``(A, B, C)`` taxon positions.

    Returns:
        A :class:`TripletSubtreeShape`, or ``None`` when a taxon is missing or
        the triplet is unresolved.
    """
    resolved = _resolve_triplet(geometry, triplet_positions)
    if resolved is None:
        return None

    a_position, b_position, c_position = triplet_positions
    if resolved.topology == TOPOLOGY_AB:
        sisters, odd_position = (a_position, b_position), c_position
    elif resolved.topology == TOPOLOGY_AC:
        sisters, odd_position = (a_position, c_position), b_position
    else:
        sisters, odd_position = (b_position, c_position), a_position

    if resolved.first_leaf > resolved.second_leaf:
        sisters = (sisters[1], sisters[0])

    return TripletSubtreeShape(
        sister_positions=sisters,
        odd_position=odd_position,
        sister_clade_first=resolved.sister_node < resolved.odd_leaf,
    )
