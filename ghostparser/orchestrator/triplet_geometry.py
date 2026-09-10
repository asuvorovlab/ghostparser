"""Triplet geometry read directly from a cached gene tree.

Produces the same ``(topology, tree_height, summary_metrics)`` observation as
extracting a triplet's subtree with DendroPy, without copying the tree. A gene
tree is cached once as parent/edge-length arrays plus a pairwise LCA table;
each triplet's geometry then follows from three table lookups and a few short
walks up the parent chain.
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
        root_edge: The tree root's own edge length, or ``None`` when the Newick
            gives it none. Held separately because ``edge_len`` carries ``0.0``
            there: no distance walk passes through the root, but writing a
            subtree out needs the real value.
    """

    parent: list
    edge_len: list
    leaf_node: list
    mrca: array
    n_taxa: int
    root_edge: float | None


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
    """Cache one gene tree's parent, edge-length, leaf and pairwise-LCA arrays.

    Walks the tree twice: pre-order to number the nodes and record the parent
    chain, then post-order to fill the LCA table. A node is the LCA of exactly
    those pairs drawn from two different children, so each pair is written once.

    Args:
        tree: A rooted DendroPy tree.
        taxon_index: Mapping of taxon label to position, from
            :func:`build_taxon_index`.

    Returns:
        A :class:`TripletGeometry` for this tree.
    """
    parent = []
    edge_len = []
    index_of = {}
    root_edge = None
    for node in tree.preorder_node_iter():
        index_of[id(node)] = len(parent)
        parent_node = node.parent_node
        if parent_node is None:
            parent.append(-1)
            edge_len.append(0.0)
            if node.edge_length is not None:
                root_edge = float(node.edge_length)
            continue
        parent.append(index_of[id(parent_node)])
        length = node.edge_length
        edge_len.append(0.0 if length is None else float(length))

    n_taxa = len(taxon_index)
    # Labels are unique: DendroPy rejects duplicate taxa when parsing, so a
    # tree that reached here has at most one leaf per label.
    leaf_node = [-1] * n_taxa
    for leaf in tree.leaf_node_iter():
        taxon = leaf.taxon
        if taxon is not None and taxon.label in taxon_index:
            leaf_node[taxon_index[taxon.label]] = index_of[id(leaf)]

    mrca = array("i", [-1]) * (n_taxa * n_taxa)
    members_of = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            taxon = node.taxon
            position = None if taxon is None else taxon_index.get(taxon.label)
            members_of[id(node)] = [] if position is None else [position]
            continue

        groups = [members_of.pop(id(child)) for child in node.child_node_iter()]
        node_position = index_of[id(node)]
        for left_group, right_group in combinations(groups, 2):
            for left in left_group:
                for right in right_group:
                    mrca[left * n_taxa + right] = node_position
                    mrca[right * n_taxa + left] = node_position
        members_of[id(node)] = [member for group in groups for member in group]

    return TripletGeometry(parent, edge_len, leaf_node, mrca, n_taxa, root_edge)


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
    """Resolve a triplet's sister pair from its pairwise LCA node indices.

    Exactly two of the three pair LCAs are the triplet's own LCA and the odd one
    out is the sister pair's, so which pair is odd *is* the topology -- no
    distance enters the call, which is what keeps it exact and what makes a
    zero-length internal branch harmless. All three equal is an unresolved
    polytomy, where the subtree path's ``find_sister_pair`` raises.

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

    # The extracted subtree puts each sister leaf two edges below its root and
    # the odd leaf one, so the sister heights are summed as two terms to match
    # how _distance_to_root accumulates them there.
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

    # _compute_triplet_tree_metrics averages in the subtree's own leaf order,
    # which is ascending pre-order index.
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
    """Name why a triplet does or does not yield an observation.

    Mirrors the two guards in :func:`geometry_observation`, which returns
    ``None`` for exactly the outcomes this reports as not resolved. It is kept
    separate so the hot path carries no reason tracking;
    ``test_triplet_resolution_agrees_with_the_observation_guards`` pins the two
    against each other so they cannot drift.

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
    """A triplet's induced subtree, as the edges and ordering needed to write it.

    The induced subtree is always ``((s1,s2),odd)``: two sister leaves under
    their own LCA, and the third leaf attached directly to the triplet's LCA.

    Attributes:
        sister_positions: The two sister taxa's positions, in the order the
            subtree lists them.
        sister_edges: Each sister leaf's edge up to the sister LCA, aligned with
            ``sister_positions``.
        internal_edge: The sister LCA's edge up to the triplet LCA.
        odd_position: The remaining taxon's position.
        odd_edge: The odd leaf's edge up to the triplet LCA.
        sister_clade_first: Whether the sister clade is listed before the odd
            leaf among the triplet LCA's children.
        root_edge: The triplet LCA's own depth, which a copied-out subtree keeps
            on its root, or ``None`` when the triplet LCA is the tree's root and
            so has no edge above it.
    """

    sister_positions: tuple
    sister_edges: tuple
    internal_edge: float
    odd_position: int
    odd_edge: float
    sister_clade_first: bool
    root_edge: float | None


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

    parent = geometry.parent
    edge_len = geometry.edge_len
    a_position, b_position, c_position = triplet_positions
    if resolved.topology == TOPOLOGY_AB:
        sisters, odd_position = (a_position, b_position), c_position
    elif resolved.topology == TOPOLOGY_AC:
        sisters, odd_position = (a_position, c_position), b_position
    else:
        sisters, odd_position = (b_position, c_position), a_position

    internal_edge = _path_sum(parent, edge_len, resolved.sister_node, resolved.root_node)
    first_edge = _path_sum(parent, edge_len, resolved.first_leaf, resolved.sister_node)
    second_edge = _path_sum(parent, edge_len, resolved.second_leaf, resolved.sister_node)
    odd_edge = _path_sum(parent, edge_len, resolved.odd_leaf, resolved.root_node)
    if None in (internal_edge, first_edge, second_edge, odd_edge):
        return None

    sister_leaves = (resolved.first_leaf, resolved.second_leaf)
    sister_edges = (first_edge, second_edge)
    if sister_leaves[0] > sister_leaves[1]:
        sisters = (sisters[1], sisters[0])
        sister_edges = (sister_edges[1], sister_edges[0])

    # A copied-out subtree keeps an edge above its root: unifurcation
    # suppression collapses the whole path from the tree's root down to the
    # triplet LCA into it, the tree root's own edge included. When the triplet
    # LCA *is* the tree root, only that root edge remains -- and if the Newick
    # gave none, the subtree has none either.
    root_edge = geometry.root_edge
    if parent[resolved.root_node] >= 0:
        depth = _path_sum(parent, edge_len, resolved.root_node, -1)
        if depth is None:
            return None
        root_edge = depth + (geometry.root_edge or 0.0)

    return TripletSubtreeShape(
        sister_positions=sisters,
        sister_edges=sister_edges,
        internal_edge=internal_edge,
        odd_position=odd_position,
        odd_edge=odd_edge,
        sister_clade_first=resolved.sister_node < resolved.odd_leaf,
        root_edge=root_edge,
    )
