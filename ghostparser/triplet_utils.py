"""Shared triplet topology utilities.

Utilities in this module keep triplet topology handling consistent across the
package; ``ghostparser.orchestrator`` uses them for both preprocessing and
inference.
"""

from __future__ import annotations


TOPOLOGY_AB = "((A,B),C)"
TOPOLOGY_BC = "((B,C),A)"
TOPOLOGY_AC = "((A,C),B)"
ALL_TOPOLOGIES = (TOPOLOGY_AB, TOPOLOGY_BC, TOPOLOGY_AC)


def triplet_taxa_labels(tree):
    """Return sorted labels for a rooted 3-tip tree.

    Raises:
        ValueError: If tree does not have exactly 3 labeled leaves.
    """
    labels = sorted(
        leaf.taxon.label
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None and leaf.taxon.label
    )
    if len(labels) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")
    return labels


def find_sister_pair(tree):
    """Return sister-pair labels for a rooted 3-tip tree as ``frozenset``."""
    tree.is_rooted = True
    labels = triplet_taxa_labels(tree)
    root = tree.seed_node

    candidate_pairs = (
        (labels[0], labels[1]),
        (labels[0], labels[2]),
        (labels[1], labels[2]),
    )

    for left, right in candidate_pairs:
        mrca = tree.mrca(taxon_labels=[left, right])
        if mrca is not None and mrca is not root:
            return frozenset((left, right))

    raise ValueError("Could not determine rooted sister pair for triplet tree")


def normalize_abc_from_sister_pair(labels, sister_pair):
    """Normalize triplet labels to ``(A, B, C)`` where ``A`` and ``B`` are sisters."""
    labels_set = set(labels)
    a_taxon, b_taxon = sorted(sister_pair)
    c_taxon = next(iter(labels_set - set(sister_pair)))
    return a_taxon, b_taxon, c_taxon


def topology_from_sister_pair(sister_pair, abc_triplet):
    """Map sister pair to one of ``((A,B),C)``, ``((B,C),A)``, ``((A,C),B)``."""
    a_taxon, b_taxon, c_taxon = abc_triplet

    if sister_pair == frozenset((a_taxon, b_taxon)):
        return TOPOLOGY_AB
    if sister_pair == frozenset((b_taxon, c_taxon)):
        return TOPOLOGY_BC
    if sister_pair == frozenset((a_taxon, c_taxon)):
        return TOPOLOGY_AC

    raise ValueError("Tree taxa do not match provided ABC triplet")


def classify_triplet_topology_string(tree, abc_triplet):
    """Classify rooted triplet topology as one of the three canonical strings."""
    sister_pair = find_sister_pair(tree)
    return topology_from_sister_pair(sister_pair, abc_triplet)


def rank_topologies_by_frequency(topology_counts, rng=None):
    """Rank topologies as (con, dis1, dis2) by descending frequency.

    Ties are broken deterministically using ``ALL_TOPOLOGIES`` order.
    """
    ranked = sorted(
        ALL_TOPOLOGIES,
        key=lambda topology: (-int(topology_counts.get(topology, 0)), ALL_TOPOLOGIES.index(topology)),
    )

    return ranked[0], ranked[1], ranked[2]
