"""Cross-library parity: DendroPy triplet extraction/collapsing agrees with a BioPython implementation.

This is the suite's only parity test. It is retained because triplet collapsing
is the one place where GhostParser relies on DendroPy's subtree semantics being
equivalent to the standard BioPython pruning behaviour; every other test derives
its expectations from definitions rather than from a second implementation.
"""

import copy
from io import StringIO

import dendropy
import pytest
from Bio import Phylo

from ghostparser.orchestrator.trees import (
    extract_triplet_subtree,
    format_newick_with_precision,
)


def _dist_to_root(node):
    """Sum edge lengths from a DendroPy node up to the root.

    Args:
        node: The DendroPy node to measure from.

    Returns:
        The accumulated root distance as a float.
    """
    dist = 0.0
    while node is not None and node.edge is not None:
        if node.edge.length is not None:
            dist += node.edge.length
        node = node.parent_node
    return dist


def _dendro_distance(tree, a, b):
    """Compute the patristic distance between two taxa in a DendroPy tree.

    Args:
        tree: The DendroPy tree.
        a: First taxon label.
        b: Second taxon label.

    Returns:
        The pairwise distance as a float.
    """
    tree.is_rooted = True
    node_a = tree.find_node_with_taxon_label(a)
    node_b = tree.find_node_with_taxon_label(b)
    mrca = tree.mrca(taxon_labels=[a, b])
    return _dist_to_root(node_a) + _dist_to_root(node_b) - 2.0 * _dist_to_root(mrca)


def _bio_distance(tree, a, b):
    """Compute the patristic distance between two taxa in a BioPython tree.

    Args:
        tree: The BioPython Phylo tree.
        a: First taxon label.
        b: Second taxon label.

    Returns:
        The pairwise distance as a float.
    """
    return tree.distance(a, b)


def _collapse_triplet_biopython(newick_str, triplet):
    """Collapse a tree to a triplet using BioPython's pruning as a reference.

    Args:
        newick_str: The source tree in Newick format.
        triplet: The three taxon labels to retain.

    Returns:
        The pruned BioPython tree, or ``None`` if the triplet is not present.
    """
    tree = Phylo.read(StringIO(newick_str), "newick")
    collapsed = copy.deepcopy(tree)

    triplet_set = set(triplet)
    terminals = {terminal.name for terminal in collapsed.get_terminals()}
    if not triplet_set.issubset(terminals):
        return None

    for terminal in list(collapsed.get_terminals()):
        if terminal.name not in triplet_set:
            collapsed.prune(target=terminal)

    return collapsed


@pytest.mark.parity
@pytest.mark.parametrize("a,b", [("A", "B"), ("A", "C"), ("B", "C")])
def test_triplet_branch_lengths_match(triplet_comparison_cases, a, b):
    """DendroPy-extracted triplet distances survive a BioPython Newick round trip."""
    for newick_str, triplet in triplet_comparison_cases:
        dendro_tree = dendropy.Tree.get(
            data=newick_str, schema="newick", preserve_underscores=True
        )

        dendro_subtree = extract_triplet_subtree(dendro_tree, triplet)
        assert dendro_subtree is not None

        dendro_newick = format_newick_with_precision(dendro_subtree)
        bio_subtree = Phylo.read(StringIO(dendro_newick), "newick")

        assert bio_subtree is not None

        bio_dist = _bio_distance(bio_subtree, a, b)
        dendro_dist = _dendro_distance(dendro_subtree, a, b)

        assert bio_dist == pytest.approx(dendro_dist)


@pytest.mark.parity
def test_triplet_collapse_consistency_dendropy_vs_biopython(triplet_comparison_cases):
    """DendroPy triplet collapsing yields the same distances as BioPython pruning."""
    for newick_str, triplet in triplet_comparison_cases:
        dendro_tree = dendropy.Tree.get(
            data=newick_str, schema="newick", preserve_underscores=True
        )
        dendro_subtree = extract_triplet_subtree(dendro_tree, triplet)
        assert dendro_subtree is not None

        bio_subtree = _collapse_triplet_biopython(newick_str, triplet)
        assert bio_subtree is not None

        for a, b in (
            (triplet[0], triplet[1]),
            (triplet[0], triplet[2]),
            (triplet[1], triplet[2]),
        ):
            dendro_dist = _dendro_distance(dendro_subtree, a, b)
            bio_dist = _bio_distance(bio_subtree, a, b)
            assert dendro_dist == pytest.approx(bio_dist, rel=0.0, abs=1e-12)
