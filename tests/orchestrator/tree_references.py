"""The two reference implementations the parity tests hold the cached geometry
to. The DendroPy path extracts each triplet's subtree and walks it; the
BioPython path measures the unpruned tree from the definitions alone. No run
uses either, so they live with the tests, and DendroPy is a ``dev`` dependency
only.
"""

from io import StringIO

import dendropy
from Bio import Phylo

from ghostparser.orchestrator.config import (
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
)
from ghostparser.triplet_utils import TOPOLOGY_AB, TOPOLOGY_AC, TOPOLOGY_BC


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract the subtree spanning only the triplet taxa.

    Args:
        tree: A DendroPy tree object.
        triplet_taxa: Iterable of the three taxon names.

    Returns:
        The extracted subtree, or ``None`` if any triplet taxon is absent.
    """
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon}
    if not set(triplet_taxa).issubset(tree_taxa):
        return None

    return tree.extract_tree_with_taxa_labels(triplet_taxa)


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


def distance_to_root(node):
    """Compute the root-to-node distance from edge lengths.

    Missing edge lengths are treated as zero.

    Args:
        node: A DendroPy node.

    Returns:
        The accumulated distance to the root as a float.
    """
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node
    return distance


def compute_triplet_tree_metrics(
    tree,
    species_triplet=None,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Compute the selected tree-height value H(T) and optional summary metrics.

    Args:
        tree: A rooted 3-tip DendroPy tree.
        species_triplet: The ``(A, B, C)`` triplet, required for the ``A``/``B``/
            ``C`` strategies.
        tree_height_calculation_strategy: One of ``AVG``/``A``/``B``/``C``/
            ``SIS``/``INT``.
        collect_summary_statistics: When ``True``, also compute the per-tree
            ``avg_tree_height``/``internal_branch``/``sister_distance`` metrics
            used for summary-statistics gathering.

    Returns:
        A tuple ``(selected_tree_height, summary_metrics)`` where
        ``summary_metrics`` is a dict of the three metrics when
        ``collect_summary_statistics`` is ``True`` (else ``None``).

    Raises:
        ValueError: If the strategy is unsupported, the tree does not have
            exactly three tips, ``species_triplet`` is missing for A/B/C, or the
            sister-pair MRCA cannot be determined.
    """
    if tree_height_calculation_strategy not in TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    leaves = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label]
    if len(leaves) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")

    leaf_by_label = {leaf.taxon.label: leaf for leaf in leaves}
    leaf_distances = {
        label: distance_to_root(leaf) for label, leaf in leaf_by_label.items()
    }
    avg_tree_height = sum(leaf_distances.values()) / 3.0

    selected_tree_height = None
    if tree_height_calculation_strategy in {"A", "B", "C"}:
        if species_triplet is None:
            raise ValueError(
                "species_triplet is required for tree height strategies A, B, and C"
            )

        strategy_index = {"A": 0, "B": 1, "C": 2}[tree_height_calculation_strategy]
        selected_taxon_label = species_triplet[strategy_index]
        if selected_taxon_label not in leaf_distances:
            raise ValueError(
                f"Selected taxon {selected_taxon_label} not found in triplet tree"
            )
        selected_tree_height = leaf_distances[selected_taxon_label]
    elif tree_height_calculation_strategy == "AVG":
        selected_tree_height = avg_tree_height

    summary_metrics = None
    if collect_summary_statistics or tree_height_calculation_strategy in {"SIS", "INT"}:
        sister_pair = find_sister_pair(tree)
        left_label, right_label = tuple(sister_pair)
        sister_mrca = tree.mrca(taxon_labels=[left_label, right_label])
        if sister_mrca is None:
            raise ValueError("Could not determine sister-pair MRCA for triplet tree")

        internal_branch = distance_to_root(sister_mrca)
        sister_distance = (
            leaf_distances[left_label]
            + leaf_distances[right_label]
            - 2.0 * internal_branch
        )

        if tree_height_calculation_strategy == "SIS":
            selected_tree_height = sister_distance
        elif tree_height_calculation_strategy == "INT":
            selected_tree_height = internal_branch

        if collect_summary_statistics:
            summary_metrics = {
                "avg_tree_height": avg_tree_height,
                "internal_branch": internal_branch,
                "sister_distance": sister_distance,
            }

    if selected_tree_height is None:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    return selected_tree_height, summary_metrics


def observation_from_subtree(
    subtree,
    triplet,
    tree_height_calculation_strategy,
    collect_summary_statistics=False,
):
    """Compute a ``(topology, tree_height, metrics)`` observation by walking a
    DendroPy subtree.

    Args:
        subtree: The extracted triplet subtree as a DendroPy tree.
        triplet: The ``(A, B, C)`` triplet.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, also compute the per-tree
            summary metrics stored as the observation's third element.

    Returns:
        A ``(topology, tree_height, summary_metrics)`` tuple where
        ``summary_metrics`` is ``None`` unless ``collect_summary_statistics`` is
        ``True``, or ``None`` if the subtree's labels do not match the triplet or
        metric computation fails.
    """
    labels = {
        leaf.taxon.label
        for leaf in subtree.leaf_node_iter()
        if leaf.taxon and leaf.taxon.label
    }
    if labels != set(triplet):
        return None

    try:
        topology = classify_triplet_topology_string(subtree, triplet)
        tree_height, summary_metrics = compute_triplet_tree_metrics(
            subtree,
            species_triplet=triplet,
            tree_height_calculation_strategy=tree_height_calculation_strategy,
            collect_summary_statistics=collect_summary_statistics,
        )
    except ValueError:
        return None

    return (topology, tree_height, summary_metrics)


def serialize_triplet_gene_trees(
    species_triplet,
    triplet_gene_trees,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Parse rooted triplet Newicks into observations.

    Trees whose leaf set does not match the triplet, or that fail metric
    computation, are skipped.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        triplet_gene_trees: Iterable of rooted triplet Newick strings.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, also compute each
            observation's summary metrics (third tuple element).

    Returns:
        A list of ``(topology, tree_height, summary_metrics)`` observation
        tuples.
    """
    observations = []
    species_set = set(species_triplet)

    for newick_str in triplet_gene_trees:
        if not str(newick_str).strip():
            continue

        tree = dendropy.Tree.get(
            data=str(newick_str).strip(), schema="newick", preserve_underscores=True
        )
        labels = {
            leaf.taxon.label
            for leaf in tree.leaf_node_iter()
            if leaf.taxon and leaf.taxon.label
        }
        if labels != species_set:
            continue

        try:
            topology = classify_triplet_topology_string(tree, species_triplet)
            tree_height, summary_metrics = compute_triplet_tree_metrics(
                tree,
                species_triplet=species_triplet,
                tree_height_calculation_strategy=tree_height_calculation_strategy,
                collect_summary_statistics=collect_summary_statistics,
            )
        except ValueError:
            continue

        observations.append((topology, tree_height, summary_metrics))

    return observations


def _parse_dendropy(newick):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def dendropy_observation(newick, triplet, strategy, collect):
    """Run the DendroPy extract-then-measure reference path."""
    subtree = extract_triplet_subtree(_parse_dendropy(newick), triplet)
    if subtree is None:
        return None
    return observation_from_subtree(subtree, triplet, strategy, collect)


def biopython_observation(newick, triplet, strategy, collect):
    """Measure a triplet with Bio.Phylo alone, straight from the definitions.

    Nothing here touches the package's tree code. The sister pair is the one
    pair of the three whose common ancestor is not the common ancestor of all
    three; a triplet whose three pairs share one ancestor is a polytomy and
    yields nothing, as does one with a taxon absent from the tree. Every
    distance is a Bio.Phylo path sum from the triplet's own common ancestor,
    which is the root of the subtree the other two paths would extract.
    """
    tree = Phylo.read(StringIO(newick), "newick")
    leaf = {terminal.name: terminal for terminal in tree.get_terminals()}
    if any(label not in leaf for label in triplet):
        return None

    a, b, c = triplet
    ancestor = tree.common_ancestor(leaf[a], leaf[b], leaf[c])
    pair_ancestor = {
        frozenset((a, b)): tree.common_ancestor(leaf[a], leaf[b]),
        frozenset((a, c)): tree.common_ancestor(leaf[a], leaf[c]),
        frozenset((b, c)): tree.common_ancestor(leaf[b], leaf[c]),
    }
    below = [pair for pair, node in pair_ancestor.items() if node is not ancestor]
    if len(below) != 1:
        return None
    sisters = below[0]
    topology = {
        frozenset((a, b)): TOPOLOGY_AB,
        frozenset((a, c)): TOPOLOGY_AC,
        frozenset((b, c)): TOPOLOGY_BC,
    }[sisters]

    depth = {label: ancestor.distance(leaf[label]) for label in triplet}
    internal_branch = ancestor.distance(pair_ancestor[sisters])
    left, right = tuple(sisters)
    sister_distance = tree.distance(leaf[left], leaf[right])
    avg_tree_height = sum(depth.values()) / 3.0

    height = {
        "AVG": avg_tree_height,
        "A": depth[a],
        "B": depth[b],
        "C": depth[c],
        "SIS": sister_distance,
        "INT": internal_branch,
    }[strategy]
    metrics = None
    if collect:
        metrics = {
            "avg_tree_height": avg_tree_height,
            "internal_branch": internal_branch,
            "sister_distance": sister_distance,
        }
    return (topology, height, metrics)
