"""The reference the parity tests hold the cached geometry to, and the
observation builders the inference tests share.

``biopython_observation`` measures an unpruned tree from the definitions alone,
with Bio.Phylo and none of the package's tree code. No run uses it, so it
lives with the tests.
"""

from io import StringIO

from Bio import Phylo

from ghostparser.orchestrator.config import DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY
from ghostparser.triplet_utils import TOPOLOGY_AB, TOPOLOGY_AC, TOPOLOGY_BC

TRIPLET = ("A", "B", "C")
SPECIES_SUBTREE = "((A,B),C);"
CON, DIS1, DIS2 = TOPOLOGY_AB, TOPOLOGY_BC, TOPOLOGY_AC
# Every gate's threshold in the inference tests, passed explicitly so the
# derivations in TEST_IO.md, written at it, hold whatever the run defaults are.
ALPHA = 0.01
ALPHAS = {"alpha_dct": ALPHA, "alpha_ks": ALPHA, "permutation_kwargs": {"alpha": ALPHA}}


def observations(con_heights, dis1_heights, dis2_heights):
    """Build ``(topology, height, None)`` observations with explicit heights.

    Args:
        con_heights: Heights for the concordant ``((A,B),C)`` topology.
        dis1_heights: Heights for the ``((B,C),A)`` topology.
        dis2_heights: Heights for the ``((A,C),B)`` topology.

    Returns:
        The observation list, grouped by topology.
    """
    return [
        (topology, float(height), None)
        for topology, heights in ((CON, con_heights), (DIS1, dis1_heights), (DIS2, dis2_heights))
        for height in heights
    ]


def biopython_observation(newick, triplet, strategy, collect):
    """Measure a triplet with Bio.Phylo alone, straight from the definitions.

    The sister pair is the one pair of the three whose common ancestor is not
    the common ancestor of all three; a triplet whose three pairs share one
    ancestor is a polytomy and yields nothing, as does one with a taxon absent
    from the tree. Every distance is a Bio.Phylo path sum from the triplet's
    own common ancestor, which is the root of the triplet's subtree.

    Args:
        newick: The tree, as a Newick string.
        triplet: The ``(A, B, C)`` labels.
        strategy: The tree-height strategy.
        collect: Whether to return the summary metrics as well.

    Returns:
        ``(topology, height, metrics)``, ``metrics`` ``None`` unless
        ``collect``; ``None`` for a polytomy or an absent taxon.
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


def serialize_triplet_gene_trees(
    triplet,
    triplet_gene_trees,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Measure rooted triplet Newicks into observations with the reference.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        triplet_gene_trees: Rooted triplet Newick strings.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: Whether to keep each tree's summary metrics.

    Returns:
        One ``(topology, tree_height, summary_metrics)`` per resolved tree.
    """
    measured = (
        biopython_observation(
            newick, triplet, tree_height_calculation_strategy, collect_summary_statistics
        )
        for newick in triplet_gene_trees
    )
    return [observation for observation in measured if observation is not None]
