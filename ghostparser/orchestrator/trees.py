"""Tree and triplet preprocessing: cleaning, rooting, pruning, triplet generation,
and species-triplet normalization.
"""

import multiprocessing as mp
import os
import re
from dataclasses import dataclass
from itertools import combinations

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from ghostparser.triplet_utils import normalize_abc_from_sister_pair

from .triplet_geometry import (
    build_taxon_index,
    build_triplet_geometry,
    triplet_subtree_shape,
)


class MetricsLogger:
    """Logger that writes each message to both stdout and a metrics file.

    Attributes:
        metrics_filepath: Path of the metrics file to write.
        metrics_file: Open file handle while used as a context manager.
        lines: List of every logged message.
    """

    def __init__(self, metrics_filepath):
        """Initialize the logger.

        Args:
            metrics_filepath: Path where log lines are persisted.
        """
        self.metrics_filepath = metrics_filepath
        self.metrics_file = None
        self.lines = []

    def __enter__(self):
        """Open the metrics file on context entry.

        Returns:
            This ``MetricsLogger`` instance.
        """
        self.metrics_file = open(self.metrics_filepath, "w")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close the metrics file on context exit.

        Args:
            exc_type: Exception type, if any.
            exc_val: Exception value, if any.
            exc_tb: Exception traceback, if any.

        Returns:
            ``False`` so any exception is not suppressed.
        """
        if self.metrics_file:
            self.metrics_file.close()
        return False

    def log(self, message):
        """Print a message to stdout and append it to the metrics file.

        Args:
            message: Text to log.
        """
        print(message)
        if self.metrics_file:
            self.metrics_file.write(message + "\n")
            self.metrics_file.flush()
        self.lines.append(message)


def _parse_outgroup_arg(outgroup_arg):
    """Parse comma-separated or iterable outgroup taxa into a name list.

    Args:
        outgroup_arg: A comma-separated string or an iterable of taxon names.

    Returns:
        A list of non-empty outgroup taxon names.
    """
    if isinstance(outgroup_arg, str):
        parts = [p.strip() for p in outgroup_arg.split(",")]
        return [p for p in parts if p]
    return [str(taxon).strip() for taxon in outgroup_arg if str(taxon).strip()]


def read_tree_file(filepath):
    """Read and validate Newick trees from a file.

    Args:
        filepath: Path to the Newick tree file.

    Returns:
        A list of ``Bio.Phylo`` tree objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains invalid Newick or a tree has no
            terminal nodes.
    """
    try:
        trees = list(Phylo.parse(filepath, "newick"))
        if not trees:
            raise ValueError(f"Invalid Newick format in {filepath}")
        for idx, tree in enumerate(trees, start=1):
            if not tree.get_terminals():
                raise ValueError(
                    f"Invalid Newick format in {filepath}: Tree {idx} has no terminal nodes"
                )
        return trees
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filepath}")
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Invalid Newick format in {filepath}: {e}")


def calculate_average_support(tree):
    """Compute the mean internal-node support value of a tree.

    Args:
        tree: A ``Bio.Phylo`` tree object.

    Returns:
        The mean support across internal nodes, or ``None`` if no support
        values are present.
    """
    support_values = []
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.confidence is not None:
            support_values.append(clade.confidence)

    if not support_values:
        return None

    return sum(support_values) / len(support_values)


def remove_support_values(tree):
    """Strip internal-node confidence/support values in place.

    Args:
        tree: A ``Bio.Phylo`` tree object.

    Returns:
        The same tree with internal-node confidence cleared.
    """
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.confidence = None
    return tree


def standardize_tree(tree):
    """Standardize a tree by removing support values.

    Args:
        tree: A ``Bio.Phylo`` tree object.

    Returns:
        The standardized tree.
    """
    tree = remove_support_values(tree)
    return tree


# Characters a Newick reader takes for structure when they appear in a bare
# label: whitespace and the punctuation Bio.Phylo or DendroPy trips on, plus
# the quote itself. A label holding any of them is written single-quoted.
_NEWICK_QUOTED_CHARS = frozenset(" \t\r\n()[]{}':;,\"\\=")

# A leaf label in the Newick this module writes: the token after ``(`` or ``,``
# and before ``:``, ``,`` or ``)``, either bare or single-quoted with inner
# quotes doubled. Branch lengths follow ``:`` and so never match.
_NEWICK_LEAF_LABEL_RE = re.compile(r"(?<=[(,])(?:'(?:[^']|'')*'|[^(),:;]+)(?=[:,)])")


def _newick_label(label):
    """Write a taxon label as a Newick token, quoting it only when required.

    Args:
        label: The taxon label.

    Returns:
        The label as-is, or single-quoted with inner quotes doubled when it
        holds a character the format would otherwise read as structure.
    """
    if _NEWICK_QUOTED_CHARS.isdisjoint(label):
        return label
    return "'" + label.replace("'", "''") + "'"


def _format_newick_with_precision_biopython(tree, decimal_places=10):
    """Serialize a BioPython tree to Newick with fixed branch-length precision.

    Args:
        tree: A ``Bio.Phylo`` tree object.
        decimal_places: Number of decimal places for branch lengths.

    Returns:
        The Newick string.
    """

    def format_branch_length(branch_length):
        formatted = f"{branch_length:.{decimal_places}f}"
        if "." in formatted:
            formatted = formatted.rstrip("0").rstrip(".")
        return formatted

    def format_clade(clade):
        if clade.is_terminal():
            result = _newick_label(clade.name) if clade.name else ""
        else:
            children = [format_clade(c) for c in clade.clades]
            result = "(" + ",".join(children) + ")"

        if clade.branch_length is not None:
            result += f":{format_branch_length(clade.branch_length)}"

        return result

    return format_clade(tree.root) + ";"


def _format_newick_with_precision_dendropy(tree, decimal_places=10):
    """Serialize a DendroPy tree to Newick with fixed branch-length precision.

    Args:
        tree: A DendroPy tree object.
        decimal_places: Number of decimal places for branch lengths.

    Returns:
        The Newick string.
    """

    def format_branch_length(branch_length):
        formatted = f"{branch_length:.{decimal_places}f}"
        if "." in formatted:
            formatted = formatted.rstrip("0").rstrip(".")
        return formatted

    def format_node(node):
        if node.is_leaf():
            result = _newick_label(node.taxon.label) if node.taxon else ""
        else:
            children = [format_node(child) for child in node.child_node_iter()]
            result = "(" + ",".join(children) + ")"

        if node.edge_length is not None:
            result += f":{format_branch_length(node.edge_length)}"

        return result

    return format_node(tree.seed_node) + ";"


def format_newick_with_precision(tree, decimal_places=10):
    """Serialize a BioPython or DendroPy tree to Newick with fixed precision.

    Args:
        tree: A ``Bio.Phylo`` or DendroPy tree object.
        decimal_places: Number of decimal places for branch lengths.

    Returns:
        The Newick string.
    """
    if hasattr(tree, "seed_node"):
        return _format_newick_with_precision_dendropy(tree, decimal_places)
    return _format_newick_with_precision_biopython(tree, decimal_places)


def write_clean_trees(trees, output_filepath, decimal_places=15):
    """Write trees to a file as one Newick per line with fixed precision.

    Args:
        trees: Iterable of tree objects.
        output_filepath: Destination file path.
        decimal_places: Number of decimal places for branch lengths.
    """
    with open(output_filepath, "w") as f:
        for tree in trees:
            newick_str = format_newick_with_precision(tree, decimal_places)
            f.write(newick_str + "\n")


def clean_and_save_trees(input_filepath, output_filepath, min_avg_support=0.5):
    """Read, support-filter, standardize, and save trees.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned trees are written.
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.

    Returns:
        A tuple ``(cleaned_trees, dropped_trees)`` where ``dropped_trees`` maps
        the 1-based input index to its average support value.
    """
    trees = read_tree_file(input_filepath)

    dropped_trees = {}
    cleaned_trees = []

    for idx, tree in enumerate(trees, start=1):
        avg_support = calculate_average_support(tree)

        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        standardized = standardize_tree(tree)
        cleaned_trees.append(standardized)

    write_clean_trees(cleaned_trees, output_filepath)

    return cleaned_trees, dropped_trees


def get_taxa_from_tree(tree):
    """Extract sorted terminal taxa names from a tree.

    Args:
        tree: A ``Bio.Phylo`` tree object.

    Returns:
        A sorted list of terminal taxon names.
    """
    taxa = [terminal.name for terminal in tree.get_terminals()]
    return sorted(taxa)


def _copy_clade_for_taxa(clade, taxa_set):
    """Copy a clade retaining only taxa in the set, collapsing unary nodes.

    Args:
        clade: The ``Bio.Phylo`` clade to copy.
        taxa_set: Set of taxon names to retain.

    Returns:
        A new ``Clade`` containing only the retained taxa, or ``None`` if none
        of the clade's taxa are retained.
    """
    if clade.is_terminal():
        if clade.name in taxa_set:
            return Clade(branch_length=clade.branch_length, name=clade.name)
        return None

    new_children = []
    for child in clade.clades:
        copied = _copy_clade_for_taxa(child, taxa_set)
        if copied is not None:
            new_children.append(copied)

    if not new_children:
        return None

    if len(new_children) == 1:
        only_child = new_children[0]
        if clade.branch_length is not None:
            if only_child.branch_length is None:
                only_child.branch_length = clade.branch_length
            else:
                only_child.branch_length += clade.branch_length
        return only_child

    new_clade = Clade(branch_length=clade.branch_length, clades=new_children)
    new_clade.confidence = clade.confidence
    return new_clade


class OutgroupRootingError(ValueError):
    """The outgroups do not root the tree.

    Attributes:
        separated_groups: When the outgroups do not branch off the rest of the
            tree at a single point, the groups the other taxa fall into,
            largest first, each a sorted tuple of names; empty for every other
            cause.
    """

    def __init__(self, message, separated_groups=()):
        super().__init__(message)
        self.separated_groups = tuple(separated_groups)


# How many of a group's taxa the rooting error spells out before "+N more".
_GROUP_NAMES_SHOWN = 20


def _describe_taxon_group(taxa):
    """Write a group of taxa as ``n taxa: a, b, c (+N more)`` for a message."""
    names = sorted(taxa)
    shown = ", ".join(names[:_GROUP_NAMES_SHOWN])
    if len(names) > _GROUP_NAMES_SHOWN:
        shown += f" (+{len(names) - _GROUP_NAMES_SHOWN} more)"
    return f"{len(names)} {'taxon' if len(names) == 1 else 'taxa'}: {shown}"


def _outgroup_attachment(tree, present):
    """Find the node where a set of outgroups parts from every other taxon,
    reading the tree unrooted.

    Args:
        tree: A ``Bio.Phylo`` tree object, in any orientation.
        present: Non-empty set of outgroup names, all of them in the tree and
            not every taxon of it.

    Returns:
        A tuple ``(host, groups)``: ``groups`` lists every outgroup-free
        subtree as ``(host_clade, taxa)``, and ``host`` is their shared host
        clade, or ``None`` when they hang off different nodes so that other
        taxa sit between the outgroups.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}

    # Outgroup leaves under each clade. A clade's branch lies on a path
    # between two outgroups exactly when it holds some but not all of them.
    outgroup_counts = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            outgroup_counts[clade] = int(clade.name in present)
        else:
            outgroup_counts[clade] = sum(
                outgroup_counts[child] for child in clade.clades
            )
    total = len(present)

    # The lowest clade holding every present outgroup. Every outgroup-free
    # subtree hangs by a single branch off that clade or off a clade below it
    # whose branch joins outgroups; the clade the branch leaves from is the
    # subtree's host. The taxa outside the lowest clade form one such subtree,
    # hosted there. The outgroups root the tree only when every host is the
    # same node: rooting there puts all of them on one side of it.
    mrca = tree.root
    while True:
        below = [child for child in mrca.clades if outgroup_counts[child] == total]
        if not below:
            break
        mrca = below[0]
    groups = []
    outside = tree_taxa - {terminal.name for terminal in mrca.get_terminals()}
    if outside:
        groups.append((mrca, outside))
    for clade in mrca.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        if clade is not mrca and not 0 < outgroup_counts[clade] < total:
            continue
        for child in clade.clades:
            if outgroup_counts[child] == 0:
                groups.append(
                    (clade, {terminal.name for terminal in child.get_terminals()})
                )

    hosts = {id(host) for host, _ in groups}
    return (groups[0][0] if len(hosts) == 1 else None), groups


def has_branch_lengths(tree):
    """Tell whether every branch of a tree, the root's own aside, has a length.

    Args:
        tree: A ``Bio.Phylo`` tree.

    Returns:
        ``True`` when no branch below the root lacks a length.
    """
    return all(
        clade.branch_length is not None
        for clade in tree.find_clades()
        if clade is not tree.root
    )


@dataclass(frozen=True)
class SpeciesTreeRooting:
    """How the species tree was rooted on its outgroups.

    Attributes:
        tree: The rooted ingroup tree with every outgroup pruned.
        ranked: The outgroups in the tree, farthest from the ingroup root
            first, or in listed order when the tree lacks some branch length.
        distances: Each outgroup in the tree mapped to its path length from
            the ingroup root, farthest first, or ``None`` when the tree lacks
            some branch length.
        missing: The outgroup names absent from the tree, in listed order.
        ingroup: The ingroup taxon names.
    """

    tree: object
    ranked: tuple
    distances: dict | None
    missing: tuple
    ingroup: frozenset

    @property
    def outgroup_order(self):
        """Every outgroup in species-tree rank, as gene trees fall back on it:
        those in the species tree first, then those it lacks in listed
        order."""
        return self.ranked + self.missing


def _outgroup_distances(tree, ingroup_taxa, outgroups):
    """Measure each outgroup's path length from the ingroup root, farthest
    first.

    Args:
        tree: A ``Bio.Phylo`` tree rooted where its outgroups branch off,
            with a length on every branch.
        ingroup_taxa: Set of ingroup taxon names.
        outgroups: The outgroups in the tree, in listed order, which breaks
            ties.

    Returns:
        A dict of outgroup name to path length, farthest first.
    """
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    ingroup_root = tree.common_ancestor([terminals[name] for name in ingroup_taxa])
    # The rooting leaves every ingroup taxon on one side of the root: either
    # the root itself, when several ingroup groups hang from it, or a single
    # child of it, whose own branch then starts every outgroup's path.
    stem = [] if ingroup_root is tree.root else [ingroup_root]
    distances = {
        outgroup: sum(
            clade.branch_length for clade in stem + tree.get_path(terminals[outgroup])
        )
        for outgroup in outgroups
    }
    listed = {outgroup: index for index, outgroup in enumerate(outgroups)}
    ranked = sorted(outgroups, key=lambda name: (-distances[name], listed[name]))
    return {outgroup: distances[outgroup] for outgroup in ranked}


def _root_tree_on_outgroup(tree, outgroup_taxa):
    """Root a tree where its outgroups branch off, in whatever orientation the
    file was written, rank the outgroups by distance and prune them.

    Args:
        tree: A ``Bio.Phylo`` tree object; rerooted in place.
        outgroup_taxa: Ordered outgroup taxon names; the order ranks them
            when the tree lacks some branch length and breaks distance ties.

    Returns:
        A :class:`SpeciesTreeRooting`.

    Raises:
        OutgroupRootingError: If no outgroup is in the tree, if every taxon
            is an outgroup, or if the outgroups branch off at more than one
            point so the other taxa fall into groups with outgroups between
            them; the message names those groups.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_list = list(dict.fromkeys(outgroup_taxa))
    missing = tuple(outgroup for outgroup in outgroup_list if outgroup not in tree_taxa)
    present = [outgroup for outgroup in outgroup_list if outgroup in tree_taxa]
    if not present:
        raise OutgroupRootingError(
            "Could not root the species tree: none of the outgroup taxa are in "
            f"it (missing: {', '.join(sorted(missing))})"
        )
    ingroup_taxa = tree_taxa - set(present)
    if not ingroup_taxa:
        raise OutgroupRootingError(
            "Could not root the species tree: every taxon in it is an outgroup"
        )

    host, groups = _outgroup_attachment(tree, set(present))
    if host is None:
        separated = sorted(
            (tuple(sorted(taxa)) for _, taxa in groups),
            key=lambda group: (-len(group), group),
        )
        raise OutgroupRootingError(
            "Could not root the species tree on the outgroups "
            f"({', '.join(sorted(present))}): they do not branch off the rest "
            f"of the tree at a single point, so the other taxa fall into "
            f"{len(separated)} groups with outgroups between them: "
            + "; ".join(f"[{_describe_taxon_group(group)}]" for group in separated)
            + ". Add every group that is an outgroup to the outgroup list, or "
            "correct the species tree, and rerun.",
            separated_groups=separated,
        )

    # Rerooting turns a missing branch length into 0, so check before it.
    has_lengths = has_branch_lengths(tree)
    lengthless = all(
        clade.branch_length is None
        for clade in tree.find_clades()
        if clade is not tree.root
    )
    tree.root_with_outgroup(host)
    if lengthless:
        for clade in tree.find_clades():
            clade.branch_length = None
    distances = (
        _outgroup_distances(tree, ingroup_taxa, present) if has_lengths else None
    )
    pruned_root = _copy_clade_for_taxa(tree.root, ingroup_taxa)
    return SpeciesTreeRooting(
        tree=Tree(root=pruned_root, rooted=True),
        ranked=tuple(distances) if has_lengths else tuple(present),
        distances=distances,
        missing=missing,
        ingroup=frozenset(ingroup_taxa),
    )


@dataclass(frozen=True)
class GeneTreeRooting:
    """How one gene tree was rooted on its outgroups.

    Attributes:
        tree: The rooted tree with every outgroup pruned, or ``None`` when the
            tree carries no outgroup or nothing but outgroups.
        present: The outgroups the tree carries, in species-tree rank.
        farthest: The outgroup farthest from the ingroup taxa, or ``None``
            when ``tree`` is ``None``.
        used: ``farthest`` and the other outgroups outside the ingroup, in
            species-tree rank; the tree is rooted at their common ancestor.
        tangled: The outgroups that sat among the ingroup taxa once the tree
            was rooted on ``farthest``, in species-tree rank; pruned unused.
        missing: The outgroups absent from the tree.
    """

    tree: object
    present: tuple
    farthest: object
    used: tuple
    tangled: tuple
    missing: frozenset


def _ingroup_distance_sums(tree, outgroups, ingroup_taxa):
    """Sum each outgroup's path lengths to every ingroup taxon, reading the
    tree unrooted.

    Args:
        tree: A ``Bio.Phylo`` tree; a missing branch length counts as 0.
        outgroups: Outgroup names in the tree.
        ingroup_taxa: Set of ingroup taxon names.

    Returns:
        A dict of outgroup name to the summed path length.
    """
    # A branch lies on the path from an outgroup to every ingroup taxon on
    # its far side: those below it, or, when the outgroup is below it, those
    # above it. Counting them once per branch sums every path in one pass.
    ingroup_counts = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            ingroup_counts[clade] = int(clade.name in ingroup_taxa)
        else:
            ingroup_counts[clade] = sum(ingroup_counts[child] for child in clade.clades)
    total = len(ingroup_taxa)
    branches = [clade for clade in ingroup_counts if clade is not tree.root]
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    sums = {}
    for outgroup in outgroups:
        above = set(tree.get_path(terminals[outgroup]))
        sums[outgroup] = sum(
            (clade.branch_length or 0.0)
            * (
                total - ingroup_counts[clade]
                if clade in above
                else ingroup_counts[clade]
            )
            for clade in branches
        )
    return sums


def root_gene_tree(tree, outgroup_taxa):
    """Root a gene tree at the common ancestor of its farthest outgroup and
    the other outgroups outside the ingroup, and prune every outgroup.

    Args:
        tree: A ``Bio.Phylo`` tree object; rerooted in place. A missing
            branch length counts as 0.
        outgroup_taxa: Outgroup taxon names in species-tree rank (see
            :attr:`SpeciesTreeRooting.outgroup_order`), which breaks exact
            ties, as between the outgroups of a tree lacking every length.

    Returns:
        A :class:`GeneTreeRooting`.
    """
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    outgroup_list = list(dict.fromkeys(outgroup_taxa))
    present = tuple(outgroup for outgroup in outgroup_list if outgroup in terminals)
    missing = frozenset(outgroup_list) - terminals.keys()
    ingroup_taxa = terminals.keys() - set(present)
    if not present or not ingroup_taxa:
        return GeneTreeRooting(None, present, None, (), (), missing)

    # The farthest outgroup has the longest mean path to the ingroup taxa,
    # which needs no root; for outgroups outside the ingroup it ranks them
    # as the species tree's distance from the ingroup root does. max() keeps
    # the first of equals, so the species-tree rank breaks ties.
    if len(present) == 1:
        farthest = present[0]
    else:
        sums = _ingroup_distance_sums(tree, present, ingroup_taxa)
        farthest = max(present, key=sums.__getitem__)

    # Rooted on the farthest, an outgroup under a child of the ingroup's
    # lowest common ancestor that also holds ingroup taxa sits among them.
    # One outside that ancestor, or under a child holding no ingroup taxon
    # (possible only at a polytomy), sits between the ingroup and the
    # farthest.
    tree.root_with_outgroup(terminals[farthest])
    ingroup_counts = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            ingroup_counts[clade] = int(clade.name in ingroup_taxa)
        else:
            ingroup_counts[clade] = sum(ingroup_counts[child] for child in clade.clades)
    total = len(ingroup_taxa)
    ingroup_mrca = tree.root
    while True:
        below = [
            child for child in ingroup_mrca.clades if ingroup_counts[child] == total
        ]
        if not below:
            break
        ingroup_mrca = below[0]
    among_ingroup = {
        terminal.name
        for child in ingroup_mrca.clades
        if ingroup_counts[child]
        for terminal in child.get_terminals()
    }
    tangled = tuple(outgroup for outgroup in present if outgroup in among_ingroup)
    used = tuple(outgroup for outgroup in present if outgroup not in among_ingroup)

    # The used outgroups part from the other taxa at one node, their common
    # ancestor read unrooted, where the species tree is rooted too. For the
    # ingroup it is the same root as the farthest outgroup's own branch.
    host, _ = _outgroup_attachment(tree, set(used))
    tree.root_with_outgroup(host)
    pruned_root = _copy_clade_for_taxa(tree.root, ingroup_taxa)
    return GeneTreeRooting(
        Tree(root=pruned_root, rooted=True),
        present,
        farthest,
        used,
        tangled,
        missing,
    )


@dataclass(frozen=True)
class GeneTreeCleaning:
    """What cleaning a gene-tree file kept, dropped and rooted.

    Attributes:
        trees: The cleaned, rooted ``Bio.Phylo`` trees, in input order, with
            their outgroups pruned.
        dropped_trees: 1-based input index of each tree dropped for low
            support, mapped to its average support.
        missing_length_indices: 1-based input indices of the kept trees
            lacking some branch length, each read as 0.
        unrootable_indices: 1-based input indices of the trees dropped for
            carrying no outgroup taxon, or nothing but outgroup taxa.
        farthest: Each outgroup label, in the order given, mapped to the
            number of kept trees in which it was the farthest outgroup.
        rooted_on: Each outgroup label, in the order given, mapped to the
            number of kept trees whose rooting used it.
        tangled: Each outgroup label, in the order given, mapped to the
            number of kept trees in which it sat among the ingroup taxa and
            was pruned without being used for rooting.
        tangled_trees: Number of kept trees with at least one tangled
            outgroup.
    """

    trees: list
    dropped_trees: dict
    missing_length_indices: list
    unrootable_indices: list
    farthest: dict
    rooted_on: dict
    tangled: dict
    tangled_trees: int

    @property
    def rooted_count(self):
        """Number of trees kept and rooted."""
        return len(self.trees)


def clean_and_save_gene_trees(
    input_filepath, output_filepath, outgroup_taxa, min_avg_support=0.5
):
    """Read, support-filter, root, standardize and save the gene trees,
    dropping any that carries no outgroup or nothing but outgroups.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned gene trees are written.
        outgroup_taxa: Outgroup taxon names in species-tree rank (see
            :attr:`SpeciesTreeRooting.outgroup_order`).
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.

    Returns:
        A :class:`GeneTreeCleaning`.
    """
    trees = read_tree_file(input_filepath)
    outgroup_list = list(outgroup_taxa)

    dropped_trees = {}
    cleaned_trees = []
    farthest = {outgroup: 0 for outgroup in outgroup_list}
    rooted_on = {outgroup: 0 for outgroup in outgroup_list}
    tangled = {outgroup: 0 for outgroup in outgroup_list}
    tangled_trees = 0
    missing_length_indices = []
    unrootable_indices = []

    for idx, tree in enumerate(trees, start=1):
        avg_support = calculate_average_support(tree)
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue
        # Checked before rooting, which turns a missing length into 0.
        lacks_length = not has_branch_lengths(tree)
        rooting = root_gene_tree(tree, outgroup_list)
        if rooting.tree is None:
            unrootable_indices.append(idx)
            continue
        if lacks_length:
            missing_length_indices.append(idx)

        farthest[rooting.farthest] += 1
        for outgroup in rooting.used:
            rooted_on[outgroup] += 1
        for outgroup in rooting.tangled:
            tangled[outgroup] += 1
        tangled_trees += bool(rooting.tangled)

        standardized = standardize_tree(rooting.tree)
        cleaned_trees.append(standardized)

    write_clean_trees(cleaned_trees, output_filepath)

    return GeneTreeCleaning(
        trees=cleaned_trees,
        dropped_trees=dropped_trees,
        missing_length_indices=missing_length_indices,
        unrootable_indices=unrootable_indices,
        farthest=farthest,
        rooted_on=rooted_on,
        tangled=tangled,
        tangled_trees=tangled_trees,
    )


# Characters the outputs put around a display name: a tab or line break ends a
# TSV cell, a comma separates the names in the ``triplet`` column, and ``=``
# and ``;`` structure ``abc_mapping``. A name holding one would corrupt the
# column it lands in, so the map is refused up front.
_DISPLAY_NAME_DELIMITERS = {
    "\t": "a tab",
    "\n": "a line break",
    "\r": "a carriage return",
    ",": "a comma",
    ";": "a semicolon",
    "=": "an equals sign",
}


def load_species_rename_map(filepath):
    """Read a display-name map keyed by the taxon labels used in the trees.

    Accepts a YAML mapping (``.yaml``/``.yml``) or a two-column TSV, where the
    first column is the label as it appears in the species and gene trees and
    the second is the name to show in the outputs.

    Args:
        filepath: Path to the rename map file.

    Returns:
        A dict mapping tree label to display name.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is malformed, maps a label more than once,
            maps two labels onto the same display name, or gives a display
            name holding a character the outputs use as a delimiter.
    """
    path = str(filepath)
    try:
        if path.lower().endswith((".yaml", ".yml")):
            import yaml

            with open(path, "r") as handle:
                payload = yaml.safe_load(handle) or {}
            if not isinstance(payload, dict):
                raise ValueError(
                    f"Species rename map {path} must be a mapping of tree label to display name"
                )
            pairs = [(str(key), str(value)) for key, value in payload.items()]
        else:
            pairs = []
            with open(path, "r") as handle:
                for line_number, line in enumerate(handle, start=1):
                    raw = line.strip()
                    if not raw or raw.startswith("#"):
                        continue
                    columns = [part.strip() for part in raw.split("\t")]
                    columns = [part for part in columns if part]
                    if len(columns) != 2:
                        raise ValueError(
                            f"Species rename map {path} line {line_number} must have "
                            f"two tab-separated columns: {raw}"
                        )
                    pairs.append((columns[0], columns[1]))
    except FileNotFoundError:
        raise FileNotFoundError(f"Species rename map not found: {path}")

    rename_map = {}
    display_names = {}
    for original, display in pairs:
        if not original or not display:
            raise ValueError(f"Species rename map {path} has an empty label or name")
        for char, what in _DISPLAY_NAME_DELIMITERS.items():
            if char in display:
                raise ValueError(
                    f"Species rename map {path} names {original} {display!r}, which "
                    f"holds {what}; a display name cannot contain a tab, line break, "
                    "comma, semicolon or equals sign because the results TSV uses "
                    "them as delimiters"
                )
        if original in rename_map:
            raise ValueError(
                f"Species rename map {path} maps {original} more than once"
            )
        if display in display_names:
            # Two labels sharing a display name would merge two taxa in every
            # output: one row for both in the consolidation matrices, an
            # ``A=X;B=X`` triplet in the results, with no trace of the map.
            raise ValueError(
                f"Species rename map {path} maps both {display_names[display]} and "
                f"{original} to {display}"
            )
        rename_map[original] = display
        display_names[display] = original
    return rename_map


def rename_newick_labels(newick, rename_map):
    """Map the leaf labels inside a Newick string written by this module,
    quoting any display name the format cannot carry bare.

    Args:
        newick: A Newick string as :func:`format_newick_with_precision` or
            :func:`_format_triplet_subtree_newick` produce it: no comments
            or whitespace, labels bare or single-quoted.
        rename_map: Mapping of tree label to display name. Labels absent from
            the map keep their name.

    Returns:
        The Newick string with its leaf labels mapped.
    """
    if not rename_map or not newick:
        return newick

    def substitute(match):
        token = match.group(0)
        label = token[1:-1].replace("''", "'") if token.startswith("'") else token
        return _newick_label(rename_map.get(label, label))

    return _NEWICK_LEAF_LABEL_RE.sub(substitute, newick)


def rename_taxon_labels(labels, rename_map):
    """Map a sequence of taxon labels through the rename map.

    Args:
        labels: Iterable of taxon labels.
        rename_map: Mapping of tree label to display name.

    Returns:
        A list of display names, leaving unmapped labels unchanged.
    """
    if not rename_map:
        return list(labels)
    return [rename_map.get(label, label) for label in labels]


def read_triplet_filter_file(filepath):
    """Read a triplet filter file with three comma-separated taxa per line.

    Args:
        filepath: Path to the triplet filter file.

    Returns:
        A tuple ``(triplets, invalid_lines)`` where ``triplets`` is a list of
        3-tuples and ``invalid_lines`` is a list of ``(line_number, text)`` for
        lines that do not parse into three taxa.
    """
    triplets = []
    invalid_lines = []
    with open(filepath, "r") as f:
        for line_number, line in enumerate(f, start=1):
            raw = line.strip()
            if not raw:
                continue
            parts = [part.strip() for part in raw.split(",")]
            parts = [part for part in parts if part]
            if len(parts) != 3:
                invalid_lines.append((line_number, raw))
                continue
            triplets.append(tuple(parts))
    return triplets, invalid_lines


def read_species_filter_file(filepath):
    """Read a species filter file naming taxa, one or more comma-separated per line.

    Args:
        filepath: Path to the species filter file.

    Returns:
        The taxon names in file order, blank entries dropped and repeats
        removed.
    """
    species = []
    seen = set()
    with open(filepath, "r") as f:
        for line in f:
            for part in line.split(","):
                name = part.strip()
                if name and name not in seen:
                    seen.add(name)
                    species.append(name)
    return species


def filter_triplets_by_taxa(triplets, taxa_set):
    """Keep triplets whose taxa are all contained in the set.

    Args:
        triplets: Iterable of triplet tuples.
        taxa_set: Set of allowed taxon names.

    Returns:
        A tuple ``(kept_triplets, skipped_triplets)`` where ``skipped_triplets``
        is a list of ``(triplet, missing_taxa)``.
    """
    kept_triplets = []
    skipped_triplets = []
    for triplet in triplets:
        missing = sorted(set(triplet) - taxa_set)
        if missing:
            skipped_triplets.append((triplet, missing))
            continue
        kept_triplets.append(triplet)
    return kept_triplets, skipped_triplets


def generate_triplets(taxa_list, outgroup):
    """Generate all unique ingroup triplet combinations.

    Args:
        taxa_list: Iterable of taxon names.
        outgroup: A comma-separated string or an iterable of outgroup names to
            exclude.

    Returns:
        A list of 3-tuples of ingroup taxa.
    """
    if isinstance(outgroup, str):
        outgroup_taxa = set(_parse_outgroup_arg(outgroup))
    else:
        outgroup_taxa = set(outgroup)
    ingroup_taxa = [taxon for taxon in taxa_list if taxon not in outgroup_taxa]
    return list(combinations(ingroup_taxa, 3))


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract the subtree spanning only the triplet taxa. This is the
    reference path the parity tests hold the cached geometry to; no run
    calls it.

    Args:
        tree: A DendroPy tree object.
        triplet_taxa: Iterable of the three taxon names.

    Returns:
        The extracted subtree, or ``None`` if any triplet taxon is absent.
    """
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon}
    if not set(triplet_taxa).issubset(tree_taxa):
        return None

    subtree = tree.extract_tree_with_taxa_labels(triplet_taxa)
    return subtree


def _format_triplet_subtree_newick(
    shape, label_of, decimal_places=10, with_lengths=True
):
    """Write the Newick for a triplet's induced subtree from its shape.

    Reproduces what serializing a copied-out subtree produces, including the
    branch-length formatting and the order the children are listed in, so the
    ``species_tree`` column is unchanged by computing the shape directly.

    Args:
        shape: The triplet's :class:`~.triplet_geometry.TripletSubtreeShape`.
        label_of: Mapping of taxon position to label.
        decimal_places: Number of decimal places for branch lengths.
        with_lengths: Write branch lengths; ``False`` writes the topology
            alone, for a tree that carries none.

    Returns:
        The Newick string, terminated with ``;``.
    """
    if not with_lengths:
        first_position, second_position = shape.sister_positions
        clade = (
            f"({_newick_label(label_of[first_position])},"
            f"{_newick_label(label_of[second_position])})"
        )
        odd = _newick_label(label_of[shape.odd_position])
        inner = f"{clade},{odd}" if shape.sister_clade_first else f"{odd},{clade}"
        return f"({inner});"

    def branch(length):
        formatted = f"{length:.{decimal_places}f}"
        if "." in formatted:
            formatted = formatted.rstrip("0").rstrip(".")
        return formatted

    first_position, second_position = shape.sister_positions
    first_edge, second_edge = shape.sister_edges
    clade = (
        f"({_newick_label(label_of[first_position])}:{branch(first_edge)},"
        f"{_newick_label(label_of[second_position])}:{branch(second_edge)})"
        f":{branch(shape.internal_edge)}"
    )
    odd = f"{_newick_label(label_of[shape.odd_position])}:{branch(shape.odd_edge)}"
    inner = f"{clade},{odd}" if shape.sister_clade_first else f"{odd},{clade}"
    root = "" if shape.root_edge is None else f":{branch(shape.root_edge)}"
    return f"({inner}){root};"


def _build_species_triplet_metadata(species_tree, triplets):
    """Normalize triplets to (A, B, C) with A and B sisters and record subtrees.

    Args:
        species_tree: Rooted species tree as a DendroPy tree.
        triplets: Iterable of triplet tuples.

    Returns:
        A tuple ``(normalized_triplets, species_triplet_trees,
        skipped_triplets)`` where ``normalized_triplets`` are ordered ``(A, B,
        C)`` with A and B as sisters, ``species_triplet_trees`` maps each
        triplet to its species subtree Newick, and ``skipped_triplets`` lists
        triplets that could not be mapped.
    """
    normalized_triplets = []
    species_triplet_trees = {}
    skipped_triplets = []

    # One cached pass over the species tree answers every triplet, the same way
    # the gene trees are handled.
    taxon_index = build_taxon_index(triplets)
    geometry = build_triplet_geometry(species_tree, taxon_index)
    with_lengths = any(
        node.edge_length is not None
        for node in species_tree.preorder_node_iter()
        if node is not species_tree.seed_node
    )

    seen = set()
    for triplet in triplets:
        positions = tuple(taxon_index[label] for label in triplet)
        shape = triplet_subtree_shape(geometry, positions)
        if shape is None:
            skipped_triplets.append(triplet)
            continue

        label_of = dict(zip(positions, triplet))
        sister_pair = frozenset(label_of[p] for p in shape.sister_positions)
        abc_triplet = normalize_abc_from_sister_pair(sorted(triplet), sister_pair)

        if abc_triplet in seen:
            continue

        seen.add(abc_triplet)
        normalized_triplets.append(abc_triplet)
        species_triplet_trees[abc_triplet] = _format_triplet_subtree_newick(
            shape, label_of, with_lengths=with_lengths
        )

    return normalized_triplets, species_triplet_trees, skipped_triplets


def _read_gene_trees_file(gene_trees_filepath):
    """Read all cleaned gene trees into memory as Newick strings.

    Args:
        gene_trees_filepath: Path to the cleaned gene trees file (one Newick per
            line).

    Returns:
        A list of Newick strings with blank lines dropped.
    """
    trees = []
    with open(str(gene_trees_filepath), "r") as f:
        for line in f:
            newick_str = line.strip()
            if newick_str:
                trees.append(newick_str)
    return trees


def _get_mp_context(prefer_fork=None):
    """Select a multiprocessing context, preferring fork then forkserver then spawn.

    Args:
        prefer_fork: Whether to prefer ``fork`` when available. Defaults to
            ``True`` on POSIX platforms when ``None``.

    Returns:
        A multiprocessing context (or the ``multiprocessing`` module itself).
    """
    if prefer_fork is None:
        prefer_fork = os.name == "posix"

    if hasattr(mp, "get_context"):
        methods = mp.get_all_start_methods()
        if prefer_fork and "fork" in methods:
            return mp.get_context("fork")
        if "forkserver" in methods:
            return mp.get_context("forkserver")
        if "spawn" in methods:
            return mp.get_context("spawn")
    return mp
