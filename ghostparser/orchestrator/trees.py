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
# label -- whitespace and the punctuation Bio.Phylo or DendroPy trips on -- plus
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
    """Find the node where a set of outgroups parts from every other taxon.

    Read unrooted, the outgroups root the tree when one branch -- or, inside
    a polytomy, one node -- parts all of them from every other taxon, in
    whatever orientation the file wrote it. Every outgroup-free subtree then
    hangs off the same node, and rooting there puts every outgroup on one
    side of it.

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


def _root_tree_on_outgroup(tree, outgroup_taxa):
    """Root a tree where the outgroups branch off and prune them.

    Read unrooted, the tree must have one branch -- or, inside a polytomy,
    one node -- where the outgroups part from every other taxon. The tree is
    rooted there, whatever orientation the file wrote it in, and the outgroup
    side is cut away.

    Args:
        tree: A ``Bio.Phylo`` tree object; rerooted in place.
        outgroup_taxa: Iterable of outgroup taxon names.

    Returns:
        A tuple ``(pruned_tree, excluded_taxa, missing_taxa, ingroup_taxa)``:
        the rooted ingroup tree, the outgroup names found and pruned, the
        outgroup names absent from the tree, and the ingroup names.

    Raises:
        OutgroupRootingError: If no outgroup is in the tree, if every taxon
            is an outgroup, or if the outgroups branch off at more than one
            point so the other taxa fall into groups with outgroups between
            them; the message names those groups.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_set = set(outgroup_taxa)
    missing = outgroup_set - tree_taxa
    present = outgroup_set & tree_taxa
    if not present:
        raise OutgroupRootingError(
            "Could not root the species tree: none of the outgroup taxa are in "
            f"it (missing: {', '.join(sorted(missing))})"
        )
    ingroup_taxa = tree_taxa - present
    if not ingroup_taxa:
        raise OutgroupRootingError(
            "Could not root the species tree: every taxon in it is an outgroup"
        )

    host, groups = _outgroup_attachment(tree, present)
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

    tree.root_with_outgroup(host)
    pruned_root = _copy_clade_for_taxa(tree.root, ingroup_taxa)
    pruned_tree = Tree(root=pruned_root, rooted=True)
    return pruned_tree, present, missing, ingroup_taxa


@dataclass(frozen=True)
class GeneTreeRooting:
    """How one gene tree was rooted on its outgroups.

    Attributes:
        tree: The rooted tree with every outgroup pruned, or ``None`` when the
            tree carries no outgroup or nothing but outgroups.
        used: The outgroups the rooting used, in listed order: the largest set
            of those present that parts from the other taxa at one point.
        tangled: The outgroups present but not used, in listed order: they
            sat among the ingroup taxa, with ingroup taxa between them and
            the outgroups used, and were pruned without rooting on them. When
            ``tree`` is ``None`` every outgroup present is listed here.
        missing: The outgroups absent from the tree.
        order_decided: ``True`` when no set of outgroups held a majority --
            another set of the same size as ``used`` also parted from the
            other taxa at one point -- so the listed order chose.
    """

    tree: object
    used: tuple
    tangled: tuple
    missing: frozenset
    order_decided: bool


def root_gene_tree(tree, outgroup_taxa):
    """Root a gene tree where the largest set of its outgroups branches off.

    Every outgroup the tree carries is tried together first. When ingroup
    taxa sit between them, the largest subset that parts from the other taxa
    at a single point roots the tree instead -- a single gene often places a
    distant outgroup somewhere inside the ingroup, and the outgroups that
    still sit together outvote it; the tangled one is pruned without being
    used. Among equally large subsets none holds a majority, and the one
    listed earliest wins, so with two outgroups apart the first listed
    decides. All outgroups are then pruned: no triplet contains one, and
    pruning a leaf changes no other taxon's rooted shape or heights.

    Args:
        tree: A ``Bio.Phylo`` tree object; rerooted in place.
        outgroup_taxa: Ordered outgroup taxon names.

    Returns:
        A :class:`GeneTreeRooting`.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_list = list(dict.fromkeys(outgroup_taxa))
    present = [outgroup for outgroup in outgroup_list if outgroup in tree_taxa]
    missing = frozenset(outgroup_list) - tree_taxa
    ingroup_taxa = tree_taxa - set(present)
    if not present or not ingroup_taxa:
        return GeneTreeRooting(None, (), tuple(present), missing, False)

    # itertools.combinations keeps the listed order, so the first subset that
    # fits at a given size is the one listed earliest.
    host = None
    for size in range(len(present), 0, -1):
        fitting = []
        for subset in combinations(present, size):
            candidate, _ = _outgroup_attachment(tree, set(subset))
            if candidate is not None:
                fitting.append((subset, candidate))
        if fitting:
            used, host = fitting[0]
            order_decided = len(fitting) > 1
            break

    tree.root_with_outgroup(host)
    pruned_root = _copy_clade_for_taxa(tree.root, ingroup_taxa)
    tangled = tuple(outgroup for outgroup in present if outgroup not in used)
    return GeneTreeRooting(
        Tree(root=pruned_root, rooted=True), used, tangled, missing, order_decided
    )


@dataclass(frozen=True)
class GeneTreeCleaning:
    """What cleaning a gene-tree file kept, dropped and rooted.

    Attributes:
        trees: The cleaned, rooted ``Bio.Phylo`` trees, in input order, with
            their outgroups pruned.
        dropped_trees: 1-based input index of each tree dropped for low
            support, mapped to its average support.
        unrootable_indices: 1-based input indices of the trees dropped for
            carrying no outgroup taxon, or nothing but outgroup taxa.
        rooted_on: Each outgroup label, in the order given, mapped to the
            number of kept trees whose rooting used it.
        tangled: Each outgroup label, in the order given, mapped to the
            number of kept trees in which it sat among the ingroup taxa and
            was pruned without being used for rooting.
        tangled_trees: Number of kept trees with at least one tangled
            outgroup.
        order_decided: Number of kept trees in which no set of outgroups held
            a majority, so the listed order chose which to root on.
    """

    trees: list
    dropped_trees: dict
    unrootable_indices: list
    rooted_on: dict
    tangled: dict
    tangled_trees: int
    order_decided: int

    @property
    def rooted_count(self):
        """Number of trees kept and rooted."""
        return len(self.trees)


def clean_and_save_gene_trees(
    input_filepath, output_filepath, outgroup_taxa, min_avg_support=0.5
):
    """Read, support-filter, root on the outgroups, standardize, and save gene trees.

    Each tree is rooted by :func:`root_gene_tree` -- where the largest set of
    its outgroups branches off, with any outgroup tangled among the ingroup
    taxa pruned unused -- and written without its outgroups. A tree carrying no
    outgroup, or nothing but outgroups, is dropped.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned gene trees are written.
        outgroup_taxa: Ordered outgroup taxon names used for rooting.
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.

    Returns:
        A :class:`GeneTreeCleaning`.
    """
    trees = read_tree_file(input_filepath)
    outgroup_list = list(outgroup_taxa)

    dropped_trees = {}
    cleaned_trees = []
    rooted_on = {outgroup: 0 for outgroup in outgroup_list}
    tangled = {outgroup: 0 for outgroup in outgroup_list}
    tangled_trees = 0
    order_decided = 0
    unrootable_indices = []

    for idx, tree in enumerate(trees, start=1):
        avg_support = calculate_average_support(tree)
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        rooting = root_gene_tree(tree, outgroup_list)
        if rooting.tree is None:
            unrootable_indices.append(idx)
            continue

        for outgroup in rooting.used:
            rooted_on[outgroup] += 1
        for outgroup in rooting.tangled:
            tangled[outgroup] += 1
        tangled_trees += bool(rooting.tangled)
        order_decided += rooting.order_decided

        standardized = standardize_tree(rooting.tree)
        cleaned_trees.append(standardized)

    write_clean_trees(cleaned_trees, output_filepath)

    return GeneTreeCleaning(
        trees=cleaned_trees,
        dropped_trees=dropped_trees,
        unrootable_indices=unrootable_indices,
        rooted_on=rooted_on,
        tangled=tangled,
        tangled_trees=tangled_trees,
        order_decided=order_decided,
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
            # output -- one row for both in the consolidation matrices, an
            # ``A=X;B=X`` triplet in the results -- with no trace of the map.
            raise ValueError(
                f"Species rename map {path} maps both {display_names[display]} and "
                f"{original} to {display}"
            )
        rename_map[original] = display
        display_names[display] = original
    return rename_map


def rename_newick_labels(newick, rename_map):
    """Map the leaf labels inside a Newick string written by this module.

    The run measures in the trees' own labels and the display names enter only
    at the outputs, so the ``species_tree`` column is renamed as text rather
    than by reparsing every triplet's subtree. Labels the display name forces
    into quotes are quoted as :func:`_newick_label` writes them.

    Args:
        newick: A Newick string as :func:`format_newick_with_precision` or
            :func:`_format_triplet_subtree_newick` produce it -- no comments
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
    """Extract the subtree spanning only the triplet taxa.

    Reference implementation. The pipeline reads triplets out of a cached
    :class:`~.triplet_geometry.TripletGeometry` instead of copying a subtree per
    (triplet, gene tree) pair; this is kept as the independent second opinion the
    parity tests check that path against.

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


def _format_triplet_subtree_newick(shape, label_of, decimal_places=10):
    """Write the Newick for a triplet's induced subtree from its shape.

    Reproduces what serializing a copied-out subtree produces, including the
    branch-length formatting and the order the children are listed in, so the
    ``species_tree`` column is unchanged by computing the shape directly.

    Args:
        shape: The triplet's :class:`~.triplet_geometry.TripletSubtreeShape`.
        label_of: Mapping of taxon position to label.
        decimal_places: Number of decimal places for branch lengths.

    Returns:
        The Newick string, terminated with ``;``.
    """

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
            shape, label_of
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
