"""Tree and triplet preprocessing for the streaming pipeline.

Provides cleaning, rooting, pruning, triplet generation, and species-triplet
normalization, without any of the intermediate-file machinery the streaming
pass makes unnecessary.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from itertools import combinations

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from ghostparser.triplet_utils import (
    find_sister_pair,
    normalize_abc_from_sister_pair,
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
            result = clade.name or ""
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
            result = node.taxon.label if node.taxon else ""
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


def _root_tree_on_outgroup(tree, outgroup_taxa):
    """Root a tree on the outgroup MRCA and prune the outgroup clade.

    Args:
        tree: A ``Bio.Phylo`` tree object.
        outgroup_taxa: Iterable of outgroup taxon names.

    Returns:
        A tuple ``(pruned_tree, excluded_taxa, missing_taxa, ingroup_taxa)``.
        ``pruned_tree`` is ``None`` when rooting/pruning is not possible.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_set = set(outgroup_taxa)
    missing = outgroup_set - tree_taxa
    present = [taxon for taxon in outgroup_taxa if taxon in tree_taxa]

    if not present:
        return None, set(), missing, set()

    present_terminals = [
        terminal for terminal in tree.get_terminals() if terminal.name in present
    ]
    if len(present_terminals) == 1:
        mrca = present_terminals[0]
    else:
        mrca = tree.common_ancestor(*present_terminals)

    if mrca is None:
        return None, set(), missing, set()

    excluded_taxa = {terminal.name for terminal in mrca.get_terminals()}
    tree.root_with_outgroup(mrca)
    ingroup_taxa = set(tree_taxa) - excluded_taxa
    if not ingroup_taxa:
        return None, excluded_taxa, missing, set()

    pruned_root = _copy_clade_for_taxa(tree.root, ingroup_taxa)
    if pruned_root is None:
        return None, excluded_taxa, missing, set()

    pruned_tree = Tree(root=pruned_root, rooted=True)

    return pruned_tree, excluded_taxa, missing, ingroup_taxa


def _root_tree_on_any_outgroup(tree, outgroup_taxa):
    """Root a tree on the first present outgroup taxon.

    Args:
        tree: A ``Bio.Phylo`` tree object.
        outgroup_taxa: Iterable of candidate outgroup taxon names.

    Returns:
        A tuple ``(rooted_tree, used_outgroup, missing_taxa)`` where
        ``used_outgroup`` is ``None`` if no outgroup taxon is present.
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_list = list(outgroup_taxa)
    missing = set(outgroup_list) - tree_taxa

    for outgroup in outgroup_list:
        if outgroup in tree_taxa:
            terminal = next(t for t in tree.get_terminals() if t.name == outgroup)
            tree.root_with_outgroup(terminal)
            return tree, outgroup, missing

    return tree, None, missing


def clean_and_save_gene_trees(
    input_filepath, output_filepath, outgroup_taxa, min_avg_support=0.5
):
    """Read, support-filter, root on outgroup, standardize, and save gene trees.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned gene trees are written.
        outgroup_taxa: Iterable of outgroup taxon names used for rooting.
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.

    Returns:
        A tuple ``(cleaned_trees, dropped_trees, rooted_count,
        missing_outgroup_indices)``. ``dropped_trees`` maps 1-based index to
        support; ``missing_outgroup_indices`` lists indices discarded for
        lacking an outgroup taxon.
    """
    trees = read_tree_file(input_filepath)

    dropped_trees = {}
    cleaned_trees = []
    rooted_count = 0
    missing_outgroup_indices = []

    for idx, tree in enumerate(trees, start=1):
        avg_support = calculate_average_support(tree)
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        rooted_tree, used_outgroup, _ = _root_tree_on_any_outgroup(tree, outgroup_taxa)
        if used_outgroup is None:
            missing_outgroup_indices.append(idx)
            continue

        rooted_count += 1

        standardized = standardize_tree(rooted_tree)
        cleaned_trees.append(standardized)

    write_clean_trees(cleaned_trees, output_filepath)

    return cleaned_trees, dropped_trees, rooted_count, missing_outgroup_indices


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

    seen = set()
    for triplet in triplets:
        subtree = extract_triplet_subtree(species_tree, triplet)
        if subtree is None:
            skipped_triplets.append(triplet)
            continue

        try:
            labels = sorted(triplet)
            sister_pair = find_sister_pair(subtree)
            abc_triplet = normalize_abc_from_sister_pair(labels, sister_pair)
        except ValueError:
            skipped_triplets.append(triplet)
            continue

        if abc_triplet in seen:
            continue

        seen.add(abc_triplet)
        normalized_triplets.append(abc_triplet)
        species_triplet_trees[abc_triplet] = format_newick_with_precision(subtree)

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
