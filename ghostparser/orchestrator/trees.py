"""Tree and triplet preprocessing: cleaning, rooting, pruning, triplet generation,
and species-triplet normalization.
"""

import multiprocessing as mp
import os
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


def clean_and_save_trees(
    input_filepath, output_filepath, min_avg_support=0.5, rename_map=None
):
    """Read, support-filter, standardize, and save trees.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned trees are written.
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.
        rename_map: Optional display-name map applied to terminal labels before
            anything else, so every later stage and output uses those names.

    Returns:
        A tuple ``(cleaned_trees, dropped_trees)`` where ``dropped_trees`` maps
        the 1-based input index to its average support value.
    """
    trees = read_tree_file(input_filepath)
    for tree in trees:
        rename_taxa_in_tree(tree, rename_map)

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
    input_filepath, output_filepath, outgroup_taxa, min_avg_support=0.5, rename_map=None
):
    """Read, support-filter, root on outgroup, standardize, and save gene trees.

    Args:
        input_filepath: Path to the input Newick file.
        output_filepath: Path where cleaned gene trees are written.
        outgroup_taxa: Iterable of outgroup taxon names used for rooting.
        min_avg_support: Minimum average support threshold; trees below it are
            dropped.
        rename_map: Optional display-name map applied to terminal labels before
            rooting, so ``outgroup_taxa`` must already be display names.

    Returns:
        A tuple ``(cleaned_trees, dropped_trees, rooted_count,
        missing_outgroup_indices)``. ``dropped_trees`` maps 1-based index to
        support; ``missing_outgroup_indices`` lists indices discarded for
        lacking an outgroup taxon.
    """
    trees = read_tree_file(input_filepath)
    for tree in trees:
        rename_taxa_in_tree(tree, rename_map)

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
        ValueError: If the file is malformed, or maps a label more than once,
            or maps two labels onto the same display name.
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
        if original in rename_map:
            raise ValueError(
                f"Species rename map {path} maps {original} more than once"
            )
        if display in display_names:
            # Two labels sharing a display name would collapse into duplicate
            # tree labels, which the Newick parser rejects further downstream.
            raise ValueError(
                f"Species rename map {path} maps both {display_names[display]} and "
                f"{original} to {display}"
            )
        rename_map[original] = display
        display_names[display] = original
    return rename_map


def rename_taxa_in_tree(tree, rename_map):
    """Rename a tree's terminal labels in place.

    Args:
        tree: A ``Bio.Phylo`` tree object.
        rename_map: Mapping of tree label to display name. Labels absent from
            the map keep their original name.

    Returns:
        The number of terminals renamed.
    """
    if not rename_map:
        return 0
    renamed = 0
    for terminal in tree.get_terminals():
        display = rename_map.get(terminal.name)
        if display is not None and display != terminal.name:
            terminal.name = display
            renamed += 1
    return renamed


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
        f"({label_of[first_position]}:{branch(first_edge)},"
        f"{label_of[second_position]}:{branch(second_edge)})"
        f":{branch(shape.internal_edge)}"
    )
    odd = f"{label_of[shape.odd_position]}:{branch(shape.odd_edge)}"
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
