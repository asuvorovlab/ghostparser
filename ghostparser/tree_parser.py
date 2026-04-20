"""Tree parsing and standardization module using BioPython."""

import argparse
from itertools import combinations
from multiprocessing import cpu_count
import multiprocessing as mp
from pathlib import Path
import shutil
import time

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import dendropy

from .config import (
    DEFAULT_MIN_SUPPORT_VALUE,
    DEFAULT_PARQUET_COMPRESSION,
    DEFAULT_PARQUET_PARTITIONS,
    DEFAULT_TRIPLET_OUTPUT_FORMAT,
    PARQUET_COMPRESSION_CHOICES,
    TRIPLET_IO_FORMAT_CHOICES,
)
from .triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    find_sister_pair,
    normalize_abc_from_sister_pair,
)


def read_tree_file(filepath):
    """Read Newick trees from a file.

    Args:
        filepath: Path to the Newick tree file.

    Returns:
        A list of Bio.Phylo tree objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains invalid Newick format.
    """
    try:
        trees = list(Phylo.parse(filepath, "newick"))
        if not trees:
            raise ValueError(f"Invalid Newick format in {filepath}")
        # Validate that each tree has at least one terminal (leaf node)
        for idx, tree in enumerate(trees, start=1):
            if not tree.get_terminals():
                raise ValueError(f"Invalid Newick format in {filepath}: Tree {idx} has no terminal nodes")
        return trees
    except FileNotFoundError:
        raise FileNotFoundError(f"Tree file not found: {filepath}")
    except ValueError:
        # Re-raise ValueError as-is (our own validation errors)
        raise
    except Exception as e:
        # Catch any parsing errors from BioPython
        raise ValueError(f"Invalid Newick format in {filepath}: {e}")


def _chunk_list(items, chunk_size):
    """Yield successive chunks from a list."""
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def calculate_average_support(tree):
    """Calculate the average support value for a tree.

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        The average support value across all internal nodes, or None if no support values found.
    """
    support_values = []
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.confidence is not None:
            support_values.append(clade.confidence)

    if not support_values:
        return None

    return sum(support_values) / len(support_values)


def remove_support_values(tree):
    """Remove support values (bootstrap/posterior probabilities) from internal nodes.

    This function removes confidence values that appear after internal nodes
    in Newick format trees like: (A,B)0.95:0.5 -> (A,B):0.5

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        The tree with support values removed from all internal nodes.
    """
    # Traverse all internal nodes (non-terminals) and remove confidence values
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.confidence = None
    return tree


def standardize_tree(tree):
    """Standardize a BioPython tree object.

    This function takes a BioPython tree object and standardizes it by:
    - Removing support values (bootstrap/posterior probabilities) from internal nodes
    - Ensuring consistent Newick format output
    - Preserving branch lengths and taxa names

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        A Bio.Phylo tree object (standardized).
    """
    tree = remove_support_values(tree)
    return tree


def _format_newick_with_precision_biopython(tree, decimal_places=10):
    """Format a BioPython tree to Newick string with custom decimal precision."""

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
    """Format a DendroPy tree to Newick string with custom decimal precision."""

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
    """Format a tree to Newick string with custom decimal precision for branch lengths."""
    if hasattr(tree, "seed_node"):
        return _format_newick_with_precision_dendropy(tree, decimal_places)
    return _format_newick_with_precision_biopython(tree, decimal_places)


def write_clean_trees(trees, output_filepath, decimal_places=15):
    """Write standardized trees to a file using BioPython with custom precision.

    Args:
        trees: A list of Bio.Phylo tree objects.
        output_filepath: Path where the cleaned trees will be written.
        decimal_places: Number of decimal places for branch lengths. Default: 10.
    """
    with open(output_filepath, "w") as f:
        for tree in trees:
            # Write Newick format with custom decimal precision
            newick_str = format_newick_with_precision(tree, decimal_places)
            f.write(newick_str + "\n")


def clean_and_save_trees(input_filepath, output_filepath, min_avg_support=0.5):
    """Read, standardize, and save trees from a file using BioPython.

    Trees with average support values below the threshold are filtered out.

    Args:
        input_filepath: Path to the input Newick format tree file.
        output_filepath: Path where the cleaned trees will be written.
        min_avg_support: Minimum average support value threshold. Default: 0.5.

    Returns:
        A tuple of (cleaned_trees_list, dropped_trees_info_dict).
        dropped_trees_info_dict maps tree index to average support value.
    """
    # Read the trees using BioPython
    trees = read_tree_file(input_filepath)

    # Track dropped trees
    dropped_trees = {}
    cleaned_trees = []

    # Process each tree
    for idx, tree in enumerate(trees, start=1):
        # Calculate average support before standardization
        avg_support = calculate_average_support(tree)

        # Check if tree meets minimum support threshold
        if avg_support is not None and avg_support < min_avg_support:
            dropped_trees[idx] = avg_support
            continue

        # Standardize the tree (removes support values)
        standardized = standardize_tree(tree)
        cleaned_trees.append(standardized)

    # Write to output file using BioPython
    write_clean_trees(cleaned_trees, output_filepath)

    return cleaned_trees, dropped_trees


def get_clean_filename(filepath):
    """Generate a processed filename with processed_ prefix.

    Args:
        filepath: The original file path.

    Returns:
        The path with processed_ prefixed to the original file name.
    """
    path = Path(filepath)
    return str(path.parent / f"processed_{path.name}")


def get_taxa_from_tree(tree):
    """Extract all terminal taxa names from a phylogenetic tree.

    Args:
        tree: A Bio.Phylo tree object.

    Returns:
        A sorted list of terminal taxa names.
    """
    taxa = [terminal.name for terminal in tree.get_terminals()]
    return sorted(taxa)


def _parse_outgroup_arg(outgroup_arg):
    """Parse comma-separated outgroup taxa into a list of names."""
    if isinstance(outgroup_arg, str):
        parts = [p.strip() for p in outgroup_arg.split(",")]
        return [p for p in parts if p]
    return [str(taxon).strip() for taxon in outgroup_arg if str(taxon).strip()]


def _copy_clade_for_taxa(clade, taxa_set):
    """Copy a clade while retaining only the specified taxa, collapsing unary nodes."""
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
    """Root a tree on the MRCA of outgroup taxa and prune the outgroup clade.

    Returns:
        Tuple of (pruned_tree, excluded_taxa_set, missing_taxa_set, ingroup_taxa_set).
    """
    tree_taxa = {terminal.name for terminal in tree.get_terminals()}
    outgroup_set = set(outgroup_taxa)
    missing = outgroup_set - tree_taxa
    present = [taxon for taxon in outgroup_taxa if taxon in tree_taxa]

    if not present:
        return None, set(), missing, set()

    present_terminals = [terminal for terminal in tree.get_terminals() if terminal.name in present]
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
    """Root a tree on the first outgroup taxon present.

    Returns:
        Tuple of (rooted_tree, used_outgroup, missing_taxa_set).
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


def clean_and_save_gene_trees(input_filepath, output_filepath, outgroup_taxa, min_avg_support=0.5):
    """Read, root, standardize, and save gene trees from a file using BioPython.

    Trees with average support values below the threshold are filtered out.

    Returns:
        Tuple of (cleaned_trees_list, dropped_trees_info_dict, rooted_count, missing_outgroup_indices).
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
    """Read triplet filter file with comma-separated taxa per line.

    Returns:
        Tuple of (triplets, invalid_lines) where invalid_lines is a list of
        (line_number, line_text) for lines that do not parse into 3 taxa.
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
    """Filter triplets to those fully contained in taxa_set.

    Returns:
        Tuple of (kept_triplets, skipped_triplets) where skipped_triplets is a
        list of (triplet, missing_taxa).
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
    """Generate all unique triplet combinations from taxa excluding the outgroup(s)."""
    if isinstance(outgroup, str):
        outgroup_taxa = set(_parse_outgroup_arg(outgroup))
    else:
        outgroup_taxa = set(outgroup)
    ingroup_taxa = [taxon for taxon in taxa_list if taxon not in outgroup_taxa]
    return list(combinations(ingroup_taxa, 3))


def write_triplets_to_file(triplets, output_filepath):
    """Write triplets to a file, one triplet per line."""
    with open(output_filepath, "w") as f:
        for triplet in triplets:
            f.write(",".join(triplet) + "\n")


def extract_triplet_subtree(tree, triplet_taxa):
    """Extract a subtree containing only the specified triplet taxa.

    Returns:
        A new tree containing only the specified taxa, or None if not all taxa are present.
    """
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon}
    if not set(triplet_taxa).issubset(tree_taxa):
        return None

    subtree = tree.extract_tree_with_taxa_labels(triplet_taxa)
    return subtree


def extract_triplet_gene_trees_for_triplet(gene_trees_filepath, triplet):
    """Extract gene-tree subtrees for one triplet by streaming a gene-tree file.

    Args:
        gene_trees_filepath: Path to cleaned gene trees file (one Newick tree per line).
        triplet: Taxon triplet tuple, expected as (A, B, C).

    Returns:
        List of Newick strings for extracted rooted triplet subtrees.
    """
    triplet_set = set(triplet)
    gene_trees = []

    with open(gene_trees_filepath, "r") as gene_f:
        for line in gene_f:
            newick_str = line.strip()
            if not newick_str:
                continue

            tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
            tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
            if not triplet_set.issubset(tree_taxa):
                continue

            subtree = extract_triplet_subtree(tree, triplet)
            if subtree is not None:
                gene_trees.append(format_newick_with_precision(subtree))

    return gene_trees


def _build_species_triplet_metadata(species_tree, triplets):
    """Normalize triplets to A/B/C and compute species triplet subtree Newick.

    Args:
        species_tree: Rooted species tree as a DendroPy tree.
        triplets: Iterable of triplet tuples.

    Returns:
        Tuple of:
            - normalized_triplets: list of triplets ordered as (A,B,C) where A and B are sisters
            - species_triplet_trees: dict triplet -> species subtree Newick string
            - skipped_triplets: list of triplets that could not be mapped on the species tree
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


def _build_topology_count_summary(triplet, newick_trees):
    """Return topology counts and discordant role summary for one triplet."""
    topology_counts = {topology: 0 for topology in ALL_TOPOLOGIES}

    for newick in newick_trees:
        tree = dendropy.Tree.get(data=str(newick).strip(), schema="newick", preserve_underscores=True)
        topology = classify_triplet_topology_string(tree, triplet)
        topology_counts[topology] += 1

    if topology_counts[TOPOLOGY_BC] >= topology_counts[TOPOLOGY_AC]:
        bc_role = "discordant1"
        ac_role = "discordant2"
    else:
        bc_role = "discordant2"
        ac_role = "discordant1"

    summary = (
        f"AB:{topology_counts[TOPOLOGY_AB]}/concordant,"
        f"BC:{topology_counts[TOPOLOGY_BC]}/{bc_role},"
        f"AC:{topology_counts[TOPOLOGY_AC]}/{ac_role}"
    )
    return topology_counts, summary


def _format_triplet_header(triplet, count, species_tree_newick=None, topology_summary=None):
    """Format a triplet section header with explicit ABC mapping and topology summary."""
    if species_tree_newick is None or not str(species_tree_newick).strip():
        raise ValueError(f"Missing species subtree for triplet header: {','.join(triplet)}")
    if topology_summary is None or not str(topology_summary).strip():
        raise ValueError(f"Missing topology summary for triplet header: {','.join(triplet)}")

    base = ",".join(triplet) + f"\t{count}"
    species_tree = str(species_tree_newick)
    label_map = f"[A={triplet[0]},B={triplet[1]},C={triplet[2]}]"
    return base + f"\t{species_tree}\t{label_map}\t{topology_summary}"


def _validate_species_triplet_trees(triplets, species_triplet_trees):
    """Validate species subtree mapping is present and non-empty for all triplets."""
    if species_triplet_trees is None:
        raise ValueError("species_triplet_trees is required and cannot be None")

    missing = [triplet for triplet in triplets if not str(species_triplet_trees.get(triplet, "")).strip()]
    if missing:
        preview = "; ".join(",".join(t) for t in missing[:5])
        raise ValueError(f"Missing species subtree mapping for triplets: {preview}")

    return species_triplet_trees


def process_gene_trees_for_triplets(gene_trees, triplets):
    """Process gene trees to extract subtrees for each triplet."""
    triplet_gene_trees = {triplet: [] for triplet in triplets}

    for gene_tree in gene_trees:
        for triplet in triplets:
            subtree = extract_triplet_subtree(gene_tree, triplet)
            if subtree:
                newick_str = format_newick_with_precision(subtree)
                triplet_gene_trees[triplet].append(newick_str)

    return triplet_gene_trees


def _get_mp_context(prefer_fork=False):
    """Get a multiprocessing context for worker pools.

    Args:
        prefer_fork: If True, prefer ``fork`` where available so workers can
            inherit large read-only data via copy-on-write.
    """
    if hasattr(mp, "get_context"):
        methods = mp.get_all_start_methods()
        if prefer_fork and "fork" in methods:
            return mp.get_context("fork")
        if "forkserver" in methods:
            return mp.get_context("forkserver")
        if "spawn" in methods:
            return mp.get_context("spawn")
    return mp


_GENE_TREES_LIST: list[str] = []
_CHUNK_DIR: str = ""
_SPECIES_TRIPLET_TREES: dict[tuple[str, str, str], str] = {}
_PARQUET_OUTPUT_DIR: str = ""
_PARQUET_PARTITIONS: int = DEFAULT_PARQUET_PARTITIONS
_PARQUET_COMPRESSION: str | None = DEFAULT_PARQUET_COMPRESSION


def _require_pyarrow():
    """Import and return pyarrow modules needed for parquet I/O."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError("Parquet output requires pyarrow to be installed") from exc
    return pa, pq


def _distance_to_root(node):
    """Compute root-to-node distance using edge lengths (missing lengths treated as 0)."""
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        if current.edge_length is not None:
            distance += float(current.edge_length)
        current = current.parent_node
    return distance


def _triplet_id(triplet):
    """Build stable triplet ID used across parquet tables."""
    return "|".join(triplet)


def _normalize_parquet_compression(parquet_compression):
    """Map CLI/config parquet compression values to pyarrow values."""
    if parquet_compression == "none":
        return None
    return parquet_compression


def _build_observation_metrics(tree, triplet):
    """Build cached topology and tree-height metrics for one triplet subtree."""
    leaves = {leaf.taxon.label: leaf for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label}
    if set(leaves.keys()) != set(triplet):
        raise ValueError("Triplet subtree labels do not match expected triplet")

    h_a = _distance_to_root(leaves[triplet[0]])
    h_b = _distance_to_root(leaves[triplet[1]])
    h_c = _distance_to_root(leaves[triplet[2]])
    h_avg = (h_a + h_b + h_c) / 3.0

    topology = classify_triplet_topology_string(tree, triplet)
    sister_pair = find_sister_pair(tree)
    left_label, right_label = tuple(sister_pair)
    sister_mrca = tree.mrca(taxon_labels=[left_label, right_label])
    if sister_mrca is None:
        raise ValueError("Could not determine sister-pair MRCA")
    h_int = _distance_to_root(sister_mrca)
    h_sis = _distance_to_root(leaves[left_label]) + _distance_to_root(leaves[right_label]) - 2.0 * h_int

    return {
        "topology": topology,
        "h_a": float(h_a),
        "h_b": float(h_b),
        "h_c": float(h_c),
        "h_avg": float(h_avg),
        "h_int": float(h_int),
        "h_sis": float(h_sis),
    }


def _read_gene_trees_file(gene_trees_filepath):
    """Read all gene trees from file into memory as Newick strings.
    
    Args:
        gene_trees_filepath: Path to the cleaned gene trees file (one Newick per line).
        
    Returns:
        List of Newick strings (empty lines filtered out).
    """
    trees = []
    with open(str(gene_trees_filepath), "r") as f:
        for line in f:
            newick_str = line.strip()
            if newick_str:
                trees.append(newick_str)
    return trees


def _init_triplet_chunk_worker(gene_trees_list, chunk_dir, species_triplet_trees=None):
    """Initializer for multiprocessing triplet chunk workers.
    
    Args:
        gene_trees_list: List of Newick strings (read once in main process).
        chunk_dir: Directory for writing chunk output files.
        species_triplet_trees: Dict mapping triplet -> species subtree Newick.
    """
    global _GENE_TREES_LIST, _CHUNK_DIR, _SPECIES_TRIPLET_TREES
    if gene_trees_list is not None:
        _GENE_TREES_LIST = gene_trees_list
    _CHUNK_DIR = str(chunk_dir)
    _SPECIES_TRIPLET_TREES = species_triplet_trees or {}


def _init_triplet_chunk_worker_parquet(
    gene_trees_list,
    output_dir,
    species_triplet_trees=None,
    parquet_partitions=DEFAULT_PARQUET_PARTITIONS,
    parquet_compression=DEFAULT_PARQUET_COMPRESSION,
):
    """Initializer for parquet-writing chunk workers."""
    global _GENE_TREES_LIST, _PARQUET_OUTPUT_DIR, _SPECIES_TRIPLET_TREES, _PARQUET_PARTITIONS, _PARQUET_COMPRESSION
    if gene_trees_list is not None:
        _GENE_TREES_LIST = gene_trees_list
    _PARQUET_OUTPUT_DIR = str(output_dir)
    _SPECIES_TRIPLET_TREES = species_triplet_trees or {}
    _PARQUET_PARTITIONS = parquet_partitions
    _PARQUET_COMPRESSION = _normalize_parquet_compression(parquet_compression)


def _process_triplet_chunk_stream(args):
    """Process a chunk of triplets against in-memory gene trees and write to a chunk file.
    
    Gene trees are read once in the main process and distributed to workers via
    the initializer, avoiding parallel disk I/O contention.
    """
    chunk_index, triplet_chunk = args
    triplet_results = {triplet: [] for triplet in triplet_chunk}

    for newick_str in _GENE_TREES_LIST:
        tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
        tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
        for triplet in triplet_chunk:
            if not set(triplet).issubset(tree_taxa):
                continue
            subtree = extract_triplet_subtree(tree, triplet)
            if subtree:
                triplet_results[triplet].append(format_newick_with_precision(subtree))

    total_subtrees = 0
    triplets_with_trees = 0
    output_lines = []
    for i, triplet in enumerate(triplet_chunk):
        newick_trees = triplet_results[triplet]
        count = len(newick_trees)
        _, topology_summary = _build_topology_count_summary(triplet, newick_trees)
        total_subtrees += count
        if count > 0:
            triplets_with_trees += 1

        species_tree_newick = _SPECIES_TRIPLET_TREES.get(triplet)
        output_lines.append(_format_triplet_header(triplet, count, species_tree_newick, topology_summary))
        output_lines.append("")
        for newick in newick_trees:
            output_lines.append(newick)

        if i < len(triplet_chunk) - 1:
            output_lines.append("")
            output_lines.append("=" * 60)
            output_lines.append("")

    chunk_path = Path(_CHUNK_DIR) / f"chunk_{chunk_index:06d}.txt"
    with open(chunk_path, "w") as out_f:
        out_f.write("\n".join(output_lines))
        if output_lines:
            out_f.write("\n")

    return total_subtrees, triplets_with_trees


def _process_triplet_chunk_parquet(args):
    """Process one chunk and write parquet triplet/observation parts."""
    pa, pq = _require_pyarrow()

    chunk_index, triplet_chunk = args
    triplet_rows = []
    obs_rows_by_partition: dict[int, list[dict]] = {}
    total_subtrees = 0
    triplets_with_trees = 0

    triplet_results = {triplet: [] for triplet in triplet_chunk}
    topology_counts_by_triplet = {
        triplet: {TOPOLOGY_AB: 0, TOPOLOGY_BC: 0, TOPOLOGY_AC: 0}
        for triplet in triplet_chunk
    }

    for newick_str in _GENE_TREES_LIST:
        tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
        tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
        for triplet in triplet_chunk:
            if not set(triplet).issubset(tree_taxa):
                continue
            subtree = extract_triplet_subtree(tree, triplet)
            if not subtree:
                continue

            subtree_newick = format_newick_with_precision(subtree)
            metrics = _build_observation_metrics(subtree, triplet)
            triplet_results[triplet].append(subtree_newick)
            topology_counts_by_triplet[triplet][metrics["topology"]] += 1

            triplet_id = _triplet_id(triplet)
            partition_id = abs(hash(triplet_id)) % max(1, _PARQUET_PARTITIONS)
            obs_rows_by_partition.setdefault(partition_id, []).append(
                {
                    "triplet_id": triplet_id,
                    "gene_tree_newick": subtree_newick,
                    "topology": metrics["topology"],
                    "h_a": metrics["h_a"],
                    "h_b": metrics["h_b"],
                    "h_c": metrics["h_c"],
                    "h_avg": metrics["h_avg"],
                    "h_int": metrics["h_int"],
                    "h_sis": metrics["h_sis"],
                    "schema_version": 1,
                }
            )

    for triplet in triplet_chunk:
        count = len(triplet_results[triplet])
        total_subtrees += count
        if count > 0:
            triplets_with_trees += 1

        topology_counts = topology_counts_by_triplet[triplet]
        n_ab = int(topology_counts.get(TOPOLOGY_AB, 0))
        n_bc = int(topology_counts.get(TOPOLOGY_BC, 0))
        n_ac = int(topology_counts.get(TOPOLOGY_AC, 0))
        if n_bc >= n_ac:
            bc_role = "discordant1"
            ac_role = "discordant2"
        else:
            bc_role = "discordant2"
            ac_role = "discordant1"

        species_tree_newick = _SPECIES_TRIPLET_TREES.get(triplet)
        triplet_rows.append(
            {
                "triplet_id": _triplet_id(triplet),
                "A": triplet[0],
                "B": triplet[1],
                "C": triplet[2],
                "count": count,
                "species_tree": species_tree_newick,
                "n_ab": n_ab,
                "n_bc": n_bc,
                "n_ac": n_ac,
                "bc_role": bc_role,
                "ac_role": ac_role,
                "schema_version": 1,
            }
        )

    triplet_dir = Path(_PARQUET_OUTPUT_DIR) / "triplets"
    triplet_dir.mkdir(parents=True, exist_ok=True)
    triplet_path = triplet_dir / f"part-{chunk_index:06d}.parquet"
    triplet_table = pa.Table.from_pylist(triplet_rows)
    pq.write_table(triplet_table, triplet_path, compression=_PARQUET_COMPRESSION)

    obs_root = Path(_PARQUET_OUTPUT_DIR) / "observations"
    for partition_id, rows in obs_rows_by_partition.items():
        partition_dir = obs_root / f"partition_id={partition_id:04d}"
        partition_dir.mkdir(parents=True, exist_ok=True)
        obs_path = partition_dir / f"part-{chunk_index:06d}.parquet"
        obs_table = pa.Table.from_pylist(rows)
        pq.write_table(obs_table, obs_path, compression=_PARQUET_COMPRESSION)

    return total_subtrees, triplets_with_trees


def write_triplet_gene_trees_parquet_multiprocess(
    triplets,
    gene_trees_filepath,
    output_path,
    species_triplet_trees=None,
    parquet_partitions=DEFAULT_PARQUET_PARTITIONS,
    parquet_compression=DEFAULT_PARQUET_COMPRESSION,
    use_multiprocessing=True,
    processes=None,
    chunksize=None,
):
    """Write triplet gene trees to a partitioned parquet dataset."""
    _require_pyarrow()

    if parquet_partitions < 1:
        raise ValueError("parquet_partitions must be >= 1")

    output_dir = Path(output_path)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not triplets:
        return 0, 0, 0

    species_triplet_trees = _validate_species_triplet_trees(triplets, species_triplet_trees)

    worker_count = _calculate_worker_count(
        len(triplets),
        use_multiprocessing=use_multiprocessing,
        processes=processes,
    )

    if chunksize is None:
        chunksize = max(1, len(triplets) // (worker_count * 4))

    triplet_chunks = list(_chunk_list(triplets, chunksize))

    temp_gene_file = None
    if not isinstance(gene_trees_filepath, (str, Path)):
        temp_gene_file = output_dir / "_gene_trees_tmp.txt"
        temp_gene_file.write_text("\n".join(gene_trees_filepath))
        gene_trees_filepath = str(temp_gene_file)

    gene_trees_list = _read_gene_trees_file(gene_trees_filepath)
    args = [(idx, chunk) for idx, chunk in enumerate(triplet_chunks)]

    if use_multiprocessing and worker_count > 1:
        ctx = _get_mp_context(prefer_fork=True)
        init_gene_trees = gene_trees_list
        if hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork":
            global _GENE_TREES_LIST
            _GENE_TREES_LIST = gene_trees_list
            init_gene_trees = None
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_chunk_worker_parquet,
            initargs=(
                init_gene_trees,
                str(output_dir),
                species_triplet_trees,
                parquet_partitions,
                parquet_compression,
            ),
        ) as pool:
            totals = pool.map(_process_triplet_chunk_parquet, args)
    else:
        _init_triplet_chunk_worker_parquet(
            gene_trees_list,
            str(output_dir),
            species_triplet_trees,
            parquet_partitions,
            parquet_compression,
        )
        totals = [_process_triplet_chunk_parquet(item) for item in args]

    if temp_gene_file and temp_gene_file.exists():
        temp_gene_file.unlink()

    total_subtrees = sum(t[0] for t in totals)
    triplets_with_trees = sum(t[1] for t in totals)
    return total_subtrees, triplets_with_trees, worker_count


def _merge_files_with_separators(input_paths, output_path):
    """Merge input files into output file with separators between file contents."""
    wrote_any = False
    with open(output_path, "w") as out_f:
        for path in input_paths:
            content = Path(path).read_text()
            if not content:
                continue
            if wrote_any:
                out_f.write("\n" + "=" * 60 + "\n")
            out_f.write(content)
            if not content.endswith("\n"):
                out_f.write("\n")
            wrote_any = True


def _merge_chunk_files(chunk_dir, output_filepath, batch_size=1000):
    """Merge chunk files in batches to limit number of temp files."""
    chunk_paths = sorted(Path(chunk_dir).glob("chunk_*.txt"))
    aggregates = []
    remaining = []

    for batch_index in range(0, len(chunk_paths), batch_size):
        batch = chunk_paths[batch_index : batch_index + batch_size]
        if len(batch) == batch_size:
            agg_path = Path(chunk_dir) / f"aggregate_{batch_index // batch_size:06d}.txt"
            _merge_files_with_separators(batch, agg_path)
            for path in batch:
                path.unlink()
            aggregates.append(agg_path)
        else:
            remaining = batch

    final_parts = aggregates + remaining
    _merge_files_with_separators(final_parts, output_filepath)

    for path in remaining:
        path.unlink()
    for path in aggregates:
        path.unlink()


def _calculate_worker_count(total_items, use_multiprocessing, processes):
    """Calculate worker count based on multiprocessing flag and available items."""
    if not use_multiprocessing:
        return 1
    available_cpus = cpu_count()
    worker_count = processes or available_cpus
    return max(1, min(worker_count, total_items))


def write_triplet_gene_trees_multiprocess(
    triplets,
    gene_trees_filepath,
    output_filepath,
    species_triplet_trees=None,
    use_multiprocessing=True,
    processes=None,
    chunksize=None,
):
    """Write triplet gene trees using multiprocessing over triplets.

    Each worker processes a chunk of triplets and streams over the gene trees file.
    Results are written to per-chunk files and merged in order to avoid high memory use.

    Args:
        triplets: List of triplet tuples to process.
        gene_trees_filepath: Path to the cleaned gene trees file.
        output_filepath: Output file path for results.
        species_triplet_trees: Dict mapping triplet -> species subtree Newick (required).
        use_multiprocessing: Whether to use multiprocessing (default: True).
                           If False, processes sequentially on single worker.
        processes: Number of worker processes when multiprocessing is enabled.
                   Defaults to cpu_count(). Only used if use_multiprocessing=True.
        chunksize: Chunk size for triplet distribution to workers.
                   Only used if use_multiprocessing=True.

    Returns:
        Tuple of (total_subtrees, triplets_with_trees, worker_count).
    """
    with open(output_filepath, "w") as f:
        f.write("")

    if not triplets:
        return 0, 0, 0

    species_triplet_trees = _validate_species_triplet_trees(triplets, species_triplet_trees)

    worker_count = _calculate_worker_count(
        len(triplets),
        use_multiprocessing=use_multiprocessing,
        processes=processes,
    )

    if chunksize is None:
        chunksize = max(1, len(triplets) // (worker_count * 4))

    triplet_chunks = list(_chunk_list(triplets, chunksize))
    chunk_dir = Path(output_filepath).with_suffix("").parent / f".{Path(output_filepath).stem}_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    temp_gene_file = None
    if not isinstance(gene_trees_filepath, (str, Path)):
        temp_gene_file = chunk_dir / "gene_trees_tmp.txt"
        temp_gene_file.write_text("\n".join(gene_trees_filepath))
        gene_trees_filepath = str(temp_gene_file)

    # Read gene trees file ONCE in the main process
    gene_trees_list = _read_gene_trees_file(gene_trees_filepath)

    args = [(idx, chunk) for idx, chunk in enumerate(triplet_chunks)]

    if use_multiprocessing and worker_count > 1:
        ctx = _get_mp_context(prefer_fork=True)
        init_gene_trees = gene_trees_list
        if hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork":
            # Under fork, workers inherit this read-only payload via copy-on-write.
            global _GENE_TREES_LIST
            _GENE_TREES_LIST = gene_trees_list
            init_gene_trees = None
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_chunk_worker,
            initargs=(init_gene_trees, str(chunk_dir), species_triplet_trees),
        ) as pool:
            totals = pool.map(_process_triplet_chunk_stream, args)
    else:
        _init_triplet_chunk_worker(gene_trees_list, str(chunk_dir), species_triplet_trees)
        totals = [_process_triplet_chunk_stream(item) for item in args]

    _merge_chunk_files(chunk_dir, output_filepath, batch_size=1000)
    if temp_gene_file and temp_gene_file.exists():
        temp_gene_file.unlink()
    chunk_dir.rmdir()

    total_subtrees = sum(t[0] for t in totals)
    triplets_with_trees = sum(t[1] for t in totals)

    return total_subtrees, triplets_with_trees, worker_count


def write_triplet_gene_trees(triplet_gene_trees, output_filepath, species_triplet_trees=None):
    """Write triplet gene trees to a file in the specified format."""
    triplets = list(triplet_gene_trees.keys())
    species_triplet_trees = _validate_species_triplet_trees(triplets, species_triplet_trees)
    with open(output_filepath, "w") as f:
        for i, (triplet, newick_trees) in enumerate(triplet_gene_trees.items()):
            count = len(newick_trees)
            _, topology_summary = _build_topology_count_summary(triplet, newick_trees)
            species_tree_newick = species_triplet_trees.get(triplet)
            f.write(_format_triplet_header(triplet, count, species_tree_newick, topology_summary) + "\n")
            f.write("\n")

            for newick in newick_trees:
                f.write(f"{newick}\n")

            if i < len(triplet_gene_trees) - 1:
                f.write("\n" + "=" * 60 + "\n")


def write_triplet_gene_trees_streaming(
    triplets,
    gene_trees_filepath,
    output_filepath,
    species_triplet_trees=None,
):
    """Write triplet gene trees by streaming one triplet at a time.

    This avoids keeping all triplet results in memory and only stores
    the current triplet's subtrees.

    Returns:
        Tuple of (total_subtrees, triplets_with_trees).
    """
    total_subtrees = 0
    triplets_with_trees = 0
    species_triplet_trees = _validate_species_triplet_trees(triplets, species_triplet_trees)

    with open(output_filepath, "w") as out_f:
        for idx, triplet in enumerate(triplets):
            newick_trees = []
            triplet_set = set(triplet)

            with open(gene_trees_filepath, "r") as gene_f:
                for line in gene_f:
                    newick_str = line.strip()
                    if not newick_str:
                        continue
                    tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
                    tree_taxa = {taxon.label for taxon in tree.taxon_namespace if taxon.label}
                    if not triplet_set.issubset(tree_taxa):
                        continue
                    subtree = extract_triplet_subtree(tree, triplet)
                    if subtree:
                        newick_trees.append(format_newick_with_precision(subtree))

            count = len(newick_trees)
            _, topology_summary = _build_topology_count_summary(triplet, newick_trees)
            total_subtrees += count
            if count > 0:
                triplets_with_trees += 1

            species_tree_newick = species_triplet_trees.get(triplet)
            out_f.write(_format_triplet_header(triplet, count, species_tree_newick, topology_summary) + "\n")
            out_f.write("\n")
            for newick in newick_trees:
                out_f.write(f"{newick}\n")

            if idx < len(triplets) - 1:
                out_f.write("\n" + "=" * 60 + "\n")

    return total_subtrees, triplets_with_trees


class MetricsLogger:
    """Logger that writes to both stdout and a metrics file."""

    def __init__(self, metrics_filepath):
        self.metrics_filepath = metrics_filepath
        self.metrics_file = None
        self.lines = []

    def __enter__(self):
        """Context manager entry: open the metrics file."""
        self.metrics_file = open(self.metrics_filepath, "w")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: close the metrics file."""
        if self.metrics_file:
            self.metrics_file.close()
        return False

    def log(self, message):
        """Print message to both stdout and metrics file."""
        print(message)
        if self.metrics_file:
            self.metrics_file.write(message + "\n")
            self.metrics_file.flush()
        self.lines.append(message)


def main():
    """Main entry point for the tree parser module."""
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except ValueError as exc:
        print(f"Error: {exc}")
        return

    species_tree_path = Path(args.species_tree)
    gene_trees_path = Path(args.gene_trees)

    if not species_tree_path.exists():
        print(f"Error: Species tree file not found: {args.species_tree}")
        return

    if not gene_trees_path.exists():
        print(f"Error: Gene trees file not found: {args.gene_trees}")
        return

    if args.output:
        output_dir = species_tree_path.parent / args.output
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = species_tree_path.parent

    species_tree_clean = str(output_dir / f"processed_{species_tree_path.name}")
    gene_trees_clean = str(output_dir / f"processed_{gene_trees_path.name}")
    metrics_filepath = str(output_dir / "metrics.txt")

    with MetricsLogger(metrics_filepath) as metrics:
        metrics.log(f"Processing species tree: {args.species_tree}")
        metrics.log(f"Processing gene trees: {args.gene_trees}")
        outgroup_taxa = _parse_outgroup_arg(args.outgroup)
        support_threshold = (
            args.min_support_value
            if args.min_support_value is not None
            else DEFAULT_MIN_SUPPORT_VALUE
        )
        metrics.log(f"Outgroup: {', '.join(outgroup_taxa)}")
        metrics.log(f"Support threshold: {support_threshold}")
        metrics.log("")

        try:
            species_start_wall = time.time()
            species_start_cpu = time.process_time()
            species_trees, dropped_species = clean_and_save_trees(
                str(species_tree_path),
                species_tree_clean,
                min_avg_support=support_threshold,
            )
            metrics.log(f"✓ Species tree cleaned and saved to: {species_tree_clean}")
            metrics.log(f"  Processed {len(species_trees)} tree(s)")
            if dropped_species:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_species)} tree(s) with avg support < {support_threshold}:"
                )
                for idx, avg_support in dropped_species.items():
                    metrics.log(f"    - Index {idx} from {args.species_tree} (avg support: {avg_support:.4f})")

            species_trees = read_tree_file(species_tree_clean)
            species_wall = time.time() - species_start_wall
            species_cpu = time.process_time() - species_start_cpu
            metrics.log(f"  Time taken (wall): {species_wall:.2f}s")
            metrics.log(f"  Time taken (CPU): {species_cpu:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing species tree: {e}")
            return

        triplets = []
        species_triplet_trees = {}

        try:
            if species_trees:
                taxa = get_taxa_from_tree(species_trees[0])
                metrics.log(f"\n✓ Found {len(taxa)} taxa in species tree")

                pruned_tree, excluded_taxa, missing_taxa, ingroup_taxa = _root_tree_on_outgroup(
                    species_trees[0], outgroup_taxa
                )

                if missing_taxa:
                    metrics.log(
                        f"⚠ Warning: Outgroup taxa not found in species tree: {', '.join(sorted(missing_taxa))}"
                    )

                if pruned_tree is None or not ingroup_taxa:
                    metrics.log("⚠ Warning: Unable to root and prune species tree on outgroups")
                    metrics.log(f"  Available taxa: {', '.join(taxa)}")
                    return

                extra_excluded = sorted(excluded_taxa - set(outgroup_taxa))
                metrics.log(f"  Outgroup taxa excluded: {len(excluded_taxa)}")
                if extra_excluded:
                    metrics.log(
                        "⚠ Warning: "
                        + ", ".join(extra_excluded)
                        + " are also removed from the tree due to rooting based on outgroup(s) "
                        + ", ".join(outgroup_taxa)
                    )

                species_trees = [pruned_tree]
                write_clean_trees(species_trees, species_tree_clean)

                if args.triplet_filter:
                    filter_path = Path(args.triplet_filter)
                    if not filter_path.exists():
                        metrics.log(
                            f"✗ Error: Triplet filter file not found: {args.triplet_filter}"
                        )
                        return

                    raw_triplets, invalid_lines = read_triplet_filter_file(str(filter_path))
                    for line_number, raw in invalid_lines:
                        metrics.log(
                            f"⚠ Warning: Skipping invalid triplet line {line_number} in {filter_path}: {raw}"
                        )

                    filtered_triplets, skipped_triplets = filter_triplets_by_taxa(raw_triplets, set(ingroup_taxa))
                    for triplet, missing in skipped_triplets:
                        metrics.log(
                            "⚠ Warning: Skipping triplet with missing taxa: "
                            f"{','.join(triplet)} (missing: {', '.join(missing)})"
                        )

                    triplets = filtered_triplets
                    metrics.log(f"✓ Using {len(triplets)} filtered triplets from: {filter_path}")
                else:
                    triplets = generate_triplets(sorted(ingroup_taxa), [])
                    metrics.log(f"✓ Generated {len(triplets)} unique triplets")
                    metrics.log(f"  ({len(ingroup_taxa)} taxa choose 3 = {len(triplets)} combinations)")

                species_tree_newick = Path(species_tree_clean).read_text().strip()
                species_dendro_tree = dendropy.Tree.get(
                    data=species_tree_newick,
                    schema="newick",
                    preserve_underscores=True,
                )

                triplets, species_triplet_trees, skipped_species_triplets = _build_species_triplet_metadata(
                    species_dendro_tree,
                    triplets,
                )

                if skipped_species_triplets:
                    metrics.log(
                        "⚠ Warning: Skipping triplets that could not be mapped on species tree: "
                        + "; ".join(",".join(triplet) for triplet in skipped_species_triplets)
                    )

                metrics.log(f"✓ Normalized {len(triplets)} triplets to A,B,C (A and B are sisters)")
        except Exception as e:
            metrics.log(f"✗ Error generating triplets: {e}")
            return

        try:
            genes_start_wall = time.time()
            genes_start_cpu = time.process_time()
            gene_trees, dropped_genes, rooted_count, missing_root_indices = clean_and_save_gene_trees(
                str(gene_trees_path),
                gene_trees_clean,
                outgroup_taxa,
                min_avg_support=support_threshold,
            )
            metrics.log(f"\n✓ Gene trees cleaned and saved to: {gene_trees_clean}")
            metrics.log(f"  Processed {len(gene_trees)} tree(s)")
            metrics.log(f"  Rooted {rooted_count} tree(s) on outgroup taxa")
            if missing_root_indices:
                metrics.log(
                    f"  Discarded {len(missing_root_indices)} gene tree(s) without outgroup taxa"
                )
                metrics.log(
                    "⚠ Warning: Discarded gene tree(s) where outgroup could not be found: "
                    + ", ".join(str(idx) for idx in missing_root_indices)
                )
            if dropped_genes:
                metrics.log(
                    f"  ⚠ Dropped {len(dropped_genes)} tree(s) with avg support < {support_threshold}:"
                )
                for idx, avg_support in dropped_genes.items():
                    metrics.log(f"    - Index {idx} from {args.gene_trees} (avg support: {avg_support:.4f})")

            genes_wall = time.time() - genes_start_wall
            genes_cpu = time.process_time() - genes_start_cpu
            metrics.log(f"  Time taken (wall): {genes_wall:.2f}s")
            metrics.log(f"  Time taken (CPU): {genes_cpu:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing gene trees: {e}")
            return

        try:
            metrics.log("\n✓ Processing gene trees for triplets...")
            triplet_start_wall = time.time()
            triplet_start_cpu = time.process_time()

            if args.triplet_output_format == "parquet":
                triplet_output_path = str(output_dir / "unique_triplets_gene_trees.parquet")
            else:
                triplet_output_path = str(output_dir / "unique_triplets_gene_trees.txt")

            worker_count = _calculate_worker_count(
                len(triplets),
                use_multiprocessing=not args.no_multiprocessing,
                processes=args.processes,
            )
            metrics.log(f"  Workers used: {worker_count}")

            if args.triplet_output_format == "parquet":
                total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_parquet_multiprocess(
                    triplets,
                    gene_trees_clean,
                    triplet_output_path,
                    species_triplet_trees=species_triplet_trees,
                    parquet_partitions=args.parquet_partitions,
                    parquet_compression=args.parquet_compression,
                    use_multiprocessing=not args.no_multiprocessing,
                    processes=args.processes,
                )
            else:
                total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
                    triplets,
                    gene_trees_clean,
                    triplet_output_path,
                    species_triplet_trees=species_triplet_trees,
                    use_multiprocessing=not args.no_multiprocessing,
                    processes=args.processes,
                )

            metrics.log(f"✓ Triplet gene trees saved to: {triplet_output_path}")
            metrics.log(f"  Total triplets: {len(triplets)}")
            metrics.log(f"  Triplets with gene trees: {triplets_with_trees}")
            metrics.log(f"  Total subtrees extracted: {total_subtrees}")
            if triplets_with_trees > 0:
                avg_trees_per_triplet = total_subtrees / triplets_with_trees
                metrics.log(f"  Average trees per triplet: {avg_trees_per_triplet:.2f}")
            triplet_wall = time.time() - triplet_start_wall
            triplet_cpu = time.process_time() - triplet_start_cpu
            metrics.log(f"  Time taken (wall): {triplet_wall:.2f}s")
            metrics.log(f"  Time taken (CPU): {triplet_cpu:.2f}s")
        except Exception as e:
            metrics.log(f"✗ Error processing triplet gene trees: {e}")
            return

        metrics.log("\nProcessing complete!")
        metrics.log(f"\nMetrics saved to: {metrics_filepath}")


def _build_argument_parser():
    """Build tree_parser CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Ghost parser for identifying ghost introgressions in phylogenetic trees."
    )

    parser.add_argument("-st", "--species-tree-path", required=True, help="Path to the species tree file in Newick format")
    parser.add_argument("-gt", "--gene-trees-path", required=True, help="Path to the gene trees file in Newick format")
    parser.add_argument("-og", "--outgroups", required=True, help="Outgroup species identifier(s), comma-separated")
    parser.add_argument(
        "--triplet-filter",
        type=str,
        default=None,
        help="Path to triplet filter file (comma-separated taxa per line)",
    )
    parser.add_argument(
        "--output-folder",
        type=str,
        default=None,
        help="Output folder relative to input data folder (default: same folder as input data)",
    )
    parser.add_argument(
        "--triplet-output-format",
        choices=TRIPLET_IO_FORMAT_CHOICES,
        default=None,
        help=f"Triplet extraction output format (default: {DEFAULT_TRIPLET_OUTPUT_FORMAT})",
    )
    parser.add_argument(
        "--parquet-partitions",
        type=int,
        default=None,
        help=f"Number of parquet hash partitions when parquet output is used (default: {DEFAULT_PARQUET_PARTITIONS})",
    )
    parser.add_argument(
        "--parquet-compression",
        choices=PARQUET_COMPRESSION_CHOICES,
        default=None,
        help=f"Parquet compression codec (default: {DEFAULT_PARQUET_COMPRESSION})",
    )
    parser.add_argument(
        "--min-support-value",
        type=float,
        default=None,
        help=f"Support threshold for tree cleaning (default: {DEFAULT_MIN_SUPPORT_VALUE})",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes for triplet extraction (default: 0 = all cores)",
    )
    parser.add_argument(
        "--no-multiprocessing",
        action="store_true",
        help="Disable multiprocessing for triplet extraction",
    )
    return parser


def _resolve_runtime_args(args):
    """Resolve runtime arguments from CLI mode."""
    species_tree = Path(args.species_tree_path).expanduser().resolve()
    gene_trees = Path(args.gene_trees_path).expanduser().resolve()
    outgroups = _parse_outgroup_arg(args.outgroups)
    if not outgroups:
        raise ValueError("Missing required CLI argument: --outgroups")

    triplet_filter = None
    if args.triplet_filter:
        triplet_filter = str(Path(args.triplet_filter).expanduser().resolve())

    output = None
    if args.output_folder:
        output = str(Path(args.output_folder).expanduser().resolve())

    processes = args.processes if args.processes is not None else 0
    if processes < 0:
        raise ValueError("CLI argument --processes must be an integer >= 0")

    return argparse.Namespace(
        species_tree=str(species_tree),
        gene_trees=str(gene_trees),
        outgroup=outgroups,
        triplet_filter=triplet_filter,
        output=output,
        triplet_output_format=args.triplet_output_format or DEFAULT_TRIPLET_OUTPUT_FORMAT,
        parquet_partitions=args.parquet_partitions if args.parquet_partitions is not None else DEFAULT_PARQUET_PARTITIONS,
        parquet_compression=args.parquet_compression or DEFAULT_PARQUET_COMPRESSION,
        processes=processes,
        no_multiprocessing=bool(args.no_multiprocessing),
        min_support_value=args.min_support_value if args.min_support_value is not None else DEFAULT_MIN_SUPPORT_VALUE,
    )


if __name__ == "__main__":
    main()
