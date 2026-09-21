"""Compute and plot triplet tree heights by topology class.
"""

import argparse
import io
from itertools import combinations
import json
from pathlib import Path

from Bio.Phylo.BaseTree import Tree
from Bio.Phylo._io import parse as phylo_parse
from Bio.Phylo._io import read as phylo_read
from Bio.Phylo._io import write as phylo_write
import dendropy
import numpy as np

from ghostparser.orchestrator.trees import _root_tree_on_any_outgroup, _root_tree_on_outgroup


TOPOLOGY_AB = "((A,B),C)"
TOPOLOGY_BC = "((B,C),A)"
TOPOLOGY_AC = "((A,C),B)"
BRANCH_STRATEGY_CHOICES = ("avg_height", "internal_branch", "sister_distance", "external_branch")


def _triplet_taxa_labels(tree: dendropy.Tree) -> list[str]:
    """Return sorted labels for a rooted 3-tip tree."""
    labels = sorted(
        leaf.taxon.label
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None and leaf.taxon.label
    )
    if len(labels) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")
    return labels


def _find_sister_pair(tree: dendropy.Tree) -> frozenset[str]:
    """Return sister-pair labels for a rooted 3-tip tree."""
    tree.is_rooted = True
    labels = _triplet_taxa_labels(tree)
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


def _normalize_abc_from_sister_pair(labels: list[str], sister_pair: frozenset[str]) -> tuple[str, str, str]:
    """Normalize labels to (A, B, C) where A and B are sisters."""
    labels_set = set(labels)
    a_taxon, b_taxon = sorted(sister_pair)
    c_taxon = next(iter(labels_set - set(sister_pair)))
    return a_taxon, b_taxon, c_taxon


def _topology_from_sister_pair(sister_pair: frozenset[str], abc_triplet: tuple[str, str, str]) -> str:
    """Map sister pair to canonical topology string."""
    a_taxon, b_taxon, c_taxon = abc_triplet

    if sister_pair == frozenset((a_taxon, b_taxon)):
        return TOPOLOGY_AB
    if sister_pair == frozenset((b_taxon, c_taxon)):
        return TOPOLOGY_BC
    if sister_pair == frozenset((a_taxon, c_taxon)):
        return TOPOLOGY_AC

    raise ValueError("Tree taxa do not match provided ABC triplet")


def classify_triplet_topology_string(tree: dendropy.Tree, abc_triplet: tuple[str, str, str]) -> str:
    """Classify rooted triplet topology as one of the canonical strings."""
    sister_pair = _find_sister_pair(tree)
    return _topology_from_sister_pair(sister_pair, abc_triplet)


def _distance_to_root(node: dendropy.Node) -> float:
    """Compute root-to-node distance using edge lengths (missing treated as 0)."""
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node
    return distance


def _distance_to_ancestor(node: dendropy.Node, ancestor: dendropy.Node) -> float:
    """Compute node-to-ancestor path distance using edge lengths."""
    distance = 0.0
    current = node
    while current is not None and current is not ancestor:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node

    if current is not ancestor:
        raise ValueError("Specified ancestor is not on node's ancestry path")
    return distance


def compute_tree_height_statistic(
    tree: dendropy.Tree,
    strategy: str = "AVG",
    species_triplet: tuple[str, str, str] | None = None,
) -> float:
    """Compute H(T) from root-to-tip distances using GhostParser's default logic."""
    leaves = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label]
    if len(leaves) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")

    if strategy in {"A", "B", "C"}:
        if species_triplet is None:
            raise ValueError("species_triplet is required for strategies A, B, and C")

        strategy_index = {"A": 0, "B": 1, "C": 2}[strategy]
        selected_taxon_label = species_triplet[strategy_index]
        for leaf in leaves:
            if leaf.taxon.label == selected_taxon_label:
                return _distance_to_root(leaf)
        raise ValueError(f"Selected taxon {selected_taxon_label} not found in triplet tree")

    if strategy != "AVG":
        raise ValueError("Unsupported tree height strategy. Use one of: AVG, A, B, C")

    total_distance = 0.0
    for leaf in leaves:
        total_distance += _distance_to_root(leaf)
    return total_distance / 3.0


def compute_branch_metric(tree: dendropy.Tree, branch_strategy: str) -> float:
    """Compute selected branch metric on a rooted 3-tip triplet tree."""
    if branch_strategy not in BRANCH_STRATEGY_CHOICES:
        raise ValueError(
            "Unsupported branch strategy. "
            f"Choose one of: {', '.join(BRANCH_STRATEGY_CHOICES)}"
        )

    labels = _triplet_taxa_labels(tree)
    label_to_leaf = {
        leaf.taxon.label: leaf
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None and leaf.taxon.label
    }
    sister_pair = _find_sister_pair(tree)
    sister_labels = sorted(sister_pair)
    sister_mrca = tree.mrca(taxon_labels=sister_labels)
    if sister_mrca is None:
        raise ValueError("Could not determine sister MRCA for triplet tree")

    root = tree.seed_node
    if root is None:
        raise ValueError("Triplet tree is missing root node")

    if branch_strategy == "avg_height":
        return compute_tree_height_statistic(tree, strategy="AVG")

    if branch_strategy == "internal_branch":
        return _distance_to_ancestor(sister_mrca, root)

    if branch_strategy == "sister_distance":
        left = label_to_leaf[sister_labels[0]]
        right = label_to_leaf[sister_labels[1]]
        return _distance_to_ancestor(left, sister_mrca) + _distance_to_ancestor(right, sister_mrca)

    if branch_strategy == "external_branch":
        external_label = next(iter(set(labels) - sister_pair))
        external_leaf = label_to_leaf[external_label]
        return _distance_to_ancestor(external_leaf, root)

    raise ValueError(
        "Unsupported branch strategy. "
        f"Choose one of: {', '.join(BRANCH_STRATEGY_CHOICES)}"
    )


def extract_triplet_subtree(tree: dendropy.Tree, triplet_taxa: tuple[str, str, str]) -> dendropy.Tree | None:
    """Extract subtree containing only the specified triplet taxa."""
    tree_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon}
    if not set(triplet_taxa).issubset(tree_taxa):
        return None
    return tree.extract_tree_with_taxa_labels(triplet_taxa)


def _build_species_triplet_metadata(
    species_tree: dendropy.Tree,
    triplets: list[tuple[str, str, str]],
) -> tuple[list[tuple[str, str, str]], dict[tuple[str, str, str], str], list[tuple[str, str, str]]]:
    """Normalize triplets to A/B/C where A and B are sister taxa."""
    normalized_triplets: list[tuple[str, str, str]] = []
    species_triplet_trees: dict[tuple[str, str, str], str] = {}
    skipped_triplets: list[tuple[str, str, str]] = []

    seen: set[tuple[str, str, str]] = set()
    for triplet in triplets:
        subtree = extract_triplet_subtree(species_tree, triplet)
        if subtree is None:
            skipped_triplets.append(triplet)
            continue

        try:
            labels = sorted(triplet)
            sister_pair = _find_sister_pair(subtree)
            abc_triplet = _normalize_abc_from_sister_pair(labels, sister_pair)
        except ValueError:
            skipped_triplets.append(triplet)
            continue

        if abc_triplet in seen:
            continue

        seen.add(abc_triplet)
        normalized_triplets.append(abc_triplet)
        species_triplet_trees[abc_triplet] = subtree.as_string(schema="newick").strip()

    return normalized_triplets, species_triplet_trees, skipped_triplets


def _read_species_tree(species_tree_path: Path) -> dendropy.Tree:
    """Read one rooted species tree from Newick file."""
    return dendropy.Tree.get(path=str(species_tree_path), schema="newick", preserve_underscores=True)


def _parse_outgroup_arg(outgroup_arg: str) -> list[str]:
    """Parse comma-separated outgroup taxa into a cleaned list."""
    parts = [part.strip() for part in outgroup_arg.split(",")]
    return [part for part in parts if part]


def _biophylo_to_dendropy(tree: Tree) -> dendropy.Tree:
    """Convert a rooted Bio.Phylo tree to DendroPy tree."""
    handle = io.StringIO()
    phylo_write([tree], handle, "newick")
    return dendropy.Tree.get(data=handle.getvalue().strip(), schema="newick", preserve_underscores=True)


def _read_and_root_species_tree(species_tree_path: Path, outgroup_taxa: list[str]) -> tuple[dendropy.Tree, dict[str, int]]:
    """Read the species tree, root and prune it on the outgroups as the orchestrator does, and convert to DendroPy."""
    with species_tree_path.open("r", encoding="utf-8") as handle:
        species_trees = list(phylo_parse(handle, "newick"))
    if not species_trees:
        raise ValueError(f"No trees found in species tree file: {species_tree_path}")

    # The rooting error names the taxa in the way when the outgroups do not
    # branch off the tree at a single point.
    rooted_tree, excluded_taxa, missing_taxa, ingroup_taxa = _root_tree_on_outgroup(species_trees[0], outgroup_taxa)

    metadata = {
        "species_outgroup_missing_count": len(missing_taxa),
        "species_outgroup_excluded_count": len(excluded_taxa),
        "species_ingroup_taxa_count": len(ingroup_taxa),
    }
    return _biophylo_to_dendropy(rooted_tree), metadata


def _read_and_root_gene_trees(gene_trees_path: Path, outgroup_taxa: list[str]) -> tuple[list[dendropy.Tree], int]:
    """Read gene trees and root each on first present outgroup; discard trees with no outgroup."""
    trees: list[dendropy.Tree] = []
    discarded_missing_outgroup = 0
    with gene_trees_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            newick = line.strip()
            if not newick:
                continue
            bio_tree = phylo_read(io.StringIO(newick), "newick")
            rooted_tree, used_outgroup, _ = _root_tree_on_any_outgroup(bio_tree, outgroup_taxa)
            if used_outgroup is None:
                discarded_missing_outgroup += 1
                continue
            trees.append(_biophylo_to_dendropy(rooted_tree))
    return trees, discarded_missing_outgroup


def _leaf_label_set(tree: dendropy.Tree) -> set[str]:
    """Return set of labeled leaf taxa for a tree."""
    return {leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label}


def compute_height_arrays(
    species_tree_path: Path,
    gene_trees_path: Path,
    outgroup_taxa: list[str],
    branch_strategy: str,
    max_triplets: int | None = None,
) -> dict[str, list[float] | dict[str, int | str]]:
    """Compute global concordant/discordant height arrays from input files."""
    species_tree, species_rooting_meta = _read_and_root_species_tree(species_tree_path, outgroup_taxa)
    gene_trees, gene_trees_discarded = _read_and_root_gene_trees(gene_trees_path, outgroup_taxa)

    species_taxa = sorted(_leaf_label_set(species_tree))
    if len(species_taxa) < 3:
        raise ValueError("Species tree must contain at least 3 taxa.")

    raw_triplets = list(combinations(species_taxa, 3))
    normalized_triplets, _, skipped_triplets = _build_species_triplet_metadata(species_tree, raw_triplets)
    if max_triplets is not None:
        normalized_triplets = normalized_triplets[: max(0, int(max_triplets))]

    concordant_heights: list[float] = []
    topology_heights = {
        TOPOLOGY_BC: [],
        TOPOLOGY_AC: [],
    }

    for triplet in normalized_triplets:
        triplet_set = set(triplet)
        by_topology = {
            TOPOLOGY_AB: [],
            TOPOLOGY_BC: [],
            TOPOLOGY_AC: [],
        }

        for gene_tree in gene_trees:
            if not triplet_set.issubset(_leaf_label_set(gene_tree)):
                continue

            subtree = extract_triplet_subtree(gene_tree, triplet)
            if subtree is None:
                continue

            try:
                topology = classify_triplet_topology_string(subtree, triplet)
                height = compute_branch_metric(subtree, branch_strategy=branch_strategy)
            except ValueError:
                continue

            by_topology[topology].append(float(height))

        concordant_heights.extend(by_topology[TOPOLOGY_AB])
        topology_heights[TOPOLOGY_BC].extend(by_topology[TOPOLOGY_BC])
        topology_heights[TOPOLOGY_AC].extend(by_topology[TOPOLOGY_AC])

    # Assign discordant1/2 as globally more/less frequent discordant topology.
    bc_heights = topology_heights[TOPOLOGY_BC]
    ac_heights = topology_heights[TOPOLOGY_AC]
    if len(ac_heights) > len(bc_heights):
        disc1_topology = "AC"
        disc2_topology = "BC"
        discordant1_heights = ac_heights
        discordant2_heights = bc_heights
    else:
        disc1_topology = "BC"
        disc2_topology = "AC"
        discordant1_heights = bc_heights
        discordant2_heights = ac_heights

    return {
        "concordant": concordant_heights,
        "discordant1": discordant1_heights,
        "discordant2": discordant2_heights,
        "metadata": {
            "n_concordant": len(concordant_heights),
            "n_discordant1": len(discordant1_heights),
            "n_discordant2": len(discordant2_heights),
            "disc1_topology": disc1_topology,
            "disc2_topology": disc2_topology,
            "species_taxa_count": len(species_taxa),
            "gene_tree_count": len(gene_trees),
            "gene_trees_discarded_missing_outgroup": gene_trees_discarded,
            "triplets_total": len(raw_triplets),
            "triplets_used": len(normalized_triplets),
            "triplets_skipped": len(skipped_triplets),
            "branch_strategy": branch_strategy,
            **species_rooting_meta,
        },
    }


def _gaussian_kde_numpy(values: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    """Estimate Gaussian KDE using NumPy only (Silverman bandwidth)."""
    n = values.size
    if n < 2:
        return np.zeros_like(x_grid)

    std = float(np.std(values, ddof=1))
    if std <= 0.0:
        return np.zeros_like(x_grid)

    bandwidth = 1.06 * std * (n ** (-1.0 / 5.0))
    if bandwidth <= 0.0:
        return np.zeros_like(x_grid)

    z = (x_grid[:, None] - values[None, :]) / bandwidth
    kernel = np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    return np.mean(kernel, axis=1) / bandwidth


def _plot_height_arrays(
    height_arrays: dict[str, list[float]],
    output_plot: Path,
    alpha: float,
    bins: int,
    branch_strategy: str,
):
    """Plot overlaid histograms with NumPy-based KDE curves for each topology class."""
    import matplotlib.pyplot as plt

    output_plot.parent.mkdir(parents=True, exist_ok=True)

    datasets = {
        "concordant": np.asarray(height_arrays["concordant"], dtype=float),
        "discordant1": np.asarray(height_arrays["discordant1"], dtype=float),
        "discordant2": np.asarray(height_arrays["discordant2"], dtype=float),
    }
    colors = {
        "concordant": "tab:blue",
        "discordant1": "tab:orange",
        "discordant2": "tab:green",
    }

    non_empty = [values for values in datasets.values() if values.size > 0]
    if not non_empty:
        raise ValueError("No tree-height values available to plot.")

    global_min = min(float(values.min()) for values in non_empty)
    global_max = max(float(values.max()) for values in non_empty)
    if np.isclose(global_min, global_max):
        padding = max(abs(global_min) * 0.1, 1e-6)
        global_min -= padding
        global_max += padding
    x_grid = np.linspace(global_min, global_max, 400)
    bin_width = (global_max - global_min) / bins if bins > 0 else 1.0

    plt.figure(figsize=(11, 6))
    for label, values in datasets.items():
        if values.size > 0:
            plt.hist(
                values,
                bins=bins,
                range=(global_min, global_max),
                density=False,
                alpha=max(0.05, alpha * 0.6),
                color=colors[label],
                label=f"{label} (hist)",
            )

        if values.size >= 2 and np.unique(values).size >= 2:
            density = _gaussian_kde_numpy(values, x_grid)
            # Scale KDE to histogram counts so both layers share the same y-axis.
            scaled_density = density * values.size * bin_width
            plt.plot(x_grid, scaled_density, color=colors[label], linewidth=2, label=f"{label} (kde)")
            continue

        if values.size == 1:
            # Fallback for one-point samples where KDE is undefined.
            x_val = float(values[0])
            plt.axvline(x=x_val, color=colors[label], alpha=0.9, linewidth=2, label=f"{label} (single point)")

    axis_label_map = {
        "avg_height": "Average Root-to-Tip Height",
        "internal_branch": "Internal Branch Length",
        "sister_distance": "Sister Pair Distance",
        "external_branch": "External Branch to Root",
    }

    plt.xlabel(axis_label_map.get(branch_strategy, "Branch Metric Value"))
    plt.ylabel("Count")
    plt.title(f"Triplet Topology Distributions ({branch_strategy})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=150)
    plt.close()


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute concordant/discordant1/discordant2 tree-height arrays from "
            "species and gene tree files, then plot them together."
        )
    )
    parser.add_argument("--species-tree-path", required=True, help="Path to species tree Newick file.")
    parser.add_argument("--gene-trees-path", required=True, help="Path to gene trees Newick file.")
    parser.add_argument(
        "--outgroup",
        required=True,
        help="Outgroup taxon name(s), comma-separated (same behavior as tree_parser).",
    )
    parser.add_argument(
        "--branch-strategy",
        default="avg_height",
        choices=BRANCH_STRATEGY_CHOICES,
        help=(
            "Metric used for concordant/discordant arrays: "
            "avg_height (default), internal_branch, sister_distance, external_branch."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Output directory for plot and text log (default: results).",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=100,
        help="Number of histogram bins (default: 100).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.45,
        help="Transparency for filled distribution curves in [0, 1] (default: 0.45).",
    )
    parser.add_argument(
        "--max-triplets",
        type=int,
        default=10,
        help="Optional cap on number of species triplets to process.",
    )
    return parser


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    species_tree_path = Path(args.species_tree_path).expanduser().resolve()
    gene_trees_path = Path(args.gene_trees_path).expanduser().resolve()
    outgroup_taxa = _parse_outgroup_arg(args.outgroup)
    branch_strategy = args.branch_strategy
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_plot = output_dir / "triplet_tree_heights_distribution.png"
    output_arrays_json = output_dir / "triplet_tree_heights_arrays.json"
    output_log = output_dir / "triplet_tree_heights_output.txt"

    if not species_tree_path.exists():
        raise FileNotFoundError(f"Species tree file not found: {species_tree_path}")
    if not gene_trees_path.exists():
        raise FileNotFoundError(f"Gene trees file not found: {gene_trees_path}")
    if not outgroup_taxa:
        raise ValueError("--outgroup must provide at least one taxon.")
    if args.bins <= 0:
        raise ValueError("--bins must be a positive integer.")
    if not (0.0 <= args.alpha <= 1.0):
        raise ValueError("--alpha must be between 0 and 1.")

    output_dir.mkdir(parents=True, exist_ok=True)

    results = compute_height_arrays(
        species_tree_path=species_tree_path,
        gene_trees_path=gene_trees_path,
        outgroup_taxa=outgroup_taxa,
        branch_strategy=branch_strategy,
        max_triplets=args.max_triplets,
    )
    concordant = results["concordant"]
    discordant1 = results["discordant1"]
    discordant2 = results["discordant2"]
    metadata = results["metadata"]

    if not isinstance(concordant, list) or not isinstance(discordant1, list) or not isinstance(discordant2, list):
        raise TypeError("Unexpected topology height array format in compute_height_arrays output.")
    if not isinstance(metadata, dict):
        raise TypeError("Unexpected metadata format in compute_height_arrays output.")

    heights: dict[str, list[float]] = {
        "concordant": concordant,
        "discordant1": discordant1,
        "discordant2": discordant2,
    }

    print("Computed tree heights and metadata.")
    print(
        "Topology counts: "
        f"concordant={metadata['n_concordant']}, "
        f"discordant1={metadata['n_discordant1']}, "
        f"discordant2={metadata['n_discordant2']}"
    )
    print(
        "Discordant topology mapping: "
        f"disc1_topology={metadata['disc1_topology']}, "
        f"disc2_topology={metadata['disc2_topology']}"
    )
    print(json.dumps(metadata, indent=2))

    with output_arrays_json.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(heights, indent=2))
        handle.write("\n")

    with output_log.open("w", encoding="utf-8") as handle:
        handle.write("Tree height metadata (JSON)\n")
        handle.write(json.dumps(metadata, indent=2))
        handle.write("\n")

    _plot_height_arrays(
        heights,
        output_plot=output_plot,
        alpha=args.alpha,
        bins=args.bins,
        branch_strategy=branch_strategy,
    )

    print(f"Saved plot: {output_plot}")
    print(f"Saved arrays JSON: {output_arrays_json}")
    print(f"Saved text output: {output_log}")


if __name__ == "__main__":
    main()