"""Preflight sanity checks for species/gene trees before GhostParser runs.

This script validates common failure modes that later surface as runtime ValueError
messages in triplet extraction and inference, including rooted-sister detection
failures and triplet shape mismatches.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations

import dendropy

from ghostparser.tree_parser import (
    _parse_outgroup_arg,
    _root_tree_on_any_outgroup,
    format_newick_with_precision,
    read_tree_file,
    read_triplet_filter_file,
    standardize_tree,
)
from ghostparser.triplet_utils import (
    find_sister_pair,
    normalize_abc_from_sister_pair,
    topology_from_sister_pair,
    triplet_taxa_labels,
)


@dataclass(frozen=True)
class Issue:
    """One detected data-quality issue."""

    category: str
    message: str


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run preflight sanity checks for species and gene trees using GhostParser "
            "triplet logic."
        )
    )
    parser.add_argument(
        "--species-tree-path",
        required=True,
        help="Path to species tree file (Newick).",
    )
    parser.add_argument(
        "--gene-trees-path",
        required=True,
        help="Path to gene trees file (Newick; one tree per line).",
    )
    parser.add_argument(
        "--outgroup",
        required=True,
        help="Comma-separated outgroup taxa, e.g. O or O1,O2.",
    )
    parser.add_argument(
        "--triplet-filter",
        default=None,
        help="Optional file with one comma-separated triplet per line.",
    )
    parser.add_argument(
        "--max-triplets",
        type=int,
        default=1500,
        help=(
            "Maximum number of triplets to analyze when --triplet-filter is not provided "
            "(default: 1500)."
        ),
    )
    parser.add_argument(
        "--max-gene-trees",
        type=int,
        default=0,
        help="Analyze at most N gene trees (0 means all).",
    )
    parser.add_argument(
        "--report-limit",
        type=int,
        default=25,
        help="Maximum example lines shown per issue category (default: 25).",
    )
    return parser.parse_args()


def _load_single_species_tree(path: str):
    trees = read_tree_file(path)
    if len(trees) != 1:
        raise ValueError(
            f"Species tree file must contain exactly one tree; found {len(trees)}"
        )
    return trees[0]


def _to_dendropy_tree(biopython_tree) -> dendropy.Tree:
    newick = format_newick_with_precision(biopython_tree)
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def _get_leaf_labels_dendropy(tree: dendropy.Tree) -> set[str]:
    return {
        leaf.taxon.label
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None and leaf.taxon.label
    }


def _load_target_triplets(
    species_labels_sorted: list[str],
    outgroups: set[str],
    triplet_filter: str | None,
    max_triplets: int,
    issues: list[Issue],
) -> list[tuple[str, str, str]]:
    ingroup = sorted(
        [taxon for taxon in species_labels_sorted if taxon not in outgroups]
    )
    if len(ingroup) < 3:
        issues.append(
            Issue(
                category="species_tree.insufficient_ingroup_taxa",
                message=(
                    f"Need at least 3 ingroup taxa after excluding outgroups; found {len(ingroup)}"
                ),
            )
        )
        return []

    if triplet_filter:
        raw_triplets, invalid_lines = read_triplet_filter_file(triplet_filter)
        for line_num, raw in invalid_lines:
            issues.append(
                Issue(
                    category="triplet_filter.invalid_line",
                    message=(
                        f"Invalid triplet-filter line {line_num}: '{raw}' "
                        "(expected exactly 3 comma-separated taxa)"
                    ),
                )
            )

        unique = []
        seen = set()
        species_label_set = set(species_labels_sorted)
        for triplet in raw_triplets:
            if triplet in seen:
                continue
            seen.add(triplet)

            missing = sorted(set(triplet) - species_label_set)
            if missing:
                issues.append(
                    Issue(
                        category="triplet_filter.taxa_missing_in_species_tree",
                        message=(
                            f"Triplet {','.join(triplet)} contains taxa absent from species tree: "
                            f"{', '.join(missing)}"
                        ),
                    )
                )
                continue

            if any(taxon in outgroups for taxon in triplet):
                issues.append(
                    Issue(
                        category="triplet_filter.includes_outgroup",
                        message=(
                            f"Triplet {','.join(triplet)} includes outgroup taxon; "
                            "triplets should use ingroup taxa only"
                        ),
                    )
                )
                continue

            unique.append(tuple(triplet))
        return unique

    all_triplets = list(combinations(ingroup, 3))
    if max_triplets > 0 and len(all_triplets) > max_triplets:
        issues.append(
            Issue(
                category="analysis.triplet_cap_applied",
                message=(
                    f"Triplet analysis capped at {max_triplets} of {len(all_triplets)} possible ingroup triplets"
                ),
            )
        )
        return all_triplets[:max_triplets]

    return all_triplets


def _normalize_species_triplets(
    species_tree_d: dendropy.Tree,
    triplets: list[tuple[str, str, str]],
    issues: list[Issue],
) -> list[tuple[str, str, str]]:
    normalized = []
    seen = set()

    for triplet in triplets:
        subtree = species_tree_d.extract_tree_with_taxa_labels(triplet)
        if subtree is None:
            issues.append(
                Issue(
                    category="species_triplet.subtree_missing",
                    message=f"Species triplet subtree missing for {','.join(triplet)}",
                )
            )
            continue

        try:
            labels = sorted(triplet)
            sister_pair = find_sister_pair(subtree)
            abc_triplet = normalize_abc_from_sister_pair(labels, sister_pair)
            topology_from_sister_pair(sister_pair, abc_triplet)
        except ValueError as exc:
            issues.append(
                Issue(
                    category="species_triplet.invalid_rooting_or_topology",
                    message=(
                        f"Species triplet {','.join(triplet)} failed rooted topology checks: {exc}"
                    ),
                )
            )
            continue

        if abc_triplet in seen:
            continue
        seen.add(abc_triplet)
        normalized.append(abc_triplet)

    return normalized


def _load_gene_tree_lines(gene_tree_path: str, max_gene_trees: int) -> list[str]:
    """Load non-empty raw gene-tree lines to report exact problematic input."""
    lines: list[str] = []
    with open(gene_tree_path, "r") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            lines.append(stripped)
            if max_gene_trees > 0 and len(lines) >= max_gene_trees:
                break
    return lines


def _check_gene_tree_triplets(
    gene_tree_path: str,
    outgroups: list[str],
    normalized_triplets: list[tuple[str, str, str]],
    max_gene_trees: int,
    issues: list[Issue],
) -> dict[str, int]:
    gene_trees = read_tree_file(gene_tree_path)
    if max_gene_trees > 0:
        gene_trees = gene_trees[:max_gene_trees]
    gene_tree_lines = _load_gene_tree_lines(gene_tree_path, max_gene_trees)

    counters: dict[str, int] = defaultdict(int)

    if len(gene_tree_lines) != len(gene_trees):
        issues.append(
            Issue(
                category="gene_tree.input_alignment_warning",
                message=(
                    "Parsed gene-tree count does not match non-empty input line count; "
                    f"parsed={len(gene_trees)}, lines={len(gene_tree_lines)}"
                ),
            )
        )

    for idx, tree in enumerate(gene_trees, start=1):
        line_preview = (
            gene_tree_lines[idx - 1]
            if idx - 1 < len(gene_tree_lines)
            else "<unavailable>"
        )
        rooted, used_outgroup, missing_outgroups = _root_tree_on_any_outgroup(
            tree, outgroups
        )
        if used_outgroup is None:
            counters["gene_tree.rooting_failed"] += 1
            issues.append(
                Issue(
                    category="gene_tree.rooting_failed",
                    message=(
                        f"Gene tree #{idx}: none of the outgroups were present; "
                        f"missing outgroups: {', '.join(sorted(missing_outgroups)) or 'all'}; "
                        f"gene_tree_line={line_preview}"
                    ),
                )
            )
            continue

        counters["gene_tree.rooted"] += 1

        try:
            rooted_std = standardize_tree(rooted)
            tree_d = _to_dendropy_tree(rooted_std)
        except Exception as exc:  # defensive parse guard
            counters["gene_tree.parse_after_rooting_failed"] += 1
            issues.append(
                Issue(
                    category="gene_tree.parse_after_rooting_failed",
                    message=(
                        f"Gene tree #{idx}: failed to parse after rooting: {exc}; "
                        f"gene_tree_line={line_preview}"
                    ),
                )
            )
            continue

        gene_labels = _get_leaf_labels_dendropy(tree_d)

        for triplet in normalized_triplets:
            if not set(triplet).issubset(gene_labels):
                continue

            subtree = tree_d.extract_tree_with_taxa_labels(triplet)
            if subtree is None:
                counters["triplet.subtree_extraction_failed"] += 1
                issues.append(
                    Issue(
                        category="triplet.subtree_extraction_failed",
                        message=(
                            f"Gene tree #{idx}, triplet {','.join(triplet)}: "
                            "subtree extraction returned None; "
                            f"gene_tree_line={line_preview}"
                        ),
                    )
                )
                continue

            try:
                labels = triplet_taxa_labels(subtree)
            except ValueError as exc:
                counters["triplet.invalid_leaf_count"] += 1
                issues.append(
                    Issue(
                        category="triplet.invalid_leaf_count",
                        message=(
                            f"Gene tree #{idx}, triplet {','.join(triplet)}: {exc}; "
                            f"gene_tree_line={line_preview}"
                        ),
                    )
                )
                continue

            try:
                sister_pair = find_sister_pair(subtree)
            except ValueError as exc:
                counters["triplet.unresolved_rooted_sister_pair"] += 1
                issues.append(
                    Issue(
                        category="triplet.unresolved_rooted_sister_pair",
                        message=(
                            f"Gene tree #{idx}, triplet {','.join(triplet)}: {exc} "
                            "(likely unresolved/polytomous triplet or ambiguous rooting); "
                            f"gene_tree_line={line_preview}"
                        ),
                    )
                )
                continue

            try:
                topology_from_sister_pair(sister_pair, triplet)
            except ValueError as exc:
                counters["triplet.abc_mapping_mismatch"] += 1
                issues.append(
                    Issue(
                        category="triplet.abc_mapping_mismatch",
                        message=(
                            f"Gene tree #{idx}, triplet {','.join(triplet)} with labels "
                            f"{','.join(labels)}: {exc}; gene_tree_line={line_preview}"
                        ),
                    )
                )
                continue

            mrca = subtree.mrca(taxon_labels=sorted(sister_pair))
            if mrca is None:
                counters["triplet.sister_mrca_missing"] += 1
                issues.append(
                    Issue(
                        category="triplet.sister_mrca_missing",
                        message=(
                            f"Gene tree #{idx}, triplet {','.join(triplet)}: "
                            "sister-pair MRCA not found; "
                            f"gene_tree_line={line_preview}"
                        ),
                    )
                )

    counters["gene_tree.total_checked"] = len(gene_trees)
    return counters


def _print_report(
    species_tree_path: str,
    gene_trees_path: str,
    outgroups: list[str],
    normalized_triplets: list[tuple[str, str, str]],
    counters: dict[str, int],
    issues: list[Issue],
    report_limit: int,
) -> int:
    print("=== GhostParser Preflight Sanity Check ===")
    print(f"species_tree_path: {species_tree_path}")
    print(f"gene_trees_path:   {gene_trees_path}")
    print(f"outgroups:         {','.join(outgroups)}")
    print(f"triplets_checked:  {len(normalized_triplets)}")
    print(f"gene_trees_checked:{counters.get('gene_tree.total_checked', 0)}")
    print()

    category_counts: dict[str, int] = defaultdict(int)
    category_examples: dict[str, list[str]] = defaultdict(list)
    for issue in issues:
        category_counts[issue.category] += 1
        if len(category_examples[issue.category]) < report_limit:
            category_examples[issue.category].append(issue.message)

    species_issue_count = sum(
        count
        for category, count in category_counts.items()
        if any(
            category.startswith(prefix)
            for prefix in ("species_tree.", "species_triplet.", "triplet_filter.")
        )
    )
    gene_issue_count = sum(
        count
        for category, count in category_counts.items()
        if category.startswith("gene_tree.") or category.startswith("triplet.")
    )

    if not category_counts:
        print("No blocking data issues detected for the checked trees/triplets.")
        print(
            "This does not guarantee biological signal quality; it checks structural sanity."
        )
        return 0

    print("Potential source assessment:")
    print(
        f"- species-tree-related issues: {'YES' if species_issue_count > 0 else 'NO'} (count={species_issue_count})"
    )
    print(
        f"- gene-tree-related issues:    {'YES' if gene_issue_count > 0 else 'NO'} (count={gene_issue_count})"
    )
    print()

    print("Detected issues:")
    for category in sorted(category_counts.keys()):
        print(f"- {category}: {category_counts[category]}")
    print()

    print("Examples:")
    for category in sorted(category_examples.keys()):
        print(f"[{category}]")
        for line in category_examples[category]:
            print(f"  - {line}")
        print()

    print("Guidance:")
    print(
        "- gene_tree.rooting_failed: outgroup labels are missing/mismatched in gene trees."
    )
    print(
        "- triplet.unresolved_rooted_sister_pair: often unresolved triplets/polytomies or ambiguous rooting; "
        "these can trigger 'Could not determine rooted sister pair for triplet tree'."
    )
    print(
        "- species_triplet.invalid_rooting_or_topology: species tree rooting/topology itself failed triplet checks."
    )

    return 1


def main() -> None:
    args = _parse_args()

    issues: list[Issue] = []
    outgroups = _parse_outgroup_arg(args.outgroup)
    if not outgroups:
        raise ValueError("No outgroup taxa were provided")

    species_tree_bio = _load_single_species_tree(args.species_tree_path)
    rooted_species_tree, used_outgroup, missing_outgroups = _root_tree_on_any_outgroup(
        species_tree_bio, outgroups
    )
    if used_outgroup is None:
        missing_str = ", ".join(sorted(missing_outgroups)) or "all"
        raise ValueError(
            "Could not root species tree: none of the provided outgroups were found "
            f"(missing outgroups: {missing_str})"
        )

    species_tree_std = standardize_tree(rooted_species_tree)
    species_tree_d = _to_dendropy_tree(species_tree_std)
    species_labels_sorted = sorted(_get_leaf_labels_dendropy(species_tree_d))

    target_triplets = _load_target_triplets(
        species_labels_sorted=species_labels_sorted,
        outgroups=set(outgroups),
        triplet_filter=args.triplet_filter,
        max_triplets=args.max_triplets,
        issues=issues,
    )

    normalized_triplets = _normalize_species_triplets(
        species_tree_d=species_tree_d,
        triplets=target_triplets,
        issues=issues,
    )

    counters = _check_gene_tree_triplets(
        gene_tree_path=args.gene_trees_path,
        outgroups=outgroups,
        normalized_triplets=normalized_triplets,
        max_gene_trees=args.max_gene_trees,
        issues=issues,
    )

    exit_code = _print_report(
        species_tree_path=args.species_tree_path,
        gene_trees_path=args.gene_trees_path,
        outgroups=outgroups,
        normalized_triplets=normalized_triplets,
        counters=counters,
        issues=issues,
        report_limit=args.report_limit,
    )
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
