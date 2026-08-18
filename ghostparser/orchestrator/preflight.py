"""Structural preflight checks for species and gene trees.

Reproduces the engine's structural logic and reports every problem at once,
rather than raising on the first from deep inside a run. The checks say whether
the data can be processed, not whether the result is biologically meaningful.
"""

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import dendropy

from ..triplet_utils import (
    find_sister_pair,
    normalize_abc_from_sister_pair,
    topology_from_sister_pair,
    triplet_taxa_labels,
)
from .trees import (
    _root_tree_on_any_outgroup,
    format_newick_with_precision,
    read_tree_file,
    read_triplet_filter_file,
    standardize_tree,
)

PREFLIGHT_REPORT_FILENAME = "preflight_data_check.txt"

# Caps keep the check fast on large inputs; both are analysis-only limits and do
# not change what a real run would process.
DEFAULT_MAX_TRIPLETS = 1500
DEFAULT_MAX_GENE_TREES = 0  # 0 means "all"
DEFAULT_REPORT_LIMIT = 25

_SPECIES_CATEGORY_PREFIXES = ("species_tree.", "species_triplet.", "triplet_filter.")
_GENE_CATEGORY_PREFIXES = ("gene_tree.", "triplet.")


@dataclass(frozen=True)
class Issue:
    """One detected data-quality issue.

    Attributes:
        category: Dotted category key used to group and count issues.
        message: Human-readable description including the offending gene-tree
            index or triplet where applicable.
    """

    category: str
    message: str


@dataclass(frozen=True)
class PreflightResult:
    """Outcome of a preflight run.

    Attributes:
        report_path: Path of the written report, or ``None`` when the caller
            requested no file.
        report_text: The full report as a single string.
        issues: Every issue detected, in discovery order.
        counters: Per-category tallies plus ``gene_tree.total_checked``.
        triplets_checked: Number of normalized species triplets analyzed.
        passed: ``True`` when no issues were detected.
    """

    report_path: str | None
    report_text: str
    issues: list[Issue]
    counters: dict[str, int]
    triplets_checked: int
    passed: bool


def _load_single_species_tree(path):
    """Read a species-tree file that must contain exactly one tree."""
    trees = read_tree_file(path)
    if len(trees) != 1:
        raise ValueError(
            f"Species tree file must contain exactly one tree; found {len(trees)}"
        )
    return trees[0]


def _to_dendropy_tree(biopython_tree) -> dendropy.Tree:
    """Convert a BioPython tree to DendroPy via a full-precision Newick string."""
    newick = format_newick_with_precision(biopython_tree)
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def _get_leaf_labels_dendropy(tree: dendropy.Tree) -> set[str]:
    """Collect the non-empty leaf labels of a DendroPy tree."""
    return {
        leaf.taxon.label
        for leaf in tree.leaf_node_iter()
        if leaf.taxon is not None and leaf.taxon.label
    }


def _load_target_triplets(
    species_labels_sorted,
    outgroups,
    triplet_filter,
    max_triplets,
    issues,
):
    """Resolve the triplets to check, from a filter file or all ingroup combinations.

    Args:
        species_labels_sorted: Sorted species-tree leaf labels.
        outgroups: Set of outgroup taxon labels to exclude from the ingroup.
        triplet_filter: Optional path to a triplet-filter file.
        max_triplets: Cap on generated triplets; ``0`` means no cap.
        issues: Mutable list that detected issues are appended to.

    Returns:
        A list of taxon triples to check. Empty when the ingroup is too small.
    """
    ingroup = sorted(
        [taxon for taxon in species_labels_sorted if taxon not in outgroups]
    )
    if len(ingroup) < 3:
        issues.append(
            Issue(
                category="species_tree.insufficient_ingroup_taxa",
                message=(
                    "Need at least 3 ingroup taxa after excluding outgroups; "
                    f"found {len(ingroup)}"
                ),
            )
        )
        return []

    if triplet_filter:
        return _load_filtered_triplets(
            triplet_filter, species_labels_sorted, outgroups, issues
        )

    all_triplets = list(combinations(ingroup, 3))
    if max_triplets > 0 and len(all_triplets) > max_triplets:
        issues.append(
            Issue(
                category="analysis.triplet_cap_applied",
                message=(
                    f"Triplet analysis capped at {max_triplets} of "
                    f"{len(all_triplets)} possible ingroup triplets"
                ),
            )
        )
        return all_triplets[:max_triplets]

    return all_triplets


def _load_filtered_triplets(triplet_filter, species_labels_sorted, outgroups, issues):
    """Read a triplet-filter file and drop entries the species tree cannot support."""
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
                        f"Triplet {','.join(triplet)} contains taxa absent from "
                        f"species tree: {', '.join(missing)}"
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


def _normalize_species_triplets(species_tree_d, triplets, issues):
    """Normalize each triplet to A/B/C order, recording those the species tree rejects.

    Args:
        species_tree_d: Rooted, standardized species tree as a DendroPy tree.
        triplets: Candidate taxon triples.
        issues: Mutable list that detected issues are appended to.

    Returns:
        Deduplicated A/B/C-normalized triplets that passed rooted topology checks.
    """
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
                        f"Species triplet {','.join(triplet)} failed rooted "
                        f"topology checks: {exc}"
                    ),
                )
            )
            continue

        if abc_triplet in seen:
            continue
        seen.add(abc_triplet)
        normalized.append(abc_triplet)

    return normalized


def _load_gene_tree_lines(gene_tree_path, max_gene_trees):
    """Load non-empty raw gene-tree lines so reports can quote the offending input."""
    lines = []
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
    gene_tree_path,
    outgroups,
    normalized_triplets,
    max_gene_trees,
    issues,
):
    """Root each gene tree and replay triplet extraction against it.

    Args:
        gene_tree_path: Path to the gene-trees Newick file.
        outgroups: Ordered outgroup taxon labels.
        normalized_triplets: A/B/C-normalized triplets from the species tree.
        max_gene_trees: Cap on gene trees to check; ``0`` means all.
        issues: Mutable list that detected issues are appended to.

    Returns:
        A mapping of counter name to count, including ``gene_tree.total_checked``.
    """
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
                    "Parsed gene-tree count does not match non-empty input line "
                    f"count; parsed={len(gene_trees)}, lines={len(gene_tree_lines)}"
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
            if set(triplet).issubset(gene_labels):
                _check_one_triplet(tree_d, triplet, idx, line_preview, counters, issues)

    counters["gene_tree.total_checked"] = len(gene_trees)
    return counters


def _check_one_triplet(tree_d, triplet, gene_index, line_preview, counters, issues):
    """Replay subtree extraction and sister-pair resolution for one triplet."""
    subtree = tree_d.extract_tree_with_taxa_labels(triplet)
    if subtree is None:
        counters["triplet.subtree_extraction_failed"] += 1
        issues.append(
            Issue(
                category="triplet.subtree_extraction_failed",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: "
                    f"subtree extraction returned None; gene_tree_line={line_preview}"
                ),
            )
        )
        return

    try:
        labels = triplet_taxa_labels(subtree)
    except ValueError as exc:
        counters["triplet.invalid_leaf_count"] += 1
        issues.append(
            Issue(
                category="triplet.invalid_leaf_count",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: {exc}; "
                    f"gene_tree_line={line_preview}"
                ),
            )
        )
        return

    try:
        sister_pair = find_sister_pair(subtree)
    except ValueError as exc:
        counters["triplet.unresolved_rooted_sister_pair"] += 1
        issues.append(
            Issue(
                category="triplet.unresolved_rooted_sister_pair",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: {exc} "
                    "(likely unresolved/polytomous triplet or ambiguous rooting); "
                    f"gene_tree_line={line_preview}"
                ),
            )
        )
        return

    try:
        topology_from_sister_pair(sister_pair, triplet)
    except ValueError as exc:
        counters["triplet.abc_mapping_mismatch"] += 1
        issues.append(
            Issue(
                category="triplet.abc_mapping_mismatch",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)} with "
                    f"labels {','.join(labels)}: {exc}; gene_tree_line={line_preview}"
                ),
            )
        )
        return

    if subtree.mrca(taxon_labels=sorted(sister_pair)) is None:
        counters["triplet.sister_mrca_missing"] += 1
        issues.append(
            Issue(
                category="triplet.sister_mrca_missing",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: "
                    f"sister-pair MRCA not found; gene_tree_line={line_preview}"
                ),
            )
        )


def _build_report(
    species_tree_path,
    gene_trees_path,
    outgroups,
    triplets_checked,
    counters,
    issues,
    report_limit,
):
    """Render the human-readable preflight report.

    Args:
        species_tree_path: Species-tree input path, echoed in the header.
        gene_trees_path: Gene-trees input path, echoed in the header.
        outgroups: Ordered outgroup taxon labels.
        triplets_checked: Number of normalized triplets analyzed.
        counters: Per-category tallies from the gene-tree pass.
        issues: Every detected issue.
        report_limit: Maximum example lines printed per category.

    Returns:
        The report as a single newline-joined string.
    """
    lines = [
        "=== GhostParser Preflight Data Check ===",
        f"species_tree_path:  {species_tree_path}",
        f"gene_trees_path:    {gene_trees_path}",
        f"outgroups:          {','.join(outgroups)}",
        f"triplets_checked:   {triplets_checked}",
        f"gene_trees_checked: {counters.get('gene_tree.total_checked', 0)}",
        "",
    ]

    category_counts: dict[str, int] = defaultdict(int)
    category_examples: dict[str, list[str]] = defaultdict(list)
    for issue in issues:
        category_counts[issue.category] += 1
        if len(category_examples[issue.category]) < report_limit:
            category_examples[issue.category].append(issue.message)

    if not category_counts:
        lines.append("No blocking data issues detected for the checked trees/triplets.")
        lines.append(
            "This does not guarantee biological signal quality; it checks "
            "structural sanity."
        )
        return "\n".join(lines) + "\n"

    species_issue_count = sum(
        count
        for category, count in category_counts.items()
        if category.startswith(_SPECIES_CATEGORY_PREFIXES)
    )
    gene_issue_count = sum(
        count
        for category, count in category_counts.items()
        if category.startswith(_GENE_CATEGORY_PREFIXES)
    )

    lines.extend(
        [
            "Potential source assessment:",
            f"- species-tree-related issues: {'YES' if species_issue_count else 'NO'} "
            f"(count={species_issue_count})",
            f"- gene-tree-related issues:    {'YES' if gene_issue_count else 'NO'} "
            f"(count={gene_issue_count})",
            "",
            "Detected issues:",
        ]
    )
    for category in sorted(category_counts):
        lines.append(f"- {category}: {category_counts[category]}")
    lines.append("")

    lines.append("Examples:")
    for category in sorted(category_examples):
        lines.append(f"[{category}]")
        lines.extend(f"  - {message}" for message in category_examples[category])
        lines.append("")

    lines.extend(
        [
            "Guidance:",
            "- gene_tree.rooting_failed: outgroup labels are missing/mismatched "
            "in gene trees.",
            "- triplet.unresolved_rooted_sister_pair: often unresolved "
            "triplets/polytomies or ambiguous rooting; these can trigger "
            "'Could not determine rooted sister pair for triplet tree'.",
            "- species_triplet.invalid_rooting_or_topology: species tree "
            "rooting/topology itself failed triplet checks.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_preflight_data_check(
    species_tree_path,
    gene_trees_path,
    outgroups,
    output_dir=None,
    triplet_filter=None,
    max_triplets=DEFAULT_MAX_TRIPLETS,
    max_gene_trees=DEFAULT_MAX_GENE_TREES,
    report_limit=DEFAULT_REPORT_LIMIT,
):
    """Run every structural check and write the report.

    Args:
        species_tree_path: Path to the species-tree Newick file.
        gene_trees_path: Path to the gene-trees Newick file.
        outgroups: Ordered outgroup taxon labels; must be non-empty.
        output_dir: Directory to write ``preflight_data_check.txt`` into. When
            ``None``, no file is written and only the text is returned.
        triplet_filter: Optional triplet-filter file restricting the triplets.
        max_triplets: Cap on generated triplets; ``0`` means no cap.
        max_gene_trees: Cap on gene trees checked; ``0`` means all.
        report_limit: Maximum example lines per issue category.

    Returns:
        A :class:`PreflightResult`.

    Raises:
        ValueError: If no outgroups are given, the species-tree file does not
            hold exactly one tree, or none of the outgroups are present in the
            species tree — all conditions that make the check itself impossible.
    """
    if not outgroups:
        raise ValueError("No outgroup taxa were provided")

    issues: list[Issue] = []

    species_tree_bio = _load_single_species_tree(species_tree_path)
    rooted_species_tree, used_outgroup, missing_outgroups = _root_tree_on_any_outgroup(
        species_tree_bio, outgroups
    )
    if used_outgroup is None:
        missing_str = ", ".join(sorted(missing_outgroups)) or "all"
        raise ValueError(
            "Could not root species tree: none of the provided outgroups were "
            f"found (missing outgroups: {missing_str})"
        )

    species_tree_d = _to_dendropy_tree(standardize_tree(rooted_species_tree))
    species_labels_sorted = sorted(_get_leaf_labels_dendropy(species_tree_d))

    target_triplets = _load_target_triplets(
        species_labels_sorted=species_labels_sorted,
        outgroups=set(outgroups),
        triplet_filter=triplet_filter,
        max_triplets=max_triplets,
        issues=issues,
    )
    normalized_triplets = _normalize_species_triplets(
        species_tree_d=species_tree_d,
        triplets=target_triplets,
        issues=issues,
    )
    counters = _check_gene_tree_triplets(
        gene_tree_path=gene_trees_path,
        outgroups=outgroups,
        normalized_triplets=normalized_triplets,
        max_gene_trees=max_gene_trees,
        issues=issues,
    )

    report_text = _build_report(
        species_tree_path=species_tree_path,
        gene_trees_path=gene_trees_path,
        outgroups=outgroups,
        triplets_checked=len(normalized_triplets),
        counters=counters,
        issues=issues,
        report_limit=report_limit,
    )

    report_path = None
    if output_dir is not None:
        destination = Path(output_dir) / PREFLIGHT_REPORT_FILENAME
        destination.write_text(report_text, encoding="utf-8")
        report_path = str(destination)

    return PreflightResult(
        report_path=report_path,
        report_text=report_text,
        issues=issues,
        counters=dict(counters),
        triplets_checked=len(normalized_triplets),
        passed=not issues,
    )
