"""Structural preflight check that replays the run's tree handling and reports
every problem at once.
"""

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import dendropy

from ..config import InputError
from ..triplet_utils import normalize_abc_from_sister_pair
from .config import DEFAULT_PREFLIGHT_TRIPLET_CAP
from .trees import (
    has_branch_lengths,
    root_gene_tree,
    _root_tree_on_outgroup,
    format_newick_with_precision,
    read_tree_file,
    read_species_filter_file,
    read_triplet_filter_file,
    standardize_tree,
)
from .triplet_geometry import (
    TRIPLET_MISSING_TAXON,
    TRIPLET_UNRESOLVED,
    build_taxon_index,
    build_triplet_geometry,
    triplet_resolution,
    triplet_subtree_shape,
)

PREFLIGHT_REPORT_FILENAME = "preflight_data_check.txt"

# The gene-tree cap is analysis-only like the triplet cap the config owns; it
# does not change what a real run would process.
DEFAULT_MAX_GENE_TREES = 0  # 0 means "all"
DEFAULT_REPORT_LIMIT = 25

_SPECIES_CATEGORY_PREFIXES = (
    "species_tree.",
    "species_triplet.",
    "triplet_filter.",
    "species_filter.",
)
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
        counters: Per-category tallies plus ``gene_tree.total_checked``,
            ``gene_tree.rooted``, ``gene_tree.farthest.<outgroup>`` (trees
            in which it was the farthest), ``gene_tree.rooted_on.<outgroup>``
            (trees rooted using it) and ``gene_tree.tangled.<outgroup>``
            (trees in which it sat among the ingroup taxa and was pruned
            unused) for each outgroup, ``gene_tree.tangled_trees`` (trees with
            at least one tangled outgroup), ``gene_tree.missing_branch_lengths``
            (rooted trees lacking some branch length, read as 0), and the
            three that account for every triplet/gene-tree pair seen:
            ``triplet.resolved``, ``triplet.unresolved_rooted_sister_pair``,
            and ``triplet.taxa_absent_from_gene_tree``.
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
        raise InputError(
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
    species_filter,
    max_triplets,
    issues,
):
    """Resolve the triplets to check, from a filter file or all ingroup combinations.

    Args:
        species_labels_sorted: Sorted species-tree leaf labels.
        outgroups: Set of outgroup taxon labels to exclude from the ingroup.
        triplet_filter: Optional path to a triplet-filter file.
        species_filter: Optional path to a species-filter file.
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
    if species_filter:
        ingroup = _load_filtered_species(
            species_filter, species_labels_sorted, outgroups, issues
        )
        if len(ingroup) < 3:
            return []

    # Named triplets are never capped; generated ones are, so the check stays
    # quick however many species there are to combine.
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


def _load_filtered_species(species_filter, species_labels_sorted, outgroups, issues):
    """Read a species-filter file and keep the sorted ingroup species it names."""
    species_label_set = set(species_labels_sorted)
    kept = []
    for taxon in read_species_filter_file(species_filter):
        if taxon not in species_label_set:
            issues.append(
                Issue(
                    category="species_filter.taxa_missing_in_species_tree",
                    message=f"Species {taxon} is absent from the species tree",
                )
            )
        elif taxon in outgroups:
            issues.append(
                Issue(
                    category="species_filter.includes_outgroup",
                    message=(
                        f"Species {taxon} is an outgroup taxon; the filter "
                        "should name ingroup taxa only"
                    ),
                )
            )
        else:
            kept.append(taxon)

    if len(kept) < 3:
        issues.append(
            Issue(
                category="species_filter.insufficient_taxa",
                message=(
                    "Need at least 3 usable species to form a triplet; "
                    f"the species filter names {len(kept)}"
                ),
            )
        )
    return sorted(kept)


def _normalize_species_triplets(species_tree_d, triplets, issues):
    """Normalize each triplet to A/B/C order from the species tree's cached
    geometry, recording the ones it cannot resolve.

    Args:
        species_tree_d: Rooted, standardized species tree as a DendroPy tree.
        triplets: Candidate taxon triples.
        issues: Mutable list that detected issues are appended to.

    Returns:
        Deduplicated A/B/C-normalized triplets that passed rooted topology checks.
    """
    normalized = []
    seen = set()

    taxon_index = build_taxon_index(triplets)
    geometry = build_triplet_geometry(species_tree_d, taxon_index)

    for triplet in triplets:
        positions = tuple(taxon_index[label] for label in triplet)
        shape = triplet_subtree_shape(geometry, positions)
        if shape is None:
            issues.append(
                Issue(
                    category="species_triplet.invalid_rooting_or_topology",
                    message=(
                        f"Species triplet {','.join(triplet)} has no resolved "
                        "rooted sister pair in the species tree"
                    ),
                )
            )
            continue

        label_of = dict(zip(positions, triplet))
        sister_pair = frozenset(
            label_of[position] for position in shape.sister_positions
        )
        abc_triplet = normalize_abc_from_sister_pair(sorted(triplet), sister_pair)

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
        outgroups: Outgroup taxon labels in species-tree rank: farthest from
            the ingroup root first, or listed order when the species tree
            lacks branch lengths.
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
    # Mirror the engine: cache each gene tree's geometry once and read every
    # triplet out of it, rather than copying a subtree per triplet.
    taxon_index = build_taxon_index(normalized_triplets)
    triplet_positions = {
        triplet: tuple(taxon_index[label] for label in triplet)
        for triplet in normalized_triplets
    }

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
        # Not a defect: the run keeps the tree and reads a missing length as
        # 0. Checked before rooting, which does the same.
        lacks_length = not has_branch_lengths(tree)
        rooting = root_gene_tree(tree, outgroups)
        if rooting.tree is None:
            counters["gene_tree.rooting_failed"] += 1
            issues.append(
                Issue(
                    category="gene_tree.rooting_failed",
                    message=(
                        f"Gene tree #{idx}: "
                        + (
                            "every taxon is an outgroup"
                            if rooting.present
                            else "none of the outgroups were present"
                        )
                        + f"; missing outgroups: {', '.join(sorted(rooting.missing)) or 'none'}; "
                        f"gene_tree_line={line_preview}"
                    ),
                )
            )
            continue

        counters["gene_tree.rooted"] += 1
        if lacks_length:
            counters["gene_tree.missing_branch_lengths"] += 1
        counters[f"gene_tree.farthest.{rooting.farthest}"] += 1
        for outgroup in rooting.used:
            counters[f"gene_tree.rooted_on.{outgroup}"] += 1
        # Not a defect: the run keeps such a tree, rooted from its farthest
        # outgroup, and prunes the tangled ones unused. Counted so the report
        # can say how often each outgroup sits among the ingroup.
        for outgroup in rooting.tangled:
            counters[f"gene_tree.tangled.{outgroup}"] += 1
        if rooting.tangled:
            counters["gene_tree.tangled_trees"] += 1

        try:
            rooted_std = standardize_tree(rooting.tree)
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

        geometry = build_triplet_geometry(tree_d, taxon_index)
        gene_labels = _get_leaf_labels_dendropy(tree_d)
        for triplet in normalized_triplets:
            if not set(triplet).issubset(gene_labels):
                # Not a defect: a gene tree need not carry every taxon, and the
                # engine skips these pairs too. Counted so the report can
                # account for every pair it looked at.
                counters["triplet.taxa_absent_from_gene_tree"] += 1
                continue
            _check_one_triplet(
                geometry,
                triplet_positions[triplet],
                triplet,
                idx,
                line_preview,
                counters,
                issues,
            )

    counters["gene_tree.total_checked"] = len(gene_trees)
    return counters


def _check_one_triplet(
    geometry, positions, triplet, gene_index, line_preview, counters, issues
):
    """Replay the engine's triplet resolution for one triplet against one tree.

    Args:
        geometry: The gene tree's cached ``TripletGeometry``.
        positions: The triplet's taxon positions in the shared taxon index.
        triplet: The ``(A, B, C)`` triplet, used for the issue message.
        gene_index: 1-based gene-tree index, used for the issue message.
        line_preview: The gene tree's input line, used for the issue message.
        counters: Mutable counter mapping, incremented in place.
        issues: Mutable list that detected issues are appended to.
    """
    status = triplet_resolution(geometry, positions)
    if status == TRIPLET_UNRESOLVED:
        counters["triplet.unresolved_rooted_sister_pair"] += 1
        issues.append(
            Issue(
                category="triplet.unresolved_rooted_sister_pair",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: "
                    "all three pairwise MRCAs coincide, so no rooted sister "
                    "pair can be determined (likely an unresolved/polytomous "
                    f"triplet or ambiguous rooting); gene_tree_line={line_preview}"
                ),
            )
        )
        return

    if status == TRIPLET_MISSING_TAXON:
        # The caller already confirmed all three labels are in this tree, so
        # this means the cache and the label set disagree.
        counters["triplet.geometry_unavailable"] += 1
        issues.append(
            Issue(
                category="triplet.geometry_unavailable",
                message=(
                    f"Gene tree #{gene_index}, triplet {','.join(triplet)}: "
                    "a triplet taxon is present in the tree's labels but has no "
                    f"leaf in its cached geometry; gene_tree_line={line_preview}"
                ),
            )
        )
        return

    counters["triplet.resolved"] += 1


def _build_report(
    species_tree_path,
    gene_trees_path,
    species_rooting,
    triplets_checked,
    counters,
    issues,
    report_limit,
):
    """Render the human-readable preflight report.

    Args:
        species_tree_path: Species-tree input path, echoed in the header.
        gene_trees_path: Gene-trees input path, echoed in the header.
        species_rooting: The species tree's
            :class:`~ghostparser.orchestrator.trees.SpeciesTreeRooting`,
            whose outgroup order the gene trees were rooted in.
        triplets_checked: Number of normalized triplets analyzed.
        counters: Per-category tallies from the gene-tree pass.
        issues: Every detected issue.
        report_limit: Maximum example lines printed per category.

    Returns:
        The report as a single newline-joined string.
    """
    outgroups = species_rooting.outgroup_order
    lines = [
        "=== GhostParser Preflight Data Check ===",
        f"species_tree_path:  {species_tree_path}",
        f"gene_trees_path:    {gene_trees_path}",
        "outgroups:          "
        + ", ".join(
            [
                outgroup
                if species_rooting.distances is None
                else f"{outgroup} {species_rooting.distances[outgroup]:g}"
                for outgroup in species_rooting.ranked
            ]
            + [
                f"{outgroup} (not in the species tree)"
                for outgroup in species_rooting.missing
            ]
        )
        + (
            "  (listed order, as the species tree lacks branch lengths"
            if species_rooting.distances is None
            else "  (farthest from the ingroup root first, by path length"
        )
        + "; it breaks exact ties between a gene tree's outgroups)",
        f"triplets_checked:   {triplets_checked}",
        f"gene_trees_checked: {counters.get('gene_tree.total_checked', 0)}",
        # Each gene tree is rooted at the common ancestor of its farthest
        # outgroup and the others outside the ingroup; an outgroup tangled
        # among the ingroup taxa is pruned without being used. The counts
        # say how often each happened.
        "farthest:           "
        + ", ".join(
            f"{outgroup}: {counters.get(f'gene_tree.farthest.{outgroup}', 0)}"
            for outgroup in outgroups
        )
        + "  (gene trees in which each outgroup was the farthest from the "
        "ingroup taxa)",
        "rooted_on:          "
        + ", ".join(
            f"{outgroup}: {counters.get(f'gene_tree.rooted_on.{outgroup}', 0)}"
            for outgroup in outgroups
        )
        + "  (gene trees rooted using each outgroup)",
        "tangled:            "
        + ", ".join(
            f"{outgroup}: {counters.get(f'gene_tree.tangled.{outgroup}', 0)}"
            for outgroup in outgroups
        )
        + "  (gene trees in which the outgroup sat among the ingroup taxa and was "
        "pruned without being used for rooting)",
        f"tangled_trees:      {counters.get('gene_tree.tangled_trees', 0)}"
        "  (gene trees with at least one tangled outgroup)",
        "missing_lengths:    "
        f"{counters.get('gene_tree.missing_branch_lengths', 0)}"
        "  (gene trees lacking some branch length; each missing length is read "
        "as 0, which may affect their tree heights and the inferences)",
        # How each triplet/gene-tree pair this check looked at would fare in a
        # run: measurable, dropped as unresolved, or skipped because the gene
        # tree does not carry all three taxa. The three sum to the pairs seen.
        f"usable_pairs:       {counters.get('triplet.resolved', 0)} usable, "
        f"{counters.get('triplet.unresolved_rooted_sister_pair', 0)} unresolved, "
        f"{counters.get('triplet.taxa_absent_from_gene_tree', 0)} with a taxon absent",
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
            "in gene trees, or a tree holds nothing but outgroups.",
            "- tangled: an outgroup that sits among the ingroup taxa once the "
            "gene tree is rooted on its farthest outgroup is pruned without "
            "being used; the run keeps the tree, rooted at the common ancestor "
            "of the farthest and the others outside the ingroup. The farthest "
            "is the outgroup with the longest mean path to the ingroup taxa in "
            "that gene tree. A distant outgroup on a long branch can attach "
            "inside the ingroup in a single gene and still be the farthest, "
            "so if the closer outgroups are tangled often, compare with a run "
            "that leaves it out.",
            "- triplet.unresolved_rooted_sister_pair: often unresolved "
            "triplets/polytomies or ambiguous rooting; the engine drops these "
            "triplet/gene-tree pairs rather than classifying them.",
            "- species_triplet.invalid_rooting_or_topology: the species tree "
            "resolves no rooted sister pair for this triplet (a polytomy, or "
            "ambiguous rooting), so it cannot be normalized to A/B/C and is "
            "dropped from the run.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_preflight_data_check(
    species_tree_path,
    gene_trees_path,
    outgroups,
    output_dir=None,
    triplet_filter=None,
    species_filter=None,
    max_triplets=DEFAULT_PREFLIGHT_TRIPLET_CAP,
    max_gene_trees=DEFAULT_MAX_GENE_TREES,
    report_limit=DEFAULT_REPORT_LIMIT,
):
    """Run every structural check and write the report.

    Args:
        species_tree_path: Path to the species-tree Newick file.
        gene_trees_path: Path to the gene-trees Newick file.
        outgroups: Outgroup taxon labels; must be non-empty. Their order
            breaks distance ties only.
        output_dir: Directory to write ``preflight_data_check.txt`` into. When
            ``None``, no file is written and only the text is returned.
        triplet_filter: Optional triplet-filter file restricting the triplets.
        species_filter: Optional species-filter file; every triplet among the
            species it names is checked.
        max_triplets: Cap on generated triplets; ``0`` means no cap.
        max_gene_trees: Cap on gene trees checked; ``0`` means all.
        report_limit: Maximum example lines per issue category.

    Returns:
        A :class:`PreflightResult`.

    Raises:
        InputError: If no outgroups are given, the species-tree file does not
            hold exactly one tree, or the outgroups do not root the species
            tree (an :class:`~ghostparser.orchestrator.trees.OutgroupRootingError`,
            whose message names the taxa in the way): all conditions that
            make the check itself impossible.
    """
    if not outgroups:
        raise InputError("No outgroup taxa were provided")

    issues: list[Issue] = []

    # The species tree takes the run's own path: standardize, then root on
    # the outgroups, rank them by distance and prune them, so the ingroup,
    # every triplet's rooted shape and the rank gene trees fall back on
    # are the ones the run would see.
    species_tree_bio = standardize_tree(_load_single_species_tree(species_tree_path))
    species_labels_sorted = sorted(
        terminal.name for terminal in species_tree_bio.get_terminals()
    )
    species_rooting = _root_tree_on_outgroup(species_tree_bio, outgroups)

    species_tree_d = _to_dendropy_tree(species_rooting.tree)

    target_triplets = _load_target_triplets(
        species_labels_sorted=species_labels_sorted,
        outgroups=set(outgroups),
        triplet_filter=triplet_filter,
        species_filter=species_filter,
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
        outgroups=species_rooting.outgroup_order,
        normalized_triplets=normalized_triplets,
        max_gene_trees=max_gene_trees,
        issues=issues,
    )

    report_text = _build_report(
        species_tree_path=species_tree_path,
        gene_trees_path=gene_trees_path,
        species_rooting=species_rooting,
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
