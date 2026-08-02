"""Per-triplet GhostParser inference for the streaming pipeline.

Covers topology classification, the discordant count test, the KS tree-height
test, bootstrap aggregation (with optional debug metrics), summary-statistics
gathering, run-wide p-value correction, and TSV writing. Ported from
``triplet_processor`` (scipy/statsmodels backend only) so the pipeline stays
self-contained.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, replace

import dendropy
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

from ghostparser.triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    find_sister_pair,
)

from .config import (
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_BOOTSTRAP,
    DEFAULT_BOOTSTRAP_DEBUG_MODE,
    DEFAULT_BOOTSTRAP_ITERATIONS,
    DEFAULT_BOOTSTRAP_SUMMARY_ONLY,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    DISCORDANT_TEST_CHOICES,
    P_VALUE_CORRECTION_CHOICES,
    SUMMARY_STATISTIC_CHOICES,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
)

Classification = str
SerializedTripletObservation = tuple[str, float, dict | None]

DISCORDANT1_TOPOLOGY_CHOICES = ("BC", "AC")

_BOOTSTRAP_CLASSES = [
    "ghost_introgression",
    "inflow_introgression",
    "outflow_introgression",
    "no_introgression",
    "unresolved",
]

# Per-iteration bootstrap-debug accumulator keys, aligned with the metric tuple
# returned by ``_iteration_full`` (excluding the leading classification).
_BOOTSTRAP_DEBUG_KEYS = (
    "dct_stats",
    "dct_p_values",
    "ks_stats",
    "ks_p_values",
    "con_summaries",
    "dis_summaries",
)

_TOPOLOGY_TO_PAIR = {
    TOPOLOGY_AB: frozenset(("A", "B")),
    TOPOLOGY_BC: frozenset(("B", "C")),
    TOPOLOGY_AC: frozenset(("A", "C")),
}
_PAIR_TO_TOPOLOGY = {pair: topology for topology, pair in _TOPOLOGY_TO_PAIR.items()}

# Integer codes for the canonical topologies, used by the vectorized bootstrap.
_TOPOLOGY_CODE = {TOPOLOGY_AB: 0, TOPOLOGY_BC: 1, TOPOLOGY_AC: 2}

# Summary-statistics layout: 3 topologies x 3 metrics x 7 statistics = 63 columns.
SUMMARY_STATISTICS = ("mean", "median", "mode", "variance", "entropy", "min", "max")
SUMMARY_TOPOLOGY_LABELS = ("concordant", "discordant1", "discordant2")
SUMMARY_METRIC_LABELS = ("avg_tree_height", "internal_branch", "sister_distance")
SUMMARY_TOPOLOGY_TO_CANONICAL = {
    "concordant": TOPOLOGY_AB,
    "discordant1": TOPOLOGY_BC,
    "discordant2": TOPOLOGY_AC,
}


@dataclass(frozen=True)
class TripletPipelineResult:
    """Result of the GhostParser pipeline for one rooted species triplet.

    Carries the topology counts, DCT/KS statistics, the final classification,
    and the bootstrap aggregates for a single triplet.
    """

    triplet: tuple[str, str, str]
    species_tree: str | None
    most_frequent_matches_concordant: bool
    n_con: int
    n_dis1: int
    n_dis2: int
    dis1_topology: str | None
    dct_statistic: float
    dct_p_value: float
    dct_p_value_corrected: float
    dct_significant: bool
    ks_p_value: float | None
    ks_p_value_corrected: float | None
    ks_statistic: float | None
    ks_significant: bool | None
    summary_con: float | None
    summary_dis: float | None
    classification: Classification
    analyzed_trees: int = 0
    topology_metric_statistics: dict[str, float | None] | None = None
    bootstrap_value: float | None = None
    all_bootstrap: dict | None = None
    bootstrap_dct_stats: dict | list | None = None
    bootstrap_dct_p_value: dict | list | None = None
    bootstrap_ks_stats: dict | list | None = None
    bootstrap_ks_p_value: dict | list | None = None
    bootstrap_con_summary: dict | list | None = None
    bootstrap_dis_summary: dict | list | None = None
    bootstrap_gene_tree_heights: list[float] | None = None

    def to_dict(self):
        """Serialize the relevant triplet statistics to a dictionary.

        Returns:
            A dict of the result's fields for downstream consumers.
        """
        return {
            "triplet": self.triplet,
            "species_tree": self.species_tree,
            "most_frequent_matches_concordant": self.most_frequent_matches_concordant,
            "n_con": self.n_con,
            "n_dis1": self.n_dis1,
            "n_dis2": self.n_dis2,
            "dis1_topology": self.dis1_topology,
            "dct_statistic": self.dct_statistic,
            "dct_p_value": self.dct_p_value,
            "dct_p_value_corrected": self.dct_p_value_corrected,
            "dct_significant": self.dct_significant,
            "ks_statistic": self.ks_statistic,
            "ks_p_value": self.ks_p_value,
            "ks_p_value_corrected": self.ks_p_value_corrected,
            "ks_significant": self.ks_significant,
            "summary_con": self.summary_con,
            "summary_dis": self.summary_dis,
            "classification": self.classification,
            "analyzed_trees": self.analyzed_trees,
            "bootstrap_value": self.bootstrap_value,
            "all_bootstrap": self.all_bootstrap,
        }


def _distance_to_root(node):
    """Compute the root-to-node distance from edge lengths.

    Missing edge lengths are treated as zero.

    Args:
        node: A DendroPy node.

    Returns:
        The accumulated distance to the root as a float.
    """
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node
    return distance


def _compute_triplet_tree_metrics(
    tree,
    species_triplet=None,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Compute the selected tree-height value H(T) and optional summary metrics.

    Args:
        tree: A rooted 3-tip DendroPy tree.
        species_triplet: The ``(A, B, C)`` triplet, required for the ``A``/``B``/
            ``C`` strategies.
        tree_height_calculation_strategy: One of ``AVG``/``A``/``B``/``C``/
            ``SIS``/``INT``.
        collect_summary_statistics: When ``True``, also compute the per-tree
            ``avg_tree_height``/``internal_branch``/``sister_distance`` metrics
            used for summary-statistics gathering.

    Returns:
        A tuple ``(selected_tree_height, summary_metrics)`` where
        ``summary_metrics`` is a dict of the three metrics when
        ``collect_summary_statistics`` is ``True`` (else ``None``).

    Raises:
        ValueError: If the strategy is unsupported, the tree does not have
            exactly three tips, ``species_triplet`` is missing for A/B/C, or the
            sister-pair MRCA cannot be determined.
    """
    if tree_height_calculation_strategy not in TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    leaves = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon and leaf.taxon.label]
    if len(leaves) != 3:
        raise ValueError("Triplet tree must contain exactly 3 terminal taxa")

    leaf_by_label = {leaf.taxon.label: leaf for leaf in leaves}
    leaf_distances = {
        label: _distance_to_root(leaf) for label, leaf in leaf_by_label.items()
    }
    avg_tree_height = sum(leaf_distances.values()) / 3.0

    selected_tree_height = None
    if tree_height_calculation_strategy in {"A", "B", "C"}:
        if species_triplet is None:
            raise ValueError(
                "species_triplet is required for tree height strategies A, B, and C"
            )

        strategy_index = {"A": 0, "B": 1, "C": 2}[tree_height_calculation_strategy]
        selected_taxon_label = species_triplet[strategy_index]
        if selected_taxon_label not in leaf_distances:
            raise ValueError(
                f"Selected taxon {selected_taxon_label} not found in triplet tree"
            )
        selected_tree_height = leaf_distances[selected_taxon_label]
    elif tree_height_calculation_strategy == "AVG":
        selected_tree_height = avg_tree_height

    summary_metrics = None
    if collect_summary_statistics or tree_height_calculation_strategy in {"SIS", "INT"}:
        sister_pair = find_sister_pair(tree)
        left_label, right_label = tuple(sister_pair)
        sister_mrca = tree.mrca(taxon_labels=[left_label, right_label])
        if sister_mrca is None:
            raise ValueError("Could not determine sister-pair MRCA for triplet tree")

        internal_branch = _distance_to_root(sister_mrca)
        sister_distance = (
            leaf_distances[left_label]
            + leaf_distances[right_label]
            - 2.0 * internal_branch
        )

        if tree_height_calculation_strategy == "SIS":
            selected_tree_height = sister_distance
        elif tree_height_calculation_strategy == "INT":
            selected_tree_height = internal_branch

        if collect_summary_statistics:
            summary_metrics = {
                "avg_tree_height": avg_tree_height,
                "internal_branch": internal_branch,
                "sister_distance": sister_distance,
            }

    if selected_tree_height is None:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    return selected_tree_height, summary_metrics


def compute_tree_height_statistic(
    tree, strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY, species_triplet=None
):
    """Compute H(T) from root-to-tip distances per the selected strategy.

    Args:
        tree: A rooted 3-tip DendroPy tree.
        strategy: One of ``AVG``/``A``/``B``/``C``/``SIS``/``INT``.
        species_triplet: The ``(A, B, C)`` triplet, required for A/B/C.

    Returns:
        The tree-height value H(T) as a float.
    """
    selected_tree_height, _ = _compute_triplet_tree_metrics(
        tree,
        species_triplet=species_triplet,
        tree_height_calculation_strategy=strategy,
    )
    return selected_tree_height


def _resolve_topology_roles(topology_counts, species_topology):
    """Resolve concordant/discordant1/discordant2 roles with deterministic ties.

    Args:
        topology_counts: Mapping of canonical topology to its count.
        species_topology: The concordant (species-tree) topology.

    Returns:
        A tuple ``(con_topology, dis1_topology, dis2_topology,
        most_frequent_matches_concordant)``.

    Raises:
        ValueError: If ``species_topology`` is not a valid topology.
    """
    if species_topology not in ALL_TOPOLOGIES:
        raise ValueError("Resolved topology roles require a valid species topology")

    all_topologies = list(ALL_TOPOLOGIES)
    all_topologies.remove(species_topology)

    first, second = all_topologies[0], all_topologies[1]
    first_count = int(topology_counts.get(first, 0))
    second_count = int(topology_counts.get(second, 0))
    if first_count >= second_count:
        dis1_topology, dis2_topology = first, second
        n_dis1, n_dis2 = first_count, second_count
    else:
        dis1_topology, dis2_topology = second, first
        n_dis1, n_dis2 = second_count, first_count

    n_con = int(topology_counts.get(species_topology, 0))
    most_frequent_matches_concordant = n_con >= n_dis1 and n_con >= n_dis2

    return (
        species_topology,
        dis1_topology,
        dis2_topology,
        most_frequent_matches_concordant,
    )


def run_discordant_count_test(n_dis1, n_dis2, method="chi-square"):
    """Run the discordant count test (SciPy/statsmodels backend).

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.
        method: ``chi-square`` (SciPy Pearson chi-square) or ``z-test``
            (statsmodels two-proportion z-test).

    Returns:
        A tuple ``(statistic, p_value)``.

    Raises:
        ValueError: If the method is unsupported.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    if method == "chi-square":
        result = stats.chisquare([n_dis1, n_dis2])
        return float(result.statistic), float(result.pvalue)
    if method == "z-test":
        z_score, p_value = proportions_ztest(
            count=[n_dis1, n_dis2],
            nobs=[total, total],
            alternative="two-sided",
        )
        return float(z_score), float(p_value)
    raise ValueError(f"Unsupported discordant test method: {method}")


def run_two_sample_ks_test(sample_a, sample_b):
    """Run the two-sample KS test (SciPy backend).

    Args:
        sample_a: First sample of numeric values.
        sample_b: Second sample of numeric values.

    Returns:
        A tuple ``(D, p_value)``.
    """
    if not sample_a or not sample_b:
        return 0.0, 1.0

    result = stats.ks_2samp(sample_a, sample_b, alternative="two-sided", method="auto")
    return float(result.statistic), float(result.pvalue)


def _median(values):
    """Compute the median of a numeric iterable.

    Args:
        values: Iterable of numeric values.

    Returns:
        The median as a float, or ``None`` if ``values`` is empty.
    """
    if not values:
        return None
    sorted_vals = sorted(float(v) for v in values)
    n = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2.0


def _mean(values):
    """Compute the mean of a numeric iterable.

    Args:
        values: Iterable of numeric values.

    Returns:
        The mean as a float, or ``None`` if ``values`` is empty.
    """
    if not values:
        return None
    values_float = [float(value) for value in values]
    return sum(values_float) / len(values_float)


def _mode_binned(values, decimals=3):
    """Compute the mode after rounding values to a fixed precision.

    When multiple modes tie, the maximum mode value is returned.

    Args:
        values: Iterable of numeric values.
        decimals: Number of decimal places used for binning.

    Returns:
        The binned mode as a float, or ``None`` if ``values`` is empty.
    """
    if not values:
        return None

    counts = {}
    for value in values:
        rounded = round(float(value), decimals)
        counts[rounded] = counts.get(rounded, 0) + 1

    max_frequency = max(counts.values())
    modes = [value for value, count in counts.items() if count == max_frequency]
    return max(modes)


def _variance(values):
    """Compute the population variance of a numeric iterable.

    Args:
        values: Iterable of numeric values.

    Returns:
        The population variance as a float, or ``None`` if ``values`` is empty.
    """
    if not values:
        return None
    values_float = [float(value) for value in values]
    mean_value = _mean(values_float)
    if mean_value is None:
        return None
    return sum((value - mean_value) ** 2 for value in values_float) / len(values_float)


def _entropy_binned(values, decimals=3):
    """Compute Shannon entropy (base 2) over rounded value frequencies.

    Args:
        values: Iterable of numeric values.
        decimals: Number of decimal places used for binning.

    Returns:
        The Shannon entropy as a float, or ``None`` if ``values`` is empty.
    """
    if not values:
        return None

    counts = {}
    for value in values:
        rounded = round(float(value), decimals)
        counts[rounded] = counts.get(rounded, 0) + 1

    total = sum(counts.values())
    if total <= 0:
        return None

    entropy = 0.0
    for count in counts.values():
        probability = count / total
        if probability > 0.0:
            entropy -= probability * math.log2(probability)
    return entropy


def _compute_summary_statistic(values, statistic_name):
    """Compute one named summary statistic over a numeric iterable.

    Args:
        values: Iterable of numeric values.
        statistic_name: One of ``mean``/``median``/``mode``/``variance``/
            ``entropy``/``min``/``max``.

    Returns:
        The computed statistic as a float, or ``None`` if undefined for the
        input.

    Raises:
        ValueError: If ``statistic_name`` is unsupported.
    """
    if statistic_name == "mean":
        return _mean(values)
    if statistic_name == "median":
        return _median(values)
    if statistic_name == "mode":
        return _mode_binned(values, decimals=3)
    if statistic_name == "variance":
        return _variance(values)
    if statistic_name == "entropy":
        return _entropy_binned(values, decimals=3)
    if statistic_name == "min":
        return min(values) if values else None
    if statistic_name == "max":
        return max(values) if values else None
    raise ValueError(f"Unsupported summary statistic name: {statistic_name}")


def _summary_statistics_column_names():
    """Return the stable ordered summary-statistics column names.

    Returns:
        A list of ``<topology>_<metric>_<statistic>`` column names covering
        every topology/metric/statistic combination.
    """
    columns = []
    for topology_label in SUMMARY_TOPOLOGY_LABELS:
        for metric_label in SUMMARY_METRIC_LABELS:
            for statistic_name in SUMMARY_STATISTICS:
                columns.append(f"{topology_label}_{metric_label}_{statistic_name}")
    return columns


def _build_empty_metric_buckets():
    """Build empty per-topology metric buckets for summary-statistics gathering.

    Returns:
        A dict keyed by canonical topology, each mapping metric label to an
        empty list.
    """
    return {
        topology: {
            "avg_tree_height": [],
            "internal_branch": [],
            "sister_distance": [],
        }
        for topology in ALL_TOPOLOGIES
    }


def _build_topology_metric_statistics(
    species_triplet, species_topology, metric_buckets
):
    """Build canonical topology/metric summary statistics for one triplet.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        species_topology: The concordant (species-tree) topology.
        metric_buckets: Per-topology metric buckets from
            :func:`_build_empty_metric_buckets`.

    Returns:
        A dict mapping each ``<topology>_<metric>_<statistic>`` column name to
        its computed value (or ``None`` when undefined).
    """
    topology_counts = {
        topology: len(metric_buckets[topology]["avg_tree_height"])
        for topology in ALL_TOPOLOGIES
    }

    _, _, canonical_to_original_topology, _ = _canonicalize_triplet_labels(
        species_triplet,
        species_topology,
        topology_counts,
    )

    statistics = {}
    for topology_label in SUMMARY_TOPOLOGY_LABELS:
        canonical_topology = SUMMARY_TOPOLOGY_TO_CANONICAL[topology_label]
        original_topology = canonical_to_original_topology[canonical_topology]
        for metric_label in SUMMARY_METRIC_LABELS:
            values = metric_buckets[original_topology][metric_label]
            for statistic_name in SUMMARY_STATISTICS:
                key = f"{topology_label}_{metric_label}_{statistic_name}"
                statistics[key] = _compute_summary_statistic(values, statistic_name)

    return statistics


def _classify_introgression(dct_significant, ks_significant, summary_con, summary_dis):
    """Apply the GhostParser decision logic to produce a classification.

    Args:
        dct_significant: Whether the discordant count test is significant.
        ks_significant: Whether the KS tree-height test is significant.
        summary_con: Concordant summary statistic value.
        summary_dis: Discordant1 summary statistic value.

    Returns:
        One of ``no_introgression``, ``inflow_introgression``,
        ``outflow_introgression``, ``ghost_introgression``, or ``unresolved``.
    """
    if not dct_significant:
        return "no_introgression"
    if not ks_significant:
        return "inflow_introgression"
    if summary_con is None or summary_dis is None:
        return "unresolved"
    if summary_con > summary_dis:
        return "outflow_introgression"
    if summary_con < summary_dis:
        return "ghost_introgression"
    return "unresolved"


def _generate_inference_description(triplet, classification, dis1_topology):
    """Render a human-readable inference direction with species names.

    Args:
        triplet: The ``(A, B, C)`` triplet where A and B are species-tree
            sisters.
        classification: The triplet classification.
        dis1_topology: The discordant1 topology label (``BC`` or ``AC``).

    Returns:
        A human-readable inference description string.

    Raises:
        ValueError: If ``dis1_topology`` is not one of the expected labels.
    """
    a_taxon, b_taxon, c_taxon = triplet

    if classification == "no_introgression":
        return "no introgression"

    if dis1_topology is None:
        return classification

    topology_label = str(dis1_topology).strip().upper()

    if topology_label == "BC":
        sisters_in_dis1 = {b_taxon, c_taxon}
        outgroup_in_dis1 = a_taxon
    elif topology_label == "AC":
        sisters_in_dis1 = {a_taxon, c_taxon}
        outgroup_in_dis1 = b_taxon
    else:
        expected = " or ".join(f"'{value}'" for value in DISCORDANT1_TOPOLOGY_CHOICES)
        raise ValueError(
            f"Invalid dis1_topology '{dis1_topology}'. Expected {expected}."
        )

    species_tree_sisters = {a_taxon, b_taxon}

    if classification == "inflow_introgression":
        if c_taxon in sisters_in_dis1:
            sister_who_moved = (sisters_in_dis1 - {c_taxon}).pop()
            return f"introgression from {c_taxon} to {sister_who_moved}"
        else:
            sisters_str = " and ".join(sorted(sisters_in_dis1))
            return f"introgression between {sisters_str} and {outgroup_in_dis1}"

    elif classification == "outflow_introgression":
        if c_taxon in sisters_in_dis1:
            sister_who_introgressed = (sisters_in_dis1 - {c_taxon}).pop()
            return f"introgression from {sister_who_introgressed} to {c_taxon}"
        else:
            sisters_str = " and ".join(sorted(species_tree_sisters))
            return f"introgression from {sisters_str} to {outgroup_in_dis1}"

    elif classification == "ghost_introgression":
        return f"introgression from ghost lineage to {outgroup_in_dis1}"

    return classification


def _build_triplet_np_rng(seed, triplet):
    """Build a deterministic per-triplet NumPy RNG when a base seed is given.

    A stable per-triplet seed is derived from ``(seed, triplet)`` via SHA-256 so
    every triplet resamples independently and a run is reproducible.

    Args:
        seed: The global base seed, or ``None`` for a non-deterministic RNG.
        triplet: The triplet used to derive a stable per-triplet seed.

    Returns:
        A ``numpy.random.Generator`` instance.
    """
    if seed is None:
        return np.random.default_rng()

    triplet_key = "|".join(triplet)
    digest = hashlib.sha256(f"{seed}|{triplet_key}".encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def observation_from_subtree(
    subtree,
    triplet,
    tree_height_calculation_strategy,
    collect_summary_statistics=False,
):
    """Compute a ``(topology, tree_height, metrics)`` observation from a subtree.

    Computes the observation directly from the extracted DendroPy subtree,
    avoiding a serialize-then-reparse round trip.

    Args:
        subtree: The extracted triplet subtree as a DendroPy tree.
        triplet: The ``(A, B, C)`` triplet.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, also compute the per-tree
            summary metrics stored as the observation's third element.

    Returns:
        A ``(topology, tree_height, summary_metrics)`` tuple where
        ``summary_metrics`` is ``None`` unless ``collect_summary_statistics`` is
        ``True``, or ``None`` if the subtree's labels do not match the triplet or
        metric computation fails.
    """
    labels = {
        leaf.taxon.label
        for leaf in subtree.leaf_node_iter()
        if leaf.taxon and leaf.taxon.label
    }
    if labels != set(triplet):
        return None

    try:
        topology = classify_triplet_topology_string(subtree, triplet)
        tree_height, summary_metrics = _compute_triplet_tree_metrics(
            subtree,
            species_triplet=triplet,
            tree_height_calculation_strategy=tree_height_calculation_strategy,
            collect_summary_statistics=collect_summary_statistics,
        )
    except ValueError:
        return None

    return (topology, tree_height, summary_metrics)


def _relabel_topology(topology, old_to_new_labels):
    """Relabel a canonical topology under an old-to-new label map.

    Args:
        topology: The canonical topology to relabel.
        old_to_new_labels: Mapping of old label to new label.

    Returns:
        The relabeled canonical topology.
    """
    old_pair = _TOPOLOGY_TO_PAIR[topology]
    new_pair = frozenset(old_to_new_labels[label] for label in old_pair)
    return _PAIR_TO_TOPOLOGY[new_pair]


def _build_topology_maps(old_to_new_labels):
    """Build forward and inverse topology relabel maps for a label permutation.

    Args:
        old_to_new_labels: Mapping of old label to new label.

    Returns:
        A tuple ``(old_to_new_topology, new_to_old_topology)``.
    """
    old_to_new_topology = {}
    new_to_old_topology = {}
    for old_topology in ALL_TOPOLOGIES:
        new_topology = _relabel_topology(old_topology, old_to_new_labels)
        old_to_new_topology[old_topology] = new_topology
        new_to_old_topology[new_topology] = old_topology
    return old_to_new_topology, new_to_old_topology


def _canonicalize_triplet_labels(species_triplet, species_topology, topology_counts):
    """Canonicalize labels so concordant is AB|C and discordant1 is BC|A.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        species_topology: The concordant (species-tree) topology.
        topology_counts: Mapping of observed topology to its count.

    Returns:
        A tuple ``(canonical_triplet, canonical_counts,
        canonical_to_original_topology, reported_dis1_topology)``.

    Raises:
        ValueError: If ``species_topology`` is invalid.
    """
    if species_topology == TOPOLOGY_AB:
        base_arrangement = ("A", "B", "C")
    elif species_topology == TOPOLOGY_BC:
        base_arrangement = ("B", "C", "A")
    elif species_topology == TOPOLOGY_AC:
        base_arrangement = ("A", "C", "B")
    else:
        raise ValueError(f"Invalid species topology: {species_topology}")

    old_taxa = {
        "A": species_triplet[0],
        "B": species_triplet[1],
        "C": species_triplet[2],
    }
    canonical_triplet = tuple(old_taxa[label] for label in base_arrangement)

    old_to_new_labels = {
        old_label: new_label
        for new_label, old_label in zip(("A", "B", "C"), base_arrangement)
    }
    _, canonical_to_original_topology = _build_topology_maps(old_to_new_labels)

    canonical_counts = {
        topology: int(topology_counts.get(canonical_to_original_topology[topology], 0))
        for topology in ALL_TOPOLOGIES
    }
    reported_dis1_topology = (
        TOPOLOGY_BC
        if canonical_counts[TOPOLOGY_BC] >= canonical_counts[TOPOLOGY_AC]
        else TOPOLOGY_AC
    )

    return (
        canonical_triplet,
        canonical_counts,
        canonical_to_original_topology,
        reported_dis1_topology,
    )


def _species_tree_topology_only_newick(species_tree_newick):
    """Return the species-subtree Newick with edge lengths suppressed.

    Args:
        species_tree_newick: The species subtree Newick, or a falsy value.

    Returns:
        The topology-only Newick string, or the input unchanged when falsy.
    """
    if not species_tree_newick:
        return species_tree_newick

    tree = dendropy.Tree.get(
        data=species_tree_newick, schema="newick", preserve_underscores=True
    )
    return tree.as_string(schema="newick", suppress_edge_lengths=True).strip()


def _serialize_triplet_gene_trees(
    species_triplet,
    triplet_gene_trees,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Parse rooted triplet Newicks into observations.

    Trees whose leaf set does not match the triplet, or that fail metric
    computation, are skipped.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        triplet_gene_trees: Iterable of rooted triplet Newick strings.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, also compute each
            observation's summary metrics (third tuple element).

    Returns:
        A list of ``(topology, tree_height, summary_metrics)`` observation
        tuples.
    """
    observations: list[SerializedTripletObservation] = []
    species_set = set(species_triplet)

    for newick_str in triplet_gene_trees:
        if not str(newick_str).strip():
            continue

        tree = dendropy.Tree.get(
            data=str(newick_str).strip(), schema="newick", preserve_underscores=True
        )
        labels = {
            leaf.taxon.label
            for leaf in tree.leaf_node_iter()
            if leaf.taxon and leaf.taxon.label
        }
        if labels != species_set:
            continue

        try:
            topology = classify_triplet_topology_string(tree, species_triplet)
            tree_height, summary_metrics = _compute_triplet_tree_metrics(
                tree,
                species_triplet=species_triplet,
                tree_height_calculation_strategy=tree_height_calculation_strategy,
                collect_summary_statistics=collect_summary_statistics,
            )
        except ValueError:
            continue

        observations.append((topology, tree_height, summary_metrics))

    return observations


def _run_triplet_pipeline_from_observations(
    species_triplet,
    observations,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    species_topology=TOPOLOGY_AB,
    species_tree_newick=None,
):
    """Run the GhostParser Figure 6 pipeline from serialized observations.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        observations: List of ``(topology, tree_height)`` tuples.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        species_topology: The concordant (species-tree) topology.
        species_tree_newick: Species subtree Newick stored on the result.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.

    Raises:
        ValueError: If ``summary_statistic`` is unsupported.
    """
    heights = {topology: [] for topology in ALL_TOPOLOGIES}
    collect_summary_statistics = bool(observations) and observations[0][2] is not None
    metric_buckets = _build_empty_metric_buckets() if collect_summary_statistics else None
    for topology, tree_height, summary_metrics in observations:
        heights[topology].append(tree_height)
        if metric_buckets is not None and summary_metrics is not None:
            for metric_label in SUMMARY_METRIC_LABELS:
                metric_buckets[topology][metric_label].append(
                    summary_metrics[metric_label]
                )

    topology_counts = {topology: len(heights[topology]) for topology in ALL_TOPOLOGIES}
    analyzed_trees = sum(topology_counts.values())

    (
        canonical_triplet,
        canonical_counts,
        canonical_to_original_topology,
        reported_dis1_topology,
    ) = _canonicalize_triplet_labels(
        species_triplet,
        species_topology,
        topology_counts,
    )
    canonical_heights = {
        topology: heights[canonical_to_original_topology[topology]]
        for topology in ALL_TOPOLOGIES
    }

    con_topology, dis1_topology, dis2_topology, most_frequent_matches_concordant = (
        _resolve_topology_roles(
            canonical_counts,
            TOPOLOGY_AB,
        )
    )

    n_con = canonical_counts[con_topology]
    n_dis1 = canonical_counts[dis1_topology]
    n_dis2 = canonical_counts[dis2_topology]

    if summary_statistic not in SUMMARY_STATISTIC_CHOICES:
        raise ValueError(
            f"Unsupported summary statistic: {summary_statistic}. "
            f"Choose one of: {', '.join(SUMMARY_STATISTIC_CHOICES)}"
        )

    dct_statistic, dct_p_value = run_discordant_count_test(
        n_dis1,
        n_dis2,
        method=discordant_test,
    )
    dct_significant = dct_p_value <= alpha_dct

    ks_statistic, ks_p_value = run_two_sample_ks_test(
        canonical_heights[dis1_topology],
        canonical_heights[con_topology],
    )
    ks_significant = ks_p_value <= alpha_ks

    if summary_statistic == "mean":
        summary_con = _mean(canonical_heights[con_topology])
        summary_dis = _mean(canonical_heights[dis1_topology])
    elif summary_statistic == "mode":
        summary_con = _mode_binned(canonical_heights[con_topology], decimals=3)
        summary_dis = _mode_binned(canonical_heights[dis1_topology], decimals=3)
    else:
        summary_con = _median(canonical_heights[con_topology])
        summary_dis = _median(canonical_heights[dis1_topology])

    classification = _classify_introgression(
        dct_significant,
        ks_significant,
        summary_con,
        summary_dis,
    )

    topology_metric_statistics = None
    if metric_buckets is not None:
        topology_metric_statistics = _build_topology_metric_statistics(
            canonical_triplet, TOPOLOGY_AB, metric_buckets
        )

    return TripletPipelineResult(
        triplet=canonical_triplet,
        species_tree=species_tree_newick,
        most_frequent_matches_concordant=most_frequent_matches_concordant,
        n_con=n_con,
        n_dis1=n_dis1,
        n_dis2=n_dis2,
        dis1_topology="BC" if reported_dis1_topology == TOPOLOGY_BC else "AC",
        dct_statistic=dct_statistic,
        dct_p_value=dct_p_value,
        dct_p_value_corrected=dct_p_value,
        dct_significant=dct_significant,
        ks_p_value=ks_p_value,
        ks_p_value_corrected=ks_p_value,
        ks_statistic=ks_statistic,
        ks_significant=ks_significant,
        summary_con=summary_con,
        summary_dis=summary_dis,
        classification=classification,
        analyzed_trees=analyzed_trees,
        topology_metric_statistics=topology_metric_statistics,
    )


def _iteration_classification(
    n_dis1,
    n_dis2,
    con_heights,
    dis1_heights,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
):
    """Classify one bootstrap iteration from its resampled per-topology heights.

    Mirrors ``_run_triplet_pipeline_from_observations`` plus the unresolved
    fallback: an iteration with an empty concordant or discordant1 sample is
    ``unresolved`` (its summary statistic would be undefined).

    Args:
        n_dis1: Discordant1 count in the resample.
        n_dis2: Discordant2 count in the resample.
        con_heights: Concordant tree heights (Python list).
        dis1_heights: Discordant1 tree heights (Python list).
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.

    Returns:
        The iteration classification string.
    """
    if not con_heights or not dis1_heights:
        return "unresolved"

    _, dct_p_value = run_discordant_count_test(
        n_dis1, n_dis2, method=discordant_test
    )
    dct_significant = dct_p_value <= alpha_dct

    _, ks_p_value = run_two_sample_ks_test(dis1_heights, con_heights)
    ks_significant = ks_p_value <= alpha_ks

    if summary_statistic == "mean":
        summary_con = _mean(con_heights)
        summary_dis = _mean(dis1_heights)
    elif summary_statistic == "mode":
        summary_con = _mode_binned(con_heights, decimals=3)
        summary_dis = _mode_binned(dis1_heights, decimals=3)
    else:
        summary_con = _median(con_heights)
        summary_dis = _median(dis1_heights)

    return _classify_introgression(
        dct_significant, ks_significant, summary_con, summary_dis
    )


def _iteration_full(
    n_dis1,
    n_dis2,
    con_heights,
    dis1_heights,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
):
    """Classify one bootstrap iteration and return its per-iteration metrics.

    Unlike :func:`_iteration_classification` (the fast path), this always
    computes the DCT/KS statistics and summary values so they can be collected
    for bootstrap-debug output.

    Args:
        n_dis1: Discordant1 count in the resample.
        n_dis2: Discordant2 count in the resample.
        con_heights: Concordant tree heights (Python list).
        dis1_heights: Discordant1 tree heights (Python list).
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.

    Returns:
        A tuple ``(classification, dct_statistic, dct_p_value, ks_statistic,
        ks_p_value, summary_con, summary_dis)``.
    """
    dct_statistic, dct_p_value = run_discordant_count_test(
        n_dis1, n_dis2, method=discordant_test
    )
    dct_significant = dct_p_value <= alpha_dct

    ks_statistic, ks_p_value = run_two_sample_ks_test(dis1_heights, con_heights)
    ks_significant = ks_p_value <= alpha_ks

    if summary_statistic == "mean":
        summary_con = _mean(con_heights)
        summary_dis = _mean(dis1_heights)
    elif summary_statistic == "mode":
        summary_con = _mode_binned(con_heights, decimals=3)
        summary_dis = _mode_binned(dis1_heights, decimals=3)
    else:
        summary_con = _median(con_heights)
        summary_dis = _median(dis1_heights)

    if not con_heights or not dis1_heights:
        classification = "unresolved"
    else:
        classification = _classify_introgression(
            dct_significant, ks_significant, summary_con, summary_dis
        )

    return (
        classification,
        dct_statistic,
        dct_p_value,
        ks_statistic,
        ks_p_value,
        summary_con,
        summary_dis,
    )


def _numeric_summary(values):
    """Build a compact summary of an iterable that may contain ``None`` values.

    Args:
        values: Iterable of numeric values (or ``None`` entries).

    Returns:
        A dict with ``count``/``non_null_count``/``mean``/``median``/``min``/
        ``max``.
    """
    cleaned = [float(value) for value in values if value is not None]
    if not cleaned:
        return {
            "count": len(values),
            "non_null_count": 0,
            "mean": None,
            "median": None,
            "min": None,
            "max": None,
        }

    return {
        "count": len(values),
        "non_null_count": len(cleaned),
        "mean": _mean(cleaned),
        "median": _median(cleaned),
        "min": min(cleaned),
        "max": max(cleaned),
    }


def _finalize_bootstrap_metric(values, summary_only):
    """Return a compact summary of a bootstrap metric, or the raw per-iteration list.

    Args:
        values: Per-iteration metric values.
        summary_only: When ``True``, return a compact numeric summary instead of
            the full list.

    Returns:
        A numeric-summary dict when ``summary_only`` is ``True``, else ``values``
        unchanged.
    """
    if summary_only:
        return _numeric_summary(values)
    return values


def _serialize_bootstrap_value(value):
    """Serialize a bootstrap-debug structure for TSV output as strict JSON.

    Args:
        value: The bootstrap-debug payload (list or dict), or ``None``.

    Returns:
        A compact JSON string, or an empty string when ``value`` is ``None``.

    Raises:
        ValueError: If ``value`` is not JSON-serializable.
    """
    if value is None:
        return ""
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Bootstrap payload is not JSON-serializable: {value!r}"
        ) from exc


def _append_iteration_debug(debug, metrics):
    """Append one iteration's debug metrics to the accumulator lists.

    Args:
        debug: Accumulator dict keyed by :data:`_BOOTSTRAP_DEBUG_KEYS`.
        metrics: The six per-iteration metric values, in
            :data:`_BOOTSTRAP_DEBUG_KEYS` order.
    """
    for key, value in zip(_BOOTSTRAP_DEBUG_KEYS, metrics):
        debug[key].append(value)


def _bootstrap_payload(bootstrap_value, fractions, debug, summary_only):
    """Assemble the bootstrap payload with optional finalized debug metrics.

    Args:
        bootstrap_value: The top class fraction.
        fractions: Per-class fractions.
        debug: The per-iteration debug accumulator, or ``None`` when disabled.
        summary_only: When ``True``, finalize debug metrics as compact summaries.

    Returns:
        A dict with ``bootstrap_value``, ``all_bootstrap``, and the six
        ``bootstrap_*`` debug entries (each ``None`` when ``debug`` is ``None``).
    """
    payload = {"bootstrap_value": bootstrap_value, "all_bootstrap": fractions}
    debug_columns = (
        "bootstrap_dct_stats",
        "bootstrap_dct_p_value",
        "bootstrap_ks_stats",
        "bootstrap_ks_p_value",
        "bootstrap_con_summary",
        "bootstrap_dis_summary",
    )
    if debug is None:
        for column in debug_columns:
            payload[column] = None
    else:
        for column, key in zip(debug_columns, _BOOTSTRAP_DEBUG_KEYS):
            payload[column] = _finalize_bootstrap_metric(debug[key], summary_only)
    return payload


def _run_bootstrap_iterations(
    observations,
    iterations,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
    rng,
    debug_mode=False,
    summary_only=False,
):
    """Resample observations and aggregate per-iteration classifications.

    The resample indices are drawn with NumPy's vectorized RNG and the topology
    counts and per-topology height groups are computed with array operations,
    replacing the per-element Python resampling loop. Assumes the concordant
    topology is ``((A,B),C)`` (always true on the pipeline path). When
    ``debug_mode`` is set, per-iteration DCT/KS statistics and con/dis summaries
    are collected too.

    Args:
        observations: List of ``(topology, tree_height, metrics)`` tuples.
        iterations: Number of bootstrap iterations.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        rng: A ``numpy.random.Generator`` used for resampling.
        debug_mode: When ``True``, collect per-iteration debug metrics.
        summary_only: When ``True`` (and debug), emit compact summaries instead
            of full per-iteration lists.

    Returns:
        A dict with ``bootstrap_value`` (top class fraction), ``all_bootstrap``
        (per-class fractions), and the six ``bootstrap_*`` debug entries (each
        ``None`` unless ``debug_mode``).
    """
    class_counts = {label: 0 for label in _BOOTSTRAP_CLASSES}
    debug = {key: [] for key in _BOOTSTRAP_DEBUG_KEYS} if debug_mode else None

    def _record(classification):
        if classification not in class_counts:
            class_counts[classification] = 0
        class_counts[classification] += 1

    if iterations <= 0:
        fractions = {label: 0.0 for label in sorted(class_counts)}
        return _bootstrap_payload(0.0, fractions, debug, summary_only)

    valid_count = len(observations)
    if valid_count == 0:
        for _ in range(iterations):
            if debug_mode:
                metrics = _iteration_full(
                    0, 0, [], [], alpha_dct, alpha_ks, discordant_test,
                    summary_statistic,
                )
                _record(metrics[0])
                _append_iteration_debug(debug, metrics[1:])
            else:
                _record("unresolved")
        fractions = {
            label: class_counts.get(label, 0) / float(iterations)
            for label in sorted(class_counts)
        }
        return _bootstrap_payload(
            max(fractions.values()), fractions, debug, summary_only
        )

    topo_codes = np.fromiter(
        (_TOPOLOGY_CODE[topology] for topology, _, _ in observations),
        dtype=np.int64,
        count=valid_count,
    )
    heights = np.fromiter(
        (height for _, height, _ in observations),
        dtype=np.float64,
        count=valid_count,
    )

    for _ in range(iterations):
        row = rng.integers(0, valid_count, size=valid_count)
        sampled_topo = topo_codes[row]
        sampled_heights = heights[row]

        counts = np.bincount(sampled_topo, minlength=3)
        n_bc = int(counts[1])
        n_ac = int(counts[2])
        if n_bc >= n_ac:
            dis1_code, n_dis1, n_dis2 = 1, n_bc, n_ac
        else:
            dis1_code, n_dis1, n_dis2 = 2, n_ac, n_bc

        con_heights = sampled_heights[sampled_topo == 0].tolist()
        dis1_heights = sampled_heights[sampled_topo == dis1_code].tolist()

        if debug_mode:
            metrics = _iteration_full(
                n_dis1,
                n_dis2,
                con_heights,
                dis1_heights,
                alpha_dct,
                alpha_ks,
                discordant_test,
                summary_statistic,
            )
            _record(metrics[0])
            _append_iteration_debug(debug, metrics[1:])
        else:
            _record(
                _iteration_classification(
                    n_dis1,
                    n_dis2,
                    con_heights,
                    dis1_heights,
                    alpha_dct,
                    alpha_ks,
                    discordant_test,
                    summary_statistic,
                )
            )

    fractions = {
        label: class_counts.get(label, 0) / float(iterations)
        for label in sorted(class_counts)
    }
    top_fraction = max(fractions.values()) if fractions else 0.0

    return _bootstrap_payload(top_fraction, fractions, debug, summary_only)


def _finalize_triplet_analysis(
    triplet,
    observations,
    species_tree_topology,
    *,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
    bootstrap_options,
    triplet_seed,
):
    """Run the base pipeline and bootstrap for one triplet's observations.

    Shared core of :func:`analyze_triplet` and
    :func:`analyze_triplet_from_observations`.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        observations: List of ``(topology, tree_height)`` tuples.
        species_tree_topology: Topology-only species subtree Newick, or ``None``.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic bootstrap.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values and bootstrap
        aggregates.
    """
    base_result = _run_triplet_pipeline_from_observations(
        triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        species_topology=TOPOLOGY_AB,
        species_tree_newick=species_tree_topology,
    )

    options = dict(bootstrap_options or {})
    iterations = int(options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    debug_mode = bool(options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE))
    summary_only = (
        bool(options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY))
        if debug_mode
        else False
    )
    rng = _build_triplet_np_rng(triplet_seed, triplet)

    bootstrap_payload = _run_bootstrap_iterations(
        observations,
        iterations=iterations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        rng=rng,
        debug_mode=debug_mode,
        summary_only=summary_only,
    )

    bootstrap_gene_tree_heights = None
    if debug_mode:
        raw_heights = [tree_height for _, tree_height, _ in observations]
        bootstrap_gene_tree_heights = _finalize_bootstrap_metric(
            raw_heights, summary_only
        )

    return replace(
        base_result,
        bootstrap_value=bootstrap_payload["bootstrap_value"],
        all_bootstrap=bootstrap_payload["all_bootstrap"],
        bootstrap_dct_stats=bootstrap_payload["bootstrap_dct_stats"],
        bootstrap_dct_p_value=bootstrap_payload["bootstrap_dct_p_value"],
        bootstrap_ks_stats=bootstrap_payload["bootstrap_ks_stats"],
        bootstrap_ks_p_value=bootstrap_payload["bootstrap_ks_p_value"],
        bootstrap_con_summary=bootstrap_payload["bootstrap_con_summary"],
        bootstrap_dis_summary=bootstrap_payload["bootstrap_dis_summary"],
        bootstrap_gene_tree_heights=bootstrap_gene_tree_heights,
    )


def analyze_triplet_from_observations(
    triplet,
    observations,
    species_subtree=None,
    *,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    bootstrap_options=None,
    triplet_seed=None,
):
    """Analyze one triplet from precomputed observations.

    This is the pipeline's per-triplet unit: observations are computed once
    during extraction (see :func:`observation_from_subtree`), so no Newick
    reparse happens here. The returned p-values are uncorrected; run-wide
    correction is applied later.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        observations: List of ``(topology, tree_height)`` tuples.
        species_subtree: The triplet's species subtree Newick (with branch
            lengths); stored topology-only on the result.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic bootstrap.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.
    """
    species_tree_topology = _species_tree_topology_only_newick(species_subtree)
    return _finalize_triplet_analysis(
        triplet,
        observations,
        species_tree_topology,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
    )


def analyze_triplet(
    triplet,
    gene_subtrees,
    species_subtree=None,
    *,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
    bootstrap_options=None,
    triplet_seed=None,
):
    """Analyze one triplet end-to-end from gene-subtree Newick strings.

    Serializes the gene subtrees into observations (reparsing each Newick), then
    runs the base pipeline and bootstrap. The returned p-values are uncorrected;
    run-wide correction is applied later.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        gene_subtrees: Iterable of the triplet's rooted gene-subtree Newicks.
        species_subtree: The triplet's species subtree Newick (with branch
            lengths); stored topology-only on the result.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, gather per-triplet
            topology/metric summary statistics onto the result.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic bootstrap.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.
    """
    species_tree_topology = _species_tree_topology_only_newick(species_subtree)
    observations = _serialize_triplet_gene_trees(
        triplet,
        gene_subtrees,
        tree_height_calculation_strategy=tree_height_calculation_strategy,
        collect_summary_statistics=collect_summary_statistics,
    )
    return _finalize_triplet_analysis(
        triplet,
        observations,
        species_tree_topology,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
    )


def _adjust_p_values(
    p_values,
    method=DEFAULT_P_VALUE_CORRECTION,
    alpha=0.05,
):
    """Adjust p-values by the selected correction method (statsmodels backend).

    Args:
        p_values: List of raw p-values.
        method: One of the supported correction methods (or ``no``).
        alpha: FDR level used by TSBH and statsmodels.

    Returns:
        The adjusted p-value list.

    Raises:
        ValueError: If the method is unsupported.
    """
    if method not in P_VALUE_CORRECTION_CHOICES:
        raise ValueError(
            f"Unsupported p-value correction method: {method}. "
            f"Choose one of: {', '.join(P_VALUE_CORRECTION_CHOICES)}"
        )

    if method == "no":
        return [float(p_value) for p_value in p_values]

    method_map = {
        "bfn": "bonferroni",
        "holm": "holm",
        "fdr_bh": "fdr_bh",
        "fdr_by": "fdr_by",
        "fdr_tsbh": "fdr_tsbh",
    }
    mapped_method = method_map[method]
    _, corrected, _, _ = multipletests(
        [float(p_value) for p_value in p_values],
        alpha=alpha,
        method=mapped_method,
    )
    return [float(p_value) for p_value in corrected]


def _apply_triplet_result_p_value_correction(
    results,
    alpha_dct,
    alpha_ks,
    method=DEFAULT_P_VALUE_CORRECTION,
):
    """Apply run-wide p-value correction to DCT and KS p-values.

    Correction is applied once across all triplets; each result's
    classification and bootstrap value are recomputed from the corrected
    significance.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        method: Correction method (or ``no``).

    Returns:
        A new list of corrected ``TripletPipelineResult`` objects.
    """
    if not results:
        return results

    dct_p_values = [result.dct_p_value for result in results]
    adjusted_dct = _adjust_p_values(
        dct_p_values,
        method=method,
        alpha=alpha_dct,
    )

    ks_indices = [
        idx for idx, result in enumerate(results) if result.ks_p_value is not None
    ]
    ks_p_values = [results[idx].ks_p_value for idx in ks_indices]
    adjusted_ks_values = _adjust_p_values(
        ks_p_values,
        method=method,
        alpha=alpha_ks,
    )
    adjusted_ks_map = {
        idx: adjusted_ks_values[pos] for pos, idx in enumerate(ks_indices)
    }

    adjusted_results = []
    for idx, result in enumerate(results):
        dct_p_value_corrected = adjusted_dct[idx]
        dct_significant = dct_p_value_corrected <= alpha_dct

        if result.ks_p_value is None:
            ks_p_value_corrected = None
            ks_significant = None
        else:
            ks_p_value_corrected = adjusted_ks_map[idx]
            ks_significant = ks_p_value_corrected <= alpha_ks

        classification = _classify_introgression(
            dct_significant,
            ks_significant,
            result.summary_con,
            result.summary_dis,
        )
        bootstrap_value = result.bootstrap_value
        if result.all_bootstrap is not None:
            bootstrap_value = float(result.all_bootstrap.get(classification, 0.0))

        adjusted_results.append(
            replace(
                result,
                dct_p_value_corrected=dct_p_value_corrected,
                dct_significant=dct_significant,
                ks_p_value_corrected=ks_p_value_corrected,
                ks_significant=ks_significant,
                classification=classification,
                bootstrap_value=bootstrap_value,
            )
        )

    return adjusted_results


def _format_all_bootstrap(value):
    """Format the all_bootstrap classification-fraction dict as key=value pairs.

    Args:
        value: The ``all_bootstrap`` dict, or ``None``.

    Returns:
        A comma-separated ``classification=fraction`` string, or an empty string
        when ``value`` is ``None``.

    Raises:
        ValueError: If ``value`` is not a dict.
    """
    if value is None:
        return ""
    if not isinstance(value, dict):
        raise ValueError(f"Bootstrap payload is not JSON-serializable: {value!r}")
    return ",".join(f"{k}={v:.12g}" for k, v in sorted(value.items()))


def write_pipeline_results(
    results,
    output_filepath,
    dct_method=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_debug_mode=DEFAULT_BOOTSTRAP_DEBUG_MODE,
):
    """Write per-triplet results to a TSV file.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        output_filepath: Destination TSV path.
        dct_method: Discordant test method, used to name the statistic column.
        summary_statistic: Summary statistic, used to name the con/dis columns.
        p_value_correction: Correction method, used to name the corrected
            columns.
        bootstrap: Whether to include the bootstrap columns.
        bootstrap_debug_mode: Whether to also include the bootstrap-debug
            columns (only meaningful when ``bootstrap`` is ``True``).

    Raises:
        ValueError: If any of the method/statistic/correction values is
            unsupported.
    """
    if dct_method not in DISCORDANT_TEST_CHOICES:
        raise ValueError(
            f"Unsupported discordant test method: {dct_method}. "
            f"Choose one of: {', '.join(DISCORDANT_TEST_CHOICES)}"
        )

    if summary_statistic not in SUMMARY_STATISTIC_CHOICES:
        raise ValueError(
            f"Unsupported summary statistic: {summary_statistic}. "
            f"Choose one of: {', '.join(SUMMARY_STATISTIC_CHOICES)}"
        )

    if p_value_correction not in P_VALUE_CORRECTION_CHOICES:
        raise ValueError(
            f"Unsupported p-value correction method: {p_value_correction}. "
            f"Choose one of: {', '.join(P_VALUE_CORRECTION_CHOICES)}"
        )

    if dct_method == "z-test":
        dct_column = "dct_z_score"
    else:
        dct_column = "dct_chi_stats"

    summary_con_column = f"{summary_statistic}_con"
    summary_dis_column = f"{summary_statistic}_dis"
    bootstrap_con_column = f"bootstrap_con_{summary_statistic}"
    bootstrap_dis_column = f"bootstrap_dis_{summary_statistic}"
    dct_corrected_column = f"dct_p_val_{p_value_correction}_corr"
    ks_corrected_column = f"ks_p_val_{p_value_correction}_corr"

    header = [
        "triplet",
        "abc_mapping",
        "species_tree",
        "n_con",
        "n_dis1",
        "n_dis2",
        "dis1_topology",
        "most_frequent_matches_concordant",
        dct_column,
        "dct_p_value",
        dct_corrected_column,
        "dct_significant",
        "ks_statistic",
        "ks_p_value",
        ks_corrected_column,
        "ks_significant",
        summary_con_column,
        summary_dis_column,
        "classification",
        "inference",
        "analyzed_trees",
    ]

    if bootstrap:
        header.extend(
            [
                "bootstrap_value",
                "all_bootstrap",
            ]
        )
        if bootstrap_debug_mode:
            header.extend(
                [
                    "bootstrap_dct_stats",
                    "bootstrap_dct_p_value",
                    "bootstrap_ks_stats",
                    "bootstrap_ks_p_value",
                    bootstrap_con_column,
                    bootstrap_dis_column,
                    "bootstrap_gene_tree_heights",
                ]
            )

    with open(output_filepath, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for result in results:
            a_taxon, b_taxon, c_taxon = result.triplet
            abc_mapping = f"A={a_taxon};B={b_taxon};C={c_taxon}"
            inference = _generate_inference_description(
                result.triplet,
                result.classification,
                result.dis1_topology,
            )
            row = [
                ",".join(result.triplet),
                abc_mapping,
                "" if result.species_tree is None else result.species_tree,
                str(result.n_con),
                str(result.n_dis1),
                str(result.n_dis2),
                "" if result.dis1_topology is None else result.dis1_topology,
                str(result.most_frequent_matches_concordant),
                f"{result.dct_statistic:.12g}",
                f"{result.dct_p_value:.12g}",
                f"{result.dct_p_value_corrected:.12g}",
                str(result.dct_significant),
                "" if result.ks_statistic is None else f"{result.ks_statistic:.12g}",
                "" if result.ks_p_value is None else f"{result.ks_p_value:.12g}",
                ""
                if result.ks_p_value_corrected is None
                else f"{result.ks_p_value_corrected:.12g}",
                "" if result.ks_significant is None else str(result.ks_significant),
                "" if result.summary_con is None else f"{result.summary_con:.12g}",
                "" if result.summary_dis is None else f"{result.summary_dis:.12g}",
                result.classification,
                inference,
                str(result.analyzed_trees),
            ]

            if bootstrap:
                row.extend(
                    [
                        ""
                        if result.bootstrap_value is None
                        else f"{result.bootstrap_value:.12g}",
                        _format_all_bootstrap(result.all_bootstrap),
                    ]
                )
                if bootstrap_debug_mode:
                    row.extend(
                        [
                            _serialize_bootstrap_value(result.bootstrap_dct_stats),
                            _serialize_bootstrap_value(result.bootstrap_dct_p_value),
                            _serialize_bootstrap_value(result.bootstrap_ks_stats),
                            _serialize_bootstrap_value(result.bootstrap_ks_p_value),
                            _serialize_bootstrap_value(result.bootstrap_con_summary),
                            _serialize_bootstrap_value(result.bootstrap_dis_summary),
                            _serialize_bootstrap_value(
                                result.bootstrap_gene_tree_heights
                            ),
                        ]
                    )

            out_f.write("\t".join(row) + "\n")


def write_summary_statistics_tsv(results, output_filepath, bootstrap=DEFAULT_BOOTSTRAP):
    """Write per-triplet topology/metric summary statistics to a TSV file.

    Emits the identity fields, the 63 topology/metric/statistic columns, the
    classification, and (when bootstrap is enabled) the bootstrap value.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        output_filepath: Destination TSV path.
        bootstrap: Whether to include the ``bootstrap_value`` column.
    """
    include_bootstrap_value = bool(bootstrap)

    header = [
        "triplet",
        "abc_mapping",
        "species_tree",
        "dis1_topology",
        "n_con",
        "n_dis1",
        "n_dis2",
    ]
    header.extend(_summary_statistics_column_names())
    header.append("classification")
    if include_bootstrap_value:
        header.append("bootstrap_value")

    with open(output_filepath, "w") as out_f:
        out_f.write("\t".join(header) + "\n")
        for result in results:
            a_taxon, b_taxon, c_taxon = result.triplet
            abc_mapping = f"A={a_taxon};B={b_taxon};C={c_taxon}"

            row = [
                ",".join(result.triplet),
                abc_mapping,
                "" if result.species_tree is None else result.species_tree,
                "" if result.dis1_topology is None else result.dis1_topology,
                str(result.n_con),
                str(result.n_dis1),
                str(result.n_dis2),
            ]

            topology_metric_statistics = result.topology_metric_statistics or {}
            for column_name in _summary_statistics_column_names():
                value = topology_metric_statistics.get(column_name)
                row.append("" if value is None else f"{value:.12g}")

            row.append(result.classification)
            if include_bootstrap_value:
                row.append(
                    ""
                    if result.bootstrap_value is None
                    else f"{result.bootstrap_value:.12g}"
                )

            out_f.write("\t".join(row) + "\n")
