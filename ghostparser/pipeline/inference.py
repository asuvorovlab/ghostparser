"""Per-triplet GhostParser inference for the streaming pipeline.

Covers topology classification, the discordant count test, the KS tree-height
test, bootstrap aggregation, run-wide p-value correction, and TSV writing.
Ported from ``triplet_processor`` (non-summary-statistics path only) so the
pipeline stays self-contained.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, replace

import dendropy
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

from ghostparser.config import (
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_BOOTSTRAP,
    DEFAULT_BOOTSTRAP_ITERATIONS,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_STATS_BACKEND,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    DISCORDANT_TEST_CHOICES,
    P_VALUE_CORRECTION_CHOICES,
    STATS_BACKEND_CHOICES,
    SUMMARY_STATISTIC_CHOICES,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
)
from ghostparser.triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    find_sister_pair,
)

Classification = str
SerializedTripletObservation = tuple[str, float]

DISCORDANT1_TOPOLOGY_CHOICES = ("BC", "AC")

_BOOTSTRAP_CLASSES = [
    "ghost_introgression",
    "inflow_introgression",
    "outflow_introgression",
    "no_introgression",
    "unresolved",
]

_TOPOLOGY_TO_PAIR = {
    TOPOLOGY_AB: frozenset(("A", "B")),
    TOPOLOGY_BC: frozenset(("B", "C")),
    TOPOLOGY_AC: frozenset(("A", "C")),
}
_PAIR_TO_TOPOLOGY = {pair: topology for topology, pair in _TOPOLOGY_TO_PAIR.items()}

# Integer codes for the canonical topologies, used by the vectorized bootstrap.
_TOPOLOGY_CODE = {TOPOLOGY_AB: 0, TOPOLOGY_BC: 1, TOPOLOGY_AC: 2}


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
):
    """Compute the selected tree-height value H(T) for one triplet tree.

    Args:
        tree: A rooted 3-tip DendroPy tree.
        species_triplet: The ``(A, B, C)`` triplet, required for the ``A``/``B``/
            ``C`` strategies.
        tree_height_calculation_strategy: One of ``AVG``/``A``/``B``/``C``/
            ``SIS``/``INT``.

    Returns:
        The selected tree-height value H(T) as a float.

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

    if tree_height_calculation_strategy in {"SIS", "INT"}:
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

    if selected_tree_height is None:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    return selected_tree_height


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
    return _compute_triplet_tree_metrics(
        tree,
        species_triplet=species_triplet,
        tree_height_calculation_strategy=strategy,
    )


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


def pearson_discordant_chi_square_test(n_dis1, n_dis2):
    """Run the custom Pearson chi-square test for discordant count imbalance.

    Uses equal expected frequencies with the analytic df=1 p-value.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.

    Returns:
        A tuple ``(statistic, p_value)``.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    expected = total / 2.0
    chi2_stat = ((n_dis1 - expected) ** 2) / expected + (
        (n_dis2 - expected) ** 2
    ) / expected

    p_value = math.erfc(math.sqrt(chi2_stat / 2.0))
    return float(chi2_stat), float(p_value)


def _pearson_discordant_chi_square_test_scipy(n_dis1, n_dis2):
    """Run the SciPy Pearson chi-square backend for discordant counts.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.

    Returns:
        A tuple ``(statistic, p_value)``.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    result = stats.chisquare([n_dis1, n_dis2])
    return float(result.statistic), float(result.pvalue)


def _two_proportion_discordant_z_test_statsmodels(n_dis1, n_dis2):
    """Run the statsmodels two-proportion z-test backend for discordant counts.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.

    Returns:
        A tuple ``(z_score, p_value)``.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    z_score, p_value = proportions_ztest(
        count=[n_dis1, n_dis2],
        nobs=[total, total],
        alternative="two-sided",
    )
    return float(z_score), float(p_value)


def two_proportion_discordant_z_test(n_dis1, n_dis2):
    """Run the custom two-proportion z-test for discordant count imbalance.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.

    Returns:
        A tuple ``(z_score, p_value)``.
    """
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    p_dis1 = n_dis1 / total
    p_dis2 = n_dis2 / total
    pooled = (n_dis1 + n_dis2) / (2 * total)
    standard_error = (pooled * (1.0 - pooled) * ((1.0 / total) + (1.0 / total))) ** 0.5
    if standard_error == 0.0:
        return 0.0, 1.0

    z_score = (p_dis1 - p_dis2) / standard_error
    p_value = 2.0 * stats.norm.sf(abs(z_score))
    return float(z_score), float(p_value)


def run_discordant_count_test(
    n_dis1, n_dis2, method="chi-square", stats_backend="custom"
):
    """Dispatch the discordant count test by method and backend.

    Args:
        n_dis1: Count of the discordant1 topology.
        n_dis2: Count of the discordant2 topology.
        method: ``chi-square`` or ``z-test``.
        stats_backend: ``custom`` or ``standard``.

    Returns:
        A tuple ``(statistic, p_value)``.

    Raises:
        ValueError: If the backend or method is unsupported.
    """
    if stats_backend not in STATS_BACKEND_CHOICES:
        raise ValueError(
            f"Unsupported stats backend: {stats_backend}. "
            f"Choose one of: {', '.join(STATS_BACKEND_CHOICES)}"
        )

    if method == "chi-square":
        if stats_backend == "standard":
            return _pearson_discordant_chi_square_test_scipy(n_dis1, n_dis2)
        return pearson_discordant_chi_square_test(n_dis1, n_dis2)
    if method == "z-test":
        if stats_backend == "standard":
            return _two_proportion_discordant_z_test_statsmodels(n_dis1, n_dis2)
        return two_proportion_discordant_z_test(n_dis1, n_dis2)
    raise ValueError(f"Unsupported discordant test method: {method}")


def two_sample_ks_test(sample_a, sample_b):
    """Run the custom two-sample KS test with an asymptotic p-value.

    Args:
        sample_a: First sample of numeric values.
        sample_b: Second sample of numeric values.

    Returns:
        A tuple ``(D, p_value)``.
    """
    if not sample_a or not sample_b:
        return 0.0, 1.0

    data_a = sorted(float(value) for value in sample_a)
    data_b = sorted(float(value) for value in sample_b)
    n1 = len(data_a)
    n2 = len(data_b)

    i = 0
    j = 0
    cdf_a = 0.0
    cdf_b = 0.0
    d_stat = 0.0

    while i < n1 and j < n2:
        a_val = data_a[i]
        b_val = data_b[j]
        if a_val <= b_val:
            while i < n1 and data_a[i] == a_val:
                i += 1
            cdf_a = i / n1
        if b_val <= a_val:
            while j < n2 and data_b[j] == b_val:
                j += 1
            cdf_b = j / n2
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    while i < n1:
        i += 1
        cdf_a = i / n1
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    while j < n2:
        j += 1
        cdf_b = j / n2
        d_stat = max(d_stat, abs(cdf_a - cdf_b))

    en = (n1 * n2) / (n1 + n2)
    if en <= 0:
        return float(d_stat), 1.0

    sqrt_en = math.sqrt(en)
    lam = (sqrt_en + 0.12 + 0.11 / sqrt_en) * d_stat

    if lam <= 0:
        p_value = 1.0
    else:
        series_sum = 0.0
        for k in range(1, 101):
            term = math.exp(-2.0 * (k**2) * (lam**2))
            if k % 2 == 1:
                series_sum += term
            else:
                series_sum -= term
            if term < 1e-12:
                break
        p_value = max(0.0, min(1.0, 2.0 * series_sum))

    return float(d_stat), float(p_value)


def _two_sample_ks_test_scipy(sample_a, sample_b):
    """Run the SciPy two-sample KS backend.

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


def run_two_sample_ks_test(sample_a, sample_b, stats_backend="custom"):
    """Dispatch the two-sample KS test by backend.

    Args:
        sample_a: First sample of numeric values.
        sample_b: Second sample of numeric values.
        stats_backend: ``custom`` or ``standard``.

    Returns:
        A tuple ``(D, p_value)``.

    Raises:
        ValueError: If the backend is unsupported.
    """
    if stats_backend not in STATS_BACKEND_CHOICES:
        raise ValueError(
            f"Unsupported stats backend: {stats_backend}. "
            f"Choose one of: {', '.join(STATS_BACKEND_CHOICES)}"
        )

    if stats_backend == "standard":
        return _two_sample_ks_test_scipy(sample_a, sample_b)
    return two_sample_ks_test(sample_a, sample_b)


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


def observation_from_subtree(subtree, triplet, tree_height_calculation_strategy):
    """Compute a ``(topology, tree_height)`` observation from a subtree object.

    Computes the observation directly from the extracted DendroPy subtree,
    avoiding a serialize-then-reparse round trip.

    Args:
        subtree: The extracted triplet subtree as a DendroPy tree.
        triplet: The ``(A, B, C)`` triplet.
        tree_height_calculation_strategy: Tree-height strategy to apply.

    Returns:
        A ``(topology, tree_height)`` tuple, or ``None`` if the subtree's labels
        do not match the triplet or metric computation fails.
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
        tree_height = _compute_triplet_tree_metrics(
            subtree,
            species_triplet=triplet,
            tree_height_calculation_strategy=tree_height_calculation_strategy,
        )
    except ValueError:
        return None

    return (topology, tree_height)


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
):
    """Parse rooted triplet Newicks into (topology, tree-height) observations.

    Trees whose leaf set does not match the triplet, or that fail metric
    computation, are skipped.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        triplet_gene_trees: Iterable of rooted triplet Newick strings.
        tree_height_calculation_strategy: Tree-height strategy to apply.

    Returns:
        A list of ``(topology, tree_height)`` observation tuples.
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
            tree_height = _compute_triplet_tree_metrics(
                tree,
                species_triplet=species_triplet,
                tree_height_calculation_strategy=tree_height_calculation_strategy,
            )
        except ValueError:
            continue

        observations.append((topology, tree_height))

    return observations


def _run_triplet_pipeline_from_observations(
    species_triplet,
    observations,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    stats_backend=DEFAULT_STATS_BACKEND,
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
        stats_backend: ``custom`` or ``standard``.
        species_topology: The concordant (species-tree) topology.
        species_tree_newick: Species subtree Newick stored on the result.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.

    Raises:
        ValueError: If ``summary_statistic`` is unsupported.
    """
    heights = {topology: [] for topology in ALL_TOPOLOGIES}
    for topology, tree_height in observations:
        heights[topology].append(tree_height)

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
        stats_backend=stats_backend,
    )
    dct_significant = dct_p_value <= alpha_dct

    ks_statistic, ks_p_value = run_two_sample_ks_test(
        canonical_heights[dis1_topology],
        canonical_heights[con_topology],
        stats_backend=stats_backend,
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
    stats_backend,
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
        stats_backend: ``custom`` or ``standard``.

    Returns:
        The iteration classification string.
    """
    if not con_heights or not dis1_heights:
        return "unresolved"

    _, dct_p_value = run_discordant_count_test(
        n_dis1, n_dis2, method=discordant_test, stats_backend=stats_backend
    )
    dct_significant = dct_p_value <= alpha_dct

    _, ks_p_value = run_two_sample_ks_test(
        dis1_heights, con_heights, stats_backend=stats_backend
    )
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


def _run_bootstrap_iterations(
    observations,
    iterations,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
    stats_backend,
    rng,
):
    """Resample observations and aggregate per-iteration classifications.

    The resample indices are drawn with NumPy's vectorized RNG and the topology
    counts and per-topology height groups are computed with array operations,
    replacing the per-element Python resampling loop. Assumes the concordant
    topology is ``((A,B),C)`` (always true on the pipeline path).

    Args:
        observations: List of ``(topology, tree_height)`` tuples.
        iterations: Number of bootstrap iterations.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        summary_statistic: ``mean``, ``median``, or ``mode``.
        stats_backend: ``custom`` or ``standard``.
        rng: A ``numpy.random.Generator`` used for resampling.

    Returns:
        A dict with ``bootstrap_value`` (top class fraction) and
        ``all_bootstrap`` (per-class fractions).
    """
    class_counts = {label: 0 for label in _BOOTSTRAP_CLASSES}

    if iterations <= 0:
        fractions = {label: 0.0 for label in sorted(class_counts)}
        return {"bootstrap_value": 0.0, "all_bootstrap": fractions}

    valid_count = len(observations)
    if valid_count == 0:
        # Every iteration has zero analyzed trees, i.e. unresolved.
        class_counts["unresolved"] = iterations
        fractions = {
            label: class_counts.get(label, 0) / float(iterations)
            for label in sorted(class_counts)
        }
        return {"bootstrap_value": max(fractions.values()), "all_bootstrap": fractions}

    topo_codes = np.fromiter(
        (_TOPOLOGY_CODE[topology] for topology, _ in observations),
        dtype=np.int64,
        count=valid_count,
    )
    heights = np.fromiter(
        (height for _, height in observations), dtype=np.float64, count=valid_count
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

        classification = _iteration_classification(
            n_dis1,
            n_dis2,
            con_heights,
            dis1_heights,
            alpha_dct,
            alpha_ks,
            discordant_test,
            summary_statistic,
            stats_backend,
        )
        if classification not in class_counts:
            class_counts[classification] = 0
        class_counts[classification] += 1

    fractions = {
        label: class_counts.get(label, 0) / float(iterations)
        for label in sorted(class_counts)
    }
    top_fraction = max(fractions.values()) if fractions else 0.0

    return {"bootstrap_value": top_fraction, "all_bootstrap": fractions}


def _finalize_triplet_analysis(
    triplet,
    observations,
    species_tree_topology,
    *,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
    stats_backend,
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
        stats_backend: ``custom`` or ``standard``.
        bootstrap_options: Optional dict; only ``iterations`` is read.
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
        stats_backend=stats_backend,
        species_topology=TOPOLOGY_AB,
        species_tree_newick=species_tree_topology,
    )

    options = dict(bootstrap_options or {})
    iterations = int(options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    rng = _build_triplet_np_rng(triplet_seed, triplet)

    bootstrap_payload = _run_bootstrap_iterations(
        observations,
        iterations=iterations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        rng=rng,
    )

    return replace(
        base_result,
        bootstrap_value=bootstrap_payload["bootstrap_value"],
        all_bootstrap=bootstrap_payload["all_bootstrap"],
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
    stats_backend=DEFAULT_STATS_BACKEND,
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
        stats_backend: ``custom`` or ``standard``.
        bootstrap_options: Optional dict; only ``iterations`` is read.
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
        stats_backend=stats_backend,
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
    stats_backend=DEFAULT_STATS_BACKEND,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
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
        stats_backend: ``custom`` or ``standard``.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        bootstrap_options: Optional dict; only ``iterations`` is read.
        triplet_seed: Optional base seed for deterministic bootstrap.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.
    """
    species_tree_topology = _species_tree_topology_only_newick(species_subtree)
    observations = _serialize_triplet_gene_trees(
        triplet,
        gene_subtrees,
        tree_height_calculation_strategy=tree_height_calculation_strategy,
    )
    return _finalize_triplet_analysis(
        triplet,
        observations,
        species_tree_topology,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
    )


def _bonferroni_adjust_p_values_custom(p_values):
    """Apply Bonferroni correction to a p-value list.

    Args:
        p_values: List of raw p-values.

    Returns:
        The adjusted p-value list.
    """
    m = len(p_values)
    if m == 0:
        return []
    return [min(1.0, float(p_value) * m) for p_value in p_values]


def _holm_adjust_p_values_custom(p_values):
    """Apply Holm step-down FWER correction to a p-value list.

    Args:
        p_values: List of raw p-values.

    Returns:
        The adjusted p-value list.
    """
    m = len(p_values)
    if m == 0:
        return []

    indexed = sorted(enumerate(float(p) for p in p_values), key=lambda item: item[1])
    adjusted_sorted = [0.0] * m
    running_max = 0.0
    for idx, (_, p_value) in enumerate(indexed):
        scaled = min(1.0, (m - idx) * p_value)
        running_max = max(running_max, scaled)
        adjusted_sorted[idx] = running_max

    adjusted = [0.0] * m
    for sorted_idx, (original_idx, _) in enumerate(indexed):
        adjusted[original_idx] = adjusted_sorted[sorted_idx]
    return adjusted


def _fdr_bh_adjust_p_values_custom(p_values):
    """Apply Benjamini-Hochberg FDR correction to a p-value list.

    Args:
        p_values: List of raw p-values.

    Returns:
        The adjusted p-value list.
    """
    m = len(p_values)
    if m == 0:
        return []

    indexed = sorted(enumerate(float(p) for p in p_values), key=lambda item: item[1])
    adjusted_sorted = [0.0] * m

    prev = 1.0
    for idx in range(m - 1, -1, -1):
        _, p_value = indexed[idx]
        rank = idx + 1
        adjusted = min(1.0, (p_value * m) / rank)
        prev = min(prev, adjusted)
        adjusted_sorted[idx] = prev

    adjusted = [0.0] * m
    for sorted_idx, (original_idx, _) in enumerate(indexed):
        adjusted[original_idx] = adjusted_sorted[sorted_idx]

    return adjusted


def _fdr_by_adjust_p_values_custom(p_values):
    """Apply Benjamini-Yekutieli FDR correction to a p-value list.

    Args:
        p_values: List of raw p-values.

    Returns:
        The adjusted p-value list.
    """
    m = len(p_values)
    if m == 0:
        return []

    c_m = sum(1.0 / j for j in range(1, m + 1))
    indexed = sorted(enumerate(float(p) for p in p_values), key=lambda item: item[1])
    adjusted_sorted = [0.0] * m

    prev = 1.0
    for idx in range(m - 1, -1, -1):
        _, p_value = indexed[idx]
        rank = idx + 1
        adjusted = min(1.0, (p_value * m * c_m) / rank)
        prev = min(prev, adjusted)
        adjusted_sorted[idx] = prev

    adjusted = [0.0] * m
    for sorted_idx, (original_idx, _) in enumerate(indexed):
        adjusted[original_idx] = adjusted_sorted[sorted_idx]
    return adjusted


def _fdr_tsbh_adjust_p_values_custom(p_values, alpha=0.05):
    """Apply two-stage Benjamini-Hochberg (TSBH) FDR correction.

    Args:
        p_values: List of raw p-values.
        alpha: FDR level used to estimate the true-null count.

    Returns:
        The adjusted p-value list.

    Raises:
        ValueError: If ``alpha`` is not in the open interval (0, 1).
    """
    m = len(p_values)
    if m == 0:
        return []

    if alpha <= 0 or alpha >= 1:
        raise ValueError("alpha must be in (0, 1) for fdr_tsbh")

    indexed = sorted(enumerate(float(p) for p in p_values), key=lambda item: item[1])
    p_sorted = [p_value for _, p_value in indexed]

    alpha_stage1 = alpha / (1.0 + alpha)
    rejects_stage1 = 0
    for idx, p_value in enumerate(p_sorted):
        rank = idx + 1
        if p_value <= (rank / m) * alpha_stage1:
            rejects_stage1 = rank

    m0_hat = max(1, m - rejects_stage1)
    adjusted_sorted = [0.0] * m
    prev = 1.0
    for idx in range(m - 1, -1, -1):
        rank = idx + 1
        adjusted = min(1.0, (p_sorted[idx] * m0_hat) / rank)
        prev = min(prev, adjusted)
        adjusted_sorted[idx] = prev

    adjusted = [0.0] * m
    for sorted_idx, (original_idx, _) in enumerate(indexed):
        adjusted[original_idx] = adjusted_sorted[sorted_idx]
    return adjusted


def _adjust_p_values(
    p_values,
    method=DEFAULT_P_VALUE_CORRECTION,
    stats_backend=DEFAULT_STATS_BACKEND,
    alpha=0.05,
):
    """Adjust p-values by the selected correction method and backend.

    Args:
        p_values: List of raw p-values.
        method: One of the supported correction methods (or ``no``).
        stats_backend: ``custom`` or ``standard``.
        alpha: FDR level used by TSBH and the statsmodels backend.

    Returns:
        The adjusted p-value list.

    Raises:
        ValueError: If the method or backend is unsupported.
    """
    if method not in P_VALUE_CORRECTION_CHOICES:
        raise ValueError(
            f"Unsupported p-value correction method: {method}. "
            f"Choose one of: {', '.join(P_VALUE_CORRECTION_CHOICES)}"
        )

    if stats_backend not in STATS_BACKEND_CHOICES:
        raise ValueError(
            f"Unsupported stats backend: {stats_backend}. "
            f"Choose one of: {', '.join(STATS_BACKEND_CHOICES)}"
        )

    if method == "no":
        return [float(p_value) for p_value in p_values]

    if stats_backend == "standard":
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

    if method == "bfn":
        return _bonferroni_adjust_p_values_custom(p_values)
    if method == "holm":
        return _holm_adjust_p_values_custom(p_values)
    if method == "fdr_bh":
        return _fdr_bh_adjust_p_values_custom(p_values)
    if method == "fdr_by":
        return _fdr_by_adjust_p_values_custom(p_values)
    return _fdr_tsbh_adjust_p_values_custom(p_values, alpha=alpha)


def _apply_triplet_result_p_value_correction(
    results,
    alpha_dct,
    alpha_ks,
    method=DEFAULT_P_VALUE_CORRECTION,
    stats_backend=DEFAULT_STATS_BACKEND,
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
        stats_backend: ``custom`` or ``standard``.

    Returns:
        A new list of corrected ``TripletPipelineResult`` objects.
    """
    if not results:
        return results

    dct_p_values = [result.dct_p_value for result in results]
    adjusted_dct = _adjust_p_values(
        dct_p_values,
        method=method,
        stats_backend=stats_backend,
        alpha=alpha_dct,
    )

    ks_indices = [
        idx for idx, result in enumerate(results) if result.ks_p_value is not None
    ]
    ks_p_values = [results[idx].ks_p_value for idx in ks_indices]
    adjusted_ks_values = _adjust_p_values(
        ks_p_values,
        method=method,
        stats_backend=stats_backend,
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

            out_f.write("\t".join(row) + "\n")
