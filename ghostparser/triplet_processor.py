"""Triplet processing and GhostParser decision pipeline.

This module implements the sequential logic for rooted species triplets:

1. Classify each triplet gene tree as concordant/dis1/dis2.
2. Compute tree height statistic ``H(T)`` using configurable strategy
    (average root-to-tip distance by default).
3. Run discordant count test (default: Pearson chi-square, alpha=0.05).
4. If significant, run tree height test (two-sample KS, alpha=0.05).
5. If significant, compare selected summary values to classify inflow vs ghost introgression.

Concordant topology is defined as the species-tree topology for the ABC triplet.
The two discordant topologies are frequency-ranked as ``dis1`` (more frequent)
and ``dis2`` (less frequent), with deterministic first-discordant tie-breaking.

The discordant count test is configurable: two-proportion z-test (default) or
Pearson chi-square.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import random
import time
from dataclasses import dataclass, replace
from multiprocessing import cpu_count
from pathlib import Path

import dendropy
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

from .config import (
    DEFAULT_ALPHA_DCT,
    DEFAULT_ALPHA_KS,
    DEFAULT_BOOTSTRAP,
    DEFAULT_BOOTSTRAP_DEBUG_MODE,
    DEFAULT_BOOTSTRAP_ITERATIONS,
    DEFAULT_BOOTSTRAP_SUMMARY_ONLY,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_GENERATE_SUMMARY_STATS,
    DEFAULT_INPUT_FORMAT,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_STATS_BACKEND,
    DEFAULT_SUMMARY_STATISTIC,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    DISCORDANT_TEST_CHOICES,
    INPUT_FORMAT_CHOICES,
    P_VALUE_CORRECTION_CHOICES,
    STATS_BACKEND_CHOICES,
    SUMMARY_STATISTIC_CHOICES,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
)
from .triplet_utils import (
    ALL_TOPOLOGIES,
    TOPOLOGY_AB,
    TOPOLOGY_AC,
    TOPOLOGY_BC,
    classify_triplet_topology_string,
    find_sister_pair,
)

Classification = str
SerializedTripletObservation = tuple[str, float]
MetricBuckets = dict[str, dict[str, list[float]]]


_TRIPLET_ENTRY_MAP: dict[tuple[str, str, str], dict] = {}

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
    """Result of running the GhostParser pipeline for one rooted species triplet."""

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
    topology_metric_statistics: dict[str, float | None] | None = None

    def to_dict(self):
        """Serialize all relevant triplet statistics to a dictionary."""
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
            "bootstrap_dct_stats": self.bootstrap_dct_stats,
            "bootstrap_dct_p_value": self.bootstrap_dct_p_value,
            "bootstrap_ks_stats": self.bootstrap_ks_stats,
            "bootstrap_ks_p_value": self.bootstrap_ks_p_value,
            "bootstrap_con_summary": self.bootstrap_con_summary,
            "bootstrap_dis_summary": self.bootstrap_dis_summary,
            "bootstrap_gene_tree_heights": self.bootstrap_gene_tree_heights,
            "topology_metric_statistics": self.topology_metric_statistics,
        }


def _distance_to_root(node):
    """Compute root-to-node distance using edge lengths (missing lengths treated as 0)."""
    distance = 0.0
    current = node
    while current is not None and current.parent_node is not None:
        edge_length = current.edge_length
        if edge_length is not None:
            distance += float(edge_length)
        current = current.parent_node
    return distance


def compute_tree_height_statistic(
    tree, strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY, species_triplet=None
):
    """Compute H(T) from root-to-tip distances according to selected strategy.

    For a rooted triplet ``((X:b2,Y:b3):b4,Z:b1)``, this equals:

    ``H(T) = (b1 + b2 + b3 + 2*b4) / 3``.

    Strategy options:
    - ``AVG``: mean over all three tip distances
    - ``A``/``B``/``C``: distance of the corresponding taxon in ``species_triplet``
    - ``SIS``: pairwise distance between the rooted sister taxa in the current topology
    - ``INT``: internal branch from sister-pair MRCA to triplet root
    """
    selected_tree_height, _ = _compute_triplet_tree_metrics(
        tree,
        species_triplet=species_triplet,
        tree_height_calculation_strategy=strategy,
        collect_summary_statistics=False,
    )
    return selected_tree_height


def _compute_triplet_tree_metrics(
    tree,
    species_triplet=None,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Compute the selected tree-height value and optional summary metrics for one tree."""
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


def classify_triplet_topology(
    tree,
    species_triplet,
    topology_counts,
    species_topology=TOPOLOGY_AB,
):
    """Classify topology as concordant/discordant1/discordant2.

    Discordant1 and discordant2 are frequency-ranked among discordant
    topologies. Returns both the per-tree label and whether concordant is
    most frequent overall.
    """
    _, dis1_topology, dis2_topology, most_frequent_matches_concordant = (
        _resolve_topology_roles(
            topology_counts,
            species_topology,
        )
    )

    topology = classify_triplet_topology_string(tree, species_triplet)
    if topology == species_topology:
        return "concordant", most_frequent_matches_concordant
    if topology == dis1_topology:
        return "discordant1", most_frequent_matches_concordant
    if topology == dis2_topology:
        return "discordant2", most_frequent_matches_concordant
    raise ValueError("Unknown triplet topology")


def _resolve_topology_roles(topology_counts, species_topology):
    """Resolve concordant/discordant roles with deterministic tie handling.

    Concordant is the topology matching ``species_topology``. Discordant roles are
    selected from the remaining two topologies by descending count; ties keep
    list order.
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
    """Custom Pearson chi-square test for discordant topology count imbalance.

    This matches the SciPy default setup for two categories with equal expected
    frequencies and computes the p-value analytically for df=1.
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
    """Backup SciPy Pearson chi-square implementation."""
    total = n_dis1 + n_dis2
    if total == 0:
        return 0.0, 1.0

    result = stats.chisquare([n_dis1, n_dis2])
    return float(result.statistic), float(result.pvalue)


def _two_proportion_discordant_z_test_statsmodels(n_dis1, n_dis2):
    """Standard-library two-proportion z-test via statsmodels."""
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
    """Two-proportion z-test for discordant count imbalance.

    Uses pooled standard error with two proportions defined over the same
    discordant-total denominator for ``dis1`` and ``dis2``.
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
    """Run selected discordant count test and return (statistic, p-value)."""
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


def run_two_sample_ks_test(sample_a, sample_b, stats_backend="custom"):
    """Run selected two-sample KS backend and return (statistic, p-value)."""
    if stats_backend not in STATS_BACKEND_CHOICES:
        raise ValueError(
            f"Unsupported stats backend: {stats_backend}. "
            f"Choose one of: {', '.join(STATS_BACKEND_CHOICES)}"
        )

    if stats_backend == "standard":
        return _two_sample_ks_test_scipy(sample_a, sample_b)
    return two_sample_ks_test(sample_a, sample_b)


def two_sample_ks_test(sample_a, sample_b):
    """Custom two-sample KS test returning (D, p-value).

    Uses the standard two-sided asymptotic Kolmogorov approximation for p-value:
    Q_K(lambda) = 2 * sum_{j>=1} (-1)^(j-1) exp(-2*j^2*lambda^2),
    with finite-sample lambda correction.
    """
    if not sample_a or not sample_b:
        return 0.0, 1.0

    # TODO: NumPy optimization: cast sample_a/sample_b to float64 arrays and use vectorized unique/sort/CDF-diff operations for KS D-stat computation.
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
    """Backup SciPy two-sample KS implementation."""
    if not sample_a or not sample_b:
        return 0.0, 1.0

    result = stats.ks_2samp(sample_a, sample_b, alternative="two-sided", method="auto")
    return float(result.statistic), float(result.pvalue)


def two_sample_ks_test_hybrid(sample_a, sample_b, alpha=0.05, borderline_margin=0.01):
    """Hybrid KS helper for future use (not used in active pipeline).

    Strategy:
    - Run custom KS first for speed.
    - If p-value is close to decision boundary (``alpha``), recompute using
      SciPy backup implementation.

    Args:
        sample_a: First sample.
        sample_b: Second sample.
        alpha: Decision threshold used to define the boundary neighborhood.
        borderline_margin: Width of boundary neighborhood around ``alpha``.

    Returns:
        Tuple ``(D, p_value)``.
    """
    if borderline_margin < 0:
        raise ValueError("borderline_margin must be >= 0")

    d_stat, p_value = two_sample_ks_test(sample_a, sample_b)
    if abs(p_value - alpha) <= borderline_margin:
        return _two_sample_ks_test_scipy(sample_a, sample_b)
    return d_stat, p_value


def _median(values):
    """Compute median of numeric iterable."""
    if not values:
        return None
    # TODO: NumPy optimization: store tree-height vectors as np.ndarray and use np.median(values) directly.
    sorted_vals = sorted(float(v) for v in values)
    n = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2.0


def _mean(values):
    """Compute mean of numeric iterable."""
    if not values:
        return None
    # TODO: NumPy optimization: keep values in np.ndarray[float64] and compute mean via np.mean(values).
    values_float = [float(value) for value in values]
    return sum(values_float) / len(values_float)


def _mode_binned(values, decimals=3):
    """Compute mode after rounding values to a fixed decimal precision.

    If multiple modes are present, return the maximum mode value.
    """
    if not values:
        return None

    # TODO: NumPy optimization: use np.round + np.unique(return_counts=True) to compute binned mode without Python dict loops.
    counts = {}
    for value in values:
        rounded = round(float(value), decimals)
        counts[rounded] = counts.get(rounded, 0) + 1

    max_frequency = max(counts.values())
    modes = [value for value, count in counts.items() if count == max_frequency]
    return max(modes)


def _variance(values):
    """Compute population variance of numeric iterable."""
    if not values:
        return None
    # TODO: NumPy optimization: compute variance using np.var on float64 arrays to reduce Python-loop overhead.
    values_float = [float(value) for value in values]
    mean_value = _mean(values_float)
    if mean_value is None:
        return None
    return sum((value - mean_value) ** 2 for value in values_float) / len(values_float)


def _entropy_binned(values, decimals=3):
    """Compute Shannon entropy (base 2) on rounded value frequencies."""
    if not values:
        return None

    # TODO: NumPy optimization: use np.round + np.unique(return_counts=True) for binning and
    # compute entropy from normalized counts (or pass probabilities to scipy.stats.entropy).
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
    """Compute one summary statistic by name."""
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
    """Return stable ordered column names for topology/metric summary statistics."""
    columns = []
    for topology_label in SUMMARY_TOPOLOGY_LABELS:
        for metric_label in SUMMARY_METRIC_LABELS:
            for statistic_name in SUMMARY_STATISTICS:
                columns.append(f"{topology_label}_{metric_label}_{statistic_name}")
    return columns


def _build_empty_metric_buckets():
    """Build empty metric buckets keyed by canonical topology and metric type."""
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
    """Build canonical topology/metric summary statistics for one triplet."""
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
    """Apply GhostParser decision logic to produce final classification."""
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
    """Generate human-readable inference direction with actual species names.

    Args:
        triplet: Tuple of (A, B, C) taxa where A and B are sisters in species tree.
        classification: String classification (no_introgression, inflow_introgression, etc).
        dis1_topology: Discordant1 topology label ("BC" or "AC").

    Returns:
        Human-readable inference description with actual species names.
    """
    a_taxon, b_taxon, c_taxon = triplet

    if classification == "no_introgression":
        return "no introgression"

    if dis1_topology is None:
        return classification  # Fallback for unresolved cases.

    # Parse dis1_topology to determine sister pairs.
    topology_label = str(dis1_topology).strip().upper()

    if topology_label == "BC":  # ((B,C),A) - B and C are sisters
        sisters_in_dis1 = {b_taxon, c_taxon}
        outgroup_in_dis1 = a_taxon
    elif topology_label == "AC":  # ((A,C),B) - A and C are sisters
        sisters_in_dis1 = {a_taxon, c_taxon}
        outgroup_in_dis1 = b_taxon
    else:
        expected = " or ".join(f"'{value}'" for value in DISCORDANT1_TOPOLOGY_CHOICES)
        raise ValueError(
            f"Invalid dis1_topology '{dis1_topology}'. Expected {expected}."
        )

    # Species tree sisters are A and B
    species_tree_sisters = {a_taxon, b_taxon}

    if classification == "inflow_introgression":
        # Introgression FROM non-sister (outgroup in species tree) TO sisters in dis1.
        # C (outgroup in species tree) moved toward A or B.
        if c_taxon in sisters_in_dis1:
            sister_who_moved = (sisters_in_dis1 - {c_taxon}).pop()
            return f"introgression from {c_taxon} to {sister_who_moved}"
        else:
            # Fallback: describe based on what we know
            sisters_str = " and ".join(sorted(sisters_in_dis1))
            return f"introgression between {sisters_str} and {outgroup_in_dis1}"

    elif classification == "outflow_introgression":
        # Introgression FROM sisters in species tree TO the non-sister (C).
        # One of A or B moved toward C.
        if c_taxon in sisters_in_dis1:
            # C paired with one of A or B
            sister_who_introgressed = (sisters_in_dis1 - {c_taxon}).pop()
            return f"introgression from {sister_who_introgressed} to {c_taxon}"
        else:
            # Fallback
            sisters_str = " and ".join(sorted(species_tree_sisters))
            return f"introgression from {sisters_str} to {outgroup_in_dis1}"

    elif classification == "ghost_introgression":
        # Use the dis1 outgroup taxon as recipient per reporting convention.
        return f"introgression from ghost lineage to {outgroup_in_dis1}"

    return classification


_BOOTSTRAP_CLASSES = [
    "ghost_introgression",
    "inflow_introgression",
    "outflow_introgression",
    "no_introgression",
    "unresolved",
]

DISCORDANT1_TOPOLOGY_CHOICES = ("BC", "AC")


def _build_triplet_rng(seed, triplet):
    """Build a deterministic per-triplet RNG when a base seed is provided."""
    if seed is None:
        return random.Random()

    triplet_key = "|".join(triplet)
    # Keep this simple: one stable string seed per (global seed, triplet).
    return random.Random(f"{seed}|{triplet_key}")


def _bootstrap_iteration_classification(iteration_result):
    """Apply unresolved fallback rule for bootstrap iterations."""
    if (
        iteration_result.analyzed_trees == 0
        or iteration_result.summary_con is None
        or iteration_result.summary_dis is None
        or iteration_result.ks_statistic is None
        or iteration_result.ks_p_value is None
    ):
        return "unresolved"
    return iteration_result.classification


def _numeric_summary(values):
    """Build compact summary statistics for an iterable with possible null values."""
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
    if summary_only:
        return _numeric_summary(values)
    return values


def _run_bootstrap_iterations(
    species_triplet,
    observations,
    iterations,
    alpha_dct,
    alpha_ks,
    discordant_test,
    summary_statistic,
    stats_backend,
    species_topology,
    species_tree_newick,
    debug_mode,
    summary_only,
    rng,
):
    """Run bootstrap iterations and return aggregate bootstrap outputs."""
    valid_count = len(observations)
    class_counts = {label: 0 for label in _BOOTSTRAP_CLASSES}

    dct_stats = [] if debug_mode else None
    dct_p_values = [] if debug_mode else None
    ks_stats = [] if debug_mode else None
    ks_p_values = [] if debug_mode else None
    con_summaries = [] if debug_mode else None
    dis_summaries = [] if debug_mode else None

    for _ in range(iterations):
        if valid_count == 0:
            sampled_observations = []
        else:
            sampled_observations = [
                observations[rng.randrange(valid_count)] for _ in range(valid_count)
            ]

        iter_result = _run_triplet_pipeline_from_observations(
            species_triplet,
            sampled_observations,
            alpha_dct=alpha_dct,
            alpha_ks=alpha_ks,
            discordant_test=discordant_test,
            summary_statistic=summary_statistic,
            stats_backend=stats_backend,
            species_topology=species_topology,
            species_tree_newick=species_tree_newick,
        )

        iter_classification = _bootstrap_iteration_classification(iter_result)
        if iter_classification not in class_counts:
            class_counts[iter_classification] = 0
        class_counts[iter_classification] += 1

        if debug_mode:
            dct_stats.append(iter_result.dct_statistic)
            dct_p_values.append(iter_result.dct_p_value)
            ks_stats.append(iter_result.ks_statistic)
            ks_p_values.append(iter_result.ks_p_value)
            con_summaries.append(iter_result.summary_con)
            dis_summaries.append(iter_result.summary_dis)

    if iterations <= 0:
        fractions = {label: 0.0 for label in sorted(class_counts)}
        return {
            "bootstrap_value": 0.0,
            "all_bootstrap": fractions,
            "bootstrap_dct_stats": _finalize_bootstrap_metric(dct_stats, summary_only)
            if debug_mode
            else None,
            "bootstrap_dct_p_value": _finalize_bootstrap_metric(
                dct_p_values, summary_only
            )
            if debug_mode
            else None,
            "bootstrap_ks_stats": _finalize_bootstrap_metric(ks_stats, summary_only)
            if debug_mode
            else None,
            "bootstrap_ks_p_value": _finalize_bootstrap_metric(
                ks_p_values, summary_only
            )
            if debug_mode
            else None,
            "bootstrap_con_summary": _finalize_bootstrap_metric(
                con_summaries, summary_only
            )
            if debug_mode
            else None,
            "bootstrap_dis_summary": _finalize_bootstrap_metric(
                dis_summaries, summary_only
            )
            if debug_mode
            else None,
        }

    fractions = {
        label: class_counts.get(label, 0) / float(iterations)
        for label in sorted(class_counts)
    }
    top_fraction = max(fractions.values()) if fractions else 0.0

    return {
        "bootstrap_value": top_fraction,
        "all_bootstrap": fractions,
        "bootstrap_dct_stats": _finalize_bootstrap_metric(dct_stats, summary_only)
        if debug_mode
        else None,
        "bootstrap_dct_p_value": _finalize_bootstrap_metric(dct_p_values, summary_only)
        if debug_mode
        else None,
        "bootstrap_ks_stats": _finalize_bootstrap_metric(ks_stats, summary_only)
        if debug_mode
        else None,
        "bootstrap_ks_p_value": _finalize_bootstrap_metric(ks_p_values, summary_only)
        if debug_mode
        else None,
        "bootstrap_con_summary": _finalize_bootstrap_metric(con_summaries, summary_only)
        if debug_mode
        else None,
        "bootstrap_dis_summary": _finalize_bootstrap_metric(dis_summaries, summary_only)
        if debug_mode
        else None,
    }


def _serialize_bootstrap_value(value):
    """Serialize bootstrap structures for TSV output as strict JSON."""
    if value is None:
        return ""

    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Bootstrap payload is not JSON-serializable: {value!r}"
        ) from exc


def _format_all_bootstrap(value):
    """Format the all_bootstrap classification-fraction dict as readable key=value pairs.

    Produces a comma-separated string of ``classification=fraction`` entries
    sorted by classification name, e.g.::

        ghost_introgression=0.42,no_introgression=0.58

    Returns an empty string when value is None.
    Raises ValueError for non-dict or non-serializable payloads (same contract
    as _serialize_bootstrap_value).
    """
    if value is None:
        return ""
    if not isinstance(value, dict):
        raise ValueError(f"Bootstrap payload is not JSON-serializable: {value!r}")
    return ",".join(f"{k}={v:.12g}" for k, v in sorted(value.items()))


def _bonferroni_adjust_p_values_custom(p_values):
    """Apply Bonferroni correction to a p-value list."""
    m = len(p_values)
    if m == 0:
        return []
    # TODO: NumPy optimization: compute corrected p-values with np.minimum(1.0, p_values_array * m).
    return [min(1.0, float(p_value) * m) for p_value in p_values]


def _holm_adjust_p_values_custom(p_values):
    """Apply Holm step-down FWER correction to a p-value list."""
    m = len(p_values)
    if m == 0:
        return []

    # TODO: NumPy optimization: use np.argsort + vectorized Holm scaling + np.maximum.accumulate for the step-down monotonic pass.
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
    """Apply Benjamini-Hochberg FDR correction to a p-value list."""
    m = len(p_values)
    if m == 0:
        return []

    # TODO: NumPy optimization: replace indexed Python loops with np.argsort + vectorized rank scaling + reverse cumulative minimum.
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
    """Apply Benjamini-Yekutieli FDR correction to a p-value list."""
    m = len(p_values)
    if m == 0:
        return []

    # TODO: NumPy optimization: compute harmonic factor and BY rank scaling with vectorized arrays instead of Python loops.
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
    """Apply two-stage Benjamini-Hochberg (TSBH) FDR correction to a p-value list."""
    m = len(p_values)
    if m == 0:
        return []

    if alpha <= 0 or alpha >= 1:
        raise ValueError("alpha must be in (0, 1) for fdr_tsbh")

    # TODO: NumPy optimization: vectorize both TSBH stages (rank-threshold rejection scan and reverse cumulative-min adjustment) on sorted arrays.
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
    """Adjust p-values using selected correction method and backend."""
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
    """Apply selected p-value correction across all triplets for DCT and KS p-values."""
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


def apply_triplet_result_p_value_correction(
    results,
    alpha_dct,
    alpha_ks,
    method=DEFAULT_P_VALUE_CORRECTION,
    stats_backend=DEFAULT_STATS_BACKEND,
):
    """Public wrapper for applying global p-value correction to triplet results."""
    return _apply_triplet_result_p_value_correction(
        results,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        method=method,
        stats_backend=stats_backend,
    )


_TOPOLOGY_TO_PAIR = {
    TOPOLOGY_AB: frozenset(("A", "B")),
    TOPOLOGY_BC: frozenset(("B", "C")),
    TOPOLOGY_AC: frozenset(("A", "C")),
}
_PAIR_TO_TOPOLOGY = {pair: topology for topology, pair in _TOPOLOGY_TO_PAIR.items()}


def _relabel_topology(topology, old_to_new_labels):
    """Relabel a canonical topology string under an old->new label mapping."""
    old_pair = _TOPOLOGY_TO_PAIR[topology]
    new_pair = frozenset(old_to_new_labels[label] for label in old_pair)
    return _PAIR_TO_TOPOLOGY[new_pair]


def _build_topology_maps(old_to_new_labels):
    """Build topology relabel maps for a label permutation.

    Returns:
        Tuple of:
        - old_to_new_topology: maps old topology key -> relabeled topology key
        - new_to_old_topology: inverse map (relabeled topology key -> old key)
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

    Discordant1 is defined as the more frequent discordant topology. If
    discordant counts tie, keep the base ordering.

    Returns:
        Tuple ``(canonical_triplet, canonical_counts, canonical_to_original_topology, reported_dis1_topology)``
        where ``canonical_to_original_topology`` maps canonical topology keys
        (AB/BC/AC) to the source topology keys in the original observation space.
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


def _species_topology_from_newick(species_tree_newick, abc_triplet):
    """Determine species topology string for an ABC triplet."""
    if not species_tree_newick:
        return TOPOLOGY_AB

    tree = dendropy.Tree.get(
        data=species_tree_newick, schema="newick", preserve_underscores=True
    )
    return classify_triplet_topology_string(tree, abc_triplet)


def _species_tree_topology_only_newick(species_tree_newick):
    """Return species tree Newick with topology only (no branch lengths)."""
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
    """Parse rooted triplet trees into observations and metric buckets."""
    observations: list[SerializedTripletObservation] = []
    metric_buckets = (
        _build_empty_metric_buckets() if collect_summary_statistics else None
    )
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

        if (
            collect_summary_statistics
            and metric_buckets is not None
            and summary_metrics is not None
        ):
            metric_buckets[topology]["avg_tree_height"].append(
                summary_metrics["avg_tree_height"]
            )
            metric_buckets[topology]["internal_branch"].append(
                summary_metrics["internal_branch"]
            )
            metric_buckets[topology]["sister_distance"].append(
                summary_metrics["sister_distance"]
            )
        observations.append((topology, tree_height))

    return observations, metric_buckets


def _serialize_triplet_observation_rows(
    observation_rows,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
):
    """Build observations from cached parquet rows."""
    observations: list[SerializedTripletObservation] = []
    metric_buckets = (
        _build_empty_metric_buckets() if collect_summary_statistics else None
    )

    for row in observation_rows:
        topology = row.get("topology")
        if topology not in ALL_TOPOLOGIES:
            continue

        avg_tree_height = float(row.get("h_avg"))
        internal_branch = float(row.get("h_int"))
        sister_distance = float(row.get("h_sis"))
        h_a = float(row.get("h_a"))
        h_b = float(row.get("h_b"))
        h_c = float(row.get("h_c"))

        if tree_height_calculation_strategy == "AVG":
            tree_height = avg_tree_height
        elif tree_height_calculation_strategy == "A":
            tree_height = h_a
        elif tree_height_calculation_strategy == "B":
            tree_height = h_b
        elif tree_height_calculation_strategy == "C":
            tree_height = h_c
        elif tree_height_calculation_strategy == "SIS":
            tree_height = sister_distance
        elif tree_height_calculation_strategy == "INT":
            tree_height = internal_branch
        else:
            raise ValueError(
                f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
                f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
            )

        if collect_summary_statistics and metric_buckets is not None:
            metric_buckets[topology]["avg_tree_height"].append(avg_tree_height)
            metric_buckets[topology]["internal_branch"].append(internal_branch)
            metric_buckets[topology]["sister_distance"].append(sister_distance)
        observations.append((topology, tree_height))

    return observations, metric_buckets


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
    rng=None,
):
    """Run GhostParser Figure 6 pipeline from lightweight serialized observations."""
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


def run_triplet_pipeline(
    species_triplet,
    triplet_gene_trees,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    stats_backend=DEFAULT_STATS_BACKEND,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    generate_summary_stats=DEFAULT_GENERATE_SUMMARY_STATS,
    species_topology=TOPOLOGY_AB,
    species_tree_newick=None,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_options=None,
    rng=None,
):
    """Run GhostParser pipeline for one rooted species triplet.

    Args:
        species_triplet: Tuple ``(A, B, C)`` where ``A`` and ``B`` are sisters
            in species tree convention.
        triplet_gene_trees: Iterable of rooted triplet gene tree Newick strings.
        alpha_dct: Significance threshold for DCT-like Z test.
        alpha_ks: Significance threshold for THT (KS test).
        discordant_test: Discordant count test method (`chi-square` or `z-test`).
        summary_statistic: Statistic used after KS for con/dis1 distributions (`mean`, `median`, or `mode`).
        stats_backend: Statistical backend for DCT/KS (`custom` or `standard`).
        tree_height_calculation_strategy: Tree-height strategy (`AVG`, `A`, `B`, or `C`).
        species_topology: Species-tree topology string for this triplet.
        species_tree_newick: Species-tree triplet Newick string for this triplet.
        rng: Optional randomizer for tie-breaking topology ranks.

    Returns:
        ``TripletPipelineResult``.
    """
    observations = _serialize_triplet_gene_trees(
        species_triplet,
        triplet_gene_trees,
        tree_height_calculation_strategy=tree_height_calculation_strategy,
        collect_summary_statistics=generate_summary_stats,
    )
    observations, metric_buckets = observations
    base_result = _run_triplet_pipeline_from_observations(
        species_triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        species_topology=species_topology,
        species_tree_newick=species_tree_newick,
        rng=rng,
    )
    if generate_summary_stats and metric_buckets is not None:
        topology_metric_statistics = _build_topology_metric_statistics(
            species_triplet,
            species_topology,
            metric_buckets,
        )
        base_result = replace(
            base_result, topology_metric_statistics=topology_metric_statistics
        )

    if not bootstrap:
        return base_result

    options = bootstrap_options or {}
    iterations = int(options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    debug_mode = bool(options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE))
    summary_only = (
        bool(options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY))
        if debug_mode
        else False
    )
    bootstrap_rng = rng if rng is not None else random.Random()

    bootstrap_payload = _run_bootstrap_iterations(
        species_triplet,
        observations,
        iterations=iterations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        species_topology=species_topology,
        species_tree_newick=None,
        debug_mode=debug_mode,
        summary_only=summary_only,
        rng=bootstrap_rng,
    )

    bootstrap_gene_tree_heights = None
    if debug_mode:
        raw_heights = [tree_height for _, tree_height in observations]
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


def parse_triplet_gene_trees_file(filepath):
    """Parse triplet gene tree file produced by ``tree_parser``.

    Required section structure per triplet:
    1. ``A,B,C<TAB>count<TAB>species_triplet_newick<TAB>[A=...,B=...,C=...]<TAB>AB:x/concordant,BC:y/discordant1|discordant2,AC:z/discordant2|discordant1``
    2. one blank line
    3. zero or more Newick gene-tree lines (each ending in ``;``)

    Sections are separated by a line of exactly 60 ``=`` characters.

    Returns:
        dict mapping triplet tuple ->
        ``{"count": int | None, "species_tree": str | None, "gene_trees": list[str]}``
    """
    section_separator = "=" * 60

    def _parse_section(section_lines):
        trimmed = list(section_lines)
        while trimmed and not trimmed[0].strip():
            trimmed.pop(0)
        while trimmed and not trimmed[-1].strip():
            trimmed.pop()

        if not trimmed:
            return None

        header = trimmed[0]
        parts = header.split("\t")
        if len(parts) != 5:
            raise ValueError(
                f"Invalid triplet header format (expected 5 tab-separated fields): {header}"
            )

        triplet_text = parts[0].strip()
        taxa = tuple(part.strip() for part in triplet_text.split(",") if part.strip())
        if len(taxa) != 3:
            raise ValueError(f"Invalid triplet header: {header}")

        if not parts[1].strip():
            raise ValueError(f"Invalid triplet count in header: {header}")
        try:
            count = int(parts[1].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid triplet count in header: {header}") from exc

        if not parts[2].strip():
            raise ValueError(f"Invalid species tree in header: {header}")
        species_tree = parts[2].strip()

        label_text = parts[3].strip()
        if not (label_text.startswith("[") and label_text.endswith("]")):
            raise ValueError(f"Invalid ABC label mapping in header: {header}")
        label_body = label_text[1:-1]
        label_pairs = [
            segment.strip() for segment in label_body.split(",") if segment.strip()
        ]
        if len(label_pairs) != 3:
            raise ValueError(f"Invalid ABC label mapping in header: {header}")

        mapped_values = []
        expected_keys = ("A", "B", "C")
        for expected_key, label_pair in zip(expected_keys, label_pairs):
            if "=" not in label_pair:
                raise ValueError(f"Invalid ABC label mapping in header: {header}")
            key, value = (piece.strip() for piece in label_pair.split("=", 1))
            if key != expected_key or not value:
                raise ValueError(f"Invalid ABC label mapping in header: {header}")
            mapped_values.append(value)
        mapped_triplet = tuple(mapped_values)
        if mapped_triplet != taxa:
            raise ValueError(
                f"Triplet/header label mapping mismatch for {','.join(taxa)}: "
                f"A,B,C mapping resolves to {','.join(mapped_triplet)}"
            )

        summary_text = parts[4].strip()
        summary_items = [
            segment.strip() for segment in summary_text.split(",") if segment.strip()
        ]
        if len(summary_items) != 3:
            raise ValueError(f"Invalid topology summary in header: {header}")

        parsed_summary = {}
        for item in summary_items:
            if ":" not in item:
                raise ValueError(f"Invalid topology summary in header: {header}")
            topology_key, payload = (piece.strip() for piece in item.split(":", 1))
            if "/" not in payload:
                raise ValueError(f"Invalid topology summary in header: {header}")
            count_text, role = (piece.strip() for piece in payload.split("/", 1))
            if topology_key in parsed_summary:
                raise ValueError(f"Invalid topology summary in header: {header}")
            try:
                parsed_count = int(count_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid topology summary in header: {header}"
                ) from exc
            parsed_summary[topology_key] = (parsed_count, role)

        if set(parsed_summary.keys()) != {"AB", "BC", "AC"}:
            raise ValueError(f"Invalid topology summary in header: {header}")

        n_ab, ab_role = parsed_summary["AB"]
        n_bc, bc_role = parsed_summary["BC"]
        n_ac, ac_role = parsed_summary["AC"]
        if ab_role != "concordant":
            raise ValueError(f"Invalid topology summary in header: {header}")
        if bc_role not in {"discordant1", "discordant2"}:
            raise ValueError(f"Invalid topology summary in header: {header}")
        if ac_role not in {"discordant1", "discordant2"}:
            raise ValueError(f"Invalid topology summary in header: {header}")
        if bc_role == ac_role:
            raise ValueError(f"Invalid discordant role assignment in header: {header}")

        summary_total = n_ab + n_bc + n_ac
        if summary_total != count:
            raise ValueError(
                f"Triplet count/header mismatch for {','.join(taxa)}: "
                f"header count={count}, topology summary total={summary_total}"
            )

        if len(trimmed) > 1 and trimmed[1].strip():
            raise ValueError(
                "Invalid triplet section format: expected blank line after header"
            )

        gene_trees = []
        for line in trimmed[2:]:
            if not line.strip():
                continue
            if not line.endswith(";"):
                raise ValueError(
                    f"Invalid gene-tree line (expected Newick ending with ';'): {line}"
                )
            gene_trees.append(line)

        if count != len(gene_trees):
            raise ValueError(
                f"Triplet count/header mismatch for {','.join(taxa)}: header count={count}, parsed trees={len(gene_trees)}"
            )

        dis1_topology = TOPOLOGY_BC if bc_role == "discordant1" else TOPOLOGY_AC

        return taxa, {
            "count": count,
            "species_tree": species_tree,
            "gene_trees": gene_trees,
            "label_map": {"A": taxa[0], "B": taxa[1], "C": taxa[2]},
            "header_topology_counts": {
                TOPOLOGY_AB: n_ab,
                TOPOLOGY_BC: n_bc,
                TOPOLOGY_AC: n_ac,
            },
            "header_dis1_topology": dis1_topology,
        }

    triplet_map = {}
    current_section = []

    with open(filepath, "r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n").rstrip("\r")
            if line == section_separator:
                parsed = _parse_section(current_section)
                if parsed is not None:
                    taxa, payload = parsed
                    if taxa in triplet_map:
                        raise ValueError(
                            f"Duplicate triplet header encountered: {','.join(taxa)}"
                        )
                    triplet_map[taxa] = payload
                current_section = []
                continue
            current_section.append(line)

    parsed = _parse_section(current_section)
    if parsed is not None:
        taxa, payload = parsed
        if taxa in triplet_map:
            raise ValueError(f"Duplicate triplet header encountered: {','.join(taxa)}")
        triplet_map[taxa] = payload

    return triplet_map


def _require_pyarrow():
    """Import and return pyarrow parquet module."""
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError("Parquet input requires pyarrow to be installed") from exc
    return pq


def parse_triplet_gene_trees_parquet(dataset_path):
    """Parse partitioned parquet triplet dataset produced by ``tree_parser``."""
    pq = _require_pyarrow()

    dataset_root = Path(dataset_path)
    triplet_dir = dataset_root / "triplets"
    observations_dir = dataset_root / "observations"
    if not triplet_dir.exists() or not observations_dir.exists():
        raise ValueError(
            "Invalid parquet triplet dataset: expected 'triplets/' and 'observations/' directories"
        )

    triplet_files = sorted(triplet_dir.glob("*.parquet"))
    if not triplet_files:
        return {}

    triplet_map = {}
    triplet_id_to_taxa = {}

    for file_path in triplet_files:
        parquet_file = pq.ParquetFile(file_path)
        for batch in parquet_file.iter_batches():
            rows = batch.to_pylist()
            for row in rows:
                taxa = (row["A"], row["B"], row["C"])
                triplet_id = row["triplet_id"]
                if taxa in triplet_map:
                    raise ValueError(
                        f"Duplicate triplet header encountered: {','.join(taxa)}"
                    )

                bc_role = row["bc_role"]
                ac_role = row["ac_role"]
                if bc_role == ac_role:
                    raise ValueError(
                        f"Invalid discordant role assignment in parquet row for {','.join(taxa)}"
                    )

                dis1_topology = TOPOLOGY_BC if bc_role == "discordant1" else TOPOLOGY_AC
                triplet_map[taxa] = {
                    "count": int(row["count"]),
                    "species_tree": row.get("species_tree"),
                    "gene_trees": [],
                    "observation_rows": [],
                    "label_map": {"A": taxa[0], "B": taxa[1], "C": taxa[2]},
                    "header_topology_counts": {
                        TOPOLOGY_AB: int(row["n_ab"]),
                        TOPOLOGY_BC: int(row["n_bc"]),
                        TOPOLOGY_AC: int(row["n_ac"]),
                    },
                    "header_dis1_topology": dis1_topology,
                }
                triplet_id_to_taxa[triplet_id] = taxa

    observation_files = sorted(observations_dir.rglob("*.parquet"))
    for file_path in observation_files:
        parquet_file = pq.ParquetFile(file_path)
        for batch in parquet_file.iter_batches(
            columns=["triplet_id", "topology", "h_avg", "h_int", "h_sis", "h_a", "h_b", "h_c"]
        ):
            rows = batch.to_pylist()
            for row in rows:
                triplet_id = row["triplet_id"]
                taxa = triplet_id_to_taxa.get(triplet_id)
                if taxa is None:
                    raise ValueError(
                        f"Observation references unknown triplet_id: {triplet_id}"
                    )
                triplet_map[taxa]["observation_rows"].append(
                    {
                        "topology": row["topology"],
                        "h_avg": row["h_avg"],
                        "h_int": row["h_int"],
                        "h_sis": row["h_sis"],
                        "h_a": row["h_a"],
                        "h_b": row["h_b"],
                        "h_c": row["h_c"],
                    }
                )

    for taxa, entry in triplet_map.items():
        if entry["count"] != len(entry["observation_rows"]):
            raise ValueError(
                f"Triplet count/header mismatch for {','.join(taxa)}: "
                f"header count={entry['count']}, parsed rows={len(entry['observation_rows'])}"
            )

    return triplet_map


def _resolve_input_format(filepath, input_format=DEFAULT_INPUT_FORMAT):
    """Resolve triplet input format from explicit value or path suffix."""
    if input_format not in INPUT_FORMAT_CHOICES:
        raise ValueError(
            f"Unsupported input format: {input_format}. "
            f"Choose one of: {', '.join(INPUT_FORMAT_CHOICES)}"
        )
    if input_format != "auto":
        return input_format

    path = Path(filepath)
    return "parquet" if path.suffix.lower() == ".parquet" else "txt"


def _get_mp_context(prefer_fork=None):
    """Get a multiprocessing context for worker pools.

    Args:
        prefer_fork: If True, prefer ``fork`` when available. If None, defaults
            to True on POSIX platforms.
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


def _resolve_processes(processes):
    """Resolve process count where 0 means all available CPU cores."""
    if processes == 0:
        return cpu_count()
    return processes


def _init_triplet_analysis_worker(entry_map=None):
    """Initialize worker-local triplet entry map.

    Under ``fork`` this map is inherited via copy-on-write and ``entry_map``
    can be None. Under other start methods, we pass the map once via initargs
    to avoid sending each large entry through every task payload.
    """
    global _TRIPLET_ENTRY_MAP
    if entry_map is not None:
        _TRIPLET_ENTRY_MAP = entry_map


def _analyze_triplet_entry(args):
    """Analyze one triplet map entry in a worker process."""
    (
        triplet,
        alpha_dct,
        alpha_ks,
        discordant_test,
        summary_statistic,
        stats_backend,
        tree_height_calculation_strategy,
        generate_summary_stats,
        bootstrap,
        bootstrap_options,
        triplet_seed,
    ) = args
    entry = _TRIPLET_ENTRY_MAP.get(triplet)
    if entry is None:
        raise ValueError(f"Missing worker triplet entry for: {','.join(triplet)}")

    start_cpu = time.process_time()
    result = analyze_triplet_entry(
        triplet,
        entry,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        tree_height_calculation_strategy=tree_height_calculation_strategy,
        generate_summary_stats=generate_summary_stats,
        bootstrap=bootstrap,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
    )
    worker_cpu_seconds = time.process_time() - start_cpu
    return result, worker_cpu_seconds


def analyze_triplet_entry(
    triplet,
    entry,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    stats_backend=DEFAULT_STATS_BACKEND,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    generate_summary_stats=DEFAULT_GENERATE_SUMMARY_STATS,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_options=None,
    triplet_seed=None,
):
    """Analyze one triplet entry payload and return a pipeline result."""
    species_tree_topology = _species_tree_topology_only_newick(
        entry.get("species_tree")
    )
    species_topology = TOPOLOGY_AB
    if entry.get("observation_rows"):
        observations, metric_buckets = _serialize_triplet_observation_rows(
            entry["observation_rows"],
            tree_height_calculation_strategy=tree_height_calculation_strategy,
            collect_summary_statistics=generate_summary_stats,
        )
    else:
        observations, metric_buckets = _serialize_triplet_gene_trees(
            triplet,
            entry["gene_trees"],
            tree_height_calculation_strategy=tree_height_calculation_strategy,
            collect_summary_statistics=generate_summary_stats,
        )

    base_result = _run_triplet_pipeline_from_observations(
        triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        species_topology=species_topology,
        species_tree_newick=species_tree_topology,
    )

    if generate_summary_stats and metric_buckets is not None:
        topology_metric_statistics = _build_topology_metric_statistics(
            triplet,
            species_topology,
            metric_buckets,
        )
        base_result = replace(
            base_result, topology_metric_statistics=topology_metric_statistics
        )

    options = dict(bootstrap_options or {})
    iterations = int(options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    debug_mode = bool(options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE))
    summary_only = (
        bool(options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY))
        if debug_mode
        else False
    )
    rng = _build_triplet_rng(triplet_seed, triplet)

    bootstrap_payload = _run_bootstrap_iterations(
        triplet,
        observations,
        iterations=iterations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        summary_statistic=summary_statistic,
        stats_backend=stats_backend,
        species_topology=species_topology,
        species_tree_newick=None,
        debug_mode=debug_mode,
        summary_only=summary_only,
        rng=rng,
    )

    bootstrap_gene_tree_heights = None
    if debug_mode:
        raw_heights = [tree_height for _, tree_height in observations]
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


def analyze_triplet_gene_tree_file(
    filepath,
    input_format="auto",
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    stats_backend=DEFAULT_STATS_BACKEND,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    generate_summary_stats=DEFAULT_GENERATE_SUMMARY_STATS,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_options=None,
    rng=None,
    use_multiprocessing=True,
    processes=None,
    return_worker_cpu=False,
):
    """Analyze all triplets from a triplet-gene-trees text or parquet input."""
    if discordant_test not in DISCORDANT_TEST_CHOICES:
        raise ValueError(
            f"Unsupported discordant test method: {discordant_test}. "
            f"Choose one of: {', '.join(DISCORDANT_TEST_CHOICES)}"
        )

    if summary_statistic not in SUMMARY_STATISTIC_CHOICES:
        raise ValueError(
            f"Unsupported summary statistic: {summary_statistic}. "
            f"Choose one of: {', '.join(SUMMARY_STATISTIC_CHOICES)}"
        )

    if stats_backend not in STATS_BACKEND_CHOICES:
        raise ValueError(
            f"Unsupported stats backend: {stats_backend}. "
            f"Choose one of: {', '.join(STATS_BACKEND_CHOICES)}"
        )

    if tree_height_calculation_strategy not in TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES:
        raise ValueError(
            f"Unsupported tree height calculation strategy: {tree_height_calculation_strategy}. "
            f"Choose one of: {', '.join(TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES)}"
        )

    if p_value_correction not in P_VALUE_CORRECTION_CHOICES:
        raise ValueError(
            f"Unsupported p-value correction method: {p_value_correction}. "
            f"Choose one of: {', '.join(P_VALUE_CORRECTION_CHOICES)}"
        )

    options = dict(bootstrap_options or {})
    if "iterations" not in options:
        options["iterations"] = DEFAULT_BOOTSTRAP_ITERATIONS
    if "summary_only" not in options:
        options["summary_only"] = DEFAULT_BOOTSTRAP_SUMMARY_ONLY
    if "debug_mode" not in options:
        options["debug_mode"] = DEFAULT_BOOTSTRAP_DEBUG_MODE
    if "seed" not in options:
        options["seed"] = None

    resolved_input_format = _resolve_input_format(filepath, input_format=input_format)
    if resolved_input_format == "parquet":
        triplet_map = parse_triplet_gene_trees_parquet(filepath)
    else:
        triplet_map = parse_triplet_gene_trees_file(filepath)
    items = list(triplet_map.items())
    if not items:
        return []

    resolved_processes = _resolve_processes(processes)
    worker_count = resolved_processes or cpu_count()
    worker_count = max(1, min(worker_count, len(items)))

    if use_multiprocessing and worker_count > 1 and rng is None:
        ctx = _get_mp_context()
        init_entry_map = triplet_map
        if hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork":
            global _TRIPLET_ENTRY_MAP
            _TRIPLET_ENTRY_MAP = triplet_map
            init_entry_map = None

        args = [
            (
                triplet,
                alpha_dct,
                alpha_ks,
                discordant_test,
                summary_statistic,
                stats_backend,
                tree_height_calculation_strategy,
                generate_summary_stats,
                bootstrap,
                options,
                options.get("seed"),
            )
            for triplet, _ in items
        ]
        chunksize = max(1, len(args) // (worker_count * 4))
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_triplet_analysis_worker,
            initargs=(init_entry_map,),
        ) as pool:
            payloads = list(
                pool.imap(_analyze_triplet_entry, args, chunksize=chunksize)
            )
        results = [payload[0] for payload in payloads]
        worker_cpu_seconds = sum(payload[1] for payload in payloads)
        corrected = _apply_triplet_result_p_value_correction(
            results,
            alpha_dct=alpha_dct,
            alpha_ks=alpha_ks,
            method=p_value_correction,
            stats_backend=stats_backend,
        )
        if return_worker_cpu:
            return corrected, worker_cpu_seconds
        return corrected

    results = []
    for triplet, entry in items:
        results.append(
            analyze_triplet_entry(
                triplet,
                entry,
                alpha_dct=alpha_dct,
                alpha_ks=alpha_ks,
                discordant_test=discordant_test,
                summary_statistic=summary_statistic,
                stats_backend=stats_backend,
                tree_height_calculation_strategy=tree_height_calculation_strategy,
                generate_summary_stats=generate_summary_stats,
                bootstrap=bootstrap,
                bootstrap_options=options,
                triplet_seed=options.get("seed"),
            )
        )
    corrected = _apply_triplet_result_p_value_correction(
        results,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        method=p_value_correction,
        stats_backend=stats_backend,
    )
    if return_worker_cpu:
        return corrected, 0.0
    return corrected


def collect_triplet_statistics(results):
    """Return triplet results as list of dictionaries for downstream use."""
    return [result.to_dict() for result in results]


def write_pipeline_statistics_json(results, output_filepath):
    """Write all triplet statistics as JSON list for downstream analysis."""
    data = collect_triplet_statistics(results)
    with open(output_filepath, "w") as out_f:
        json.dump(data, out_f, indent=2)


def write_pipeline_results(
    results,
    output_filepath,
    dct_method=DEFAULT_DISCORDANT_TEST,
    summary_statistic=DEFAULT_SUMMARY_STATISTIC,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_debug_mode=DEFAULT_BOOTSTRAP_DEBUG_MODE,
):
    """Write pipeline results to TSV file."""
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
    """Write per-triplet topology/metric summary statistics to TSV."""
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


def main():
    """CLI entry point for triplet processing pipeline."""
    parser = _build_argument_parser()
    parsed_args = parser.parse_args()

    try:
        args = _resolve_runtime_args(parsed_args)
    except ValueError as exc:
        print(f"Error: {exc}")
        return

    input_path = Path(args.input)
    output_path = (
        Path(args.output)
        if args.output
        else input_path.parent / "triplet_introgression_results.tsv"
    )

    results = analyze_triplet_gene_tree_file(
        str(input_path),
        input_format=args.input_format,
        alpha_dct=args.alpha_dct,
        alpha_ks=args.alpha_ks,
        discordant_test=args.discordant_test,
        summary_statistic=args.summary_statistic,
        stats_backend=args.stats_backend,
        tree_height_calculation_strategy=args.tree_height_calculation_strategy,
        p_value_correction=args.p_value_correction,
        generate_summary_stats=args.generate_summary_stats,
        bootstrap=args.bootstrap,
        bootstrap_options=args.bootstrap_options,
        use_multiprocessing=not args.no_multiprocessing,
        processes=args.processes,
    )
    write_pipeline_results(
        results,
        str(output_path),
        dct_method=args.discordant_test,
        summary_statistic=args.summary_statistic,
        p_value_correction=args.p_value_correction,
        bootstrap=args.bootstrap,
        bootstrap_debug_mode=args.bootstrap_options["debug_mode"],
    )

    summary_output = output_path.parent / "summary_statistics.tsv"
    if args.generate_summary_stats:
        write_summary_statistics_tsv(
            results,
            str(summary_output),
            bootstrap=args.bootstrap,
        )

    stats_output = (
        Path(args.stats_output)
        if args.stats_output
        else output_path.with_suffix(".json")
    )
    write_pipeline_statistics_json(results, str(stats_output))

    print(f"Processed {len(results)} triplets")
    print(f"Results written to: {output_path}")
    if args.generate_summary_stats:
        print(f"Summary statistics written to: {summary_output}")
    else:
        print("Summary statistics written to: skipped")
    print(f"Statistics JSON written to: {stats_output}")


def _build_argument_parser():
    """Build triplet_processor CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Run GhostParser triplet processing pipeline (Fig. 6)."
    )
    parser.add_argument(
        "--input-path",
        required=True,
        help="Path to unique_triplets_gene_trees.txt or parquet dataset",
    )
    parser.add_argument(
        "--input-format",
        choices=INPUT_FORMAT_CHOICES,
        default=None,
        help=f"Input format (default: {DEFAULT_INPUT_FORMAT}; auto infers from input path suffix)",
    )
    parser.add_argument(
        "--output-path", default=None, help="Output TSV path (default: alongside input)"
    )
    parser.add_argument(
        "--stats-output",
        dest="stats_output",
        default=None,
        help="Optional JSON output path for full per-triplet statistics",
    )
    parser.add_argument(
        "--alpha-dct",
        type=float,
        default=None,
        help=f"DCT significance threshold (default: {DEFAULT_ALPHA_DCT})",
    )
    parser.add_argument(
        "--alpha-ks",
        type=float,
        default=None,
        help=f"KS significance threshold (default: {DEFAULT_ALPHA_KS})",
    )
    parser.add_argument(
        "--discordant-test",
        choices=DISCORDANT_TEST_CHOICES,
        default=None,
        help=f"Discordant count test to use (default: {DEFAULT_DISCORDANT_TEST})",
    )
    parser.add_argument(
        "--summary-statistic",
        choices=SUMMARY_STATISTIC_CHOICES,
        default=None,
        help=f"Statistic used for con/dis1 distributions after KS test (default: {DEFAULT_SUMMARY_STATISTIC})",
    )
    parser.add_argument(
        "--stats-backend",
        choices=STATS_BACKEND_CHOICES,
        default=None,
        help=f"Statistical backend for DCT/KS calculations (default: {DEFAULT_STATS_BACKEND})",
    )
    parser.add_argument(
        "--tree-height-calculation-strategy",
        choices=TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
        default=None,
        help=(
            "Tree-height strategy: AVG uses mean root-to-tip distance, "
            "A/B/C use the selected taxon's root-to-tip distance, "
            "SIS uses sister-taxon distance, and INT uses internal branch length "
            f"(default: {DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY})"
        ),
    )
    parser.add_argument(
        "--p-value-correction",
        choices=P_VALUE_CORRECTION_CHOICES,
        default=None,
        help=f"Multiple-testing correction for triplet p-values (default: {DEFAULT_P_VALUE_CORRECTION})",
    )
    parser.add_argument(
        "--no-bootstrap",
        dest="bootstrap",
        action="store_false",
        help="Disable bootstrap sampling-with-replacement during triplet analysis",
    )
    parser.set_defaults(bootstrap=None)
    parser.add_argument(
        "--bootstrap-iterations",
        type=int,
        default=None,
        help=f"Number of bootstrap iterations per triplet (default: {DEFAULT_BOOTSTRAP_ITERATIONS})",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=None,
        help="Optional bootstrap random seed for reproducibility",
    )
    parser.add_argument(
        "--bootstrap-debug-mode",
        dest="bootstrap_debug_mode",
        action="store_true",
        help="Enable bootstrap debug output columns",
    )
    parser.add_argument(
        "--bootstrap-summary-only",
        dest="bootstrap_summary_only",
        action="store_true",
        help="When bootstrap debug mode is enabled, emit compact summaries instead of full per-iteration lists",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes for triplet analysis (0 = all cores)",
    )
    parser.add_argument(
        "--generate-summary-stats",
        dest="generate_summary_stats",
        action="store_true",
        default=None,
        help="Generate summary_statistics.tsv output (default: False)",
    )
    parser.add_argument(
        "--no-multiprocessing",
        action="store_true",
        help="Disable multiprocessing for triplet analysis",
    )
    return parser


def _resolve_runtime_args(args):
    """Resolve runtime arguments from CLI mode."""
    input_path = Path(args.input_path).expanduser().resolve()

    output = None
    if args.output_path:
        output = str(Path(args.output_path).expanduser().resolve())

    stats_output = None
    if args.stats_output:
        stats_output = str(Path(args.stats_output).expanduser().resolve())

    processes = args.processes if args.processes is not None else 0
    if processes < 0:
        raise ValueError("CLI argument --processes must be an integer >= 0")

    bootstrap_iterations = (
        args.bootstrap_iterations
        if args.bootstrap_iterations is not None
        else DEFAULT_BOOTSTRAP_ITERATIONS
    )
    if bootstrap_iterations < 1:
        raise ValueError("CLI argument --bootstrap-iterations must be an integer >= 1")

    return argparse.Namespace(
        input=str(input_path),
        input_format=args.input_format or DEFAULT_INPUT_FORMAT,
        output=output,
        stats_output=stats_output,
        alpha_dct=args.alpha_dct if args.alpha_dct is not None else DEFAULT_ALPHA_DCT,
        alpha_ks=args.alpha_ks if args.alpha_ks is not None else DEFAULT_ALPHA_KS,
        discordant_test=args.discordant_test or DEFAULT_DISCORDANT_TEST,
        summary_statistic=args.summary_statistic or DEFAULT_SUMMARY_STATISTIC,
        stats_backend=args.stats_backend or DEFAULT_STATS_BACKEND,
        tree_height_calculation_strategy=args.tree_height_calculation_strategy
        or DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
        p_value_correction=args.p_value_correction or DEFAULT_P_VALUE_CORRECTION,
        bootstrap=DEFAULT_BOOTSTRAP if args.bootstrap is None else bool(args.bootstrap),
        bootstrap_options={
            "iterations": bootstrap_iterations,
            "seed": args.bootstrap_seed,
            "debug_mode": bool(args.bootstrap_debug_mode),
            "summary_only": bool(args.bootstrap_summary_only),
        },
        processes=processes,
        generate_summary_stats=bool(args.generate_summary_stats),
        no_multiprocessing=bool(args.no_multiprocessing),
    )


if __name__ == "__main__":
    main()
