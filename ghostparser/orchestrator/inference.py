"""Per-triplet GhostParser inference for the orchestrator.

Topology classification, the three decision gates, bootstrap aggregation,
summary statistics, run-wide p-value correction, and TSV writing. All statistics
use the scipy/statsmodels backend.
"""

import hashlib
import json
import math
from dataclasses import dataclass, replace

import dendropy
import numpy as np
from scipy import stats
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
    DEFAULT_ALPHA_PERM,
    DEFAULT_BOOTSTRAP,
    DEFAULT_BOOTSTRAP_DEBUG_MODE,
    DEFAULT_BOOTSTRAP_ITERATIONS,
    DEFAULT_BOOTSTRAP_SUMMARY_ONLY,
    DEFAULT_DISCORDANT_TEST,
    DEFAULT_P_VALUE_CORRECTION,
    DEFAULT_PERMUTATION_MAX_RESAMPLES,
    DEFAULT_PERMUTATION_MIN_RESAMPLES,
    DEFAULT_PIPELINE_MODE,
    DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    DISCORDANT_TEST_CHOICES,
    P_VALUE_CORRECTION_CHOICES,
    PIPELINE_MODE_EFFICIENT,
    TREE_HEIGHT_CALCULATION_STRATEGY_CHOICES,
)
from .correction import (
    adjust_p_value_inline as _adjust_p_value_inline,
    adjust_p_values as _adjust_p_values,
    is_inline_correction,
)
from .shape import SHAPE_FIELD_NAMES, describe_shape
from .permutation import (
    DECISION_GREATER,
    DECISION_INCONCLUSIVE,
    DECISION_LESS,
    bootstrap_resample_budget,
    run_studentized_permutation_test,
    studentized_mean_diff,
)

Classification = str
SerializedTripletObservation = tuple[str, float, dict | None]

DISCORDANT1_TOPOLOGY_CHOICES = ("BC", "AC")

# Height groups the optional shape diagnostics describe, in column order. They
# are written to the results TSV only: the summary-statistics TSV is a feature
# matrix, and these are undefined for small groups, so including them there
# would punch holes in it.
SHAPE_GROUP_LABELS = ("con", "dis1", "dis2")

# Which of the three tests produced a triplet's classification, named after the
# test itself. ``DECISION_GATE_PERMUTATION`` is the only value for which the
# permutation columns took part in the call; on the other two nothing consulted
# the direction test, and under the efficient pipeline mode it did not run.
DECISION_GATE_DCT = "DCT"
DECISION_GATE_THT = "THT"
DECISION_GATE_PERMUTATION = "PERM"
DECISION_GATE_CHOICES = (
    DECISION_GATE_DCT,
    DECISION_GATE_THT,
    DECISION_GATE_PERMUTATION,
)

# Integer codes for a bootstrap iteration's direction outcome, so the deferred
# per-iteration record stays a compact numeric array.
_DIRECTION_SKIPPED = 0
_DIRECTION_GREATER = 1
_DIRECTION_LESS = 2
_DIRECTION_CODES = {DECISION_GREATER: _DIRECTION_GREATER, DECISION_LESS: _DIRECTION_LESS}
_DIRECTION_LABELS = {
    _DIRECTION_GREATER: DECISION_GREATER,
    _DIRECTION_LESS: DECISION_LESS,
    _DIRECTION_SKIPPED: DECISION_INCONCLUSIVE,
}

# Recorded in ``perm_note`` when the efficient pipeline mode skipped the
# direction test because the cascade had already been settled by an earlier
# gate. It marks an empty ``perm_*`` block as deliberate rather than as a test
# that ran and failed, which is what the guard slugs in ``permutation`` mean.
PERM_NOTE_NOT_CONSULTED = "direction_test_not_consulted"

_BOOTSTRAP_CLASSES = [
    "ghost_introgression",
    "inflow_introgression",
    "outflow_introgression",
    "no_introgression",
    "ambiguous",
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
# Which canonical topology each summary label refers to is resolved per triplet
# from the observed counts, not fixed here: `discordant1` means the more
# frequent discordant topology, matching the `dis1_topology` column and the
# group the DCT, KS, and permutation tests all operate on.


@dataclass(frozen=True)
class DeferredBootstrapRecord:
    """One triplet's raw per-iteration bootstrap p-values, awaiting correction.

    Rank-based corrections need every triplet's p-value for the same iteration
    index before any of them can be adjusted, so iterations that cannot be
    settled during the stream park their raw values here instead.

    Attributes:
        dct_p_values: Raw DCT p-value per iteration.
        ks_p_values: Raw KS p-value per iteration.
        directions: Per-iteration direction code, ``_DIRECTION_SKIPPED`` where
            the permutation test was short-circuited or found no direction.
    """

    dct_p_values: np.ndarray
    ks_p_values: np.ndarray
    directions: np.ndarray


@dataclass(frozen=True)
class _BootstrapPolicy:
    """How bootstrap iterations apply the run's correction method.

    Attributes:
        correction: The configured p-value correction method.
        family_size: Triplet count, the size of each per-iteration family.
        inline: Whether the method can be applied to one p-value from the family
            size alone, which lets an iteration tighten its short-circuit test.
    """

    correction: str
    family_size: int
    inline: bool


@dataclass(frozen=True)
class TripletPipelineResult:
    """Result of the GhostParser orchestrator for one rooted species triplet.

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
    classification: Classification
    decision_gate: str | None = None
    perm_decision: str | None = None
    perm_statistic: float | None = None
    perm_p_greater: float | None = None
    perm_p_less: float | None = None
    perm_null_skew: float | None = None
    perm_p_greater_corrected: float | None = None
    perm_p_less_corrected: float | None = None
    perm_p_tost: float | None = None
    bootstrap_stat_ci_low: float | None = None
    bootstrap_stat_ci_high: float | None = None
    perm_n_resamples: int | None = None
    perm_converged: bool | None = None
    perm_n_resamples_skew: int | None = None
    perm_note: str | None = None
    shape_statistics: dict[str, float | None] | None = None
    analyzed_trees: int = 0
    topology_metric_statistics: dict[str, float | None] | None = None
    bootstrap_value: float | None = None
    all_bootstrap: dict | None = None
    bootstrap_deferred: "DeferredBootstrapRecord | None" = None
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
            "classification": self.classification,
            "decision_gate": self.decision_gate,
            "perm_decision": self.perm_decision,
            "perm_statistic": self.perm_statistic,
            "perm_p_greater": self.perm_p_greater,
            "perm_p_less": self.perm_p_less,
            "perm_null_skew": self.perm_null_skew,
            "perm_p_greater_corrected": self.perm_p_greater_corrected,
            "perm_p_less_corrected": self.perm_p_less_corrected,
            "perm_p_tost": self.perm_p_tost,
            "bootstrap_stat_ci_low": self.bootstrap_stat_ci_low,
            "bootstrap_stat_ci_high": self.bootstrap_stat_ci_high,
            "perm_n_resamples": self.perm_n_resamples,
            "perm_converged": self.perm_converged,
            "perm_n_resamples_skew": self.perm_n_resamples_skew,
            "perm_note": self.perm_note,
            "shape_statistics": self.shape_statistics,
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

    Reached only from the DendroPy-tree reference path
    (:func:`observation_from_subtree`, :func:`_serialize_triplet_gene_trees`);
    a run computes these values from a cached geometry.

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


def shape_column_names(group_labels=SHAPE_GROUP_LABELS):
    """Return the shape-diagnostic column names in output order.

    Args:
        group_labels: Per-group prefixes. Only the results TSV carries these
            columns, so the short forms are the only ones in use.

    Returns:
        A list of ``<group>_<field>`` column names.
    """
    return [
        f"{label}_{field}" for label in group_labels for field in SHAPE_FIELD_NAMES
    ]


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
    species_triplet, species_topology, metric_buckets, dis1_topology, dis2_topology
):
    """Build topology/metric summary statistics for one triplet.

    The ``discordant1`` columns describe whichever discordant topology is more
    frequent, so they name the same set of gene trees as the ``dis1_topology``
    column and as the groups the DCT, KS, and permutation tests use. The roles
    are passed in rather than recomputed here so they cannot drift from the ones
    the decision logic resolved.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        species_topology: The concordant (species-tree) topology.
        metric_buckets: Per-topology metric buckets from
            :func:`_build_empty_metric_buckets`.
        dis1_topology: Canonical topology of the more frequent discordant.
        dis2_topology: Canonical topology of the less frequent discordant.

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

    label_to_canonical = {
        "concordant": species_topology,
        "discordant1": dis1_topology,
        "discordant2": dis2_topology,
    }

    statistics = {}
    for topology_label in SUMMARY_TOPOLOGY_LABELS:
        canonical_topology = label_to_canonical[topology_label]
        original_topology = canonical_to_original_topology[canonical_topology]
        for metric_label in SUMMARY_METRIC_LABELS:
            values = metric_buckets[original_topology][metric_label]
            for statistic_name in SUMMARY_STATISTICS:
                key = f"{topology_label}_{metric_label}_{statistic_name}"
                statistics[key] = _compute_summary_statistic(values, statistic_name)

    return statistics


def _classify_introgression(dct_significant, ks_significant, direction):
    """Apply the GhostParser decision cascade to produce a classification.

    Each branch returns the classification and the gate that produced it in one
    pass, so the two cannot drift apart. What each gate means biologically is
    described under "The statistical tests" in the orchestrator guide.

    Args:
        dct_significant: Whether the discordant count test is significant.
        ks_significant: Whether the tree-height (KS) test is significant.
            ``None`` (no test ran) takes the same branch as non-significant.
        direction: The permutation decision; only ``greater`` and ``less`` name
            a direction.

    Returns:
        A ``(classification, decision_gate)`` tuple.
    """
    if not dct_significant:
        return "no_introgression", DECISION_GATE_DCT
    if not ks_significant:
        return "inflow_introgression", DECISION_GATE_THT
    if direction == DECISION_GREATER:
        return "outflow_introgression", DECISION_GATE_PERMUTATION
    if direction == DECISION_LESS:
        return "ghost_introgression", DECISION_GATE_PERMUTATION
    # The heights differ in distribution but not detectably in mean, so the
    # difference is in shape rather than location and no direction is defensible.
    return "ambiguous", DECISION_GATE_PERMUTATION


def _permutation_result_fields(permutation_result, note=None):
    """Flatten a permutation result into ``TripletPipelineResult`` field values.

    Args:
        permutation_result: A ``PermutationTestResult``, or ``None`` when the
            permutation test was disabled or skipped.
        note: Slug recording why no test ran, reported in ``perm_note`` so an
            empty block is distinguishable from one a guard emptied. Ignored
            when ``permutation_result`` is present, which carries its own note.

    Returns:
        A dict of ``perm_*`` keyword arguments, all ``None`` when no test ran.
    """
    if permutation_result is None:
        return {
            "perm_statistic": None,
            "perm_p_greater": None,
            "perm_p_less": None,
            "perm_null_skew": None,
            "perm_p_greater_corrected": None,
            "perm_p_less_corrected": None,
            "perm_p_tost": None,
            "perm_n_resamples": None,
            "perm_converged": None,
            "perm_n_resamples_skew": None,
            "perm_note": note,
        }

    return {
        "perm_statistic": permutation_result.statistic,
        "perm_p_greater": permutation_result.p_greater,
        "perm_p_less": permutation_result.p_less,
        "perm_null_skew": permutation_result.null_skew,
        "perm_p_greater_corrected": permutation_result.p_greater_corrected,
        "perm_p_less_corrected": permutation_result.p_less_corrected,
        "perm_p_tost": permutation_result.p_tost,
        "perm_n_resamples": permutation_result.n_resamples,
        "perm_converged": permutation_result.converged,
        "perm_n_resamples_skew": permutation_result.n_resamples_skew,
        "perm_note": permutation_result.note,
    }


def _decide_direction(
    con_heights,
    dis1_heights,
    *,
    permutation_kwargs,
    rng,
):
    """Decide the concordant-vs-discordant1 direction for one sample pair.

    Args:
        con_heights: Concordant tree heights.
        dis1_heights: Discordant1 tree heights.
        permutation_kwargs: Keyword arguments forwarded to
            :func:`~ghostparser.orchestrator.permutation.run_studentized_permutation_test`.
        rng: A ``numpy.random.Generator`` for the permutation resampling.

    Returns:
        A tuple ``(direction, permutation_result)`` where ``permutation_result``
        is ``None`` when either group is empty and the test cannot run.
    """
    if not len(con_heights) or not len(dis1_heights):
        return DECISION_INCONCLUSIVE, None

    result = run_studentized_permutation_test(
        con_heights,
        dis1_heights,
        rng=rng,
        **(permutation_kwargs or {}),
    )
    return result.decision, result


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


def _build_triplet_seed_sequence(seed, triplet):
    """Build a deterministic per-triplet seed sequence from the run's base seed.

    The run carries a single configured seed; a stable per-triplet sequence is
    derived from ``(seed, triplet)`` via SHA-256. Deriving per triplet rather
    than drawing from one shared stream is what makes a run reproducible
    independently of how triplets are chunked across workers, so the ``taxon``,
    ``gene``, and serial paths all produce identical results.

    Args:
        seed: The global base seed, or ``None`` for a non-deterministic sequence.
        triplet: The triplet used to derive a stable per-triplet seed.

    Returns:
        A ``numpy.random.SeedSequence`` instance.
    """
    if seed is None:
        return np.random.SeedSequence()

    triplet_key = "|".join(triplet)
    digest = hashlib.sha256(f"{seed}|{triplet_key}".encode("utf-8")).digest()
    return np.random.SeedSequence(int.from_bytes(digest[:8], "little"))


def observation_from_subtree(
    subtree,
    triplet,
    tree_height_calculation_strategy,
    collect_summary_statistics=False,
):
    """Compute a ``(topology, tree_height, metrics)`` observation from a subtree.

    Reference implementation, paired with
    :func:`~.trees.extract_triplet_subtree`. The pipeline derives observations
    from a cached :class:`~.triplet_geometry.TripletGeometry`
    (:func:`~.triplet_geometry.geometry_observation`); this path walks a real
    DendroPy subtree instead, and the parity tests hold the two to agreement.

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


def _build_shape_statistics(canonical_heights, role_topologies, rng):
    """Describe the shape of each height group under its decision role.

    Args:
        canonical_heights: Heights keyed by canonical topology.
        role_topologies: The ``(con, dis1, dis2)`` topologies, in that order.
        rng: A ``numpy.random.Generator`` for the modality bootstrap.

    Returns:
        A flat dict keyed ``<group>_<field>`` over :data:`SHAPE_GROUP_LABELS`
        and :data:`~ghostparser.orchestrator.shape.SHAPE_FIELD_NAMES`.
    """
    statistics = {}
    for label, topology in zip(SHAPE_GROUP_LABELS, role_topologies):
        described = describe_shape(canonical_heights[topology], rng=rng)
        for field in SHAPE_FIELD_NAMES:
            statistics[f"{label}_{field}"] = described[field]
    return statistics


def _run_triplet_pipeline_from_observations(
    species_triplet,
    observations,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    permutation_kwargs=None,
    rng=None,
    species_topology=TOPOLOGY_AB,
    species_tree_newick=None,
    shape_diagnostics=False,
    shape_rng=None,
    skip_settled_gates=False,
):
    """Run the GhostParser Figure 6 pipeline from serialized observations.

    Args:
        species_triplet: The ``(A, B, C)`` triplet.
        observations: List of ``(topology, tree_height)`` tuples.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        permutation_kwargs: Keyword arguments forwarded to the permutation test.
        rng: A ``numpy.random.Generator`` for the permutation resampling.
        species_topology: The concordant (species-tree) topology.
        species_tree_newick: Species subtree Newick stored on the result.
        shape_diagnostics: Whether to describe each height group's shape.
        shape_rng: A ``numpy.random.Generator`` for the modality bootstrap.
        skip_settled_gates: When ``True``, run each test only where the cascade
            can still consult it.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values.
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

    dct_statistic, dct_p_value = run_discordant_count_test(
        n_dis1,
        n_dis2,
        method=discordant_test,
    )
    dct_significant = dct_p_value <= alpha_dct

    # Gate one settles the cascade on its own when it fails, so the tree-height
    # test below it decides nothing. The run-wide correction enrols only the
    # triplets that cleared gate one in the tree-height family, so a value
    # measured below a failed gate would be discarded anyway; the efficient mode
    # declines to compute it. The detailed mode reports it, and it still takes no
    # part in the family, which is what keeps the two modes' results identical.
    ks_statistic = None
    ks_p_value = None
    ks_significant = None
    if dct_significant or not skip_settled_gates:
        ks_statistic, ks_p_value = run_two_sample_ks_test(
            canonical_heights[dis1_topology],
            canonical_heights[con_topology],
        )
        ks_significant = ks_p_value <= alpha_ks

    con_heights = canonical_heights[con_topology]
    dis1_heights = canonical_heights[dis1_topology]

    shape_statistics = None
    if shape_diagnostics:
        shape_statistics = _build_shape_statistics(
            canonical_heights,
            (con_topology, dis1_topology, dis2_topology),
            shape_rng,
        )

    # Gates one and two settle the cascade by themselves whenever either is
    # non-significant, and no supported correction can revive a failed raw gate,
    # so the direction test cannot change the classification for those triplets.
    # The efficient mode skips it there and marks the empty block in
    # ``perm_note``; the detailed mode runs it for every triplet so the columns
    # stay populated for debugging.
    direction = None
    permutation_result = None
    perm_note = None
    if (dct_significant and ks_significant) or not skip_settled_gates:
        direction, permutation_result = _decide_direction(
            con_heights,
            dis1_heights,
            permutation_kwargs=permutation_kwargs,
            rng=rng,
        )
    else:
        perm_note = PERM_NOTE_NOT_CONSULTED

    classification, decision_gate = _classify_introgression(
        dct_significant,
        ks_significant,
        direction,
    )

    topology_metric_statistics = None
    if metric_buckets is not None:
        topology_metric_statistics = _build_topology_metric_statistics(
            canonical_triplet,
            TOPOLOGY_AB,
            metric_buckets,
            dis1_topology,
            dis2_topology,
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
        classification=classification,
        decision_gate=decision_gate,
        perm_decision=direction,
        shape_statistics=shape_statistics,
        analyzed_trees=analyzed_trees,
        topology_metric_statistics=topology_metric_statistics,
        **_permutation_result_fields(permutation_result, note=perm_note),
    )


def _bootstrap_policy(correction, family_size):
    """Decide how bootstrap iterations handle correction for one run.

    Args:
        correction: The configured p-value correction method.
        family_size: Number of triplets in the run, which is the size of each
            per-iteration testing family.

    Returns:
        A :class:`_BootstrapPolicy`.
    """
    return _BootstrapPolicy(
        correction=correction,
        family_size=max(1, int(family_size)),
        inline=is_inline_correction(correction),
    )


def _iteration_outcome(
    n_dis1,
    n_dis2,
    con_heights,
    dis1_heights,
    alpha_dct,
    alpha_ks,
    discordant_test,
    permutation_kwargs,
    rng,
    policy,
    collect_metrics=False,
):
    """Run one bootstrap iteration's cascade over its resampled heights.

    Uses the same thresholds as the point estimate; only the direction test is
    skipped, and only when an upstream gate has already failed under a
    correction that cannot lower a p-value.

    Args:
        n_dis1: Discordant1 count in the resample.
        n_dis2: Discordant2 count in the resample.
        con_heights: Concordant tree heights (Python list).
        dis1_heights: Discordant1 tree heights (Python list).
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        permutation_kwargs: Keyword arguments forwarded to the permutation test.
        rng: A ``numpy.random.Generator`` for the permutation resampling.
        policy: The run's :class:`_BootstrapPolicy`.
        collect_metrics: When ``True``, compute every statistic even where the
            cascade no longer needs it, for bootstrap-debug output.


    Returns:
        A ``(dct_p_value, ks_p_value, direction, metrics)`` tuple; each of the
        last three is ``None`` when that step did not run.
    """
    dct_statistic, dct_p_value = run_discordant_count_test(
        n_dis1, n_dis2, method=discordant_test
    )
    dct_failed = _iteration_corrected(dct_p_value, policy) > alpha_dct

    ks_statistic = None
    ks_p_value = None
    # An iteration measures the tree-height test on the same rule the point
    # estimate uses, so its per-iteration family is shaped like the run-wide one
    # it is compared against: only triplets that cleared the count gate.
    if collect_metrics or not dct_failed:
        ks_statistic, ks_p_value = run_two_sample_ks_test(dis1_heights, con_heights)

    ks_failed = (
        ks_p_value is not None
        and _iteration_corrected(ks_p_value, policy) > alpha_ks
    )

    direction = None
    if not dct_failed and not ks_failed:
        direction, _ = _decide_direction(
            con_heights,
            dis1_heights,
            permutation_kwargs=permutation_kwargs,
            rng=rng,
        )

    metrics = None
    if collect_metrics:
        metrics = (
            dct_statistic,
            dct_p_value,
            ks_statistic,
            ks_p_value,
            _mean(con_heights),
            _mean(dis1_heights),
        )
    return dct_p_value, ks_p_value, direction, metrics


def _iteration_corrected(p_value, policy):
    """Correct one iteration p-value as far as the policy allows.

    Args:
        p_value: The raw p-value.
        policy: The run's :class:`_BootstrapPolicy`.

    Only the short-circuit test reads this. An inline method can tighten the
    gate with the family size alone, which lets an iteration skip more work; a
    rank-based one falls back to the raw p-value, which is the conservative
    direction and so still safe to skip on.

    Returns:
        The corrected p-value on the inline path, else the raw one.
    """
    if policy.inline:
        return _adjust_p_value_inline(p_value, policy.correction, policy.family_size)
    return p_value


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


def _bootstrap_payload(
    bootstrap_value,
    fractions,
    debug,
    summary_only,
    deferred=None,
    studentized_ci=(None, None),
):
    """Assemble the bootstrap payload with optional finalized debug metrics.

    Args:
        bootstrap_value: The top class fraction, or ``None`` when the votes are
            still deferred.
        fractions: Per-class fractions, or ``None`` when still deferred.
        debug: The per-iteration debug accumulator, or ``None`` when disabled.
        summary_only: When ``True``, finalize debug metrics as compact summaries.
        deferred: A :class:`DeferredBootstrapRecord` when the votes cannot be
            tallied until every triplet's p-values are in hand.
        studentized_ci: The ``(low, high)`` percentile interval on the
            studentized mean difference.

    Returns:
        A dict with ``bootstrap_value``, ``all_bootstrap``, ``deferred``,
        ``studentized_ci``, and the six ``bootstrap_*`` debug entries (each
        ``None`` when ``debug`` is ``None``).
    """
    payload = {
        "bootstrap_value": bootstrap_value,
        "all_bootstrap": fractions,
        "deferred": deferred,
        "studentized_ci": studentized_ci,
    }
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


def _bootstrap_fractions(class_counts, iterations):
    """Turn per-class iteration counts into fractions.

    Args:
        class_counts: Mapping of classification to iteration count.
        iterations: Total iterations run.

    Returns:
        A dict of classification to fraction, ordered by classification.
    """
    counts = {label: 0 for label in _BOOTSTRAP_CLASSES}
    counts.update(class_counts)
    return {
        label: counts[label] / float(iterations) for label in sorted(counts)
    }


def _studentized_percentile_interval(values, alpha):
    """Compute the bootstrap-percentile interval on the studentized difference.

    Taken at the ``1 - 2 * alpha`` level, where an interval and a one-sided test
    at ``alpha`` agree.

    Args:
        values: Per-iteration studentized mean differences, possibly with
            ``nan`` entries for iterations that had too few observations.
        alpha: The permutation test's significance threshold.

    Returns:
        A ``(low, high)`` tuple, or ``(None, None)`` when no iteration produced
        a finite statistic.
    """
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (None, None)
    low, high = np.percentile(finite, [100.0 * alpha, 100.0 * (1.0 - alpha)])
    return (float(low), float(high))


def _run_bootstrap_iterations(
    observations,
    iterations,
    alpha_dct,
    alpha_ks,
    discordant_test,
    permutation_kwargs,
    rng,
    permutation_rng,
    policy,
    alpha_perm,
    debug_mode=False,
    summary_only=False,
):
    """Resample observations and aggregate per-iteration classifications.

    Resampling is vectorized and assumes the concordant topology is ``((A,B),C)``
    (always true on the orchestrator path). Every iteration parks its raw
    p-values in a :class:`DeferredBootstrapRecord` for the run-wide pass to
    correct: the tree-height family holds only the triplets that cleared the
    count gate, and how many those are is not known until the stream ends, so no
    method can classify an iteration while it runs.

    Args:
        observations: List of ``(topology, tree_height, metrics)`` tuples.
        iterations: Number of bootstrap iterations.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        permutation_kwargs: Keyword arguments forwarded to the permutation test,
            already scaled to the reduced bootstrap resample budget.
        rng: A ``numpy.random.Generator`` used for resampling.
        permutation_rng: A separate generator for the permutation tests, so the
            resample sequence stays independent of how many permutations run.
        policy: The run's :class:`_BootstrapPolicy`.
        alpha_perm: Sets the level of the studentized-difference interval.
        debug_mode: When ``True``, collect per-iteration debug metrics.
        summary_only: When ``True`` (and debug), emit compact summaries instead
            of full per-iteration lists.

    Returns:
        The payload dict built by :func:`_bootstrap_payload`.
    """
    debug = {key: [] for key in _BOOTSTRAP_DEBUG_KEYS} if debug_mode else None

    if iterations <= 0:
        return _bootstrap_payload(0.0, _bootstrap_fractions({}, 1), debug, summary_only)

    dct_p_values = np.ones(iterations, dtype=np.float64)
    ks_p_values = np.full(iterations, np.nan, dtype=np.float64)
    directions = np.zeros(iterations, dtype=np.int8)
    studentized = np.full(iterations, np.nan, dtype=np.float64)

    valid_count = len(observations)
    if valid_count:
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

    for index in range(iterations):
        if valid_count:
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
        else:
            n_dis1 = n_dis2 = 0
            con_heights = dis1_heights = []

        dct_p_value, ks_p_value, direction, metrics = _iteration_outcome(
            n_dis1,
            n_dis2,
            con_heights,
            dis1_heights,
            alpha_dct,
            alpha_ks,
            discordant_test,
            permutation_kwargs,
            permutation_rng,
            policy,
            collect_metrics=debug_mode,
        )

        studentized[index] = studentized_mean_diff(con_heights, dis1_heights)

        dct_p_values[index] = dct_p_value
        if ks_p_value is not None:
            ks_p_values[index] = ks_p_value
        directions[index] = _DIRECTION_CODES.get(direction, _DIRECTION_SKIPPED)

        if debug is not None:
            _append_iteration_debug(debug, metrics)

    studentized_ci = _studentized_percentile_interval(studentized, alpha_perm)

    return _bootstrap_payload(
        None,
        None,
        debug,
        summary_only,
        deferred=DeferredBootstrapRecord(
            dct_p_values=dct_p_values,
            ks_p_values=ks_p_values,
            directions=directions,
        ),
        studentized_ci=studentized_ci,
    )


def _finalize_triplet_analysis(
    triplet,
    observations,
    species_tree_topology,
    *,
    alpha_dct,
    alpha_ks,
    discordant_test,
    permutation_kwargs,
    bootstrap_options,
    triplet_seed,
    p_value_correction,
    family_size,
    shape_diagnostics=False,
    pipeline_mode=DEFAULT_PIPELINE_MODE,
    bootstrap=DEFAULT_BOOTSTRAP,
):
    """Run the base pipeline and bootstrap for one triplet's observations.

    Shared core of the two public entry points.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        observations: List of ``(topology, tree_height)`` tuples.
        species_tree_topology: Topology-only species subtree Newick, or ``None``.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        permutation_kwargs: Keyword arguments forwarded to the permutation test.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic resampling.
        p_value_correction: The run's correction method.
        family_size: Number of triplets in the run.
        shape_diagnostics: Whether to describe each height group's shape.
        pipeline_mode: ``efficient`` skips tests the cascade cannot consult;
            ``detailed`` runs every test for every triplet.
        bootstrap: Whether to run the bootstrap iterations at all.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values and bootstrap
        aggregates.
    """
    # One configured seed drives the whole run, but each consumer gets its own
    # child stream. That keeps the bootstrap resample sequence identical whether
    # or not the permutation test runs, and independent of how many resamples
    # any single adaptive test happens to draw.
    point_seed, bootstrap_seed, bootstrap_perm_seed, shape_seed = (
        _build_triplet_seed_sequence(triplet_seed, triplet).spawn(4)
    )

    permutation_kwargs = dict(permutation_kwargs or {})
    efficient = pipeline_mode == PIPELINE_MODE_EFFICIENT
    base_result = _run_triplet_pipeline_from_observations(
        triplet,
        observations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        permutation_kwargs=permutation_kwargs,
        rng=np.random.default_rng(point_seed),
        species_topology=TOPOLOGY_AB,
        species_tree_newick=species_tree_topology,
        shape_diagnostics=shape_diagnostics,
        shape_rng=np.random.default_rng(shape_seed),
        skip_settled_gates=efficient,
    )

    # Switching the bootstrap off is an instruction about what to compute, so it
    # holds in both pipeline modes: ``detailed`` declines to skip work the
    # cascade cannot consult, which is not the same as reinstating work the user
    # turned off. The bootstrap-sourced fields keep their ``None`` defaults.
    if not bootstrap:
        return base_result

    options = dict(bootstrap_options or {})
    iterations = int(options.get("iterations", DEFAULT_BOOTSTRAP_ITERATIONS))
    debug_mode = bool(options.get("debug_mode", DEFAULT_BOOTSTRAP_DEBUG_MODE))
    summary_only = (
        bool(options.get("summary_only", DEFAULT_BOOTSTRAP_SUMMARY_ONLY))
        if debug_mode
        else False
    )

    # Each bootstrap iteration re-runs the direction test, so it uses a fraction
    # of the configured resample budget.
    bootstrap_permutation_kwargs = dict(permutation_kwargs)
    scaled_min, scaled_max = bootstrap_resample_budget(
        permutation_kwargs.get("min_resamples", DEFAULT_PERMUTATION_MIN_RESAMPLES),
        permutation_kwargs.get("max_resamples", DEFAULT_PERMUTATION_MAX_RESAMPLES),
    )
    bootstrap_permutation_kwargs["min_resamples"] = scaled_min
    bootstrap_permutation_kwargs["max_resamples"] = scaled_max
    # Equivalence and inconclusiveness classify a triplet the same way, so an
    # iteration's vote never depends on which of the two it is.
    bootstrap_permutation_kwargs["equivalence_test"] = False

    bootstrap_payload = _run_bootstrap_iterations(
        observations,
        iterations=iterations,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        permutation_kwargs=bootstrap_permutation_kwargs,
        rng=np.random.default_rng(bootstrap_seed),
        permutation_rng=np.random.default_rng(bootstrap_perm_seed),
        policy=_bootstrap_policy(p_value_correction, family_size),
        alpha_perm=permutation_kwargs.get("alpha", DEFAULT_ALPHA_PERM),
        debug_mode=debug_mode,
        summary_only=summary_only,
    )

    bootstrap_gene_tree_heights = None
    if debug_mode:
        raw_heights = [tree_height for _, tree_height, _ in observations]
        bootstrap_gene_tree_heights = _finalize_bootstrap_metric(
            raw_heights, summary_only
        )

    ci_low, ci_high = bootstrap_payload["studentized_ci"]
    result = replace(
        base_result,
        bootstrap_stat_ci_low=ci_low,
        bootstrap_stat_ci_high=ci_high,
        bootstrap_value=bootstrap_payload["bootstrap_value"],
        all_bootstrap=bootstrap_payload["all_bootstrap"],
        bootstrap_deferred=bootstrap_payload["deferred"],
        bootstrap_dct_stats=bootstrap_payload["bootstrap_dct_stats"],
        bootstrap_dct_p_value=bootstrap_payload["bootstrap_dct_p_value"],
        bootstrap_ks_stats=bootstrap_payload["bootstrap_ks_stats"],
        bootstrap_ks_p_value=bootstrap_payload["bootstrap_ks_p_value"],
        bootstrap_con_summary=bootstrap_payload["bootstrap_con_summary"],
        bootstrap_dis_summary=bootstrap_payload["bootstrap_dis_summary"],
        bootstrap_gene_tree_heights=bootstrap_gene_tree_heights,
    )

    # Bootstrap votes are resolved by the run-wide pass, because the tree-height
    # family spans triplets. A family of one is the exception: it is fully known
    # here, so a caller analyzing a single triplet gets a finished result rather
    # than a deferred record. The pass is idempotent, so a one-triplet run that
    # reaches it again is unaffected.
    if family_size == 1:
        result = _apply_triplet_result_p_value_correction(
            [result], alpha_dct=alpha_dct, alpha_ks=alpha_ks, method=p_value_correction
        )[0]
    return result


def analyze_triplet_from_observations(
    triplet,
    observations,
    species_subtree=None,
    *,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    permutation_kwargs=None,
    bootstrap_options=None,
    triplet_seed=None,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    family_size=1,
    shape_diagnostics=False,
    pipeline_mode=DEFAULT_PIPELINE_MODE,
    bootstrap=DEFAULT_BOOTSTRAP,
):
    """Analyze one triplet from precomputed observations.

    This is the orchestrator's per-triplet unit: observations are computed once
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
        permutation_kwargs: Keyword arguments forwarded to the permutation test.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic resampling.
        p_value_correction: The run's correction method.
        family_size: Number of triplets in the run, which sizes the bootstrap's
            per-iteration testing families. Left at ``1`` this analyzes the
            triplet as a family of one.
        shape_diagnostics: Whether to describe each height group's shape.
        pipeline_mode: ``efficient`` skips tests the cascade cannot consult;
            ``detailed`` runs every test for every triplet.
        bootstrap: Whether to run the bootstrap iterations at all.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values. Under a
        rank-based correction its bootstrap votes are still deferred, and
        :func:`_apply_triplet_result_p_value_correction` fills them in.
    """
    species_tree_topology = _species_tree_topology_only_newick(species_subtree)
    return _finalize_triplet_analysis(
        triplet,
        observations,
        species_tree_topology,
        alpha_dct=alpha_dct,
        alpha_ks=alpha_ks,
        discordant_test=discordant_test,
        permutation_kwargs=permutation_kwargs,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
        p_value_correction=p_value_correction,
        family_size=family_size,
        shape_diagnostics=shape_diagnostics,
        pipeline_mode=pipeline_mode,
        bootstrap=bootstrap,
    )


def analyze_triplet(
    triplet,
    gene_subtrees,
    species_subtree=None,
    *,
    alpha_dct=DEFAULT_ALPHA_DCT,
    alpha_ks=DEFAULT_ALPHA_KS,
    discordant_test=DEFAULT_DISCORDANT_TEST,
    permutation_kwargs=None,
    tree_height_calculation_strategy=DEFAULT_TREE_HEIGHT_CALCULATION_STRATEGY,
    collect_summary_statistics=False,
    bootstrap_options=None,
    triplet_seed=None,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    family_size=1,
    shape_diagnostics=False,
    pipeline_mode=DEFAULT_PIPELINE_MODE,
    bootstrap=DEFAULT_BOOTSTRAP,
):
    """Analyze one triplet end-to-end from gene-subtree Newick strings.

    Serializes the gene subtrees into observations (reparsing each Newick), then
    runs the base pipeline and bootstrap. The returned p-values are uncorrected;
    run-wide correction is applied later.

    A run reaches the same pipeline through
    :func:`analyze_triplet_from_observations`, which takes observations the
    stream already has; this Newick-taking entry point is the reference the tests
    hold that one against.

    Args:
        triplet: The ``(A, B, C)`` triplet.
        gene_subtrees: Iterable of the triplet's rooted gene-subtree Newicks.
        species_subtree: The triplet's species subtree Newick (with branch
            lengths); stored topology-only on the result.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        discordant_test: ``chi-square`` or ``z-test``.
        permutation_kwargs: Keyword arguments forwarded to the permutation test.
        tree_height_calculation_strategy: Tree-height strategy to apply.
        collect_summary_statistics: When ``True``, gather per-triplet
            topology/metric summary statistics onto the result.
        bootstrap_options: Optional dict; ``iterations``/``debug_mode``/
            ``summary_only`` are read.
        triplet_seed: Optional base seed for deterministic resampling.
        p_value_correction: The run's correction method.
        family_size: Number of triplets in the run, which sizes the bootstrap's
            per-iteration testing families. Left at ``1`` this analyzes the
            triplet as a family of one.
        shape_diagnostics: Whether to describe each height group's shape.
        pipeline_mode: ``efficient`` skips tests the cascade cannot consult;
            ``detailed`` runs every test for every triplet.
        bootstrap: Whether to run the bootstrap iterations at all.

    Returns:
        A ``TripletPipelineResult`` with uncorrected p-values. Under a
        rank-based correction its bootstrap votes are still deferred, and
        :func:`_apply_triplet_result_p_value_correction` fills them in.
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
        permutation_kwargs=permutation_kwargs,
        bootstrap_options=bootstrap_options,
        triplet_seed=triplet_seed,
        p_value_correction=p_value_correction,
        family_size=family_size,
        shape_diagnostics=shape_diagnostics,
        pipeline_mode=pipeline_mode,
        bootstrap=bootstrap,
    )


# Classification codes for the vectorized bootstrap tally, in the same order as
# ``_BOOTSTRAP_CLASSES``.
_CLASS_CODES = {label: code for code, label in enumerate(_BOOTSTRAP_CLASSES)}


def _classification_codes(dct_significant, ks_significant, directions):
    """Apply the classification cascade to whole arrays of gate outcomes.

    The array form of :func:`_classify_introgression`; the two are pinned
    against each other by the truth-table test.

    Args:
        dct_significant: Boolean array of DCT outcomes.
        ks_significant: Boolean array of KS outcomes, aligned with it.
        directions: Integer array of direction codes, aligned with both.

    Returns:
        An integer array of ``_CLASS_CODES`` values.
    """
    codes = np.full(dct_significant.shape, _CLASS_CODES["ambiguous"], dtype=np.int8)
    both = dct_significant & ks_significant
    codes[~dct_significant] = _CLASS_CODES["no_introgression"]
    codes[dct_significant & ~ks_significant] = _CLASS_CODES["inflow_introgression"]
    codes[both & (directions == _DIRECTION_GREATER)] = _CLASS_CODES[
        "outflow_introgression"
    ]
    codes[both & (directions == _DIRECTION_LESS)] = _CLASS_CODES["ghost_introgression"]
    return codes


def _resolve_deferred_bootstrap(results, alpha_dct, alpha_ks, method):
    """Correct the parked per-iteration bootstrap p-values and tally the votes.

    Iteration ``i`` forms one testing family across triplets and is corrected as
    a whole before any of its members is classified.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        alpha_dct: Significance threshold for the discordant count test.
        alpha_ks: Significance threshold for the KS test.
        method: Correction method.

    Returns:
        A dict mapping each deferred result's index to its ``all_bootstrap``
        fractions, empty when nothing was deferred.
    """
    indices = [
        idx for idx, result in enumerate(results) if result.bootstrap_deferred is not None
    ]
    if not indices:
        return {}

    records = [results[idx].bootstrap_deferred for idx in indices]
    dct = np.stack([record.dct_p_values for record in records])
    ks = np.stack([record.ks_p_values for record in records])
    directions = np.stack([record.directions for record in records])
    n_triplets, iterations = dct.shape

    dct_significant = np.empty(dct.shape, dtype=bool)
    ks_significant = np.zeros(ks.shape, dtype=bool)
    for column in range(iterations):
        dct_significant[:, column] = (
            np.asarray(
                _adjust_p_values(dct[:, column], method=method, alpha=alpha_dct)
            )
            <= alpha_dct
        )
        # The tree-height family holds only the triplets this iteration's count
        # gate let through, matching the run-wide pass. Its size therefore varies
        # from iteration to iteration, which is why no method can correct it
        # while the stream is still running.
        present = np.isfinite(ks[:, column]) & dct_significant[:, column]
        if present.any():
            adjusted = _adjust_p_values(
                ks[present, column], method=method, alpha=alpha_ks
            )
            ks_significant[present, column] = np.asarray(adjusted) <= alpha_ks

    codes = _classification_codes(dct_significant, ks_significant, directions)

    resolved = {}
    for position, idx in enumerate(indices):
        counts = np.bincount(codes[position], minlength=len(_BOOTSTRAP_CLASSES))
        resolved[idx] = _bootstrap_fractions(
            {label: int(counts[code]) for label, code in _CLASS_CODES.items()},
            iterations,
        )
    return resolved


def _apply_triplet_result_p_value_correction(
    results,
    alpha_dct,
    alpha_ks,
    method=DEFAULT_P_VALUE_CORRECTION,
):
    """Apply run-wide p-value correction to DCT and KS p-values.

    Runs once across all triplets, recomputing each classification and bootstrap
    value from the corrected significance, and resolves any deferred bootstrap
    votes in the same pass.

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

    resolved_bootstrap = _resolve_deferred_bootstrap(
        results, alpha_dct, alpha_ks, method
    )

    dct_p_values = [result.dct_p_value for result in results]
    adjusted_dct = _adjust_p_values(
        dct_p_values,
        method=method,
        alpha=alpha_dct,
    )
    dct_significant_flags = [value <= alpha_dct for value in adjusted_dct]

    # The tree-height test only decides anything for a triplet whose count gate
    # cleared, so only those form its correction family. Enrolling the rest would
    # inflate the family with p-values nothing reads and push the corrected
    # values of the triplets that do matter towards non-significance. This is
    # also what makes the two pipeline modes agree: correction can only raise a
    # p-value, so the corrected survivors are a subset of the raw ones, and the
    # efficient mode measured the tree-height test for every raw survivor.
    ks_indices = [
        idx
        for idx, result in enumerate(results)
        if result.ks_p_value is not None and dct_significant_flags[idx]
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
        dct_significant = dct_significant_flags[idx]

        if idx in adjusted_ks_map:
            ks_p_value_corrected = adjusted_ks_map[idx]
            ks_significant = ks_p_value_corrected <= alpha_ks
        else:
            # Outside the family: either the test never ran, or the count gate
            # settled the triplet and its raw value takes no part in the
            # correction. Either way there is no corrected value to report.
            ks_p_value_corrected = None
            ks_significant = None

        # The permutation p-values are corrected inside the test, across the
        # one-tailed family, so they take no part in this run-wide pass. Its
        # stored decision is replayed rather than recomputed, which is sound
        # because no supported correction can lower a p-value: a triplet whose
        # raw gate failed cannot clear the corrected one, so the efficient mode
        # only ever skips a direction test the cascade could not have reached.
        # A skipped test leaves ``perm_decision`` as ``None`` -- guards and empty
        # groups still record a string -- so reaching gate three without one is
        # an invariant violation rather than an ambiguous call, and says so.
        if dct_significant and ks_significant and result.perm_decision is None:
            raise RuntimeError(
                f"Triplet {result.triplet} reached the direction gate with no "
                "permutation decision. Correction moved a gate the efficient "
                "pipeline mode assumed settled; rerun with pipeline_mode: "
                "detailed and report this."
            )
        classification, decision_gate = _classify_introgression(
            dct_significant,
            ks_significant,
            result.perm_decision,
        )
        all_bootstrap = resolved_bootstrap.get(idx, result.all_bootstrap)
        bootstrap_value = result.bootstrap_value
        if all_bootstrap is not None:
            bootstrap_value = float(all_bootstrap.get(classification, 0.0))

        adjusted_results.append(
            replace(
                result,
                dct_p_value_corrected=dct_p_value_corrected,
                dct_significant=dct_significant,
                ks_p_value_corrected=ks_p_value_corrected,
                ks_significant=ks_significant,
                classification=classification,
                decision_gate=decision_gate,
                all_bootstrap=all_bootstrap,
                bootstrap_value=bootstrap_value,
                bootstrap_deferred=None,
            )
        )

    return adjusted_results


def _format_optional_float(value):
    """Format an optional float for TSV output.

    Args:
        value: A float, or ``None``.

    Returns:
        The value at 12 significant digits, or an empty string for ``None``.
    """
    return "" if value is None else f"{value:.12g}"


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


def _format_shape_values(shape_statistics, columns):
    """Render one triplet's shape diagnostics in column order.

    Args:
        shape_statistics: The result's shape dict, or ``None``.
        columns: The ``<group>_<field>`` names, in output order.

    Returns:
        A list of strings; empty entries where a group was too small to
        describe, and ``n_modes`` rendered as an integer.
    """
    statistics = shape_statistics or {}
    values = []
    for column in columns:
        value = statistics.get(column)
        if value is None:
            values.append("")
        elif column.endswith("_n_modes"):
            values.append(str(int(value)))
        else:
            values.append(f"{float(value):.12g}")
    return values


def write_pipeline_results(
    results,
    output_filepath,
    dct_method=DEFAULT_DISCORDANT_TEST,
    p_value_correction=DEFAULT_P_VALUE_CORRECTION,
    bootstrap=DEFAULT_BOOTSTRAP,
    bootstrap_debug_mode=DEFAULT_BOOTSTRAP_DEBUG_MODE,
):
    """Write per-triplet results to a TSV file.

    Args:
        results: List of ``TripletPipelineResult`` objects.
        output_filepath: Destination TSV path.
        dct_method: Discordant test method, used to name the statistic column.
        p_value_correction: Correction method; names the corrected columns, and
            omits them under ``no``.
        bootstrap: Whether to include the bootstrap columns.
        bootstrap_debug_mode: Whether to also include the bootstrap-debug
            columns (only meaningful when ``bootstrap`` is ``True``).

    Raises:
        ValueError: If the method or correction value is unsupported.
    """
    if dct_method not in DISCORDANT_TEST_CHOICES:
        raise ValueError(
            f"Unsupported discordant test method: {dct_method}. "
            f"Choose one of: {', '.join(DISCORDANT_TEST_CHOICES)}"
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

    # Under ``no`` a corrected column would repeat its raw neighbour exactly.
    corrected = p_value_correction != "no"
    # Read off the results rather than re-plumbed, so the columns cannot
    # disagree with what the run actually measured.
    shape = any(result.shape_statistics is not None for result in results)
    shape_columns = shape_column_names()
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
        *([dct_corrected_column] if corrected else []),
        "dct_significant",
        "ks_statistic",
        "ks_p_value",
        *([ks_corrected_column] if corrected else []),
        "ks_significant",
        "perm_statistic",
        "perm_p_greater",
        "perm_p_less",
        *(
            [
                f"perm_p_greater_{p_value_correction}_corr",
                f"perm_p_less_{p_value_correction}_corr",
            ]
            if corrected
            else []
        ),
        "perm_p_tost",
        "bootstrap_stat_ci_low",
        "bootstrap_stat_ci_high",
        "perm_n_resamples",
        "perm_converged",
        "perm_null_skew",
        "perm_note",
    ]

    header.extend(
        [
            "perm_decision",
            "decision_gate",
            "classification",
            "inference",
            "analyzed_trees",
            *(shape_columns if shape else []),
        ]
    )

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
                    "bootstrap_con_mean",
                    "bootstrap_dis_mean",
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
                *(
                    [f"{result.dct_p_value_corrected:.12g}"]
                    if corrected
                    else []
                ),
                str(result.dct_significant),
                "" if result.ks_statistic is None else f"{result.ks_statistic:.12g}",
                "" if result.ks_p_value is None else f"{result.ks_p_value:.12g}",
                *(
                    [_format_optional_float(result.ks_p_value_corrected)]
                    if corrected
                    else []
                ),
                "" if result.ks_significant is None else str(result.ks_significant),
                _format_optional_float(result.perm_statistic),
                _format_optional_float(result.perm_p_greater),
                _format_optional_float(result.perm_p_less),
                *(
                    [
                        _format_optional_float(result.perm_p_greater_corrected),
                        _format_optional_float(result.perm_p_less_corrected),
                    ]
                    if corrected
                    else []
                ),
                _format_optional_float(result.perm_p_tost),
                _format_optional_float(result.bootstrap_stat_ci_low),
                _format_optional_float(result.bootstrap_stat_ci_high),
                "" if result.perm_n_resamples is None else str(result.perm_n_resamples),
                "" if result.perm_converged is None else str(result.perm_converged),
                _format_optional_float(result.perm_null_skew),
                result.perm_note or "",
            ]

            row.extend(
                [
                    result.perm_decision or "",
                    result.decision_gate or "",
                    result.classification,
                    inference,
                    str(result.analyzed_trees),
                    *(
                        _format_shape_values(result.shape_statistics, shape_columns)
                        if shape
                        else []
                    ),
                ]
            )

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
    classification, and (when bootstrap is enabled) the bootstrap value. Shape
    diagnostics are not part of this file; they appear only in the results TSV.

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
