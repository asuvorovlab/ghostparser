"""Triplet geometry read from a cached gene tree matches two independent references.

The parity tests hold the cached-geometry path to agreement with the DendroPy
subtree extraction it replaced and with a BioPython implementation written
from the definitions alone: root distances, pairwise distances and common
ancestors as ``Bio.Phylo`` computes them on the unpruned tree. Both live in
``tests/orchestrator/tree_references.py``. A hand-derived
table pins absolute values, so the three implementations cannot agree on a
shared mistake.
"""

import itertools
from io import StringIO

import pytest
from Bio import Phylo

from ghostparser.orchestrator.triplet_geometry import (
    TRIPLET_MISSING_TAXON,
    TRIPLET_RESOLVED,
    TRIPLET_UNRESOLVED,
    build_taxon_index,
    build_triplet_geometry,
    geometry_observation,
    triplet_resolution,
)
from ghostparser.triplet_utils import TOPOLOGY_AB
from tests.orchestrator.tree_references import (
    biopython_observation,
    dendropy_observation,
)

_STRATEGIES = ("AVG", "A", "B", "C", "SIS", "INT")

# (label, newick, triplet) covering nested pairs, pairs spanning the root,
# extra taxa pruned away, an unbalanced ladder, and absent edge lengths.
_PARITY_CASES = [
    ("nested_pair", "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);", ("P", "Q", "R")),
    ("spans_root", "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);", ("P", "R", "S")),
    ("both_sides", "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);", ("P", "S", "T")),
    ("odd_first", "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);", ("S", "P", "Q")),
    ("ladder", "((((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6):0.7,E:0.8);", ("A", "C", "D")),
    ("deep_ladder", "((((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6):0.7,E:0.8);", ("B", "D", "E")),
    ("pruned_taxa", "(((A:0.1,X:0.1):0.2,(B:0.1,Y:0.1):0.2):0.3,(C:0.1,Z:0.1):0.4);", ("A", "B", "C")),
    ("no_lengths", "(((A,B),C),(D,E));", ("A", "B", "C")),
    ("mixed_lengths", "(((A:0.5,B):0.25,C:1.0):0.5,D:2.0);", ("A", "B", "C")),
    ("zero_internal", "((A:1.0,B:1.0):0.0,C:1.0);", ("A", "B", "C")),
    ("tiny_internal", "(((A:1e-12,B:1e-12):1e-13,C:1.0000000000001):2.5,D:3.0);", ("A", "B", "C")),
]

# A 9-taxon tree with an outgroup, nested clades, a long ladder, uneven branch
# lengths and one zero-length internal branch: 84 triplets over one topology.
_LARGE_TREE = (
    "((((T1:0.11,T2:0.19):0.23,(T3:0.07,T4:0.31):0.0):0.17,"
    "((T5:0.29,T6:0.13):0.41,T7:0.53):0.09):0.37,(T8:0.61,OUT:0.71):0.43);"
)
_LARGE_LABELS = ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "OUT")


def _parse_for_geometry(newick):
    return Phylo.read(StringIO(newick), "newick")


_REFERENCES = {
    "dendropy": dendropy_observation,
    "biopython": biopython_observation,
}


def _geometry_observation(newick, triplet, strategy, collect):
    """Run the cached-geometry path."""
    taxon_index = build_taxon_index([triplet])
    geometry = build_triplet_geometry(_parse_for_geometry(newick), taxon_index)
    positions = tuple(taxon_index[label] for label in triplet)
    return geometry_observation(geometry, positions, strategy, collect)


def _assert_same_observation(actual, expected, context):
    """Compare a cached-geometry observation with a reference one."""
    assert actual[0] == expected[0], context
    assert actual[1] == pytest.approx(expected[1], rel=1e-12, abs=1e-15), context
    for metric, value in expected[2].items():
        assert actual[2][metric] == pytest.approx(value, rel=1e-12, abs=1e-15), (
            context,
            metric,
        )


@pytest.mark.parity
@pytest.mark.parametrize(
    "newick,triplet", [case[1:] for case in _PARITY_CASES],
    ids=[case[0] for case in _PARITY_CASES],
)
def test_geometry_matches_each_reference(newick, triplet):
    """Cached geometry reproduces both references under every height strategy."""
    for reference, strategy in itertools.product(sorted(_REFERENCES), _STRATEGIES):
        expected = _REFERENCES[reference](newick, triplet, strategy, True)
        actual = _geometry_observation(newick, triplet, strategy, True)

        assert expected is not None, "fixture should produce an observation"
        assert actual is not None
        _assert_same_observation(actual, expected, (reference, strategy))


@pytest.mark.parity
def test_geometry_matches_each_reference_across_a_nine_taxon_tree():
    """Every triplet of a 9-taxon tree agrees with each reference under each strategy.

    Sweeps all 84 triplets rather than hand-picked shapes, so sister pairs on
    either side of the root, across the zero-length internal branch, and down
    the ladder are all exercised. One cache built over the whole tree answers
    every triplet, and the diagnostic ``triplet_resolution`` must call each of
    them resolved exactly where the observation path yields one.
    """
    triplets = list(itertools.combinations(_LARGE_LABELS, 3))
    assert len(triplets) == 84

    taxon_index = build_taxon_index([_LARGE_LABELS])
    geometry = build_triplet_geometry(_parse_for_geometry(_LARGE_TREE), taxon_index)

    for triplet in triplets:
        positions = tuple(taxon_index[label] for label in triplet)
        assert triplet_resolution(geometry, positions) == TRIPLET_RESOLVED, triplet
        for reference, strategy in itertools.product(sorted(_REFERENCES), _STRATEGIES):
            expected = _REFERENCES[reference](_LARGE_TREE, triplet, strategy, True)
            actual = geometry_observation(geometry, positions, strategy, True)

            assert expected is not None, triplet
            assert actual is not None, triplet
            _assert_same_observation(actual, expected, (triplet, reference, strategy))


@pytest.mark.parametrize(
    "newick, triplet, strategy, topology, height, metrics",
    [
        # (P,Q) are sisters below node X. Subtree ((P:1,Q:1):2,R:3) puts every
        # tip 3.0 from its root, the internal branch at 2.0, and P..Q at 1+1.
        (
            "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);",
            ("P", "Q", "R"),
            "AVG",
            TOPOLOGY_AB,
            3.0,
            {"avg_tree_height": 3.0, "internal_branch": 2.0, "sister_distance": 2.0},
        ),
        # (P,R) are sisters below X while S sits across the root, so every tip
        # is 4.0 down, the internal branch is X's own 1.0, and P..R is 1+2+3.
        (
            "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);",
            ("P", "R", "S"),
            "AVG",
            TOPOLOGY_AB,
            4.0,
            {"avg_tree_height": 4.0, "internal_branch": 1.0, "sister_distance": 6.0},
        ),
        # A zero-length internal branch is resolved from LCA identity, not from
        # a depth comparison that would tie: INT is exactly 0 and A..B is 1+1.
        (
            "((A:1.0,B:1.0):0.0,C:1.0);",
            ("A", "B", "C"),
            "INT",
            TOPOLOGY_AB,
            0.0,
            {"avg_tree_height": 1.0, "internal_branch": 0.0, "sister_distance": 2.0},
        ),
        # Absent Newick lengths contribute nothing, so every distance is 0.
        (
            "(((A,B),C),(D,E));",
            ("A", "B", "C"),
            "AVG",
            TOPOLOGY_AB,
            0.0,
            {"avg_tree_height": 0.0, "internal_branch": 0.0, "sister_distance": 0.0},
        ),
    ],
    ids=["nested_pair", "spans_root", "zero_internal", "no_lengths"],
)
def test_geometry_matches_hand_derived_values(
    newick, triplet, strategy, topology, height, metrics
):
    """Geometry values match distances computed by hand from the Newick."""
    observed_topology, observed_height, observed_metrics = _geometry_observation(
        newick, triplet, strategy, True
    )
    assert observed_topology == topology
    assert observed_height == pytest.approx(height)
    assert observed_metrics == pytest.approx(metrics)
    # The metrics slot stays empty unless summary statistics are requested.
    assert _geometry_observation(newick, triplet, strategy, False)[2] is None


@pytest.mark.parity
@pytest.mark.parametrize(
    "newick, triplet, status",
    [
        ("(A:1.0,B:1.0,C:1.0);", ("A", "B", "C"), TRIPLET_UNRESOLVED),
        ("((A:1.0,B:1.0):1.0,D:2.0);", ("A", "B", "C"), TRIPLET_MISSING_TAXON),
        ("((A:1.0,B:1.0):0.0,C:1.0);", ("A", "B", "C"), TRIPLET_RESOLVED),
    ],
    ids=["root_polytomy", "absent_taxon", "zero_internal"],
)
def test_geometry_skips_exactly_what_each_reference_skips(newick, triplet, status):
    """The cached path, each reference and the diagnostic decline the same triplets.

    ``triplet_resolution`` restates ``geometry_observation``'s guards so the
    preflight can name a reason without the hot path tracking one; nothing
    shares code between them, so this is what stops them drifting apart.
    """
    taxon_index = build_taxon_index([triplet])
    geometry = build_triplet_geometry(_parse_for_geometry(newick), taxon_index)
    positions = tuple(taxon_index[label] for label in triplet)
    resolved = status == TRIPLET_RESOLVED

    assert triplet_resolution(geometry, positions) == status
    assert (geometry_observation(geometry, positions, "AVG", True) is not None) is resolved
    for reference in sorted(_REFERENCES):
        assert (
            _REFERENCES[reference](newick, triplet, "AVG", True) is not None
        ) is resolved, reference
