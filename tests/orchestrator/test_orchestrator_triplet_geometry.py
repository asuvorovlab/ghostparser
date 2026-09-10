"""Triplet geometry read from a cached gene tree matches subtree extraction."""

import itertools

import dendropy
import pytest

from ghostparser.orchestrator.inference import observation_from_subtree
from ghostparser.orchestrator.trees import extract_triplet_subtree
from ghostparser.orchestrator.triplet_geometry import (
    TRIPLET_MISSING_TAXON,
    TRIPLET_RESOLVED,
    TRIPLET_UNRESOLVED,
    build_taxon_index,
    build_triplet_geometry,
    geometry_observation,
    triplet_resolution,
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


def _parse(newick):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def _dendropy_observation(newick, triplet, strategy, collect):
    """Run the current extract-then-measure path."""
    subtree = extract_triplet_subtree(_parse(newick), triplet)
    if subtree is None:
        return None
    return observation_from_subtree(subtree, triplet, strategy, collect)


def _geometry_observation(newick, triplet, strategy, collect):
    """Run the cached-geometry path."""
    taxon_index = build_taxon_index([triplet])
    geometry = build_triplet_geometry(_parse(newick), taxon_index)
    positions = tuple(taxon_index[label] for label in triplet)
    return geometry_observation(geometry, positions, strategy, collect)


@pytest.mark.parametrize("strategy", _STRATEGIES)
@pytest.mark.parametrize(
    "newick,triplet", [case[1:] for case in _PARITY_CASES],
    ids=[case[0] for case in _PARITY_CASES],
)
def test_geometry_matches_subtree_extraction(newick, triplet, strategy):
    """Cached geometry reproduces the extracted subtree's observation."""
    expected = _dendropy_observation(newick, triplet, strategy, True)
    actual = _geometry_observation(newick, triplet, strategy, True)

    assert expected is not None, "fixture should produce an observation"
    assert actual is not None

    assert actual[0] == expected[0]
    assert actual[1] == pytest.approx(expected[1], rel=1e-12, abs=1e-15)
    for metric, value in expected[2].items():
        assert actual[2][metric] == pytest.approx(value, rel=1e-12, abs=1e-15)


def test_geometry_omits_summary_metrics_when_not_collecting():
    """The metrics slot stays empty unless summary statistics are requested."""
    newick, triplet = _PARITY_CASES[0][1:]
    assert _geometry_observation(newick, triplet, "AVG", False)[2] is None
    assert _dendropy_observation(newick, triplet, "AVG", False)[2] is None


def test_geometry_on_a_hand_derived_tree():
    """Geometry values match a tree whose distances are computed by hand."""
    newick = "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);"

    # (P,Q) are sisters below node X. Subtree ((P:1,Q:1):2,R:3) puts every tip
    # 3.0 from its root, the internal branch at 2.0, and P..Q at 1.0+1.0.
    topology, height, metrics = _geometry_observation(newick, ("P", "Q", "R"), "AVG", True)
    assert topology == "((A,B),C)"
    assert height == pytest.approx(3.0)
    assert metrics["avg_tree_height"] == pytest.approx(3.0)
    assert metrics["internal_branch"] == pytest.approx(2.0)
    assert metrics["sister_distance"] == pytest.approx(2.0)

    # (P,R) are sisters below X while S sits across the root, so every tip is
    # 4.0 down, the internal branch is X's own 1.0, and P..R is 1.0+2.0+3.0.
    topology, height, metrics = _geometry_observation(newick, ("P", "R", "S"), "AVG", True)
    assert topology == "((A,B),C)"
    assert height == pytest.approx(4.0)
    assert metrics["internal_branch"] == pytest.approx(1.0)
    assert metrics["sister_distance"] == pytest.approx(6.0)


@pytest.mark.parametrize(
    "label,newick,triplet",
    [
        ("root_polytomy", "(A:1.0,B:1.0,C:1.0);", ("A", "B", "C")),
        ("absent_taxon", "((A:1.0,B:1.0):1.0,D:2.0);", ("A", "B", "C")),
    ],
)
def test_geometry_skips_exactly_what_extraction_skips(label, newick, triplet):
    """Both paths decline the same unusable triplets."""
    assert _dendropy_observation(newick, triplet, "AVG", True) is None
    assert _geometry_observation(newick, triplet, "AVG", True) is None


def test_zero_length_internal_branch_still_resolves():
    """A zero-length internal branch is resolved, not treated as a polytomy.

    Comparing LCA depths would tie here and drop the observation; comparing LCA
    node identity keeps it, which is what the subtree path does.
    """
    newick = "((A:1.0,B:1.0):0.0,C:1.0);"
    expected = _dendropy_observation(newick, ("A", "B", "C"), "INT", True)
    actual = _geometry_observation(newick, ("A", "B", "C"), "INT", True)

    assert expected is not None and actual is not None
    assert actual[0] == expected[0] == "((A,B),C)"
    assert actual[1] == pytest.approx(0.0)
    assert actual[2]["sister_distance"] == pytest.approx(expected[2]["sister_distance"])


def test_missing_edge_lengths_count_as_zero():
    """Absent Newick lengths contribute nothing, as _distance_to_root assumes."""
    topology, height, metrics = _geometry_observation(
        "(((A,B),C),(D,E));", ("A", "B", "C"), "AVG", True
    )
    assert topology == "((A,B),C)"
    assert height == pytest.approx(0.0)
    assert metrics["internal_branch"] == pytest.approx(0.0)


def test_one_cache_serves_every_triplet_in_the_tree():
    """A single cache answers all triplets, which is what makes it worth building."""
    newick = "(((P:1.0,Q:1.0):2.0,R:3.0):1.0,(S:2.0,T:2.0):2.0);"
    labels = ("P", "Q", "R", "S", "T")
    taxon_index = build_taxon_index([labels])
    geometry = build_triplet_geometry(_parse(newick), taxon_index)

    for triplet in itertools.combinations(labels, 3):
        positions = tuple(taxon_index[label] for label in triplet)
        actual = geometry_observation(geometry, positions, "AVG", True)
        expected = _dendropy_observation(newick, triplet, "AVG", True)
        assert (actual is None) == (expected is None)
        if expected is None:
            continue
        assert actual[0] == expected[0]
        assert actual[1] == pytest.approx(expected[1], rel=1e-12, abs=1e-15)


# A 9-taxon tree with an outgroup, nested clades, a long ladder, uneven branch
# lengths and one zero-length internal branch: 84 triplets over one topology.
_LARGE_TREE = (
    "((((T1:0.11,T2:0.19):0.23,(T3:0.07,T4:0.31):0.0):0.17,"
    "((T5:0.29,T6:0.13):0.41,T7:0.53):0.09):0.37,(T8:0.61,OUT:0.71):0.43);"
)


@pytest.mark.parametrize("strategy", _STRATEGIES)
def test_geometry_matches_dendropy_across_a_nine_taxon_tree(strategy):
    """Every triplet of a 9-taxon tree agrees with DendroPy extraction.

    Sweeps all 84 triplets rather than hand-picked shapes, so sister pairs on
    either side of the root, across the zero-length internal branch, and down
    the ladder are all exercised for each tree-height strategy.
    """
    labels = ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "OUT")
    triplets = list(itertools.combinations(labels, 3))
    assert len(triplets) == 84

    taxon_index = build_taxon_index([labels])
    geometry = build_triplet_geometry(_parse(_LARGE_TREE), taxon_index)

    compared = 0
    for triplet in triplets:
        expected = _dendropy_observation(_LARGE_TREE, triplet, strategy, True)
        actual = geometry_observation(
            geometry, tuple(taxon_index[label] for label in triplet), strategy, True
        )

        assert (actual is None) == (expected is None), triplet
        if expected is None:
            continue

        assert actual[0] == expected[0], triplet
        assert actual[1] == pytest.approx(expected[1], rel=1e-12, abs=1e-15), triplet
        for metric, value in expected[2].items():
            assert actual[2][metric] == pytest.approx(
                value, rel=1e-12, abs=1e-15
            ), (triplet, metric)
        compared += 1

    assert compared == len(triplets), "every triplet should yield an observation"


def test_triplet_resolution_agrees_with_the_observation_guards():
    """The diagnostic and the hot path decline exactly the same triplets.

    `triplet_resolution` restates `geometry_observation`'s guards so preflight
    can name a reason without the hot path tracking one. Nothing shares code
    between them, so this is what stops them drifting apart.
    """
    cases = [
        (_LARGE_TREE, ("T1", "T2", "T3"), TRIPLET_RESOLVED),
        ("((A:1.0,B:1.0):0.0,C:1.0);", ("A", "B", "C"), TRIPLET_RESOLVED),
        ("(A:1.0,B:1.0,C:1.0);", ("A", "B", "C"), TRIPLET_UNRESOLVED),
        ("((A:1.0,B:1.0):1.0,D:2.0);", ("A", "B", "C"), TRIPLET_MISSING_TAXON),
    ]
    for newick, triplet, expected in cases:
        taxon_index = build_taxon_index([triplet])
        geometry = build_triplet_geometry(_parse(newick), taxon_index)
        positions = tuple(taxon_index[label] for label in triplet)

        status = triplet_resolution(geometry, positions)
        observation = geometry_observation(geometry, positions, "AVG", False)

        assert status == expected, (newick, triplet)
        assert (status == TRIPLET_RESOLVED) == (observation is not None), (
            newick,
            triplet,
        )

    # And over a whole tree, the two must agree on every triplet.
    labels = ("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "OUT")
    taxon_index = build_taxon_index([labels])
    geometry = build_triplet_geometry(_parse(_LARGE_TREE), taxon_index)
    for triplet in itertools.combinations(labels, 3):
        positions = tuple(taxon_index[label] for label in triplet)
        resolved = triplet_resolution(geometry, positions) == TRIPLET_RESOLVED
        assert resolved == (
            geometry_observation(geometry, positions, "AVG", False) is not None
        ), triplet
