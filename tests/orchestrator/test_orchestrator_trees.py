"""Verify orchestrator.trees preprocessing: support filtering, outgroup rooting, triplet generation, and species-subtree construction.

Expected values are stated as explicit Newick literals derived by hand from each
fixture (see ``tests/TEST_IO.md``), not by comparison with another GhostParser
module.
"""

import re

import dendropy
import pytest

from ghostparser.orchestrator import trees as ptrees

_OUTGROUP = ["OUT"]


def test_clean_and_save_trees_keeps_well_supported_trees_and_drops_the_rest(
    orchestrator_species_tree, low_support_tree_file, tmp_path
):
    """Trees at or above the mean-support threshold pass unchanged; the rest are dropped."""
    out_path = tmp_path / "clean_species.tree"
    ptrees.clean_and_save_trees(
        str(orchestrator_species_tree), str(out_path), min_avg_support=0.5
    )
    # No internal support labels -> nothing to filter; the topology and branch
    # lengths round-trip verbatim.
    assert out_path.read_text() == "(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);\n"
    assert len(ptrees.read_tree_file(str(out_path))) == 1

    out_path = tmp_path / "clean.tree"
    ptrees.clean_and_save_trees(
        str(low_support_tree_file), str(out_path), min_avg_support=0.5
    )
    # Tree 0 supports (0.95, 0.99, 0.98) -> mean 0.9733 >= 0.5, kept.
    # Tree 1 supports (0.30, 0.20, 0.40) -> mean 0.3000 <  0.5, dropped.
    kept = ptrees.read_tree_file(str(out_path))
    assert len(kept) == 1
    labels = {terminal.name for terminal in kept[0].get_terminals()}
    assert labels == {"TaxaC", "TaxaD", "TaxaF", "TaxaG", "OutGroup"}
    # Support values are stripped from the cleaned output.
    assert "0.95" not in out_path.read_text()


@pytest.mark.output
def test_clean_and_save_trees_quotes_labels_the_format_needs(tmp_path):
    """Labels a bare Newick token cannot hold are written quoted and read back intact.

    The processed trees are read back by the run itself, so a label with a
    space, a dot or a quote must survive the write in both parsers the run uses.
    """
    in_path = tmp_path / "quoted.tree"
    in_path.write_text("(('Homo sapiens':0.1,'Pan sp.':0.2):0.3,'O''Brien':0.4,Mus_musculus:0.5);\n")
    out_path = tmp_path / "clean.tree"
    ptrees.clean_and_save_trees(str(in_path), str(out_path), min_avg_support=0.5)

    # Only the labels that need it are quoted; the underscore label stays bare.
    assert out_path.read_text() == (
        "(('Homo sapiens':0.1,'Pan sp.':0.2):0.3,'O''Brien':0.4,Mus_musculus:0.5);\n"
    )
    expected = ["Homo sapiens", "Pan sp.", "O'Brien", "Mus_musculus"]
    biopython = ptrees.read_tree_file(str(out_path))[0]
    assert [t.name for t in biopython.get_terminals()] == expected
    dendro = dendropy.Tree.get(
        path=str(out_path), schema="newick", preserve_underscores=True
    )
    assert [leaf.taxon.label for leaf in dendro.leaf_node_iter()] == expected


@pytest.mark.parametrize(
    "newick, outgroups, expected_pruned, expected_distances, expected_order",
    [
        # A single outgroup: removing OUT dissolves the (D,OUT) node, so D
        # keeps 0.1 and the ((A,B),C) clade absorbs 0.1 + 0.2. That node is
        # the ingroup root, 0.5 from OUT.
        (
            "(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);",
            ["OUT"],
            "(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;",
            {"OUT": 0.5},
            ("OUT",),
        ),
        # OUT1 and OUT2 are not a clade as written (they sit on either side
        # of the file's root), but both branch off the ingroup at the node
        # joining OUT2 to it. Rooting there and pruning leaves the ingroup
        # clade with its own 0.4 edge. From the ingroup root OUT2 is
        # 0.4 + 0.2 away and OUT1 0.4 + 0.5 + 0.3, so OUT1 ranks first
        # though listed second.
        (
            "(OUT1:0.3,(OUT2:0.2,(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4):0.5);",
            ["OUT2", "OUT1"],
            "(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.1):0.4;",
            {"OUT1": 1.2, "OUT2": 0.6},
            ("OUT1", "OUT2"),
        ),
        # A root polytomy: both outgroups and both ingroup clades meet at the
        # same node, so that node is the ingroup's root and keeps both clades;
        # each outgroup's distance is its own edge.
        (
            "(OUT1:0.3,OUT2:0.2,(A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6);",
            ["OUT2", "OUT1"],
            "((A:0.1,B:0.1):0.4,(C:0.1,D:0.1):0.6);",
            {"OUT1": 0.3, "OUT2": 0.2},
            ("OUT1", "OUT2"),
        ),
        # No branch lengths, so there are no distances and the outgroups
        # keep their listed order, with OUTX, absent from the tree, last.
        (
            "(OUT1,(OUT2,(OUT3,((A,B),(C,D)))));",
            ["OUT3", "OUT2", "OUT1", "OUTX"],
            "((A,B),(C,D));",
            None,
            ("OUT3", "OUT2", "OUT1", "OUTX"),
        ),
    ],
    ids=["single_outgroup", "either_side_of_root", "root_polytomy", "no_lengths"],
)
def test_root_tree_on_outgroup_roots_where_the_outgroups_branch_off(
    newick, outgroups, expected_pruned, expected_distances, expected_order, tmp_path
):
    """The tree is rooted where the outgroups branch off, pruned of them, and
    the outgroups ranked farthest from the ingroup root first, or kept in
    listed order when the tree has no branch lengths.

    The outgroups need not be a clade as written, only branch off at one
    point; the removed edges fold into the ingroup's root edge. The ranking
    is the order that breaks gene-tree rooting ties.
    """
    path = tmp_path / "species.tree"
    path.write_text(newick + "\n")
    species_tree = ptrees.read_tree_file(str(path))[0]

    rooting = ptrees._root_tree_on_outgroup(species_tree, outgroups)
    pruned = rooting.tree

    if expected_distances is None:
        assert rooting.distances is None
    else:
        assert list(rooting.distances) == list(expected_distances)
        assert rooting.distances == pytest.approx(expected_distances)
    assert rooting.outgroup_order == expected_order
    assert sorted(rooting.ingroup) == ["A", "B", "C", "D"]
    assert ptrees.format_newick_with_precision(pruned) == expected_pruned


@pytest.mark.parametrize(
    "newick, expected_groups",
    [
        # X nests among the outgroups: OUT2 branches off between X and the
        # rest, so the taxa fall into two groups.
        (
            "(((A:1,B:1):1,C:1):1,(D:1,(OUT1:1,(OUT2:1,X:1):1):1):1);",
            (("A", "B", "C", "D"), ("X",)),
        ),
        # Each outgroup carries its own pocket of ingroup taxa.
        (
            "((OUT1:1,(A:1,B:1):1):1,(OUT2:1,(C:1,D:1):1):1);",
            (("A", "B"), ("C", "D")),
        ),
    ],
)
def test_root_tree_on_outgroup_rejects_outgroups_that_branch_off_twice(
    newick, expected_groups, tmp_path
):
    """Outgroups with other taxa between them raise, naming the separated groups."""
    path = tmp_path / "species.tree"
    path.write_text(newick + "\n")
    species_tree = ptrees.read_tree_file(str(path))[0]

    with pytest.raises(ptrees.OutgroupRootingError) as excinfo:
        ptrees._root_tree_on_outgroup(species_tree, ["OUT1", "OUT2"])

    assert excinfo.value.separated_groups == expected_groups


@pytest.mark.parametrize(
    "with_lengths, expected_subtrees",
    [
        (
            True,
            {
                ("A", "B", "C"): "((A:0.1,B:0.1):0.1,C:0.2):0.8;",
                ("A", "B", "D"): "((A:0.1,B:0.1):0.4,D:0.1):0.5;",
                ("A", "C", "D"): "((A:0.2,C:0.2):0.3,D:0.1):0.5;",
                ("B", "C", "D"): "((B:0.2,C:0.2):0.3,D:0.1):0.5;",
            },
        ),
        # The same tree with its branch lengths removed keeps them removed.
        (
            False,
            {
                ("A", "B", "C"): "((A,B),C);",
                ("A", "B", "D"): "((A,B),D);",
                ("A", "C", "D"): "((A,C),D);",
                ("B", "C", "D"): "((B,C),D);",
            },
        ),
    ],
    ids=["lengths", "no_lengths"],
)
def test_generate_triplets_and_species_subtrees(
    with_lengths, expected_subtrees, orchestrator_species_tree, tmp_path
):
    """Triplet enumeration and per-triplet species subtrees match hand
    derivation, with or without branch lengths in the species tree."""
    if not with_lengths:
        orchestrator_species_tree.write_text(
            re.sub(r":[0-9.]+", "", orchestrator_species_tree.read_text())
        )
    out_path = tmp_path / "clean_species.tree"
    ptrees.clean_and_save_trees(
        str(orchestrator_species_tree), str(out_path), min_avg_support=0.5
    )
    rooting = ptrees._root_tree_on_outgroup(
        ptrees.read_tree_file(str(out_path))[0], _OUTGROUP
    )
    ingroup = rooting.ingroup
    pruned_newick = ptrees.format_newick_with_precision(rooting.tree)
    dendro = dendropy.Tree.get(
        data=pruned_newick, schema="newick", preserve_underscores=True
    )

    # 4 ingroup taxa -> C(4,3) = 4 triplets, enumerated in sorted order.
    raw_triplets = ptrees.generate_triplets(sorted(ingroup), [])
    assert raw_triplets == [
        ("A", "B", "C"),
        ("A", "B", "D"),
        ("A", "C", "D"),
        ("B", "C", "D"),
    ]

    triplets, species_map, skipped = ptrees._build_species_triplet_metadata(
        dendro, raw_triplets
    )
    assert skipped == []
    # Species tree (((A,B),C),D): in every triplet the first two taxa are already
    # the sister pair, so A/B/C normalization is the identity here.
    assert triplets == raw_triplets
    assert species_map == expected_subtrees


def test_read_species_filter_file_collects_names_in_order(tmp_path):
    """Names may share a line or take one each; blanks and repeats are dropped."""
    path = tmp_path / "species.txt"
    path.write_text(" A, B \n\nC\nB,,D\n  \nA\n")

    assert ptrees.read_species_filter_file(str(path)) == ["A", "B", "C", "D"]


def test_clean_and_save_gene_trees_roots_each_tree_from_its_farthest_outgroup(
    orchestrator_gene_trees, tmp_path
):
    """Each tree is rooted at the common ancestor of its farthest outgroup and the others outside the ingroup, and pruned of them.

    The farthest is the outgroup with the longest mean path to the ingroup
    taxa in that gene tree, a missing branch length counting as 0 and the
    species-tree rank breaking ties. Rooting folds the outgroups' edge into the ingroup
    clade's edge; the outgroups are then cut away, so the written tree is the
    rooted ingroup. An outgroup that sits among the ingroup taxa once the
    tree is rooted on the farthest is pruned without being used, and the
    counts say which and how often.
    """
    out_path = tmp_path / "clean_genes.tree"
    cleaning = ptrees.clean_and_save_gene_trees(
        str(orchestrator_gene_trees), str(out_path), _OUTGROUP, min_avg_support=0.5
    )
    cleaned = ptrees._read_gene_trees_file(str(out_path))

    # All 12 fixture trees carry OUT and have no support labels, so all survive.
    assert len(cleaned) == 12
    assert cleaning.rooted_count == 12
    assert cleaning.farthest == {"OUT": 12}
    assert cleaning.rooted_on == {"OUT": 12}
    assert cleaning.tangled == {"OUT": 0}
    assert cleaning.tangled_trees == 0
    # Tree 0: ingroup edge 0.10 + OUT 0.50 = 0.60, and OUT is gone.
    assert cleaned[0] == "(((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6;"
    # Tree 3: ingroup edge 0.10 + OUT 0.55 = 0.65.
    assert cleaned[3] == "(((B:0.4,C:0.4):0.1,A:0.6):0.1,D:0.35):0.65;"
    assert all("OUT" not in line for line in cleaned)

    in_path = tmp_path / "genes.tree"
    in_path.write_text(
        # OUT1 alone.
        "((((A:1,B:1):1,C:1):1,D:1):1,OUT1:1);\n"
        # OUT2 nests among the ingroup; OUT1 is farther (mean path 4 against
        # 10/3), so the tree roots on OUT1 and OUT2 is pruned unused.
        "(((A:1,B:1):1,(C:1,OUT2:1):1):1,OUT1:1);\n"
        # OUT1 nests beside C, closer to the ingroup (mean path 3.5) than the
        # OUT2,OUT3 pair (5.5 each), so the pair roots the tree though the
        # species tree ranks OUT1 first.
        "((((A:1,B:1):1,(C:1,OUT1:1):1):1,D:1):1,(OUT2:1,OUT3:1):1);\n"
        # The same tree with OUT1 on a long branch (mean path 11.5): OUT1 is
        # the farthest, roots the tree, and the pair sits among the ingroup.
        "((((A:1,B:1):1,(C:1,OUT1:9):1):1,D:1):1,(OUT2:1,OUT3:1):1);\n"
        # OUT2 alone.
        "((((A:1,B:1):1,C:1):1,D:1):1,OUT2:1);\n"
        # No branch lengths: every mean path is 0, so the species-tree rank
        # makes OUT1 the farthest, and OUT2 sits with C.
        "(((A,B),(C,OUT2)),OUT1);\n"
        # The (A,B) edge lacks a length, read as 0: OUT1's mean path is 10/3
        # against OUT2's 8/3, so OUT1 roots it and OUT2 is tangled.
        "(((A:1,B:1),(C:1,OUT2:1):1):1,OUT1:1);\n"
        # No outgroup at all: dropped.
        "(((A:1,B:1):1,C:1):1,D:1);\n"
        # Nothing but outgroups: dropped.
        "(OUT1:1,(OUT2:1,OUT3:1):1);\n"
    )
    cleaning = ptrees.clean_and_save_gene_trees(
        str(in_path), str(out_path), ["OUT1", "OUT2", "OUT3"], min_avg_support=0.5
    )

    assert cleaning.rooted_count == 7
    assert cleaning.farthest == {"OUT1": 5, "OUT2": 2, "OUT3": 0}
    assert cleaning.rooted_on == {"OUT1": 5, "OUT2": 2, "OUT3": 1}
    assert cleaning.tangled == {"OUT1": 1, "OUT2": 4, "OUT3": 1}
    assert cleaning.tangled_trees == 5
    assert cleaning.missing_length_indices == [6, 7]
    assert cleaning.unrootable_indices == [8, 9]
    assert cleaning.dropped_trees == {}
    assert ptrees._read_gene_trees_file(str(out_path)) == [
        "(((A:1,B:1):1,C:1):1,D:1):2;",
        # OUT2 pruned from (C:1,OUT2:1):1 leaves C on a 2-long edge.
        "((A:1,B:1):1,C:2):2;",
        # OUT1 pruned from (C:1,OUT1:1):1 likewise; the root is where the
        # pair joins the ingroup, whose edge takes the pair's edge 1.
        "(((A:1,B:1):1,C:2):1,D:1):2;",
        # The root falls on OUT1's branch: C against the rest, where D takes
        # the edge the pair's removal leaves above it.
        "((D:2,(A:1,B:1):1):1,C:1):9;",
        "(((A:1,B:1):1,C:1):1,D:1):2;",
        # Rerooting gives the missing root edge 0; the rest stay unwritten.
        "((A,B),C):0;",
        # OUT1's edge joins the ingroup's, and C takes OUT2's parent edge.
        "((A:1,B:1),C:2):2;",
    ]
