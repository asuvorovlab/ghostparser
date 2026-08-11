"""Verify orchestrator.trees preprocessing: support filtering, outgroup rooting, triplet generation, and species-subtree construction.

Expected values are stated as explicit Newick literals derived by hand from each
fixture (see ``tests/TEST_IO.md``), not by comparison with another GhostParser
module.
"""

import dendropy

from ghostparser.orchestrator import trees as ptrees

_OUTGROUP = ["OUT"]


def test_clean_and_save_trees_preserves_a_well_supported_tree(
    orchestrator_species_tree, tmp_path
):
    """A support-free species tree passes the filter unchanged."""
    out_path = tmp_path / "clean_species.tree"
    ptrees.clean_and_save_trees(
        str(orchestrator_species_tree), str(out_path), min_avg_support=0.5
    )

    # No internal support labels -> nothing to filter; the topology and branch
    # lengths round-trip verbatim.
    assert out_path.read_text() == "(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);\n"
    assert len(ptrees.read_tree_file(str(out_path))) == 1


def test_clean_and_save_trees_drops_low_average_support(low_support_tree_file, tmp_path):
    """Trees whose mean internal support is below the threshold are dropped."""
    out_path = tmp_path / "clean.tree"
    ptrees.clean_and_save_trees(
        str(low_support_tree_file), str(out_path), min_avg_support=0.5
    )

    # Tree 0 supports (0.95, 0.99, 0.98) -> mean 0.9733 >= 0.5, kept.
    # Tree 1 supports (0.30, 0.20, 0.40) -> mean 0.3000 <  0.5, dropped.
    # read_tree_file returns BioPython Phylo trees.
    kept = ptrees.read_tree_file(str(out_path))
    assert len(kept) == 1
    labels = {terminal.name for terminal in kept[0].get_terminals()}
    assert labels == {"TaxaC", "TaxaD", "TaxaF", "TaxaG", "OutGroup"}
    # Support values are stripped from the cleaned output.
    assert "0.95" not in out_path.read_text()


def test_root_tree_on_outgroup_prunes_and_reports_ingroup(
    orchestrator_species_tree, tmp_path
):
    """Rooting on the outgroup removes it and folds its edge into the ingroup."""
    out_path = tmp_path / "clean_species.tree"
    ptrees.clean_and_save_trees(
        str(orchestrator_species_tree), str(out_path), min_avg_support=0.5
    )
    species_tree = ptrees.read_tree_file(str(out_path))[0]

    pruned, excluded, missing, ingroup = ptrees._root_tree_on_outgroup(
        species_tree, _OUTGROUP
    )

    assert excluded == {"OUT"}
    assert missing == set()
    assert sorted(ingroup) == ["A", "B", "C", "D"]
    # Input ((( A,B ),C):0.1,(D:0.1,OUT:0.5):0.2). Removing OUT dissolves the
    # (D,OUT) node, so D keeps 0.1 and the ((A,B),C) clade absorbs 0.1 + 0.2.
    assert ptrees.format_newick_with_precision(pruned) == (
        "(((A:0.1,B:0.1):0.1,C:0.2):0.3,D:0.1):0.5;"
    )


def test_generate_triplets_and_species_subtrees(orchestrator_species_tree, tmp_path):
    """Triplet enumeration and per-triplet species subtrees match hand derivation."""
    out_path = tmp_path / "clean_species.tree"
    ptrees.clean_and_save_trees(
        str(orchestrator_species_tree), str(out_path), min_avg_support=0.5
    )
    pruned, _excluded, _missing, ingroup = ptrees._root_tree_on_outgroup(
        ptrees.read_tree_file(str(out_path))[0], _OUTGROUP
    )
    pruned_newick = ptrees.format_newick_with_precision(pruned)
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
    assert species_map == {
        ("A", "B", "C"): "((A:0.1,B:0.1):0.1,C:0.2):0.8;",
        ("A", "B", "D"): "((A:0.1,B:0.1):0.4,D:0.1):0.5;",
        ("A", "C", "D"): "((A:0.2,C:0.2):0.3,D:0.1):0.5;",
        ("B", "C", "D"): "((B:0.2,C:0.2):0.3,D:0.1):0.5;",
    }


def test_clean_and_save_gene_trees_roots_every_tree_on_the_outgroup(
    orchestrator_gene_trees, tmp_path
):
    """Gene-tree cleaning reroots on the outgroup, zeroing its edge."""
    out_path = tmp_path / "clean_genes.tree"
    ptrees.clean_and_save_gene_trees(
        str(orchestrator_gene_trees), str(out_path), _OUTGROUP, min_avg_support=0.5
    )
    cleaned = ptrees._read_gene_trees_file(str(out_path))

    # All 12 fixture trees carry OUT and have no support labels, so all survive.
    assert len(cleaned) == 12
    # Rerooting at the outgroup attachment point: OUT's original edge is folded
    # into the ingroup clade's edge and OUT itself is left at length 0.
    # Tree 0: ingroup edge 0.10 + OUT 0.50 = 0.60.
    assert cleaned[0] == "((((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3):0.6,OUT:0);"
    # Tree 3: ingroup edge 0.10 + OUT 0.55 = 0.65.
    assert cleaned[3] == "((((B:0.4,C:0.4):0.1,A:0.6):0.1,D:0.35):0.65,OUT:0);"
    for line in cleaned:
        assert line.endswith("OUT:0);")


def test_extract_triplet_subtree_selects_the_triplet_taxa(orchestrator_gene_trees, tmp_path):
    """Subtree extraction keeps exactly the requested taxa with their distances."""
    out_path = tmp_path / "clean_genes.tree"
    ptrees.clean_and_save_gene_trees(
        str(orchestrator_gene_trees), str(out_path), _OUTGROUP, min_avg_support=0.5
    )
    first = ptrees._read_gene_trees_file(str(out_path))[0]
    tree = dendropy.Tree.get(data=first, schema="newick", preserve_underscores=True)

    subtree = ptrees.extract_triplet_subtree(tree, ("A", "B", "C"))
    assert subtree is not None
    labels = {
        leaf.taxon.label
        for leaf in subtree.leaf_node_iter()
        if leaf.taxon is not None
    }
    assert labels == {"A", "B", "C"}
