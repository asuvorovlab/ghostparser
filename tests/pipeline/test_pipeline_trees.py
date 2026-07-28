"""Assert pipeline.trees preprocessing (cleaning, rooting/pruning, triplet normalization, gene cleaning) matches tree_parser on the same inputs."""

import dendropy

from ghostparser import tree_parser as tpz
from ghostparser.pipeline import trees as ptrees

_OUTGROUP = ["OUT"]


def _species_pipeline(species_path, out_path):
    """Run pipeline preprocessing to normalized triplets and a species subtree map.

    Args:
        species_path: Path to the input species tree.
        out_path: Path for the cleaned species tree.

    Returns:
        A tuple ``(triplets, species_map, ingroup_set, pruned_newick)``.
    """
    ptrees.clean_and_save_trees(str(species_path), str(out_path), min_avg_support=0.5)
    species_trees = ptrees.read_tree_file(str(out_path))
    pruned, _excluded, _missing, ingroup = ptrees._root_tree_on_outgroup(
        species_trees[0], _OUTGROUP
    )
    pruned_newick = ptrees.format_newick_with_precision(pruned)
    dendro = dendropy.Tree.get(
        data=pruned_newick, schema="newick", preserve_underscores=True
    )
    raw_triplets = ptrees.generate_triplets(sorted(ingroup), [])
    triplets, species_map, _skipped = ptrees._build_species_triplet_metadata(
        dendro, raw_triplets
    )
    return triplets, species_map, ingroup, pruned_newick


def _species_reference(species_path, out_path):
    """Run tree_parser preprocessing to normalized triplets and a species subtree map.

    Args:
        species_path: Path to the input species tree.
        out_path: Path for the cleaned species tree.

    Returns:
        A tuple ``(triplets, species_map, ingroup_set, pruned_newick)``.
    """
    tpz.clean_and_save_trees(str(species_path), str(out_path), min_avg_support=0.5)
    species_trees = tpz.read_tree_file(str(out_path))
    pruned, _excluded, _missing, ingroup = tpz._root_tree_on_outgroup(
        species_trees[0], _OUTGROUP
    )
    pruned_newick = tpz.format_newick_with_precision(pruned)
    dendro = dendropy.Tree.get(
        data=pruned_newick, schema="newick", preserve_underscores=True
    )
    raw_triplets = tpz.generate_triplets(sorted(ingroup), [])
    triplets, species_map, _skipped = tpz._build_species_triplet_metadata(
        dendro, raw_triplets
    )
    return triplets, species_map, ingroup, pruned_newick


def test_species_preprocessing_matches_tree_parser(pipeline_species_tree, tmp_path):
    """Pipeline species cleaning, rooting, and triplet normalization equal tree_parser."""
    p_triplets, p_map, p_ingroup, p_newick = _species_pipeline(
        pipeline_species_tree, tmp_path / "p_species.tree"
    )
    r_triplets, r_map, r_ingroup, r_newick = _species_reference(
        pipeline_species_tree, tmp_path / "r_species.tree"
    )

    assert p_newick == r_newick
    assert p_ingroup == r_ingroup
    assert p_triplets == r_triplets
    assert p_map == r_map


def test_gene_tree_cleaning_matches_tree_parser(pipeline_gene_trees, tmp_path):
    """Pipeline gene-tree cleaning and rooting equal tree_parser output."""
    p_out = tmp_path / "p_genes.tree"
    r_out = tmp_path / "r_genes.tree"
    ptrees.clean_and_save_gene_trees(
        str(pipeline_gene_trees), str(p_out), _OUTGROUP, min_avg_support=0.5
    )
    tpz.clean_and_save_gene_trees(
        str(pipeline_gene_trees), str(r_out), _OUTGROUP, min_avg_support=0.5
    )
    assert p_out.read_text() == r_out.read_text()
    assert ptrees._read_gene_trees_file(str(p_out)) == tpz._read_gene_trees_file(
        str(r_out)
    )
