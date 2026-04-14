"""Test cases for the tree_parser module."""

import argparse
import copy
from io import StringIO
import tempfile
from pathlib import Path

import dendropy
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from ghostparser.tree_parser import (
    calculate_average_support,
    clean_and_save_gene_trees,
    clean_and_save_trees,
    extract_triplet_subtree,
    filter_triplets_by_taxa,
    format_newick_with_precision,
    generate_triplets,
    get_clean_filename,
    get_taxa_from_tree,
    MetricsLogger,
    process_gene_trees_for_triplets,
    read_triplet_filter_file,
    read_tree_file,
    remove_support_values,
    standardize_tree,
    _build_species_triplet_metadata,
    _resolve_runtime_args,
    write_clean_trees,
    write_triplet_gene_trees,
    write_triplet_gene_trees_multiprocess,
    write_triplet_gene_trees_streaming,
    write_triplets_to_file,
)


def _species_triplet_map(triplets):
    return {triplet: f"(({triplet[0]}:1,{triplet[1]}:1):1,{triplet[2]}:1);" for triplet in triplets}


def _tree_parser_args(**overrides):
    base = {
        "config_file": None,
        "species_tree_path": "species.nwk",
        "gene_trees_path": "genes.nwk",
        "outgroups": "Out1,Out2",
        "triplet_filter": None,
        "output_folder": None,
        "processes": None,
        "no_multiprocessing": False,
    }
    base.update(overrides)
    return argparse.Namespace(**base)


@pytest.mark.parametrize(
    "processes,expected_processes",
    [
        (None, 0),
        (6, 6),
    ],
)
def test_resolve_runtime_args_tree_parser_cli_processes(processes, expected_processes):
    resolved = _resolve_runtime_args(_tree_parser_args(processes=processes))
    # Paths are resolved to absolute paths
    assert resolved.species_tree == str(Path("species.nwk").resolve())
    assert resolved.gene_trees == str(Path("genes.nwk").resolve())
    assert resolved.outgroup == ["Out1", "Out2"]
    assert resolved.min_support_value == 0.5
    assert resolved.processes == expected_processes


def test_resolve_runtime_args_tree_parser_config_warns_and_ignores(tmp_path, capsys):
    config_path = tmp_path / "tree_parser_config.json"
    config_path.write_text(
        """
{
  "species_tree_path": "s.nwk",
  "gene_trees_path": "g.nwk",
  "outgroup": "OutA",
  "processes": 3,
  "no_multiprocessing": true,
  "min_support_value": 0.7
}
""".strip()
    )

    args = argparse.Namespace(
        config_file=str(config_path),
        species_tree_path="species.nwk",
        gene_trees_path=None,
        outgroups=None,
        triplet_filter=None,
        output_folder=None,
        processes=None,
        no_multiprocessing=False,
    )

    resolved = _resolve_runtime_args(args)
    captured = capsys.readouterr()

    assert "Warning: --config-file provided; CLI arguments not in config will be ignored" in captured.out
    # Paths are resolved to absolute paths
    assert resolved.species_tree == str(Path("s.nwk").resolve())
    assert resolved.gene_trees == str(Path("g.nwk").resolve())
    assert resolved.outgroup == ["OutA"]
    assert resolved.processes == 3
    assert resolved.no_multiprocessing is True
    assert resolved.min_support_value == 0.7


@pytest.mark.parametrize(
    "config_body,expected_processes",
    [
        (
            """
{
  "species_tree_path": "s.nwk",
  "gene_trees_path": "g.nwk",
  "outgroup": "OutA",
  "processes": 3,
  "no_multiprocessing": true,
  "min_support_value": 0.7
}
""".strip(),
            3,
        ),
        (
            """
{
  "species_tree_path": "s.nwk",
  "gene_trees_path": "g.nwk",
  "outgroup": "OutA"
}
""".strip(),
            0,
        ),
    ],
)
def test_resolve_runtime_args_tree_parser_config_processes_behavior(tmp_path, config_body, expected_processes):
    config_path = tmp_path / "tree_parser_config_processes.json"
    config_path.write_text(config_body)

    resolved = _resolve_runtime_args(
        _tree_parser_args(
            config_file=str(config_path),
            species_tree_path=None,
            gene_trees_path=None,
            outgroups=None,
            processes=9,
        )
    )
    assert resolved.processes == expected_processes

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def gene_trees_missing_outgroup_file(tmp_path):
    """Create a temporary file with one tree missing the outgroup taxon."""
    newick_lines = [
        "((TaxaA:0.1,TaxaB:0.2):0.3,OutGroup:0.4);",
        "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);",
    ]
    tree_file = tmp_path / "gene_trees_missing_outgroup.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


# ============================================================================
# Tests for read_tree_file
# ============================================================================


def test_read_tree_file_single_tree(simple_newick_file):
    """Test reading a single tree from a Newick file."""
    trees = read_tree_file(str(simple_newick_file))
    assert len(trees) == 1
    assert isinstance(trees[0], Tree)


def test_read_tree_file_multiple_trees(multiple_trees_file):
    """Test reading multiple trees from a Newick file."""
    trees = read_tree_file(str(multiple_trees_file))
    assert len(trees) == 3
    assert all(isinstance(tree, Tree) for tree in trees)


def test_read_tree_file_not_found():
    """Test error handling for non-existent file."""
    with pytest.raises(FileNotFoundError):
        read_tree_file("nonexistent_file.nwk")


@pytest.mark.parametrize(
    "filename,content",
    [
        ("invalid.nwk", "((TaxaA,TaxaB),TaxaC"),
        ("random.nwk", "(((,("),
        ("empty.nwk", ""),
    ],
)
def test_read_tree_file_invalid_inputs_raise_value_error(tmp_path, filename, content):
    invalid_file = tmp_path / filename
    invalid_file.write_text(content)

    with pytest.raises(ValueError, match="Invalid Newick format"):
        read_tree_file(str(invalid_file))


# ============================================================================
# Tests for calculate_average_support
# ============================================================================


def test_calculate_average_support_with_values(newick_with_support_file):
    """Test calculating average support from a tree with support values."""
    trees = read_tree_file(str(newick_with_support_file))
    avg_support = calculate_average_support(trees[0])

    # Average of 0.95, 0.99, 0.98, 0.85 should be ~0.9425
    assert avg_support is not None
    assert 0.93 < avg_support < 0.95


def test_calculate_average_support_no_values(simple_newick_file):
    """Test calculating average support from a tree without support values."""
    trees = read_tree_file(str(simple_newick_file))
    avg_support = calculate_average_support(trees[0])

    # Tree has no support values, should return None
    assert avg_support is None


# ============================================================================
# Tests for remove_support_values
# ============================================================================


def test_remove_support_values(newick_with_support_file):
    """Test removing support values from a tree."""
    trees = read_tree_file(str(newick_with_support_file))
    tree = trees[0]

    # Verify tree has support values before
    original_support = calculate_average_support(tree)
    assert original_support is not None

    # Remove support values
    cleaned_tree = remove_support_values(tree)

    # Verify support values are gone
    new_support = calculate_average_support(cleaned_tree)
    assert new_support is None


# ============================================================================
# Tests for standardize_tree
# ============================================================================


def test_standardize_tree_removes_support(newick_with_support_file):
    """Test that standardize_tree removes support values."""
    trees = read_tree_file(str(newick_with_support_file))
    tree = trees[0]

    standardized = standardize_tree(tree)
    avg_support = calculate_average_support(standardized)

    assert avg_support is None


def test_standardize_tree_preserves_branch_lengths(simple_newick_file):
    """Test that standardize_tree preserves branch lengths."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    # Get original branch lengths
    original_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]

    standardized = standardize_tree(tree)

    # Get standardized branch lengths
    standardized_lengths = [
        clade.branch_length for clade in standardized.find_clades() if clade.branch_length is not None
    ]

    assert len(original_lengths) == len(standardized_lengths)
    for orig, std in zip(original_lengths, standardized_lengths):
        assert abs(orig - std) < 1e-10


# ============================================================================
# Tests for format_newick_with_precision
# ============================================================================


def test_format_newick_with_precision_trailing_zeros(simple_newick_file):
    """Test that trailing zeros are removed in Newick output."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree, decimal_places=10)

    # Should not have many zeros in a row
    assert "0000000000" not in newick_str
    # Should end with semicolon
    assert newick_str.endswith(";")


def test_format_newick_with_precision_default_places(simple_newick_file):
    """Test formatting with default decimal places."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree)

    # Should be valid Newick format
    assert newick_str.endswith(";")
    assert "(" in newick_str
    assert ")" in newick_str


def test_format_newick_with_custom_precision(simple_newick_file):
    """Test formatting with custom decimal places."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    newick_str = format_newick_with_precision(tree, decimal_places=5)

    # Should be valid Newick format with custom precision
    assert newick_str.endswith(";")
    # Check that no more than 5 decimal places are used
    import re

    decimals = re.findall(r":0\.(\d+)", newick_str)
    for decimal_part in decimals:
        assert len(decimal_part) <= 5


# ============================================================================
# Tests for write_clean_trees
# ============================================================================


def test_write_clean_trees(simple_newick_file, tmp_path):
    """Test writing cleaned trees to a file."""
    trees = read_tree_file(str(simple_newick_file))
    output_file = tmp_path / "output_trees.nwk"

    write_clean_trees(trees, str(output_file))

    assert output_file.exists()

    # Verify the output file contains valid Newick
    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == len(trees)


def test_write_clean_trees_multiple(multiple_trees_file, tmp_path):
    """Test writing multiple cleaned trees to a file."""
    trees = read_tree_file(str(multiple_trees_file))
    output_file = tmp_path / "output_trees.nwk"

    write_clean_trees(trees, str(output_file))

    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == 3


# ============================================================================
# Tests for clean_and_save_trees
# ============================================================================


def test_clean_and_save_trees_filters_low_support(low_support_tree_file, tmp_path):
    """Test that clean_and_save_trees filters trees with low average support."""
    output_file = tmp_path / "cleaned_trees.nwk"

    cleaned, dropped = clean_and_save_trees(str(low_support_tree_file), str(output_file), min_avg_support=0.5)

    # Should keep high support tree and drop low support tree
    assert len(cleaned) == 1
    assert len(dropped) == 1
    assert 2 in dropped  # Second tree should be dropped


def test_clean_and_save_trees_no_filters(simple_newick_file, tmp_path):
    """Test clean_and_save_trees with trees that pass filter."""
    output_file = tmp_path / "cleaned_trees.nwk"

    cleaned, dropped = clean_and_save_trees(str(simple_newick_file), str(output_file), min_avg_support=0.0)

    assert len(cleaned) == 1
    assert len(dropped) == 0


def test_clean_and_save_trees_creates_output_file(simple_newick_file, tmp_path):
    """Test that clean_and_save_trees creates the output file."""
    output_file = tmp_path / "cleaned_trees.nwk"

    clean_and_save_trees(str(simple_newick_file), str(output_file))

    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_clean_and_save_gene_trees_discards_missing_outgroup(gene_trees_missing_outgroup_file, tmp_path):
    """Test that gene trees missing an outgroup taxon are discarded."""
    output_file = tmp_path / "cleaned_gene_trees.nwk"

    cleaned, dropped, rooted_count, missing_indices = clean_and_save_gene_trees(
        str(gene_trees_missing_outgroup_file),
        str(output_file),
        outgroup_taxa=["OutGroup"],
        min_avg_support=0.0,
    )

    assert len(cleaned) == 1
    assert dropped == {}
    assert rooted_count == 1
    assert missing_indices == [2]

    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == 1


# ============================================================================
# Tests for get_taxa_from_tree
# ============================================================================


def test_get_taxa_from_tree_correct_names(simple_newick_file):
    """Test that correct taxa names are extracted."""
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    taxa = get_taxa_from_tree(tree)
    expected_taxa = sorted(["TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG", "OutGroup"])

    assert taxa == expected_taxa


# ============================================================================
# Tests for generate_triplets
# ============================================================================


def test_generate_triplets_count():
    """Test that correct number of triplets are generated."""
    taxa = ["A", "B", "C", "D", "E"]
    outgroup = "E"

    triplets = generate_triplets(taxa, outgroup)

    # 4C3 = 4 triplets
    assert len(triplets) == 4


def test_generate_triplets_excludes_outgroup():
    """Test that outgroup is excluded from triplets."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup"]
    outgroup = "OutGroup"

    triplets = generate_triplets(taxa, outgroup)

    # Verify no triplet contains the outgroup
    for triplet in triplets:
        assert "OutGroup" not in triplet


def test_generate_triplets_content():
    """Test the content of generated triplets."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD"]
    outgroup = "TaxaD"

    triplets = generate_triplets(taxa, outgroup)
    triplets_set = set(triplets)

    expected = {
        ("TaxaA", "TaxaB", "TaxaC"),
    }

    assert len(triplets_set) == 1
    assert triplets_set == expected


def test_generate_triplets_large_set():
    """Test triplet generation with larger taxa set."""
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "TaxaE", "TaxaF", "TaxaG", "OutGroup"]  # 8 taxa
    outgroup = "OutGroup"

    triplets = generate_triplets(taxa, outgroup)

    # 7C3 = 35
    assert len(triplets) == 35
    # All triplets should be unique
    assert len(set(triplets)) == 35


# ============================================================================
# Tests for write_triplets_to_file
# ============================================================================


def test_write_triplets_to_file(tmp_path):
    """Test writing triplets to a file."""
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaB", "TaxaD"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplets.txt"

    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()

    # Verify file content
    lines = output_file.read_text().strip().split("\n")
    assert len(lines) == 3
    assert lines[0] == "TaxaA,TaxaB,TaxaC"
    assert lines[1] == "TaxaA,TaxaB,TaxaD"
    assert lines[2] == "TaxaA,TaxaC,TaxaD"


def test_write_triplets_to_file_empty(tmp_path):
    """Test writing empty triplet list."""
    triplets = []
    output_file = tmp_path / "triplets.txt"

    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()
    assert output_file.read_text() == ""


# ============================================================================
# Tests for get_clean_filename
# ============================================================================


def test_get_clean_filename_simple():
    """Test generating clean filename."""
    filepath = "/path/to/tree.nwk"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/processed_tree.nwk"


def test_get_clean_filename_different_extension():
    """Test clean filename with different extension."""
    filepath = "/path/to/mytrees.txt"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/processed_mytrees.txt"


def test_get_clean_filename_no_extension():
    """Test clean filename for file without extension."""
    filepath = "/path/to/treefile"
    clean_filepath = get_clean_filename(filepath)

    assert clean_filepath == "/path/to/processed_treefile"


# ============================================================================
# Integration Tests
# ============================================================================


def test_integration_full_workflow(simple_newick_file, tmp_path):
    """Test the full workflow: read, clean, and write trees."""
    # Read tree
    trees = read_tree_file(str(simple_newick_file))
    assert len(trees) == 1

    # Standardize
    standardized = standardize_tree(trees[0])

    # Write
    output_file = tmp_path / "output.nwk"
    write_clean_trees([standardized], str(output_file))

    # Verify output
    output_trees = read_tree_file(str(output_file))
    assert len(output_trees) == 1


def test_integration_triplets_workflow(simple_newick_file, tmp_path):
    """Test the workflow: read tree, extract taxa, and generate triplets."""
    # Read tree
    trees = read_tree_file(str(simple_newick_file))
    tree = trees[0]

    # Get taxa
    taxa = get_taxa_from_tree(tree)
    assert len(taxa) == 8

    # Generate triplets
    triplets = generate_triplets(taxa, "OutGroup")
    assert len(triplets) == 35  # 7C3

    # Write triplets
    output_file = tmp_path / "triplets.txt"
    write_triplets_to_file(triplets, str(output_file))

    assert output_file.exists()
    lines = output_file.read_text().strip().split("\n")
    assert len(lines) == 35


# ============================================================================
# Triplet Extraction Tests
# ============================================================================


def test_extract_triplet_subtree_all_taxa_present():
    """Test extracting a triplet subtree when all taxa are present."""
    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,(TaxaC:0.15,TaxaD:0.25):0.35);"
    tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)

    triplet = ("TaxaA", "TaxaB", "TaxaC")
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is not None
    subtree_taxa = {leaf.taxon.label for leaf in subtree.leaf_nodes() if leaf.taxon}
    assert subtree_taxa == {"TaxaA", "TaxaB", "TaxaC"}
    assert len(subtree.leaf_nodes()) == 3


def test_extract_triplet_subtree_missing_taxa():
    """Test extracting a triplet when not all taxa are present."""
    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.15);"
    tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)

    triplet = ("TaxaA", "TaxaB", "TaxaD")  # TaxaD is not in tree
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is None


def test_extract_triplet_subtree_preserves_branch_lengths():
    """Test that branch lengths are preserved/adjusted correctly."""
    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)

    triplet = ("TaxaA", "TaxaB", "TaxaC")
    subtree = extract_triplet_subtree(tree, triplet)

    assert subtree is not None

    # Convert to Newick to verify structure
    newick = format_newick_with_precision(subtree)
    # Should have the three taxa with branch lengths
    assert "TaxaA" in newick
    assert "TaxaB" in newick
    assert "TaxaC" in newick


def test_process_gene_trees_for_triplets():
    """Test processing multiple gene trees for triplet extraction."""
    # Create gene trees
    tree1_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree2_str = "((TaxaA:0.15,TaxaC:0.25):0.35,TaxaD:0.45);"
    tree3_str = "((TaxaB:0.12,TaxaC:0.22):0.32,TaxaD:0.42);"

    tree1 = dendropy.Tree.get(data=tree1_str, schema="newick", preserve_underscores=True)
    tree2 = dendropy.Tree.get(data=tree2_str, schema="newick", preserve_underscores=True)
    tree3 = dendropy.Tree.get(data=tree3_str, schema="newick", preserve_underscores=True)

    gene_trees = [tree1, tree2, tree3]

    # Define triplets
    triplets = [
        ("TaxaA", "TaxaB", "TaxaC"),
        ("TaxaA", "TaxaC", "TaxaD"),
        ("TaxaB", "TaxaC", "TaxaD"),
    ]

    # Process
    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    # Verify results
    assert len(triplet_gene_trees) == 3

    # First triplet should have 1 tree (from tree1 only)
    assert len(triplet_gene_trees[("TaxaA", "TaxaB", "TaxaC")]) == 1

    # Second triplet should have 1 tree (from tree2 only)
    assert len(triplet_gene_trees[("TaxaA", "TaxaC", "TaxaD")]) == 1

    # Third triplet should have 1 tree (from tree3 only)
    assert len(triplet_gene_trees[("TaxaB", "TaxaC", "TaxaD")]) == 1


def test_process_gene_trees_for_triplets_empty():
    """Test processing when no gene trees match triplets."""
    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)

    gene_trees = [tree]
    triplets = [("TaxaX", "TaxaY", "TaxaZ")]  # None of these taxa exist

    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    assert len(triplet_gene_trees[("TaxaX", "TaxaY", "TaxaZ")]) == 0


def test_write_triplet_gene_trees(tmp_path):
    """Test writing triplet gene trees to file in the specified format."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): [
            "(TaxaA:0.1,TaxaB:0.2,TaxaC:0.3);",
            "(TaxaA:0.15,TaxaB:0.25,TaxaC:0.35);",
        ],
        ("TaxaD", "TaxaE", "TaxaF"): [
            "(TaxaD:0.4,TaxaE:0.5,TaxaF:0.6);",
        ],
        ("TaxaG", "TaxaH", "TaxaI"): [],  # Empty triplet
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(
        triplet_gene_trees,
        str(output_file),
        species_triplet_trees=_species_triplet_map(list(triplet_gene_trees.keys())),
    )

    assert output_file.exists()

    content = output_file.read_text()
    lines = content.split("\n")

    # Check first triplet header
    assert lines[0] == "TaxaA,TaxaB,TaxaC\t2\t((TaxaA:1,TaxaB:1):1,TaxaC:1);"

    # Check blank line after header
    assert lines[1] == ""

    # Check first tree
    assert "TaxaA:0.1,TaxaB:0.2,TaxaC:0.3" in lines[2]

    # Check second tree
    assert "TaxaA:0.15,TaxaB:0.25,TaxaC:0.35" in lines[3]

    # Check separator line format: blank line then separator
    assert lines[4] == ""
    assert lines[5] == "=" * 60

    # Check second triplet header
    assert lines[6] == "TaxaD,TaxaE,TaxaF\t1\t((TaxaD:1,TaxaE:1):1,TaxaF:1);"


def test_write_triplet_gene_trees_includes_species_tree_header(tmp_path):
    """Test writing species triplet subtree in header with count."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): ["(TaxaA:0.1,TaxaB:0.2,TaxaC:0.3);"]
    }
    species_triplet_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): "((TaxaA:1,TaxaB:1):1,TaxaC:1);"
    }

    output_file = tmp_path / "triplet_gene_trees_species_header.txt"
    write_triplet_gene_trees(triplet_gene_trees, str(output_file), species_triplet_trees=species_triplet_trees)

    header = output_file.read_text().splitlines()[0]
    assert header == "TaxaA,TaxaB,TaxaC\t1\t((TaxaA:1,TaxaB:1):1,TaxaC:1);"


def test_write_triplet_gene_trees_empty_triplet(tmp_path):
    """Test writing when a triplet has no gene trees."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): [],
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(
        triplet_gene_trees,
        str(output_file),
        species_triplet_trees=_species_triplet_map(list(triplet_gene_trees.keys())),
    )

    assert output_file.exists()

    content = output_file.read_text()
    lines = content.split("\n")

    # Should have header with count 0 and blank line
    assert lines[0] == "TaxaA,TaxaB,TaxaC\t0\t((TaxaA:1,TaxaB:1):1,TaxaC:1);"
    assert lines[1] == ""
    # No trees after blank line


def test_integration_full_triplet_extraction_workflow(tmp_path):
    """Test the complete workflow from species tree to triplet gene trees."""
    # Create species tree file
    species_tree_str = "(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);"
    species_file = tmp_path / "species.tree"
    species_file.write_text(species_tree_str)

    # Create gene trees file
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    # Read species tree and generate triplets
    species_trees = read_tree_file(str(species_file))
    taxa = get_taxa_from_tree(species_trees[0])
    triplets = generate_triplets(taxa, "OutGroup")

    # Read gene trees with DendroPy
    gene_trees = []
    for line in gene_file.read_text().splitlines():
        newick_str = line.strip()
        if newick_str:
            gene_trees.append(
                dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
            )

    # Process gene trees for triplets
    triplet_gene_trees = process_gene_trees_for_triplets(gene_trees, triplets)

    # Write output
    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(
        triplet_gene_trees,
        str(output_file),
        species_triplet_trees=_species_triplet_map(list(triplet_gene_trees.keys())),
    )

    # Verify output
    assert output_file.exists()
    content = output_file.read_text()

    # Check that we have triplets in the output
    assert "TaxaA,TaxaB,TaxaC" in content
    assert "TaxaA,TaxaC,TaxaD" in content
    assert "TaxaB,TaxaC,TaxaD" in content

    # Check that counts are present
    assert "\t" in content  # Tab-separated counts


def test_write_triplet_gene_trees_multiprocess_with_workers(tmp_path):
    """Test multiprocessing triplet writer with multiple workers."""
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaD", "TaxaE", "TaxaF")]
    gene_trees_newick = [
        "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);",
        "((TaxaD:0.5,TaxaE:0.6):0.7,TaxaF:0.8);",
    ]

    output_file = tmp_path / "triplet_gene_trees_mp.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        gene_trees_newick,
        str(output_file),
        species_triplet_trees=_species_triplet_map(triplets),
        use_multiprocessing=True,
        processes=2,
    )

    # Verify output file exists
    assert output_file.exists()

    # Verify worker count is 2
    assert worker_count == 2

    # Verify statistics
    assert triplets_with_trees > 0
    assert total_subtrees > 0

    # Verify file content
    content = output_file.read_text()
    assert "TaxaA,TaxaB,TaxaC" in content
    assert "TaxaD,TaxaE,TaxaF" in content


def test_write_triplet_gene_trees_multiprocess_includes_species_header(tmp_path):
    """Test multiprocess writer includes species subtree in each triplet header."""
    triplets = [("TaxaA", "TaxaB", "TaxaC")]
    gene_trees_newick = ["((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"]
    species_triplet_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): "((TaxaA:1,TaxaB:1):1,TaxaC:1);"
    }

    output_file = tmp_path / "triplet_gene_trees_mp_species_header.txt"
    write_triplet_gene_trees_multiprocess(
        triplets,
        gene_trees_newick,
        str(output_file),
        species_triplet_trees=species_triplet_trees,
        use_multiprocessing=False,
    )

    first_line = output_file.read_text().splitlines()[0]
    assert first_line == "TaxaA,TaxaB,TaxaC\t1\t((TaxaA:1,TaxaB:1):1,TaxaC:1);"


def test_metrics_logger_context_manager(tmp_path):
    """Test MetricsLogger as a context manager."""
    metrics_file = tmp_path / "metrics.txt"

    with MetricsLogger(str(metrics_file)) as metrics:
        metrics.log("Test message 1")
        metrics.log("Test message 2")
        assert metrics.metrics_file is not None  # File should be open

    # File should be closed after exiting context
    assert metrics.metrics_file.closed

    # Verify file content
    content = metrics_file.read_text()
    assert "Test message 1" in content
    assert "Test message 2" in content


def test_metrics_logger_file_not_opened_before_enter(tmp_path):
    """Test that MetricsLogger doesn't open file until __enter__."""
    metrics_file = tmp_path / "metrics.txt"

    logger = MetricsLogger(str(metrics_file))
    # File should not exist yet
    assert not metrics_file.exists()
    assert logger.metrics_file is None

    # After entering context, file should exist
    with logger:
        assert logger.metrics_file is not None


def test_triplet_gene_trees_separator_format(tmp_path):
    """Test that separator lines are correctly formatted between triplets."""
    triplet_gene_trees = {
        ("TaxaA", "TaxaB", "TaxaC"): ["(TaxaA:0.1,TaxaB:0.2,TaxaC:0.3);"],
        ("TaxaD", "TaxaE", "TaxaF"): ["(TaxaD:0.4,TaxaE:0.5,TaxaF:0.6);"],
        ("TaxaG", "TaxaH", "TaxaI"): ["(TaxaG:0.7,TaxaH:0.8,TaxaI:0.9);"],
    }

    output_file = tmp_path / "triplet_gene_trees.txt"
    write_triplet_gene_trees(
        triplet_gene_trees,
        str(output_file),
        species_triplet_trees=_species_triplet_map(list(triplet_gene_trees.keys())),
    )

    content = output_file.read_text()
    lines = content.split("\n")

    # Count separator lines
    separator_count = sum(1 for line in lines if line == "=" * 60)

    # Should have 2 separators between 3 triplets
    assert separator_count == 2

    # Verify separator format: blank line, separator, next triplet
    for i, line in enumerate(lines):
        if line == "=" * 60:
            # Before separator should be blank line or tree content
            # After separator should be next triplet header or content
            if i > 0:
                assert lines[i - 1] == "" or ")" in lines[i - 1]  # Blank or tree line
            if i + 1 < len(lines):
                # Next line could be another header or blank
                assert lines[i + 1] == "" or "," in lines[i + 1]  # Blank or triplet header


def test_multiprocessing_triplet_writer_handles_empty_triplets(tmp_path):
    """Test that multiprocessing writer handles empty triplets correctly."""
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaX", "TaxaY", "TaxaZ")]
    gene_trees_newick = [
        "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);",
        # TaxaX, TaxaY, TaxaZ not present in any gene tree
    ]

    output_file = tmp_path / "triplet_gene_trees_empty.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        gene_trees_newick,
        str(output_file),
        species_triplet_trees=_species_triplet_map(triplets),
        use_multiprocessing=False,
    )

    # Verify that triplets with no matching trees are still written
    content = output_file.read_text()
    assert "TaxaX,TaxaY,TaxaZ\t0" in content

    # Only one triplet should have trees
    assert triplets_with_trees == 1
    assert total_subtrees > 0


# ============================================================================
# Triplet Utilities
# ============================================================================


def test_generate_triplets_multiple_outgroups():
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup1", "OutGroup2"]
    triplets = generate_triplets(taxa, ["OutGroup1", "OutGroup2"])

    for triplet in triplets:
        assert "OutGroup1" not in triplet
        assert "OutGroup2" not in triplet


def test_generate_triplets_outgroup_comma_separated_with_spaces():
    taxa = ["TaxaA", "TaxaB", "TaxaC", "TaxaD", "OutGroup1", "OutGroup2"]
    triplets = generate_triplets(taxa, "OutGroup1, OutGroup2")

    for triplet in triplets:
        assert "OutGroup1" not in triplet
        assert "OutGroup2" not in triplet


def test_read_triplet_filter_file_parses_valid_and_skips_invalid(tmp_path):
    filter_file = tmp_path / "triplets.txt"
    filter_file.write_text("TaxaA,TaxaB,TaxaC\nTaxaA,TaxaB\n\n")

    triplets, invalid_lines = read_triplet_filter_file(str(filter_file))

    assert triplets == [("TaxaA", "TaxaB", "TaxaC")]
    assert invalid_lines == [(2, "TaxaA,TaxaB")]


def test_filter_triplets_by_taxa_skips_missing_taxa():
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaB", "TaxaX")]
    taxa_set = {"TaxaA", "TaxaB", "TaxaC"}

    kept, skipped = filter_triplets_by_taxa(triplets, taxa_set)

    assert kept == [("TaxaA", "TaxaB", "TaxaC")]
    assert skipped == [(("TaxaA", "TaxaB", "TaxaX"), ["TaxaX"])]


def test_write_triplet_gene_trees_streaming(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees = write_triplet_gene_trees_streaming(
        triplets,
        str(gene_file),
        str(output_file),
        species_triplet_trees=_species_triplet_map(triplets),
    )

    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0
    content = output_file.read_text()
    assert "TaxaA,TaxaB,TaxaC" in content


def test_write_triplet_gene_trees_multiprocess_triplets_single_worker(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)

    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        str(gene_file),
        str(output_file),
        species_triplet_trees=_species_triplet_map(triplets),
        use_multiprocessing=False,
    )

    assert worker_count == 1
    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0


def test_write_triplet_gene_trees_multiprocess_accepts_list(tmp_path):
    gene_trees_newick = [
        "((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);",
        "((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);",
    ]
    triplets = [("TaxaA", "TaxaB", "TaxaC"), ("TaxaA", "TaxaC", "TaxaD")]
    output_file = tmp_path / "triplet_gene_trees.txt"

    total_subtrees, triplets_with_trees, worker_count = write_triplet_gene_trees_multiprocess(
        triplets,
        gene_trees_newick,
        str(output_file),
        species_triplet_trees=_species_triplet_map(triplets),
        use_multiprocessing=False,
    )

    assert worker_count == 1
    assert output_file.exists()
    assert total_subtrees > 0
    assert triplets_with_trees > 0
    chunk_dirs = list(tmp_path.glob(".*_chunks"))
    assert not chunk_dirs


def test_build_species_triplet_metadata_normalizes_abc(tmp_path):
    """Species triplets should be normalized so A and B are sisters."""
    species_tree = dendropy.Tree.get(
        data="((TaxaB:0.1,TaxaA:0.2):0.3,TaxaC:0.4);",
        schema="newick",
        preserve_underscores=True,
    )

    normalized_triplets, species_triplet_trees, skipped = _build_species_triplet_metadata(
        species_tree,
        [("TaxaA", "TaxaB", "TaxaC")],
    )

    assert skipped == []
    assert normalized_triplets == [("TaxaA", "TaxaB", "TaxaC")]
    assert ("TaxaA", "TaxaB", "TaxaC") in species_triplet_trees


def test_format_newick_with_precision_triplet_parser():
    tree_str = "((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4);"
    tree = dendropy.Tree.get(data=tree_str, schema="newick", preserve_underscores=True)

    newick = format_newick_with_precision(tree, decimal_places=5)
    assert newick.endswith(";")


# ============================================================================
# Triplet Equivalence
# ============================================================================


def _dist_to_root(node):
    dist = 0.0
    while node is not None and node.edge is not None:
        if node.edge.length is not None:
            dist += node.edge.length
        node = node.parent_node
    return dist


def _dendro_distance(tree, a, b):
    tree.is_rooted = True
    node_a = tree.find_node_with_taxon_label(a)
    node_b = tree.find_node_with_taxon_label(b)
    mrca = tree.mrca(taxon_labels=[a, b])
    return _dist_to_root(node_a) + _dist_to_root(node_b) - 2.0 * _dist_to_root(mrca)


def _bio_distance(tree, a, b):
    return tree.distance(a, b)


def _collapse_triplet_biopython(newick_str, triplet):
    tree = Phylo.read(StringIO(newick_str), "newick")
    collapsed = copy.deepcopy(tree)

    triplet_set = set(triplet)
    terminals = {terminal.name for terminal in collapsed.get_terminals()}
    if not triplet_set.issubset(terminals):
        return None

    for terminal in list(collapsed.get_terminals()):
        if terminal.name not in triplet_set:
            collapsed.prune(target=terminal)

    return collapsed


@pytest.mark.parametrize("a,b", [("A", "B"), ("A", "C"), ("B", "C")])
def test_triplet_branch_lengths_match(triplet_comparison_cases, a, b):
    for newick_str, triplet in triplet_comparison_cases:
        dendro_tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)

        dendro_subtree = extract_triplet_subtree(dendro_tree, triplet)
        assert dendro_subtree is not None

        dendro_newick = format_newick_with_precision(dendro_subtree)
        bio_subtree = Phylo.read(StringIO(dendro_newick), "newick")

        assert bio_subtree is not None

        bio_dist = _bio_distance(bio_subtree, a, b)
        dendro_dist = _dendro_distance(dendro_subtree, a, b)

        assert bio_dist == pytest.approx(dendro_dist)


def test_triplet_collapse_consistency_dendropy_vs_biopython(triplet_comparison_cases):
    for newick_str, triplet in triplet_comparison_cases:
        dendro_tree = dendropy.Tree.get(data=newick_str, schema="newick", preserve_underscores=True)
        dendro_subtree = extract_triplet_subtree(dendro_tree, triplet)
        assert dendro_subtree is not None

        bio_subtree = _collapse_triplet_biopython(newick_str, triplet)
        assert bio_subtree is not None

        for a, b in ((triplet[0], triplet[1]), (triplet[0], triplet[2]), (triplet[1], triplet[2])):
            dendro_dist = _dendro_distance(dendro_subtree, a, b)
            bio_dist = _bio_distance(bio_subtree, a, b)
            assert dendro_dist == pytest.approx(bio_dist, rel=0.0, abs=1e-12)
