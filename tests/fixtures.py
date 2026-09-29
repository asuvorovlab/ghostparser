import csv

import pytest


@pytest.fixture
def simple_newick_file(tmp_path):
    """Create a temporary file with a simple Newick tree."""
    newick_str = (
        "(TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,"
        "TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,"
        "OutGroup:0.558);"
    )
    tree_file = tmp_path / "simple_tree.nwk"
    tree_file.write_text(newick_str)
    return tree_file


@pytest.fixture
def newick_with_support_file(tmp_path):
    """Create a temporary file with Newick tree containing support values."""
    newick_str = "(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup)0.85;"
    tree_file = tmp_path / "tree_with_support.nwk"
    tree_file.write_text(newick_str)
    return tree_file


@pytest.fixture
def multiple_trees_file(tmp_path):
    """Create a temporary file with multiple Newick trees."""
    newick_lines = [
        "(TaxaA:0.001,(TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001):0.012,OutGroup:0.558);",
        "(TaxaB:0.098,(TaxaC:0.001,TaxaD:0.001):0.001,OutGroup:0.558);",
        "((TaxaC:0.001,TaxaD:0.001):0.001,(TaxaB:0.098,TaxaA:0.001):0.012,OutGroup:0.558);",
    ]
    tree_file = tmp_path / "multiple_trees.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


@pytest.fixture
def low_support_tree_file(tmp_path):
    """Create a temporary file with trees having low average support."""
    newick_lines = [
        "(((TaxaC,TaxaD)0.95:0.110599,(TaxaF,TaxaG)0.99:1.860334)0.98:0.500000,OutGroup);",
        "(((TaxaC,TaxaD)0.3:0.110599,(TaxaF,TaxaG)0.2:1.860334)0.4:0.500000,OutGroup);",
    ]
    tree_file = tmp_path / "low_support_trees.nwk"
    tree_file.write_text("\n".join(newick_lines))
    return tree_file


@pytest.fixture()
def simple_species_tree(tmp_path):
    species_tree_str = (
        "(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);"
    )
    species_file = tmp_path / "species.tree"
    species_file.write_text(species_tree_str)
    return species_file


@pytest.fixture()
def simple_gene_trees(tmp_path):
    gene_trees_str = """((TaxaA:0.15,TaxaB:0.25):0.35,TaxaC:0.45);
((TaxaA:0.11,TaxaC:0.22):0.33,TaxaD:0.44);
((TaxaB:0.12,TaxaC:0.23):0.34,TaxaD:0.45);
"""
    gene_file = tmp_path / "genes.tree"
    gene_file.write_text(gene_trees_str)
    return gene_file


@pytest.fixture()
def summary_statistics_tsv(tmp_path):
    rows = [
        {
            "class": "100001",
            "dis1_topology": "BC",
            "feature_1": "1.0",
            "feature_2": "0.6",
            "feature_3": "0.4",
            "feature_4": "0.7",
        },
        {
            "class": "100001",
            "dis1_topology": "AC",
            "feature_1": "1.1",
            "feature_2": "0.7",
            "feature_3": "0.5",
            "feature_4": "0.75",
        },
        {
            "class": "010010",
            "dis1_topology": "BC",
            "feature_1": "1.2",
            "feature_2": "0.8",
            "feature_3": "0.6",
            "feature_4": "0.8",
        },
        {
            "class": "010010",
            "dis1_topology": "AC",
            "feature_1": "1.3",
            "feature_2": "0.85",
            "feature_3": "0.65",
            "feature_4": "0.82",
        },
        {
            "class": "001100",
            "dis1_topology": "BC",
            "feature_1": "1.4",
            "feature_2": "0.9",
            "feature_3": "0.7",
            "feature_4": "0.85",
        },
        {
            "class": "001100",
            "dis1_topology": "AC",
            "feature_1": "1.5",
            "feature_2": "0.95",
            "feature_3": "0.75",
            "feature_4": "0.9",
        },
        {
            "class": "110000",
            "dis1_topology": "BC",
            "feature_1": "1.6",
            "feature_2": "1.0",
            "feature_3": "0.8",
            "feature_4": "0.92",
        },
        {
            "class": "110000",
            "dis1_topology": "AC",
            "feature_1": "1.7",
            "feature_2": "1.05",
            "feature_3": "0.85",
            "feature_4": "0.95",
        },
        {
            "class": "110000",
            "dis1_topology": "BC",
            "feature_1": "1.8",
            "feature_2": "1.1",
            "feature_3": "0.9",
            "feature_4": "1.0",
        },
    ]
    path = tmp_path / "summary_statistics.tsv"
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture()
def summary_statistics_tsv_tuning(tmp_path):
    rows = []
    class_rows = [
        ("100001", "BC", 1.0, 0.6, 0.4, 0.7),
        ("010010", "AC", 1.2, 0.8, 0.6, 0.8),
        ("001100", "BC", 1.4, 0.9, 0.7, 0.85),
        ("110000", "AC", 1.6, 1.0, 0.8, 0.92),
    ]
    for repeat in range(4):
        for (
            class_label,
            topology,
            feature_1,
            feature_2,
            feature_3,
            feature_4,
        ) in class_rows:
            rows.append(
                {
                    "class": class_label,
                    "dis1_topology": topology,
                    "feature_1": f"{feature_1 + 0.01 * repeat:.2f}",
                    "feature_2": f"{feature_2 + 0.01 * repeat:.2f}",
                    "feature_3": f"{feature_3 + 0.01 * repeat:.2f}",
                    "feature_4": f"{feature_4 + 0.01 * repeat:.2f}",
                }
            )

    path = tmp_path / "summary_statistics_tuning.tsv"
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path
