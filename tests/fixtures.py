import pytest
import csv


@pytest.fixture
def simple_newick_file(tmp_path):
    """Create a temporary file with a simple Newick tree."""
    newick_str = "(TaxaA:0.001,(TaxaB:0.098,(((TaxaC:0.001,TaxaD:0.001):0.001,TaxaE:0.001):0.086,(TaxaF:0.001,TaxaG:0.001):0.032):0.001):0.012,OutGroup:0.558);"
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
    species_tree_str = "(((TaxaA:0.1,TaxaB:0.2):0.3,TaxaC:0.4):0.5,(TaxaD:0.6,OutGroup:0.7):0.8);"
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
def triplet_comparison_cases():
    return [
        (
            "((A:1.0,B:1.0):2.0,C:3.0,D:4.0);",
            ("A", "B", "C"),
        ),
        (
            "(((A:0.1,X:0.1):0.2,(B:0.1,Y:0.1):0.2):0.3,(C:0.1,Z:0.1):0.4);",
            ("A", "B", "C"),
        ),
        (
            "((A:0.5,B:0.5):0.5,(C:0.2,D:0.2):0.8);",
            ("A", "B", "C"),
        ),
    ]


@pytest.fixture()
def summary_statistics_tsv(tmp_path):
    rows = [
        {
            "triplet": "A,B,C",
            "abc_mapping": "A,B,C",
            "species_tree": "((A,B),C)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.9",
            "source_folder": "run_a",
            "classes": "100001",
            "concordant_avg_tree_height_mean": "1.0",
            "discordant1_avg_tree_height_mean": "0.6",
            "discordant2_avg_tree_height_mean": "0.4",
            "concordant_internal_branch_mean": "0.7",
            "discordant1_internal_branch_mean": "0.5",
            "discordant2_internal_branch_mean": "0.3",
            "concordant_sister_distance_mean": "0.8",
            "discordant1_sister_distance_mean": "0.45",
            "discordant2_sister_distance_mean": "0.25",
        },
        {
            "triplet": "D,E,F",
            "abc_mapping": "D,E,F",
            "species_tree": "((D,E),F)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.7",
            "source_folder": "run_a",
            "classes": "100001",
            "concordant_avg_tree_height_mean": "1.1",
            "discordant1_avg_tree_height_mean": "0.7",
            "discordant2_avg_tree_height_mean": "0.5",
            "concordant_internal_branch_mean": "0.75",
            "discordant1_internal_branch_mean": "0.55",
            "discordant2_internal_branch_mean": "0.35",
            "concordant_sister_distance_mean": "0.85",
            "discordant1_sister_distance_mean": "0.5",
            "discordant2_sister_distance_mean": "0.3",
        },
        {
            "triplet": "G,H,I",
            "abc_mapping": "G,H,I",
            "species_tree": "((G,H),I)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.8",
            "source_folder": "run_a",
            "classes": "010010",
            "concordant_avg_tree_height_mean": "1.2",
            "discordant1_avg_tree_height_mean": "0.8",
            "discordant2_avg_tree_height_mean": "0.6",
            "concordant_internal_branch_mean": "0.8",
            "discordant1_internal_branch_mean": "0.6",
            "discordant2_internal_branch_mean": "0.4",
            "concordant_sister_distance_mean": "0.9",
            "discordant1_sister_distance_mean": "0.55",
            "discordant2_sister_distance_mean": "0.35",
        },
        {
            "triplet": "J,K,L",
            "abc_mapping": "J,K,L",
            "species_tree": "((J,K),L)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.85",
            "source_folder": "run_a",
            "classes": "010010",
            "concordant_avg_tree_height_mean": "1.3",
            "discordant1_avg_tree_height_mean": "0.85",
            "discordant2_avg_tree_height_mean": "0.65",
            "concordant_internal_branch_mean": "0.82",
            "discordant1_internal_branch_mean": "0.62",
            "discordant2_internal_branch_mean": "0.42",
            "concordant_sister_distance_mean": "0.95",
            "discordant1_sister_distance_mean": "0.57",
            "discordant2_sister_distance_mean": "0.37",
        },
        {
            "triplet": "M,N,O",
            "abc_mapping": "M,N,O",
            "species_tree": "((M,N),O)",
            "classification": "outflow_introgression",
            "bootstrap_value": "0.65",
            "source_folder": "run_a",
            "classes": "001100",
            "concordant_avg_tree_height_mean": "1.4",
            "discordant1_avg_tree_height_mean": "0.9",
            "discordant2_avg_tree_height_mean": "0.7",
            "concordant_internal_branch_mean": "0.85",
            "discordant1_internal_branch_mean": "0.65",
            "discordant2_internal_branch_mean": "0.45",
            "concordant_sister_distance_mean": "0.98",
            "discordant1_sister_distance_mean": "0.6",
            "discordant2_sister_distance_mean": "0.4",
        },
        {
            "triplet": "P,Q,R",
            "abc_mapping": "P,Q,R",
            "species_tree": "((P,Q),R)",
            "classification": "outflow_introgression",
            "bootstrap_value": "0.6",
            "source_folder": "run_a",
            "classes": "001100",
            "concordant_avg_tree_height_mean": "1.5",
            "discordant1_avg_tree_height_mean": "0.95",
            "discordant2_avg_tree_height_mean": "0.75",
            "concordant_internal_branch_mean": "0.9",
            "discordant1_internal_branch_mean": "0.7",
            "discordant2_internal_branch_mean": "0.5",
            "concordant_sister_distance_mean": "1.0",
            "discordant1_sister_distance_mean": "0.62",
            "discordant2_sister_distance_mean": "0.42",
        },
        {
            "triplet": "S,T,U",
            "abc_mapping": "S,T,U",
            "species_tree": "((S,T),U)",
            "classification": "outflow_introgression",
            "bootstrap_value": "0.55",
            "source_folder": "run_a",
            "classes": "110000",
            "concordant_avg_tree_height_mean": "1.6",
            "discordant1_avg_tree_height_mean": "1.0",
            "discordant2_avg_tree_height_mean": "0.8",
            "concordant_internal_branch_mean": "0.92",
            "discordant1_internal_branch_mean": "0.72",
            "discordant2_internal_branch_mean": "0.52",
            "concordant_sister_distance_mean": "1.02",
            "discordant1_sister_distance_mean": "0.64",
            "discordant2_sister_distance_mean": "0.44",
        },
        {
            "triplet": "V,W,X",
            "abc_mapping": "V,W,X",
            "species_tree": "((V,W),X)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.5",
            "source_folder": "run_a",
            "classes": "110000",
            "concordant_avg_tree_height_mean": "1.7",
            "discordant1_avg_tree_height_mean": "1.05",
            "discordant2_avg_tree_height_mean": "0.85",
            "concordant_internal_branch_mean": "0.95",
            "discordant1_internal_branch_mean": "0.75",
            "discordant2_internal_branch_mean": "0.55",
            "concordant_sister_distance_mean": "1.05",
            "discordant1_sister_distance_mean": "0.66",
            "discordant2_sister_distance_mean": "0.46",
        },
        {
            "triplet": "Y,Z,AA",
            "abc_mapping": "Y,Z,AA",
            "species_tree": "((Y,Z),AA)",
            "classification": "ghost_introgression",
            "bootstrap_value": "0.45",
            "source_folder": "run_a",
            "classes": "110000",
            "concordant_avg_tree_height_mean": "1.8",
            "discordant1_avg_tree_height_mean": "1.1",
            "discordant2_avg_tree_height_mean": "0.9",
            "concordant_internal_branch_mean": "1.0",
            "discordant1_internal_branch_mean": "0.8",
            "discordant2_internal_branch_mean": "0.6",
            "concordant_sister_distance_mean": "1.1",
            "discordant1_sister_distance_mean": "0.68",
            "discordant2_sister_distance_mean": "0.48",
        },
    ]
    path = tmp_path / "summary_statistics.tsv"
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path
