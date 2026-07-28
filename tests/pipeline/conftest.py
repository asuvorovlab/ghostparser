"""Shared fixtures for pipeline parity tests: a species tree and gene trees that all carry the outgroup so rooting yields real observations."""

import pytest

_SPECIES_TREE = "(((A:0.1,B:0.1):0.1,C:0.2):0.1,(D:0.1,OUT:0.5):0.2);"

# Gene trees include A,B,C,D,OUT so every ingroup triplet extracts a subtree and
# rooting on OUT always succeeds. Topologies and branch lengths vary across lines.
_GENE_TREES = [
    "((((A:0.10,B:0.10):0.10,C:0.20):0.10,D:0.30):0.10,OUT:0.50);",
    "((((A:0.12,B:0.11):0.09,C:0.22):0.11,D:0.31):0.10,OUT:0.52);",
    "((((A:0.30,C:0.30):0.10,B:0.20):0.10,D:0.30):0.10,OUT:0.50);",
    "((((B:0.40,C:0.40):0.10,A:0.60):0.10,D:0.35):0.10,OUT:0.55);",
    "((((A:0.15,B:0.15):0.12,C:0.25):0.10,D:0.28):0.10,OUT:0.49);",
    "((((B:0.35,C:0.35):0.10,A:0.55):0.10,D:0.33):0.10,OUT:0.51);",
    "((((A:0.20,B:0.20):0.10,C:0.40):0.10,D:0.30):0.10,OUT:0.50);",
    "((((A:0.25,C:0.25):0.10,B:0.45):0.10,D:0.32):0.10,OUT:0.53);",
    "((((A:0.11,B:0.13):0.10,C:0.21):0.11,D:0.29):0.10,OUT:0.50);",
    "((((B:0.42,C:0.44):0.10,A:0.62):0.10,D:0.34):0.10,OUT:0.54);",
    "((((A:0.18,B:0.16):0.10,C:0.26):0.10,D:0.30):0.10,OUT:0.50);",
    "((((A:0.14,B:0.12):0.10,C:0.24):0.10,D:0.27):0.10,OUT:0.48);",
]


@pytest.fixture()
def pipeline_species_tree(tmp_path):
    """Write the shared species tree to a temp file.

    Args:
        tmp_path: Pytest temporary directory fixture.

    Returns:
        The path to the written species tree file.
    """
    path = tmp_path / "species.tree"
    path.write_text(_SPECIES_TREE + "\n")
    return path


@pytest.fixture()
def pipeline_gene_trees(tmp_path):
    """Write the shared gene trees to a temp file.

    Args:
        tmp_path: Pytest temporary directory fixture.

    Returns:
        The path to the written gene trees file.
    """
    path = tmp_path / "genes.tree"
    path.write_text("\n".join(_GENE_TREES) + "\n")
    return path
