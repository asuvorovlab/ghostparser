"""Species rename map: loading, validation, and the label helpers the outputs use."""

from io import StringIO

import pytest
from Bio import Phylo

from ghostparser.orchestrator.trees import (
    load_species_rename_map,
    rename_newick_labels,
    rename_taxon_labels,
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.mark.config
@pytest.mark.parametrize(
    "name, text",
    [
        ("names.tsv", "# tree label\tdisplay name\nT1\tHomo sapiens\n\nT2\tPan troglodytes\n"),
        ("names.yaml", "T1: Homo sapiens\nT2: Pan troglodytes\n"),
    ],
    ids=["tsv", "yaml"],
)
def test_rename_map_reads_a_tsv_or_a_yaml_mapping(tmp_path, name, text):
    """A TSV maps its first column onto its second, skipping blanks and comments; YAML maps directly."""
    path = _write(tmp_path, name, text)
    assert load_species_rename_map(str(path)) == {
        "T1": "Homo sapiens",
        "T2": "Pan troglodytes",
    }


@pytest.mark.config
@pytest.mark.parametrize(
    "name, text, message",
    [
        ("bad.tsv", "T1\tA\textra\n", "two tab-separated columns"),
        ("bad.tsv", "T1\n", "two tab-separated columns"),
        ("dup_key.tsv", "T1\tA\nT1\tB\n", "more than once"),
        ("dup_value.tsv", "T1\tA\nT2\tA\n", "both T1 and T2 to A"),
        ("dup_value.yaml", "T1: Ape\nT2: Ape\n", "both T1 and T2 to Ape"),
        ("bad.yaml", "- T1\n- T2\n", "must be a mapping"),
        ("tab.yaml", "T1: \"Homo\\tsapiens\"\n", "holds a tab"),
        ("comma.tsv", "T1\tHomo, sapiens\n", "holds a comma"),
        ("absent.tsv", None, "Species rename map not found"),
    ],
)
def test_rename_map_rejects_malformed_or_missing_files(tmp_path, name, text, message):
    """Malformed or missing maps fail loudly rather than silently renaming nothing."""
    path = tmp_path / name
    if text is not None:
        path.write_text(text)
    with pytest.raises((ValueError, FileNotFoundError), match=message):
        load_species_rename_map(str(path))


_RENAME_MAP = {"T1": "Alpha", "T2": "Beta sp."}


@pytest.mark.parametrize(
    "newick, rename_map, expected",
    [
        # Mapped labels renamed, the unmapped one kept, lengths untouched.
        ("((T1:0.1,T2:0.2):0.3,T3:0.4);", _RENAME_MAP, "((Alpha:0.1,'Beta sp.':0.2):0.3,T3:0.4);"),
        # The odd taxon listed first, and a root edge.
        ("(T3:0.4,(T1:0.1,T2:0.2):0.3):0.5;", _RENAME_MAP, "(T3:0.4,(Alpha:0.1,'Beta sp.':0.2):0.3):0.5;"),
        # Whole-token matching: labels that merely contain a key are left alone.
        ("((T1:0.1,T10:0.2):0.3,XT1:0.4);", _RENAME_MAP, "((Alpha:0.1,T10:0.2):0.3,XT1:0.4);"),
        # A quoted input label is unquoted before lookup and requoted as needed.
        ("(('O''Brien':0.1,T2:0.2):0.3,T3:0.4);", _RENAME_MAP, "(('O''Brien':0.1,'Beta sp.':0.2):0.3,T3:0.4);"),
        # An empty map is a no-op.
        ("((T1:0.1,T2:0.2):0.3,T3:0.4);", {}, "((T1:0.1,T2:0.2):0.3,T3:0.4);"),
    ],
    ids=["mapped", "odd_first", "whole_token", "quoted_input", "empty_map"],
)
def test_renaming_labels_maps_leaves_and_quotes_as_needed(newick, rename_map, expected):
    """Leaf labels are mapped in place, quoting display names the format needs.

    Branch lengths follow ``:`` so they are never mistaken for labels, and the
    result must still parse to the display names with its lengths intact. The
    plain label helper maps the same leaves to the same names.
    """
    renamed = rename_newick_labels(newick, rename_map)
    assert renamed == expected

    before = Phylo.read(StringIO(newick), "newick").get_terminals()
    after = Phylo.read(StringIO(renamed), "newick").get_terminals()
    before_labels = [leaf.name for leaf in before]
    after_labels = [leaf.name for leaf in after]
    assert after_labels == [rename_map.get(label, label) for label in before_labels]
    assert rename_taxon_labels(before_labels, rename_map) == after_labels
    assert [leaf.branch_length for leaf in after] == [leaf.branch_length for leaf in before]
