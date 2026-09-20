"""Species rename map: loading, validation, and the label helpers the outputs use."""

import dendropy
import pytest

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
def test_rename_map_reads_a_two_column_tsv(tmp_path):
    """A TSV maps its first column onto its second, skipping blanks and comments."""
    path = _write(
        tmp_path,
        "names.tsv",
        "# tree label\tdisplay name\nT1\tHomo sapiens\n\nT2\tPan troglodytes\n",
    )
    assert load_species_rename_map(str(path)) == {
        "T1": "Homo sapiens",
        "T2": "Pan troglodytes",
    }


@pytest.mark.config
def test_rename_map_reads_a_yaml_mapping(tmp_path):
    """A YAML mapping is accepted in place of a TSV."""
    path = _write(tmp_path, "names.yaml", "T1: Homo sapiens\nT2: Pan troglodytes\n")
    assert load_species_rename_map(str(path)) == {
        "T1": "Homo sapiens",
        "T2": "Pan troglodytes",
    }


@pytest.mark.config
@pytest.mark.parametrize(
    "name,text,message",
    [
        ("bad.tsv", "T1\tA\textra\n", "two tab-separated columns"),
        ("bad.tsv", "T1\n", "two tab-separated columns"),
        ("dup_key.tsv", "T1\tA\nT1\tB\n", "more than once"),
        ("dup_value.tsv", "T1\tA\nT2\tA\n", "both T1 and T2 to A"),
        ("dup_value.yaml", "T1: Ape\nT2: Ape\n", "both T1 and T2 to Ape"),
        ("bad.yaml", "- T1\n- T2\n", "must be a mapping"),
        ("tab.yaml", "T1: \"Homo\\tsapiens\"\n", "holds a tab"),
        ("comma.tsv", "T1\tHomo, sapiens\n", "holds a comma"),
        ("semicolon.tsv", "T1\tA;B\n", "holds a semicolon"),
        ("equals.tsv", "T1\tA=B\n", "holds an equals sign"),
    ],
)
def test_rename_map_rejects_malformed_files(tmp_path, name, text, message):
    """Malformed maps fail loudly rather than silently renaming nothing."""
    path = _write(tmp_path, name, text)
    with pytest.raises(ValueError, match=message):
        load_species_rename_map(str(path))


@pytest.mark.config
def test_rename_map_rejects_a_missing_file(tmp_path):
    """A missing map is reported by path."""
    with pytest.raises(FileNotFoundError, match="Species rename map not found"):
        load_species_rename_map(str(tmp_path / "absent.tsv"))


@pytest.mark.parametrize(
    "newick,expected",
    [
        # Mapped labels renamed, the unmapped one kept, lengths untouched.
        ("((T1:0.1,T2:0.2):0.3,T3:0.4);", "((Alpha:0.1,'Beta sp.':0.2):0.3,T3:0.4);"),
        # The odd taxon listed first, and a root edge.
        ("(T3:0.4,(T1:0.1,T2:0.2):0.3):0.5;", "(T3:0.4,(Alpha:0.1,'Beta sp.':0.2):0.3):0.5;"),
        # Whole-token matching: labels that merely contain a key are left alone.
        ("((T1:0.1,T10:0.2):0.3,XT1:0.4);", "((Alpha:0.1,T10:0.2):0.3,XT1:0.4);"),
        # A quoted input label is unquoted before lookup and requoted as needed.
        ("(('O''Brien':0.1,T2:0.2):0.3,T3:0.4);", "(('O''Brien':0.1,'Beta sp.':0.2):0.3,T3:0.4);"),
    ],
)
def test_renaming_newick_labels_maps_leaves_and_quotes_as_needed(newick, expected):
    """Leaf labels are mapped in place, quoting display names the format needs.

    Branch lengths follow ``:`` so they are never mistaken for labels, and the
    result must still parse to the display names with its lengths intact.
    """
    rename_map = {"T1": "Alpha", "T2": "Beta sp."}
    renamed = rename_newick_labels(newick, rename_map)
    assert renamed == expected

    before = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    after = dendropy.Tree.get(data=renamed, schema="newick", preserve_underscores=True)
    assert [leaf.taxon.label for leaf in after.leaf_node_iter()] == [
        rename_map.get(leaf.taxon.label, leaf.taxon.label)
        for leaf in before.leaf_node_iter()
    ]
    assert [leaf.edge_length for leaf in after.leaf_node_iter()] == [
        leaf.edge_length for leaf in before.leaf_node_iter()
    ]


def test_renaming_labels_leaves_unmapped_names_alone():
    """The label helper maps a plain name list, and an empty map is a no-op."""
    assert rename_taxon_labels(("T1", "T3"), {"T1": "Alpha"}) == ["Alpha", "T3"]
    assert rename_taxon_labels(("T1",), {}) == ["T1"]
    assert rename_newick_labels("((T1:0.1,T2:0.2):0.3,T3:0.4);", {}) == (
        "((T1:0.1,T2:0.2):0.3,T3:0.4);"
    )
