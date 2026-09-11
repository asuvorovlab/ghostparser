"""Species rename map: loading, validation, and its effect on the outputs."""

import pytest

from ghostparser.orchestrator.trees import (
    load_species_rename_map,
    read_tree_file,
    rename_taxa_in_tree,
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


def test_renaming_a_tree_touches_only_mapped_terminals(tmp_path):
    """Terminals absent from the map keep their original labels."""
    tree_path = _write(tmp_path, "t.tre", "((T1:0.1,T2:0.2):0.3,T3:0.4);\n")
    tree = read_tree_file(str(tree_path))[0]

    renamed = rename_taxa_in_tree(tree, {"T1": "Alpha", "T2": "Beta"})

    assert renamed == 2
    assert sorted(t.name for t in tree.get_terminals()) == ["Alpha", "Beta", "T3"]


def test_renaming_labels_leaves_unmapped_names_alone():
    """The label helper mirrors the tree helper for plain name lists."""
    assert rename_taxon_labels(("T1", "T3"), {"T1": "Alpha"}) == ["Alpha", "T3"]
    assert rename_taxon_labels(("T1",), {}) == ["T1"]
