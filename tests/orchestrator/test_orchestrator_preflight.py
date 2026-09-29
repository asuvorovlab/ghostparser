"""Tests for the orchestrator's structural preflight data check."""

import pytest

from ghostparser.orchestrator import runner
from ghostparser.orchestrator.preflight import (
    PREFLIGHT_REPORT_FILENAME,
    run_preflight_data_check,
)

# A clean 5-taxon species tree: ingroup A,B,C,D gives C(4,3) = 4 triplets.
_CLEAN_SPECIES = "(((A,B),C),(D,OUT));\n"

# Gene tree 1 is well formed. Gene tree 2 is a polytomy over A,B,C, so the
# rooted sister pair for triplet A,B,C is undeterminable. Gene tree 3 carries no
# outgroup label at all, so rooting fails outright. Gene tree 4 has no branch
# lengths, which is counted, not reported as a defect.
_DIRTY_GENES = (
    "((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);\n"
    "(((A:1,B:1,C:1):1,D:1):1,OUT:1);\n"
    "(((A:1,B:1):1,C:1):1,MISSING:1);\n"
    "((((A,B),C),D),OUT);\n"
)

_CLEAN_GENES = (
    "((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);\n"
    "((((A:1,C:1):1,B:1):1,D:1):1,OUT:1);\n"
)


@pytest.fixture()
def clean_inputs(tmp_path):
    """Write a species tree and gene trees with no structural defects."""
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text(_CLEAN_SPECIES)
    genes.write_text(_CLEAN_GENES)
    return species, genes


@pytest.fixture()
def dirty_inputs(tmp_path):
    """Write gene trees carrying one polytomy and one missing-outgroup defect."""
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text(_CLEAN_SPECIES)
    genes.write_text(_DIRTY_GENES)
    return species, genes


@pytest.mark.output
def test_clean_inputs_pass_and_the_report_lands_where_documented(clean_inputs, tmp_path):
    """Well-formed trees give an empty issue list, and the report is written only on request.

    With an output directory the report lands under the documented filename
    and equals the returned text; without one nothing is written and the same
    text comes back, so the check is usable without touching disk.
    """
    species, genes = clean_inputs
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(output_dir),
    )

    assert result.passed is True
    assert result.issues == []
    # Ingroup A,B,C,D -> C(4,3) = 4 triplets; both gene trees are checked.
    assert result.triplets_checked == 4
    assert result.counters["gene_tree.total_checked"] == 2
    assert result.counters["gene_tree.rooted"] == 2
    written = output_dir / PREFLIGHT_REPORT_FILENAME
    assert result.report_path == str(written)
    assert written.read_text() == result.report_text

    unwritten = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=None,
    )
    assert unwritten.report_path is None
    assert unwritten.report_text == result.report_text


def test_species_tree_is_rooted_where_the_outgroups_branch_off(tmp_path):
    """The check roots both trees as the run does, whatever the file's orientation.

    OUT1 and OUT2 are not a clade as written (they sit on either side of the
    file's root), but both branch off the ingroup at one node, so the ingroup
    is A,B,C,D and 4 triplets are checked. Each gene tree roots from its
    farthest outgroup, read from its own branch lengths, and the counters
    say which outgroup was the farthest, which were used, and which were
    tangled and pruned unused.
    """
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text("(OUT1:1,(OUT2:1,(((A:1,B:1):1,C:1):1,D:1):1):1);\n")
    genes.write_text(
        # OUT1 only; OUT1 with OUT2 nested among the ingroup; OUT2 only.
        "((((A:1,B:1):1,C:1):1,D:1):1,OUT1:1);\n"
        "(((A:1,B:1):1,(C:1,OUT2:1):1):1,OUT1:1);\n"
        "((((A:1,C:1):1,B:1):1,D:1):1,OUT2:1);\n"
    )

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT2", "OUT1"],
        output_dir=str(tmp_path),
    )

    assert result.passed is True
    assert result.triplets_checked == 4
    assert result.counters["gene_tree.rooted"] == 3
    # In tree 2 OUT1's mean path to A, B, C is 4 against OUT2's 10/3, so
    # OUT1 is the farthest and roots it; OUT2, C's sister, is tangled and
    # pruned unused, counted, not reported as a defect.
    assert result.counters["gene_tree.farthest.OUT1"] == 2
    assert result.counters["gene_tree.farthest.OUT2"] == 1
    assert result.counters["gene_tree.rooted_on.OUT1"] == 2
    assert result.counters["gene_tree.rooted_on.OUT2"] == 1
    assert result.counters["gene_tree.tangled.OUT2"] == 1
    assert result.counters["gene_tree.tangled_trees"] == 1
    # Tree 2 lacks D, so its three D triplets are skipped, not failed.
    assert result.counters["triplet.resolved"] == 9
    assert result.counters["triplet.taxa_absent_from_gene_tree"] == 3


def test_detects_polytomy_and_missing_outgroup(dirty_inputs, tmp_path):
    """Each planted defect is reported under its own category exactly once."""
    species, genes = dirty_inputs

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(tmp_path),
    )

    categories = [issue.category for issue in result.issues]
    assert result.passed is False
    # Gene tree 3 has no OUT label; gene tree 4, lacking branch lengths, is
    # counted but raises no issue.
    assert categories.count("gene_tree.rooting_failed") == 1
    assert result.counters["gene_tree.missing_branch_lengths"] == 1
    # Gene tree 2 is a polytomy, and only triplet A,B,C is affected by it.
    assert categories.count("triplet.unresolved_rooted_sister_pair") == 1
    # Gene tree 3 never gets rooted, so trees 1, 2 and 4 reach triplet checks.
    assert result.counters["gene_tree.rooted"] == 3
    assert result.counters["gene_tree.total_checked"] == 4
    # The offending gene-tree index and the failing triplet are both named.
    polytomy_message = next(
        issue.message
        for issue in result.issues
        if issue.category == "triplet.unresolved_rooted_sister_pair"
    )
    assert "Gene tree #2" in polytomy_message
    assert "A,B,C" in polytomy_message
    # Every triplet/gene-tree pair the check looked at is accounted for: 4
    # triplets x 3 rooted trees = 12 pairs, of which tree 2's A,B,C is the
    # only one that cannot be measured.
    assert result.counters["triplet.resolved"] == 11
    assert result.counters.get("triplet.taxa_absent_from_gene_tree", 0) == 0
    assert (
        result.counters["triplet.resolved"]
        + result.counters["triplet.unresolved_rooted_sister_pair"]
        + result.counters.get("triplet.taxa_absent_from_gene_tree", 0)
        == result.triplets_checked * result.counters["gene_tree.rooted"]
        == 12
    )


@pytest.mark.parametrize(
    "kind, text",
    [
        # One valid line, one naming an unknown taxon, one naming the outgroup.
        ("triplet_filter", "A,B,C\nA,B,NOPE\nA,B,OUT\n"),
        # Three usable species, an unknown name and the outgroup.
        ("species_filter", "A,B\nC\nNOPE\nOUT\n"),
    ],
)
def test_filter_entries_are_validated(clean_inputs, tmp_path, kind, text):
    """Filter entries naming unknown taxa or the outgroup are reported, not run.

    Each is reported under its own category and the rest are checked: the one
    valid triplet line, or the one triplet the three usable species form out
    of the tree's four.
    """
    species, genes = clean_inputs
    filter_path = tmp_path / f"{kind}.txt"
    filter_path.write_text(text)

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(tmp_path),
        **{kind: str(filter_path)},
    )

    categories = [issue.category for issue in result.issues]
    assert categories.count(f"{kind}.taxa_missing_in_species_tree") == 1
    assert categories.count(f"{kind}.includes_outgroup") == 1
    assert result.triplets_checked == 1


@pytest.mark.parametrize(
    "species_text, outgroups, expected_message",
    [
        (_CLEAN_SPECIES, [], "No outgroup taxa were provided"),
        (_CLEAN_SPECIES, ["NOT_PRESENT"], "none of the outgroup taxa"),
        (_CLEAN_SPECIES, ["A", "B", "C", "D", "OUT"], "every taxon"),
        # OUT and C branch off at different points: A,B and D end up on
        # different sides of the outgroups, and the message names both groups.
        (_CLEAN_SPECIES, ["OUT", "C"], "1 taxon: D"),
        (_CLEAN_SPECIES + _CLEAN_SPECIES, ["OUT"], "exactly one tree"),
    ],
    ids=["no_outgroups", "outgroup_absent", "all_outgroups", "outgroups_tangled", "two_trees"],
)
def test_impossible_checks_raise(tmp_path, species_text, outgroups, expected_message):
    """Conditions that make the check itself impossible raise ValueError."""
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text(species_text)
    genes.write_text(_CLEAN_GENES)

    with pytest.raises(ValueError, match=expected_message):
        run_preflight_data_check(
            species_tree_path=str(species),
            gene_trees_path=str(genes),
            outgroups=outgroups,
            output_dir=str(tmp_path),
        )


@pytest.mark.integration
@pytest.mark.output
@pytest.mark.parametrize("outgroup", ["OUT", "NOT_PRESENT"], ids=["runs", "cannot_run"])
def test_runner_preflight_mode_skips_analysis(dirty_inputs, tmp_path, outgroup):
    """The flag short-circuits the run: only the report is produced, or nothing.

    The cap is set below the four ingroup triplets so its effect is visible in
    the result, proving the config value reaches the check rather than the
    module default. An impossible check (an outgroup the species tree lacks)
    is reported without raising out of the runner, which returns ``None``
    and writes no report.
    """
    species, genes = dirty_inputs
    output_dir = tmp_path / "results"

    result = runner.run_orchestrator(
        {
            "species_tree": str(species),
            "gene_trees": str(genes),
            "outgroup": outgroup,
            "output": str(output_dir),
            "overwrite": True,
            "triplet_filter": None,
            "species_filter": None,
            "preflight_data_check": True,
            "preflight_triplet_cap": 3,
        }
    )

    if outgroup == "NOT_PRESENT":
        assert result is None
        assert list(output_dir.iterdir()) == []
        return
    assert result.passed is False
    assert result.triplets_checked == 3
    assert any(
        issue.category == "analysis.triplet_cap_applied" for issue in result.issues
    )
    assert sorted(path.name for path in output_dir.iterdir()) == [PREFLIGHT_REPORT_FILENAME]
