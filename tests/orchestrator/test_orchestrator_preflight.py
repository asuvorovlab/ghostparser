"""Tests for the orchestrator's structural preflight data check."""

import pytest

from ghostparser.orchestrator import runner
from ghostparser.orchestrator.preflight import (
    PREFLIGHT_REPORT_FILENAME,
    run_preflight_data_check,
)

# A clean 5-taxon species tree: ingroup A,B,C,D gives C(4,3) = 4 triplets.
_CLEAN_SPECIES = "(((A:1,B:1):1,C:1):1,(D:1,OUT:1):1);\n"

# Gene tree 1 is well formed. Gene tree 2 is a polytomy over A,B,C, so the
# rooted sister pair for triplet A,B,C is undeterminable. Gene tree 3 carries no
# outgroup label at all, so rooting fails outright.
_DIRTY_GENES = (
    "((((A:1,B:1):1,C:1):1,D:1):1,OUT:1);\n"
    "(((A:1,B:1,C:1):1,D:1):1,OUT:1);\n"
    "(((A:1,B:1):1,C:1):1,MISSING:1);\n"
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


def test_clean_inputs_pass_with_no_issues(clean_inputs, tmp_path):
    """Well-formed trees produce an empty issue list and a passing report."""
    species, genes = clean_inputs

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(tmp_path),
    )

    assert result.passed is True
    assert result.issues == []
    # Ingroup A,B,C,D -> C(4,3) = 4 triplets; both gene trees are checked.
    assert result.triplets_checked == 4
    assert result.counters["gene_tree.total_checked"] == 2
    assert result.counters["gene_tree.rooted"] == 2


@pytest.mark.output
def test_report_is_written_only_when_an_output_dir_is_given(clean_inputs, tmp_path):
    """The report lands under the documented filename, or nowhere without a dir."""
    species, genes = clean_inputs
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    written_result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(output_dir),
    )
    written = output_dir / PREFLIGHT_REPORT_FILENAME
    assert written_result.report_path == str(written)
    assert written.read_text() == written_result.report_text

    unwritten_result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=None,
    )
    assert unwritten_result.report_path is None
    assert unwritten_result.report_text == written_result.report_text


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
    # Gene tree 3 has no OUT label.
    assert categories.count("gene_tree.rooting_failed") == 1
    # Gene tree 2 is a polytomy, and only triplet A,B,C is affected by it.
    assert categories.count("triplet.unresolved_rooted_sister_pair") == 1
    # Gene tree 3 never gets rooted, so only trees 1 and 2 reach triplet checks.
    assert result.counters["gene_tree.rooted"] == 2
    assert result.counters["gene_tree.total_checked"] == 3
    # The offending gene-tree index and the failing triplet are both named.
    polytomy_message = next(
        issue.message
        for issue in result.issues
        if issue.category == "triplet.unresolved_rooted_sister_pair"
    )
    assert "Gene tree #2" in polytomy_message
    assert "A,B,C" in polytomy_message
    # Every triplet/gene-tree pair the check looked at is accounted for: 4
    # triplets x 2 rooted trees = 8 pairs, of which tree 2's A,B,C is the only
    # one that cannot be measured.
    assert result.counters["triplet.resolved"] == 7
    assert result.counters.get("triplet.taxa_absent_from_gene_tree", 0) == 0
    assert (
        result.counters["triplet.resolved"]
        + result.counters["triplet.unresolved_rooted_sister_pair"]
        + result.counters.get("triplet.taxa_absent_from_gene_tree", 0)
        == result.triplets_checked * result.counters["gene_tree.rooted"]
        == 8
    )


def test_triplet_filter_entries_are_validated(clean_inputs, tmp_path):
    """Filter lines naming unknown taxa or the outgroup are rejected, not run."""
    species, genes = clean_inputs
    triplet_filter = tmp_path / "triplets.txt"
    triplet_filter.write_text("A,B,C\nA,B,NOPE\nA,B,OUT\n")

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(tmp_path),
        triplet_filter=str(triplet_filter),
    )

    categories = [issue.category for issue in result.issues]
    assert categories.count("triplet_filter.taxa_missing_in_species_tree") == 1
    assert categories.count("triplet_filter.includes_outgroup") == 1
    # Only the one valid line survives to be checked.
    assert result.triplets_checked == 1


def test_species_filter_entries_are_validated(clean_inputs, tmp_path):
    """Filter names absent from the tree or in the outgroup are reported, not combined."""
    species, genes = clean_inputs
    species_filter = tmp_path / "species.txt"
    species_filter.write_text("A,B\nC\nNOPE\nOUT\n")

    result = run_preflight_data_check(
        species_tree_path=str(species),
        gene_trees_path=str(genes),
        outgroups=["OUT"],
        output_dir=str(tmp_path),
        species_filter=str(species_filter),
    )

    categories = [issue.category for issue in result.issues]
    assert categories.count("species_filter.taxa_missing_in_species_tree") == 1
    assert categories.count("species_filter.includes_outgroup") == 1
    # The three usable species form one triplet, out of the tree's four.
    assert result.triplets_checked == 1


@pytest.mark.parametrize(
    "outgroups, expected_message",
    [
        ([], "No outgroup taxa were provided"),
        (["NOT_PRESENT"], "Could not root species tree"),
    ],
)
def test_impossible_checks_raise(clean_inputs, tmp_path, outgroups, expected_message):
    """Conditions that make the check itself impossible raise ValueError."""
    species, genes = clean_inputs

    with pytest.raises(ValueError, match=expected_message):
        run_preflight_data_check(
            species_tree_path=str(species),
            gene_trees_path=str(genes),
            outgroups=outgroups,
            output_dir=str(tmp_path),
        )


def test_multi_tree_species_file_raises(tmp_path):
    """A species-tree file holding more than one tree is rejected."""
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text(_CLEAN_SPECIES + _CLEAN_SPECIES)
    genes.write_text(_CLEAN_GENES)

    with pytest.raises(ValueError, match="exactly one tree"):
        run_preflight_data_check(
            species_tree_path=str(species),
            gene_trees_path=str(genes),
            outgroups=["OUT"],
            output_dir=str(tmp_path),
        )


@pytest.mark.integration
@pytest.mark.output
def test_runner_preflight_mode_skips_analysis(dirty_inputs, tmp_path):
    """The flag short-circuits the run: only the report is produced.

    The cap is set below the four ingroup triplets so its effect is visible in
    the result, proving the config value reaches the check rather than the
    module default.
    """
    species, genes = dirty_inputs
    output_dir = tmp_path / "results"

    result = runner.run_orchestrator(
        {
            "species_tree": str(species),
            "gene_trees": str(genes),
            "outgroup": "OUT",
            "output": str(output_dir),
            "overwrite": True,
            "triplet_filter": None,
            "species_filter": None,
            "preflight_data_check": True,
            "preflight_triplet_cap": 3,
        }
    )

    assert result.passed is False
    assert result.triplets_checked == 3
    assert any(
        issue.category == "analysis.triplet_cap_applied" for issue in result.issues
    )
    written = sorted(path.name for path in output_dir.iterdir())
    assert written == [PREFLIGHT_REPORT_FILENAME]


@pytest.mark.integration
def test_runner_returns_none_when_preflight_cannot_run(tmp_path):
    """An impossible check is reported without raising out of the runner."""
    species = tmp_path / "species.tree"
    genes = tmp_path / "genes.tre"
    species.write_text(_CLEAN_SPECIES)
    genes.write_text(_CLEAN_GENES)

    result = runner.run_orchestrator(
        {
            "species_tree": str(species),
            "gene_trees": str(genes),
            "outgroup": "NOT_PRESENT",
            "output": str(tmp_path / "results"),
            "overwrite": True,
            "triplet_filter": None,
            "species_filter": None,
            "preflight_data_check": True,
            "preflight_triplet_cap": 15000,
        }
    )

    assert result is None
