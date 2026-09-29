"""Tests for the shared command-line wrapper in ``ghostparser.cli_config``.

Every entry point runs through ``run_cli``, which maps how a run ends onto an
exit status and, with ``--debug``, lets the error through with its traceback.
"""

import argparse

import pytest

from ghostparser.cli_config import (
    EXIT_CHECK_FAILED,
    EXIT_CONFIG_ERROR,
    EXIT_INTERNAL_ERROR,
    EXIT_INTERRUPTED,
    EXIT_OK,
    EXIT_RUN_FAILED,
    run_cli,
)
from ghostparser.config import ConfigError, InputError
from ghostparser.orchestrator.trees import OutgroupRootingError


@pytest.mark.parametrize(
    "outcome, expected_status",
    [
        (None, EXIT_OK),
        (EXIT_CHECK_FAILED, EXIT_CHECK_FAILED),
        (ConfigError("Missing required config field: species_tree_path"), EXIT_CONFIG_ERROR),
        (InputError("Tree file not found: genes.tree"), EXIT_RUN_FAILED),
        # A subclass of the input error reports as one.
        (OutgroupRootingError("Could not root the species tree"), EXIT_RUN_FAILED),
        # Anything outside the package's own errors is a bug.
        (TypeError("unsupported operand"), EXIT_INTERNAL_ERROR),
        (KeyboardInterrupt(), EXIT_INTERRUPTED),
    ],
    ids=["success", "check_failed", "config", "input", "rooting", "bug", "interrupt"],
)
def test_run_cli_maps_each_outcome_to_its_exit_status(outcome, expected_status, capsys):
    """Each way a run can end gets its own exit status, a failure reported on
    stderr without a traceback; ``--debug`` lets the error itself through."""

    def run(_args):
        if isinstance(outcome, BaseException):
            raise outcome
        return outcome

    assert run_cli(argparse.ArgumentParser(), run, argv=[]) == expected_status
    if isinstance(outcome, BaseException):
        assert "Traceback" not in capsys.readouterr().err
        with pytest.raises(type(outcome)):
            run_cli(argparse.ArgumentParser(), run, argv=["--debug"])
