"""Shared fixtures, plus the rule that sorts every test into a category.

Four explicit markers -- ``config``, ``output``, ``integration`` and ``parity``
-- are applied in the test files. Anything carrying none of the first three is
the pipeline's statistics and decisions, and is marked ``core`` here so
``pytest -m core`` selects it without every such test naming itself. A test
may carry several: consolidation tests that check a computed average by reading
the TSV it was written to are both ``core`` and ``output``.
"""

import pytest

from tests.fixtures import *

_EXPLICIT_CATEGORIES = {"config", "output", "integration"}


def pytest_collection_modifyitems(items):
    """Add ``core`` to every test that carries none of the explicit categories."""
    for item in items:
        if not any(marker.name in _EXPLICIT_CATEGORIES for marker in item.iter_markers()):
            item.add_marker(pytest.mark.core)
