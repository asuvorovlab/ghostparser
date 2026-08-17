"""Streaming introgression orchestrator.

Fuses triplet subtree extraction and per-triplet inference into a single
streaming pass so the global intermediate triplet-gene-trees structure is never
materialized. The subpackage owns its tree preprocessing, inference, and
configuration; it imports only the shared ``triplet_utils`` topology helpers,
``introgression_mapper`` for consolidation, and the ``config`` trunk.
"""

from .runner import run_orchestrator

__all__ = ["run_orchestrator"]
