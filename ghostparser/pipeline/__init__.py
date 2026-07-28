"""Self-contained streaming introgression pipeline.

Fuses triplet subtree extraction and per-triplet inference into a single
streaming pass so the global intermediate triplet-gene-trees structure is never
materialized. This subpackage ports (copies + cleans) the logic it needs from
``tree_parser``/``triplet_processor`` and does not import from them, from
``orchestrator``, so those modules can eventually be deleted.
"""

from __future__ import annotations

from .runner import run_pipeline

__all__ = ["run_pipeline"]
