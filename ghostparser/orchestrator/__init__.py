"""Streaming introgression orchestrator.

Fuses triplet subtree extraction and per-triplet inference into a single pass.
Owns its tree preprocessing, inference, consolidation, and configuration.
"""

from .runner import run_orchestrator

__all__ = ["run_orchestrator"]
