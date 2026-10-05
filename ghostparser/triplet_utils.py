"""Shared triplet topology utilities, used for preprocessing and inference."""


TOPOLOGY_AB = "((A,B),C)"
TOPOLOGY_BC = "((B,C),A)"
TOPOLOGY_AC = "((A,C),B)"
ALL_TOPOLOGIES = (TOPOLOGY_AB, TOPOLOGY_BC, TOPOLOGY_AC)


def normalize_abc_from_sister_pair(labels, sister_pair):
    """Normalize triplet labels to ``(A, B, C)`` where ``A`` and ``B`` are sisters."""
    labels_set = set(labels)
    a_taxon, b_taxon = sorted(sister_pair)
    c_taxon = next(iter(labels_set - set(sister_pair)))
    return a_taxon, b_taxon, c_taxon
