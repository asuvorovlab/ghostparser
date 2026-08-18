"""Multiple-testing correction shared by the orchestrator's statistical tests.

Lives in its own module so ``inference`` and ``permutation`` can both import it
without a circular dependency.
"""

from statsmodels.stats.multitest import multipletests

from .config import DEFAULT_P_VALUE_CORRECTION, P_VALUE_CORRECTION_CHOICES

__all__ = [
    "INLINE_CORRECTION_CHOICES",
    "MONOTONE_CORRECTION_CHOICES",
    "adjust_p_value_inline",
    "adjust_p_values",
    "is_inline_correction",
    "is_monotone_correction",
]

# GhostParser config names mapped onto the statsmodels ``multipletests`` methods.
_METHOD_MAP = {
    "bfn": "bonferroni",
    "holm": "holm",
    "fdr_bh": "fdr_bh",
    "fdr_by": "fdr_by",
    "fdr_tsbh": "fdr_tsbh",
}

# Corrections that never return an adjusted p-value below the raw one, so
# ``raw > alpha`` already implies ``adjusted > alpha``. That implication is what
# licenses skipping a downstream test once an upstream raw gate has failed; the
# per-method argument is under "Correction inside the bootstrap" in the
# orchestrator guide.
MONOTONE_CORRECTION_CHOICES = ("no", "bfn", "holm", "fdr_bh", "fdr_by")

# Corrections that can be applied to a single p-value knowing only the family
# size, without the other members. Everything else is rank-based and needs the
# whole family in hand.
INLINE_CORRECTION_CHOICES = ("no", "bfn")


def is_monotone_correction(method):
    """Report whether a correction method can only raise p-values.

    Args:
        method: One of :data:`~ghostparser.orchestrator.config.P_VALUE_CORRECTION_CHOICES`.

    Returns:
        ``True`` when the method is in :data:`MONOTONE_CORRECTION_CHOICES`.
    """
    return method in MONOTONE_CORRECTION_CHOICES


def is_inline_correction(method):
    """Report whether a correction method needs only the family size.

    Args:
        method: One of :data:`~ghostparser.orchestrator.config.P_VALUE_CORRECTION_CHOICES`.

    Returns:
        ``True`` when the method is in :data:`INLINE_CORRECTION_CHOICES`.
    """
    return method in INLINE_CORRECTION_CHOICES


def adjust_p_value_inline(p_value, method, family_size):
    """Adjust one p-value from the family size alone.

    Args:
        p_value: The raw p-value.
        method: One of :data:`INLINE_CORRECTION_CHOICES`.
        family_size: Number of tests in the family.

    Returns:
        The adjusted p-value, clipped to 1.0.

    Raises:
        ValueError: If ``method`` cannot be applied without the whole family.
    """
    if method == "no":
        return float(p_value)
    if method == "bfn":
        return min(1.0, float(p_value) * max(1, int(family_size)))
    raise ValueError(
        f"Correction method {method} needs the whole family. "
        f"Choose one of: {', '.join(INLINE_CORRECTION_CHOICES)}"
    )


def adjust_p_values(p_values, method=DEFAULT_P_VALUE_CORRECTION, alpha=0.05):
    """Adjust a family of p-values by the selected correction method.

    Uses the statsmodels ``multipletests`` backend. The ``no`` method returns
    the inputs unchanged as floats.

    Args:
        p_values: Sequence of raw p-values forming one testing family.
        method: One of :data:`~ghostparser.orchestrator.config.P_VALUE_CORRECTION_CHOICES`.
        alpha: FDR level, used by the two-stage BH method and by statsmodels
            when deciding rejections.

    Returns:
        A list of adjusted p-values in the input order.

    Raises:
        ValueError: If ``method`` is not a supported correction method.
    """
    if method not in P_VALUE_CORRECTION_CHOICES:
        raise ValueError(
            f"Unsupported p-value correction method: {method}. "
            f"Choose one of: {', '.join(P_VALUE_CORRECTION_CHOICES)}"
        )

    values = [float(p_value) for p_value in p_values]
    if method == "no" or not values:
        return values

    _, corrected, _, _ = multipletests(
        values,
        alpha=alpha,
        method=_METHOD_MAP[method],
    )
    return [float(p_value) for p_value in corrected]
