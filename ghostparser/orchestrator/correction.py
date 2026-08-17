"""Multiple-testing correction shared by the orchestrator's statistical tests.

Two consumers apply the same correction machinery at different scopes:
:mod:`ghostparser.orchestrator.inference` corrects the DCT and KS p-values once
across every triplet in a run, while
:mod:`ghostparser.orchestrator.permutation` corrects the pair of one-tailed
p-values within a single permutation test. Keeping :func:`adjust_p_values` here
lets both import it without a circular dependency.
"""

from statsmodels.stats.multitest import multipletests

from .config import DEFAULT_P_VALUE_CORRECTION, P_VALUE_CORRECTION_CHOICES

__all__ = ["adjust_p_values"]

# GhostParser config names mapped onto the statsmodels ``multipletests`` methods.
_METHOD_MAP = {
    "bfn": "bonferroni",
    "holm": "holm",
    "fdr_bh": "fdr_bh",
    "fdr_by": "fdr_by",
    "fdr_tsbh": "fdr_tsbh",
}


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
