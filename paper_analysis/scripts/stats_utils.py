"""Shared statistics helpers for the paper_analysis figure scripts."""

import numpy as np
import pandas as pd


def bh_fdr_mask(pvals: pd.Series, q: float) -> pd.Series:
    """Benjamini-Hochberg FDR significance mask at threshold q.

    Returns a boolean Series aligned to `pvals`' index: True where the
    cell type survives BH-FDR correction at the given q.
    """
    p_sorted = pvals.sort_values()
    n = len(p_sorted)
    thresh = q * (np.arange(1, n + 1) / n)
    below = p_sorted.values <= thresh
    mask = pd.Series(False, index=pvals.index)
    if below.any():
        k = np.max(np.where(below))
        mask.loc[p_sorted.index[: k + 1]] = True
    return mask
