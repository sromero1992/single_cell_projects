"""
utils.py — TimeSCape utility functions
=======================================
Benjamini-Hochberg correction, acrophase wrapping, tmeta helpers.
"""

from __future__ import annotations

import re
import numpy as np
import pandas as pd


# ── Benjamini-Hochberg multiple testing correction ────────────────────────────

def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.

    Translated directly from the MATLAB bh_adjust_pvalues() helper in
    sce_circ_phase_estimation_stattest.m — results are bit-for-bit identical.

    Parameters
    ----------
    pvalues : array-like of float
        Raw p-values (NaN entries are preserved).

    Returns
    -------
    p_adj : np.ndarray
        BH-adjusted p-values, same shape as input.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    original_shape = pvalues.shape
    pvalues = pvalues.ravel()

    p_adj = np.full_like(pvalues, np.nan)
    finite_mask = np.isfinite(pvalues)
    p_finite = pvalues[finite_mask]

    if p_finite.size == 0:
        return p_adj.reshape(original_shape)

    m = len(p_finite)
    sort_idx = np.argsort(p_finite)
    p_sorted = p_finite[sort_idx]
    ranks = np.arange(1, m + 1)

    p_adj_sorted = np.minimum(1.0, p_sorted * m / ranks)

    # Enforce monotonicity right-to-left
    for k in range(m - 2, -1, -1):
        p_adj_sorted[k] = min(p_adj_sorted[k], p_adj_sorted[k + 1])

    p_adj_finite = np.empty(m)
    p_adj_finite[sort_idx] = p_adj_sorted

    p_adj[finite_mask] = p_adj_finite
    return p_adj.reshape(original_shape)


# ── Acrophase wrapping ────────────────────────────────────────────────────────

def wrap_acrophase(acro: float, period: float = 24.0) -> float:
    """
    Wrap an acrophase value into [0, period).

    Equivalent to the MATLAB post-processing:
        acro_fmt(acro <  0) = acro + 24;
        acro_fmt(acro > 24) = acro - 24;
    """
    if np.isnan(acro):
        return np.nan
    return float(acro % period)


# ── ZT string parsing ─────────────────────────────────────────────────────────

_ZT_PATTERN = re.compile(r"ZT\s*(\d+\.?\d*)", re.IGNORECASE)


def parse_zt_string(s: str) -> float | None:
    """
    Parse a ZT label string to a numeric hour.

    Examples
    --------
    >>> parse_zt_string("ZT06")   # 6.0
    >>> parse_zt_string("ZT03")   # 3.0
    >>> parse_zt_string("zt18")   # 18.0
    >>> parse_zt_string("SampleA") # None  (no ZT pattern)
    """
    m = _ZT_PATTERN.search(str(s))
    return float(m.group(1)) if m else None


def build_tmeta(zt_strings: list[str], exclude: list[str] | None = None) -> pd.DataFrame:
    """
    Build a tmeta DataFrame from a list of ZT label strings.

    Attempts to auto-parse numeric hours from strings like 'ZT00', 'ZT03',
    'zt12', 'ZT 06', etc. Labels that cannot be parsed are dropped (or you
    can pass ``exclude`` to drop specific labels explicitly).

    Parameters
    ----------
    zt_strings : list of str
        Unique values of the ZT metadata column in ``adata.obs``.
    exclude : list of str, optional
        Labels to exclude from the analysis (set their numeric time to NaN
        so they are dropped automatically).

    Returns
    -------
    tmeta : pd.DataFrame
        Columns: ``old_labels`` (str), ``new_labels`` (str),
        ``ZT_times`` (float, numeric hour).

    Examples
    --------
    >>> unique_zt = sorted(adata.obs["ZT_time_str"].unique())
    >>> tmeta = build_tmeta(unique_zt)
    >>> print(tmeta)
      old_labels new_labels  ZT_times
    0       ZT00       ZT00       0.0
    1       ZT03       ZT03       3.0
    ...
    """
    exclude = set(exclude or [])
    rows = []
    for label in zt_strings:
        if label in exclude:
            continue
        numeric = parse_zt_string(label)
        if numeric is None:
            continue
        rows.append({"old_labels": label, "new_labels": label, "ZT_times": numeric})

    if not rows:
        raise ValueError(
            "Could not parse any ZT hours from the provided labels. "
            "Expected labels like 'ZT00', 'ZT03', 'ZT06', ... "
            "Use build_tmeta() with a manual mapping if your labels differ."
        )

    tmeta = pd.DataFrame(rows).sort_values("ZT_times").reset_index(drop=True)
    return tmeta
