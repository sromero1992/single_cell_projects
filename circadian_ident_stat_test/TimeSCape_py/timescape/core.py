"""
core.py — TimeSCape core statistical engine
=============================================
estimate_phase_r(): cosinor fitting + F-test + Pearson correlation.

This is the Python equivalent of estimate_phaseR.m (MATLAB v0.2) and
estimate_phaseR.R (TimeSCape_R).  Statistical results are numerically
equivalent to < 0.1 % across all three platforms.

Key difference vs. original notebook
--------------------------------------
The notebook used ``time_step * np.arange(n_zts)`` as time values, which
assumes evenly-spaced ZT points.  This function requires callers to pass
``actual_times`` — the true numeric ZT hours for each time slot — so that
missing time points (e.g. a cell type absent at ZT12) are handled correctly
without imputation.
"""

from __future__ import annotations

import warnings
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats


def estimate_phase_r(
    Xg_zts: list[np.ndarray],
    actual_times: np.ndarray,
    period: float = 24.0,
    test_type: str = "Ftest",
) -> dict:
    """
    Fit a cosinor model and test significance for a single gene.

    Parameters
    ----------
    Xg_zts : list of np.ndarray
        Each element holds the log-normalised expression values of every cell
        collected at that ZT time point.  Empty arrays are automatically
        removed (robust to missing time points).
    actual_times : np.ndarray
        Numeric ZT hours corresponding to each slot in ``Xg_zts``,
        e.g. ``[0, 3, 6, 9, 15, 18, 21]`` when ZT12 is absent.
    period : float
        Circadian period in hours.  24 (default) or 12 for ultradian.
    test_type : str
        ``"Ftest"`` (default) or ``"LRT"`` (likelihood-ratio test).

    Returns
    -------
    dict with keys:
        acrophase, amp, period, mesor   — cosinor parameters
        p_value                         — F-test or LRT p-value (single-cell)
        rho                             — Pearson correlation (per-ZT means vs fit)
        p_value_macro                   — p-value for that correlation
        R0                              — per-ZT mean expression vector

    All values are ``np.nan`` on fitting failure.
    """
    _nan_result = dict(
        acrophase=np.nan, amp=np.nan, period=period, mesor=np.nan,
        p_value=np.nan, rho=np.nan, p_value_macro=np.nan, R0=[]
    )

    actual_times = np.asarray(actual_times, dtype=float)

    # ── 1. Drop empty time slots ──────────────────────────────────────────────
    non_empty = [i for i, x in enumerate(Xg_zts) if len(x) > 0]
    if len(non_empty) < 4:          # need ≥ 4 time points to fit 3-param cosine
        return _nan_result

    Xg_zts      = [np.asarray(Xg_zts[i], dtype=float) for i in non_empty]
    actual_times = actual_times[non_empty]
    nzts         = len(Xg_zts)

    # ── 2. Flatten cells → arrays ─────────────────────────────────────────────
    icells    = [len(x) for x in Xg_zts]
    num_cells = sum(icells)

    R         = np.concatenate(Xg_zts)
    time_grid = np.repeat(actual_times, icells)   # each cell gets its ZT hour
    R0        = np.array([x.mean() for x in Xg_zts])

    if not np.any(np.isfinite(R)):
        return _nan_result

    # ── 3. Null model ─────────────────────────────────────────────────────────
    mean_model = np.nanmean(R)
    if not np.isfinite(mean_model):
        return _nan_result
    SSR_null = np.nansum((R - mean_model) ** 2)

    # ── 4. Cosinor model ──────────────────────────────────────────────────────
    def cos_model(t, acro, amp, mesor):
        return amp * np.cos(2 * np.pi * (t - acro) / period) + mesor

    max_amp_guess = float(R0.max() - mean_model)
    max_peak_t    = float(actual_times[int(np.argmax(R0))])

    p0     = [max_peak_t, max_amp_guess, mean_model]
    bounds = ([0.0,    -np.inf, 0.0   ],
              [period,  np.inf, np.inf])

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, _ = curve_fit(
                cos_model, time_grid, R,
                p0=p0, bounds=bounds,
                maxfev=10_000,
                method="trf",   # Trust-Region Reflective — matches MATLAB
            )
    except (RuntimeError, ValueError):
        return _nan_result

    acrophase_fit, amp_fit, mesor_fit = popt
    fitted_vals = cos_model(time_grid, *popt)
    SSR_sine    = float(np.sum((R - fitted_vals) ** 2))

    # ── 5. Significance test ──────────────────────────────────────────────────
    p_value = np.nan
    d2 = num_cells - 3

    if test_type == "Ftest":
        if d2 > 0 and np.isfinite(SSR_sine):
            if SSR_null == 0 and SSR_sine == 0:
                p_value = 1.0           # perfectly flat gene
            elif SSR_sine > 0:
                d1     = 2
                F_stat = ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2)
                p_value = float(1.0 - stats.f.cdf(F_stat, d1, d2))

    elif test_type == "LRT":
        n = num_cells
        if n > 0 and SSR_null > 0 and SSR_sine > 0:
            logL_null = -0.5 * n * (np.log(2 * np.pi * SSR_null / n) + 1)
            logL_sine = -0.5 * n * (np.log(2 * np.pi * SSR_sine / n) + 1)
            lrt_stat  = -2.0 * (logL_null - logL_sine)
            p_value   = float(1.0 - stats.chi2.cdf(lrt_stat, df=2))

    # ── Canonical form: positive amplitude ───────────────────────────────────
    # If the optimizer converges with amp < 0 the trough is at acrophase_fit
    # and the true peak is half a period away. Flip both so the reported
    # amplitude is always positive and acrophase is the true peak time.
    if np.isfinite(amp_fit) and amp_fit < 0:
        amp_fit       = -amp_fit
        acrophase_fit = (acrophase_fit + period / 2) % period
    acrophase_fit = float(acrophase_fit % period)   # robust wrap into [0, period)

    # ── 6. Pearson correlation: per-ZT means vs fitted cosine ─────────────────
    fval = cos_model(actual_times, acrophase_fit, amp_fit, mesor_fit)

    rho          = np.nan
    p_value_macro = np.nan
    if np.std(R0) > 1e-12 and np.std(fval) > 1e-12 and nzts >= 4:
        try:
            r, p = stats.pearsonr(R0, fval)
            rho          = float(r)
            p_value_macro = float(p)
        except Exception:
            pass

    return dict(
        acrophase     = float(acrophase_fit),
        amp           = float(amp_fit),
        period        = float(period),
        mesor         = float(mesor_fit),
        p_value       = p_value,
        rho           = rho,
        p_value_macro = p_value_macro,
        R0            = R0.tolist(),
    )
