"""
pipeline.py — TimeSCape main analysis pipeline
===============================================
run_timescape(): full per-cell-type circadian rhythm detection.

Mirrors sce_circ_phase_estimation_stattest.m (MATLAB v0.2) and
run_timescape.R (TimeSCape_R) in both API and statistical behaviour.

Output files per cell type (6 CSVs + optional heatmap PNG):
    {ct}_period_{T}_circadian_analysis_all.csv
    {ct}_period_{T}_circadian_analysis_confident.csv
    {ct}_period_{T}_circadian_ZTs_mean.csv
    {ct}_period_{T}_circadian_ZTs_mean_normalized.csv
    {ct}_period_{T}_circadian_ZTs_mean_confident.csv
    {ct}_period_{T}_circadian_ZTs_mean_normalized_confident.csv
    summary_results.csv  (cross-cell-type)
"""

from __future__ import annotations

import os
import re
import time
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed

from .core      import estimate_phase_r
from .normalize import normalize_lib_size
from .utils     import bh_adjust, wrap_acrophase, build_tmeta


# ── Internal worker (called by joblib) ───────────────────────────────────────

def _fit_gene(
    gene_idx: int,
    X_ct: np.ndarray,          # cells × genes dense float32
    zt_v: np.ndarray,          # per-cell numeric ZT hours
    actual_times: np.ndarray,  # unique ZT hours present for this cell type
    period: float,
    test_type: str,
) -> dict:
    """Fit one gene and return a flat result dict."""
    Xg_zts = [X_ct[zt_v == z, gene_idx] for z in actual_times]
    return estimate_phase_r(Xg_zts, actual_times, period=period, test_type=test_type)


# ── Public API ───────────────────────────────────────────────────────────────

def run_timescape(
    adata,
    tmeta: pd.DataFrame,
    *,
    celltype_col: str  = "cell_type",
    zt_col: str        = "ZT_time",
    rm_low_conf: bool  = True,
    period: float      = 24.0,
    custom_genelist: list[str] | None = None,
    custom_celltype: list[str] | None = None,
    norm_str: str      = "lib_size",
    plot_heat: bool    = True,
    outdir: str        = ".",
    test_type: str     = "Ftest",
    n_jobs: int        = -1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run the TimeSCape circadian rhythm detection pipeline.

    Parameters
    ----------
    adata : AnnData
        Single-cell object.  ``adata.X`` should contain raw counts
        (when ``norm_str="lib_size"``) or pre-normalized log counts.
    tmeta : pd.DataFrame
        ZT metadata table with columns:
            ``old_labels``  — values that appear in ``adata.obs[zt_col]``
            ``new_labels``  — human-readable ZT label (e.g. "ZT06")
            ``ZT_times``    — numeric ZT hour (float)
        Build with :func:`build_tmeta` or supply manually.
    celltype_col : str
        ``adata.obs`` column holding cell type labels.
    zt_col : str
        ``adata.obs`` column holding ZT time strings (matched to
        ``tmeta.old_labels``).
    rm_low_conf : bool
        Write confident-only output files (default True).
    period : float
        Circadian period in hours: 24 (default) or 12.
    custom_genelist : list of str, optional
        Restrict analysis to these genes.  None = all genes.
    custom_celltype : list of str, optional
        Restrict analysis to these cell types.  None = all cell types.
    norm_str : str
        Normalization strategy: ``"lib_size"`` (default), ``"logcounts"``,
        ``"none"``.
    plot_heat : bool
        Generate a heatmap PNG after each cell type (default True).
    outdir : str
        Root output directory.  A sub-folder is created per cell type.
    test_type : str
        ``"Ftest"`` (default) or ``"LRT"``.
    n_jobs : int
        Number of parallel workers for joblib (default -1 = all CPUs).

    Returns
    -------
    T1 : pd.DataFrame
        Circadian statistics for the last processed cell type.
    T2 : pd.DataFrame
        Per-ZT mean expression for the last processed cell type.

    Notes
    -----
    Six CSV files are written per cell type to ``outdir/<CellType>/``.
    A cross-cell-type summary CSV is written to ``outdir/``.
    """
    t_start = time.time()
    os.makedirs(outdir, exist_ok=True)

    per_label = f"_period_{int(period)}_"

    # ── Validate tmeta ────────────────────────────────────────────────────────
    required_cols = {"old_labels", "new_labels", "ZT_times"}
    missing = required_cols - set(tmeta.columns)
    if missing:
        raise ValueError(f"tmeta is missing columns: {missing}")

    tmeta = tmeta.copy().sort_values("ZT_times").reset_index(drop=True)

    # Map old_labels → numeric ZT hours
    label_to_zt: dict[str, float] = dict(
        zip(tmeta["old_labels"].astype(str), tmeta["ZT_times"].astype(float))
    )
    label_to_new: dict[str, str] = dict(
        zip(tmeta["old_labels"].astype(str), tmeta["new_labels"].astype(str))
    )

    # Add numeric ZT column to obs
    adata.obs["_ts_zt_num"] = (
        adata.obs[zt_col].astype(str).map(label_to_zt)
    )
    adata.obs["_ts_zt_new"] = (
        adata.obs[zt_col].astype(str).map(label_to_new)
    )

    # ── Normalize full object once if using lib_size ──────────────────────────
    if norm_str == "lib_size":
        raw_X = adata.X
        if sp.issparse(raw_X):
            raw_X_dense = raw_X.toarray().astype(np.float32)
        else:
            raw_X_dense = np.asarray(raw_X, dtype=np.float32)
    elif norm_str in ("logcounts", "none"):
        raw_X = adata.X
        if sp.issparse(raw_X):
            raw_X_dense = raw_X.toarray().astype(np.float32)
        else:
            raw_X_dense = np.asarray(raw_X, dtype=np.float32)
    else:
        raise ValueError(f"Unknown norm_str '{norm_str}'")

    gene_names = list(adata.var_names)

    # ── Cell type loop ────────────────────────────────────────────────────────
    all_cell_types = sorted(adata.obs[celltype_col].dropna().unique())
    if custom_celltype:
        all_cell_types = [ct for ct in all_cell_types if ct in custom_celltype]

    nztps_expected = tmeta.shape[0]
    summary_rows   = []

    T1_last = pd.DataFrame()
    T2_last = pd.DataFrame()

    for cell_type in all_cell_types:
        ct_safe = re.sub(r"[^a-zA-Z0-9_]", "_", str(cell_type).strip())
        ct_safe = re.sub(r"_+", "_", ct_safe).strip("_")
        ct_outdir = os.path.join(outdir, ct_safe)
        os.makedirs(ct_outdir, exist_ok=True)

        print(f"\nProcessing cell type: {cell_type}")
        print(f"  Output directory  : {ct_outdir}")

        # ── Subset cells ──────────────────────────────────────────────────────
        ct_mask = adata.obs[celltype_col] == cell_type
        n_cells_ct = int(ct_mask.sum())
        if n_cells_ct == 0:
            print(f"  ⚠ No cells — skipping.")
            continue

        zt_num_ct = adata.obs.loc[ct_mask, "_ts_zt_num"].values.astype(float)

        # ── Per-cell-type normalization ───────────────────────────────────────
        X_ct_raw = raw_X_dense[ct_mask.values, :]          # cells × genes
        if norm_str == "lib_size":
            X_ct = normalize_lib_size(X_ct_raw)
        else:
            X_ct = X_ct_raw                                 # already normalized

        # ── Gene list ─────────────────────────────────────────────────────────
        if custom_genelist:
            gene_idx = [i for i, g in enumerate(gene_names) if g in set(custom_genelist)]
            gl       = [gene_names[i] for i in gene_idx]
            X_ct     = X_ct[:, gene_idx]
        else:
            gl       = gene_names
        num_genes = len(gl)

        # ── Time-point discovery ──────────────────────────────────────────────
        present_times = np.array(sorted(set(zt_num_ct[np.isfinite(zt_num_ct)])))
        nzts = len(present_times)
        print(f"  Time points found : {nzts} / {nztps_expected} expected")

        if nzts < 4:
            print(f"  ⚠ Fewer than 4 time points — skipping {cell_type}.")
            continue
        if nzts < nztps_expected:
            print(f"  ⚠ Missing {nztps_expected - nzts} time point(s) — "
                  f"fitting on available points only.")

        # ── Parallel gene fitting ─────────────────────────────────────────────
        print(f"  Fitting {num_genes:,} genes ...")
        t0 = time.time()

        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_fit_gene)(
                i, X_ct, zt_num_ct, present_times, period, test_type
            )
            for i in range(num_genes)
        )

        print(f"  Fitting complete in {time.time() - t0:.1f} s")

        # ── Assemble result tables ────────────────────────────────────────────
        acro_raw   = np.array([r["acrophase"]     for r in results])
        amp        = np.array([r["amp"]           for r in results])
        mesor      = np.array([r["mesor"]         for r in results])
        pvalue     = np.array([r["p_value"]       for r in results])
        rho        = np.array([r["rho"]           for r in results])
        pvalue_rho = np.array([r["p_value_macro"] for r in results])
        R0_mat     = np.array([r["R0"]            for r in results],
                               dtype=object)    # jagged; pad below

        acro_fmt = np.array([wrap_acrophase(a, period) for a in acro_raw])
        p_adj     = bh_adjust(pvalue)
        p_adj_rho = bh_adjust(pvalue_rho)

        # Build T1
        T1 = pd.DataFrame({
            "Genes"          : gl,
            "Amp"            : amp,
            "Abs_Amp"        : np.abs(amp),
            "Mesor"          : mesor,
            "Acrophase"      : acro_raw,
            "Acrophase_24"   : acro_fmt,
            "Period"         : period,
            "pvalue"         : pvalue,
            "pvalue_adj"     : p_adj,
            "Sine_corr"      : rho,
            "pvalue_corr"    : pvalue_rho,
            "pvalue_adj_corr": p_adj_rho,
        })

        # Build T2: per-ZT means (time columns named by new_labels)
        zt_col_names = []
        R0_cols      = {}
        for j, zt_hr in enumerate(present_times):
            # Find new_label for this numeric hour
            row = tmeta[np.isclose(tmeta["ZT_times"], zt_hr)]
            col_name = row["new_labels"].iloc[0] if len(row) else f"ZT{int(zt_hr):02d}"
            zt_col_names.append(col_name)
            R0_cols[col_name] = np.array([
                r["R0"][j] if len(r["R0"]) > j else np.nan
                for r in results
            ])

        T2 = pd.DataFrame({"Genes": gl, **R0_cols})

        # Drop genes where both p-values are NaN (fit failed completely)
        valid = ~(np.isnan(pvalue) & np.isnan(pvalue_rho))
        T1 = T1[valid].reset_index(drop=True)
        T2 = T2[valid].reset_index(drop=True)

        # ── Confidence ───────────────────────────────────────────────────────
        conf_ftest = T1["pvalue"]     < 0.05
        conf_corr  = T1["pvalue_corr"] < 0.05
        conf_both  = conf_ftest & conf_corr

        n_conf_both  = int(conf_both.sum())
        n_ftest_only = int(conf_ftest.sum())
        n_corr_only  = int(conf_corr.sum())
        n_adj_ftest  = int((T1["pvalue_adj"]      < 0.05).sum())
        n_adj_corr   = int((T1["pvalue_adj_corr"] < 0.05).sum())

        print(f"  Total genes tested  : {len(T1):,}")
        print(f"  Confident (F+corr)  : {n_conf_both:,}")
        print(f"  F-test p<0.05       : {n_ftest_only:,}")
        print(f"  Corr-test p<0.05    : {n_corr_only:,}")

        # ── Sort ──────────────────────────────────────────────────────────────
        T1 = T1.sort_values(
            ["pvalue_adj_corr", "pvalue_adj", "Acrophase_24", "Abs_Amp"],
            ascending=[True, True, True, False],
        ).reset_index(drop=True)
        T2 = T2.set_index("Genes").loc[T1["Genes"]].reset_index()

        # ── T3: ZT-normalized means ───────────────────────────────────────────
        T3 = T2.copy()
        zt_data_cols = [c for c in T3.columns if c != "Genes"]
        ref = T3[zt_data_cols[0]].values.copy()
        zero_mask = ref == 0
        if zero_mask.any() and len(zt_data_cols) > 1:
            ref[zero_mask] = T3[zt_data_cols[1]].values[zero_mask]
        ref[ref == 0] = 1.0
        T3[zt_data_cols] = T3[zt_data_cols].div(ref, axis=0)

        # ── Write all-gene files ──────────────────────────────────────────────
        fbase = os.path.join(ct_outdir, f"{ct_safe}{per_label}")
        T1.to_csv(f"{fbase}circadian_analysis_all.csv",         index=False)
        T2.to_csv(f"{fbase}circadian_ZTs_mean.csv",             index=False)
        T3.to_csv(f"{fbase}circadian_ZTs_mean_normalized.csv",  index=False)

        # ── Write confident-only files ────────────────────────────────────────
        if rm_low_conf:
            conf_mask_sorted = T1["Genes"].isin(T1.loc[conf_both, "Genes"])
            T1c = T1[conf_mask_sorted].reset_index(drop=True)
            T2c = T2.set_index("Genes").loc[T1c["Genes"]].reset_index() if len(T1c) else T2.iloc[0:0].copy()
            T3c = T3.set_index("Genes").loc[T1c["Genes"]].reset_index() if len(T1c) else T3.iloc[0:0].copy()
            T1c.to_csv(f"{fbase}circadian_analysis_confident.csv",        index=False)
            T2c.to_csv(f"{fbase}circadian_ZTs_mean_confident.csv",        index=False)
            T3c.to_csv(f"{fbase}circadian_ZTs_mean_normalized_confident.csv", index=False)

        # ── Optional heatmap ──────────────────────────────────────────────────
        if plot_heat:
            try:
                from .visualize import generate_heatmap
                generate_heatmap(
                    celltype=cell_type,
                    outdir=ct_outdir,
                    period=period,
                    strict=True,
                )
            except Exception as e:
                print(f"  ⚠ Heatmap failed: {e}")

        summary_rows.append({
            "CellType"           : cell_type,
            "NumCells"           : n_cells_ct,
            "NumGenes"           : num_genes,
            "NumTested"          : len(T1),
            "NumConfident_Both"  : n_conf_both,
            "NumNonConfident"    : len(T1) - n_conf_both,
            "NumPvalConf_Ftest"  : n_ftest_only,
            "NumPvalConf_Corr"   : n_corr_only,
            "NumAdjConf_Ftest"   : n_adj_ftest,
            "NumAdjConf_Corr"    : n_adj_corr,
        })

        T1_last, T2_last = T1, T2

    # ── Summary CSV ───────────────────────────────────────────────────────────
    if summary_rows:
        T0 = pd.DataFrame(summary_rows)
        ct_label = custom_celltype[0] if (custom_celltype and len(custom_celltype) == 1) \
            else "all_cell_types"
        ct_label_safe = re.sub(r"[^a-zA-Z0-9_]", "_", ct_label)
        T0.to_csv(os.path.join(outdir, f"{ct_label_safe}{per_label}summary_results.csv"),
                  index=False)

    print(f"\n=== TimeSCape complete.  Total time: {time.time() - t_start:.1f} s ===")
    return T1_last, T2_last
