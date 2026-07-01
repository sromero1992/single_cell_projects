"""
visualize.py — TimeSCape visualization functions
=================================================
generate_heatmap()  — Z-score heatmap of confident circadian genes.
plot_gene_single()  — cosine fit + per-ZT means for one gene.
save_batch_plots()  — batch PNG export of gene plots.
"""

from __future__ import annotations

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
# Only force the non-interactive backend when no GUI backend is already active.
# This allows app_gui.py (TkAgg) to import this module without conflict.
if matplotlib.get_backend().lower() in ("agg", ""):
    try:
        matplotlib.use("Agg")
    except Exception:
        pass
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import zscore

# Classic circadian clock genes (case-insensitive match)
_CLOCK_GENES = {
    "arntl", "bmal1", "clock", "cry1", "cry2",
    "per1", "per2", "per3", "nr1d1", "nr1d2",
    "rev-erba", "rev-erbb", "dbp", "tef", "hlf",
    "rora", "rorb", "rorc", "timeless", "npas2",
    "ciart", "nfil3", "dec1", "dec2",
}


# ── Heatmap ───────────────────────────────────────────────────────────────────

def generate_heatmap(
    celltype: str,
    outdir: str = ".",
    period: float = 24.0,
    strict: bool = True,
    circ_only: bool = False,
    custom_name: str = "",
    figsize: tuple = (10, 8),
    cmap: str = "RdBu_r",
    max_genes: int = 100,
) -> str | None:
    """
    Generate a Z-score heatmap from TimeSCape output CSVs.

    Reads ``*circadian_analysis_confident.csv`` and
    ``*circadian_ZTs_mean_confident.csv`` from ``outdir`` and plots a
    heatmap with genes ordered by acrophase.

    Parameters
    ----------
    celltype : str
        Cell type name (used to locate the CSV files in ``outdir``).
    outdir : str
        Directory containing the TimeSCape output CSVs for this cell type.
    period : float
        Period used in the analysis (24 or 12).
    strict : bool
        If True, additionally filter to BH-adjusted p < 0.05 for both tests.
    circ_only : bool
        If True, restrict to classic circadian clock genes only.
    custom_name : str
        Suffix added to the output PNG filename.
    figsize : tuple
        Figure size.
    cmap : str
        Matplotlib colormap.
    max_genes : int
        Maximum number of genes to show (top ranked).

    Returns
    -------
    str | None
        Path to the saved PNG, or None if no genes passed the filters.
    """
    ct_safe = re.sub(r"[^a-zA-Z0-9_]", "_", str(celltype).strip())
    ct_safe = re.sub(r"_+", "_", ct_safe).strip("_")
    per_label = f"_period_{int(period)}_"
    fbase = os.path.join(outdir, f"{ct_safe}{per_label}")

    stats_file = f"{fbase}circadian_analysis_confident.csv"
    means_file = f"{fbase}circadian_ZTs_mean_confident.csv"

    if not os.path.isfile(stats_file) or not os.path.isfile(means_file):
        # Fall back to all-genes file
        stats_file = f"{fbase}circadian_analysis_all.csv"
        means_file = f"{fbase}circadian_ZTs_mean.csv"

    if not os.path.isfile(stats_file):
        warnings.warn(f"Stats file not found: {stats_file}")
        return None

    T1 = pd.read_csv(stats_file)
    T2 = pd.read_csv(means_file)

    # Additional filters
    if strict and "pvalue_adj_corr" in T1.columns:
        T1 = T1[
            (T1["pvalue_adj_corr"] < 0.05) & (T1["pvalue_adj"] < 0.05)
        ].reset_index(drop=True)

    if circ_only:
        T1 = T1[T1["Genes"].str.lower().isin(_CLOCK_GENES)].reset_index(drop=True)

    if len(T1) == 0:
        print(f"  No genes passed filters for {celltype} — heatmap skipped.")
        return None

    # Clamp to max_genes (already sorted by p-value in pipeline)
    T1 = T1.iloc[:max_genes]
    T2 = T2.set_index("Genes").loc[T1["Genes"]].copy()

    zt_cols = [c for c in T2.columns if c != "Genes"]
    data = T2[zt_cols].values.astype(float)

    # Z-score each gene (row)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data_z = np.apply_along_axis(
            lambda x: zscore(x, nan_policy="omit") if np.nanstd(x) > 0 else x,
            axis=1, arr=data
        )

    # Sort by acrophase
    acro_vals = T1["Acrophase_24"].values
    sort_order = np.argsort(acro_vals)
    data_z = data_z[sort_order, :]
    gene_labels = T1["Genes"].values[sort_order]

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data_z, aspect="auto", cmap=cmap, vmin=-2, vmax=2)

    ax.set_xticks(range(len(zt_cols)))
    ax.set_xticklabels(zt_cols, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(gene_labels)))
    ax.set_yticklabels(gene_labels, fontsize=max(5, 9 - len(gene_labels) // 20))
    ax.set_xlabel("Zeitgeber Time", fontsize=11)
    ax.set_ylabel("Gene", fontsize=11)

    title = f"{celltype} — Circadian Genes (n={len(gene_labels)})"
    if strict:
        title += " [BH-adjusted p<0.05]"
    if circ_only:
        title += " [clock genes only]"
    ax.set_title(title, fontsize=12)

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=9)

    plt.tight_layout()

    suffix = f"_{custom_name}" if custom_name else ""
    png_path = f"{fbase}heatmap{'_strict' if strict else ''}{'_circ' if circ_only else ''}{suffix}.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Heatmap saved: {png_path}")
    return png_path


# ── Single-gene plot ──────────────────────────────────────────────────────────

def plot_gene_single(
    adata,
    gene: str,
    celltype: str,
    tmeta: pd.DataFrame,
    *,
    celltype_col: str = "cell_type",
    zt_col: str       = "ZT_time",
    period: float     = 24.0,
    norm_str: str     = "lib_size",
    show_cells: bool  = True,
    use_violin: bool  = False,
    figsize: tuple    = (9, 5),
    save_path: str | None = None,
    ax=None,
) -> "matplotlib.figure.Figure":
    """
    Plot cosine fit + per-ZT means for a single gene in a single cell type.

    Reproduces the three-layer visualization from plot_gene.R and the MATLAB
    sce_circ_plot_gene.m:
      * Blue line  — fitted cosine curve (evaluated densely over [0, period])
      * Orange circles + line — per-ZT-point mean expression
      * Red dashed vertical — acrophase
      * Violin / strip overlay (optional) — single-cell distribution per ZT

    Parameters
    ----------
    adata : AnnData
        Full or subsetted object containing raw or normalized data.
    gene : str
        Gene to plot.
    celltype : str
        Cell type to subset.
    tmeta : pd.DataFrame
        ZT metadata (old_labels, new_labels, ZT_times).
    celltype_col, zt_col : str
        Metadata column names.
    period : float
        Period (24 or 12).
    norm_str : str
        Normalization applied before plotting.
    show_cells : bool
        Overlay single-cell data points.
    use_violin : bool
        Use violins instead of strip (scatter) for single-cell overlay.
    figsize : tuple
        Figure size.
    save_path : str, optional
        If given, save the figure to this path.
    ax : matplotlib Axes, optional
        Axes to plot into.  If None, a new figure is created.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import scipy.sparse as sp
    from .core      import estimate_phase_r
    from .normalize import normalize_lib_size

    # ── Data extraction ───────────────────────────────────────────────────────
    ct_mask = adata.obs[celltype_col] == celltype
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in adata.var_names.")
    g_idx = list(adata.var_names).index(gene)

    X_ct = adata[ct_mask, :].X
    if sp.issparse(X_ct):
        X_ct = X_ct.toarray()
    X_ct = np.asarray(X_ct, dtype=np.float32)

    if norm_str == "lib_size":
        X_ct = normalize_lib_size(X_ct)

    x_gene = X_ct[:, g_idx]

    label_to_zt = dict(zip(tmeta["old_labels"].astype(str),
                           tmeta["ZT_times"].astype(float)))
    label_to_new = dict(zip(tmeta["old_labels"].astype(str),
                            tmeta["new_labels"].astype(str)))

    zt_num = adata.obs.loc[ct_mask, zt_col].astype(str).map(label_to_zt).values.astype(float)
    actual_times = np.array(sorted(set(zt_num[np.isfinite(zt_num)])))
    nzts = len(actual_times)

    if nzts < 4:
        raise ValueError(f"Only {nzts} time points for {celltype} — cannot fit.")

    Xg_zts = [x_gene[zt_num == z] for z in actual_times]
    fit = estimate_phase_r(Xg_zts, actual_times, period=period)

    R0 = np.array([x.mean() for x in Xg_zts])
    zt_labels = [label_to_new.get(
        next((k for k, v in label_to_zt.items() if np.isclose(v, t)), ""),
        f"ZT{int(t):02d}"
    ) for t in actual_times]

    # ── Plot ──────────────────────────────────────────────────────────────────
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Single-cell overlay
    if show_cells:
        if use_violin:
            parts = ax.violinplot(
                [x_gene[zt_num == z] for z in actual_times],
                positions=actual_times,
                widths=1.5, showmedians=True,
            )
            for pc in parts["bodies"]:
                pc.set_facecolor("#AACCEE")
                pc.set_alpha(0.4)
        else:
            jitter = np.random.default_rng(42).uniform(-0.5, 0.5, size=len(x_gene))
            ax.scatter(zt_num + jitter, x_gene,
                       c="#AACCEE", alpha=0.25, s=8, linewidths=0, zorder=1)

    # Per-ZT means
    ax.plot(actual_times, R0, "o-", color="#FF8C00", lw=2, ms=7,
            zorder=3, label="Per-ZT mean")

    # Cosine fit
    if not np.isnan(fit["acrophase"]):
        t_dense = np.linspace(0, period, 500)
        y_dense = (fit["amp"] * np.cos(2 * np.pi * (t_dense - fit["acrophase"]) / period)
                   + fit["mesor"])
        ax.plot(t_dense, y_dense, "-", color="#1F6DB5", lw=2.5,
                zorder=4, label="Cosine fit")

        # Acrophase line
        ax.axvline(x=fit["acrophase"], color="red", linestyle="--",
                   lw=1.5, alpha=0.8, zorder=5,
                   label=f"Acrophase: ZT{fit['acrophase']:.1f}")

    ax.set_xticks(actual_times)
    ax.set_xticklabels(zt_labels, rotation=45, ha="right")
    ax.set_xlabel("Zeitgeber Time", fontsize=11)
    ax.set_ylabel("log-normalized expression", fontsize=11)

    p_str  = f"p={fit['p_value']:.2e}"  if np.isfinite(fit['p_value'])  else "p=NA"
    pc_str = f"p_corr={fit['p_value_macro']:.2e}" if np.isfinite(fit['p_value_macro']) else "p_corr=NA"
    ax.set_title(
        f"{gene}  [{celltype}]\n"
        f"Acrophase: ZT{fit['acrophase']:.1f}  Amp: {fit['amp']:.3f}  {p_str}  {pc_str}",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")

    return fig


# ── Batch gene plots ──────────────────────────────────────────────────────────

def save_batch_plots(
    adata,
    tmeta: pd.DataFrame,
    celltype: str,
    outdir: str,
    *,
    celltype_col: str = "cell_type",
    zt_col: str       = "ZT_time",
    period: float     = 24.0,
    norm_str: str     = "lib_size",
    plot_type: int    = 1,
    max_genes: int    = 200,
    show_cells: bool  = True,
    use_violin: bool  = False,
    n_cols: int       = 4,
) -> None:
    """
    Export individual PNG plots for a set of genes.

    Parameters
    ----------
    adata : AnnData
    tmeta : pd.DataFrame
    celltype : str
    outdir : str
        Root TimeSCape output directory (must contain the cell-type sub-folder).
    celltype_col, zt_col : str
    period : float
    norm_str : str
    plot_type : int
        1 = confident genes only (default)
        2 = non-confident genes
        3 = classic circadian clock genes
    max_genes : int
        Maximum number of plots to generate.
    show_cells : bool
        Overlay single-cell points.
    use_violin : bool
        Use violin overlay instead of scatter.
    n_cols : int
        Columns per panel figure.
    """
    ct_safe   = re.sub(r"[^a-zA-Z0-9_]", "_", str(celltype).strip()).strip("_")
    per_label = f"_period_{int(period)}_"
    ct_outdir = os.path.join(outdir, ct_safe)
    fbase     = os.path.join(ct_outdir, f"{ct_safe}{per_label}")

    if plot_type == 1:
        stats_csv = f"{fbase}circadian_analysis_confident.csv"
        subdir    = os.path.join(ct_outdir, "plots_confident")
    elif plot_type == 2:
        stats_csv = f"{fbase}circadian_analysis_all.csv"
        subdir    = os.path.join(ct_outdir, "plots_non_confident")
    else:
        stats_csv = f"{fbase}circadian_analysis_confident.csv"
        subdir    = os.path.join(ct_outdir, "plots_clock_genes")

    if not os.path.isfile(stats_csv):
        stats_csv = f"{fbase}circadian_analysis_all.csv"

    if not os.path.isfile(stats_csv):
        warnings.warn(f"Stats CSV not found: {stats_csv}")
        return

    T1 = pd.read_csv(stats_csv)
    if plot_type == 2:
        T1 = T1[~((T1["pvalue"] < 0.05) & (T1["pvalue_corr"] < 0.05))]
    elif plot_type == 3:
        T1 = T1[T1["Genes"].str.lower().isin(_CLOCK_GENES)]

    T1 = T1.iloc[:max_genes].reset_index(drop=True)
    if len(T1) == 0:
        print(f"  No genes to plot for {celltype} (type={plot_type}).")
        return

    os.makedirs(subdir, exist_ok=True)
    print(f"  Saving {len(T1)} gene plots to {subdir}")

    for _, row in T1.iterrows():
        gene = row["Genes"]
        try:
            fig = plot_gene_single(
                adata, gene, celltype, tmeta,
                celltype_col=celltype_col, zt_col=zt_col,
                period=period, norm_str=norm_str,
                show_cells=show_cells, use_violin=use_violin,
            )
            safe_gene = re.sub(r"[^a-zA-Z0-9_]", "_", gene)
            fig.savefig(os.path.join(subdir, f"{safe_gene}.png"),
                        dpi=120, bbox_inches="tight")
            plt.close(fig)
        except Exception as e:
            warnings.warn(f"  Gene {gene} plot failed: {e}")

    print(f"  Batch plots complete.")
