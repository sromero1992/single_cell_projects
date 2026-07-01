"""
grn.py — TimeSCape gene regulatory network (GRN) module
=========================================================
Python equivalent of the GRN section in pathway_circadian.R.

Pipeline:
  1. select_hub_genes()    — Pearson correlation on pooled cells, top % by degree
  2. plot_grn_timeseries() — per-ZT correlation network panels (matplotlib)

Dependencies
------------
  pip install networkx matplotlib

networkx mirrors igraph; matplotlib mirrors R's ggraph/ggplot2.
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats as scipy_stats


# ── 1. Hub gene selection ─────────────────────────────────────────────────────

def select_hub_genes(
    X_ct: np.ndarray,
    gene_names: list[str],
    cor_thresh: float = 0.30,
    p_thresh: float = 0.05,
    hub_pct: float = 0.10,
    min_hub: int = 5,
) -> dict:
    """
    Select hub genes from a circadian gene pool by co-expression degree.

    Mirrors R's ``select_hub_genes()`` in pathway_circadian.R.

    Approach:
      1. Z-score each gene across cells (standardise).
      2. Compute the gene × gene Pearson correlation matrix.
      3. Convert to t-statistics to obtain p-values (df = n_cells − 2).
      4. Build adjacency: edge exists if |r| >= cor_thresh AND p < p_thresh.
      5. Rank genes by degree (number of significant edges); return top hub_pct.

    Why pool all ZT time points?
    The gene pool consists of confirmed circadian oscillators. Their pooled
    correlation is dominated by phase relationships — co-phased genes are
    positively correlated; anti-phased genes (e.g. clock activators vs.
    repressors) are negatively correlated. Pooling gives statistical power
    that per-ZT slices (small n) cannot. Time-varying biology is recovered
    downstream by ``plot_grn_timeseries()``, which recomputes per-ZT edge
    weights independently.

    Both positive and negative correlations count toward degree because
    anti-phase coupling (e.g. Arntl vs Nr1d1) is the core feedback loop.
    Signed r is preserved as edge weight in the output graph for directional
    interpretation.

    Parameters
    ----------
    X_ct : np.ndarray  (n_cells × n_genes) OR (n_genes × n_cells)
        Dense expression matrix for ONE cell type, ALL ZT time points pooled.
        If shape[0] > shape[1] it is assumed cells × genes and transposed.
        Pass only the genes in ``gene_names`` (already filtered).
    gene_names : list of str
        Gene labels for columns of X_ct (or rows if genes × cells).
    cor_thresh : float
        Minimum |Pearson r| for an edge to count toward degree (default 0.30).
    p_thresh : float
        Maximum p-value for an edge (default 0.05).
    hub_pct : float
        Top fraction of genes by degree to call hubs (default 0.10 = top 10%).
    min_hub : int
        Minimum number of hubs to return; cor_thresh is relaxed step-wise
        (−0.05 per step) if too few genes pass (default 5).

    Returns
    -------
    dict with keys:
        ``genes``      — list of gene names actually used
        ``cor_mat``    — np.ndarray (n_genes × n_genes) Pearson correlation matrix
        ``adj_mat``    — np.ndarray bool, True where edge exists
        ``degree``     — np.ndarray int, edge count per gene
        ``hub_genes``  — list of hub gene names (sorted by degree descending)
        ``n_cells``    — int, number of cells used
    """
    X = np.asarray(X_ct, dtype=np.float64)

    # Ensure shape is genes × cells
    if X.shape[0] > X.shape[1]:
        X = X.T   # now genes × cells

    n_genes, n_cells = X.shape

    if n_cells < 10:
        raise ValueError(
            f"select_hub_genes: only {n_cells} cells — too few for reliable correlation. "
            "Use at least 10 cells (ideally >50)."
        )
    if n_genes < 5:
        raise ValueError(
            f"select_hub_genes: only {n_genes} genes in pool — too few to build a network."
        )

    # Z-score each gene (row) across cells
    mu  = X.mean(axis=1, keepdims=True)
    sig = X.std(axis=1, keepdims=True)
    sig[sig < 1e-12] = 1.0          # constant genes → zero correlation after z-scoring
    X_z = (X - mu) / sig

    # Gene × gene Pearson correlation
    # r[i,j] = (X_z[i] · X_z[j]) / (n_cells - 1)
    cor_mat = (X_z @ X_z.T) / (n_cells - 1)
    np.clip(cor_mat, -1.0, 1.0, out=cor_mat)
    np.fill_diagonal(cor_mat, 0.0)

    # t-statistic → p-value (df = n_cells - 2)
    df  = n_cells - 2
    eps = 1e-10
    denom = np.sqrt(np.maximum(1.0 - cor_mat**2, eps))
    t_mat = cor_mat * np.sqrt(df) / denom
    # Two-sided p-value
    p_mat = 2 * scipy_stats.t.sf(np.abs(t_mat), df=df)
    np.fill_diagonal(p_mat, 1.0)

    # Adjacency: |r| >= cor_thresh AND p < p_thresh
    # Use only upper triangle to avoid double-counting, then symmetrise
    def _build_adj(cth):
        adj = (np.abs(cor_mat) >= cth) & (p_mat < p_thresh)
        np.fill_diagonal(adj, False)
        return adj

    adj_mat = _build_adj(cor_thresh)
    degree  = adj_mat.sum(axis=1).astype(int)

    # Relax threshold if too few hubs would be selected
    current_thresh = cor_thresh
    while True:
        n_top_needed = max(min_hub, int(np.ceil(hub_pct * n_genes)))
        n_with_edges = int((degree > 0).sum())
        if n_with_edges >= min_hub:
            break
        current_thresh -= 0.05
        if current_thresh < 0.05:
            warnings.warn(
                "select_hub_genes: very few edges even at low threshold. "
                "Check that the gene pool has sufficient expression variation.",
                RuntimeWarning,
                stacklevel=2,
            )
            break
        adj_mat = _build_adj(current_thresh)
        degree  = adj_mat.sum(axis=1).astype(int)

    if current_thresh != cor_thresh:
        print(f"  select_hub_genes: cor_thresh relaxed to {current_thresh:.2f} "
              f"to find ≥{min_hub} genes with edges")

    # Rank by degree and select top hub_pct
    n_hubs     = max(min_hub, int(np.ceil(hub_pct * n_genes)))
    sorted_idx = np.argsort(-degree)          # descending
    hub_idx    = sorted_idx[:n_hubs]
    hub_idx    = hub_idx[degree[hub_idx] > 0]  # must have at least one edge

    hub_genes = [gene_names[i] for i in hub_idx]

    print(f"  select_hub_genes: {n_genes} genes, {n_cells} cells, "
          f"threshold |r|≥{current_thresh:.2f} → {len(hub_genes)} hubs")

    return {
        "genes":     gene_names,
        "cor_mat":   cor_mat,
        "adj_mat":   adj_mat,
        "degree":    degree,
        "hub_genes": hub_genes,
        "n_cells":   n_cells,
    }


# ── 2. Per-ZT network builder ─────────────────────────────────────────────────

def _build_zt_network(
    X_gc: np.ndarray,
    gene_names: list[str],
    gene_subset: list[str],
    zt_v: np.ndarray,
    zt_target: float,
    cor_thresh: float = 0.20,
    p_thresh: float = 0.05,
):
    """
    Build a correlation network for ``gene_subset`` at one ZT time point.

    Returns
    -------
    (nodes, edges) where:
        nodes : list of str — gene names present in data
        edges : list of (gene_i, gene_j, r) tuples
    """
    cell_idx = np.where(zt_v == zt_target)[0]
    if len(cell_idx) < 3:
        return [], []

    gset = set(gene_subset)
    g_idx = [i for i, g in enumerate(gene_names) if g in gset]
    g_present = [gene_names[i] for i in g_idx]
    if len(g_present) < 2:
        return g_present, []

    sub = X_gc[np.ix_(g_idx, cell_idx)].astype(np.float64)   # genes × cells

    ng, nc = sub.shape
    if nc < 3:
        return g_present, []

    df  = nc - 2
    mu  = sub.mean(axis=1, keepdims=True)
    sig = sub.std(axis=1, keepdims=True)
    sig[sig < 1e-12] = 1.0
    z   = (sub - mu) / sig
    r   = (z @ z.T) / (nc - 1)
    np.clip(r, -1.0, 1.0, out=r)

    eps = 1e-10
    t   = r * np.sqrt(df) / np.sqrt(np.maximum(1 - r**2, eps))
    p   = 2 * scipy_stats.t.sf(np.abs(t), df=df)

    edges = []
    for i in range(ng):
        for j in range(i + 1, ng):
            if abs(r[i, j]) >= cor_thresh and p[i, j] < p_thresh:
                edges.append((g_present[i], g_present[j], float(r[i, j])))

    return g_present, edges


# ── 3. Per-ZT GRN time series ─────────────────────────────────────────────────

def plot_grn_timeseries(
    X_gc: np.ndarray,
    gene_names: list[str],
    zt_v: np.ndarray,
    hub_genes: list[str],
    actual_times: np.ndarray,
    title: str = "",
    cor_thresh: float = 0.20,
    p_thresh: float = 0.05,
    figsize_per_panel: tuple[float, float] = (3.5, 3.5),
    n_cols: int = 4,
    outpath: str | None = None,
):
    """
    Build one co-expression network panel per ZT time point and arrange in a grid.

    A fixed Fruchterman-Reingold layout is computed from the pooled (all-ZT)
    correlation network so node positions are stable across panels.  Per-ZT
    edge weights are recomputed independently using only cells from each slice.

    Parameters
    ----------
    X_gc : np.ndarray  (genes × cells)
        Dense expression matrix for the target cell type, ALL ZT pooled.
    gene_names : list of str
        Gene labels for rows of ``X_gc``.
    zt_v : np.ndarray
        Per-cell numeric ZT hour (aligned to columns of ``X_gc``).
    hub_genes : list of str
        Genes to display in the network (output of ``select_hub_genes()``).
        Should be 3–15 genes for readable plots.
    actual_times : np.ndarray
        Unique ZT hours to display (one panel each).
    title : str
        Overall figure title.
    cor_thresh : float
        Minimum |r| for an edge to be drawn (default 0.20).
    p_thresh : float
        Maximum p-value for an edge (default 0.05).
    figsize_per_panel : tuple
        (width, height) in inches per panel (default (3.5, 3.5)).
    n_cols : int
        Number of columns in the panel grid (default 4).
    outpath : str, optional
        If given, save the figure to this path (PNG/PDF).

    Returns
    -------
    matplotlib.figure.Figure
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
    except ImportError:
        raise ImportError(
            "networkx and matplotlib are required for plot_grn_timeseries(). "
            "Install with:  pip install networkx matplotlib"
        )

    n_panels = len(actual_times)
    n_rows   = int(np.ceil(n_panels / n_cols))
    fig_w    = figsize_per_panel[0] * n_cols
    fig_h    = figsize_per_panel[1] * n_rows

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(fig_w, fig_h),
                             squeeze=False)
    if title:
        fig.suptitle(title, fontsize=12, fontweight="bold", y=1.01)

    # Compute pooled (fixed) layout using all-ZT network
    _, pooled_edges = _build_zt_network(
        X_gc, gene_names, hub_genes,
        zt_v, zt_v[0],          # use first ZT just to get gene subset
        cor_thresh=0.05,         # low threshold for layout stability
        p_thresh=1.0,
    )
    # Use all hub genes for layout even if some have no edges
    gset = set(hub_genes) & set(gene_names)
    G_layout = nx.Graph()
    G_layout.add_nodes_from(sorted(gset))
    _, all_edges = _build_zt_network(
        X_gc, gene_names, hub_genes,
        zt_v, zt_v[0],
        cor_thresh=0.0, p_thresh=1.0,
    )
    pos = nx.spring_layout(G_layout, seed=42, k=1.5 / max(len(gset)**0.5, 1))

    cmap = cm.RdBu_r
    norm = mcolors.Normalize(vmin=-1, vmax=1)

    for panel_i, zt in enumerate(actual_times):
        row = panel_i // n_cols
        col = panel_i % n_cols
        ax  = axes[row][col]

        nodes, edges = _build_zt_network(
            X_gc, gene_names, hub_genes,
            zt_v, zt,
            cor_thresh=cor_thresh,
            p_thresh=p_thresh,
        )

        G = nx.Graph()
        G.add_nodes_from(sorted(gset))
        for (u, v, r) in edges:
            G.add_edge(u, v, weight=r, color=cmap(norm(r)))

        node_pos = {g: pos[g] for g in G.nodes() if g in pos}

        edge_colors = [G[u][v]["color"] for u, v in G.edges()]
        edge_widths = [abs(G[u][v]["weight"]) * 3 for u, v in G.edges()]

        nx.draw_networkx_nodes(G, node_pos, ax=ax,
                                node_size=300, node_color="#4E79A7",
                                alpha=0.9)
        nx.draw_networkx_labels(G, node_pos, ax=ax,
                                 font_size=6, font_color="white",
                                 font_weight="bold")
        if edges:
            nx.draw_networkx_edges(G, node_pos, ax=ax,
                                    edge_color=edge_colors,
                                    width=edge_widths, alpha=0.75)

        ax.set_title(f"ZT{int(zt):02d}  ({len(edges)} edges)", fontsize=9)
        ax.axis("off")

    # Hide unused panels
    for panel_i in range(n_panels, n_rows * n_cols):
        axes[panel_i // n_cols][panel_i % n_cols].axis("off")

    # Colorbar for edge r values
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.04, label="Pearson r")

    fig.tight_layout()

    if outpath:
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        print(f"  GRN time series saved → {outpath}")

    return fig
