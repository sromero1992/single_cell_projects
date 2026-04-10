# Author: Selim Romero, Texas A&M University
"""
Visualization module for co-expression networks and heatmaps.

Generates publication-quality figures for differential co-expression matrices,
network graphs with condition-specific edge coloring, and comparative heatmaps.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import TwoSlopeNorm


def plot_coexpr_heatmap(mat, gene_names, title='Differential Co-expression',
                        cmap='RdBu_r', ax=None, vmin=None, vmax=None):
    """
    Plot co-expression heatmap with divergent colormap.

    Parameters
    ----------
    mat : np.ndarray, shape (n_genes, n_genes)
        Co-expression (or differential co-expression) matrix.
    gene_names : list of str
        Gene names for axis labels, length n_genes.
    title : str, optional
        Figure title (default 'Differential Co-expression').
    cmap : str, optional
        Colormap name (default 'RdBu_r' for divergent red-blue-red).
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, creates new figure.
    vmin, vmax : float, optional
        Min/max values for colormap. If None, uses data range.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axis object.

    Notes
    -----
    Uses divergent colormap centered at 0 for better visualization of
    positive and negative co-expression changes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.figure

    # Set color normalization
    if vmin is None:
        vmin = np.min(mat)
    if vmax is None:
        vmax = np.max(mat)

    # Center colormap at 0 for divergent colormaps
    if 'RdBu' in cmap or 'diverging' in cmap.lower():
        max_abs = np.max(np.abs([vmin, vmax]))
        vmin, vmax = -max_abs, max_abs

    # Plot heatmap
    im = ax.imshow(mat, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)

    # Set ticks and labels
    ax.set_xticks(np.arange(len(gene_names)))
    ax.set_yticks(np.arange(len(gene_names)))
    ax.set_xticklabels(gene_names, rotation=90, fontsize=8)
    ax.set_yticklabels(gene_names, fontsize=8)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Co-expression', rotation=270, labelpad=15)

    ax.set_title(title, fontsize=12, fontweight='bold')

    return fig, ax


def plot_gene_network(dS, gene_names, selected_mask=None,
                      top_pct=10,
                      ax=None,
                      title='Gene Co-expression Network',
                      test_label='test', reference_label='reference',
                      perm_pval=None, perm_zscore=None):
    """
    Plot force-directed differential co-expression network.

    Displays all G pathway genes as nodes.  QUBO-selected genes are drawn
    larger and in orange; background pathway genes are small and light blue.

    Edge colour encodes the direction of co-expression change:

    * **Red**  — condition A (test/disease) co-expression gain:
      dS[i,j] < 0  →  S_A[i,j] > S_B[i,j]
    * **Blue** — condition B (reference/control) co-expression gain:
      dS[i,j] > 0  →  S_B[i,j] > S_A[i,j]

    Edge width is proportional to |dS[i,j]|.  Only the top ``top_pct``%
    of gene pairs by |dS| are drawn so the graph stays readable — this
    is purely a **display threshold**, not a statistical filter.

    Hub-level significance from ``test_differential_hub()`` can be passed
    as ``perm_pval`` / ``perm_zscore`` to annotate the plot.

    Parameters
    ----------
    dS : np.ndarray, shape (G, G)
        Differential co-expression matrix passed to the plot.
        Convention: positive entry = condition B (reference) gain.
        Pass Q0 (= raw dS) for the biological signal, or Q1
        (= dS − MNN − penalty terms) for the optimisation landscape.
    gene_names : list of str, length G
        Gene names for all pathway genes.
    selected_mask : np.ndarray, dtype bool, shape (G,), optional
        Boolean mask of QUBO-selected genes.
    top_pct : float, optional
        Percentage of gene pairs (by |dS|) to draw as edges (default 10).
        Set to 100 to show every edge (hairball for large pathways).
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.  Creates a new figure if None.
    title : str, optional
        Figure title.
    test_label : str, optional
        Label for condition A (test / disease) — used in legend.
    reference_label : str, optional
        Label for condition B (reference / control) — used in legend.
    perm_pval : float, optional
        Hub permutation p-value from ``test_differential_hub()``.
        If provided, annotated below the plot.
    perm_zscore : float, optional
        Hub z-score from ``test_differential_hub()``.

    Returns
    -------
    fig, ax, G_nx : (Figure, Axes, networkx.Graph)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    else:
        fig = ax.figure

    # Negate so positive = condition A (test/disease) co-expression gain
    # dS = S_B − S_A  →  mat = S_A − S_B  (positive = disease gain)
    mat     = -dS
    n_genes = mat.shape[0]

    # ── Edge display threshold: top top_pct% by |mat| ─────────────────────
    triu_idx  = np.triu_indices(n_genes, k=1)
    triu_vals = np.abs(mat[triu_idx])
    max_w     = triu_vals.max() if triu_vals.size > 0 else 1.0

    if triu_vals.size > 0 and top_pct < 100:
        threshold = np.percentile(triu_vals, 100 - top_pct)
    else:
        threshold = 0.0                    # show all edges (top_pct=100)

    # ── Build graph ────────────────────────────────────────────────────────
    G_nx = nx.Graph()
    for i in range(n_genes):
        G_nx.add_node(i, label=gene_names[i])

    edges_pos, edges_neg             = [], []
    edge_widths_pos, edge_widths_neg = [], []

    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            w = mat[i, j]
            if np.abs(w) >= threshold:
                G_nx.add_edge(i, j, weight=w)
                width = 0.3 + 2.7 * (np.abs(w) / max_w)
                if w > 0:
                    edges_pos.append((i, j))
                    edge_widths_pos.append(width)
                else:
                    edges_neg.append((i, j))
                    edge_widths_neg.append(width)

    # ── Layout: weight by edge strength so hub genes cluster ──────────────
    pos = nx.spring_layout(G_nx, k=1.2, iterations=80, seed=42,
                           weight='weight')

    # ── Draw nodes ─────────────────────────────────────────────────────────
    if selected_mask is not None:
        sel_nodes = [i for i in range(n_genes) if selected_mask[i]]
        bg_nodes  = [i for i in range(n_genes) if not selected_mask[i]]
        nx.draw_networkx_nodes(G_nx, pos, nodelist=bg_nodes,
                               node_color='lightblue', node_size=150,
                               alpha=0.6, ax=ax)
        nx.draw_networkx_nodes(G_nx, pos, nodelist=sel_nodes,
                               node_color='orange', node_size=450,
                               alpha=1.0, ax=ax)
    else:
        nx.draw_networkx_nodes(G_nx, pos, node_color='lightblue',
                               node_size=200, ax=ax)

    # ── Draw edges ─────────────────────────────────────────────────────────
    if edges_pos:
        nx.draw_networkx_edges(G_nx, pos, edgelist=edges_pos,
                               width=edge_widths_pos,
                               edge_color='red', alpha=0.65, ax=ax)
    if edges_neg:
        nx.draw_networkx_edges(G_nx, pos, edgelist=edges_neg,
                               width=edge_widths_neg,
                               edge_color='royalblue', alpha=0.65, ax=ax)

    # ── Labels ─────────────────────────────────────────────────────────────
    if selected_mask is not None:
        sel_labels = {i: gene_names[i] for i in range(n_genes) if selected_mask[i]}
        bg_labels  = {i: gene_names[i] for i in range(n_genes) if not selected_mask[i]}
        nx.draw_networkx_labels(G_nx, pos, bg_labels,
                                font_size=6, font_color='grey', ax=ax)
        nx.draw_networkx_labels(G_nx, pos, sel_labels,
                                font_size=8, font_weight='bold', ax=ax)
    else:
        nx.draw_networkx_labels(G_nx, pos,
                                {i: gene_names[i] for i in range(n_genes)},
                                font_size=7, ax=ax)

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.axis('off')

    # ── Hub significance annotation (permutation test result) ─────────────
    if perm_pval is not None:
        if   perm_pval < 0.001: stars = '***'
        elif perm_pval < 0.01:  stars = '**'
        elif perm_pval < 0.05:  stars = '*'
        else:                   stars = 'n.s.'
        ann = f'Hub significance:  p = {perm_pval:.3f}  {stars}'
        if perm_zscore is not None:
            ann += f'   z = {perm_zscore:.2f}'
        ax.text(0.5, -0.03, ann, transform=ax.transAxes,
                ha='center', va='top', fontsize=9, style='italic',
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='lightyellow', alpha=0.85))

    # ── Legend ─────────────────────────────────────────────────────────────
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red',       lw=2,
               label=f'{test_label} co-expr gain  (dS < 0)'),
        Line2D([0], [0], color='royalblue', lw=2,
               label=f'{reference_label} co-expr gain  (dS > 0)'),
    ]
    if selected_mask is not None:
        legend_elements += [
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor='orange',   markersize=10,
                   label='QUBO selected (hub)'),
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor='lightblue', markersize=8,
                   label='pathway gene'),
        ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8,
              framealpha=0.8)

    return fig, ax, G_nx


def plot_hub_network(results,
                     top_pct=None,
                     edge_threshold=None,
                     n_sigma=None,
                     layout='spring',
                     ax=None,
                     title=None,
                     test_label=None,
                     reference_label=None,
                     figsize=(10, 9)):
    """
    Plot the QUBO-selected hub subnetwork from run_pipeline() output.

    Works directly on the dict returned by ``run_pipeline()`` or
    ``run_pipeline_classical()``.  Only the K selected hub genes are drawn
    (no background pathway genes), giving a much cleaner view than the
    full-pathway network produced during the pipeline run.

    Edge display is controlled by **exactly one** of three parameters:

    * ``top_pct``        – keep the top X% of hub edges by |dS|.
    * ``edge_threshold`` – keep edges with |dS(i,j)| ≥ threshold.
    * ``n_sigma``        – keep edges with |dS(i,j)| ≥ n_sigma × σ_edge,
      where σ_edge is estimated from the permutation null distribution.
      This is the statistically principled mode: each displayed edge must
      rise above the per-edge noise level implied by random-label shuffling.
      Requires ``run_pipeline(run_perm_test=True)``.

    If none of the three is specified, defaults to ``top_pct=25``.

    Parameters
    ----------
    results : dict
        Output of ``run_pipeline()`` or ``run_pipeline_classical()``.
        Must contain ``'sub_Q0'`` (K×K dS submatrix) and ``'sub_g_net'``
        (list of K selected gene names).
    top_pct : float, optional
        Percentage of hub-internal edges to display (by |dS| strength).
    edge_threshold : float, optional
        Absolute |dS| cutoff for edge display.
    n_sigma : float, optional
        Multiples of the implied per-edge null noise floor required for an
        edge to be displayed.  Derived from the permutation null as::

            σ_edge = σ_null / (2 × √(K(K−1)/2))

        where σ_null = std(E_null) and E = 2·Σ_{i<j} dS(i,j) sums
        n_pairs = K(K−1)/2 terms.  Under the assumption that edge terms
        are approximately independent under random label shuffling (isotropic
        null noise), per-edge variance ≈ σ²_null / (4·n_pairs).
        This is a display filter — NOT a per-edge statistical test.
        σ_edge is a uniform noise floor; n_sigma sets how many multiples of
        it an edge must exceed to be shown.  The independence assumption is
        approximate (hub genes share co-expression partners); interpret
        σ_edge as a principled heuristic, not a rigorous standard error.
        Requires permutation test to have been run (run_perm_test=True).
    layout : {'spring', 'circular', 'kamada_kawai'}, optional
        NetworkX layout algorithm (default 'spring').
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  Creates a new figure if None.
    title : str, optional
        Figure title.  Auto-generated (with edge count and filter info) if None.
    test_label : str, optional
        Display name for the test / disease condition.
    reference_label : str, optional
        Display name for the reference / control condition.
    figsize : tuple, optional
        Figure size in inches (default (10, 9)).

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax  : matplotlib.axes.Axes
    G_nx : networkx.Graph
        Graph built from the displayed hub genes and edges.

    Examples
    --------
    >>> results = run_pipeline(X, g, batch, ...)
    >>> # Statistically principled: show edges above 1-sigma null noise
    >>> fig, ax, G = plot_hub_network(results, n_sigma=1.0)
    >>> # Quick exploration: top 30% of edges
    >>> fig, ax, G = plot_hub_network(results, top_pct=30)
    >>> # Absolute threshold
    >>> fig, ax, G = plot_hub_network(results, edge_threshold=0.05)
    """
    # ── Unpack results dict ───────────────────────────────────────────────
    sub_Q0    = results.get('sub_Q0')
    sub_g_net = results.get('sub_g_net')

    if sub_Q0 is None or sub_g_net is None:
        raise KeyError(
            "'results' must contain 'sub_Q0' and 'sub_g_net'. "
            "These are produced by run_pipeline() and run_pipeline_classical()."
        )

    perm_result = results.get('hub_perm_result')
    perm_pval   = perm_result['pval']    if perm_result else None
    perm_zscore = perm_result['zscore']  if perm_result else None
    E_null      = perm_result.get('E_null') if perm_result else None

    if test_label is None:
        test_label      = results.get('test_label',      'test')
    if reference_label is None:
        reference_label = results.get('reference_label', 'reference')

    K             = len(sub_g_net)
    n_pairs       = K * (K - 1) // 2          # unique off-diagonal pairs

    # mat: positive value → test/disease co-expression gain (negated dS)
    mat = -sub_Q0

    triu_idx  = np.triu_indices(K, k=1)
    triu_vals = np.abs(mat[triu_idx])
    max_w     = triu_vals.max() if triu_vals.size > 0 else 1.0

    # ── Edge display threshold ────────────────────────────────────────────
    if n_sigma is not None:
        if E_null is None:
            raise ValueError(
                "n_sigma requires permutation-test results. "
                "Re-run run_pipeline() with run_perm_test=True."
            )
        null_std = float(np.std(E_null))

        # --- σ_edge derivation -------------------------------------------
        # The hub energy is a sum over n_pairs = K(K-1)/2 edge terms:
        #   E = 2 × Σ_{i<j} dS(i,j)
        # so σ²_null = Var(E_null) aggregates variance from all n_pairs edges.
        #
        # Working assumption: under random label shuffling, individual edge
        # contributions are approximately INDEPENDENT with equal variance.
        # Justification: permutation destroys systematic co-expression
        # structure, leaving approximately isotropic finite-sample noise —
        # no particular gene pair is systematically noisier than any other.
        #
        # Under this assumption:
        #   σ²_null ≈ 4 × n_pairs × σ²_per_edge
        #   ⟹  σ_edge = σ_null / (2 × √n_pairs)
        #
        # IMPORTANT CAVEAT: the independence assumption is approximate.
        # Within a hub, gene i appears in K-1 edges simultaneously, so
        # edges (i,j) and (i,k) share gene i and are correlated under
        # permutation. σ_edge is therefore an implied per-edge noise floor
        # (the per-edge noise that *would* explain the total hub variance if
        # edges were independent), NOT a rigorous per-edge standard error.
        #
        # Practical interpretation:
        #   - σ_edge is a single fixed threshold applied uniformly to all edges.
        #   - |dS(i,j)| ≥ n_sigma × σ_edge keeps edges visibly above null noise.
        #   - This is a display filter calibrated by the permutation null,
        #     NOT a per-edge statistical test (no per-edge p-value is computed).
        # -----------------------------------------------------------------
        sigma_edge        = null_std / (2.0 * np.sqrt(max(n_pairs, 1)))
        display_threshold = n_sigma * sigma_edge
        mode_label = (f'|dS| ≥ {n_sigma}σ_edge'
                      f'  (σ_edge = {sigma_edge:.4f}, '
                      f'σ_null = {null_std:.4f})')
    elif edge_threshold is not None:
        display_threshold = float(edge_threshold)
        mode_label = f'|dS| ≥ {edge_threshold:.4f}'
    else:
        pct = float(top_pct) if top_pct is not None else 25.0
        display_threshold = (np.percentile(triu_vals, 100.0 - pct)
                             if triu_vals.size > 0 and pct < 100
                             else 0.0)
        mode_label = f'top {pct:.4g}% edges'

    # ── Build graph ───────────────────────────────────────────────────────
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    G_nx = nx.Graph()
    for i, name in enumerate(sub_g_net):
        G_nx.add_node(i, label=name)

    edges_pos, edges_neg             = [], []
    edge_widths_pos, edge_widths_neg = [], []

    for i in range(K):
        for j in range(i + 1, K):
            w = mat[i, j]
            if np.abs(w) >= display_threshold:
                G_nx.add_edge(i, j, weight=float(np.abs(w)))
                width = 0.5 + 3.5 * (np.abs(w) / max_w)
                if w > 0:
                    edges_pos.append((i, j))
                    edge_widths_pos.append(width)
                else:
                    edges_neg.append((i, j))
                    edge_widths_neg.append(width)

    n_shown = len(edges_pos) + len(edges_neg)

    # ── Layout ────────────────────────────────────────────────────────────
    if layout == 'circular':
        pos = nx.circular_layout(G_nx)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G_nx, weight='weight')
    else:
        pos = nx.spring_layout(G_nx, k=1.5, iterations=100,
                               seed=42, weight='weight')

    # ── Draw nodes & labels (all hub genes in orange) ─────────────────────
    nx.draw_networkx_nodes(G_nx, pos,
                           node_color='orange', node_size=800,
                           alpha=1.0, ax=ax)
    nx.draw_networkx_labels(G_nx, pos,
                            {i: sub_g_net[i] for i in range(K)},
                            font_size=12, font_weight='bold', ax=ax)

    # ── Draw edges ────────────────────────────────────────────────────────
    if edges_pos:
        nx.draw_networkx_edges(G_nx, pos, edgelist=edges_pos,
                               width=edge_widths_pos,
                               edge_color='red', alpha=0.70, ax=ax)
    if edges_neg:
        nx.draw_networkx_edges(G_nx, pos, edgelist=edges_neg,
                               width=edge_widths_neg,
                               edge_color='royalblue', alpha=0.70, ax=ax)

    # ── Title ─────────────────────────────────────────────────────────────
    if title is None:
        title = (f'Q0 Hub  (K={K})  —  {n_shown} / {n_pairs} edges shown'
                 f'\nFilter: {mode_label}')
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axis('off')

    # ── Significance annotation ───────────────────────────────────────────
    if perm_pval is not None:
        if   perm_pval < 0.001: stars = '***'
        elif perm_pval < 0.01:  stars = '**'
        elif perm_pval < 0.05:  stars = '*'
        else:                   stars = 'n.s.'
        ann = f'Hub significance:  p = {perm_pval:.3f}  {stars}'
        if perm_zscore is not None:
            ann += f'   z = {perm_zscore:.2f}'
        ax.text(0.5, -0.03, ann, transform=ax.transAxes,
                ha='center', va='top', fontsize=9, style='italic',
                bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='lightyellow', alpha=0.85))

    # ── Legend ────────────────────────────────────────────────────────────
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='red',       lw=2,
               label=f'{test_label} co-expr gain  (dS < 0)'),
        Line2D([0], [0], color='royalblue', lw=2,
               label=f'{reference_label} co-expr gain  (dS > 0)'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='orange', markersize=10,
               label='QUBO hub gene'),
    ]
    ax.legend(handles=legend_elements, loc='upper left',
              fontsize=8, framealpha=0.8)

    return fig, ax, G_nx


def plot_condition_heatmaps(Xa2norm, Xb2norm, dS, gene_names,
                            test_label='test', reference_label='control'):
    """
    Plot 3-panel comparison: test, reference, and Differential co-expression heatmaps.

    Parameters
    ----------
    Xa2norm : np.ndarray, shape (genes, genes)
        Condition A (test) Gram matrix (cosine similarity).
    Xb2norm : np.ndarray, shape (genes, genes)
        Condition B (reference) Gram matrix (cosine similarity).
    dS : np.ndarray, shape (genes, genes)
        Differential co-expression matrix (reference - test).
    gene_names : list of str
        Gene names for labels.
    test_label : str, optional
        Display label for condition A / test (default 'test').
    reference_label : str, optional
        Display label for condition B / reference (default 'control').

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object with 3 subplots.

    Notes
    -----
    Left:   test co-expression (Blues colormap, 0–1).
    Middle: reference co-expression (Blues colormap, 0–1).
    Right:  Differential = test − reference (RdBu_r, centered at 0).
            Red  → test co-expression gain  (positive = test > reference).
            Blue → reference co-expression gain (negative = reference > test).
            Internally computed as −dS because dS = reference − test.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Test (condition A) heatmap
    im0 = axes[0].imshow(Xa2norm, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    axes[0].set_xticks(np.arange(len(gene_names)))
    axes[0].set_yticks(np.arange(len(gene_names)))
    axes[0].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[0].set_yticklabels(gene_names, fontsize=8)
    axes[0].set_title(f'{test_label} Co-expression', fontweight='bold')
    plt.colorbar(im0, ax=axes[0], label='Similarity')

    # Reference (condition B) heatmap
    im1 = axes[1].imshow(Xb2norm, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    axes[1].set_xticks(np.arange(len(gene_names)))
    axes[1].set_yticks(np.arange(len(gene_names)))
    axes[1].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[1].set_yticklabels(gene_names, fontsize=8)
    axes[1].set_title(f'{reference_label} Co-expression', fontweight='bold')
    plt.colorbar(im1, ax=axes[1], label='Similarity')

    # Differential heatmap
    # Show  (test − reference) = −dS  so that positive values (red) mean
    # the test condition gained co-expression — biologically natural.
    diff_display = -dS
    max_abs = np.max(np.abs(diff_display))
    im2 = axes[2].imshow(diff_display, cmap='RdBu_r', aspect='auto',
                         vmin=-max_abs, vmax=max_abs)
    axes[2].set_xticks(np.arange(len(gene_names)))
    axes[2].set_yticks(np.arange(len(gene_names)))
    axes[2].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[2].set_yticklabels(gene_names, fontsize=8)
    axes[2].set_title(f'Differential ({test_label} − {reference_label})', fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=axes[2], label=f'{test_label} − {reference_label}')

    plt.tight_layout()

    return fig
