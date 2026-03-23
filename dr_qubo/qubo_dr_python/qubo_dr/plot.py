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


def plot_gene_network(sub_Q_net, gene_names, edge_pct=95, ax=None,
                      title='Gene Co-expression Network'):
    """
    Plot force-directed network with condition-specific edge coloring.

    Edges are colored based on sign: blue for WT co-expression gain (positive),
    red for KO co-expression gain (negative). Edge width is proportional to
    absolute co-expression strength.

    Parameters
    ----------
    sub_Q_net : np.ndarray, shape (n_genes, n_genes)
        Subnetwork co-expression matrix.
    gene_names : list of str
        Gene names, length n_genes.
    edge_pct : float, optional
        Percentile threshold for edge inclusion (default 95).
        Edges with |weight| < percentile(95) are not drawn.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, creates new figure.
    title : str, optional
        Figure title (default 'Gene Co-expression Network').

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axis object.
    G : networkx.Graph
        NetworkX graph object used for layout.

    Notes
    -----
    Uses spring_layout for force-directed positioning.
    Positive edges (blue): WT condition shows higher co-expression.
    Negative edges (red): KO condition shows higher co-expression.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.figure

    n_genes = sub_Q_net.shape[0]

    # Create graph
    G = nx.Graph()
    for i in range(n_genes):
        G.add_node(i, label=gene_names[i])

    # Add edges, filtering by percentile
    edge_weights = []
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            weight = sub_Q_net[i, j]
            if weight != 0:
                edge_weights.append(np.abs(weight))

    if edge_weights:
        threshold = np.percentile(edge_weights, edge_pct)
    else:
        threshold = 0

    edges_pos = []
    edges_neg = []
    edge_widths_pos = []
    edge_widths_neg = []

    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            weight = sub_Q_net[i, j]
            if np.abs(weight) >= threshold:
                G.add_edge(i, j, weight=weight)
                # Normalize edge width
                width = 0.5 + 2.0 * (np.abs(weight) / np.max(np.abs(edge_weights)))
                if weight > 0:
                    edges_pos.append((i, j))
                    edge_widths_pos.append(width)
                else:
                    edges_neg.append((i, j))
                    edge_widths_neg.append(width)

    # Compute layout
    pos = nx.spring_layout(G, k=1, iterations=50, seed=42)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=300, ax=ax)

    # Draw positive edges (blue = WT gain)
    if edges_pos:
        nx.draw_networkx_edges(
            G, pos, edgelist=edges_pos, width=edge_widths_pos,
            edge_color='blue', alpha=0.6, ax=ax
        )

    # Draw negative edges (red = KO gain)
    if edges_neg:
        nx.draw_networkx_edges(
            G, pos, edgelist=edges_neg, width=edge_widths_neg,
            edge_color='red', alpha=0.6, ax=ax
        )

    # Draw labels
    labels = {i: gene_names[i] for i in range(n_genes)}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold', ax=ax)

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.axis('off')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', lw=2, label='WT co-expr gain'),
        Line2D([0], [0], color='red', lw=2, label='KO co-expr gain'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

    return fig, ax, G


def plot_condition_heatmaps(Xko2norm, Xwt2norm, Xdiff, gene_names):
    """
    Plot 3-panel comparison: KO, WT, and Differential co-expression heatmaps.

    Parameters
    ----------
    Xko2norm : np.ndarray, shape (genes, genes)
        KO condition Gram matrix (cosine similarity).
    Xwt2norm : np.ndarray, shape (genes, genes)
        WT condition Gram matrix (cosine similarity).
    Xdiff : np.ndarray, shape (genes, genes)
        Differential co-expression matrix (WT - KO).
    gene_names : list of str
        Gene names for labels.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object with 3 subplots.

    Notes
    -----
    Left: KO co-expression (blue colormap)
    Middle: WT co-expression (blue colormap)
    Right: Differential (red-blue divergent, centered at 0)
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # KO heatmap
    im0 = axes[0].imshow(Xko2norm, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    axes[0].set_xticks(np.arange(len(gene_names)))
    axes[0].set_yticks(np.arange(len(gene_names)))
    axes[0].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[0].set_yticklabels(gene_names, fontsize=8)
    axes[0].set_title('KO Co-expression', fontweight='bold')
    plt.colorbar(im0, ax=axes[0], label='Similarity')

    # WT heatmap
    im1 = axes[1].imshow(Xwt2norm, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    axes[1].set_xticks(np.arange(len(gene_names)))
    axes[1].set_yticks(np.arange(len(gene_names)))
    axes[1].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[1].set_yticklabels(gene_names, fontsize=8)
    axes[1].set_title('WT Co-expression', fontweight='bold')
    plt.colorbar(im1, ax=axes[1], label='Similarity')

    # Differential heatmap
    max_abs = np.max(np.abs(Xdiff))
    im2 = axes[2].imshow(Xdiff, cmap='RdBu_r', aspect='auto',
                         vmin=-max_abs, vmax=max_abs)
    axes[2].set_xticks(np.arange(len(gene_names)))
    axes[2].set_yticks(np.arange(len(gene_names)))
    axes[2].set_xticklabels(gene_names, rotation=90, fontsize=8)
    axes[2].set_yticklabels(gene_names, fontsize=8)
    axes[2].set_title('Differential (WT - KO)', fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=axes[2], label='Difference')

    plt.tight_layout()

    return fig
