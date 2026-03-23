"""
Visualization functions for scReservoir.

All functions accept numpy arrays directly — no class instances needed.
Each function creates its own figure/axes if none are provided, or draws
onto existing axes if `ax` is passed.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401  (registers 3D projection)


# ---------------------------------------------------------------------------
# GRN heatmap
# ---------------------------------------------------------------------------

def plot_grn_heatmap(GRN, gene_names=None, top_n=50, ax=None, cmap='YlOrRd', figsize=(12, 10)):
    """
    Plot the gene regulatory network as a heatmap.

    Only the top `top_n` genes by out-degree (row sum) are shown so the
    plot remains readable.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        GRN[i, j] = influence of gene i on gene j.
    gene_names : np.ndarray of shape (n_genes,), optional
        Gene labels.  Defaults to integer indices.
    top_n : int, default=50
        Number of highest-activity genes to include.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  Creates a new figure if None.
    cmap : str, default='YlOrRd'
        Matplotlib colormap name.
    figsize : tuple, default=(12, 10)
        Figure size (used only when ax is None).

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    n_genes = GRN.shape[0]
    if gene_names is None:
        gene_names = np.arange(n_genes).astype(str)

    # Select top genes by out-degree (row sums)
    activity = GRN.sum(axis=1)
    top_idx  = np.argsort(activity)[-top_n:][::-1]

    GRN_top       = GRN[np.ix_(top_idx, top_idx)]
    names_top     = gene_names[top_idx]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        GRN_top,
        xticklabels=names_top,
        yticklabels=names_top,
        cmap=cmap,
        cbar_kws={'label': 'Influence strength'},
        ax=ax,
        square=True
    )
    ax.set_xlabel('Target gene (j)', fontsize=12)
    ax.set_ylabel('Regulator gene (i)', fontsize=12)
    ax.set_title(f'Gene Regulatory Network (top {top_n} genes)', fontsize=14, fontweight='bold')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0,  fontsize=8)

    return ax


# ---------------------------------------------------------------------------
# Energy landscape  (2D and 3D)
# ---------------------------------------------------------------------------

def plot_energy_landscape(energy, embedding, pseudotime=None, ax=None,
                          figsize=(12, 10), cmap='RdYlBu_r'):
    """
    Plot the Waddington energy landscape over a 2D or 3D embedding.

    Parameters
    ----------
    energy : np.ndarray of shape (n_cells,)
        Energy per cell (low = attractor / valley, high = transitional).
    embedding : np.ndarray of shape (n_cells, 2) or (n_cells, 3)
        2D or 3D cell coordinates (e.g., UMAP or PCA output).
    pseudotime : np.ndarray of shape (n_cells,), optional
        If provided, draws a dashed pseudotime trajectory on the plot.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  For 3D plots, pass a 3D axes.
    figsize : tuple, default=(12, 10)
    cmap : str, default='RdYlBu_r'

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if embedding.shape[1] == 2:
        return _plot_energy_2d(energy, embedding, pseudotime, ax, figsize, cmap)
    elif embedding.shape[1] == 3:
        return _plot_energy_3d(energy, embedding, pseudotime, ax, figsize, cmap)
    else:
        raise ValueError("embedding must have 2 or 3 columns.")


def _plot_energy_2d(energy, embedding, pseudotime, ax, figsize, cmap):
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sc = ax.scatter(
        embedding[:, 0], embedding[:, 1],
        c=energy, s=50, alpha=0.6, cmap=cmap,
        edgecolors='black', linewidth=0.5
    )
    if pseudotime is not None:
        idx = np.argsort(pseudotime)
        ax.plot(embedding[idx, 0], embedding[idx, 1],
                'k--', alpha=0.3, linewidth=1, label='Pseudotime trajectory')
        ax.legend(fontsize=10)

    ax.set_xlabel('UMAP1 / PC1', fontsize=12)
    ax.set_ylabel('UMAP2 / PC2', fontsize=12)
    ax.set_title('Waddington Energy Landscape', fontsize=14, fontweight='bold')
    plt.colorbar(sc, ax=ax).set_label('Energy', fontsize=11)
    return ax


def _plot_energy_3d(energy, embedding, pseudotime, ax, figsize, cmap):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(
        embedding[:, 0], embedding[:, 1], embedding[:, 2],
        c=energy, s=50, alpha=0.6, cmap=cmap,
        edgecolors='black', linewidth=0.5
    )
    if pseudotime is not None:
        idx = np.argsort(pseudotime)
        ax.plot(embedding[idx, 0], embedding[idx, 1], embedding[idx, 2],
                'k--', alpha=0.3, linewidth=1)

    ax.set_xlabel('PC1', fontsize=11)
    ax.set_ylabel('PC2', fontsize=11)
    ax.set_zlabel('PC3', fontsize=11)
    ax.set_title('Waddington Energy Landscape (3D)', fontsize=14, fontweight='bold')
    plt.colorbar(sc, ax=ax, pad=0.1).set_label('Energy', fontsize=11)
    return ax


# ---------------------------------------------------------------------------
# Fate probabilities
# ---------------------------------------------------------------------------

def plot_fate_probabilities(fate_probs, embedding, attractor_names=None, figsize=(16, 4)):
    """
    Plot one subplot per attractor showing each cell's probability of reaching it.

    Parameters
    ----------
    fate_probs : np.ndarray of shape (n_cells, n_attractors)
        Output of compute_fate_probabilities().
    embedding : np.ndarray of shape (n_cells, 2)
        2D cell coordinates (UMAP or PCA).
    attractor_names : list of str, optional
        Labels for each attractor.  Defaults to 'Attractor 1', 'Attractor 2', ...
    figsize : tuple, default=(16, 4)

    Returns
    -------
    axes : np.ndarray of matplotlib.axes.Axes
    """
    n_attractors = fate_probs.shape[1]
    if attractor_names is None:
        attractor_names = [f'Attractor {i + 1}' for i in range(n_attractors)]

    fig, axes = plt.subplots(1, n_attractors, figsize=figsize)
    if n_attractors == 1:
        axes = np.array([axes])

    for k, ax in enumerate(axes):
        sc = ax.scatter(
            embedding[:, 0], embedding[:, 1],
            c=fate_probs[:, k], s=50, alpha=0.7, cmap='viridis',
            edgecolors='black', linewidth=0.5, vmin=0, vmax=1
        )
        ax.set_xlabel('UMAP1 / PC1', fontsize=11)
        ax.set_ylabel('UMAP2 / PC2', fontsize=11)
        ax.set_title(attractor_names[k], fontsize=12, fontweight='bold')
        plt.colorbar(sc, ax=ax).set_label('Fate probability', fontsize=10)

    plt.tight_layout()
    return axes


# ---------------------------------------------------------------------------
# Attractor gene programs
# ---------------------------------------------------------------------------

def plot_attractor_genes(attractor_gene_programs, gene_names=None, top_n=20,
                         attractor_names=None, figsize=(14, 8), cmap='RdBu_r'):
    """
    Heatmap of the top variable genes across attractor states.

    Parameters
    ----------
    attractor_gene_programs : np.ndarray of shape (n_attractors, n_genes)
        Mean gene expression per attractor.
        Compute this as: X_sorted[cells_per_attractor[k]].mean(axis=0) for each k.
    gene_names : np.ndarray of shape (n_genes,), optional
    top_n : int, default=20
        Number of genes with highest variance across attractors to show.
    attractor_names : list of str, optional
    figsize : tuple, default=(14, 8)
    cmap : str, default='RdBu_r'

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    n_attractors, n_genes = attractor_gene_programs.shape
    if gene_names is None:
        gene_names = np.arange(n_genes).astype(str)
    if attractor_names is None:
        attractor_names = [f'Attractor {i + 1}' for i in range(n_attractors)]

    gene_var = attractor_gene_programs.std(axis=0)
    top_idx  = np.argsort(gene_var)[-top_n:][::-1]

    _, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        attractor_gene_programs[:, top_idx].T,
        xticklabels=attractor_names,
        yticklabels=gene_names[top_idx],
        cmap=cmap,
        cbar_kws={'label': 'Mean expression'},
        ax=ax,
        center=0
    )
    ax.set_xlabel('Attractor state', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)
    ax.set_title(f'Attractor Gene Programs (top {top_n} variable genes)',
                 fontsize=14, fontweight='bold')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)
    plt.tight_layout()
    return ax


# ---------------------------------------------------------------------------
# Reservoir state dynamics
# ---------------------------------------------------------------------------

def plot_reservoir_dynamics(H, pseudotime=None, n_dims=3, figsize=(14, 10)):
    """
    Visualize reservoir state trajectory after PCA reduction to 2 or 3 dimensions.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix.
    pseudotime : np.ndarray of shape (n_cells,), optional
        Colors points along the pseudotime gradient.
    n_dims : int, default=3
        2 for a 2D scatter, 3 for a 3D scatter.
    figsize : tuple, default=(14, 10)

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    from sklearn.decomposition import PCA

    pca   = PCA(n_components=n_dims)
    H_pca = pca.fit_transform(H)
    evr   = pca.explained_variance_ratio_

    if n_dims == 2:
        _, ax = plt.subplots(figsize=figsize)
        kw = dict(s=50, alpha=0.6, edgecolors='black', linewidth=0.5)

        if pseudotime is not None:
            sc = ax.scatter(H_pca[:, 0], H_pca[:, 1], c=pseudotime, cmap='viridis', **kw)
            plt.colorbar(sc, ax=ax).set_label('Pseudotime', fontsize=11)
            idx = np.argsort(pseudotime)
            ax.plot(H_pca[idx, 0], H_pca[idx, 1], 'k--', alpha=0.3, linewidth=1)
        else:
            ax.scatter(H_pca[:, 0], H_pca[:, 1], **kw)

        ax.set_xlabel(f'PC1 ({evr[0]:.1%} var)', fontsize=11)
        ax.set_ylabel(f'PC2 ({evr[1]:.1%} var)', fontsize=11)

    elif n_dims == 3:
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111, projection='3d')
        kw  = dict(s=50, alpha=0.6, edgecolors='black', linewidth=0.5)

        if pseudotime is not None:
            sc = ax.scatter(H_pca[:, 0], H_pca[:, 1], H_pca[:, 2],
                            c=pseudotime, cmap='viridis', **kw)
            plt.colorbar(sc, ax=ax, pad=0.1).set_label('Pseudotime', fontsize=11)
            idx = np.argsort(pseudotime)
            ax.plot(H_pca[idx, 0], H_pca[idx, 1], H_pca[idx, 2],
                    'k--', alpha=0.3, linewidth=1)
        else:
            ax.scatter(H_pca[:, 0], H_pca[:, 1], H_pca[:, 2], **kw)

        ax.set_xlabel(f'PC1 ({evr[0]:.1%} var)', fontsize=11)
        ax.set_ylabel(f'PC2 ({evr[1]:.1%} var)', fontsize=11)
        ax.set_zlabel(f'PC3 ({evr[2]:.1%} var)', fontsize=11)
    else:
        raise ValueError("n_dims must be 2 or 3.")

    ax.set_title('Reservoir State Dynamics', fontsize=14, fontweight='bold')
    plt.tight_layout()
    return ax
