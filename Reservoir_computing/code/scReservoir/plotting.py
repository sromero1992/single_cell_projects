"""
Visualization module for scReservoir.

Provides plotting functions for GRN heatmaps, energy landscapes,
fate probabilities, and attractor programs.
"""

from typing import Optional, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns


def plot_grn_heatmap(
    GRN: np.ndarray,
    gene_names: Optional[np.ndarray] = None,
    top_n: int = 50,
    ax: Optional[plt.Axes] = None,
    cmap: str = 'YlOrRd',
    figsize: Tuple[int, int] = (12, 10)
) -> plt.Axes:
    """
    Plot gene regulatory network as heatmap.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene regulatory network matrix.
    gene_names : np.ndarray of shape (n_genes,), optional
        Gene names/IDs. If None, uses integer indices.
    top_n : int, default=50
        Number of top genes by mean activity to include.
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. Creates new figure if None.
    cmap : str, default='YlOrRd'
        Colormap name.
    figsize : tuple, default=(12, 10)
        Figure size if creating new axes.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes object with plot.
    """
    n_genes = GRN.shape[0]

    if gene_names is None:
        gene_names = np.arange(n_genes).astype(str)

    # Select top genes by activity (out-degree)
    activity = GRN.sum(axis=1)
    top_idx = np.argsort(activity)[-top_n:][::-1]

    GRN_top = GRN[np.ix_(top_idx, top_idx)]
    gene_names_top = gene_names[top_idx]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    sns.heatmap(
        GRN_top,
        xticklabels=gene_names_top,
        yticklabels=gene_names_top,
        cmap=cmap,
        cbar_kws={'label': 'Influence strength'},
        ax=ax,
        square=True
    )

    ax.set_xlabel('Target gene (j)', fontsize=12)
    ax.set_ylabel('Regulator gene (i)', fontsize=12)
    ax.set_title(f'Gene Regulatory Network (top {top_n} genes)', fontsize=14, fontweight='bold')

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    return ax


def plot_energy_landscape(
    energy: np.ndarray,
    embedding: np.ndarray,
    pseudotime: Optional[np.ndarray] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[int, int] = (12, 10),
    cmap: str = 'RdYlBu_r'
) -> plt.Axes:
    """
    Plot energy landscape over 2D or 3D embedding.

    Parameters
    ----------
    energy : np.ndarray of shape (n_cells,)
        Energy values per cell (lower = attractors).
    embedding : np.ndarray of shape (n_cells, 2) or (n_cells, 3)
        2D or 3D cell coordinates (e.g., from PCA or UMAP).
    pseudotime : np.ndarray of shape (n_cells,), optional
        Pseudotime values for ordering cells on trajectory.
    ax : matplotlib.axes.Axes, optional
        Existing axes. Creates new figure if None.
    figsize : tuple, default=(12, 10)
        Figure size.
    cmap : str, default='RdYlBu_r'
        Colormap for energy values.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes object with plot.
    """
    if embedding.shape[1] == 2:
        return _plot_energy_2d(energy, embedding, pseudotime, ax, figsize, cmap)
    elif embedding.shape[1] == 3:
        return _plot_energy_3d(energy, embedding, pseudotime, ax, figsize, cmap)
    else:
        raise ValueError("Embedding must be 2D or 3D")


def _plot_energy_2d(
    energy: np.ndarray,
    embedding: np.ndarray,
    pseudotime: Optional[np.ndarray] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[int, int] = (12, 10),
    cmap: str = 'RdYlBu_r'
) -> plt.Axes:
    """Plot 2D energy landscape."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=energy,
        s=50,
        alpha=0.6,
        cmap=cmap,
        edgecolors='black',
        linewidth=0.5
    )

    # Plot trajectory if pseudotime provided
    if pseudotime is not None:
        sort_idx = np.argsort(pseudotime)
        ax.plot(
            embedding[sort_idx, 0],
            embedding[sort_idx, 1],
            'k--',
            alpha=0.3,
            linewidth=1,
            label='Pseudotime trajectory'
        )
        ax.legend(fontsize=10)

    ax.set_xlabel('PC1 / UMAP1', fontsize=12)
    ax.set_ylabel('PC2 / UMAP2', fontsize=12)
    ax.set_title('Waddington Energy Landscape', fontsize=14, fontweight='bold')

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Energy', fontsize=11)

    return ax


def _plot_energy_3d(
    energy: np.ndarray,
    embedding: np.ndarray,
    pseudotime: Optional[np.ndarray] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[int, int] = (14, 10),
    cmap: str = 'RdYlBu_r'
) -> plt.Axes:
    """Plot 3D energy landscape."""
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

    scatter = ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        embedding[:, 2],
        c=energy,
        s=50,
        alpha=0.6,
        cmap=cmap,
        edgecolors='black',
        linewidth=0.5
    )

    # Plot trajectory if pseudotime provided
    if pseudotime is not None:
        sort_idx = np.argsort(pseudotime)
        ax.plot(
            embedding[sort_idx, 0],
            embedding[sort_idx, 1],
            embedding[sort_idx, 2],
            'k--',
            alpha=0.3,
            linewidth=1
        )

    ax.set_xlabel('PC1', fontsize=11)
    ax.set_ylabel('PC2', fontsize=11)
    ax.set_zlabel('PC3', fontsize=11)
    ax.set_title('Waddington Energy Landscape (3D)', fontsize=14, fontweight='bold')

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Energy', fontsize=11)

    return ax


def plot_fate_probabilities(
    fate_probs: np.ndarray,
    embedding: np.ndarray,
    attractor_names: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (16, 4)
) -> np.ndarray:
    """
    Plot fate probabilities for each attractor.

    Parameters
    ----------
    fate_probs : np.ndarray of shape (n_cells, n_attractors)
        Fate probability matrix.
    embedding : np.ndarray of shape (n_cells, 2)
        2D cell coordinates (e.g., PCA or UMAP).
    attractor_names : list of str, optional
        Names for attractors. If None, uses "Attractor 1", etc.
    figsize : tuple, default=(16, 4)
        Figure size.

    Returns
    -------
    axes : np.ndarray of plt.Axes
        Array of subplot axes.
    """
    n_attractors = fate_probs.shape[1]

    if attractor_names is None:
        attractor_names = [f'Attractor {i+1}' for i in range(n_attractors)]

    fig, axes = plt.subplots(1, n_attractors, figsize=figsize)
    if n_attractors == 1:
        axes = np.array([axes])

    for ax_idx, ax in enumerate(axes):
        scatter = ax.scatter(
            embedding[:, 0],
            embedding[:, 1],
            c=fate_probs[:, ax_idx],
            s=50,
            alpha=0.7,
            cmap='viridis',
            edgecolors='black',
            linewidth=0.5,
            vmin=0,
            vmax=1
        )

        ax.set_xlabel('PC1 / UMAP1', fontsize=11)
        ax.set_ylabel('PC2 / UMAP2', fontsize=11)
        ax.set_title(attractor_names[ax_idx], fontsize=12, fontweight='bold')

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Fate probability', fontsize=10)

    plt.tight_layout()
    return axes


def plot_attractor_genes(
    attractor_gene_programs: np.ndarray,
    gene_names: Optional[np.ndarray] = None,
    top_n: int = 20,
    attractor_names: Optional[List[str]] = None,
    figsize: Tuple[int, int] = (14, 8),
    cmap: str = 'RdBu_r'
) -> plt.Axes:
    """
    Plot marker genes of attractor states.

    Parameters
    ----------
    attractor_gene_programs : np.ndarray of shape (n_attractors, n_genes)
        Mean gene expression per attractor.
    gene_names : np.ndarray of shape (n_genes,), optional
        Gene names. If None, uses integer indices.
    top_n : int, default=20
        Number of top variable genes to plot.
    attractor_names : list of str, optional
        Names for attractors.
    figsize : tuple, default=(14, 8)
        Figure size.
    cmap : str, default='RdBu_r'
        Colormap.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes object with plot.
    """
    n_attractors, n_genes = attractor_gene_programs.shape

    if gene_names is None:
        gene_names = np.arange(n_genes).astype(str)

    if attractor_names is None:
        attractor_names = [f'Attractor {i+1}' for i in range(n_attractors)]

    # Select top variable genes
    gene_var = attractor_gene_programs.std(axis=0)
    top_idx = np.argsort(gene_var)[-top_n:][::-1]

    program_top = attractor_gene_programs[:, top_idx]
    gene_names_top = gene_names[top_idx]

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        program_top.T,
        xticklabels=attractor_names,
        yticklabels=gene_names_top,
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


def plot_reservoir_dynamics(
    H: np.ndarray,
    pseudotime: Optional[np.ndarray] = None,
    n_dims: int = 3,
    figsize: Tuple[int, int] = (14, 10)
) -> plt.Axes:
    """
    Plot reservoir state dynamics.

    Parameters
    ----------
    H : np.ndarray of shape (n_cells, n_reservoir)
        Reservoir state matrix.
    pseudotime : np.ndarray of shape (n_cells,), optional
        Pseudotime for coloring trajectory.
    n_dims : int, default=3
        Number of dimensions to plot (2 or 3).
    figsize : tuple, default=(14, 10)
        Figure size.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes object with plot.
    """
    from sklearn.decomposition import PCA

    # Reduce to n_dims via PCA
    pca = PCA(n_components=n_dims)
    H_pca = pca.fit_transform(H)

    if n_dims == 2:
        fig, ax = plt.subplots(figsize=figsize)
        if pseudotime is not None:
            scatter = ax.scatter(
                H_pca[:, 0],
                H_pca[:, 1],
                c=pseudotime,
                s=50,
                alpha=0.6,
                cmap='viridis',
                edgecolors='black',
                linewidth=0.5
            )
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Pseudotime', fontsize=11)

            # Plot trajectory
            sort_idx = np.argsort(pseudotime)
            ax.plot(H_pca[sort_idx, 0], H_pca[sort_idx, 1], 'k--', alpha=0.3, linewidth=1)
        else:
            ax.scatter(
                H_pca[:, 0],
                H_pca[:, 1],
                s=50,
                alpha=0.6,
                edgecolors='black',
                linewidth=0.5
            )

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} var)', fontsize=11)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} var)', fontsize=11)

    elif n_dims == 3:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        if pseudotime is not None:
            scatter = ax.scatter(
                H_pca[:, 0],
                H_pca[:, 1],
                H_pca[:, 2],
                c=pseudotime,
                s=50,
                alpha=0.6,
                cmap='viridis',
                edgecolors='black',
                linewidth=0.5
            )
            cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
            cbar.set_label('Pseudotime', fontsize=11)

            # Plot trajectory
            sort_idx = np.argsort(pseudotime)
            ax.plot(H_pca[sort_idx, 0], H_pca[sort_idx, 1], H_pca[sort_idx, 2],
                   'k--', alpha=0.3, linewidth=1)
        else:
            ax.scatter(
                H_pca[:, 0],
                H_pca[:, 1],
                H_pca[:, 2],
                s=50,
                alpha=0.6,
                edgecolors='black',
                linewidth=0.5
            )

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} var)', fontsize=11)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} var)', fontsize=11)
        ax.set_zlabel(f'PC3 ({pca.explained_variance_ratio_[2]:.1%} var)', fontsize=11)

    ax.set_title('Reservoir State Dynamics', fontsize=14, fontweight='bold')
    plt.tight_layout()

    return ax
