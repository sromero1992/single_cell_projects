"""
Example 3: Attractor Landscape and Developmental Fates

This example demonstrates landscape-based analysis:
1. Load developmental single-cell data
2. Identify attractor states and developmental fates
3. Compute Waddington energy landscape
4. Analyze gene programs of attractor states
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import scReservoir as scr


def generate_developmental_data(n_cells=500, n_genes=100, n_attractors=3, seed=42):
    """
    Generate synthetic developmental single-cell data with multiple fates.

    Parameters
    ----------
    n_cells : int
        Number of cells.
    n_genes : int
        Number of genes.
    n_attractors : int
        Number of attractor states (fates).
    seed : int
        Random seed.

    Returns
    -------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime values.
    true_fates : np.ndarray of shape (n_cells,)
        True attractor assignment.
    gene_names : np.ndarray of shape (n_genes,)
        Gene names.
    """
    rng = np.random.RandomState(seed)

    # Create bifurcating trajectory
    X = np.zeros((n_cells, n_genes))
    pseudotime = np.linspace(0, 10, n_cells)
    true_fates = np.zeros(n_cells, dtype=int)

    # Attractor centers
    attractor_centers = rng.randn(n_attractors, n_genes)

    for i in range(n_cells):
        t = pseudotime[i]

        if t < 5:
            # Common progenitor phase
            base = rng.randn(n_genes) * 2
            X[i] = base
            true_fates[i] = 0
        else:
            # Differentiation: cells commit to different fates
            which_fate = int((i - n_cells * 0.5) / (n_cells * 0.5) * (n_attractors - 1)) + 1
            which_fate = min(which_fate, n_attractors - 1)
            true_fates[i] = which_fate

            # Cell state interpolates toward attractor
            progress = (t - 5) / 5  # 0 to 1
            X[i] = attractor_centers[which_fate] * progress + rng.randn(n_genes) * 0.5

    # Normalize
    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-6)
    X = X * 3  # Scale for better dynamics

    gene_names = np.array([f'Gene_{i}' for i in range(n_genes)])

    return X, pseudotime, true_fates, gene_names


def main():
    """Main workflow."""
    print("=" * 70)
    print("scReservoir Example 3: Attractor Landscape and Developmental Fates")
    print("=" * 70)

    # Generate developmental data
    print("\n1. Generating developmental single-cell data...")
    n_attractors = 3
    X, pseudotime, true_fates, gene_names = generate_developmental_data(
        n_cells=500, n_genes=100, n_attractors=n_attractors, seed=42
    )
    print(f"   Data shape: {X.shape}")
    print(f"   Number of fates: {n_attractors}")
    print(f"   Pseudotime range: [{pseudotime.min():.2f}, {pseudotime.max():.2f}]")
    print(f"   Cells per fate: {np.bincount(true_fates)}")

    # Initialize reservoir
    print("\n2. Initializing reservoir computing model...")
    reservoir = scr.ScReservoir(
        n_reservoir=200,
        spectral_radius=0.9,
        input_scale=0.5,
        leak_rate=0.3,
        density=0.01,
        random_state=42
    )

    # Compute reservoir states
    print("   Computing reservoir states...")
    H = reservoir.compute_states(X, pseudotime=pseudotime, washout=50)
    print(f"   Reservoir state shape: {H.shape}")

    # Fit attractor landscape
    print("\n3. Fitting attractor landscape...")
    landscape = scr.ScReservoirLandscape(
        n_attractors=n_attractors,
        svd_rank=30,
        lambda_reg=1e-3,
        random_state=42
    )

    landscape.fit(X, H, pseudotime)
    print("   Landscape fitting complete")

    # Get results
    energy = landscape.get_energy_landscape()
    fate_probs = landscape.get_fate_probabilities()
    H_latent = landscape.get_latent_states()
    attractors = landscape.get_attractor_states()

    print(f"   Energy landscape: {energy.shape}")
    print(f"   Fate probabilities: {fate_probs.shape}")
    print(f"   Latent states: {H_latent.shape}")
    print(f"   Attractors identified: {attractors.shape[0]}")

    # Analyze attractor programs
    print("\n4. Analyzing attractor gene programs...")
    attractor_programs = landscape.get_attractor_gene_programs(X)
    print(f"   Attractor programs shape: {attractor_programs.shape}")

    for i in range(n_attractors):
        top_genes_idx = np.argsort(attractor_programs[i])[-5:][::-1]
        top_genes = gene_names[top_genes_idx]
        print(f"   Fate {i}: top genes = {', '.join(top_genes)}")

    # Compute PCA for visualization
    print("\n5. Computing PCA for visualization...")
    pca = PCA(n_components=3, random_state=42)
    X_pca = pca.fit_transform(X)
    print(f"   PCA explained variance: {pca.explained_variance_ratio_[:2].sum():.1%}")

    # Predicted fates
    pred_fates = np.argmax(fate_probs, axis=1)
    accuracy = np.mean(pred_fates == true_fates)
    print(f"\n6. Prediction accuracy: {accuracy:.1%}")

    # Visualizations
    print("\n7. Generating landscape visualizations...")

    fig = plt.figure(figsize=(18, 12))

    # Energy landscape in PCA space
    ax1 = plt.subplot(2, 3, 1)
    scatter = ax1.scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=energy, s=50, alpha=0.6, cmap='RdYlBu_r', edgecolors='black', linewidth=0.5
    )
    sort_idx = np.argsort(pseudotime)
    ax1.plot(X_pca[sort_idx, 0], X_pca[sort_idx, 1], 'k--', alpha=0.3, linewidth=1)
    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=11)
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=11)
    ax1.set_title('Energy Landscape', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax1, label='Energy')

    # Pseudotime
    ax2 = plt.subplot(2, 3, 2)
    scatter = ax2.scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=pseudotime, s=50, alpha=0.6, cmap='viridis', edgecolors='black', linewidth=0.5
    )
    ax2.set_xlabel(f'PC1', fontsize=11)
    ax2.set_ylabel(f'PC2', fontsize=11)
    ax2.set_title('Pseudotime', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax2, label='Pseudotime')

    # True fates
    ax3 = plt.subplot(2, 3, 3)
    colors = ['red', 'blue', 'green', 'orange', 'purple'][:n_attractors]
    for fate in range(n_attractors):
        mask = true_fates == fate
        ax3.scatter(
            X_pca[mask, 0], X_pca[mask, 1],
            c=colors[fate], label=f'Fate {fate}', s=50, alpha=0.6, edgecolors='black', linewidth=0.5
        )
    ax3.set_xlabel(f'PC1', fontsize=11)
    ax3.set_ylabel(f'PC2', fontsize=11)
    ax3.set_title('True Fates', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)

    # Fate probabilities (Fate 0)
    ax4 = plt.subplot(2, 3, 4)
    scatter = ax4.scatter(
        X_pca[:, 0], X_pca[:, 1],
        c=fate_probs[:, 0], s=50, alpha=0.6, cmap='Reds', edgecolors='black', linewidth=0.5,
        vmin=0, vmax=1
    )
    ax4.set_xlabel(f'PC1', fontsize=11)
    ax4.set_ylabel(f'PC2', fontsize=11)
    ax4.set_title('Fate 0 Probability', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax4, label='Probability')

    # Predicted fates
    ax5 = plt.subplot(2, 3, 5)
    for fate in range(n_attractors):
        mask = pred_fates == fate
        ax5.scatter(
            X_pca[mask, 0], X_pca[mask, 1],
            c=colors[fate], label=f'Fate {fate}', s=50, alpha=0.6, edgecolors='black', linewidth=0.5
        )
    ax5.set_xlabel(f'PC1', fontsize=11)
    ax5.set_ylabel(f'PC2', fontsize=11)
    ax5.set_title('Predicted Fates', fontsize=12, fontweight='bold')
    ax5.legend(fontsize=10)

    # 3D energy landscape
    ax6 = plt.subplot(2, 3, 6, projection='3d')
    scatter = ax6.scatter(
        X_pca[:, 0], X_pca[:, 1], X_pca[:, 2],
        c=energy, s=30, alpha=0.6, cmap='RdYlBu_r', edgecolors='black', linewidth=0.3
    )
    ax6.set_xlabel('PC1', fontsize=10)
    ax6.set_ylabel('PC2', fontsize=10)
    ax6.set_zlabel('PC3', fontsize=10)
    ax6.set_title('3D Energy Landscape', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig('landscape_analysis.png', dpi=150, bbox_inches='tight')
    print("   Saved: landscape_analysis.png")

    # Attractor programs heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    scr.plot_attractor_genes(
        attractor_programs, gene_names, top_n=25,
        attractor_names=[f'Fate {i}' for i in range(n_attractors)], ax=ax
    )
    plt.tight_layout()
    plt.savefig('attractor_genes.png', dpi=150, bbox_inches='tight')
    print("   Saved: attractor_genes.png")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Cells analyzed: {X.shape[0]} along developmental trajectory")
    print(f"Genes: {X.shape[1]}")
    print(f"Attractor states identified: {n_attractors}")
    print(f"Fate prediction accuracy: {accuracy:.1%}")
    print(f"Mean energy: {energy.mean():.3f}")
    print(f"Energy range: [{energy.min():.3f}, {energy.max():.3f}]")
    print("=" * 70)

    plt.show()


if __name__ == '__main__':
    main()
