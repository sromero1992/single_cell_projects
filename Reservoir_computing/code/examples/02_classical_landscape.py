"""
Example 2 — Classical Reservoir Attractor Landscape Analysis

Shows the full landscape pipeline: energy landscape, attractor detection,
and fate probabilities using the classical reservoir module.

  classical/
    reservoir.py  → build_reservoir, compute_reservoir_states
    landscape.py  → run_landscape_pipeline
    plotting.py   → plot_energy_landscape, plot_fate_probabilities, plot_attractor_genes
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'classical'))
import utils
import reservoir
import landscape
import plotting


def make_branching_data(n_cells=600, n_genes=80, n_attractors=3, seed=42):
    """Synthetic data with branching developmental trajectory."""
    rng = np.random.RandomState(seed)

    # Three attractor states
    attractor_centers = rng.randn(n_attractors, n_genes)
    pseudotime = np.sort(rng.rand(n_cells))

    X = np.zeros((n_cells, n_genes))
    for t, pt in enumerate(pseudotime):
        # Cells gradually drift toward one of the attractors
        attractor_idx = rng.choice(n_attractors, p=[0.4, 0.35, 0.25])
        target = attractor_centers[attractor_idx]
        # Linear interpolation + noise
        X[t] = (1 - pt) * rng.randn(n_genes) * 0.5 + pt * target + rng.randn(n_genes) * 0.1

    return X, pseudotime


def main():
    print("=" * 60)
    print("Example 2 — Attractor Landscape Analysis")
    print("=" * 60)

    X, pseudotime = make_branching_data()

    # Sort by pseudotime
    X_s, pt_s, _ = utils.order_by_pseudotime(X, pseudotime)

    # Build and run classical reservoir
    weights = reservoir.build_reservoir(n_reservoir=300, n_genes=X_s.shape[1])
    H = reservoir.compute_reservoir_states(X_s, weights, leak_rate=0.3, washout=80)
    print(f"H matrix: {H.shape}")

    # Full landscape pipeline
    result = landscape.run_landscape_pipeline(
        X_s, H, pt_s,
        n_attractors=3,
        svd_rank=30,
        lambda_reg=1e-3
    )
    print(f"H_latent: {result['H_latent'].shape}")
    print(f"Attractors found: {len(result['attractors'])}")
    print(f"Energy range: [{result['energy'].min():.2f}, {result['energy'].max():.2f}]")

    # 2D PCA for visualization
    pca = PCA(n_components=2)
    emb = pca.fit_transform(X_s)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    plotting.plot_energy_landscape(result['energy'], emb, pseudotime=pt_s, ax=axes[0])
    axes[1].set_visible(False)
    plt.tight_layout()
    plt.savefig('ex02_energy_landscape.png', dpi=150, bbox_inches='tight')
    print("Saved: ex02_energy_landscape.png")

    axes_fate = plotting.plot_fate_probabilities(result['fate_probs'], emb)
    plt.savefig('ex02_fate_probs.png', dpi=150, bbox_inches='tight')
    print("Saved: ex02_fate_probs.png")

    # Attractor gene programs
    gene_names = np.array([f'Gene_{i}' for i in range(X.shape[1])])
    attractor_programs = np.array([
        X_s[cells].mean(axis=0) if len(cells) > 0 else np.zeros(X.shape[1])
        for cells in result['cells_per_attractor']
    ])
    plotting.plot_attractor_genes(attractor_programs, gene_names, top_n=20)
    plt.savefig('ex02_attractor_genes.png', dpi=150, bbox_inches='tight')
    print("Saved: ex02_attractor_genes.png")
    plt.show()


if __name__ == '__main__':
    main()
