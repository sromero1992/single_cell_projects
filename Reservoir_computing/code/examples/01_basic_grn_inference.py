"""
Example 1: Basic Gene Regulatory Network Inference

This example demonstrates the core workflow of scReservoir:
1. Load and preprocess single-cell RNA-seq data
2. Initialize and run reservoir computing
3. Infer gene regulatory network
4. Visualize top regulators and GRN structure
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import scReservoir as scr


def generate_synthetic_data(n_cells=500, n_genes=100, n_tf=10, seed=42):
    """
    Generate synthetic single-cell RNA-seq data with known GRN structure.

    Parameters
    ----------
    n_cells : int
        Number of cells.
    n_genes : int
        Number of genes.
    n_tf : int
        Number of transcription factors (true regulators).
    seed : int
        Random seed.

    Returns
    -------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    gene_names : np.ndarray of shape (n_genes,)
        Gene identifiers.
    true_grn : np.ndarray of shape (n_genes, n_genes)
        True underlying GRN (sparse).
    """
    rng = np.random.RandomState(seed)

    # Create sparse ground truth GRN
    # TFs (first n_tf genes) regulate target genes
    true_grn = np.zeros((n_genes, n_genes))

    # TF self-regulation
    for i in range(n_tf):
        true_grn[i, i] = 0.8

    # TF -> target gene regulations
    for i in range(n_tf):
        n_targets = rng.randint(3, 10)
        targets = rng.choice(n_genes, n_targets, replace=False)
        for j in targets:
            true_grn[i, j] = rng.uniform(0.3, 1.0)

    # Generate expression data from linear dynamics
    X = np.zeros((n_cells, n_genes))
    x = rng.randn(n_genes) * 0.5

    for t in range(n_cells):
        # Update dynamics: x' = tanh(GRN @ x) + noise
        x = np.tanh(true_grn @ x) + rng.randn(n_genes) * 0.3
        X[t] = x

    # Normalize to [0, 1]
    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-6)

    gene_names = np.array([f'Gene_{i}' for i in range(n_genes)])
    gene_names[:n_tf] = np.array([f'TF_{i}' for i in range(n_tf)])

    return X, gene_names, true_grn


def main():
    """Main workflow."""
    print("=" * 70)
    print("scReservoir Example 1: Basic GRN Inference")
    print("=" * 70)

    # Generate synthetic data
    print("\n1. Generating synthetic single-cell RNA-seq data...")
    X, gene_names, true_grn = generate_synthetic_data(
        n_cells=500, n_genes=100, n_tf=10, seed=42
    )
    print(f"   Data shape: {X.shape}")
    print(f"   Gene names: {gene_names[:15]} ...")

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
    H = reservoir.compute_states(X, washout=50)
    print(f"   Reservoir state shape: {H.shape}")

    # Initialize GRN inference
    print("\n3. Inferring gene regulatory network...")
    grn_model = scr.ScReservoirGRN(
        reservoir=reservoir,
        lambda_reg=1e-3,
        n_top_regulators=10
    )

    # Fit GRN
    grn_model.fit(X, H, mode='standard')
    print("   GRN inference complete")

    # Get inferred GRN
    GRN_inferred = grn_model.get_grn()
    print(f"   Inferred GRN shape: {GRN_inferred.shape}")
    print(f"   GRN sparsity: {(GRN_inferred == 0).sum() / GRN_inferred.size:.1%}")

    # Get top regulators
    print("\n4. Top regulators per gene:")
    regulators = grn_model.get_top_regulators(gene_names, k=5)

    for gene_idx in range(5):
        gene = gene_names[gene_idx]
        top_regs = regulators[gene]
        print(f"   {gene}: {', '.join([f'{r[0]} ({r[1]:.3f})' for r in top_regs])}")

    # Compute prediction error
    print("\n5. Model performance:")
    mse = grn_model.compute_prediction_error(X, H)
    print(f"   Mean squared prediction error: {mse:.4f}")

    # Plot inferred GRN
    print("\n6. Generating visualizations...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot inferred GRN heatmap
    ax = axes[0]
    scr.plot_grn_heatmap(GRN_inferred, gene_names, top_n=30, ax=ax)

    # Plot true GRN heatmap
    ax = axes[1]
    scr.plot_grn_heatmap(true_grn, gene_names, top_n=30, ax=ax)
    ax.set_title('True GRN (ground truth)', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('grn_comparison.png', dpi=150, bbox_inches='tight')
    print("   Saved: grn_comparison.png")

    # Plot reservoir dynamics
    print("\n7. Visualizing reservoir dynamics...")
    fig, ax = plt.subplots(figsize=(12, 8))
    scr.plot_reservoir_dynamics(H, n_dims=2, figsize=(10, 8))
    plt.savefig('reservoir_dynamics.png', dpi=150, bbox_inches='tight')
    print("   Saved: reservoir_dynamics.png")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Cells analyzed: {X.shape[0]}")
    print(f"Genes profiled: {X.shape[1]}")
    print(f"Reservoir neurons: {reservoir.n_reservoir}")
    print(f"Regularization (lambda): {grn_model.lambda_reg}")
    print(f"Prediction MSE: {mse:.4f}")
    print(f"Top regulator: {gene_names[np.argsort(GRN_inferred.sum(axis=1))[-1]]}")
    print("=" * 70)

    plt.show()


if __name__ == '__main__':
    main()
