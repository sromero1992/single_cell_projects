"""
Example 2: Pseudotime-based Causal GRN Inference

This example demonstrates causal GRN inference using pseudotime:
1. Load data with pseudotime information
2. Fit causal model using time-lagged reservoir states
3. Compare causal vs. standard GRN inference
4. Visualize directional regulatory relationships
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import scReservoir as scr


def generate_pseudotime_data(n_cells=300, n_genes=50, seed=42):
    """
    Generate pseudo-temporal single-cell data with causal structure.

    Parameters
    ----------
    n_cells : int
        Number of cells.
    n_genes : int
        Number of genes.
    seed : int
        Random seed.

    Returns
    -------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression along developmental trajectory.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime values.
    gene_names : np.ndarray of shape (n_genes,)
        Gene identifiers.
    """
    rng = np.random.RandomState(seed)

    # Create linear pseudotime
    pseudotime = np.sort(rng.uniform(0, 10, n_cells))

    # Generate data with causal structure
    X = np.zeros((n_cells, n_genes))

    # TF genes drive developmental progression
    for t_idx in range(n_cells):
        t = pseudotime[t_idx]

        # TF1: early activator
        X[t_idx, 0] = 5 * np.exp(-(t - 2)**2 / 2) + rng.randn() * 0.2

        # TF2: intermediate regulator
        X[t_idx, 1] = 4 * np.exp(-(t - 5)**2 / 2) + rng.randn() * 0.2

        # TF3: late activator
        X[t_idx, 2] = 5 * np.exp(-(t - 8)**2 / 2) + rng.randn() * 0.2

        # Target genes regulated by TFs
        for i in range(3, n_genes):
            # Each target is regulated by combination of TFs
            tf_contrib = (
                0.3 * X[t_idx, 0] +
                0.4 * X[t_idx, 1] +
                0.3 * X[t_idx, 2]
            )
            X[t_idx, i] = np.maximum(tf_contrib + rng.randn() * 0.3, 0)

    # Normalize
    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-6)

    gene_names = np.array([f'Gene_{i}' for i in range(n_genes)])
    gene_names[:3] = np.array(['TF1_early', 'TF2_mid', 'TF3_late'])

    return X, pseudotime, gene_names


def main():
    """Main workflow."""
    print("=" * 70)
    print("scReservoir Example 2: Pseudotime-based Causal GRN Inference")
    print("=" * 70)

    # Generate pseudo-temporal data
    print("\n1. Generating pseudo-temporal single-cell data...")
    X, pseudotime, gene_names = generate_pseudotime_data(
        n_cells=300, n_genes=50, seed=42
    )
    print(f"   Data shape: {X.shape}")
    print(f"   Pseudotime range: [{pseudotime.min():.2f}, {pseudotime.max():.2f}]")

    # Initialize reservoir
    print("\n2. Initializing reservoir computing model...")
    reservoir = scr.ScReservoir(
        n_reservoir=150,
        spectral_radius=0.9,
        input_scale=0.5,
        leak_rate=0.3,
        density=0.01,
        random_state=42
    )

    # Compute reservoir states
    print("   Computing reservoir states...")
    H = reservoir.compute_states(X, pseudotime=pseudotime, washout=30)
    print(f"   Reservoir state shape: {H.shape}")

    # Initialize GRN inference models
    print("\n3. Inferring gene regulatory networks...")

    # Standard GRN
    print("   Fitting STANDARD GRN model...")
    grn_standard = scr.ScReservoirGRN(
        reservoir=reservoir,
        lambda_reg=1e-3,
        n_top_regulators=10
    )
    grn_standard.fit(X, H, mode='standard')
    GRN_standard = grn_standard.get_grn()
    mse_standard = grn_standard.compute_prediction_error(X, H)

    # Causal GRN
    print("   Fitting CAUSAL GRN model...")
    grn_causal = scr.ScReservoirGRN(
        reservoir=reservoir,
        lambda_reg=1e-3,
        n_top_regulators=10
    )
    grn_causal.fit(X, H, pseudotime=pseudotime, mode='causal')
    GRN_causal = grn_causal.get_grn()
    mse_causal = grn_causal.compute_prediction_error(X[1:], H[:-1])

    # Velocity-based GRN
    print("   Fitting VELOCITY-based GRN model...")
    velocity = scr.compute_gene_velocity(X, pseudotime)

    grn_velocity = scr.ScReservoirGRN(
        reservoir=reservoir,
        lambda_reg=1e-3,
        n_top_regulators=10
    )
    grn_velocity.fit(X, H, pseudotime=pseudotime, mode='velocity')
    GRN_velocity = grn_velocity.get_grn()

    print("   GRN inference complete")

    # Compare models
    print("\n4. Model comparison:")
    print(f"   Standard MSE:  {mse_standard:.4f}")
    print(f"   Causal MSE:    {mse_causal:.4f}")

    # Top regulators
    print("\n5. Top TF regulators (causal model):")
    regulators_causal = grn_causal.get_top_regulators(gene_names, k=8)

    for tf_idx in range(3):
        tf_name = gene_names[tf_idx]
        targets = regulators_causal[tf_name]
        print(f"   {tf_name} regulates:")
        for target_name, influence in targets[:5]:
            print(f"      {target_name}: {influence:.3f}")

    # Analyze causal direction
    print("\n6. Causal directionality analysis:")
    print("   GRN sparsity comparison:")
    print(f"   Standard: {(GRN_standard == 0).sum() / GRN_standard.size:.1%} zeros")
    print(f"   Causal:   {(GRN_causal == 0).sum() / GRN_causal.size:.1%} zeros")
    print(f"   Velocity: {(GRN_velocity == 0).sum() / GRN_velocity.size:.1%} zeros")

    # Visualizations
    print("\n7. Generating visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Standard GRN
    scr.plot_grn_heatmap(GRN_standard, gene_names, top_n=25, ax=axes[0, 0])
    axes[0, 0].set_title('Standard GRN', fontsize=13, fontweight='bold')

    # Causal GRN
    scr.plot_grn_heatmap(GRN_causal, gene_names, top_n=25, ax=axes[0, 1])
    axes[0, 1].set_title('Causal GRN (time-lagged)', fontsize=13, fontweight='bold')

    # Velocity GRN
    scr.plot_grn_heatmap(GRN_velocity, gene_names, top_n=25, ax=axes[1, 0])
    axes[1, 0].set_title('Velocity GRN (dX/dt)', fontsize=13, fontweight='bold')

    # Gene expression along pseudotime
    ax = axes[1, 1]
    sort_idx = np.argsort(pseudotime)
    for i in range(3):  # Plot TFs
        ax.plot(pseudotime[sort_idx], X[sort_idx, i], label=gene_names[i], linewidth=2)
    ax.set_xlabel('Pseudotime', fontsize=12)
    ax.set_ylabel('Gene expression', fontsize=12)
    ax.set_title('TF Expression Dynamics', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('causal_grn_comparison.png', dpi=150, bbox_inches='tight')
    print("   Saved: causal_grn_comparison.png")

    # Plot velocity vs standard GRN differences
    fig, ax = plt.subplots(figsize=(10, 8))
    grn_diff = np.abs(GRN_causal - GRN_standard)
    ax.imshow(grn_diff[:30, :30], cmap='hot', aspect='auto')
    ax.set_xlabel('Target gene', fontsize=12)
    ax.set_ylabel('Regulator gene', fontsize=12)
    ax.set_title('Difference: Causal - Standard GRN', fontsize=14, fontweight='bold')
    plt.colorbar(ax.images[0], ax=ax, label='Absolute difference')
    plt.tight_layout()
    plt.savefig('grn_difference.png', dpi=150, bbox_inches='tight')
    print("   Saved: grn_difference.png")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Cells profiled: {X.shape[0]} along pseudotime trajectory")
    print(f"Genes: {X.shape[1]} (3 TFs + {X.shape[1]-3} targets)")
    print(f"Inference modes: Standard | Causal | Velocity")
    print(f"Causal improves prediction: {mse_standard > mse_causal}")
    print("=" * 70)

    plt.show()


if __name__ == '__main__':
    main()
