"""
Example 1 — Classical Reservoir GRN Inference

Shows the minimal end-to-end workflow using the classical reservoir module.

  classical/
    utils.py      → preprocess, order_by_pseudotime, get_top_regulators
    reservoir.py  → build_reservoir, compute_reservoir_states
    grn.py        → fit_readout_causal, infer_grn
    plotting.py   → plot_grn_heatmap, plot_reservoir_dynamics
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add classical module to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'classical'))
import utils
import reservoir
import grn
import plotting


# ---------------------------------------------------------------------------
# Synthetic data with known GRN
# ---------------------------------------------------------------------------

def make_synthetic_data(n_cells=500, n_genes=100, n_tf=10, seed=42):
    rng = np.random.RandomState(seed)
    true_grn = np.zeros((n_genes, n_genes))
    for i in range(n_tf):
        true_grn[i, i] = 0.8
        for j in rng.choice(n_genes, rng.randint(3, 10), replace=False):
            true_grn[i, j] = rng.uniform(0.3, 1.0)

    X = np.zeros((n_cells, n_genes))
    x = rng.randn(n_genes) * 0.5
    for t in range(n_cells):
        x    = np.tanh(true_grn @ x) + rng.randn(n_genes) * 0.3
        X[t] = x

    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-6)
    gene_names         = np.array([f'Gene_{i}' for i in range(n_genes)])
    gene_names[:n_tf]  = [f'TF_{i}' for i in range(n_tf)]
    return X, gene_names, true_grn


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Example 1 — Classical ESN GRN Inference")
    print("=" * 60)

    # Data
    X, gene_names, true_grn = make_synthetic_data()
    pseudotime = np.linspace(0, 1, X.shape[0])   # synthetic linear pseudotime

    # Sort by pseudotime
    X_s, pt_s, _ = utils.order_by_pseudotime(X, pseudotime)

    # Build reservoir  (fixed random weights, never trained)
    weights = reservoir.build_reservoir(
        n_reservoir=200, n_genes=X_s.shape[1],
        spectral_radius=0.9, input_scale=0.5, density=0.01
    )
    print(f"W_res: {weights['W_res'].shape}   W_in: {weights['W_in'].shape}")

    # Compute H matrix
    H = reservoir.compute_reservoir_states(X_s, weights, leak_rate=0.3, washout=50)
    print(f"H matrix: {H.shape}   (n_cells × n_reservoir)")

    # Fit readout weights — causal mode (H[t-1] → X[t])
    W_out = grn.fit_readout_causal(H, X_s, pt_s, lambda_reg=1e-3)
    print(f"W_out: {W_out.shape}   (n_reservoir × n_genes)")

    # Infer GRN
    GRN = grn.infer_grn(weights['W_in'], W_out)
    print(f"GRN: {GRN.shape}   sparsity: {(GRN == 0).mean():.1%}")

    # Top regulators
    top_regs = utils.get_top_regulators(GRN, gene_names, k=5)
    print("\nTop regulators for first 3 genes:")
    for g in gene_names[:3]:
        print(f"  {g:12s} ← {', '.join(f'{r[0]}({r[1]:.3f})' for r in top_regs[g])}")

    # Prediction error
    mse = grn.compute_prediction_mse(X_s, H, W_out)
    print(f"\nPrediction MSE: {mse:.4f}")

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    plotting.plot_grn_heatmap(GRN,      gene_names, top_n=30, ax=axes[0])
    plotting.plot_grn_heatmap(true_grn, gene_names, top_n=30, ax=axes[1])
    axes[1].set_title('True GRN', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('ex01_grn_comparison.png', dpi=150, bbox_inches='tight')
    print("\nSaved: ex01_grn_comparison.png")

    plotting.plot_reservoir_dynamics(H, pseudotime=pt_s, n_dims=2)
    plt.savefig('ex01_reservoir_dynamics.png', dpi=150, bbox_inches='tight')
    print("Saved: ex01_reservoir_dynamics.png")
    plt.show()


if __name__ == '__main__':
    main()
