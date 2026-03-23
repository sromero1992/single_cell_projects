"""
Example 3 — Quantum Reservoir GRN Inference (two approaches)

Shows both quantum approaches:
  A) Simulated quantum reservoir (transverse-field Ising model)
  B) NG-RC: Next-Generation quantum-inspired features (no quantum hardware)

  quantum/
    quantum_reservoir.py → build_quantum_reservoir, compute_quantum_reservoir_states,
                           compute_ngrc_reservoir_states, run_quantum_grn_pipeline
  classical/
    grn.py      → fit_readout_causal, infer_grn
    utils.py    → order_by_pseudotime, get_top_regulators
    plotting.py → plot_grn_heatmap

Key insight:
  The only difference from the classical pipeline is HOW the H matrix is built.
  Once you have H_q (from the quantum reservoir) or H_ng (from NG-RC),
  all downstream analysis is identical to the classical case.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'classical'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'quantum'))

import utils
import grn
import plotting
import quantum_reservoir as qr


def make_synthetic_data(n_cells=300, n_genes=40, n_tf=6, seed=42):
    rng = np.random.RandomState(seed)
    true_grn = np.zeros((n_genes, n_genes))
    for i in range(n_tf):
        true_grn[i, i] = 0.8
        for j in rng.choice(n_genes, rng.randint(3, 8), replace=False):
            true_grn[i, j] = rng.uniform(0.3, 1.0)

    X = np.zeros((n_cells, n_genes))
    x = rng.randn(n_genes) * 0.5
    for t in range(n_cells):
        x    = np.tanh(true_grn @ x) + rng.randn(n_genes) * 0.2
        X[t] = x

    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-6)
    gene_names        = np.array([f'Gene_{i}' for i in range(n_genes)])
    gene_names[:n_tf] = [f'TF_{i}' for i in range(n_tf)]
    return X, gene_names, true_grn


def main():
    print("=" * 65)
    print("Example 3 — Quantum Reservoir GRN Inference")
    print("=" * 65)

    X, gene_names, true_grn = make_synthetic_data()
    pseudotime = np.linspace(0, 1, X.shape[0])

    X_s, pt_s, _ = utils.order_by_pseudotime(X, pseudotime)

    # ==================================================================
    # Approach A: Simulated quantum reservoir (Ising model)
    # ==================================================================
    # NOTE: n_qubits=6 → 64-dimensional Hilbert space.
    # The quantum reservoir has 6+15=21 observables per step
    # (6 single-qubit + 15 two-qubit correlators).
    # This is the "H_q" matrix — directly replaces the classical H.
    # ==================================================================
    print("\n[A] Simulated quantum reservoir (n_qubits=6, Ising model)")

    q_weights = qr.build_quantum_reservoir(
        n_qubits=6,
        J=1.0,               # ZZ coupling
        h=1.0,               # transverse X field  (h/J ≈ 1 → near critical point)
        dt=0.5,              # time step: J·dt = 0.5 (within optimal range from Paper 2)
        hamiltonian_type='ising',
    )
    print(f"  Hilbert space dim: {2**q_weights['n_qubits']}")
    print(f"  Observables per step: {q_weights['n_obs']}")

    H_q = qr.compute_quantum_reservoir_states(
        X_s, q_weights,
        encoding_strength=np.pi / 4,
        washout=30
    )
    print(f"  H_q shape: {H_q.shape}   (n_cells × n_observables)")

    # Fit readout (causal mode)
    W_out_q = grn.fit_readout_causal(H_q, X_s, pt_s, lambda_reg=1e-3)

    # For GRN back-projection we build a random encoding matrix W_enc
    # that approximates how genes map to quantum observables.
    # (In a real quantum device, W_enc would be derived analytically
    # from the encoding gate derivatives.)
    rng    = np.random.RandomState(42)
    W_enc  = rng.randn(q_weights['n_obs'], X.shape[1]) * (np.pi / 4) / np.sqrt(q_weights['n_obs'])
    GRN_q  = utils.normalize_grn(np.abs(W_enc.T @ W_out_q))
    print(f"  GRN shape: {GRN_q.shape}")

    mse_q = grn.compute_prediction_mse(X_s, H_q, W_out_q)
    print(f"  Prediction MSE: {mse_q:.4f}")

    # ==================================================================
    # Approach B: NG-RC (Next-Generation quantum-inspired, no hardware)
    # ==================================================================
    # Feature vector φ_k = [x(t), x(t-1),  x(t)⊗x(t-1)]
    # The outer product ⊗ introduces ALL quadratic interactions between
    # current and past gene expression — mimicking quantum entanglement.
    # Feature dimension: 2*n_genes + (2*n_genes)^2 = 2*40 + 6400 = 6480
    # (reduce with fewer genes or lower poly_degree for large datasets)
    # ==================================================================
    print("\n[B] NG-RC quantum-inspired features (delay=1, poly_degree=2)")

    H_ng, n_feat = qr.compute_ngrc_reservoir_states(X_s, delay_steps=1, poly_degree=2)
    print(f"  H_ng shape: {H_ng.shape}   (n_valid_cells × n_features)")
    print(f"  Features: {n_feat}")

    # Match expression to valid cells
    X_ng = X_s[1:]     # drop first cell (no prior context)
    pt_ng = pt_s[1:]

    W_out_ng  = grn.fit_readout_standard(H_ng, X_ng, lambda_reg=1e-3)
    # Simple GRN from readout only (no W_in analog for NG-RC)
    W_ng_slice = W_out_ng[:X.shape[1], :]   # use first n_genes rows as proxy for W_in
    GRN_ng     = utils.normalize_grn(np.abs(W_ng_slice.T @ W_out_ng))

    mse_ng = grn.compute_prediction_mse(X_ng, H_ng, W_out_ng)
    print(f"  Prediction MSE: {mse_ng:.4f}")

    # ==================================================================
    # Comparison plot
    # ==================================================================
    top_q  = utils.get_top_regulators(GRN_q,  gene_names, k=5)
    top_ng = utils.get_top_regulators(GRN_ng, gene_names, k=5)
    print("\nTop regulators for TF_0 (first TF):")
    print(f"  Quantum:  {[r[0] for r in top_q['TF_0']]}")
    print(f"  NG-RC:    {[r[0] for r in top_ng['TF_0']]}")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plotting.plot_grn_heatmap(true_grn, gene_names, top_n=25, ax=axes[0])
    axes[0].set_title('True GRN', fontsize=13, fontweight='bold')
    plotting.plot_grn_heatmap(GRN_q,   gene_names, top_n=25, ax=axes[1])
    axes[1].set_title('Quantum ESN (Ising)', fontsize=13, fontweight='bold')
    plotting.plot_grn_heatmap(GRN_ng,  gene_names, top_n=25, ax=axes[2])
    axes[2].set_title('NG-RC (polynomial)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig('ex03_quantum_grn.png', dpi=150, bbox_inches='tight')
    print("\nSaved: ex03_quantum_grn.png")
    plt.show()


if __name__ == '__main__':
    main()
