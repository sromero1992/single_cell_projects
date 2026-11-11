"""
scReservoir — Quantum Reservoir Computing module.

All functions are plain Python functions — no classes.

Two approaches are available:

  1. SIMULATED QUANTUM RESERVOIR
     Classically simulates a transverse-field Ising model or random regular
     graph model on a CPU.  Practical for n_qubits ≤ 10.

  2. NEXT-GENERATION (NG-RC) QUANTUM-INSPIRED
     No simulation needed.  Polynomial time-delay features approximate
     the nonlinear mixing of a quantum reservoir.  Scales to any dataset.

Typical usage — simulated quantum reservoir
--------------------------------------------
import sys; sys.path.insert(0, 'path/to/code/quantum')
from quantum_reservoir import (
    build_quantum_reservoir,
    compute_quantum_reservoir_states,
    run_quantum_grn_pipeline,
)
sys.path.insert(0, 'path/to/code/classical')
from grn import fit_readout_causal, infer_grn
from utils import order_by_pseudotime, get_top_regulators

X_s, pt_s, _ = order_by_pseudotime(X, pseudotime)
q_weights     = build_quantum_reservoir(n_qubits=6, J=1.0, h=1.0, dt=0.5)
H_q           = compute_quantum_reservoir_states(X_s, q_weights)
W_out         = fit_readout_causal(H_q, X_s, pt_s)
GRN           = infer_grn(q_weights['W_enc_approx'], W_out)   # see run_quantum_grn_pipeline

Typical usage — NG-RC (no quantum hardware)
---------------------------------------------
from quantum_reservoir import compute_ngrc_reservoir_states
from classical.grn import fit_readout_causal, infer_grn

H_ng, _ = compute_ngrc_reservoir_states(X_s, delay_steps=1, poly_degree=2)
# Use H_ng[1:] paired with X_s[1:] for regression (first cell has no delay)
"""
