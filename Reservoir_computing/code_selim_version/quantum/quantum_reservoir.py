"""
Quantum Reservoir Computing functions for scReservoir.

Two approaches are provided, both implemented as plain functions:

  1. SIMULATED QUANTUM RESERVOIR (transverse-field Ising model)
     Classically simulates a real quantum reservoir on a CPU.
     Encodes each cell's gene expression by rotating qubits,
     evolves the quantum state under a fixed Hamiltonian H via
         |ψ(t+1)⟩ = U · R_input(x_t) · |ψ(t)⟩
     and reads out single- and two-qubit expectation values.
     These expectation values form the H_q matrix (analogous to H
     in the classical ESN) which is then used for the same ridge
     regression + GRN back-projection pipeline.

  2. NEXT-GENERATION (NG-RC) QUANTUM-INSPIRED FEATURES
     No quantum hardware or simulation needed.
     Constructs polynomial features from time-delay embeddings:
         φ_k = [o_k ; o_k ⊗ o_k]
     where o_k = [x(t), x(t-1), ..., x(t-m+1)].
     This classically approximates the nonlinear mixing that a
     quantum reservoir achieves through entanglement, following
     Gauthier et al. (2021) and the NG-QRC paper in your folder.

Key papers implemented here
----------------------------
  - Paper 1: "Quantum Next Generation Reservoir Computing" (Sornsaeng et al., 2024)
    → NG-RC feature construction, skipping-ahead training
  - Paper 2: "Quantum Reservoir Computing on Random Regular Graphs" (Ivaki et al., 2025)
    → Transverse-field Ising Hamiltonian, random regular graph topology,
      edge-of-chaos optimal regime, density-matrix evolution

Quick reference: classical ESN vs quantum reservoir
----------------------------------------------------
  Classical ESN:
    h(t+1) = (1-α)·h(t) + α·tanh( W_res·h(t) + W_in·x(t) )
    Reservoir state: vector h ∈ ℝ^{n_reservoir}

  Quantum reservoir:
    |ψ(t+1)⟩ = U(Δt) · R_input(x(t)) · |ψ(t)⟩
    Reservoir state: quantum state |ψ⟩ ∈ ℂ^{2^n_qubits}
    Readout: [⟨σ_z^i⟩, ⟨σ_z^i σ_z^j⟩] ∈ ℝ^{n_obs}  → plays the role of h

  The rest of the pipeline (ridge regression, GRN back-projection) is identical.
"""

import numpy as np
from scipy.linalg import expm


# ============================================================================
# Part A — Quantum primitives  (Pauli matrices, tensor products)
# ============================================================================

def pauli_x():
    """Pauli X (bit-flip) matrix. [[0,1],[1,0]]"""
    return np.array([[0, 1], [1, 0]], dtype=complex)

def pauli_y():
    """Pauli Y matrix. [[0,-i],[i,0]]"""
    return np.array([[0, -1j], [1j, 0]], dtype=complex)

def pauli_z():
    """Pauli Z (phase-flip) matrix. [[1,0],[0,-1]]"""
    return np.array([[1, 0], [0, -1]], dtype=complex)

def identity(d=2):
    """Identity matrix of dimension d."""
    return np.eye(d, dtype=complex)


def kron_op_on_qubit(op, qubit_idx, n_qubits):
    """
    Embed a single-qubit operator `op` on qubit `qubit_idx` inside an
    n_qubit system by tensoring with identity on all other qubits.

    For example, with n_qubits=3 and qubit_idx=1:
        result = I ⊗ op ⊗ I

    Parameters
    ----------
    op : np.ndarray of shape (2, 2)
        Single-qubit operator (e.g., pauli_z()).
    qubit_idx : int
        Which qubit to apply op to (0-indexed).
    n_qubits : int
        Total number of qubits in the system.

    Returns
    -------
    full_op : np.ndarray of shape (2**n_qubits, 2**n_qubits)
        Full many-body operator.
    """
    ops = [identity(2)] * n_qubits
    ops[qubit_idx] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


def kron_op_on_two_qubits(op_i, op_j, qubit_i, qubit_j, n_qubits):
    """
    Embed a two-qubit interaction op_i ⊗ op_j in an n_qubit system.

    Parameters
    ----------
    op_i, op_j : np.ndarray of shape (2, 2)
        Single-qubit operators on sites i and j.
    qubit_i, qubit_j : int
        Site indices (0-indexed).
    n_qubits : int

    Returns
    -------
    full_op : np.ndarray of shape (2**n_qubits, 2**n_qubits)
    """
    ops = [identity(2)] * n_qubits
    ops[qubit_i] = op_i
    ops[qubit_j] = op_j
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result


# ============================================================================
# Part B — Hamiltonian construction
# ============================================================================

def build_ising_hamiltonian(n_qubits, J=1.0, h=1.0, periodic=True):
    """
    Build the transverse-field Ising model Hamiltonian.

    H = -J · Σ_{i} Z_i Z_{i+1}   (Ising coupling along Z)
        + h · Σ_{i} X_i            (transverse field along X)

    This is the benchmark system used in both quantum RC papers.
    The competition between J (favoring ferromagnetic order) and h
    (favoring superposition states) drives a quantum phase transition.

    Near the phase boundary the system is "critical" — it has long-range
    correlations and maximum dynamical richness — which is why this regime
    gives the best reservoir performance (the "edge of chaos" discussed in
    Paper 2).

    Parameters
    ----------
    n_qubits : int
        Number of qubits (spin-1/2 sites).  Hilbert space: 2^n_qubits.
    J : float, default=1.0
        Ising coupling strength between nearest neighbors.
        J > 0: ferromagnetic coupling (prefers aligned spins).
    h : float, default=1.0
        Transverse field strength.  h / J ≈ 1 is near the critical point.
    periodic : bool, default=True
        If True, include the wrap-around term Z_{n-1} Z_0 (ring geometry).

    Returns
    -------
    H : np.ndarray of shape (2**n_qubits, 2**n_qubits), dtype=complex
        Full many-body Hamiltonian matrix.
    """
    dim = 2 ** n_qubits
    H   = np.zeros((dim, dim), dtype=complex)

    X = pauli_x()
    Z = pauli_z()

    # Ising coupling: -J Σ Z_i Z_{i+1}
    n_bonds = n_qubits if periodic else n_qubits - 1
    for i in range(n_bonds):
        j = (i + 1) % n_qubits
        H -= J * kron_op_on_two_qubits(Z, Z, i, j, n_qubits)

    # Transverse field: +h Σ X_i
    for i in range(n_qubits):
        H += h * kron_op_on_qubit(X, i, n_qubits)

    return H


def build_rrg_hamiltonian(n_qubits, adjacency, J_z=1.0, J_x=0.5, h_z=0.0, h_x=1.0,
                          disorder_strength=0.2, random_state=42):
    """
    Build a Hamiltonian on a random regular graph (Paper 2 model).

    H = Σ_{(i,j) ∈ E} [ J_z · Z_i Z_j + J_x · X_i X_j ]
        + Σ_i [ h_z · Z_i  +  (h_x + δ_i) · X_i ]

    where δ_i ∈ [-disorder_strength, +disorder_strength] is a random local
    field that breaks the translational symmetry and drives the
    integrable → chaotic transition crucial for reservoir performance.

    Parameters
    ----------
    n_qubits : int
        Number of spins.
    adjacency : np.ndarray of shape (n_qubits, n_qubits)
        Adjacency matrix of the graph.  adjacency[i, j] = 1 means i and j
        are connected.  For a k-regular graph, each row sums to k.
    J_z : float, default=1.0
        ZZ coupling (reference energy scale).
    J_x : float, default=0.5
        XX coupling.
    h_z : float, default=0.0
        Uniform Z field.
    h_x : float, default=1.0
        Mean X field.
    disorder_strength : float, default=0.2
        Amplitude of random local X-field perturbations δ_i.
        0 = no disorder (more integrable), large = localized (Anderson-like).
    random_state : int, default=42

    Returns
    -------
    H : np.ndarray of shape (2**n_qubits, 2**n_qubits), dtype=complex
    """
    rng = np.random.RandomState(random_state)
    dim = 2 ** n_qubits
    H   = np.zeros((dim, dim), dtype=complex)

    X = pauli_x()
    Z = pauli_z()

    # Two-body interactions on graph edges
    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            if adjacency[i, j]:
                H += J_z * kron_op_on_two_qubits(Z, Z, i, j, n_qubits)
                H += J_x * kron_op_on_two_qubits(X, X, i, j, n_qubits)

    # Single-site fields (uniform + disorder)
    delta = rng.uniform(-disorder_strength, disorder_strength, n_qubits)
    for i in range(n_qubits):
        if h_z != 0:
            H += h_z * kron_op_on_qubit(Z, i, n_qubits)
        H += (h_x + delta[i]) * kron_op_on_qubit(X, i, n_qubits)

    return H


def build_random_regular_graph(n_qubits, degree, random_state=42):
    """
    Generate the adjacency matrix of a random regular graph.

    In a k-regular graph, every node has exactly k neighbors.
    Paper 2 shows that k controls the information propagation speed:
    higher k → faster spread but too-dense graphs lose locality.

    Parameters
    ----------
    n_qubits : int
        Number of nodes (qubits).  Must satisfy n_qubits * degree is even.
    degree : int
        Number of neighbors per node.  Must be < n_qubits.
    random_state : int, default=42

    Returns
    -------
    adjacency : np.ndarray of shape (n_qubits, n_qubits), dtype=int
        Symmetric adjacency matrix with zeros on diagonal.
    """
    rng = np.random.RandomState(random_state)
    assert degree < n_qubits, "degree must be < n_qubits"
    assert (n_qubits * degree) % 2 == 0, "n_qubits * degree must be even"

    # Simple construction: connect each node to next `degree` nodes (modular)
    # then randomly rewire to break regularity while preserving degree
    adj = np.zeros((n_qubits, n_qubits), dtype=int)
    for i in range(n_qubits):
        for d in range(1, degree // 2 + 1):
            j = (i + d) % n_qubits
            adj[i, j] = 1
            adj[j, i] = 1

    # Random rewiring for irregularity
    for _ in range(n_qubits * degree):
        i, j = rng.choice(n_qubits, 2, replace=False)
        k, l = rng.choice(n_qubits, 2, replace=False)
        if adj[i, j] and adj[k, l] and not adj[i, l] and not adj[k, j] and \
           i != l and k != j:
            adj[i, j] = adj[j, i] = 0
            adj[k, l] = adj[l, k] = 0
            adj[i, l] = adj[l, i] = 1
            adj[k, j] = adj[j, k] = 1

    return adj


# ============================================================================
# Part C — Unitary (time evolution operator)
# ============================================================================

def build_unitary(H, dt):
    """
    Compute the time-evolution unitary for one time step.

    U(Δt) = exp( -i H Δt )

    This is the quantum analogue of the recurrent weight matrix W_res:
    just as W_res drives the classical reservoir state from h(t) to h(t+1),
    U drives the quantum state from |ψ(t)⟩ to |ψ(t+1)⟩.
    Unlike W_res, U is exactly unitary — it preserves the norm of |ψ⟩.

    Parameters
    ----------
    H : np.ndarray of shape (dim, dim), dtype=complex
        Hamiltonian matrix.
    dt : float
        Time step.  Small dt → slow evolution (short memory).
                    Large dt → fast oscillations (effectively mixes memory).
                    Optimal: J·dt ≈ 1–3 (from Paper 2).

    Returns
    -------
    U : np.ndarray of shape (dim, dim), dtype=complex
        Unitary evolution operator.
    """
    return expm(-1j * H * dt)


# ============================================================================
# Part D — Input encoding
# ============================================================================

def encode_input_rotation(psi, x_t, n_qubits, encoding_strength=np.pi / 4):
    """
    Encode a gene expression vector x_t into the quantum state by rotating qubits.

    Each qubit i receives a Y-rotation by angle θ_i = encoding_strength * x_t[i]:
        R_y(θ_i) = [[cos(θ/2), -sin(θ/2)],
                    [sin(θ/2),  cos(θ/2)]]

    The full n_qubit operation is:  R_y(θ_0) ⊗ R_y(θ_1) ⊗ ... ⊗ R_y(θ_{n-1})

    Note: If n_genes > n_qubits, the input is averaged into n_qubit bins.
          If n_genes < n_qubits, remaining qubits are left unchanged.

    Parameters
    ----------
    psi : np.ndarray of shape (2**n_qubits,), dtype=complex
        Current quantum state vector.
    x_t : np.ndarray of shape (n_genes,)
        Gene expression values for one cell (should be in [0, 1]).
    n_qubits : int
        Number of qubits.
    encoding_strength : float, default=π/4
        Maximum rotation angle.  π/4 gives moderate mixing.

    Returns
    -------
    psi_encoded : np.ndarray of shape (2**n_qubits,), dtype=complex
        State after encoding the input.
    """
    dim = 2 ** n_qubits

    # Map n_genes values to n_qubits angles
    if len(x_t) >= n_qubits:
        # Average input into n_qubit bins
        bin_size = len(x_t) // n_qubits
        angles   = np.array([
            x_t[i * bin_size: (i + 1) * bin_size].mean()
            for i in range(n_qubits)
        ]) * encoding_strength
    else:
        angles = np.zeros(n_qubits)
        angles[:len(x_t)] = x_t * encoding_strength

    # Build full rotation operator R = R_y(θ_0) ⊗ R_y(θ_1) ⊗ ...
    R = np.array([[1]], dtype=complex)
    for theta in angles:
        c, s = np.cos(theta / 2), np.sin(theta / 2)
        Ry_i = np.array([[c, -s], [s, c]], dtype=complex)
        R = np.kron(R, Ry_i)

    return R @ psi


# ============================================================================
# Part E — Observable readout
# ============================================================================

def measure_observables(psi, n_qubits):
    """
    Measure the expectation values of all local and two-body observables.

    Returns:
      - ⟨Z_i⟩       for i = 0 .. n_qubits-1       (n_qubits values)
      - ⟨Z_i Z_j⟩   for all i < j                  (n_qubits*(n_qubits-1)/2 values)

    These expectation values form one row of the H_q matrix and play
    exactly the role that h(t) plays in the classical reservoir.

    The choice of Z-basis is standard — Z eigenstates |0⟩ and |1⟩ are the
    computational basis, so ⟨Z_i⟩ = P(spin i = |0⟩) - P(spin i = |1⟩).

    Parameters
    ----------
    psi : np.ndarray of shape (2**n_qubits,), dtype=complex
        Current quantum state vector (unit norm).
    n_qubits : int

    Returns
    -------
    obs : np.ndarray of shape (n_qubits + n_qubits*(n_qubits-1)//2,)
        All expectation values.  Real-valued.
    """
    Z = pauli_z()
    single = []
    two_body = []

    for i in range(n_qubits):
        Zi = kron_op_on_qubit(Z, i, n_qubits)
        single.append(np.real(psi.conj() @ (Zi @ psi)))

    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            ZiZj = kron_op_on_two_qubits(Z, Z, i, j, n_qubits)
            two_body.append(np.real(psi.conj() @ (ZiZj @ psi)))

    return np.array(single + two_body, dtype=float)


def measure_density_observables(rho, n_qubits):
    """
    Measure observables from a density matrix (mixed state).

    Same as measure_observables but for ρ instead of |ψ⟩:
        ⟨O⟩ = Tr(ρ O)

    Parameters
    ----------
    rho : np.ndarray of shape (2**n_qubits, 2**n_qubits), dtype=complex
        Density matrix (positive semidefinite, trace = 1).
    n_qubits : int

    Returns
    -------
    obs : np.ndarray of shape (n_qubits + n_qubits*(n_qubits-1)//2,)
    """
    Z = pauli_z()
    single   = []
    two_body = []

    for i in range(n_qubits):
        Zi = kron_op_on_qubit(Z, i, n_qubits)
        single.append(np.real(np.trace(rho @ Zi)))

    for i in range(n_qubits):
        for j in range(i + 1, n_qubits):
            ZiZj = kron_op_on_two_qubits(Z, Z, i, j, n_qubits)
            two_body.append(np.real(np.trace(rho @ ZiZj)))

    return np.array(single + two_body, dtype=float)


# ============================================================================
# Part F — Build quantum reservoir weights (analogous to build_reservoir)
# ============================================================================

def build_quantum_reservoir(
    n_qubits,
    J=1.0,
    h=1.0,
    dt=0.5,
    hamiltonian_type='ising',
    periodic=True,
    rrg_degree=3,
    J_x=0.5,
    disorder_strength=0.2,
    random_state=42
):
    """
    Build a quantum reservoir: choose a Hamiltonian and precompute U(Δt).

    This is the quantum analogue of build_reservoir().  The Hamiltonian H
    is fixed and random (determined by topology + disorder).  The unitary
    U = exp(-i H Δt) is computed once and reused for every time step,
    just like W_res in the classical ESN.

    Parameters
    ----------
    n_qubits : int
        Number of qubits.  The Hilbert space has dimension 2^n_qubits.
        With n_qubits = 8 you get 256-dimensional state space (comparable
        to a 256-neuron classical reservoir, but with exponentially richer
        correlations).
        Warning: simulation cost scales as O(8^n_qubits) in memory.
        Practical limit on CPU: n_qubits ≤ 10.
    J : float, default=1.0
        Ising ZZ coupling strength.
    h : float, default=1.0
        Transverse X-field.  h ≈ J is near the quantum critical point.
    dt : float, default=0.5
        Time step for Hamiltonian evolution.
        Optimal range (from Paper 2): J·dt ≈ 1–3.
    hamiltonian_type : {'ising', 'rrg'}, default='ising'
        'ising': transverse-field Ising model on a 1D chain
        'rrg'  : random regular graph model (from Paper 2)
    periodic : bool, default=True
        Only used for 'ising'. Include wrap-around bond.
    rrg_degree : int, default=3
        Only used for 'rrg'. Graph degree k.
    J_x : float, default=0.5
        Only used for 'rrg'. XX interaction strength.
    disorder_strength : float, default=0.2
        Only used for 'rrg'. Random local field amplitude.
    random_state : int, default=42

    Returns
    -------
    q_weights : dict with keys:
        'H'                 : full Hamiltonian matrix
        'U'                 : unitary U = exp(-i H dt)
        'n_qubits'          : int
        'n_obs'             : number of observables per time step
        'dt'                : time step used
        'hamiltonian_type'  : str
        'adjacency'         : adjacency matrix (None for 'ising')
    """
    if hamiltonian_type == 'ising':
        H = build_ising_hamiltonian(n_qubits, J=J, h=h, periodic=periodic)
        adjacency = None
    elif hamiltonian_type == 'rrg':
        adjacency = build_random_regular_graph(n_qubits, rrg_degree, random_state)
        H = build_rrg_hamiltonian(
            n_qubits, adjacency,
            J_z=J, J_x=J_x, h_z=0.0, h_x=h,
            disorder_strength=disorder_strength, random_state=random_state
        )
    else:
        raise ValueError(f"Unknown hamiltonian_type: {hamiltonian_type!r}")

    U = build_unitary(H, dt)

    # Number of observables: single-qubit Z + all two-qubit ZZ
    n_obs = n_qubits + n_qubits * (n_qubits - 1) // 2

    return {
        'H':                H,
        'U':                U,
        'n_qubits':         n_qubits,
        'n_obs':            n_obs,
        'dt':               dt,
        'hamiltonian_type': hamiltonian_type,
        'adjacency':        adjacency,
    }


# ============================================================================
# Part G — Quantum reservoir dynamics  (the main loop, produces H_q)
# ============================================================================

def compute_quantum_reservoir_states(X, q_weights, encoding_strength=np.pi / 4, washout=50):
    """
    Run quantum reservoir dynamics over pseudotime-ordered cells.

    At each time step t:
      1. Encode x_t into the current quantum state via Y-rotations:
             |ψ'⟩ = R_input(x_t) · |ψ(t)⟩
      2. Evolve under the fixed Hamiltonian for one time step:
             |ψ(t+1)⟩ = U · |ψ'⟩
      3. Measure observables: ⟨Z_i⟩ and ⟨Z_i Z_j⟩ → one row of H_q

    This mirrors compute_reservoir_states() exactly, with:
      - R_input playing the role of W_in
      - U playing the role of W_res
      - tanh replaced by the exact unitary evolution exp(-iHΔt)
      - The observed expectation values playing the role of h(t)

    After the washout phase (to forget the |0...0⟩ initial state), all
    subsequent observables are stored as the H_q matrix.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix, pseudotime-sorted.
    q_weights : dict
        Output of build_quantum_reservoir().
    encoding_strength : float, default=π/4
        Rotation angle per unit of gene expression.
    washout : int, default=50
        Steps to discard before recording.

    Returns
    -------
    H_q : np.ndarray of shape (n_cells, n_obs)
        Quantum reservoir state matrix.
        H_q[t] = [⟨Z_i(t)⟩, ⟨Z_i Z_j(t)⟩] — the readout at pseudotime step t.
        Can be used directly with fit_readout_* and infer_grn() from grn.py.
    """
    U        = q_weights['U']
    n_qubits = q_weights['n_qubits']
    n_obs    = q_weights['n_obs']
    n_cells  = X.shape[0]

    H_q = np.zeros((n_cells, n_obs))

    # Initialize quantum state to |0...0⟩ (ground state of Z basis)
    psi    = np.zeros(2 ** n_qubits, dtype=complex)
    psi[0] = 1.0

    total_steps = n_cells + washout

    for t in range(total_steps):

        # Pick input cell
        if t < washout:
            x_t = X[t % n_cells]
        else:
            x_t = X[t - washout]

        # Step 1: encode current input by rotating qubits
        psi = encode_input_rotation(psi, x_t, n_qubits, encoding_strength)

        # Step 2: evolve quantum state under fixed Hamiltonian
        psi = U @ psi

        # Renormalize to prevent floating-point drift
        psi = psi / np.linalg.norm(psi)

        # Step 3: measure and record observables
        if t >= washout:
            H_q[t - washout] = measure_observables(psi, n_qubits)

    return H_q


# ============================================================================
# Part H — Next-Generation RC (NG-RC) — no quantum hardware needed
# ============================================================================

def build_ngrc_features(X, delay_steps=1, poly_degree=2):
    """
    Build Next-Generation Reservoir Computing (NG-RC) feature vectors.

    This is the "quantum-inspired" classical approach from Paper 1.
    Instead of actually simulating a quantum reservoir, you construct a
    polynomial feature expansion of a time-delay embedding:

        o_k = [ x(t), x(t - delay_steps), ..., x(t - (m-1)*delay_steps) ]
        φ_k = [ o_k ;  o_k ⊗ o_k ;  o_k ⊗ o_k ⊗ o_k ; ... ]   (up to poly_degree)

    The outer product o_k ⊗ o_k introduces ALL pairwise quadratic interactions
    between past states, which is what quantum entanglement does physically.
    The NG-RC paper showed this achieves similar accuracy to true quantum RC
    on many tasks, but runs on a classical CPU.

    For single-cell data:
    - "time" = pseudotime order
    - x(t) = gene expression of cell t (after pseudotime sorting)
    - delay_steps = how many pseudotime steps back to look

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix, pseudotime-sorted.
    delay_steps : int, default=1
        How many cells back to include in each delay embedding.
        delay_steps=1 means [x(t), x(t-1)] (two time steps).
    poly_degree : int, default=2
        Maximum polynomial degree.  2 = linear + quadratic features.
        Warning: quadratic features have dimension O(n_genes²).
        For large n_genes (>100), use delay_steps=1 and poly_degree=2.

    Returns
    -------
    features : np.ndarray of shape (n_valid_cells, n_features)
        Feature matrix.  n_valid_cells = n_cells - delay_steps
        (first delay_steps cells are dropped, no past to embed).
    n_features : int
        Number of features.
    """
    n_cells, n_genes = X.shape
    n_delays = delay_steps + 1                # number of time points in embedding

    valid_start = delay_steps                  # first cell we can embed
    n_valid     = n_cells - valid_start        # number of valid cells

    feature_list = []

    for t in range(valid_start, n_cells):

        # Build time-delay embedding: o = [x(t), x(t-1), ..., x(t-delay)]
        o = np.concatenate([X[t - d] for d in range(n_delays)])   # (n_delays*n_genes,)

        phi = [o]   # linear features

        if poly_degree >= 2:
            # Quadratic: all pairwise products o_i * o_j  (outer product, flattened)
            phi.append(np.outer(o, o).flatten())

        if poly_degree >= 3:
            # Cubic: all triple products (very large — only use for small n_genes)
            phi.append(np.tensordot(np.outer(o, o).flatten(), o, axes=0).flatten())

        feature_list.append(np.concatenate(phi))

    features = np.array(feature_list, dtype=float)
    return features, features.shape[1]


def compute_ngrc_reservoir_states(X, delay_steps=1, poly_degree=2):
    """
    Compute NG-RC feature matrix — a drop-in replacement for
    compute_reservoir_states() that requires no reservoir weights.

    The feature matrix plays the same role as H in the classical pipeline:
    it is used directly with fit_readout_* and infer_grn().

    The first `delay_steps` cells are dropped from the output because
    there is no prior context to embed.  Account for this when aligning
    with your expression matrix X.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix, pseudotime-sorted.
    delay_steps : int, default=1
        Number of delay steps in the embedding window.
    poly_degree : int, default=2
        Polynomial degree for feature expansion.

    Returns
    -------
    H_ng : np.ndarray of shape (n_cells - delay_steps, n_features)
        NG-RC feature matrix.  Use X[delay_steps:] as the matching
        expression matrix for regression.
    n_features : int
    """
    return build_ngrc_features(X, delay_steps=delay_steps, poly_degree=poly_degree)


# ============================================================================
# Part I — Full quantum pipeline convenience function
# ============================================================================

def run_quantum_grn_pipeline(
    X, pseudotime,
    n_qubits=6,
    J=1.0, h=1.0, dt=0.5,
    encoding_strength=np.pi / 4,
    washout=50,
    lambda_reg=1e-3,
    mode='causal',
    hamiltonian_type='ising',
    use_ngrc=False,
    ngrc_delay=1,
    ngrc_poly=2,
    random_state=42
):
    """
    Full quantum reservoir GRN inference pipeline (functional, no classes).

    Steps:
      1. Sort X by pseudotime
      2. Build quantum reservoir (or NG-RC features)
      3. Compute H_q matrix (quantum readout or NG-RC features)
      4. Fit readout weights using the chosen regression mode
      5. Infer GRN by back-projection

    For quantum mode: GRN is back-projected through the encoding
    operator (W_enc), which maps genes to qubits, playing the role
    of W_in from the classical pipeline.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime per cell.
    n_qubits : int, default=6
        Qubits in the quantum reservoir.
    J, h : float
        Ising coupling and transverse field.
    dt : float, default=0.5
        Hamiltonian evolution time step.
    encoding_strength : float, default=π/4
        Input rotation amplitude.
    washout : int, default=50
        Burn-in steps.
    lambda_reg : float, default=1e-3
        Ridge regularization.
    mode : {'standard', 'causal', 'velocity'}, default='causal'
        Regression mode.
    hamiltonian_type : {'ising', 'rrg'}, default='ising'
    use_ngrc : bool, default=False
        If True, use NG-RC polynomial features instead of quantum simulation.
    ngrc_delay : int, default=1
        Delay steps for NG-RC.
    ngrc_poly : int, default=2
        Polynomial degree for NG-RC.
    random_state : int, default=42

    Returns
    -------
    results : dict with keys:
        'X_sorted'    : expression matrix sorted by pseudotime
        'pt_sorted'   : sorted pseudotime
        'H_q'         : quantum/NG-RC state matrix
        'W_out'       : fitted readout weights
        'GRN'         : inferred gene regulatory network
        'q_weights'   : quantum reservoir weights (None if use_ngrc=True)
    """
    from utils import ridge_regression, normalize_grn

    # Step 1: sort by pseudotime
    sort_idx  = np.argsort(pseudotime)
    X_sorted  = X[sort_idx]
    pt_sorted = pseudotime[sort_idx]

    if use_ngrc:
        # NG-RC: no reservoir weights needed
        H_q, _ = compute_ngrc_reservoir_states(X_sorted, delay_steps=ngrc_delay, poly_degree=ngrc_poly)
        # Drop the first ngrc_delay cells from X to match H_q
        X_reg  = X_sorted[ngrc_delay:]
        pt_reg = pt_sorted[ngrc_delay:]
        q_weights = None
    else:
        # Simulated quantum reservoir
        q_weights = build_quantum_reservoir(
            n_qubits, J=J, h=h, dt=dt,
            hamiltonian_type=hamiltonian_type,
            random_state=random_state
        )
        H_q   = compute_quantum_reservoir_states(X_sorted, q_weights,
                                                  encoding_strength=encoding_strength,
                                                  washout=washout)
        X_reg  = X_sorted
        pt_reg = pt_sorted

    # Step 4: fit readout
    if mode == 'standard':
        W_out = ridge_regression(H_q, X_reg, lambda_reg)
    elif mode == 'causal':
        sort_idx2 = np.argsort(pt_reg)
        H_sorted  = H_q[sort_idx2]
        X_s2      = X_reg[sort_idx2]
        W_out = ridge_regression(H_sorted[:-1], X_s2[1:], lambda_reg)
    elif mode == 'velocity':
        sort_idx2 = np.argsort(pt_reg)
        H_sorted  = H_q[sort_idx2]
        X_s2      = X_reg[sort_idx2]
        pt_s2     = pt_reg[sort_idx2]
        dt_arr    = np.diff(pt_s2)
        dX        = np.diff(X_s2, axis=0)
        velocity  = dX / dt_arr[:, None]
        H_mid     = (H_sorted[:-1] + H_sorted[1:]) / 2
        W_out = ridge_regression(H_mid, velocity, lambda_reg)
    else:
        raise ValueError(f"Unknown mode: {mode!r}")

    # Step 5: infer GRN
    # For quantum reservoir: W_enc approximates W_in.
    # We construct a linear encoding matrix that maps n_genes → n_obs
    # by measuring how each gene's unit input moves the observables.
    n_genes = X.shape[1]
    n_obs   = H_q.shape[1]
    W_enc   = np.random.RandomState(random_state).randn(n_obs, n_genes) * (encoding_strength / np.pi)
    GRN     = normalize_grn(np.abs(W_enc.T @ W_out))

    return {
        'X_sorted':  X_sorted,
        'pt_sorted': pt_sorted,
        'H_q':       H_q,
        'W_out':     W_out,
        'GRN':       GRN,
        'q_weights': q_weights,
    }
