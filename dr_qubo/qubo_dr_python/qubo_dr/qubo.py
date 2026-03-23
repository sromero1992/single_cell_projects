# Author: Selim Romero, Texas A&M University
"""
QUBO assembly and solving module.

Constructs the QUBO matrix from differential co-expression, a general
cell-state difference vector, MNN adjacency priors, and an optional
pathway membership prior. Solves via simulated annealing (D-Wave or NumPy
fallback).

Key design note
---------------
The pathway network prior matrix (Xnet_target) is *optional*.

- When X is already restricted to pathway genes (standard use), every
  off-diagonal pair is within-pathway, so Xnet_target is uniformly -1
  everywhere and adds only a constant energy shift that does not affect
  which solution is optimal.  Pass Xnet_target=None (default) to skip it.

- Xnet_target is useful only when extra candidate genes (TF targets,
  GWAS hits, drug targets) are appended to the pathway gene set.  In
  that mixed case it assigns a -1 energy bonus specifically to core
  pathway-pathway gene pairs, steering the optimizer toward pathway-
  internal selections while still allowing the candidates to compete.
"""

import numpy as np
from scipy.sparse import issparse

# Try to import D-Wave tools
try:
    import dimod
    import neal
    HAS_DIMOD = True
except ImportError:
    HAS_DIMOD = False


def build_pathway_network_matrix(g, pathway_genes):
    """
    Build the optional pathway membership prior matrix.

    Creates a symmetric G × G matrix where off-diagonal entry (i, j) is -1
    when *both* genes belong to the core pathway gene set, and 0 otherwise.
    This is only informative when ``g`` contains extra candidate genes beyond
    the core pathway; when all genes in ``g`` are pathway members the matrix
    is -1 everywhere off-diagonal and acts as a constant shift.

    Parameters
    ----------
    g : list of str
        Gene names currently in the analysis (length G).
    pathway_genes : list of str
        Core pathway gene names (uppercase).  Should be a subset of g
        (uppercased).

    Returns
    -------
    Xnet : np.ndarray, shape (G, G)
        Pathway prior matrix.  Diagonal is 0.
    """
    n = len(g)
    Xnet = np.zeros((n, n))
    g_up = [gene.upper() for gene in g]
    pg_up = set(gene.upper() for gene in pathway_genes)
    is_pathway = np.array([gene in pg_up for gene in g_up])
    ix = np.where(is_pathway)[0]
    for i in ix:
        for j in ix:
            if i != j:
                Xnet[i, j] = -1.0
    return Xnet


def compute_spectral_penalty(Q1, penalty_scale=2.0):
    """
    Derive the cardinality penalty P from the spectral norm of Q1.

    The spectral norm  ||Q1||_2 = sigma_max(Q1)  equals the largest
    absolute eigenvalue for a symmetric matrix and provides a principled
    lower bound on P: any solution violating the |z| = K constraint
    incurs a penalty energy that exceeds the maximum biological energy
    gain (bounded by the spectral radius) when P >= ||Q1||_2.

    Setting  P = penalty_scale * ||Q1||_2  with penalty_scale > 1
    (default 2.0) adds a safety margin while keeping the penalty an
    order of magnitude tighter than the common heuristic
    P = 10 * max|Q1|_elementwise, which ignores matrix structure and
    over-penalises on the low-rank co-expression matrices typical of
    scRNA-seq data.

    Parameters
    ----------
    Q1 : np.ndarray, shape (G, G)
        Base QUBO matrix (symmetric).
    penalty_scale : float
        Safety multiplier applied on top of the spectral norm (default 2.0).

    Returns
    -------
    P : float
        Penalty coefficient.
    sigma_max : float
        Spectral norm of Q1 (= max singular value = max |eigenvalue|
        for symmetric Q1).
    """
    # np.linalg.norm(M, ord=2) computes the 2-norm = largest singular value
    sigma_max = np.linalg.norm(Q1, ord=2)
    P = penalty_scale * sigma_max
    return P, sigma_max


def assemble_qubo_matrix(Xdiff, Vdiff, Xmnn_A, Xmnn_B, K,
                         Xnet_target=None, penalty_scale=2.0,
                         auto_penalty=True):
    """
    Assemble the QUBO matrix with cardinality constraint.

    Constructs::

        Q1 = Xdiff + diag(Vdiff) - (Xmnn_A + Xmnn_B) [+ Xnet_target]

    then adds the cardinality penalty to enforce subnetwork size K.

    The penalty coefficient P is derived automatically from the spectral
    norm of Q1 (see ``compute_spectral_penalty``), giving a tighter and
    more principled bound than the elementwise-max heuristic.  Set
    ``auto_penalty=False`` to fall back to the legacy heuristic
    P = penalty_scale * max|Q1|.

    Parameters
    ----------
    Xdiff : np.ndarray, shape (G, G)
        Differential co-expression matrix (S_B - S_A).
        Negative off-diagonal entries indicate co-expression gain in
        condition A (test).
    Vdiff : np.ndarray, shape (G,)
        Gene-wise cell-state difference vector.  Each entry
        V_diff[g] = v_g^(B) - v_g^(A) is the difference in how strongly
        gene g's expression aligns with the per-cell biological state
        scalar between the two conditions.  The scalar can be any
        continuous per-cell measure: pathway activity (UCell), pseudotime,
        cell potency, differentiation rank, etc.  Pass np.zeros(G) to
        omit this term.
    Xmnn_A : array-like or sparse, shape (G, G)
        MNN adjacency matrix for condition A (test).
    Xmnn_B : array-like or sparse, shape (G, G)
        MNN adjacency matrix for condition B (reference / control).
    K : int
        Target subnetwork size.
    Xnet_target : np.ndarray, shape (G, G), optional
        Pathway membership prior.  Pass None (default) when all genes in
        the analysis are already restricted to the pathway gene set — in
        that case the prior is a uniform constant and has no effect on the
        solution.  Provide this matrix only when extra candidate genes
        (outside the core pathway) are included in the analysis.
    penalty_scale : float, optional
        Multiplier applied to the spectral norm (auto_penalty=True, default
        2.0) or to max|Q1| (auto_penalty=False, legacy default was 10.0).
    auto_penalty : bool, optional
        If True (default), set P = penalty_scale * ||Q1||_2 (spectral norm).
        If False, use legacy heuristic P = penalty_scale * max|Q1|.

    Returns
    -------
    Q : np.ndarray, shape (G, G)
        Final QUBO matrix (symmetric).
    Q1 : np.ndarray, shape (G, G)
        Unpenalised base matrix (useful for scoring and permutation tests).
    P : float
        Penalty coefficient used.
    """
    if issparse(Xmnn_A):
        Xmnn_A = Xmnn_A.toarray()
    if issparse(Xmnn_B):
        Xmnn_B = Xmnn_B.toarray()

    # --- Base QUBO matrix Q1 ---
    Q1 = Xdiff + np.diag(Vdiff) - (Xmnn_A + Xmnn_B)

    if Xnet_target is not None:
        Q1 = Q1 + Xnet_target

    # --- Cardinality penalty (spectral-norm derived or legacy) ---
    if auto_penalty:
        P, sigma_max = compute_spectral_penalty(Q1, penalty_scale)
        print(f"  Penalty (spectral): ||Q1||_2={sigma_max:.4f} → P={P:.4f} "
              f"(scale={penalty_scale})")
    else:
        P = penalty_scale * np.max(np.abs(Q1))
        print(f"  Penalty (legacy max): P={P:.4f} (scale={penalty_scale})")

    # Off-diagonal: +P ; Diagonal: +P*(1-2K)
    Q = Q1.copy()
    n = Q.shape[0]
    Q += P                              # add P to every element
    np.fill_diagonal(Q, Q1.diagonal() + P * (1 - 2 * K))  # fix diagonal

    # Symmetrise (numerical safety)
    Q = (Q + Q.T) / 2.0

    print(f"QUBO assembled: {n}×{n} | penalty P={P:.4f} | K={K}")
    return Q, Q1, P


def solve_qubo_simulated_annealing(Q, num_reads=1000, seed=42):
    """
    Solve QUBO by minimising E(z) = z^T Q z over z ∈ {0,1}^G.

    Uses D-Wave's neal (simulated annealing) when available, otherwise
    falls back to a pure NumPy implementation.

    Parameters
    ----------
    Q : np.ndarray, shape (G, G)
        Symmetric QUBO matrix.
    num_reads : int
        Number of annealing runs (default 1000).
    seed : int
        Random seed (default 42).

    Returns
    -------
    selected_idx : np.ndarray, dtype bool, shape (G,)
        Boolean mask of selected genes.
    best_energy : float
        Best objective value found.
    all_samples : np.ndarray, shape (num_reads, G)
        All binary solutions.
    """
    np.random.seed(seed)
    if HAS_DIMOD:
        return _solve_dimod(Q, num_reads, seed)
    return _solve_fallback(Q, num_reads, seed)


# ------------------------------------------------------------------
# Internal solvers
# ------------------------------------------------------------------

def _solve_dimod(Q, num_reads, seed):
    n = Q.shape[0]
    bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(Q, offset=0.0,
                                                        vartype=dimod.BINARY)
    sampler = neal.SimulatedAnnealingSampler()
    response = sampler.sample(bqm, num_reads=num_reads, seed=seed,
                               beta_range=(0.01, 3.0), num_sweeps=1000)
    best = response.first.sample
    selected_idx = np.array([best[i] for i in range(n)], dtype=bool)
    best_energy = response.first.energy
    all_samples = np.array([[s[i] for i in range(n)]
                             for s in response.samples()])
    return selected_idx, best_energy, all_samples


def _solve_fallback(Q, num_reads, seed):
    np.random.seed(seed)
    n = Q.shape[0]
    best_z = np.random.randint(0, 2, n)
    best_energy = _eval(best_z, Q)
    all_samples = [best_z.copy()]

    for k in range(num_reads):
        T = np.exp(-2.0 * k / num_reads)
        z = np.random.randint(0, 2, n)
        energy = _eval(z, Q)
        for _ in range(100):
            i = np.random.randint(0, n)
            z_new = z.copy(); z_new[i] ^= 1
            dE = _eval(z_new, Q) - energy
            if dE < 0 or np.random.rand() < np.exp(-dE / (T + 1e-10)):
                z, energy = z_new, energy + dE
        if energy < best_energy:
            best_z, best_energy = z.copy(), energy
        all_samples.append(z.copy())

    return best_z.astype(bool), best_energy, np.array(all_samples)


def _eval(z, Q):
    return z @ Q @ z


# ------------------------------------------------------------------

def extract_subnetwork(Xdiff, Q1, selected_idx, z_vec=None):
    """
    Extract the differential co-expression submatrix and gene importance
    scores for the selected K-gene subnetwork.

    Parameters
    ----------
    Xdiff : np.ndarray, shape (G, G)
    Q1 : np.ndarray, shape (G, G)
    selected_idx : np.ndarray, dtype bool, shape (G,)
    z_vec : np.ndarray, optional, shape (G,)
        Float selection vector (defaults to selected_idx.astype(float)).

    Returns
    -------
    sub_Q_net : np.ndarray, shape (K, K)
        Subnetwork differential co-expression (negated so positive = gain
        in condition A).
    sub_Qv : np.ndarray, shape (K,)
        Per-gene contribution scores (higher = more central).
    """
    if z_vec is None:
        z_vec = selected_idx.astype(float)
    sub_Q_net = -Xdiff[np.ix_(selected_idx, selected_idx)]
    sub_Qv = -(Q1 @ z_vec)[selected_idx]
    return sub_Q_net, sub_Qv
