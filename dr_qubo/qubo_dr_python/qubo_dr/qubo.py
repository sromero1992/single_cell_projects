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

# SteepestDescentSolver lives in dwave.samplers (classical, no QPU required)
try:
    from dwave.samplers import SteepestDescentSolver
    HAS_DWAVE_SAMPLERS = True
except ImportError:
    HAS_DWAVE_SAMPLERS = False


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


def assemble_qubo_matrix(dS, Vdiff, A_mnn_test, A_mnn_ref, K,
                         Xnet_target=None, penalty_scale=2.0,
                         auto_penalty=True):
    """
    Assemble the QUBO matrix with cardinality constraint.

    Constructs three nested matrices::

        Q0 = dS                                          (pure differential signal)
        Q1 = Q0 + diag(Vdiff) - (A_mnn_test + A_mnn_ref)         (+ Xnet_target if provided)
        Q  = Q1 + P * cardinality_constraint                (full penalised QUBO)

    then adds the cardinality penalty to enforce subnetwork size K.

    The penalty coefficient P is derived automatically from the spectral
    norm of Q1 (see ``compute_spectral_penalty``), giving a tighter and
    more principled bound than the elementwise-max heuristic.  Set
    ``auto_penalty=False`` to fall back to the legacy heuristic
    P = penalty_scale * max|Q1|.

    Parameters
    ----------
    dS : np.ndarray, shape (G, G)
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
    A_mnn_test : array-like or sparse, shape (G, G)
        MNN adjacency matrix for condition A (test).
    A_mnn_ref : array-like or sparse, shape (G, G)
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
        Final QUBO matrix with cardinality penalty (symmetric).
        This is what the solver minimises.
    Q1 : np.ndarray, shape (G, G)
        Unpenalised model matrix (dS + MNN + Vdiff [+ Xnet]).
        Useful for gene scoring, permutation tests, and comparison to Q0.
    Q0 : np.ndarray, shape (G, G)
        Pure differential co-expression (= dS, S_B - S_A).
        Contains only the biological signal — no graph adjacency,
        no pathway prior, no cardinality penalty.
        Use this matrix for direct visualisation and as the baseline
        signal before the other model terms are added.
    P : float
        Penalty coefficient used.
    """
    if issparse(A_mnn_test):
        A_mnn_test = A_mnn_test.toarray()
    if issparse(A_mnn_ref):
        A_mnn_ref = A_mnn_ref.toarray()

    # --- Q0: pure biological differential signal ---
    Q0 = dS.copy()

    # --- Q1: base QUBO matrix (no penalty yet) ---
    Q1 = dS + np.diag(Vdiff) - (A_mnn_test + A_mnn_ref)

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
    return Q, Q1, Q0, P


def solve_qubo_simulated_annealing(Q, num_reads=1000, seed=42, solver='sa'):
    """
    Solve QUBO by minimising E(z) = z^T Q z over z ∈ {0,1}^G.

    Parameters
    ----------
    Q : np.ndarray, shape (G, G)
        Symmetric QUBO matrix.
    num_reads : int
        Number of annealing runs (default 1000).
    seed : int
        Random seed (default 42).
    solver : str, {'sa', 'sd', 'fallback'}
        Solver backend:

        'sa'       — D-Wave neal SimulatedAnnealingSampler (default).
                     Stochastic, temperature-scheduled.  Best global search.
                     Requires: dwave-neal.
        'sd'       — SteepestDescentSolver (classical greedy descent).
                     Deterministic, very fast.  Fine for small K (≤ 30).
                     Equivalent to the classical optimizer in other QUBO
                     pipelines; useful when porting to R or non-D-Wave envs.
                     Requires: dwave-samplers.
        'fallback' — Pure NumPy SA.  No D-Wave dependency at all.
                     Use when neither dwave-neal nor dwave-samplers is
                     available, or to match R/MATLAB implementations exactly.

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
    solver_lower = solver.lower()

    if solver_lower == 'sa':
        if HAS_DIMOD:
            return _solve_dimod(Q, num_reads, seed)
        print("  [solver] dwave-neal not found — falling back to NumPy SA")
        return _solve_fallback(Q, num_reads, seed)

    elif solver_lower == 'sd':
        if HAS_DWAVE_SAMPLERS:
            return _solve_steepest_descent(Q, num_reads, seed)
        print("  [solver] dwave-samplers not found — falling back to NumPy SA")
        return _solve_fallback(Q, num_reads, seed)

    elif solver_lower == 'fallback':
        return _solve_fallback(Q, num_reads, seed)

    else:
        raise ValueError(
            f"Unknown solver '{solver}'. Choose 'sa', 'sd', or 'fallback'."
        )


# ------------------------------------------------------------------
# Internal solvers
# ------------------------------------------------------------------

def _qubo_to_dict(Q):
    """Convert a dense QUBO matrix to the {(i,j): value} dict that dimod samplers expect."""
    n = Q.shape[0]
    return {(i, j): float(Q[i, j])
            for i in range(n) for j in range(n)
            if Q[i, j] != 0}


def _solve_dimod(Q, num_reads, seed):
    """Solve via D-Wave neal SimulatedAnnealingSampler using sample_qubo (dict API).

    Uses sample_qubo() with a plain dict — no BinaryQuadraticModel construction,
    so it's version-agnostic and avoids the deprecated from_numpy_matrix/vartype API.
    """
    n = Q.shape[0]
    qubo = _qubo_to_dict(Q)
    sampler  = neal.SimulatedAnnealingSampler()
    response = sampler.sample_qubo(qubo, num_reads=num_reads, seed=seed,
                                   beta_range=(0.01, 3.0), num_sweeps=1000)
    best         = response.first.sample
    selected_idx = np.array([best[i] for i in range(n)], dtype=bool)
    best_energy  = response.first.energy
    all_samples  = np.array([[s[i] for i in range(n)] for s in response.samples()])
    return selected_idx, best_energy, all_samples


def _solve_steepest_descent(Q, num_reads, seed):
    """Solve via SteepestDescentSolver (classical greedy, dwave.samplers).

    Deterministic greedy descent — fast and dependency-light.
    Equivalent to the classical optimizer used in other QUBO pipelines.
    Works well for small K (≤ 30) where the landscape is not too rugged.
    num_reads initializations are run from random starts; the best is returned.
    """
    n    = Q.shape[0]
    qubo = _qubo_to_dict(Q)
    np.random.seed(seed)
    sampler  = SteepestDescentSolver()
    response = sampler.sample_qubo(qubo, num_reads=num_reads,
                                   initial_states_generator='random')
    best         = response.first.sample
    selected_idx = np.array([best[i] for i in range(n)], dtype=bool)
    best_energy  = response.first.energy
    all_samples  = np.array([[s[i] for i in range(n)] for s in response.samples()])
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

def extract_subnetwork(dS, Q1, selected_idx, z_vec=None):
    """
    Extract the differential co-expression submatrix and gene importance
    scores for the selected K-gene subnetwork.

    Parameters
    ----------
    dS : np.ndarray, shape (G, G)
    Q1 : np.ndarray, shape (G, G)
    selected_idx : np.ndarray, dtype bool, shape (G,)
        Boolean mask.  Internally converted to integer indices so behaviour
        is consistent across NumPy versions (np.ix_ with bool arrays is
        version-specific).
    z_vec : np.ndarray, optional, shape (G,)
        Float selection vector (defaults to selected_idx.astype(float)).

    Returns
    -------
    sub_Q0 : np.ndarray, shape (K, K)
        K×K submatrix of dS for the selected genes (Q0 restricted to hub).
        Sign is preserved (positive = reference gain, negative = disease gain)
        — same convention as dS itself, suitable for plot_gene_network.
    sub_Q1 : np.ndarray, shape (K, K)
        K×K submatrix of Q1 for the selected genes (optimisation landscape
        restricted to hub).
    sub_Qv : np.ndarray, shape (K,)
        Per-gene contribution scores (higher = more central in the hub).
    """
    if z_vec is None:
        z_vec = selected_idx.astype(float)
    # Use integer indices — avoids np.ix_ bool-array ambiguity across numpy versions
    sel_int   = np.where(selected_idx)[0]
    sub_Q0    = dS[np.ix_(sel_int, sel_int)]     # K×K dS (Q0)
    sub_Q1    = Q1[np.ix_(sel_int, sel_int)]         # K×K Q1
    sub_Qv    = -(Q1 @ z_vec)[selected_idx]          # per-gene score
    return sub_Q0, sub_Q1, sub_Qv


# ===========================================================================
#  Classical baseline solver
# ===========================================================================

def solve_classical_topk(dS, K, adjacency=None, method='rowsum'):
    """
    Select K genes classically as a baseline for comparison with the QUBO solution.

    No penalty terms, no annealing — pure algebraic ranking based on the
    differential co-expression matrix dS.

    Parameters
    ----------
    dS : np.ndarray, shape (G, G)
        Differential co-expression matrix (S_B - S_A).
        Negative values mean S_A > S_B, i.e. test/disease gains.
    K : int
        Number of genes to select.
    adjacency : np.ndarray, shape (G, G), optional
        MNN adjacency matrix.  When provided, scores are restricted to
        adjacent pairs (same as Q1 in the QUBO), giving a fairer comparison.
        If None, all gene pairs are used.
    method : str, {'rowsum', 'absrowsum', 'spectral'}
        Ranking criterion:

        'rowsum'    — row sum of –dS (= Q1 diagonal contribution):
                      score(i) = Σ_j (–dS[i,j]).  Selects genes whose
                      net co-expression increased in condition A.  Directly
                      comparable to the QUBO linear term.
        'absrowsum' — row sum of |dS|: score(i) = Σ_j |dS[i,j]|.
                      Selects genes with the most total co-expression change
                      (direction-agnostic).
        'spectral'  — leading eigenvector of –dS (or Q1 when adjacency
                      is provided).  Captures the dominant differential
                      co-expression pattern; equivalent to PCA-based ranking.

    Returns
    -------
    selected_idx : np.ndarray, dtype bool, shape (G,)
        Boolean mask of the K selected genes.
    scores : np.ndarray, shape (G,)
        Per-gene scores in the same order as the input genes.

    Notes
    -----
    This is the natural classical reference for the QUBO problem:
    instead of minimising z^T Q z subject to ||z||_0 = K (NP-hard),
    we rank genes by a simple linear or spectral score derived from Q1
    and take the top K.  The comparison highlights what the MNN graph
    adjacency constraint and the joint combinatorial optimisation add
    over a greedy ranking.
    """
    G = dS.shape[0]
    Q1 = -dS                             # Q1 = -dS: positive = disease gain

    if adjacency is not None:
        Q1_adj = Q1 * adjacency             # restrict to MNN-connected pairs
    else:
        Q1_adj = Q1

    if method == 'rowsum':
        scores = Q1_adj.sum(axis=1)         # net disease co-expression gain per gene
    elif method == 'absrowsum':
        scores = np.abs(Q1_adj).sum(axis=1) # total change magnitude per gene
    elif method == 'spectral':
        # Leading eigenvector of Q1_adj (symmetric) — sign convention: positive = gain
        eigvals, eigvecs = np.linalg.eigh(Q1_adj)
        scores = eigvecs[:, -1]             # eigenvector of largest eigenvalue
        if scores.sum() < 0:
            scores = -scores                # flip sign so higher = more disease gain
    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'rowsum', 'absrowsum', or 'spectral'.")

    K_safe = min(K, G)
    top_idx = np.argsort(scores)[-K_safe:][::-1]   # top-K indices, descending score
    selected_idx = np.zeros(G, dtype=bool)
    selected_idx[top_idx] = True
    return selected_idx, scores
