# Author: Selim Romero, Texas A&M University
"""
Hub-level permutation test for QUBO differential co-expression significance.
"""

import numpy as np


def test_differential_hub(X_A, X_B, dS, selected_idx,
                          n_perm=200, seed=42, verbose=True):
    """
    Permutation test for differential co-expression hub significance.

    Tests whether the QUBO-selected K-gene hub shows more extreme aggregate
    differential co-expression than expected under random condition-label
    permutation.

    Test statistic — hub energy::

        E = z^T · dS · z  =  Σ_{i,j ∈ hub} dS[i,j]

    where z is the binary selection vector (z_i = 1 if gene i is in the hub).
    A very negative E means hub genes gained co-expression in the test/disease
    condition relative to reference/control.

    The null is built by randomly reassigning cells between conditions
    (preserving N_A and N_B), recomputing dS, and re-evaluating E.
    The p-value is two-sided (captures hubs differential in either direction).

    Note on z-score vs. per-gene stats
    -----------------------------------
    This produces **one z-score for the entire hub** — not a per-gene score.
    The hub is selected jointly (all K genes at once), so there is no
    independent per-gene null. The z-score is:

        z = (E_real − mean(E_null)) / std(E_null)

    It measures how many standard deviations the hub's aggregate differential
    signal is from random expectation.

    The ``n_sigma`` parameter in ``plot_hub_network`` derives a *display
    threshold* for individual edges from this same null:
        σ_edge = std(E_null) / (2 · √(K(K−1)/2))
    This estimates per-edge noise level assuming the hub energy variance
    distributes evenly across the K(K−1)/2 edge terms.  It is **not** a
    statistical test per edge — just a principled edge-visibility cutoff.

    Parameters
    ----------
    X_A : np.ndarray, shape (G, N_A)
        Library-size normalised expression for condition A (test).
    X_B : np.ndarray, shape (G, N_B)
        Same for condition B (reference/control).
    dS : np.ndarray, shape (G, G)
        Observed differential co-expression matrix (S_ref − S_test).
    selected_idx : np.ndarray, dtype bool, shape (G,)
        Boolean mask of QUBO-selected hub genes.
    n_perm : int
        Number of label permutations (default 200). Min resolvable p = 1/n_perm.
    seed : int
        Random seed (default 42).
    verbose : bool
        Print result summary (default True).

    Returns
    -------
    dict with keys:
        'pval'               – two-sided permutation p-value
        'zscore'             – z-score of E_real vs. null
        'E_real'             – observed hub energy (negative → test/disease hub)
        'E_null'             – null energies, shape (n_perm,)
        'null_mean'          – mean of null
        'null_std'           – std of null
        'significance_label' – '***' / '**' / '*' / 'n.s.'
    """
    from sklearn.preprocessing import normalize as _normalize

    X_A = np.asarray(X_A, dtype=float)
    X_B = np.asarray(X_B, dtype=float)
    z   = selected_idx.astype(float)
    G, N_A = X_A.shape
    N_B    = X_B.shape[1]
    N      = N_A + N_B

    E_real = float(z @ dS @ z)

    X_pool = np.concatenate([X_A, X_B], axis=1)
    rng    = np.random.default_rng(seed)
    E_null = np.zeros(n_perm)

    if verbose:
        print(f"  Hub permutation test: {n_perm} perms, "
              f"K={int(selected_idx.sum())} genes, N={N} cells…")

    for k in range(n_perm):
        perm  = rng.permutation(N)
        XAp   = X_pool[:, perm[:N_A]]
        XBp   = X_pool[:, perm[N_A:]]
        XAp_n = _normalize(XAp, norm='l2', axis=1)
        XBp_n = _normalize(XBp, norm='l2', axis=1)
        SAp   = np.clip(XAp_n @ XAp_n.T, 0, None)
        SBp   = np.clip(XBp_n @ XBp_n.T, 0, None)
        E_null[k] = float(z @ (SBp - SAp) @ z)

    null_mean = float(np.mean(E_null))
    null_std  = float(np.std(E_null))
    pval      = float(np.mean(np.abs(E_null) >= np.abs(E_real)))
    zscore    = (E_real - null_mean) / null_std if null_std > 0 else np.nan

    if   pval < 0.001: sig = '***'
    elif pval < 0.01:  sig = '**'
    elif pval < 0.05:  sig = '*'
    else:              sig = 'n.s.'

    if verbose:
        print(f"  E_real={E_real:.4f}  null: mean={null_mean:.4f}, std={null_std:.4f}")
        print(f"  z={zscore:.2f}  p={pval:.4f}  {sig}")

    return {
        'pval':               pval,
        'zscore':             zscore,
        'E_real':             E_real,
        'E_null':             E_null,
        'null_mean':          null_mean,
        'null_std':           null_std,
        'significance_label': sig,
    }
