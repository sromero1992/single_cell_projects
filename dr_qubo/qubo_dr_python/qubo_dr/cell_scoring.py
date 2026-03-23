# Author: Selim Romero, Texas A&M University
# UCell reference: Andreatta & Carmona, Comput Struct Biotechnol J, 2021
# doi: 10.1016/j.csbj.2021.06.043

"""
Cell-level pathway activity scoring using UCell rank-based methodology.

This module provides functions for computing UCell scores (a rank-based gene set
activity scoring method) and differential activity analysis between conditions.
UCell is particularly suited for scRNA-seq data with sparse expression matrices.
"""

import warnings
import numpy as np
from scipy import stats
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import decoupler as dc
    HAS_DECOUPLER = True
except ImportError:
    HAS_DECOUPLER = False


def compute_ucell_score(X, g, geneset, n_max_rank=1500, chunk_size=500):
    """
    Compute UCell rank-based gene set activity scores for each cell.

    UCell is a rank-based, non-parametric single-cell gene set scoring method
    that performs well on sparse expression matrices. The algorithm ranks genes
    within each cell and computes a Mann-Whitney U statistic based on the ranks
    of genes in the input set.

    Parameters
    ----------
    X : array-like, shape (N, G)
        Expression matrix (cells × genes). Can be sparse (scipy.sparse).
    g : array-like, shape (G,)
        Gene names (strings), must match input genes.
    geneset : array-like or set
        Gene names in the set (strings). Matching is case-insensitive.
    n_max_rank : int, default=1500
        Maximum rank to use for normalization. Reduces the effect of highly
        expressed genes and dropout. UCell parameter.
    chunk_size : int, default=500
        Number of cells to process in each chunk (for memory efficiency).

    Returns
    -------
    scores : np.ndarray, shape (N,)
        UCell score for each cell in [0, 1].

    Warnings
    --------
    UserWarning
        If no genes from geneset are found in g, or if n_set >= n_max_rank.

    Notes
    -----
    Algorithm outline:
    1. For each cell, rank genes by expression (1 = lowest, G = highest).
       Ties are handled with average rank (scipy.stats.rankdata, method='average').
    2. Flip ranks so highest expression → highest rank.
    3. Cap ranks at n_max_rank.
    4. Compute Mann-Whitney U statistic: U = R - n_set*(n_set+1)/2
       where R is sum of ranks for genes in geneset.
    5. Normalize: score = U / (n_set * (n_max_rank - n_set)), clamp to [0, 1].

    References
    ----------
    Andreatta, M., & Carmona, S. J. (2021). joint contribution of TCGA and
    single-cell RNA-seq enables separation of tumor- and stromal-derived
    signatures in human glioblastoma. Computational and Structural Biotechnology
    Journal, 19, 3798-3809.

    Examples
    --------
    >>> X = np.random.poisson(2, (1000, 2000))  # 1000 cells, 2000 genes
    >>> gene_names = [f'Gene_{i}' for i in range(2000)]
    >>> pathway = ['Gene_10', 'Gene_23', 'Gene_99', 'Gene_500']
    >>> scores = compute_ucell_score(X, gene_names, pathway, n_max_rank=1500)
    >>> print(scores.shape)
    (1000,)
    >>> print(scores.min(), scores.max())
    0.0 1.0
    """
    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    # Case-insensitive gene matching
    geneset_lower = set(s.lower() for s in geneset)
    g_lower = np.array([gg.lower() for gg in g])
    set_mask = np.isin(g_lower, geneset_lower)

    n_set = np.sum(set_mask)
    if n_set == 0:
        warnings.warn(f"No genes from geneset found in gene list. Returning zero scores.", UserWarning)
        return np.zeros(X.shape[0])

    if n_set >= n_max_rank:
        warnings.warn(f"Gene set size ({n_set}) >= n_max_rank ({n_max_rank}). "
                      "Results may be unreliable. Consider increasing n_max_rank.", UserWarning)

    N = X.shape[0]
    scores = np.zeros(N)

    # Process in chunks
    for start_idx in range(0, N, chunk_size):
        end_idx = min(start_idx + chunk_size, N)
        X_chunk = X[start_idx:end_idx, :]

        # Rank genes for each cell (scipy.stats.rankdata with average for ties)
        # Shape: (chunk_size, G)
        ranks = np.apply_along_axis(
            lambda row: stats.rankdata(row, method='average'),
            axis=1,
            arr=X_chunk
        )

        # Flip so highest expression = highest rank
        G = X_chunk.shape[1]
        ranks = (G + 1) - ranks

        # Cap at n_max_rank
        ranks = np.minimum(ranks, n_max_rank)

        # Extract ranks for genes in set
        set_ranks = ranks[:, set_mask]  # (chunk_size, n_set)

        # Sum of ranks: R
        R = np.sum(set_ranks, axis=1)

        # Mann-Whitney U statistic
        U = R - n_set * (n_set + 1) / 2

        # Normalize
        denom = n_set * (n_max_rank - n_set)
        if denom > 0:
            ucell_score_chunk = U / denom
        else:
            ucell_score_chunk = np.zeros_like(U)

        # Clamp to [0, 1]
        ucell_score_chunk = np.clip(ucell_score_chunk, 0.0, 1.0)

        scores[start_idx:end_idx] = ucell_score_chunk

    return scores


def compute_differential_activity(X, g, batch_id, geneset, condition_A='KO', condition_B='WT',
                                  n_max_rank=1500, return_stats=True):
    """
    Compute differential UCell pathway activity between two conditions.

    Computes UCell scores per cell within each condition, performs statistical
    testing, and derives a differential activity vector (Vdiff) suitable as a
    linear term in QUBO optimization.

    Parameters
    ----------
    X : array-like, shape (N, G)
        Expression matrix (cells × genes). Can be sparse.
    g : array-like, shape (G,)
        Gene names (strings).
    batch_id : array-like, shape (N,)
        Condition labels for each cell (strings).
    geneset : array-like or set
        Gene names in the pathway (strings). Case-insensitive matching.
    condition_A : str, default='KO'
        Condition A label (substring search in batch_id).
    condition_B : str, default='WT'
        Condition B label (substring search in batch_id).
    n_max_rank : int, default=1500
        UCell parameter (max rank cap).
    return_stats : bool, default=True
        If True, compute Mann-Whitney U test and AUROC.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'scores_A': np.ndarray, shape (N_A,)
            Per-cell UCell scores for condition A cells.
        - 'scores_B': np.ndarray, shape (N_B,)
            Per-cell UCell scores for condition B cells.
        - 'Vdiff': np.ndarray, shape (G,)
            Differential activity vector: v_B - v_A. Suitable for QUBO
            linear term. Genes with positive Vdiff are more active in B.
        - 'ucell_pval': float or None
            Two-tailed Mann-Whitney U test p-value (condition A vs B scores).
        - 'ucell_auc': float or None
            AUROC score (B > A): probability that a random B score exceeds
            a random A score.

    Notes
    -----
    Vdiff computation:
    1. Normalize UCell score vectors to unit L2 norm.
    2. Row-normalize expression within each condition (L2 norm per cell).
    3. Project normalized expression onto normalized activity: v = Z @ psi.
    4. Vdiff = v_B - v_A.

    This vector captures both the direction and magnitude of differential
    expression weighted by pathway activity in each condition.

    Examples
    --------
    >>> X = np.random.poisson(2, (500, 1000))
    >>> genes = [f'Gene_{i}' for i in range(1000)]
    >>> batch = ['KO'] * 250 + ['WT'] * 250
    >>> pathway = ['Gene_10', 'Gene_50', 'Gene_100']
    >>> result = compute_differential_activity(X, genes, batch, pathway)
    >>> print(result.keys())
    dict_keys(['scores_A', 'scores_B', 'Vdiff', 'ucell_pval', 'ucell_auc'])
    """
    batch_id = np.asarray(batch_id)
    X = np.asarray(X)

    # Extract cells for each condition (substring matching, case-insensitive)
    condition_A_mask = np.array([condition_A.lower() in str(b).lower() for b in batch_id])
    condition_B_mask = np.array([condition_B.lower() in str(b).lower() for b in batch_id])

    Xa = X[condition_A_mask, :]
    Xb = X[condition_B_mask, :]

    # Compute UCell scores
    scores_A = compute_ucell_score(Xa, g, geneset, n_max_rank=n_max_rank)
    scores_B = compute_ucell_score(Xb, g, geneset, n_max_rank=n_max_rank)

    # Compute Vdiff
    # Normalize score vectors
    norm_a = np.linalg.norm(scores_A)
    norm_b = np.linalg.norm(scores_B)

    if norm_a > 0:
        psi_a = scores_A / norm_a
    else:
        psi_a = np.zeros_like(scores_A)

    if norm_b > 0:
        psi_b = scores_B / norm_b
    else:
        psi_b = np.zeros_like(scores_B)

    # Row-normalize expression (L2 norm per cell, axis=1)
    Za = normalize(Xa, norm='l2', axis=1)
    Zb = normalize(Xb, norm='l2', axis=1)

    # Project onto normalized activity
    v_a = Za.T @ psi_a  # (G,)
    v_b = Zb.T @ psi_b  # (G,)

    Vdiff = v_b - v_a

    # Statistical tests
    ucell_pval = None
    ucell_auc = None

    if return_stats and len(scores_A) > 0 and len(scores_B) > 0:
        # Mann-Whitney U test (two-sided)
        stat, pval = mannwhitneyu(scores_A, scores_B, alternative='two-sided')
        ucell_pval = pval

        # AUROC: probability that B > A
        try:
            labels = np.concatenate([np.zeros(len(scores_A)), np.ones(len(scores_B))])
            scores = np.concatenate([scores_A, scores_B])
            ucell_auc = roc_auc_score(labels, scores)
        except Exception:
            ucell_auc = None

    return {
        'scores_A': scores_A,
        'scores_B': scores_B,
        'Vdiff': Vdiff,
        'ucell_pval': ucell_pval,
        'ucell_auc': ucell_auc,
    }


def plot_ucell_violin(scores_A, scores_B, condition_A='KO', condition_B='WT',
                      title='Pathway Activity (UCell)', ax=None):
    """
    Create a violin plot comparing UCell score distributions between conditions.

    Parameters
    ----------
    scores_A : array-like, shape (N_A,)
        UCell scores for condition A cells.
    scores_B : array-like, shape (N_B,)
        UCell scores for condition B cells.
    condition_A : str, default='KO'
        Label for condition A.
    condition_B : str, default='WT'
        Label for condition B.
    title : str, default='Pathway Activity (UCell)'
        Plot title.
    ax : matplotlib.axes.Axes, optional
        Axis object to plot on. If None, creates new figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axes object.

    Notes
    -----
    The plot includes:
    - Violin plots for each condition (violin color: coral/red for A, steelblue for B)
    - Individual points with jitter (stripplot style)
    - Mann-Whitney U test p-value annotation
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib and seaborn are required for plotting. "
                          "Install with: pip install matplotlib seaborn")

    scores_A = np.asarray(scores_A)
    scores_B = np.asarray(scores_B)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    # Prepare data for seaborn
    data_list = []
    for score in scores_A:
        data_list.append({'Score': score, 'Condition': condition_A})
    for score in scores_B:
        data_list.append({'Score': score, 'Condition': condition_B})

    import pandas as pd
    df = pd.DataFrame(data_list)

    # Violin plot
    parts = ax.violinplot([scores_A, scores_B], positions=[0, 1], widths=0.7,
                          showmeans=False, showmedians=False, showextrema=False)

    # Color violins
    for pc in parts['bodies']:
        pc.set_facecolor('coral')
        pc.set_alpha(0.6)
    parts['bodies'][0].set_facecolor('coral')
    parts['bodies'][1].set_facecolor('steelblue')

    # Add individual points with jitter
    jitter_a = np.random.normal(0, 0.04, len(scores_A))
    jitter_b = np.random.normal(1, 0.04, len(scores_B))

    ax.scatter(jitter_a, scores_A, alpha=0.4, s=30, color='darkred', label=condition_A)
    ax.scatter(jitter_b, scores_B, alpha=0.4, s=30, color='darkblue', label=condition_B)

    # Mann-Whitney U test
    if len(scores_A) > 0 and len(scores_B) > 0:
        stat, pval = mannwhitneyu(scores_A, scores_B, alternative='two-sided')

        # Format p-value
        if pval < 0.001:
            pval_text = '***'
        elif pval < 0.01:
            pval_text = '**'
        elif pval < 0.05:
            pval_text = '*'
        else:
            pval_text = 'n.s.'

        # Add annotation
        y_max = max(scores_A.max(), scores_B.max())
        ax.text(0.5, y_max * 1.05, pval_text, ha='center', fontsize=14, fontweight='bold')

    ax.set_xticks([0, 1])
    ax.set_xticklabels([condition_A, condition_B])
    ax.set_ylabel('UCell Score')
    ax.set_title(title)
    ax.set_ylim(bottom=-0.05)
    ax.grid(axis='y', alpha=0.3)

    return fig, ax
