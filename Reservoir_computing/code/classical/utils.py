"""
Utility functions for scReservoir.

Provides data preprocessing, velocity computation, and linear algebra routines.
All functions are stateless: they take inputs and return outputs, with no side effects.
"""

import numpy as np
from scipy import sparse


# ---------------------------------------------------------------------------
# Data loading / preprocessing
# ---------------------------------------------------------------------------

def preprocess(adata, n_hvg=2000, log_transform=True):
    """
    Preprocess single-cell data using scanpy conventions.

    Steps:
      1. Normalize each cell to 10,000 total counts
      2. log1p transform (if log_transform=True)
      3. Select top n_hvg highly variable genes
      4. Scale each gene to zero mean / unit variance (clip at 10)

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix with raw counts in .X
    n_hvg : int, default=2000
        Number of highly variable genes to retain.
    log_transform : bool, default=True
        Whether to apply log1p normalization (assume raw counts as input).

    Returns
    -------
    X : np.ndarray of shape (n_cells, n_hvg)
        Preprocessed, normalized expression matrix (float32).
    """
    import scanpy as sc

    adata_pp = adata.copy()

    if log_transform:
        sc.pp.normalize_total(adata_pp, target_sum=1e4)
        sc.pp.log1p(adata_pp)

    sc.pp.highly_variable_genes(adata_pp, n_top_genes=n_hvg)
    adata_pp = adata_pp[:, adata_pp.var.highly_variable]

    sc.pp.scale(adata_pp, max_value=10)

    X = adata_pp.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    return X.astype(np.float32)


def validate_expression_matrix(X):
    """
    Validate that an expression matrix is well-formed.

    Raises ValueError for:
      - Non-2D arrays
      - Zero-sized dimensions
      - NaN or Inf values

    Parameters
    ----------
    X : np.ndarray
        Expression matrix to validate.

    Returns
    -------
    True if valid (raises on failure).
    """
    if X.ndim != 2:
        raise ValueError("Expression matrix must be 2D (n_cells, n_genes).")
    if X.shape[0] == 0 or X.shape[1] == 0:
        raise ValueError("Expression matrix cannot have empty dimensions.")
    if not np.isfinite(X).all():
        raise ValueError("Expression matrix contains NaN or Inf values.")
    return True


# ---------------------------------------------------------------------------
# Pseudotime ordering
# ---------------------------------------------------------------------------

def order_by_pseudotime(X, pseudotime):
    """
    Sort cells by pseudotime values.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime value per cell.

    Returns
    -------
    X_sorted : np.ndarray of shape (n_cells, n_genes)
        Expression matrix with rows sorted by pseudotime (ascending).
    pt_sorted : np.ndarray of shape (n_cells,)
        Sorted pseudotime array.
    sort_idx : np.ndarray of shape (n_cells,)
        Integer indices used for sorting (useful to reorder other arrays).
    """
    sort_idx = np.argsort(pseudotime)
    return X[sort_idx], pseudotime[sort_idx], sort_idx


# ---------------------------------------------------------------------------
# Velocity computation
# ---------------------------------------------------------------------------

def compute_gene_velocity(X, pseudotime, window=5):
    """
    Estimate gene velocity (dX/dt) via finite differences with smoothing.

    Assumes X and pseudotime are already sorted in ascending pseudotime order.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix, pseudotime-sorted.
    pseudotime : np.ndarray of shape (n_cells,)
        Sorted pseudotime values.
    window : int, default=5
        Moving-average window for smoothing raw finite-difference estimates.

    Returns
    -------
    velocity : np.ndarray of shape (n_cells, n_genes)
        Estimated rate of change dX/dt per cell per gene.
    """
    n_cells, n_genes = X.shape

    dt = np.diff(pseudotime)          # shape (n_cells-1,)
    dX = np.diff(X, axis=0)           # shape (n_cells-1, n_genes)
    vel_raw = dX / dt[:, None]        # finite differences

    velocity = np.zeros((n_cells, n_genes))
    velocity[1:] = vel_raw            # first cell gets zero (no predecessor)

    # Smooth with a moving average
    if window > 1:
        for i in range(1, n_cells - 1):
            w_start = max(0, i - window // 2)
            w_end   = min(n_cells, i + window // 2 + 1)
            velocity[i] = velocity[w_start:w_end].mean(axis=0)

    return velocity


# ---------------------------------------------------------------------------
# Linear algebra helpers
# ---------------------------------------------------------------------------

def ridge_regression(H, y, lambda_reg=1e-3):
    """
    Solve ridge regression analytically.

    Minimizes:  ||H @ W - y||^2  +  lambda * ||W||^2

    Closed-form solution:
        W = (H^T H + lambda * I)^{-1}  H^T y

    Parameters
    ----------
    H : np.ndarray of shape (n_samples, n_features)
        Design matrix (typically the reservoir state matrix).
    y : np.ndarray of shape (n_samples,) or (n_samples, n_outputs)
        Target values (typically gene expression).
    lambda_reg : float, default=1e-3
        Regularization strength (larger = more shrinkage).

    Returns
    -------
    W : np.ndarray of shape (n_features,) or (n_features, n_outputs)
        Ridge regression coefficients.
    """
    n_features = H.shape[1]
    HTH = H.T @ H
    HTy = H.T @ y
    HTH_reg = HTH + lambda_reg * np.eye(n_features)
    W = np.linalg.solve(HTH_reg, HTy)
    return W


# ---------------------------------------------------------------------------
# GRN normalization helpers
# ---------------------------------------------------------------------------

def normalize_grn(GRN, method='minmax'):
    """
    Normalize a gene regulatory network matrix to [0, 1].

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Raw gene influence matrix (non-negative).
    method : {'minmax', 'zscore', 'log'}, default='minmax'
        Normalization strategy:
          - 'minmax' : scale to [0, 1] by (x - min) / (max - min)
          - 'zscore' : standardize then clip to [0, 1]
          - 'log'    : log1p then scale to [0, 1]

    Returns
    -------
    GRN_norm : np.ndarray of shape (n_genes, n_genes)
        Normalized values, clipped to [0, 1].
    """
    G = GRN.copy().astype(float)

    if method == 'minmax':
        lo, hi = G.min(), G.max()
        G = (G - lo) / (hi - lo) if hi > lo else np.zeros_like(G)

    elif method == 'zscore':
        mu, sigma = G.mean(), G.std()
        G = (G - mu) / sigma if sigma > 0 else np.zeros_like(G)

    elif method == 'log':
        G = np.log1p(G)
        hi = G.max()
        G = G / hi if hi > 0 else G

    return np.clip(G, 0, 1)


def threshold_grn(GRN, threshold=0.1):
    """
    Binarize a GRN by zeroing values below a threshold.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Normalized gene regulatory network.
    threshold : float, default=0.1
        Interactions with strength < threshold are set to 0.

    Returns
    -------
    GRN_sparse : np.ndarray of shape (n_genes, n_genes), dtype int
        Binary matrix (1 = regulation above threshold, 0 = otherwise).
    """
    return (GRN > threshold).astype(int)


# ---------------------------------------------------------------------------
# GRN summary statistics
# ---------------------------------------------------------------------------

def compute_sensitivity(GRN):
    """
    Compute per-gene in-degree (how strongly each gene is regulated by others).

    This is the normalized column sum of GRN.
    A high sensitivity score means many genes regulate this gene.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)

    Returns
    -------
    sensitivity : np.ndarray of shape (n_genes,)
    """
    n = GRN.shape[0]
    return GRN.sum(axis=0) / (n - 1) if n > 1 else GRN.sum(axis=0)


def compute_activity(GRN):
    """
    Compute per-gene out-degree (how strongly each gene regulates others).

    This is the normalized row sum of GRN.
    A high activity score suggests the gene is a master regulator / TF.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)

    Returns
    -------
    activity : np.ndarray of shape (n_genes,)
    """
    n = GRN.shape[1]
    return GRN.sum(axis=1) / (n - 1) if n > 1 else GRN.sum(axis=1)


def get_top_regulators(GRN, gene_names=None, k=10):
    """
    For each gene, return the top k regulators ranked by influence strength.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        GRN[i, j] = influence of gene i on gene j.
    gene_names : np.ndarray of shape (n_genes,), optional
        Gene labels. If None, uses integer indices as strings.
    k : int, default=10
        Number of top regulators to return per gene.

    Returns
    -------
    regulators : dict
        Keys are gene names (targets). Values are lists of
        (regulator_name, influence_score) tuples, sorted descending.
    """
    n_genes = GRN.shape[0]
    if gene_names is None:
        gene_names = np.arange(n_genes).astype(str)

    regulators = {}
    for j in range(n_genes):
        col = GRN[:, j]                         # influences flowing INTO gene j
        top_idx = np.argsort(col)[-k:][::-1]    # top-k descending
        regulators[gene_names[j]] = [(gene_names[i], col[i]) for i in top_idx]

    return regulators
