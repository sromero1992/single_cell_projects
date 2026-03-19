"""
Utility functions for scReservoir.

Provides data preprocessing, velocity computation, and linear algebra routines.
"""

from typing import Optional, Tuple
import numpy as np
from scipy import sparse


def preprocess(
    adata,
    n_hvg: int = 2000,
    log_transform: bool = True
) -> np.ndarray:
    """
    Preprocess single-cell data using scanpy conventions.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix with gene expression in .X
    n_hvg : int, default=2000
        Number of highly variable genes to select.
    log_transform : bool, default=True
        If True, applies log1p transformation (assumes raw counts).

    Returns
    -------
    X : np.ndarray of shape (n_cells, n_hvg)
        Preprocessed and normalized expression matrix.
    """
    import scanpy as sc

    # Make a copy
    adata_pp = adata.copy()

    # Log transform if needed
    if log_transform:
        sc.pp.normalize_total(adata_pp, target_sum=1e4)
        sc.pp.log1p(adata_pp)

    # Select highly variable genes
    sc.pp.highly_variable_genes(adata_pp, n_top_genes=n_hvg)
    adata_pp = adata_pp[:, adata_pp.var.highly_variable]

    # Standardize
    sc.pp.scale(adata_pp, max_value=10)

    X = adata_pp.X
    if hasattr(X, 'toarray'):  # sparse matrix
        X = X.toarray()

    return X.astype(np.float32)


def order_by_pseudotime(
    X: np.ndarray,
    pseudotime: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sort cells by pseudotime.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix.
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime values.

    Returns
    -------
    X_sorted : np.ndarray of shape (n_cells, n_genes)
        Expression matrix sorted by pseudotime.
    pseudotime_sorted : np.ndarray of shape (n_cells,)
        Sorted pseudotime values.
    """
    sort_idx = np.argsort(pseudotime)
    return X[sort_idx], pseudotime[sort_idx]


def compute_gene_velocity(
    X: np.ndarray,
    pseudotime: np.ndarray,
    window: int = 5
) -> np.ndarray:
    """
    Compute gene velocity (dX/dt) via finite differences.

    Parameters
    ----------
    X : np.ndarray of shape (n_cells, n_genes)
        Gene expression matrix (should be sorted by pseudotime).
    pseudotime : np.ndarray of shape (n_cells,)
        Pseudotime values (should be sorted).
    window : int, default=5
        Window size for smoothing velocity estimates.

    Returns
    -------
    velocity : np.ndarray of shape (n_cells, n_genes)
        Gene velocity estimates.
    """
    n_cells, n_genes = X.shape

    # Compute raw velocity via finite differences
    dt = np.diff(pseudotime)
    dX = np.diff(X, axis=0)
    vel_raw = dX / dt[:, None]

    # Pad first cell
    velocity = np.zeros((n_cells, n_genes))
    velocity[1:] = vel_raw

    # Smooth with moving average
    if window > 1:
        for i in range(1, n_cells - 1):
            w_start = max(0, i - window // 2)
            w_end = min(n_cells, i + window // 2 + 1)
            velocity[i] = velocity[w_start:w_end].mean(axis=0)

    return velocity


def sparse_reservoir(
    n_res: int,
    n_genes: int,
    density: float = 0.01,
    input_scale: float = 0.5,
    random_state: int = 42
) -> Tuple[sparse.csr_matrix, np.ndarray]:
    """
    Initialize sparse reservoir and input weight matrices.

    Parameters
    ----------
    n_res : int
        Number of reservoir neurons.
    n_genes : int
        Number of input genes.
    density : float, default=0.01
        Sparsity of W_res (fraction of non-zero entries).
    input_scale : float, default=0.5
        Scaling factor for input weights.
    random_state : int, default=42
        Random seed.

    Returns
    -------
    W_res : scipy.sparse.csr_matrix of shape (n_res, n_res)
        Sparse recurrent weight matrix.
    W_in : np.ndarray of shape (n_res, n_genes)
        Dense input weight matrix.
    """
    rng = np.random.RandomState(random_state)

    # Sparse recurrent matrix
    W_res = sparse.random(
        n_res, n_res, density=density, random_state=rng, format='csr'
    )

    # Normalize by spectral radius
    eigenvalues = sparse.linalg.eigsh(
        W_res.T @ W_res, k=1, which='LM', return_eigenvectors=False
    )
    rho = np.sqrt(eigenvalues[0])
    W_res *= 0.9 / (rho + 1e-10)

    # Input weights
    W_in = rng.randn(n_res, n_genes) * input_scale

    return W_res, W_in


def ridge_regression(
    H: np.ndarray,
    y: np.ndarray,
    lambda_reg: float = 1e-3
) -> np.ndarray:
    """
    Solve ridge regression: argmin_w ||H @ w - y||^2 + lambda * ||w||^2

    Parameters
    ----------
    H : np.ndarray of shape (n_samples, n_features)
        Design matrix (reservoir states).
    y : np.ndarray of shape (n_samples, n_outputs)
        Target matrix (gene expression).
    lambda_reg : float, default=1e-3
        Regularization strength.

    Returns
    -------
    w : np.ndarray of shape (n_features, n_outputs)
        Ridge regression coefficients.
    """
    n_features = H.shape[1]

    # Solve: (H.T @ H + lambda*I) @ w = H.T @ y
    HTH = H.T @ H
    HTy = H.T @ y

    # Add regularization
    HTH_reg = HTH + lambda_reg * np.eye(n_features)

    # Solve linear system
    w = np.linalg.solve(HTH_reg, HTy)

    return w


def normalize_grn(
    GRN: np.ndarray,
    method: str = 'minmax'
) -> np.ndarray:
    """
    Normalize gene regulatory network.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene influence matrix.
    method : {'minmax', 'zscore', 'log'}, default='minmax'
        Normalization method.

    Returns
    -------
    GRN_norm : np.ndarray of shape (n_genes, n_genes)
        Normalized GRN.
    """
    GRN_norm = GRN.copy()

    if method == 'minmax':
        GRN_min = GRN_norm.min()
        GRN_max = GRN_norm.max()
        if GRN_max > GRN_min:
            GRN_norm = (GRN_norm - GRN_min) / (GRN_max - GRN_min)
        else:
            GRN_norm = np.zeros_like(GRN_norm)

    elif method == 'zscore':
        mu = GRN_norm.mean()
        sigma = GRN_norm.std()
        if sigma > 0:
            GRN_norm = (GRN_norm - mu) / sigma
        else:
            GRN_norm = np.zeros_like(GRN_norm)

    elif method == 'log':
        GRN_norm = np.log1p(GRN_norm)
        GRN_max = GRN_norm.max()
        if GRN_max > 0:
            GRN_norm = GRN_norm / GRN_max

    return np.clip(GRN_norm, 0, 1)


def threshold_grn(
    GRN: np.ndarray,
    threshold: float = 0.1
) -> np.ndarray:
    """
    Binarize GRN by thresholding.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene regulatory network.
    threshold : float, default=0.1
        Threshold value.

    Returns
    -------
    GRN_binary : np.ndarray of shape (n_genes, n_genes)
        Binarized GRN (0 or 1).
    """
    return (GRN > threshold).astype(int)


def compute_sensitivity(
    GRN: np.ndarray,
    n_genes: Optional[int] = None
) -> np.ndarray:
    """
    Compute in-degree (sensitivity) of each gene.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene regulatory network.
    n_genes : int, optional
        Number of genes (defaults to GRN.shape[0]).

    Returns
    -------
    sensitivity : np.ndarray of shape (n_genes,)
        In-degree (number/strength of incoming regulations).
    """
    if n_genes is None:
        n_genes = GRN.shape[0]

    sensitivity = GRN.sum(axis=0)
    return sensitivity / (GRN.shape[0] - 1) if GRN.shape[0] > 1 else sensitivity


def compute_activity(
    GRN: np.ndarray,
    n_genes: Optional[int] = None
) -> np.ndarray:
    """
    Compute out-degree (activity) of each gene.

    Parameters
    ----------
    GRN : np.ndarray of shape (n_genes, n_genes)
        Gene regulatory network.
    n_genes : int, optional
        Number of genes (defaults to GRN.shape[0]).

    Returns
    -------
    activity : np.ndarray of shape (n_genes,)
        Out-degree (number/strength of outgoing regulations).
    """
    if n_genes is None:
        n_genes = GRN.shape[0]

    activity = GRN.sum(axis=1)
    return activity / (GRN.shape[1] - 1) if GRN.shape[1] > 1 else activity


def validate_expression_matrix(X: np.ndarray) -> bool:
    """
    Validate gene expression matrix format.

    Parameters
    ----------
    X : np.ndarray
        Expression matrix.

    Returns
    -------
    valid : bool
        True if matrix is valid, raises ValueError otherwise.
    """
    if len(X.shape) != 2:
        raise ValueError("Expression matrix must be 2D (n_cells, n_genes)")

    if X.shape[0] == 0 or X.shape[1] == 0:
        raise ValueError("Expression matrix cannot have empty dimensions")

    if not np.isfinite(X).all():
        raise ValueError("Expression matrix contains NaN or Inf values")

    return True
