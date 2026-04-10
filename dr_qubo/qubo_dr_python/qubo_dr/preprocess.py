# Author: Selim Romero, Texas A&M University
"""
Preprocessing module for scRNA-seq data.

Handles library-size normalization, log transformation, subsetting by condition,
and computation of gene similarity matrices (Gram matrices via cosine similarity).
"""

import numpy as np
from sklearn.preprocessing import normalize


def normalize_libsize(X, target_sum=1e4):
    """
    Library-size normalize and log-transform scRNA-seq data.

    Scales each cell to a target library size, then applies log1p transformation.

    Parameters
    ----------
    X : np.ndarray, shape (genes, cells)
        Raw count matrix (genes x cells). Typically UMI or read counts.
    target_sum : float, optional
        Target sum per cell (default 1e4). Each cell is scaled to this total,
        then log1p is applied.

    Returns
    -------
    X_norm : np.ndarray, shape (genes, cells)
        Normalized and log-transformed count matrix, same shape as input.

    Notes
    -----
    This is equivalent to Seurat or scanpy's default normalization:
    X_norm = log1p(X / X.sum(axis=0) * target_sum)
    """
    if X.shape[0] == 0 or X.shape[1] == 0:
        return X.copy()

    # Library size normalization per cell
    lib_size = X.sum(axis=0)
    lib_size = np.where(lib_size == 0, 1, lib_size)  # Avoid division by zero
    X_norm = X / lib_size * target_sum

    # Log transformation
    X_norm = np.log1p(X_norm)

    return X_norm


def subset_by_condition(X, cell_ids, condition, label):
    """
    Extract submatrix for cells matching a condition label.

    Parameters
    ----------
    X : np.ndarray, shape (genes, cells)
        Count matrix (genes x cells).
    cell_ids : array-like, shape (cells,)
        Cell identifier array, same length as X.shape[1].
    condition : array-like, shape (cells,)
        Condition label per cell (e.g., array of strings like ['KO', 'WT', 'KO', ...]).
    label : str
        Condition label to match (e.g., 'KO' or 'WT'). Uses string exact match.

    Returns
    -------
    X_sub : np.ndarray, shape (genes, num_matching_cells)
        Submatrix containing only cells where condition == label.
    idx : np.ndarray, dtype bool, shape (cells,)
        Boolean index array indicating which cells were selected.

    Raises
    ------
    ValueError
        If no cells match the given label.
    """
    condition = np.asarray(condition, dtype=str)
    idx = condition == label

    if idx.sum() == 0:
        raise ValueError(f"No cells found matching condition label '{label}'")

    X_sub = X[:, idx]

    return X_sub, idx


def compute_gene_similarity(X):
    """
    Compute row-wise L2 normalization and cosine similarity (Gram matrix).

    Normalizes each gene's expression vector to unit L2 norm, then computes
    the Gram matrix (dot product of normalized genes), which yields cosine similarity.

    Parameters
    ----------
    X : np.ndarray, shape (genes, cells)
        Gene expression matrix (genes x cells).

    Returns
    -------
    Xnorm2 : np.ndarray, shape (genes, genes)
        Gram matrix (cosine similarity). Values in [0, 1].
    Xnorm : np.ndarray, shape (genes, cells)
        Row-normalized gene expression matrix (L2 norm per gene = 1).

    Notes
    -----
    Xnorm2[i, j] = cos(angle between gene i and gene j across cells).
    Diagonal elements are typically 1.0 if genes are normalized.
    """
    # Row-wise L2 normalization
    Xnorm = normalize(X, norm='l2', axis=1)

    # Gram matrix: cosine similarity
    Xnorm2 = Xnorm @ Xnorm.T

    return Xnorm2, Xnorm


import numpy as np
import warnings

def compute_differential(Xref2norm, Xtest2norm, Xref_norm=None, Xtest_norm=None,
                         cs_ref=None, cs_test=None):
    """
    Compute differential co-expression matrix and optional cell state vector.

    dS = Xref2norm - Xtest2norm, with diagonal set to zero.
    Optionally computes Vdiff from cell state projections if provided.

    Parameters
    ----------
    Xref2norm : np.ndarray, shape (genes, genes)
        Reference (WT) Gram matrix (cosine similarity).
    Xtest2norm : np.ndarray, shape (genes, genes)
        Test (KO) Gram matrix (cosine similarity).
    Xref_norm : np.ndarray, shape (genes, cells_ref), optional
        Reference normalized expression matrix.
    Xtest_norm : np.ndarray, shape (genes, cells_test), optional
        Test normalized expression matrix.
    cs_ref : np.ndarray, shape (cells_ref,), optional
        Reference cell state vector.
    cs_test : np.ndarray, shape (cells_test,), optional
        Test cell state vector.

    Returns
    -------
    dS : np.ndarray, shape (genes, genes)
        Differential co-expression matrix.
        Positive values: Reference co-expression gain.
        Negative values: Test co-expression gain.
    Vdiff : np.ndarray, shape (genes,)
        Differential cell state projection vector.
    """
    # dS represents Reference - Test
    dS = Xref2norm - Xtest2norm

    # ── Diagonal sanity check ─────────────────────────────────────────────
    diag_abs  = np.abs(np.diag(dS))
    mean_diag = diag_abs.mean()
    n_nonzero = int((diag_abs > 0.1).sum())

    if n_nonzero > 0:
        print(f"  Diagonal check: {n_nonzero}/{len(diag_abs)} genes have "
              f"|dS diagonal| > 0.1  (mean |diag| = {mean_diag:.3f}). "
              f"Likely genes with near-zero expression in one condition.")

    if mean_diag > 0.15:
        warnings.warn(
            f"compute_differential: mean |diag(dS)| = {mean_diag:.3f}. "
            f"This may indicate systematic normalization issues between conditions.",
            RuntimeWarning, stacklevel=2,
        )

    # Compute Vdiff from cell states if provided
    n_genes = Xref2norm.shape[0]
    Vdiff = np.zeros(n_genes)

    if cs_ref is not None and cs_test is not None:
        # Convert to arrays and normalize
        cs_ref = np.asarray(cs_ref)
        cs_test = np.asarray(cs_test)

        norm_ref = np.linalg.norm(cs_ref)
        norm_test = np.linalg.norm(cs_test)

        if norm_ref > 0:
            cs_ref = cs_ref / norm_ref
        if norm_test > 0:
            cs_test = cs_test / norm_test

        # Project onto gene space
        if Xref_norm is not None and Xtest_norm is not None:
            proj_ref = Xref_norm @ cs_ref
            proj_test = Xtest_norm @ cs_test
            Vdiff = proj_ref - proj_test

    return dS, Vdiff