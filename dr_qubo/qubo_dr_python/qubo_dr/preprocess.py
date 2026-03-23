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


def compute_differential(Xwt2norm, Xko2norm, Xwt_norm=None, Xko_norm=None,
                        cs_wt=None, cs_ko=None):
    """
    Compute differential co-expression matrix and optional cell state vector.

    Xdiff = Xwt2norm - Xko2norm, with diagonal set to zero.
    Optionally computes Vdiff from cell state projections if provided.

    Parameters
    ----------
    Xwt2norm : np.ndarray, shape (genes, genes)
        WT Gram matrix (cosine similarity).
    Xko2norm : np.ndarray, shape (genes, genes)
        KO Gram matrix (cosine similarity).
    Xwt_norm : np.ndarray, shape (genes, cells_wt), optional
        WT normalized expression matrix (needed if cs_wt is provided).
    Xko_norm : np.ndarray, shape (genes, cells_ko), optional
        KO normalized expression matrix (needed if cs_ko is provided).
    cs_wt : np.ndarray, shape (cells_wt,), optional
        WT cell state vector (e.g., from trajectory or differentiation stage).
        If provided, Vdiff is computed from projection onto gene space.
    cs_ko : np.ndarray, shape (cells_ko,), optional
        KO cell state vector.

    Returns
    -------
    Xdiff : np.ndarray, shape (genes, genes)
        Differential co-expression matrix. Diagonal set to 0.
        Positive values: WT co-expression gain.
        Negative values: KO co-expression gain.
    Vdiff : np.ndarray, shape (genes,)
        Differential cell state projection vector. All zeros if cell states not provided.

    Notes
    -----
    If cell states are provided:
    - They are L2-normalized per condition.
    - Projected onto normalized gene space: proj_wt = Xwt_norm @ cs_wt / ||cs_wt||
    - Vdiff = proj_wt - proj_ko
    """
    Xdiff = Xwt2norm - Xko2norm
    np.fill_diagonal(Xdiff, 0)

    # Compute Vdiff from cell states if provided
    n_genes = Xwt2norm.shape[0]
    Vdiff = np.zeros(n_genes)

    if cs_wt is not None and cs_ko is not None:
        cs_wt = np.asarray(cs_wt)
        cs_ko = np.asarray(cs_ko)

        # Normalize cell state vectors
        norm_wt = np.linalg.norm(cs_wt)
        norm_ko = np.linalg.norm(cs_ko)

        if norm_wt > 0:
            cs_wt = cs_wt / norm_wt
        if norm_ko > 0:
            cs_ko = cs_ko / norm_ko

        # Project onto gene space
        if Xwt_norm is not None and Xko_norm is not None:
            proj_wt = Xwt_norm @ cs_wt
            proj_ko = Xko_norm @ cs_ko
            Vdiff = proj_wt - proj_ko

    return Xdiff, Vdiff
