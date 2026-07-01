"""
normalize.py — Normalization utilities
=======================================
Library-size normalization + log1p, matching MATLAB pkg.norm_libsize(X, 1e4)
and Scanpy sc.pp.normalize_total() + sc.pp.log1p().
"""

from __future__ import annotations

import warnings
import numpy as np
import scipy.sparse as sp


def normalize_lib_size(
    X: "np.ndarray | sp.spmatrix",
    target_sum: float = 1e4,
) -> np.ndarray:
    """
    Library-size normalize a genes × cells (or cells × genes) matrix.

    Each cell is scaled so that its total count equals ``target_sum``, then
    log1p-transformed.  Cells with zero total counts are left as zeros.

    Parameters
    ----------
    X : array or sparse matrix  (cells × genes)
        Raw count matrix.
    target_sum : float
        Target library size (default 10 000).

    Returns
    -------
    X_norm : np.ndarray  (cells × genes), dense float32
    """
    if sp.issparse(X):
        X = X.toarray()
    X = np.array(X, dtype=np.float32)

    col_sums = X.sum(axis=1, keepdims=True)          # per-cell totals
    col_sums[col_sums == 0] = 1.0                     # avoid div-by-zero
    X_norm = (X / col_sums) * target_sum
    return np.log1p(X_norm)


def normalize_adata(
    adata,
    norm_str: str = "lib_size",
    target_sum: float = 1e4,
    layer: str | None = None,
    inplace: bool = True,
):
    """
    Normalize an AnnData object for TimeSCape.

    Parameters
    ----------
    adata : AnnData
        Input object.  ``adata.X`` must contain raw counts when
        ``norm_str="lib_size"``.
    norm_str : str
        ``"lib_size"``  — normalize to ``target_sum`` + log1p (default).
        ``"logcounts"`` — use existing ``adata.layers["logcounts"]`` or
                         ``adata.X`` if already log-normalized by Scanpy.
        ``"none"``      — use ``adata.X`` as-is.
    target_sum : float
        Library size target (default 10 000).
    layer : str, optional
        If given, read from ``adata.layers[layer]`` instead of ``adata.X``.
    inplace : bool
        If True, store result back into ``adata.X`` and return ``adata``.
        If False, return the normalized dense matrix without modifying adata.

    Returns
    -------
    adata (inplace=True) or np.ndarray (inplace=False)
    """
    import scanpy as sc  # optional import — keeps package lightweight

    X_src = adata.layers[layer] if layer else adata.X

    if norm_str == "lib_size":
        X_norm = normalize_lib_size(X_src, target_sum=target_sum)
        if inplace:
            adata.X = sp.csr_matrix(X_norm)
            return adata
        return X_norm

    elif norm_str == "logcounts":
        # Assume data is already log-normalized
        if sp.issparse(X_src):
            X_norm = X_src.toarray().astype(np.float32)
        else:
            X_norm = np.asarray(X_src, dtype=np.float32)
        if inplace:
            adata.X = sp.csr_matrix(X_norm)
            return adata
        return X_norm

    elif norm_str == "none":
        if sp.issparse(X_src):
            X_norm = X_src.toarray().astype(np.float32)
        else:
            X_norm = np.asarray(X_src, dtype=np.float32)
        if inplace:
            return adata
        return X_norm

    else:
        raise ValueError(
            f"Unknown norm_str '{norm_str}'. "
            "Choose 'lib_size', 'logcounts', or 'none'."
        )
