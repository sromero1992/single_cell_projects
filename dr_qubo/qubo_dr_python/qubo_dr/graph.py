# Author: Selim Romero, Texas A&M University
"""
Graph construction module for gene adjacency matrices.

Builds sparse K-nearest neighbor (KNN) and Mutual Nearest Neighbor (MNN)
adjacency matrices for genes after dimensionality reduction via SVD.
"""

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors


def build_mnn_adjacency(X, method='mnn', K=15, n_svd=50):
    """
    Build sparse KNN or Mutual Nearest Neighbor adjacency matrix for genes.

    Applies SVD dimensionality reduction to the gene space (if n_genes > n_svd),
    then constructs a K-nearest neighbor graph. For MNN, only retains edges where
    both genes are in each other's K-nearest neighbors.

    Parameters
    ----------
    X : np.ndarray, shape (genes, cells)
        Gene expression matrix (genes x cells).
    method : str, {'knn', 'mnn'}, optional
        'knn': standard K-nearest neighbors.
        'mnn': Mutual Nearest Neighbors (default).
    K : int, optional
        Number of neighbors to use (default 15).
    n_svd : int, optional
        SVD components to retain. If n_genes <= n_svd, no reduction applied (default 50).

    Returns
    -------
    adj : scipy.sparse.csr_matrix, shape (genes, genes)
        Sparse symmetric adjacency matrix. adj[i,j] = 1 if j is in i's K-neighbors
        (and vice versa if MNN).

    Notes
    -----
    SVD is applied to reduce the gene feature space (cells) to n_svd components.
    Then KNN is computed in this reduced space using Euclidean distance.
    The resulting adjacency is symmetrized.
    """
    n_genes, n_cells = X.shape

    # Apply SVD reduction if needed
    if n_cells > n_svd:
        svd = TruncatedSVD(n_components=n_svd, random_state=42)
        X_reduced = svd.fit_transform(X.T).T  # Transpose to get gene features
    else:
        X_reduced = X.copy()

    # Build KNN graph
    nbrs = NearestNeighbors(n_neighbors=K + 1, metric='euclidean', n_jobs=-1)
    nbrs.fit(X_reduced.T)  # Fit on cells; find neighbors in gene space
    distances, indices = nbrs.kneighbors(X_reduced.T)

    # Build adjacency matrix from neighbor lists
    row_indices = []
    col_indices = []

    for i in range(n_genes):
        # Add edges from gene i to its K neighbors
        neighbors = indices[i, 1:]  # Skip self (first neighbor)
        for neighbor in neighbors:
            row_indices.append(i)
            col_indices.append(neighbor)

    # Create sparse matrix
    data = np.ones(len(row_indices))
    adj = csr_matrix((data, (row_indices, col_indices)), shape=(n_genes, n_genes))

    # For MNN, intersect with reverse adjacency
    if method == 'mnn':
        adj_t = adj.T
        adj = adj.multiply(adj_t)  # Element-wise multiplication (symmetric)

    # Symmetrize (union of both directions)
    adj = adj + adj.T
    adj.data[adj.data > 0] = 1  # Binarize (remove duplicate weights)

    return adj.tocsr()
