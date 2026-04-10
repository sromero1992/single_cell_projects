# Author: Selim Romero, Texas A&M University
# Graph construction: KNN / MNN adjacency matrices for gene co-expression networks.

#' Build KNN or MNN adjacency matrix for genes
#'
#' Computes a binary (or weighted) adjacency matrix between genes using
#' K-Nearest Neighbors (KNN) or Mutual Nearest Neighbors (MNN).
#'
#' For MNN: gene i and gene j are connected if i is among the K nearest
#' neighbors of j AND j is among the K nearest neighbors of i — i.e.,
#' the relationship is mutual. This is more conservative and reduces
#' spurious edges from asymmetric similarity.
#'
#' @param X Numeric matrix, shape (genes x cells). Gene expression matrix.
#' @param method Character. Either \code{"knn"} or \code{"mnn"} (default
#'   \code{"mnn"}).
#' @param K Integer. Number of neighbors (default 30).
#' @param n_svd Integer. Number of SVD components for dimensionality reduction
#'   before building the neighbor graph (default 50). Reduces noise and
#'   computational cost for large gene sets.
#'
#' @return A symmetric binary adjacency matrix (genes x genes), class
#'   \code{\link[Matrix]{sparseMatrix}}. Entry (i,j) = 1 if genes i and j
#'   are neighbors (KNN: at least one direction; MNN: both directions).
#'
#' @details
#'   SVD is performed on the cell dimension (columns) so each gene obtains a
#'   \code{n_svd}-dimensional coordinate in cell-space; KNN is then computed
#'   between genes in that reduced space using Euclidean distance.
#'
#' @importFrom Matrix sparseMatrix
#' @export
build_mnn_adjacency <- function(X, method = "mnn", K = 30, n_svd = 50) {
  n_genes <- nrow(X)
  n_cells <- ncol(X)

  K <- min(K, n_genes - 1)

  # ── Dimensionality reduction via truncated SVD ───────────────────────────
  if (n_svd > 0 && min(n_genes, n_cells) > n_svd) {
    # irlba for memory-efficient truncated SVD
    if (requireNamespace("irlba", quietly = TRUE)) {
      sv  <- irlba::irlba(X, nv = n_svd)
      emb <- sv$u %*% diag(sv$d)    # genes × n_svd embedding
    } else {
      # Fall back to base SVD (slower)
      sv  <- svd(X, nu = n_svd, nv = 0)
      emb <- sv$u %*% diag(sv$d[seq_len(n_svd)])
    }
  } else {
    emb <- X  # use raw expression
  }

  # ── Pairwise Euclidean distances between genes in embedding space ────────
  dist_mat <- as.matrix(dist(emb, method = "euclidean"))
  diag(dist_mat) <- Inf   # exclude self

  # ── Build KNN adjacency ───────────────────────────────────────────────────
  knn_adj <- matrix(0L, nrow = n_genes, ncol = n_genes)

  for (i in seq_len(n_genes)) {
    nn_idx <- order(dist_mat[i, ])[seq_len(K)]
    knn_adj[i, nn_idx] <- 1L
  }

  if (method == "mnn") {
    # Mutual: both i→j and j→i required
    adj <- knn_adj & t(knn_adj)
  } else {
    # KNN: either direction suffices (union)
    adj <- (knn_adj | t(knn_adj)) * 1L
  }

  # Symmetrize and return as sparse matrix
  adj <- (adj + t(adj)) > 0
  adj <- Matrix::sparseMatrix(
    i = which(adj, arr.ind = TRUE)[, 1],
    j = which(adj, arr.ind = TRUE)[, 2],
    x = 1,
    dims = c(n_genes, n_genes)
  )

  return(adj)
}
