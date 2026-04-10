# Author: Selim Romero, Texas A&M University
# Preprocessing functions for scRNA-seq data.

#' Library-size normalize and log-transform scRNA-seq data
#'
#' Scales each cell to a target library size, then applies log1p transformation.
#' Equivalent to Seurat / scanpy default normalization.
#'
#' @param X Numeric matrix, shape (genes x cells). Raw count matrix.
#' @param target_sum Numeric. Target sum per cell (default 1e4).
#'
#' @return Numeric matrix, same shape as X, library-size normalized and log1p
#'   transformed.
#'
#' @export
normalize_libsize <- function(X, target_sum = 1e4) {
  if (ncol(X) == 0 || nrow(X) == 0) return(X)

  # Library size per cell (column sums)
  lib_size <- colSums(X)
  lib_size[lib_size == 0] <- 1  # avoid division by zero

  # Scale and log-transform
  X_norm <- sweep(X, 2, lib_size, FUN = "/") * target_sum
  X_norm <- log1p(X_norm)

  return(X_norm)
}


#' Compute row-wise L2 normalization and cosine similarity (Gram matrix)
#'
#' Normalizes each gene's expression vector to unit L2 norm, then computes
#' the Gram matrix (dot product of normalized genes), yielding cosine similarity.
#'
#' @param X Numeric matrix, shape (genes x cells).
#'
#' @return A list with:
#'   \item{S}{Numeric matrix, shape (genes x genes). Gram matrix (cosine
#'     similarity), values in [0, 1].}
#'   \item{X_norm}{Numeric matrix, shape (genes x cells). Row-normalized
#'     gene expression (L2 norm per gene = 1).}
#'
#' @export
compute_gene_similarity <- function(X) {
  # Row-wise L2 norms
  row_norms <- sqrt(rowSums(X^2))
  row_norms[row_norms == 0] <- 1  # avoid division by zero

  # Normalize rows to unit L2 norm
  X_norm <- sweep(X, 1, row_norms, FUN = "/")

  # Gram matrix: cosine similarity (genes x genes)
  S <- X_norm %*% t(X_norm)
  S <- pmax(S, 0)  # non-negative cosine (clip negatives to 0)

  return(list(S = S, X_norm = X_norm))
}


#' Compute differential co-expression matrix and optional cell-state vector
#'
#' Computes \code{dS = S_ref - S_test}, the differential Gram matrix.
#' Positive entries indicate reference co-expression gain; negative entries
#' indicate test/disease co-expression gain.
#'
#' Optionally computes \code{Vdiff} from per-cell biological state scalars
#' (e.g., pseudotime, pathway activity score, potency).
#'
#' @param S_ref Numeric matrix (genes x genes). Reference Gram matrix.
#' @param S_test Numeric matrix (genes x genes). Test Gram matrix.
#' @param X_ref_norm Numeric matrix (genes x ref_cells), optional.
#'   Reference normalized expression. Required for Vdiff computation.
#' @param X_test_norm Numeric matrix (genes x test_cells), optional.
#'   Test normalized expression. Required for Vdiff computation.
#' @param cs_ref Numeric vector (ref_cells), optional. Per-cell state scalar
#'   for reference condition.
#' @param cs_test Numeric vector (test_cells), optional. Per-cell state scalar
#'   for test condition.
#'
#' @return A list with:
#'   \item{dS}{Numeric matrix (genes x genes). Differential co-expression.
#'     Positive = reference gain, negative = test/disease gain.}
#'   \item{Vdiff}{Numeric vector (genes). Differential cell-state projection.
#'     Zero vector if cell-state scalars not provided.}
#'
#' @details
#'   Sign convention: \code{dS = S_ref - S_test}. When plotting, negate
#'   \code{dS} so that positive values (red) represent test/disease co-expression
#'   gain — the biologically natural direction.
#'
#' @export
compute_differential <- function(S_ref, S_test,
                                 X_ref_norm = NULL, X_test_norm = NULL,
                                 cs_ref = NULL, cs_test = NULL) {
  dS <- S_ref - S_test

  # Zero the diagonal (by construction S(i,i)=1 for both, so dS(i,i)=0,
  # but guard against numerical drift)
  diag(dS) <- 0

  # Vdiff: differential cell-state projection
  n_genes <- nrow(dS)
  Vdiff <- rep(0, n_genes)

  if (!is.null(cs_ref) && !is.null(cs_test) &&
      !is.null(X_ref_norm) && !is.null(X_test_norm)) {

    cs_ref  <- as.numeric(cs_ref)
    cs_test <- as.numeric(cs_test)

    norm_ref  <- sqrt(sum(cs_ref^2))
    norm_test <- sqrt(sum(cs_test^2))

    if (norm_ref  > 0) cs_ref  <- cs_ref  / norm_ref
    if (norm_test > 0) cs_test <- cs_test / norm_test

    proj_ref  <- X_ref_norm  %*% cs_ref
    proj_test <- X_test_norm %*% cs_test
    Vdiff <- as.numeric(proj_ref - proj_test)
  }

  return(list(dS = dS, Vdiff = Vdiff))
}
