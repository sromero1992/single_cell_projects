# Author: Selim Romero, Texas A&M University
# Permutation test for hub-level differential co-expression significance.

#' Hub-level permutation test for differential co-expression
#'
#' Tests whether the QUBO-selected K-gene hub subnetwork is significantly more
#' differential than expected by chance. The null distribution is built by
#' randomly permuting condition labels (shuffling which cells are "test" and
#' which are "reference"), recomputing \eqn{\Delta S}, and measuring the hub
#' energy on the permuted matrix.
#'
#' The test statistic is:
#' \deqn{E_\mathrm{real} = (\mathbf{z}^*)^T \Delta S\, \mathbf{z}^* = \sum_{i,j \in \mathrm{hub}} \Delta S(i,j)}
#'
#' where \eqn{\mathbf{z}^*} is the fixed QUBO solution (binary indicator of
#' selected hub genes). Note that \eqn{\mathbf{z}^*} is NOT re-optimized for each
#' permutation — this avoids circularity and measures only the biological signal.
#'
#' @param X_test Numeric matrix (genes x test_cells). Test condition expression.
#' @param X_ref Numeric matrix (genes x ref_cells). Reference condition expression.
#' @param dS Numeric matrix (G x G). Observed differential co-expression matrix.
#' @param selected_idx Integer vector. Indices (1-based) of selected hub genes.
#' @param n_perm Integer. Number of label permutations (default 200).
#' @param seed Integer. Random seed (default 42).
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return A list with:
#'   \item{pval}{Numeric. Permutation p-value (fraction of null energies <= real
#'     energy). Small p-value means the hub is significantly more differential
#'     than random — i.e., the selected genes gained co-expression in the test
#'     condition.}
#'   \item{zscore}{Numeric. Z-score: \code{(E_real - mean(E_null)) / sd(E_null)}.
#'     Negative z-score = hub energy is more negative than null (test gained).}
#'   \item{E_real}{Numeric. Hub energy on the observed data.}
#'   \item{E_null}{Numeric vector (n_perm). Null distribution of hub energies.}
#'   \item{null_mean}{Numeric. Mean of null distribution.}
#'   \item{null_sd}{Numeric. Standard deviation of null distribution.}
#'
#' @details
#'   Sign: because \code{dS = S_ref - S_test}, a negative E_real means the hub
#'   genes collectively have LESS co-expression in the reference than in the test
#'   — i.e., the test/disease condition gained co-expression. A significantly
#'   negative z-score therefore indicates a disease-associated hub.
#'
#' @export
test_differential_hub <- function(X_test, X_ref, dS, selected_idx,
                                  n_perm = 200, seed = 42, verbose = TRUE) {
  set.seed(seed)

  n_test <- ncol(X_test)
  n_ref  <- ncol(X_ref)
  n_tot  <- n_test + n_ref

  # Combined expression matrix
  X_all <- cbind(X_test, X_ref)

  # Binary hub indicator vector
  bz <- rep(0, nrow(dS))
  bz[selected_idx] <- 1

  # Real hub energy
  E_real <- as.numeric(bz %*% dS %*% bz)

  # Permutation null
  E_null <- numeric(n_perm)

  for (p in seq_len(n_perm)) {
    # Shuffle condition labels
    perm    <- sample(n_tot)
    idx_a_p <- perm[seq_len(n_test)]
    idx_b_p <- perm[(n_test + 1):n_tot]

    Xa_p <- X_all[, idx_a_p]
    Xb_p <- X_all[, idx_b_p]

    # L2-normalize rows
    norm_a <- sqrt(rowSums(Xa_p^2)); norm_a[norm_a == 0] <- 1
    norm_b <- sqrt(rowSums(Xb_p^2)); norm_b[norm_b == 0] <- 1
    Xa_p_n <- Xa_p / norm_a
    Xb_p_n <- Xb_p / norm_b

    # Gram matrices
    Sa_p <- Xa_p_n %*% t(Xa_p_n)
    Sb_p <- Xb_p_n %*% t(Xb_p_n)

    dS_p <- Sb_p - Sa_p

    E_null[p] <- as.numeric(bz %*% dS_p %*% bz)

    if (verbose && p %% 50 == 0) {
      message(sprintf("  Permutation %d / %d", p, n_perm))
    }
  }

  null_mean <- mean(E_null)
  null_sd   <- sd(E_null)

  zscore <- if (null_sd > 0) (E_real - null_mean) / null_sd else 0
  pval   <- mean(E_null <= E_real)

  if (verbose) {
    stars <- if      (pval < 0.001) "***"
              else if (pval < 0.01)  "**"
              else if (pval < 0.05)  "*"
              else                   "n.s."
    message(sprintf(
      "Hub permutation test: E_real=%.4f  null_mean=%.4f  z=%.2f  p=%.4f  %s",
      E_real, null_mean, zscore, pval, stars))
  }

  return(list(pval      = pval,
              zscore    = zscore,
              E_real    = E_real,
              E_null    = E_null,
              null_mean = null_mean,
              null_sd   = null_sd))
}
