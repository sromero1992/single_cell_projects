# Author: Selim Romero, Texas A&M University
# QUBO matrix assembly, simulated-annealing solver, and subnetwork extraction.

#' Assemble the QUBO matrix with cardinality penalty
#'
#' Assembles the base QUBO matrix
#' \deqn{Q_1 = \Delta S + \mathrm{diag}(V_\mathrm{diff}) - (A_{\mathrm{MNN,test}} + A_{\mathrm{MNN,ref}})}
#' then adds the cardinality penalty to enforce exactly K gene selections.
#'
#' @param dS Numeric matrix (G x G). Differential co-expression matrix,
#'   \eqn{\Delta S = S_\mathrm{ref} - S_\mathrm{test}}. Negative off-diagonal
#'   entries indicate co-expression gained in the test condition.
#' @param Vdiff Numeric vector (G). Gene-wise cell-state difference vector.
#'   Pass \code{rep(0, G)} to omit this term.
#' @param A_mnn_test Sparse or dense matrix (G x G). MNN adjacency for the
#'   test condition.
#' @param A_mnn_ref Sparse or dense matrix (G x G). MNN adjacency for the
#'   reference condition.
#' @param K Integer. Target subnetwork size.
#' @param Xnet_target Numeric matrix (G x G), optional. Pathway membership
#'   prior: off-diagonal entry (i,j) = -1 if both genes are core pathway
#'   members, 0 otherwise. Pass \code{NULL} to skip (default).
#' @param penalty_scale Numeric. Safety multiplier for the cardinality
#'   penalty (default 2).
#' @param auto_penalty Logical. If \code{TRUE} (default), use spectral-norm-
#'   derived penalty \eqn{P = penalty\_scale \times \|Q_1\|_2}. If
#'   \code{FALSE}, use legacy elementwise maximum.
#'
#' @return A list with:
#'   \item{Q}{Final QUBO matrix with cardinality penalty (G x G).}
#'   \item{Q1}{Base QUBO matrix before penalty (G x G).}
#'   \item{P}{Penalty coefficient used.}
#'
#' @export
assemble_qubo_matrix <- function(dS, Vdiff, A_mnn_test, A_mnn_ref, K,
                                 Xnet_target = NULL,
                                 penalty_scale = 2,
                                 auto_penalty  = TRUE) {
  ng <- nrow(dS)

  stopifnot(length(Vdiff)     == ng,
            nrow(A_mnn_test)  == ng,
            nrow(A_mnn_ref)   == ng)

  # Convert sparse to dense if needed
  if (methods::is(A_mnn_test, "sparseMatrix"))
    A_mnn_test <- as.matrix(A_mnn_test)
  if (methods::is(A_mnn_ref, "sparseMatrix"))
    A_mnn_ref <- as.matrix(A_mnn_ref)

  # Base QUBO matrix
  Q1 <- dS + diag(Vdiff) - (A_mnn_test + A_mnn_ref)

  if (!is.null(Xnet_target)) {
    Q1 <- Q1 + Xnet_target
    message("  Pathway prior Xnet: included")
  } else {
    message("  Pathway prior Xnet: skipped (pure pathway-gene analysis)")
  }

  # Cardinality penalty
  if (auto_penalty) {
    sigma_max <- max(abs(eigen(Q1, only.values = TRUE)$values))
    P <- penalty_scale * sigma_max
    message(sprintf("  Penalty (spectral): ||Q1||_2 = %.4f  ->  P = %.4f  (scale = %.1f)",
                    sigma_max, P, penalty_scale))
  } else {
    P <- penalty_scale * max(abs(Q1))
    message(sprintf("  Penalty (legacy max): P = %.4f  (scale = %.1f)", P, penalty_scale))
  }

  # Build full penalty matrix: all entries = P, then fix diagonal
  Q_penalty         <- matrix(P, nrow = ng, ncol = ng)
  diag(Q_penalty)   <- P * (1 - 2 * K)

  Q <- Q1 + Q_penalty
  Q <- (Q + t(Q)) / 2   # symmetrise

  message(sprintf("QUBO assembled: %dx%d | K=%d | P=%.4f", ng, ng, K, P))

  return(list(Q = Q, Q1 = Q1, P = P))
}


#' Solve a QUBO problem via Simulated Annealing
#'
#' Minimises \eqn{E(\mathbf{z}) = \mathbf{z}^T Q \mathbf{z}} over binary vectors
#' \eqn{\mathbf{z} \in \{0,1\}^G} by simulated annealing.  The cardinality constraint
#' \eqn{|\mathbf{z}| = K} is enforced via swap moves (flip one 1→0 and one 0→1
#' simultaneously), so the constraint is never violated.
#'
#' @param Q Numeric matrix (G x G). QUBO matrix.
#' @param K Integer. Exact cardinality constraint: exactly K ones in solution.
#' @param num_reads Integer. Number of independent SA restarts (default 1000).
#'   The restart with lowest energy is returned.
#' @param seed Integer. Random seed (default 42).
#' @param T_start Numeric. Initial temperature (default auto: 0.1 * |Q|_inf).
#' @param T_end Numeric. Final temperature (default 1e-6).
#' @param n_steps Integer. SA steps per restart (default 5000).
#'
#' @return A list with:
#'   \item{z}{Integer vector (G). Best binary solution (0/1).}
#'   \item{energy}{Numeric. Best energy achieved.}
#'   \item{selected_idx}{Integer vector. Indices (1-based) of selected genes.}
#'
#' @export
solve_qubo_sa <- function(Q, K, num_reads = 1000, seed = 42,
                          T_start = NULL, T_end = 1e-6, n_steps = 5000) {
  set.seed(seed)
  ng <- nrow(Q)

  if (is.null(T_start)) T_start <- 0.1 * max(abs(Q))
  if (T_start <= 0) T_start <- 0.1

  best_z      <- NULL
  best_energy <- Inf

  for (read in seq_len(num_reads)) {
    # Random initial solution with exactly K ones
    z <- integer(ng)
    z[sample(ng, K)] <- 1L

    E <- as.numeric(t(z) %*% Q %*% z)

    T_curr <- T_start

    for (step in seq_len(n_steps)) {
      # Exponential cooling schedule
      T_curr <- T_start * (T_end / T_start)^(step / n_steps)

      # Swap move: pick one 1 and one 0, swap them
      ones  <- which(z == 1L)
      zeros <- which(z == 0L)
      if (length(ones) == 0 || length(zeros) == 0) break

      i <- ones[sample(length(ones),  1)]
      j <- zeros[sample(length(zeros), 1)]

      # Compute delta energy for swap (i: 1->0, j: 0->1)
      # dE = E_new - E_old
      # Efficient O(G) update: dE = Q[j,j]-Q[i,i] + 2*(col_j - col_i)^T z
      #   after excluding i and j contributions
      z_copy    <- z
      z_copy[i] <- 0L
      z_copy[j] <- 1L
      E_new <- as.numeric(t(z_copy) %*% Q %*% z_copy)
      dE    <- E_new - E

      if (dE < 0 || runif(1) < exp(-dE / T_curr)) {
        z <- z_copy
        E <- E_new
      }
    }

    if (E < best_energy) {
      best_energy <- E
      best_z      <- z
    }
  }

  selected_idx <- which(best_z == 1L)

  return(list(z = best_z, energy = best_energy, selected_idx = selected_idx))
}


#' Extract the K-gene subnetwork from QUBO solution
#'
#' Given the QUBO solution binary vector, extracts the K×K submatrix of
#' \code{dS} (Q0) and computes per-gene contribution scores.
#'
#' @param dS Numeric matrix (G x G). Differential co-expression matrix.
#' @param Q1 Numeric matrix (G x G). Base QUBO matrix (before penalty).
#' @param selected_idx Integer vector. Indices (1-based) of selected genes.
#'
#' @return A list with:
#'   \item{sub_Q0}{Numeric matrix (K x K). dS submatrix for selected genes.}
#'   \item{sub_Q1}{Numeric matrix (K x K). Q1 submatrix for selected genes.}
#'   \item{sub_Qv}{Numeric vector (K). Per-gene contribution: row sums of
#'     \code{sub_Q0} (total differential co-expression with other hub genes).}
#'   \item{sub_Q_net}{Numeric matrix (K x K). = -sub_Q0.  Positive entries
#'     represent test/disease co-expression gain (biologically natural sign
#'     for plotting: red = test gain).}
#'
#' @export
extract_subnetwork <- function(dS, Q1, selected_idx) {
  sub_Q0    <- dS[selected_idx, selected_idx]
  sub_Q1    <- Q1[selected_idx, selected_idx]
  sub_Qv    <- rowSums(sub_Q0)      # per-gene sum of differential co-expression
  sub_Q_net <- -sub_Q0              # negate: positive = test gain (for plotting)

  return(list(sub_Q0    = sub_Q0,
              sub_Q1    = sub_Q1,
              sub_Qv    = sub_Qv,
              sub_Q_net = sub_Q_net))
}
