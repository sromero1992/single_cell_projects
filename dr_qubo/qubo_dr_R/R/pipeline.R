# Author: Selim Romero, Texas A&M University
# End-to-end DR-QUBO pipeline for differential co-expression hub discovery.

#' End-to-end QUBO differential co-expression pipeline
#'
#' Identifies the K-gene subnetwork within a pathway whose pairwise
#' co-expression relationships have changed most between a test condition
#' (e.g., knockout, disease) and a reference condition (e.g., wild-type,
#' healthy).
#'
#' The method formulates hub selection as a Quadratic Unconstrained Binary
#' Optimization (QUBO) problem and solves it via simulated annealing.
#'
#' @param X Numeric matrix (genes x cells). Library-size normalized + log1p
#'   scRNA-seq count matrix. Use \code{\link{normalize_libsize}} upstream.
#' @param g Character vector (genes). Gene names, length = nrow(X).
#' @param batch_id Character vector (cells). Condition label per cell.
#'   Must contain entries matching \code{cond_a_label} and \code{cond_b_label}.
#' @param pathway_genelist Character vector or NULL. Explicit pathway gene list.
#'   If NULL, the gene set is fetched from Enrichr using \code{pathway_name}.
#' @param K Integer. Target subnetwork size (default 30).
#' @param pathway_name Character or NULL. Pathway keyword for Enrichr lookup.
#'   Required when \code{pathway_genelist} is NULL.
#' @param library Character or NULL. Enrichr library to restrict to. See
#'   \code{\link{fetch_pathway_genes}}.
#' @param organism Character. \code{"human"} or \code{"mouse"} (default
#'   \code{"human"}).
#' @param cond_a_label Character. Label for condition A (test). Default
#'   \code{"A"}.
#' @param cond_b_label Character. Label for condition B (reference). Default
#'   \code{"B"}.
#' @param test_label Character. Display name for the test condition (default
#'   \code{"disease"}).
#' @param reference_label Character. Display name for the reference condition
#'   (default \code{"control"}).
#' @param method Character. Graph construction: \code{"mnn"} (default) or
#'   \code{"knn"}.
#' @param n_neighbors Integer. Number of neighbors for KNN/MNN (default 30).
#' @param n_svd Integer. SVD components for gene embedding before MNN
#'   (default 50).
#' @param cell_state Numeric vector (N cells) or NULL. Per-cell biological
#'   state scalar (e.g., UCell pathway activity, pseudotime, potency score).
#'   The pipeline splits this by condition automatically.
#'   Pass NULL (default) to omit the V_diff linear bias term.
#' @param extra_genes Character vector or NULL. Additional gene names to
#'   include alongside pathway genes (e.g., GWAS hits, drug targets).
#' @param use_pathway_prior Logical or NULL. Include pathway membership prior
#'   matrix Xnet. NULL (default) = auto: TRUE when extra_genes provided,
#'   FALSE otherwise.
#' @param penalty_scale Numeric. QUBO cardinality penalty multiplier (default 2).
#' @param num_reads Integer. Simulated annealing restarts (default 1000).
#' @param plotit Logical. Generate diagnostic plots (default TRUE).
#' @param top_pct Numeric. Percentage of gene pairs shown as edges in network
#'   plots (default 10).
#' @param run_perm_test Logical. Run hub-level permutation test (default TRUE).
#' @param n_perm Integer. Number of permutations for hub test (default 200).
#'
#' @return A named list with:
#'   \item{sub_g_net}{Character vector (K). Selected hub gene names.}
#'   \item{sub_Q0}{Numeric matrix (K x K). dS submatrix for hub genes.}
#'   \item{sub_Q_net}{Numeric matrix (K x K). = -sub_Q0 (test − reference,
#'     positive = test gain, for heatmap display).}
#'   \item{sub_Qv}{Numeric vector (K). Per-gene contribution scores.}
#'   \item{selected_idx}{Integer vector. Indices (1-based) of selected genes
#'     in the pathway gene set.}
#'   \item{gene_names_pathway}{Character vector (G). All pathway gene names.}
#'   \item{dS}{Numeric matrix (G x G). Full differential co-expression matrix.}
#'   \item{Q0}{Numeric matrix (G x G). = dS (alias, before QUBO assembly).}
#'   \item{Q1}{Numeric matrix (G x G). QUBO matrix before cardinality penalty.}
#'   \item{Q}{Numeric matrix (G x G). Full QUBO matrix with penalty.}
#'   \item{best_energy}{Numeric. Best QUBO objective value.}
#'   \item{hub_perm_result}{List or NULL. Permutation test result from
#'     \code{\link{test_differential_hub}}.}
#'   \item{test_label}{Character. Test condition display label.}
#'   \item{reference_label}{Character. Reference condition display label.}
#'
#' @seealso \code{\link{plot_hub_network}}, \code{\link{test_differential_hub}},
#'   \code{\link{fetch_pathway_genes}}
#'
#' @export
run_pipeline <- function(X, g, batch_id,
                         pathway_genelist  = NULL,
                         K                 = 30,
                         pathway_name      = NULL,
                         library           = NULL,
                         organism          = "human",
                         cond_a_label      = "A",
                         cond_b_label      = "B",
                         test_label        = "disease",
                         reference_label   = "control",
                         method            = "mnn",
                         n_neighbors       = 30,
                         n_svd             = 50,
                         cell_state        = NULL,
                         extra_genes       = NULL,
                         use_pathway_prior = NULL,
                         penalty_scale     = 2,
                         num_reads         = 1000,
                         plotit            = TRUE,
                         top_pct           = 10,
                         run_perm_test     = TRUE,
                         n_perm            = 200) {

  # ─────────────────────────────────────────────────────────────────────────
  # Step 1 — Partition cells by condition
  # ─────────────────────────────────────────────────────────────────────────
  batch_id <- as.character(batch_id)
  idx_a    <- which(batch_id == cond_a_label)
  idx_b    <- which(batch_id == cond_b_label)

  if (length(idx_a) == 0) stop(sprintf("No cells found with label '%s'.", cond_a_label))
  if (length(idx_b) == 0) stop(sprintf("No cells found with label '%s'.", cond_b_label))

  message(sprintf("Condition A (%s): %d cells | Condition B (%s): %d cells",
                  cond_a_label, length(idx_a), cond_b_label, length(idx_b)))

  X_A <- X[, idx_a, drop = FALSE]
  X_B <- X[, idx_b, drop = FALSE]

  # Split cell_state if provided
  cs_A <- if (!is.null(cell_state)) cell_state[idx_a] else NULL
  cs_B <- if (!is.null(cell_state)) cell_state[idx_b] else NULL

  # ─────────────────────────────────────────────────────────────────────────
  # Step 2 — Pathway gene subsetting
  # ─────────────────────────────────────────────────────────────────────────
  if (is.null(pathway_genelist)) {
    if (is.null(pathway_name)) {
      stop("Provide either 'pathway_genelist' or 'pathway_name'.")
    }
    message(sprintf("Fetching pathway genes for '%s' ...", pathway_name))
    pathway_genelist <- fetch_pathway_genes(pathway_name,
                                            organism = organism,
                                            library  = library)
    if (length(pathway_genelist) == 0) {
      stop(sprintf("Pathway '%s' not found. Try a different keyword.", pathway_name))
    }
  }

  # Add extra_genes if provided
  if (!is.null(extra_genes)) {
    pathway_genelist <- unique(c(pathway_genelist, extra_genes))
  }

  # Intersect with genes in data
  g_path <- intersect(pathway_genelist, g)
  if (length(g_path) < K) {
    stop(sprintf("Only %d pathway genes found in data (need >= K=%d).",
                 length(g_path), K))
  }
  message(sprintf("Pathway genes in data: %d / %d", length(g_path), length(pathway_genelist)))

  # Subset X to pathway genes
  path_idx <- match(g_path, g)
  X_A_path <- X_A[path_idx, , drop = FALSE]
  X_B_path <- X_B[path_idx, , drop = FALSE]
  G <- length(g_path)

  # ─────────────────────────────────────────────────────────────────────────
  # Step 3 — Gram matrices and differential co-expression (dS)
  # ─────────────────────────────────────────────────────────────────────────
  message("Computing gene similarity matrices...")
  sim_A <- compute_gene_similarity(X_A_path)
  sim_B <- compute_gene_similarity(X_B_path)

  message("Computing differential co-expression (dS = S_ref - S_test)...")
  diff_res <- compute_differential(sim_B$S, sim_A$S,
                                   X_ref_norm = sim_B$X_norm,
                                   X_test_norm = sim_A$X_norm,
                                   cs_ref  = cs_B,
                                   cs_test = cs_A)
  dS    <- diff_res$dS
  Vdiff <- diff_res$Vdiff

  # ─────────────────────────────────────────────────────────────────────────
  # Step 4 — MNN adjacency matrices
  # ─────────────────────────────────────────────────────────────────────────
  message(sprintf("Building %s adjacency (K=%d, n_svd=%d)...", toupper(method), n_neighbors, n_svd))
  A_mnn_test <- build_mnn_adjacency(X_A_path, method = method,
                                    K = n_neighbors, n_svd = n_svd)
  A_mnn_ref  <- build_mnn_adjacency(X_B_path, method = method,
                                    K = n_neighbors, n_svd = n_svd)

  # ─────────────────────────────────────────────────────────────────────────
  # Step 5 — Pathway prior (Xnet)
  # ─────────────────────────────────────────────────────────────────────────
  if (is.null(use_pathway_prior)) {
    use_pathway_prior <- !is.null(extra_genes)
  }

  Xnet_target <- NULL
  if (use_pathway_prior) {
    message("Including pathway membership prior (Xnet)...")
    core_genes   <- intersect(pathway_genelist, g_path)
    is_core      <- g_path %in% core_genes
    # Off-diagonal entry = -1 if both genes are core pathway members
    Xnet_target  <- outer(is_core * 1, is_core * 1)
    diag(Xnet_target) <- 0
    Xnet_target  <- -Xnet_target
  }

  # ─────────────────────────────────────────────────────────────────────────
  # Step 6 — Assemble and solve QUBO
  # ─────────────────────────────────────────────────────────────────────────
  message("Assembling QUBO matrix...")
  qubo_res <- assemble_qubo_matrix(dS, Vdiff, A_mnn_test, A_mnn_ref, K,
                                   Xnet_target   = Xnet_target,
                                   penalty_scale = penalty_scale)
  Q  <- qubo_res$Q
  Q1 <- qubo_res$Q1

  message(sprintf("Solving QUBO (SA, %d reads)...", num_reads))
  sol <- solve_qubo_sa(Q, K = K, num_reads = num_reads)
  message(sprintf("Best energy: %.4f | K selected: %d", sol$energy, length(sol$selected_idx)))

  selected_idx <- sol$selected_idx
  sub_res      <- extract_subnetwork(dS, Q1, selected_idx)

  sub_g_net <- g_path[selected_idx]
  message("Selected genes: ", paste(sub_g_net, collapse = ", "))

  # ─────────────────────────────────────────────────────────────────────────
  # Step 7 — Hub permutation test
  # ─────────────────────────────────────────────────────────────────────────
  hub_perm_result <- NULL
  if (run_perm_test) {
    message(sprintf("Running hub permutation test (%d permutations)...", n_perm))
    hub_perm_result <- test_differential_hub(
      X_test       = X_A_path,
      X_ref        = X_B_path,
      dS           = dS,
      selected_idx = selected_idx,
      n_perm       = n_perm,
      verbose      = TRUE
    )
  }

  # ─────────────────────────────────────────────────────────────────────────
  # Step 8 — Plots
  # ─────────────────────────────────────────────────────────────────────────
  if (plotit) {
    message("Generating plots...")

    # 3-panel co-expression heatmaps
    p_heatmaps <- plot_condition_heatmaps(
      sim_A$S, sim_B$S, dS, g_path,
      test_label = test_label, reference_label = reference_label
    )
    print(p_heatmaps)

    # Hub network
    tmp_results <- list(
      sub_Q0          = sub_res$sub_Q0,
      sub_g_net       = sub_g_net,
      hub_perm_result = hub_perm_result,
      test_label      = test_label,
      reference_label = reference_label
    )
    p_hub <- plot_hub_network(tmp_results, top_pct = top_pct,
                              test_label = test_label,
                              reference_label = reference_label)
    print(p_hub)
  }

  message("Pipeline complete.")

  return(list(
    sub_g_net        = sub_g_net,
    sub_Q0           = sub_res$sub_Q0,
    sub_Q_net        = sub_res$sub_Q_net,
    sub_Q1           = sub_res$sub_Q1,
    sub_Qv           = sub_res$sub_Qv,
    selected_idx     = selected_idx,
    gene_names_pathway = g_path,
    dS               = dS,
    Q0               = dS,          # alias: Q0 = dS (pure differential signal)
    Q1               = Q1,
    Q                = Q,
    best_energy      = sol$energy,
    hub_perm_result  = hub_perm_result,
    test_label       = test_label,
    reference_label  = reference_label
  ))
}
