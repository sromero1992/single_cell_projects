# Author: Selim Romero, Texas A&M University
# Visualization functions for DR-QUBO results.

#' Plot co-expression heatmap
#'
#' Plots a symmetric co-expression (or differential co-expression) matrix as
#' a heatmap with a divergent colormap centered at zero.
#'
#' @param mat Numeric matrix (n x n). Co-expression matrix to visualize.
#' @param gene_names Character vector (n). Gene names for axis labels.
#' @param title Character. Plot title.
#' @param cmap Character. Color scale to use: \code{"diverging"} (default,
#'   red-white-blue) or \code{"sequential"} (Blues).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 scale_fill_distiller
#'   theme element_text element_blank labs coord_fixed
#' @export
plot_coexpr_heatmap <- function(mat, gene_names,
                                title = "Differential Co-expression",
                                cmap  = "diverging") {
  n <- nrow(mat)
  df <- expand.grid(x = seq_len(n), y = seq_len(n))
  df$value <- as.vector(mat)
  df$x_name <- factor(gene_names[df$x], levels = gene_names)
  df$y_name <- factor(gene_names[df$y], levels = rev(gene_names))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_name, y = y_name, fill = value)) +
    ggplot2::geom_tile()

  if (cmap == "diverging") {
    max_abs <- max(abs(df$value), na.rm = TRUE)
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs, max_abs),
      name = "Value"
    )
  } else {
    p <- p + ggplot2::scale_fill_distiller(palette = "Blues", direction = 1,
                                           limits = c(0, 1), name = "Similarity")
  }

  p <- p +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
      axis.text.y     = ggplot2::element_text(size = 7),
      plot.title      = ggplot2::element_text(face = "bold", size = 11),
      panel.background = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank()
    )

  return(p)
}


#' Plot 3-panel condition comparison heatmaps
#'
#' Creates side-by-side heatmaps: test co-expression, reference co-expression,
#' and the differential (test − reference). Positive values in the differential
#' panel (red) indicate test/disease co-expression gain.
#'
#' @param S_test Numeric matrix (G x G). Test condition Gram matrix.
#' @param S_ref Numeric matrix (G x G). Reference condition Gram matrix.
#' @param dS Numeric matrix (G x G). Differential co-expression (reference − test).
#'   The function internally negates this to display test − reference.
#' @param gene_names Character vector (G). Gene names.
#' @param test_label Character. Label for the test condition (default \code{"test"}).
#' @param reference_label Character. Label for the reference condition
#'   (default \code{"reference"}).
#'
#' @return A \code{ggplot} object (3-panel via \code{patchwork} if available,
#'   otherwise prints three separate plots).
#'
#' @export
plot_condition_heatmaps <- function(S_test, S_ref, dS, gene_names,
                                    test_label = "test",
                                    reference_label = "reference") {
  p1 <- plot_coexpr_heatmap(S_test, gene_names,
                             title = paste(test_label, "Co-expression"),
                             cmap  = "sequential")
  p2 <- plot_coexpr_heatmap(S_ref, gene_names,
                             title = paste(reference_label, "Co-expression"),
                             cmap  = "sequential")

  # Negate dS so positive = test gain (red)
  diff_display <- -dS
  p3 <- plot_coexpr_heatmap(diff_display, gene_names,
                             title = paste0("Differential (", test_label,
                                            " \u2212 ", reference_label, ")"),
                             cmap  = "diverging")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p1 + p2 + p3)
  } else {
    message("Install 'patchwork' for combined 3-panel plot. Printing separately.")
    print(p1); print(p2); print(p3)
    invisible(list(test = p1, ref = p2, diff = p3))
  }
}


#' Plot hub gene co-expression network
#'
#' Visualizes the QUBO-selected hub subnetwork as a force-directed graph.
#' Edge color encodes direction of co-expression change:
#' \itemize{
#'   \item \strong{Red}: test/disease co-expression gain (mat > 0, i.e. -dS > 0)
#'   \item \strong{Blue}: reference co-expression gain (mat < 0)
#' }
#'
#' @param results List. Output of \code{\link{run_pipeline}}. Must contain
#'   \code{sub_Q0} (K×K dS submatrix) and \code{sub_g_net} (K gene names).
#' @param top_pct Numeric or NULL. Keep top X% of hub edges by |dS| strength.
#'   If \code{NULL} and \code{edge_threshold} is also \code{NULL}, defaults to
#'   top 25%.
#' @param edge_threshold Numeric or NULL. Absolute |dS| cutoff for edges.
#' @param n_sigma Numeric or NULL. Display edges with |dS(i,j)| >= n_sigma *
#'   sigma_edge (statistically principled; requires permutation test results).
#' @param layout Character. Layout algorithm: \code{"fr"} (Fruchterman-Reingold,
#'   default), \code{"kk"} (Kamada-Kawai), or \code{"circle"}.
#' @param title Character or NULL. Plot title. Auto-generated if NULL.
#' @param test_label Character. Label for test condition.
#' @param reference_label Character. Label for reference condition.
#'
#' @return A \code{ggplot} object (using \pkg{ggraph} if available) or base R
#'   igraph plot.
#'
#' @importFrom igraph graph_from_data_frame V E
#' @export
plot_hub_network <- function(results,
                             top_pct         = NULL,
                             edge_threshold  = NULL,
                             n_sigma         = NULL,
                             layout          = "fr",
                             title           = NULL,
                             test_label      = NULL,
                             reference_label = NULL) {
  sub_Q0    <- results[["sub_Q0"]]
  sub_g_net <- results[["sub_g_net"]]

  if (is.null(sub_Q0) || is.null(sub_g_net)) {
    stop("'results' must contain 'sub_Q0' and 'sub_g_net' (output of run_pipeline).")
  }

  perm_result <- results[["hub_perm_result"]]
  perm_pval   <- perm_result[["pval"]]
  perm_zscore <- perm_result[["zscore"]]
  E_null      <- perm_result[["E_null"]]

  if (is.null(test_label))      test_label      <- results[["test_label"]]      %||% "test"
  if (is.null(reference_label)) reference_label <- results[["reference_label"]] %||% "reference"

  K       <- length(sub_g_net)
  n_pairs <- K * (K - 1) / 2

  # mat: positive = test/disease co-expression gain
  mat <- -sub_Q0

  # Upper-triangle values
  triu_idx  <- which(upper.tri(mat), arr.ind = TRUE)
  triu_vals <- abs(mat[triu_idx])
  max_w     <- max(triu_vals, na.rm = TRUE)
  if (max_w == 0) max_w <- 1

  # ── Edge threshold ─────────────────────────────────────────────────────
  if (!is.null(n_sigma)) {
    if (is.null(E_null)) {
      stop("n_sigma requires permutation test results. Run run_pipeline with run_perm_test=TRUE.")
    }
    sigma_null <- sd(E_null)
    sigma_edge <- sigma_null / (2 * sqrt(max(n_pairs, 1)))
    display_threshold <- n_sigma * sigma_edge
    mode_label <- sprintf("|dS| >= %.2f * sigma_edge (sigma_edge=%.4f)", n_sigma, sigma_edge)
  } else if (!is.null(edge_threshold)) {
    display_threshold <- edge_threshold
    mode_label <- sprintf("|dS| >= %.4f", edge_threshold)
  } else {
    pct <- if (!is.null(top_pct)) top_pct else 25
    display_threshold <- if (length(triu_vals) > 0 && pct < 100)
      quantile(triu_vals, 1 - pct / 100) else 0
    mode_label <- sprintf("top %.4g%% edges", pct)
  }

  # ── Build edge data frame ──────────────────────────────────────────────
  edges_df <- data.frame(
    from   = character(0), to = character(0),
    weight = numeric(0),   color = character(0),
    stringsAsFactors = FALSE
  )

  for (k in seq_len(nrow(triu_idx))) {
    i <- triu_idx[k, 1]; j <- triu_idx[k, 2]
    w <- mat[i, j]
    if (abs(w) >= display_threshold) {
      edges_df <- rbind(edges_df, data.frame(
        from   = sub_g_net[i],
        to     = sub_g_net[j],
        weight = abs(w),
        color  = if (w > 0) "red" else "royalblue",
        stringsAsFactors = FALSE
      ))
    }
  }

  n_shown <- nrow(edges_df)

  if (is.null(title)) {
    title <- sprintf("Hub (K=%d) — %d/%d edges\nFilter: %s",
                     K, n_shown, n_pairs, mode_label)
  }

  # ── igraph / ggraph plot ───────────────────────────────────────────────
  nodes_df <- data.frame(name = sub_g_net, stringsAsFactors = FALSE)

  G_ig <- igraph::graph_from_data_frame(
    d        = edges_df[, c("from","to","weight","color")],
    directed = FALSE,
    vertices = nodes_df
  )

  if (requireNamespace("ggraph", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE)) {

    layout_obj <- ggraph::create_layout(G_ig, layout = layout)
    p <- ggraph::ggraph(layout_obj) +
      ggraph::geom_edge_link(ggplot2::aes(width = weight, color = color),
                             alpha = 0.7) +
      ggraph::scale_edge_color_identity() +
      ggraph::scale_edge_width(range = c(0.5, 4)) +
      ggraph::geom_node_point(color = "orange", size = 6) +
      ggraph::geom_node_text(ggplot2::aes(label = name),
                             fontface = "bold", size = 4, repel = TRUE) +
      ggplot2::labs(title = title) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11))

    # Significance annotation
    if (!is.null(perm_pval)) {
      stars <- if      (perm_pval < 0.001) "***"
                else if (perm_pval < 0.01)  "**"
                else if (perm_pval < 0.05)  "*"
                else                        "n.s."
      ann <- sprintf("Hub p=%.3f %s  z=%.2f", perm_pval, stars, perm_zscore %||% NA)
      p <- p + ggplot2::labs(subtitle = ann)
    }

    return(p)

  } else {
    # Fallback: base R igraph plot
    message("Install 'ggraph' for nicer plots. Using base igraph.")
    edge_colors <- igraph::E(G_ig)$color
    igraph::plot.igraph(G_ig,
                        vertex.color = "orange",
                        vertex.size  = 18,
                        vertex.label.font = 2,
                        edge.color   = edge_colors,
                        edge.width   = igraph::E(G_ig)$weight * 3,
                        layout       = igraph::layout_with_fr(G_ig),
                        main         = title)
    invisible(G_ig)
  }
}

# Helper: %||% (null-coalescing operator)
`%||%` <- function(a, b) if (!is.null(a)) a else b
