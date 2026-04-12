# ============================================================================
# TimeSCape v0.2 — Gene Plot Functions (ggplot2, Seurat-native)
# ============================================================================

#' Plot Single Gene with Circadian Fit
#'
#' Reads pre-computed CSVs and returns a ggplot2 object showing the cosine
#' fit, per-ZT means, acrophase marker, and optionally single-cell overlay.
#'
#' @param tmeta        data.frame from build_tmeta_from_seurat(): columns
#'                     zt_str (e.g. "ZT03") and ZT_times (numeric, e.g. 3).
#' @param cust_cells   Cell type name (must match the directory name produced
#'                     by run_timescape, i.e. the sanitised ct_safe string).
#' @param period12     Logical. 12-hr (TRUE) or 24-hr (FALSE, default).
#' @param cust_gene    Gene name string.
#' @param print_scdata Logical. Overlay individual cell expression (default FALSE).
#' @param sce          Seurat object — required when print_scdata = TRUE.
#' @param celltype_col Seurat metadata column for cell type.
#' @param zt_col       Seurat metadata column for ZT string (e.g. "ZT_str").
#' @param use_violin   Logical. Violin (TRUE) or dots scatter (FALSE, default).
#' @param outdir       Directory containing the cell-type subdirectory.
#'
#' @return A ggplot2 object.
#' @export
plot_gene_single <- function(
  tmeta,
  cust_cells,
  period12     = FALSE,
  cust_gene,
  print_scdata = FALSE,
  sce          = NULL,
  celltype_col = "cell_type",
  zt_col       = "ZT_str",
  use_violin   = FALSE,
  outdir       = getwd()
) {
  per_label  <- if (period12) "_period_12_" else "_period_24_"
  period_val <- if (period12) 12 else 24

  # ── Locate CSV files ───────────────────────────────────────────────────────
  ct_safe   <- gsub("[^[:alnum:]_]", "_", trimws(cust_cells))
  ct_outdir <- file.path(outdir)   # outdir already points to the cell-type dir
  fbase     <- paste0(ct_safe, per_label)

  f_stats   <- file.path(ct_outdir, paste0(fbase, "circadian_analysis_all.csv"))
  f_zts     <- file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean.csv"))

  if (!file.exists(f_stats))
    stop("Stats CSV not found:\n  ", f_stats,
         "\nRun the analysis first (Step ④).")
  if (!file.exists(f_zts))
    stop("ZT-means CSV not found:\n  ", f_zts)

  T1 <- read.csv(f_stats, stringsAsFactors = FALSE)
  T2 <- read.csv(f_zts,   stringsAsFactors = FALSE)

  # ── Find gene ─────────────────────────────────────────────────────────────
  gi <- which(T1$Genes == cust_gene)
  if (length(gi) == 0)
    stop(sprintf("Gene '%s' not found in results for '%s'.", cust_gene, ct_safe))
  gi <- gi[1]

  amp_g   <- T1$Amp[gi]
  acro_g  <- T1$Acrophase[gi]
  mesor_g <- T1$Mesor[gi]
  acro_24 <- T1$Acrophase_24[gi]
  pval    <- T1$pvalue[gi]
  pcorr   <- T1$pvalue_corr[gi]

  # ── Identify ZT columns in T2 ──────────────────────────────────────────────
  # T2 columns: Genes, ZT00, ZT03, ...  (zt_str labels)
  zt_cols     <- setdiff(colnames(T2), "Genes")
  # Match column names back to numeric ZT_times via tmeta
  zt_numeric  <- tmeta$ZT_times[match(zt_cols, tmeta$zt_str)]
  # Fall back: try stripping "ZT" prefix
  if (any(is.na(zt_numeric))) {
    zt_numeric <- suppressWarnings(as.numeric(gsub("(?i)^ZT[_ ]?", "", zt_cols, perl=TRUE)))
  }
  Rzts <- as.numeric(T2[gi, zt_cols])

  t0   <- min(zt_numeric, na.rm=TRUE)
  tf   <- max(zt_numeric, na.rm=TRUE)
  tval <- seq(t0, tf, length.out = 300)
  fval <- amp_g * cos(2*pi*(tval - acro_g) / period_val) + mesor_g

  # ── Always count total cells for N in title ───────────────────────────────
  n_cells_type <- 0
  if (!is.null(sce)) {
    meta_ct <- as.character(sce@meta.data[[celltype_col]])
    n_cells_type <- sum(meta_ct == cust_cells, na.rm=TRUE)
    # Also try the sanitised name
    if (n_cells_type == 0) n_cells_type <- sum(meta_ct == ct_safe, na.rm=TRUE)
  }

  # ── Collect per-cell expression when overlay requested ────────────────────
  nc_cum   <- 0
  sc_data  <- NULL

  if (print_scdata && !is.null(sce)) {
    meta <- sce@meta.data
    ct_mask <- (as.character(meta[[celltype_col]]) == cust_cells) |
               (as.character(meta[[celltype_col]]) == ct_safe)

    if (sum(ct_mask) > 0) {
      # Normalize
      X_raw  <- as.matrix(.get_seurat_matrix(sce, use_normalized = FALSE))
      csumsX <- colSums(X_raw); csumsX[csumsX == 0] <- 1
      X_norm <- log1p(t(t(X_raw) / csumsX * 1e4))
      rm(X_raw)

      ig_sc  <- which(rownames(X_norm) == cust_gene)
      if (length(ig_sc) > 0) {
        zt_str_cells <- as.character(meta[[zt_col]])[ct_mask]
        zt_num_cells <- tmeta$ZT_times[match(zt_str_cells, tmeta$zt_str)]
        expr_cells   <- as.numeric(X_norm[ig_sc[1], ct_mask])

        ok   <- !is.na(zt_num_cells)
        sc_data <- data.frame(
          ZT   = zt_num_cells[ok],
          Expr = expr_cells[ok]
        )
        nc_cum <- nrow(sc_data)
      }
    }
  }

  # ── White plot theme ───────────────────────────────────────────────────────
  wt <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill="white", colour=NA),
      plot.background   = ggplot2::element_rect(fill="white", colour=NA),
      panel.grid.major  = ggplot2::element_line(colour="#d4d4d4"),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.text         = ggplot2::element_text(colour="black", size=11),
      axis.title        = ggplot2::element_text(colour="black", size=11),
      plot.title        = ggplot2::element_text(colour="black", size=12, face="bold"),
      legend.background = ggplot2::element_rect(fill="white", colour="#b0b0b0"),
      legend.text       = ggplot2::element_text(colour="black", size=10),
      legend.title      = ggplot2::element_blank(),
      legend.position   = "top",
      legend.direction  = "horizontal"
    )

  # ── Build plot ─────────────────────────────────────────────────────────────
  cos_df  <- data.frame(ZT = tval, Expr = fval)
  mean_df <- data.frame(ZT = zt_numeric, Expr = Rzts)
  col_cos  <- "#2278B5"
  col_mean <- "#D65A0C"
  col_dots <- "#4D7FD1"
  col_acro <- "#7A0000"
  col_viol <- "#99CCEE"

  p <- ggplot2::ggplot() + wt +
    ggplot2::scale_x_continuous(breaks = zt_numeric) +
    ggplot2::labs(
      x     = "Zeitgeber Time (hrs)",
      y     = "Expression (log-normalised)",
      title = sprintf("%s  |  F-test p=%.3g  |  Corr p=%.3g  |  Phase ZT%.1f  |  N=%d",
                      cust_gene, pval, pcorr, acro_24, n_cells_type)
    )

  # ── VIOLIN MODE ────────────────────────────────────────────────────────────
  if (print_scdata && !is.null(sc_data) && nc_cum > 0 && use_violin) {
    y_max <- max(stats::quantile(sc_data$Expr, 0.99) * 1.15,
                 max(Rzts, na.rm=TRUE) * 3)
    sc_data$ZT_f <- factor(sc_data$ZT)
    viol_cols <- setNames(rep(col_viol, length(unique(sc_data$ZT))),
                          as.character(sort(unique(sc_data$ZT))))
    p <- p +
      ggplot2::geom_violin(
        data = sc_data,
        ggplot2::aes(x = ZT, y = Expr, group = ZT_f, fill = "Violin"),
        colour = col_cos, alpha = 0.50, linewidth = 0.4
      ) +
      ggplot2::scale_fill_manual(values = c("Violin" = col_viol)) +
      ggplot2::geom_line(data = cos_df,
                         ggplot2::aes(x=ZT, y=Expr, colour="Cosine fit"),
                         linewidth = 2.2, inherit.aes = FALSE) +
      ggplot2::geom_line(data = mean_df,
                         ggplot2::aes(x=ZT, y=Expr, colour="Mean per ZT"),
                         linewidth = 1.8, inherit.aes = FALSE) +
      ggplot2::geom_point(data = mean_df,
                          ggplot2::aes(x=ZT, y=Expr, colour="Mean per ZT"),
                          size=3.5, inherit.aes = FALSE) +
      ggplot2::geom_vline(xintercept = acro_24,
                          linetype="dashed", colour=col_acro, linewidth=1.0) +
      ggplot2::annotate("text", x=acro_24, y=-Inf,
                        label=sprintf("ZT%.1f", acro_24),
                        vjust=-0.4, hjust=-0.1, colour=col_acro, size=3.5) +
      ggplot2::scale_colour_manual(
        values = c("Cosine fit" = col_cos, "Mean per ZT" = col_mean),
        breaks = c("Violin", "Cosine fit", "Mean per ZT")
      ) +
      ggplot2::ylim(0, y_max)

  # ── DOTS / NORMAL MODE ─────────────────────────────────────────────────────
  } else {
    if (print_scdata && !is.null(sc_data) && nc_cum > 0) {
      p <- p +
        ggplot2::geom_jitter(
          data = sc_data,
          ggplot2::aes(x=ZT, y=Expr, colour="Single cells"),
          width=0.25, size=0.8, alpha=0.22, inherit.aes=FALSE
        )
    }
    p <- p +
      ggplot2::geom_line(data = cos_df,
                         ggplot2::aes(x=ZT, y=Expr, colour="Cosine fit"),
                         linewidth=1.6, inherit.aes=FALSE) +
      ggplot2::geom_line(data = mean_df,
                         ggplot2::aes(x=ZT, y=Expr, colour="Mean per ZT"),
                         linewidth=1.3, inherit.aes=FALSE) +
      ggplot2::geom_point(data = mean_df,
                          ggplot2::aes(x=ZT, y=Expr, colour="Mean per ZT"),
                          size=2.8, inherit.aes=FALSE) +
      ggplot2::geom_vline(xintercept=acro_24,
                          linetype="dashed", colour=col_acro, linewidth=0.9) +
      ggplot2::annotate("text", x=acro_24, y=-Inf,
                        label=sprintf("ZT%.1f", acro_24),
                        vjust=-0.4, hjust=-0.1, colour=col_acro, size=3.5) +
      ggplot2::scale_colour_manual(
        values = c("Cosine fit"   = col_cos,
                   "Mean per ZT"  = col_mean,
                   "Single cells" = col_dots),
        breaks = c("Cosine fit", "Mean per ZT",
                   if (print_scdata && nc_cum > 0) "Single cells",
                   "Acrophase")
      )
  }

  # Dummy legend entry for Acrophase (dashed line swatch)
  p <- p +
    ggplot2::geom_blank(
      ggplot2::aes(colour = "Acrophase"),
      data = data.frame(ZT=NA_real_, Expr=NA_real_)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(
          linetype = c("solid","solid",
                       if (print_scdata && nc_cum > 0) "blank",
                       "dashed"),
          shape    = c(NA, 16,
                       if (print_scdata && nc_cum > 0) 16,
                       NA),
          linewidth= c(1.6,1.3,
                       if (print_scdata && nc_cum > 0) 0,
                       0.9)
        )
      )
    )

  p
}


#' Save Batch Gene Plots as PNG Files
#'
#' Generates one PNG per gene and saves them to a subdirectory of outdir.
#'
#' @param tmeta       data.frame from build_tmeta_from_seurat().
#' @param cust_cells  Sanitised cell type name (directory name).
#' @param plot_type   1 = confident (F-test AND corr p<0.05),
#'                    2 = non-confident, 3 = classical circadian genes.
#' @param period12    Logical.
#' @param outdir      Cell-type output directory (contains the CSVs).
#'
#' @export
save_batch_plots <- function(
  tmeta,
  cust_cells,
  plot_type  = 1,
  period12   = FALSE,
  outdir     = getwd()
) {
  per_label  <- if (period12) "_period_12_" else "_period_24_"
  period_val <- if (period12) 12 else 24

  ct_safe   <- gsub("[^[:alnum:]_]", "_", trimws(cust_cells))
  fbase     <- paste0(ct_safe, per_label)
  f_stats   <- file.path(outdir, paste0(fbase, "circadian_analysis_all.csv"))
  f_zts     <- file.path(outdir, paste0(fbase, "circadian_ZTs_mean.csv"))

  if (!file.exists(f_stats)) stop("Stats CSV not found: ", f_stats)
  T1 <- read.csv(f_stats, stringsAsFactors=FALSE)
  T2 <- read.csv(f_zts,   stringsAsFactors=FALSE)

  zt_cols    <- setdiff(colnames(T2), "Genes")
  zt_numeric <- tmeta$ZT_times[match(zt_cols, tmeta$zt_str)]
  if (any(is.na(zt_numeric)))
    zt_numeric <- suppressWarnings(as.numeric(gsub("(?i)^ZT[_ ]?","",zt_cols,perl=TRUE)))

  classic_circ <- c("arnt","bhlh","clock","cry","dbp","tef","hlf","raf","erk",
                    "mek","ras","mtor","map","ral","akt","hif","kras","myc",
                    "nfkb","per","wnt","nr1d","rev","pik","ror")

  conf   <- (T1$pvalue < 0.05) & (T1$pvalue_corr < 0.05)
  is_cl  <- grepl(paste(classic_circ, collapse="|"), tolower(T1$Genes))

  if (plot_type == 1) {
    gene_idx <- which(conf)
    sub_dir  <- file.path(outdir, paste0(fbase, "plots_confident"))
  } else if (plot_type == 2) {
    gene_idx <- which(!conf)
    sub_dir  <- file.path(outdir, paste0(fbase, "plots_non_confident"))
  } else {
    gene_idx <- which(is_cl)
    sub_dir  <- file.path(outdir, paste0(fbase, "plots_classic_circadian"))
  }

  if (!dir.exists(sub_dir)) dir.create(sub_dir, showWarnings=FALSE, recursive=TRUE)
  message("  Saving ", length(gene_idx), " gene plots → ", sub_dir)

  col_cos  <- "#2278B5"; col_mean <- "#D65A0C"; col_acro <- "#7A0000"
  wt <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill="white", colour=NA),
      plot.background  = ggplot2::element_rect(fill="white", colour=NA),
      panel.grid.major = ggplot2::element_line(colour="#d4d4d4"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(colour="black", size=10),
      axis.title       = ggplot2::element_text(colour="black", size=10),
      plot.title       = ggplot2::element_text(colour="black", size=11, face="bold"),
      legend.position  = "top",
      legend.title     = ggplot2::element_blank()
    )

  for (gi in gene_idx) {
    gene   <- T1$Genes[gi]
    Rzts   <- as.numeric(T2[gi, zt_cols])
    t0     <- min(zt_numeric, na.rm=TRUE)
    tf     <- max(zt_numeric, na.rm=TRUE)
    tval   <- seq(t0, tf, length.out=300)
    fval   <- T1$Amp[gi] * cos(2*pi*(tval - T1$Acrophase[gi]) / period_val) + T1$Mesor[gi]
    acro24 <- T1$Acrophase_24[gi]

    cos_df  <- data.frame(ZT=tval, Expr=fval)
    mean_df <- data.frame(ZT=zt_numeric, Expr=Rzts)

    p <- ggplot2::ggplot() + wt +
      ggplot2::scale_x_continuous(breaks=zt_numeric) +
      ggplot2::labs(
        x     = "Zeitgeber Time (hrs)",
        y     = "Expression (log-normalised)",
        title = sprintf("%s  |  F-test p=%.3g  |  Corr p=%.3g",
                        gene, T1$pvalue[gi], T1$pvalue_corr[gi])
      ) +
      ggplot2::geom_line(data=cos_df,  ggplot2::aes(x=ZT,y=Expr,colour="Cosine fit"),  linewidth=1.6, inherit.aes=FALSE) +
      ggplot2::geom_line(data=mean_df, ggplot2::aes(x=ZT,y=Expr,colour="Mean per ZT"), linewidth=1.3, inherit.aes=FALSE) +
      ggplot2::geom_point(data=mean_df,ggplot2::aes(x=ZT,y=Expr,colour="Mean per ZT"), size=2.8,      inherit.aes=FALSE) +
      ggplot2::geom_vline(xintercept=acro24, linetype="dashed", colour=col_acro, linewidth=0.9) +
      ggplot2::annotate("text", x=acro24, y=-Inf,
                        label=sprintf("ZT%.1f",acro24),
                        vjust=-0.4, hjust=-0.1, colour=col_acro, size=3.2) +
      ggplot2::scale_colour_manual(
        values = c("Cosine fit"="$col_cos", "Mean per ZT"=col_mean),
        breaks = c("Cosine fit","Mean per ZT","Acrophase")
      )

    # Fix colour literal
    p <- p + ggplot2::scale_colour_manual(
        values = c("Cosine fit"  = col_cos,
                   "Mean per ZT" = col_mean),
        breaks = c("Cosine fit", "Mean per ZT")
      )

    out_png <- file.path(sub_dir, paste0("plot_", gene, ".png"))
    ggplot2::ggsave(out_png, p, width=8, height=5, dpi=150, bg="white")
  }
  message("  Done.")
}
