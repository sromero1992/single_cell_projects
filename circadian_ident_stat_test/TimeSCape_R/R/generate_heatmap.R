# ============================================================================
# TimeSCape v0.2 — Heatmap Generation
# ============================================================================

#' Generate Circadian Heatmap from Analysis Results
#'
#' Reads CSVs written by run_timescape(), filters/sorts genes, and renders a
#' blue-white-red z-score heatmap (genes × ZT time points) equivalent to the
#' MATLAB generateHeatmap_circ_simple().
#'
#' @param celltype    Sanitised cell type name (directory + file-name prefix).
#' @param strict      Logical. TRUE (default) = keep only genes with both
#'                    F-test p<0.05 AND corr p<0.05.
#' @param custom_name Optional suffix appended to the PNG filename.
#' @param circ        Logical. Restrict to classical circadian gene prefixes.
#' @param period12    Logical. 12-hr (TRUE) or 24-hr (FALSE, default).
#' @param outdir      Directory containing the CSVs (cell-type subdirectory).
#' @param return_obj  Logical. If TRUE, return the pheatmap object invisibly
#'                    (for inline rendering in Shiny). Default FALSE.
#'
#' @return Invisibly: the pheatmap list (when return_obj=TRUE), else NULL.
#' @export
generate_heatmap <- function(
  celltype,
  strict      = TRUE,
  custom_name = "",
  circ        = FALSE,
  period12    = FALSE,
  outdir      = getwd(),
  return_obj  = FALSE
) {
  per_label  <- if (period12) "_period_12_" else "_period_24_"
  fbase      <- paste0(celltype, per_label)

  f_stats <- file.path(outdir, paste0(fbase, "circadian_analysis_all.csv"))
  f_zts   <- file.path(outdir, paste0(fbase, "circadian_ZTs_mean.csv"))

  if (!file.exists(f_stats) || !file.exists(f_zts)) {
    warning("Analysis files not found:\n  ", f_stats, "\nRun the analysis first.")
    return(invisible(NULL))
  }

  D    <- read.csv(f_stats, stringsAsFactors = FALSE)
  Dzts <- read.csv(f_zts,   stringsAsFactors = FALSE)

  # ── Sort by acrophase then amplitude ──────────────────────────────────────
  ord  <- order(D$Acrophase_24, -D$Abs_Amp)
  D    <- D[ord, ]
  Dzts <- Dzts[ord, ]

  # ── Filter to confident genes ──────────────────────────────────────────────
  if (strict) {
    keep <- (D$pvalue < 0.05) & (D$pvalue_corr < 0.05)
    D    <- D[keep, ]
    Dzts <- Dzts[keep, ]
  }

  # ── Correct negative amplitudes ───────────────────────────────────────────
  neg <- D$Amp < 0
  D$Amp[neg]          <- -D$Amp[neg]
  D$Acrophase_24[neg] <- (D$Acrophase_24[neg] + 12) %% 24
  # Re-sort after correction
  ord2 <- order(D$Acrophase_24, -D$Abs_Amp)
  D    <- D[ord2, ]
  Dzts <- Dzts[ord2, ]

  # ── Restrict to classical circadian genes ────────────────────────────────
  if (circ) {
    circ_pats <- paste0("(?i)(arnt|bhlh|clock|cry|dbp|tef|hlf|raf|erk|mek|",
                        "ras|mtor|map|ral|akt|hif|kras|myc|nfkb|per|wnt|",
                        "nr1d|rev|pik|ror)")
    keep_c <- grepl(circ_pats, D$Genes, perl=TRUE)
    D    <- D[keep_c, ]
    Dzts <- Dzts[keep_c, ]
  }

  if (nrow(D) <= 1) {
    message("Too few genes after filtering (", nrow(D), ") – skipping heatmap.")
    return(invisible(NULL))
  }

  # ── Build expression matrix: rows = genes, cols = ZT labels ──────────────
  zt_cols    <- setdiff(colnames(Dzts), "Genes")   # e.g. "ZT00","ZT03",...
  mat        <- as.matrix(Dzts[, zt_cols])          # genes × ZT
  rownames(mat) <- D$Genes

  # Z-score each gene across ZT time points (row-wise)
  mat_z <- t(scale(t(mat)))
  mat_z[is.nan(mat_z)] <- 0

  # ── Colour palette: blue → white → red ───────────────────────────────────
  cmap <- grDevices::colorRampPalette(c("#2166AC","white","#B2182B"))(256)

  # ── Row annotation: acrophase bar ────────────────────────────────────────
  row_ann <- data.frame(
    Acrophase = D$Acrophase_24,
    row.names = D$Genes
  )
  ann_colors <- list(
    Acrophase = grDevices::colorRampPalette(c("#3a1f72","#e05c2f"))(50)
  )

  # ── pheatmap arguments ────────────────────────────────────────────────────
  fontsize_row <- if (nrow(mat_z) <= 60) max(6, min(9, floor(280/nrow(mat_z)))) else 0
  show_rownames <- nrow(mat_z) <= 60

  title_str <- sprintf("%s  –  %d genes  (z-score, sorted by acrophase)",
                       celltype, nrow(mat_z))

  ph_args <- list(
    mat              = mat_z,
    color            = cmap,
    cluster_rows     = FALSE,
    cluster_cols     = FALSE,
    scale            = "none",         # already z-scored
    show_rownames    = show_rownames,
    fontsize_row     = fontsize_row,
    fontsize_col     = 9,
    border_color     = NA,
    main             = title_str,
    annotation_row   = row_ann,
    annotation_colors= ann_colors,
    silent           = TRUE            # suppress console output
  )

  # ── Save PNG to outdir ────────────────────────────────────────────────────
  sfx  <- if (nchar(custom_name) > 0) paste0("_", custom_name) else ""
  sfx2 <- if (circ) "_core_circ" else ""
  sfx3 <- if (strict) "_strict" else "_all"
  out_png <- file.path(outdir, paste0(fbase, "heatmap", sfx3, sfx2, sfx, ".png"))
  out_csv <- file.path(outdir, paste0(fbase, "heatmap", sfx3, sfx2, sfx, "_genes.csv"))

  h <- do.call(pheatmap::pheatmap, c(ph_args,
    list(filename = out_png,
         width    = 10,
         height   = max(4, min(18, nrow(mat_z) * 0.14)))))

  write.csv(D, out_csv, row.names=FALSE)
  message("Heatmap saved → ", out_png)

  # ── Return object for inline Shiny rendering ──────────────────────────────
  if (return_obj) {
    ph_args[["silent"]] <- FALSE
    ph_out <- do.call(pheatmap::pheatmap, ph_args)
    return(invisible(ph_out))
  }

  invisible(NULL)
}
