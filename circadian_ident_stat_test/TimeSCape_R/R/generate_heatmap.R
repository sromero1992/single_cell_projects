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


# ============================================================================
#' Load Stage-Level Cosinor Results from Disk (Post-Processing Helper)
#'
#' Reads the \code{*_circadian_analysis_all.csv} files written by
#' \code{run_timescape()} for every cell type in one cancer stage and returns
#' a named list compatible with \code{plot_clock_acrophase()}.
#'
#' Use this after the pipeline has completed to regenerate or customise the
#' acrophase polar plot without re-running the analysis. The function recovers
#' original cell type names (with spaces and special characters) from the
#' \code{run_log_{stage}.csv} file so that \code{cell_group_rules} in
#' \code{plot_clock_acrophase()} can use the same names you see in RStudio.
#'
#' @param out_dir   Root output directory for the stage (the folder that
#'   contains per-cell-type subdirectories, e.g.
#'   \code{"Z:/data/TimeSCape_output_early"}).
#' @param period12  Logical. Must match the \code{period12} value used during
#'   the original pipeline run. Default \code{FALSE} (24-hr).
#' @param skip_failed Logical. If \code{TRUE} (default), cell types marked
#'   \code{Skipped = TRUE} in the run log are excluded. Set to \code{FALSE}
#'   to attempt loading all directories regardless.
#'
#' @return A named list: original_cell_type_name -> data.frame of cosinor
#'   stats (same columns as T1 from \code{run_timescape()}).
#'   Returns an empty list with a warning if no files are found.
#'
#' @examples
#' \dontrun{
#' # After run_full_pipeline_early.R has finished:
#' clock_results <- load_stage_results(
#'   out_dir  = "Z:/data/TimeSCape_output_early",
#'   period12 = FALSE
#' )
#'
#' # Inspect which genes are available across all cell types
#' all_genes <- sort(unique(unlist(lapply(clock_results, `[[`, "Genes"))))
#' head(all_genes, 20)
#'
#' # Custom plot: pick your own genes + cell-type grouping
#' plot_clock_acrophase(
#'   results_list     = clock_results,
#'   stage            = "early",
#'   outfile          = "Z:/data/TimeSCape_output_early/clock_acrophase_custom.png",
#'   gene_list        = c("Per2", "Nr1d1", "Cry1", "Dbp", "Arntl"),
#'   cell_group_rules = list(
#'     "T cells"  = c("CD8+ T cells", "CD4+ T cells", "T-Reg cells"),
#'     "Myeloid"  = c("M1 macrophages", "M2 macrophages", "Monocytes"),
#'     "Stromal"  = c("Endothelial", "Fibroblasts"),
#'     "Tumor"    = c("Tumor cells")
#'   ),
#'   strict = TRUE,
#'   dpi    = 300
#' )
#' }
#' @export
load_stage_results <- function(out_dir,
                               period12     = FALSE,
                               skip_failed  = TRUE) {

  per_label <- if (period12) "_period_12_" else "_period_24_"

  if (!dir.exists(out_dir))
    stop("out_dir does not exist: ", out_dir)

  results_list <- list()

  # ── Strategy 1: use run_log to recover original cell type names ──────────────
  log_files <- list.files(out_dir, pattern = "^run_log_.*\\.csv$",
                          full.names = TRUE)

  if (length(log_files) > 0) {
    run_log <- read.csv(log_files[1], stringsAsFactors = FALSE)

    if (skip_failed && "Skipped" %in% names(run_log))
      run_log <- run_log[!run_log$Skipped, ]

    for (i in seq_len(nrow(run_log))) {
      ct_name <- run_log$CellType[i]
      ct_safe <- gsub("[^[:alnum:]_]", "_", trimws(ct_name))
      ct_dir  <- file.path(out_dir, ct_safe)
      f_csv   <- file.path(ct_dir,
                   paste0(ct_safe, per_label, "circadian_analysis_all.csv"))

      if (!file.exists(f_csv)) {
        message("  load_stage_results: CSV not found for '", ct_name,
                "' -- skipping.")
        next
      }

      df <- read.csv(f_csv, stringsAsFactors = FALSE)
      if (nrow(df) == 0) next

      results_list[[ct_name]] <- df
      message(sprintf("  Loaded %-30s  %d genes", ct_name, nrow(df)))
    }

  } else {
    # ── Strategy 2: no run_log — scan subdirectories directly ─────────────────
    message("  load_stage_results: run_log not found; scanning subdirectories.")

    ct_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = FALSE)

    for (ct_safe in ct_dirs) {
      ct_dir <- file.path(out_dir, ct_safe)
      f_csv  <- file.path(ct_dir,
                  paste0(ct_safe, per_label, "circadian_analysis_all.csv"))

      if (!file.exists(f_csv)) next

      df <- read.csv(f_csv, stringsAsFactors = FALSE)
      if (nrow(df) == 0) next

      # Use directory name as key (spaces not recoverable without run_log)
      results_list[[ct_safe]] <- df
      message(sprintf("  Loaded %-30s  %d genes", ct_safe, nrow(df)))
    }
  }

  if (length(results_list) == 0)
    warning("load_stage_results: no analysis CSV files found in ", out_dir)
  else
    message(sprintf("\n  Total: %d cell types loaded from %s",
                    length(results_list), out_dir))

  results_list
}


# ============================================================================
#' Plot Gene Acrophase — Polar Chart (All Cell Types × Stage)
#'
#' Takes accumulated gene-level cosinor results from all cell types in one
#' cancer stage and renders a polar / circular plot showing where in the 24-hr
#' cycle each gene peaks, per cell type and cell group.
#'
#' Call this AFTER the cell-type loop completes in the production scripts.
#'
#' @param results_list   Named list: cell_type_label -> T1 data.frame from
#'   run_timescape(). Each T1 must contain columns: Genes, pvalue,
#'   pvalue_corr, Acrophase_24.
#' @param stage          Character label used in the plot title
#'   (e.g. "early", "intermediate", "advanced").
#' @param outfile        Full path for the saved PNG (including .png extension).
#' @param dpi            Resolution for PNG output. Default 300.
#' @param strict         If TRUE (default), only genes with pvalue < 0.05 AND
#'   pvalue_corr < 0.05 are plotted.
#' @param gene_list      Character vector of gene symbols to plot.
#'   NULL (default) uses the canonical core clock gene set:
#'   Clock, Arntl, Npas2, Per1, Per2, Cry1, Cry2, Nr1d1, Nr1d2, Rorc, Dbp.
#' @param cell_group_rules Named list mapping group label -> character vector
#'   of cell type names belonging to that group. Cell types not matched by any
#'   rule are placed in "Other". NULL (default) uses the built-in immune /
#'   tumour / structural classification.
#'   Example:
#'   \code{list(
#'     "T cells"    = c("CD8+ T cells", "CD4+ T cells", "T-Reg cells"),
#'     "Myeloid"    = c("M1 macrophages", "M2 macrophages", "Monocytes"),
#'     "Epithelial" = c("Tumor cells", "Epithelial cells")
#'   )}
#'
#' @return Invisibly: the ggplot object.
#' @export
plot_clock_acrophase <- function(
  results_list,
  stage,
  outfile,
  dpi              = 300,
  strict           = TRUE,
  gene_list        = NULL,
  cell_group_rules = NULL
) {

  # ── 0. Packages ──────────────────────────────────────────────────────────────
  for (pkg in c("ggplot2", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Install %s:  install.packages('%s')", pkg, pkg))
  }

  # ── 1. Resolve gene list ──────────────────────────────────────────────────────
  default_clock_genes <- c("Clock", "Arntl", "Npas2",
                            "Per1",  "Per2",
                            "Cry1",  "Cry2",
                            "Nr1d1", "Nr1d2",
                            "Rorc",  "Dbp")

  target_genes <- if (is.null(gene_list)) default_clock_genes else unique(gene_list)

  plot_title_genes <- if (is.null(gene_list)) "Core Clock Gene" else "Gene"

  # ── 2. Resolve cell group rules ───────────────────────────────────────────────
  default_cell_groups <- list(
    "Adaptive Immune"       = c("CD8+ T cells", "CD4+ T cells", "T-Reg cells",
                                 "NK cells", "B cells", "γδ T cells", "DN T cells"),
    "Innate Immune"         = c("M0 macrophages", "M1 macrophages", "M2 macrophages",
                                 "Mast cells", "Neutrophils", "Monocytes", "pDCs", "DCs"),
    "Tumor"                 = c("Tumor cells"),
    "Structural/Connective" = c("Endothelial", "Epithelial cells", "Fibroblasts",
                                 "Adipocytes", "Myocytes", "Basal KC-like")
  )

  group_rules <- if (is.null(cell_group_rules)) default_cell_groups else cell_group_rules

  # Validate: must be a named list of character vectors
  if (!is.list(group_rules) || is.null(names(group_rules)) ||
      any(names(group_rules) == "")) {
    stop("cell_group_rules must be a named list, e.g. list('T cells' = c('CD8+ T cells', ...))")
  }

  group_labels <- names(group_rules)

  # ── 3. Combine all cell-type T1 tables ───────────────────────────────────────
  if (length(results_list) == 0) {
    message("plot_clock_acrophase: results_list is empty -- skipping.")
    return(invisible(NULL))
  }

  combined_df <- dplyr::bind_rows(
    lapply(names(results_list), function(ct) {
      df <- results_list[[ct]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$Cell_Type <- ct
      df
    })
  )

  if (nrow(combined_df) == 0) {
    message("plot_clock_acrophase: no data after combining -- skipping.")
    return(invisible(NULL))
  }

  # ── 4. Filter to target genes ─────────────────────────────────────────────────
  genes_present <- intersect(combined_df$Genes, target_genes)
  if (length(genes_present) == 0) {
    message("plot_clock_acrophase: none of the requested genes found in results -- skipping.")
    message("  Genes requested: ", paste(target_genes, collapse = ", "))
    return(invisible(NULL))
  }

  df_clock <- dplyr::filter(combined_df, Genes %in% genes_present)

  # ── 5. Confidence filter ──────────────────────────────────────────────────────
  if (strict) {
    df_clock <- dplyr::filter(df_clock, pvalue < 0.05, pvalue_corr < 0.05)
  }

  if (nrow(df_clock) == 0) {
    message("plot_clock_acrophase: no ", if (strict) "confident " else "",
            "hits for the requested genes -- skipping.")
    if (strict)
      message("  Tip: set strict = FALSE to include all fitted genes.")
    return(invisible(NULL))
  }

  # ── 6. Classify cell types into groups ───────────────────────────────────────
  classify_group <- function(ct) {
    for (grp in group_labels) {
      if (ct %in% group_rules[[grp]]) return(grp)
    }
    "Other"
  }

  df_clock$Acrophase_24 <- as.numeric(df_clock$Acrophase_24)
  df_clock$Cell_Group   <- factor(
    vapply(df_clock$Cell_Type, classify_group, character(1L)),
    levels = c(group_labels, "Other")
  )

  # Radial position: one ring per gene, sorted alphabetically
  df_clock$Base_Radius <- as.numeric(factor(
    df_clock$Genes, levels = sort(unique(df_clock$Genes))
  ))

  # ── 7. Colours — one per gene ─────────────────────────────────────────────────
  gene_order   <- sort(unique(df_clock$Genes))
  n_genes      <- length(gene_order)

  base_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
                    "#b4e946","#8c564b","#e377c2","#bcbd22","#17becf","#c5b0d5")
  if (n_genes > length(base_palette))
    base_palette <- grDevices::colorRampPalette(base_palette)(n_genes)

  gene_colors <- stats::setNames(base_palette[seq_len(n_genes)], gene_order)

  # ── 8. Shapes — one per group, auto-assigned from a fixed pool ───────────────
  shape_pool  <- c(16, 17, 15, 18, 4, 8, 10, 11, 12, 13, 14)
  all_groups  <- c(group_labels, "Other")
  # Only keep groups that actually appear in the data
  groups_used <- intersect(all_groups, unique(as.character(df_clock$Cell_Group)))
  group_shapes <- stats::setNames(
    shape_pool[seq_along(groups_used)],
    groups_used
  )

  max_radius <- max(df_clock$Base_Radius) + 1.0
  all_breaks <- sort(unique(df_clock$Base_Radius))

  # ── 9. Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df_clock,
         ggplot2::aes(x     = Acrophase_24,
                      y     = Base_Radius,
                      color = Genes,
                      shape = Cell_Group)) +

    ggplot2::geom_hline(ggplot2::aes(yintercept = max_radius),
                        color = "black", linewidth = 0.8, linetype = "solid") +

    ggplot2::geom_point(size = 3, alpha = 0.85) +

    ggplot2::coord_polar(theta = "x", start = 0) +

    ggplot2::scale_x_continuous(
      limits = c(0, 24),
      breaks = seq(0, 22, by = 2),
      labels = as.character(seq(0, 22, by = 2))
    ) +

    ggplot2::scale_y_continuous(
      limits = c(0, max_radius),
      breaks = all_breaks,
      labels = rep("", length(all_breaks))
    ) +

    ggplot2::scale_color_manual(values = gene_colors,  name = "Gene") +
    ggplot2::scale_shape_manual(values = group_shapes,  name = "Cell Group") +

    ggplot2::labs(
      title    = sprintf("%s Acrophase — %s stage",
                         plot_title_genes, tools::toTitleCase(stage)),
      subtitle = sprintf(
        "%d gene × cell-type combinations  |  n genes: %d  |  confident only: %s",
        nrow(df_clock), n_genes, if (strict) "yes" else "no"
      ),
      x = "Zeitgeber Time (hrs)",
      y = ""
    ) +

    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "white", color = NA),
      plot.background    = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major.x = ggplot2::element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey80", linetype = "dotted"),
      axis.text.x        = ggplot2::element_text(size = 11, colour = "black"),
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(size = 10, colour = "grey40", hjust = 0.5),
      legend.position    = "right",
      legend.title       = ggplot2::element_text(size = 11, face = "bold"),
      legend.text        = ggplot2::element_text(size = 10)
    )

  # ── 10. Save ──────────────────────────────────────────────────────────────────
  if (!dir.exists(dirname(outfile)))
    dir.create(dirname(outfile), recursive = TRUE)

  ggplot2::ggsave(outfile, plot = p,
                  width = 8, height = 7, units = "in",
                  dpi = dpi, bg = "white")

  message(sprintf("Acrophase plot saved -> %s", outfile))
  invisible(p)
}
