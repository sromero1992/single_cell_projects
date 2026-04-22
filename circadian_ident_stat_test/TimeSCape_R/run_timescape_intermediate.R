# ============================================================================
# TimeSCape v0.2  --  Multi-Cell-Type Analysis Script  [Intermediate stage]
#
# Runs the full circadian + pathway + GRN pipeline for each target cell type.
# Gene sets are pulled ONCE; the main pipeline loops over each target.
#
# Targets
#   Cell type        Metadata column     Notes
#    -- 
#   CD8+ T cells     sub_cell_types2     fine-grained annotation
#   Endothelial      broad_cell_types    broad annotation
#   Tumor            broad_cell_types    broad annotation
#
# Output layout:
#   TimeSCape_output/
#     CD8__T_cells/          <- circadian results, plots, heatmap
#       CD8__T_cells_GRN/    <- per-pathway GRN PNGs
#     Endothelial/
#       Endothelial_GRN/
#     Tumor/
#       Tumor_GRN/
#
# Workflow per cell type:
#   5.  Gene-level cosinor (confident circadian genes)
#   6.  Top-gene plots + heatmap
#   7.  [Gene sets pulled once before loop]
#   8.  Phase-bin ORA  (phyper  --  EnrichR consensus)
#   9.  AUCell scoring on phase-restricted gene sets
#   10. Pathway cosinor + Excel + violin plots
#   11. Per-pathway GRN time series
# ============================================================================

library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)

set.seed(123)

# =============================================================================
# 1. PATHS
# =============================================================================

base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm\\r_pre_process\\TimeSCape_R_tes"
src_path  <- "C:\\Users\\selim\\Documentos\\vscode_working_dir\\single_cell_projects\\circadian_ident_stat_test\\TimeSCape_R\\R"
out_dir   <- file.path(base_dir, "TimeSCape_output_intermediate")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =============================================================================
# 2. SOURCE
# =============================================================================

source(file.path(src_path, "estimate_phaseR.R"))
source(file.path(src_path, "run_timescape.R"))
source(file.path(src_path, "generate_heatmap.R"))
source(file.path(src_path, "plot_gene.R"))
source(file.path(src_path, "pathway_circadian.R"))

# =============================================================================
# 3. LOAD DATA
# =============================================================================

cat("Loading Seurat object...\n")
data <- readRDS(file.path(base_dir,
                          "circadian_mouse_breast_decontX_QC_2_23_2026_subann_intermediate.rds"))
cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(data), nrow(data)))

# =============================================================================
# 4. GLOBAL PARAMETERS  (shared across all cell types)
# =============================================================================

zt_col       <- "ZT_time_str"   # ZT string column ("ZT00", "ZT06", ...)
group_col    <- NULL
custom_zt    <- NULL
custom_group <- NULL
period12     <- FALSE           # 24-hr period
norm_str     <- "logcounts"     # decontX cleaned data in the 'data' slot
n_workers    <- 2L

# Build tmeta once (same ZT structure for all cell types)
tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)

# =============================================================================
# 5. CELL-TYPE TARGETS
#    All unique cell types from sub_cell_types2 are run automatically.
#    To restrict to a subset, set ct_subset to a character vector of labels;
#    leave as NULL to run every cell type found in the column.
# =============================================================================

celltype_col_global <- "sub_cell_types2"
ct_subset           <- NULL   # e.g. c("CD8+ T cells", "NK cells") or NULL = all

all_labels <- sort(unique(na.omit(
  as.character(data@meta.data[[celltype_col_global]])
)))

ct_targets <- if (!is.null(ct_subset)) {
  intersect(ct_subset, all_labels)
} else {
  all_labels
}

cat(sprintf("\nFound %d cell types in '%s':\n", length(ct_targets), celltype_col_global))
for (lbl in ct_targets) {
  n <- sum(as.character(data@meta.data[[celltype_col_global]]) == lbl, na.rm = TRUE)
  cat(sprintf("  %-35s  n=%d\n", lbl, n))
}

# =============================================================================
# 6. PULL GENE SETS  --  done ONCE before the loop
#    KEGG + Reactome + GO:BP pulled together so phyper and EnrichR
#    test the same pathway universe (enables consensus hits).
# =============================================================================

cat("\n-- Pulling gene sets (KEGG + Reactome + GO:BP) --\n")

gs_kegg <- pull_genesets(
  collection  = "C2",
  subcategory = "CP:KEGG_LEGACY",
  species     = "Mus musculus",
  min_size    = 10L,
  max_size    = 500L,
  deduplicate = FALSE
)
cat(sprintf("  KEGG     : %d gene sets\n", length(gs_kegg)))

gs_reactome <- pull_genesets(
  collection  = "C2",
  subcategory = "CP:REACTOME",
  species     = "Mus musculus",
  min_size    = 30L,          # removes small viral/disease-specific sets
  max_size    = 500L,
  deduplicate = FALSE
)
cat(sprintf("  Reactome : %d gene sets\n", length(gs_reactome)))

gs_gobp <- pull_genesets(
  collection  = "C5",
  subcategory = "GO:BP",
  species     = "Mus musculus",
  min_size    = 20L,
  max_size    = 200L,
  deduplicate = TRUE
)
cat(sprintf("  GO:BP    : %d gene sets (after dedup)\n", length(gs_gobp)))

genesets <- c(gs_kegg, gs_reactome, gs_gobp)
cat(sprintf("  Combined : %d gene sets total\n\n", length(genesets)))

# (GRN appearance and clock gene list removed -- GRN not run in this script)

# =============================================================================
# 8. MAIN LOOP  --  runs the full pipeline for every cell-type target
# =============================================================================

for (focus_ct in ct_targets) {

  celltype_col <- celltype_col_global   # all targets use the same column
  focus_safe   <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))
  ct_dir       <- file.path(out_dir, focus_safe)

  cat(sprintf(
    "\n%s\n=== [%s]  col='%s' ===\n%s\n",
    strrep("=", 70), focus_ct, celltype_col, strrep("=", 70)
  ))
  
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)
  
  #  --  8-A. Gene-level circadian analysis  -- 
  cat(sprintf("\n--- 5. Gene-level cosinor: '%s' ---\n", focus_ct))
  
  results <- run_timescape(
    obj             = data,
    celltype_col    = celltype_col,
    zt_col          = zt_col,
    tmeta           = tmeta,
    rm_low_conf     = TRUE,
    period12        = period12,
    custom_celltype = focus_ct,
    plot_heat       = TRUE,
    norm_str        = norm_str,
    outdir          = out_dir,
    group_col       = group_col,
    custom_group    = custom_group,
    custom_zt       = custom_zt,
    n_workers       = n_workers
  )
  
  key       <- names(results)[1]
  T1        <- results[[key]]$T1
  conf_mask <- T1$pvalue < 0.05 & T1$pvalue_corr < 0.05
  conf_genes <- T1$Genes[conf_mask]
  
  cat(sprintf("\n  Confident circadian genes: %d / %d\n",
              length(conf_genes), nrow(T1)))
  cat("  Top 10:\n")
  print(head(T1[conf_mask,
                c("Genes","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))
  
  if (length(conf_genes) == 0) {
    cat(sprintf(
      "  [!] No confident genes found for '%s'  --  skipping downstream steps.\n",
      focus_ct))
    next
  }
  
  #  --  8-B. Top-6 gene grid + heatmap  --
  cat("\n--- 6. Gene plots + heatmap ---\n")

  top_n_show    <- 6L
  genes_to_show <- head(conf_genes, top_n_show)
  gene_plots <- lapply(genes_to_show, function(g) {
    tryCatch(
      plot_gene_single(
        tmeta        = tmeta,
        cust_cells   = key,
        period12     = period12,
        cust_gene    = g,
        print_scdata = TRUE,
        sce          = data,
        celltype_col = celltype_col,
        zt_col       = zt_col,
        use_violin   = FALSE,
        outdir       = ct_dir
      ),
      error = function(e) { message("  Skip ", g, ": ", e$message); NULL }
    )
  })
  gene_plots <- Filter(Negate(is.null), gene_plots)

  if (length(gene_plots) > 0) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
    n_col  <- min(3L, length(gene_plots))
    n_row  <- ceiling(length(gene_plots) / n_col)
    grid_p <- gridExtra::arrangeGrob(grobs = gene_plots, ncol = n_col)
    ggplot2::ggsave(
      filename = file.path(ct_dir,
                           paste0(focus_safe, "_top", top_n_show, "_genes.png")),
      plot   = grid_p,
      width  = 6*n_col, height = 5*n_row, dpi = 150, bg = "white"
    )
    cat(sprintf("  Top-%d gene grid saved.\n", top_n_show))
  }

  # No per-gene batch saving -- use the explorer block at the bottom for
  # individual gene inspection. Batch saving is impractical for large cell types.

  cat("  Generating heatmap...\n")
  tryCatch(
    generate_heatmap(
      celltype = key,
      strict   = TRUE,
      circ     = FALSE,
      period12 = period12,
      outdir   = ct_dir
    ),
    error = function(e) message("  Heatmap skipped: ", e$message)
  )
  
  #  --  8-C. Phase-bin ORA  -- 
  cat(sprintf("\n--- 8. Phase-bin ORA: '%s' ---\n", focus_ct))
  
  phase_results <- phase_bin_analysis(
    T1            = T1,
    conf_mask     = conf_mask,
    genesets      = genesets,
    bin_width     = 3,
    n_top         = 5L,
    min_overlap   = 3L,
    min_bin_genes = 3L,
    p_thresh      = 0.05,
    use_padj      = FALSE,
    use_enrichr   = FALSE,   # set TRUE when online for consensus validation
    enrichr_dbs   = c(
      "KEGG_2026",
      "Reactome_Pathways_2024",
      "GO_Biological_Process_2025"
    )
  )
  
  # Print ORA summary
  cat("\n  Enriched pathways per phase bin:\n")
  for (bin in names(phase_results$ora_results)) {
    df <- phase_results$ora_results[[bin]]
    for (i in seq_len(nrow(df)))
      cat(sprintf("    %s [%s]  %s\n", bin, df$source[i], df$Pathway[i]))
  }
  
  # Save ORA Excel
  ora_xlsx <- file.path(ct_dir, paste0(focus_safe, "_phase_bin_ORA.xlsx"))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    style_consensus <- openxlsx::createStyle(fgFill = "#C6EFCE")
    style_phyper    <- openxlsx::createStyle(fgFill = "#FFEB9C")
    style_enrichr   <- openxlsx::createStyle(fgFill = "#BDD7EE")
    hdr_style       <- openxlsx::createStyle(
      fontColour = "#FFFFFF", fgFill = "#2F4F8F",
      halign = "center", textDecoration = "bold")
    
    all_rows <- do.call(rbind, lapply(names(phase_results$ora_results), function(b)
      cbind(Bin = b, phase_results$ora_results[[b]])))
    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", all_rows, rowNames = FALSE)
    openxlsx::setColWidths(wb, "Summary", cols = seq_len(ncol(all_rows)), widths = "auto")
    openxlsx::addStyle(wb, "Summary", hdr_style, rows = 1L,
                       cols = seq_len(ncol(all_rows)), gridExpand = TRUE)
    for (ri in seq_len(nrow(all_rows))) {
      sty <- switch(all_rows$source[ri],
                    "consensus"    = style_consensus,
                    "phyper_only"  = style_phyper,
                    "enrichr_only" = style_enrichr, NULL)
      if (!is.null(sty))
        openxlsx::addStyle(wb, "Summary", sty, rows = ri + 1L,
                           cols = seq_len(ncol(all_rows)),
                           gridExpand = TRUE, stack = TRUE)
    }
    for (bin in names(phase_results$ora_results)) {
      sname <- gsub("[^[:alnum:]_]", "_", bin)
      df    <- phase_results$ora_results[[bin]]
      openxlsx::addWorksheet(wb, sname)
      openxlsx::writeData(wb, sname, df, rowNames = FALSE)
      openxlsx::setColWidths(wb, sname, cols = seq_len(ncol(df)), widths = "auto")
      openxlsx::addStyle(wb, sname, hdr_style, rows = 1L,
                         cols = seq_len(ncol(df)), gridExpand = TRUE)
      for (ri in seq_len(nrow(df))) {
        sty <- switch(df$source[ri],
                      "consensus"    = style_consensus,
                      "phyper_only"  = style_phyper,
                      "enrichr_only" = style_enrichr, NULL)
        if (!is.null(sty))
          openxlsx::addStyle(wb, sname, sty, rows = ri + 1L,
                             cols = seq_len(ncol(df)),
                             gridExpand = TRUE, stack = TRUE)
      }
    }
    openxlsx::saveWorkbook(wb, ora_xlsx, overwrite = TRUE)
    cat(sprintf("  ORA saved -> %s\n", ora_xlsx))
  }
  
  #  --  8-D. AUCell scoring on phase-restricted gene sets  -- 
  cat(sprintf("\n--- 9. AUCell scoring: '%s' ---\n", focus_ct))
  
  # Cache is per-cell-type (different confident genes each time)
  auc_cache <- file.path(ct_dir,
                         paste0(focus_safe, "_auc_matrix_phase_KEGG_Reactome_GOBP.rds"))
  
  if (file.exists(auc_cache)) {
    cat("  Loading cached AUC matrix...\n")
    auc_phase <- readRDS(auc_cache)
    cat(sprintf("  %d phase gene sets x %d cells\n",
                nrow(auc_phase), ncol(auc_phase)))
  } else {
    cat("  Running AUCell...\n")
    auc_phase <- auc_score_cells(
      obj          = data,
      genesets     = phase_results$phase_gs,
      use_norm     = TRUE,
      auc_max_rank = 0.05,
      n_cores      = 1L,
      min_gs_size  = 3L
    )
    saveRDS(auc_phase, auc_cache)
    cat(sprintf("  Cached -> %s\n", auc_cache))
  }
  
  #  --  8-E. Pathway cosinor  -- 
  cat(sprintf("\n--- 10. Pathway cosinor: '%s' ---\n", focus_ct))
  
  path_phase <- pathway_cosinor(
    auc_mat      = auc_phase,
    meta         = data@meta.data,
    celltype_col = celltype_col,
    zt_col       = zt_col,
    tmeta        = tmeta,
    target_ct    = focus_ct,
    period12     = period12,
    custom_zt    = custom_zt
  )
  
  conf_phase <- path_phase$stats[
    path_phase$stats$pvalue < 0.05 & path_phase$stats$pvalue_corr < 0.05, ]
  cat(sprintf("  Confident phase-pathway sets: %d / %d\n",
              nrow(conf_phase), nrow(path_phase$stats)))
  if (nrow(conf_phase) > 0) {
    cat("  Top 10:\n")
    print(head(conf_phase[,
                          c("Pathway","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))
  }
  
  # Save Excel
  xlsx_phase <- file.path(ct_dir,
                          paste0(focus_safe, "_phase_pathway_circadian.xlsx"))
  write_pathway_results(path_phase, xlsx_phase, celltype = focus_ct)
  cat(sprintf("  Results saved -> %s\n", xlsx_phase))
  
  # Top-6 pathway grid
  path_plot_dir <- file.path(ct_dir, paste0(focus_safe, "_phase_pathway_plots"))
  if (!dir.exists(path_plot_dir)) dir.create(path_plot_dir, recursive = TRUE)
  
  if (nrow(conf_phase) > 0) {
    top_pw   <- head(conf_phase$Pathway, 6L)
    pw_plots <- lapply(top_pw, function(pw) {
      tryCatch(
        plot_pathway_single(
          auc_mat        = auc_phase,
          path_results   = path_phase,
          meta           = data@meta.data,
          celltype_col   = celltype_col,
          zt_col         = zt_col,
          tmeta          = tmeta,
          target_ct      = focus_ct,
          target_pathway = pw,
          period12       = period12,
          use_violin     = TRUE
        ),
        error = function(e) { message("  Skip ", pw, ": ", e$message); NULL }
      )
    })
    pw_plots <- Filter(Negate(is.null), pw_plots)
    
    if (length(pw_plots) > 0) {
      n_col   <- min(3L, length(pw_plots))
      n_row   <- ceiling(length(pw_plots) / n_col)
      grid_pw <- gridExtra::arrangeGrob(grobs = pw_plots, ncol = n_col)
      ggplot2::ggsave(
        filename = file.path(path_plot_dir,
                             paste0(focus_safe, "_top_phase_pathways_grid.png")),
        plot   = grid_pw,
        width  = 6*n_col, height = 5*n_row, dpi = 150, bg = "white"
      )
    }
  }
  
  # Batch save all confident pathway plots
  cat("  Batch saving pathway plots...\n")
  save_batch_pathway_plots(
    auc_mat      = auc_phase,
    path_results = path_phase,
    meta         = data@meta.data,
    celltype_col = celltype_col,
    zt_col       = zt_col,
    tmeta        = tmeta,
    target_ct    = focus_ct,
    n_top        = 20L,
    period12     = period12,
    use_violin   = TRUE,
    outdir       = path_plot_dir
  )
  
  cat(sprintf("\n  [DONE] '%s'\n", focus_ct))
  
} # end loop over ct_targets

cat("\n", strrep("=", 70), "\n")
cat("All cell types complete.\n")
cat(sprintf("Output root: %s\n", out_dir))

# =============================================================================
# INTERACTIVE SINGLE-GENE EXPLORER
# Run this block manually after the loop for any cell type you want to inspect.
# Set focus_ct / celltype_col / ct_dir / key / T1 / conf_mask / conf_genes
# to the values for the cell type of interest, then run the block.
#
# Example  --  switch to Endothelial:
#   focus_ct     <- "Endothelial"
#   celltype_col <- "broad_cell_types"
#   focus_safe   <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))
#   ct_dir       <- file.path(out_dir, focus_safe)
#   key          <- focus_safe
#   T1           <- read.csv(file.path(ct_dir,
#                     paste0(key, "_period_24_circadian_analysis_all.csv")))
#   conf_mask    <- T1$pvalue < 0.05 & T1$pvalue_corr < 0.05
#   conf_genes   <- T1$Genes[conf_mask]
# =============================================================================

#  --  GENE & DATA  -- 
# explore_gene   <- "Per2"
# explore_violin <- FALSE
# explore_sc     <- TRUE
# explore_save   <- TRUE
#
#  --  APPEARANCE  -- 
# title_size     <- 13
# title_bold     <- TRUE
# title_hjust    <- 0.5
# axis_title_sz  <- 12
# axis_text_sz   <- 11
# legend_text_sz <- 10
# legend_pos     <- "top"
# y_min          <- NULL
# y_max          <- NULL
# fig_width      <- 9
# fig_height     <- 5.5
# fig_dpi        <- 180
#
# if (!explore_gene %in% T1$Genes) {
#   cat(sprintf("  '%s' not found in T1.\n", explore_gene))
# } else {
#   gr <- T1[T1$Genes == explore_gene, ]
#   cat(sprintf("\n[%s]  Amp=%.3f  Mesor=%.3f  Acrophase=%.1fh  p=%.4g  padj=%.4g  %s\n",
#     explore_gene, gr$Abs_Amp, gr$Mesor, gr$Acrophase_24, gr$pvalue, gr$pvalue_corr,
#     ifelse(explore_gene %in% conf_genes, "** CONFIDENT **", "(not significant)")))
#   explore_plot <- tryCatch(
#     plot_gene_single(
#       tmeta=tmeta, cust_cells=key, period12=period12, cust_gene=explore_gene,
#       print_scdata=explore_sc, sce=data, celltype_col=celltype_col, zt_col=zt_col,
#       use_violin=explore_violin, outdir=ct_dir),
#     error=function(e){message("  Plot error: ",conditionMessage(e)); NULL})
#   if (!is.null(explore_plot)) {
#     explore_plot <- explore_plot +
#       ggplot2::theme(
#         plot.title=ggplot2::element_text(size=title_size,
#           face=if(title_bold)"bold"else"plain", hjust=title_hjust),
#         axis.title=ggplot2::element_text(size=axis_title_sz),
#         axis.text=ggplot2::element_text(size=axis_text_sz),
#         legend.text=ggplot2::element_text(size=legend_text_sz),
#         legend.position=legend_pos)
#     if (!is.null(y_min)||!is.null(y_max))
#       explore_plot <- explore_plot +
#         ggplot2::coord_cartesian(ylim=c(
#           if(!is.null(y_min))y_min else NA_real_,
#           if(!is.null(y_max))y_max else NA_real_))
#     suppressWarnings(print(explore_plot))
#     if (explore_save) {
#       out_file <- file.path(ct_dir,
#         sprintf("%s_%s_explorer.png", focus_safe, explore_gene))
#       ggplot2::ggsave(out_file, explore_plot,
#         width=fig_width, height=fig_height, dpi=fig_dpi, bg="white")
#       cat(sprintf("  Saved -> %s\n", out_file))
#     }
#   }
# }