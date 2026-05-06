# ============================================================================
# TimeSCape v0.2 — Complete Analysis Script
# Run interactively in RStudio, or: Rscript run_timescape_test.R
#
# Workflow:
#   1-2.  Paths + source
#   3.    Load data & inspect
#   4.    Set parameters
#   5.    Run gene-level circadian (all cell types or one for testing)
#   6.    Plot top genes + heatmap
#   7.    Choose + pull gene sets (used as ORA background)
#   8.    Phase-bin enrichment (bin genes by acrophase -> ORA -> phase-restricted sets)
#   9.    AUCell scoring on phase-restricted gene sets
#   10.   Circadian cosinor on phase scores -> Excel + plots
#   11.   GRN time series using top confident phase-pathway
# ============================================================================

library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)

set.seed(123)

# -- 1. PATHS ------------------------------------------------------------------

base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm_organized\\r_pre_process\\TimeSCape_R_tes"
src_path  <- "C:\\Users\\selim\\Documentos\\vscode_working_dir\\single_cell_projects\\circadian_ident_stat_test\\TimeSCape_R\\R"
out_dir   <- file.path(base_dir, "TimeSCape_output")

# -- 2. SOURCE ----------------------------------------------------------------

source(file.path(src_path, "estimate_phaseR.R"))
source(file.path(src_path, "run_timescape.R"))
source(file.path(src_path, "generate_heatmap.R"))
source(file.path(src_path, "plot_gene.R"))
source(file.path(src_path, "pathway_circadian.R"))

# -- 3. LOAD DATA -------------------------------------------------------------

cat("Loading Seurat object...\n")
data <- readRDS(file.path(base_dir,
  "circadian_mouse_breast_decontX_QC_2_23_2026_subann_early.rds"))
cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(data), nrow(data)))

# -- 4. PARAMETERS ------------------------------------------------------------

celltype_col <- "CellType"      # broad/fine cell type column
zt_col       <- "ZT_time_str"   # ZT string column ("ZT00", "ZT06", ...)
group_col    <- NULL            # secondary grouping, e.g. "tumor_stage"; NULL = off
custom_zt    <- NULL            # filter ZT points; NULL = all
custom_group <- NULL            # filter group values; NULL = all
period12     <- FALSE           # FALSE = 24-hr circadian; TRUE = 12-hr

# norm_str:
#   "logcounts" -> use 'data' slot (NormalizeData / SCTransform / decontX output)
#   "lib_size"  -> re-normalise from raw 'counts' slot
#   "none"      -> use counts as-is
norm_str  <- "logcounts"   # decontX stores cleaned data in the 'data' slot

# Workers:
#   1 = sequential (safest for dense matrices)
#   2 = parallel via future_lapply (only sparse slice sent to workers)
n_workers <- 2L

# Cell type to focus on for pathway + GRN analysis
focus_ct  <- "CD8+ T cells"

# Build tmeta (ZT strings -> numeric hours)
tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)

cat(sprintf("\nCell types in '%s':\n", celltype_col))
print(sort(unique(data@meta.data[[celltype_col]])))

# -- 5. GENE-LEVEL CIRCADIAN ANALYSIS -----------------------------------------
# Set custom_celltype = focus_ct for a quick test.
# Set to NULL to run all cell types.

cat(sprintf("\n=== Running circadian analysis: '%s' ===\n", focus_ct))
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

# Results for focus cell type
key    <- names(results)[1]
T1     <- results[[key]]$T1
ct_dir <- file.path(out_dir, key)

# ── Output subfolders (mirrors pipeline structure) ─────────────────────────────
dir_gene <- file.path(ct_dir, "01_gene_circadian")   # gene plots, heatmap, ORA
dir_pwA  <- file.path(ct_dir, "02_pathway_A")         # phase-bin ORA + AUCell A
dir_pwB  <- file.path(ct_dir, "03_pathway_B")         # full-pathway AUCell B
for (d in c(ct_dir, dir_gene, dir_pwA, dir_pwB))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ── Core clock gene reference (used by GRN hub selection) ─────────────────────
clock_genes_ref <- c(
  "Arntl", "Bmal1", "Clock", "Npas2",        # positive arm
  "Per1",  "Per2",  "Per3",                    # period genes
  "Cry1",  "Cry2",                             # cryptochrome genes
  "Nr1d1", "Nr1d2",                            # REV-ERB alpha/beta
  "Rora",  "Rorb",  "Rorc",                   # ROR alpha/beta/gamma
  "Dbp",   "Tef",   "Hlf",   "Nfil3",        # PAR/D-box TFs
  "Ciart", "Bhlhe40", "Bhlhe41",              # DEC1/DEC2
  "Timeless", "Csnk1d", "Csnk1e"             # CK1 delta/epsilon
)

conf_mask  <- T1$pvalue < 0.05 & T1$pvalue_corr < 0.05
conf_genes <- T1$Genes[conf_mask]
focus_safe <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))

cat(sprintf("\nConfident circadian genes in '%s': %d / %d\n",
            focus_ct, length(conf_genes), nrow(T1)))
cat("Top 10:\n")
print(head(T1[conf_mask, c("Genes","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))

# -- 6. GENE PLOTS + HEATMAP --------------------------------------------------
# 6a. Interactive grid of top 6 confident genes
top_n_show <- 6L
sc_style   <- "violin"    # "dots" = jitter  |  "violin" = violin density

if (length(conf_genes) > 0) {
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
        use_violin   = (sc_style == "violin"),
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
    print(gridExtra::grid.arrange(grid_p))
    ggplot2::ggsave(
      filename = file.path(dir_gene, paste0(key, "_top", top_n_show, "_genes.png")),
      plot     = grid_p,
      width    = 6*n_col, height = 5*n_row, dpi = 150, bg = "white"
    )
  }
}

# 6b. Heatmap (confident genes sorted by acrophase)
# NOTE: batch gene plots removed — use the 6d explorer to inspect individual
# genes interactively. Batch-saving hundreds/thousands of PNGs is too slow
# for large cell types and the heatmap + explorer cover the same need.
cat("\nGenerating heatmap...\n")
generate_heatmap(
  celltype = key,
  strict   = TRUE,
  circ     = FALSE,
  period12 = period12,
  outdir   = ct_dir
)

# 6c. Clock gene acrophase polar plot (single cell type)
# Uses T1 already in memory — no re-running needed.
#
# ── Customisation options ────────────────────────────────────────────────────
# gene_list (default NULL = built-in core clock genes):
#   Pass a character vector to plot any genes instead of the clock set, e.g.:
#     gene_list = c("Per2", "Nr1d1", "Dbp", "Tsc22d3", "Cxcr4")
#
# cell_group_rules (default NULL = built-in immune/tumour/structural bins):
#   Pass a named list to define your own cell-type groupings, e.g.:
#     cell_group_rules = list(
#       "T cells"   = c("CD8+ T cells", "CD4+ T cells", "T-Reg cells"),
#       "Myeloid"   = c("M1 macrophages", "M2 macrophages", "Monocytes"),
#       "Stromal"   = c("Endothelial", "Fibroblasts"),
#       "Malignant" = c("Tumor cells")
#     )
#   Cell types not matched by any rule are placed in "Other" automatically.
# ─────────────────────────────────────────────────────────────────────────────
cat("\nGenerating clock gene acrophase polar plot...\n")
tryCatch(
  plot_clock_acrophase(
    results_list = stats::setNames(list(T1), focus_ct),
    stage        = "test",
    outfile      = file.path(dir_gene,
                     paste0(focus_safe, "_clock_acrophase_test.png")),
    dpi          = 300,
    strict       = TRUE
    # gene_list        = NULL   # NULL = core clock genes; or c("Per2", "Nr1d1", ...)
    # cell_group_rules = NULL   # NULL = built-in bins;    or list("T cells" = c(...), ...)
  ),
  error = function(e) message("  Clock acrophase plot skipped: ", e$message)
)

# 6d. ── CUSTOM SINGLE-GENE EXPLORER ─────────────────────────────────────────
# Run this block any time after Section 5 to inspect a specific gene.
# Tweak the three sections below (gene, display, appearance) then run the block.
#
#   Useful browsing commands (run in console):
#     View(T1[conf_mask, ])                                             # all confident
#     T1[T1$Genes == "Per2", ]                                         # one gene row
#     head(T1[conf_mask,][order(-T1[conf_mask,]$Abs_Amp),     ], 20)  # top amplitude
#     head(T1[conf_mask,][order( T1[conf_mask,]$Acrophase_24),], 20)  # by peak time

# ── GENE & DATA ───────────────────────────────────────────────────────────────
explore_gene   <- "Per2"   # gene symbol to plot
explore_violin <- TRUE    # FALSE = jitter dots  |  TRUE = violin density
explore_sc     <- TRUE     # TRUE  = overlay single-cell data
explore_save   <- TRUE     # TRUE  = save PNG to ct_dir

# ── APPEARANCE ────────────────────────────────────────────────────────────────
title_size     <- 13       # plot title font size (pts)
title_bold     <- TRUE     # TRUE = bold title
title_hjust    <- 0.5      # 0 = left  |  0.5 = centred  |  1 = right
axis_title_sz  <- 12       # x / y axis label size
axis_text_sz   <- 11       # x / y tick label size
legend_text_sz <- 10       # legend entry size
legend_pos     <- "top"    # "top" | "bottom" | "right" | "none"
y_min          <- NULL     # NULL = auto  |  e.g. 0   to force y-axis floor
y_max          <- NULL     # NULL = auto  |  e.g. 2.5 to cap y-axis ceiling
fig_width      <- 9        # saved PNG width  (inches)
fig_height     <- 5.5      # saved PNG height (inches)
fig_dpi        <- 180      # saved PNG resolution

# ─────────────────────────────────────────────────────────────────────────────
if (!explore_gene %in% T1$Genes) {
  cat(sprintf("  '%s' not found in T1 — check spelling.\n", explore_gene))
} else {
  gr <- T1[T1$Genes == explore_gene, ]
  cat(sprintf(
    "\n[%s]  Amp=%.3f  Mesor=%.3f  Acrophase=%.1fh  p=%.4g  padj=%.4g  %s\n",
    explore_gene,
    gr$Abs_Amp, gr$Mesor, gr$Acrophase_24,
    gr$pvalue, gr$pvalue_corr,
    ifelse(explore_gene %in% conf_genes, "** CONFIDENT **", "(not significant)")
  ))

  # NOTE: outdir = dir_gene is required — the function reads its input CSVs from
  # this folder. Passing NULL would break file.path() → "argument is of length zero".
  explore_plot <- tryCatch(
    plot_gene_single(
      tmeta        = tmeta,
      cust_cells   = key,
      period12     = period12,
      cust_gene    = explore_gene,
      print_scdata = explore_sc,
      sce          = data,
      celltype_col = celltype_col,
      zt_col       = zt_col,
      use_violin   = explore_violin,
      outdir       = ct_dir         # reads CSVs from here; saving is done below
    ),
    error = function(e) {
      message("  Plot error: ", conditionMessage(e))
      message("  Tip: make sure Section 5 has been run and ct_dir exists.")
      NULL
    }
  )

  if (!is.null(explore_plot)) {
    # Apply appearance overrides on top of the returned ggplot object
    explore_plot <- explore_plot +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(
                            size  = title_size,
                            face  = if (title_bold) "bold" else "plain",
                            hjust = title_hjust),
        axis.title      = ggplot2::element_text(size = axis_title_sz),
        axis.text       = ggplot2::element_text(size = axis_text_sz),
        legend.text     = ggplot2::element_text(size = legend_text_sz),
        legend.position = legend_pos
      )

    # Y-axis range: coord_cartesian zooms without removing data (safe for
    # violins whose tails may extend beyond the visible window).
    # Only applied when at least one of y_min / y_max is set.
    if (!is.null(y_min) || !is.null(y_max)) {
      explore_plot <- explore_plot +
        ggplot2::coord_cartesian(
          ylim = c(
            if (!is.null(y_min)) y_min else NA_real_,
            if (!is.null(y_max)) y_max else NA_real_
          )
        )
    }

    # Force render to the RStudio plot pane
    suppressWarnings(print(explore_plot))

    if (explore_save) {
      out_file <- file.path(dir_gene,
        sprintf("%s_%s_explorer.png", focus_safe, explore_gene))
      ggplot2::ggsave(out_file, explore_plot,
                      width = fig_width, height = fig_height,
                      dpi = fig_dpi, bg = "white")
      cat(sprintf("  Saved -> %s\n", out_file))
    }
  } else {
    cat("  Plot could not be generated — see error above.\n")
  }
}

# =============================================================================
# -- Step 7: Pull gene sets (KEGG + Reactome + GO:BP) --
# =============================================================================
# Using the local pull_genesets function (msigdbr wrapper)
cat("\n-- Step 7: Pulling MSigDB gene sets --\n")

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
  min_size    = 30L,
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

# Combine all lists into one background set for clusterProfiler
genesets <- c(gs_kegg, gs_reactome, gs_gobp)
cat(sprintf("  Combined : %d gene sets total\n\n", length(genesets)))


# =============================================================================
# -- Step 8: Phase-bin enrichment (Local clusterProfiler) --
# =============================================================================
# KEY IDEA: Instead of scoring whole pathways, we:
#   1. Bin confident circadian genes into narrow acrophase windows (e.g., ZT2-4).
#   2. Run ORA per bin using clusterProfiler to find time-specific biology.
#   3. Build PHASE-RESTRICTED gene sets: pathway genes INTERSECT bin genes.
#      This ensures all genes in the set co-peak, avoiding AUCell cancellation.

cat(sprintf("\n-- Step 8: Phase-bin enrichment for '%s' --\n", focus_ct))

phase_results <- phase_bin_analysis(
  T1               = T1,
  conf_mask        = conf_mask,
  genesets         = genesets,
  bin_width        = 3,      # Recommended 2hr bins for robust ORA
  n_top            = 20L,     # Top enriched pathways to keep per bin
  min_overlap      = 3L,     # Minimum genes in overlap to count
  min_bin_genes    = 3L,     # Minimum genes in a bin to attempt ORA
  p_thresh         = 0.05,
  use_padj         = FALSE,   # Using BH-correction (standard for clusterProfiler)
  exclude_patterns = c( "VIRAL", "INFECTION") 
)

# Print enrichment summary to console
cat("\nTop Enriched pathways per phase bin (via clusterProfiler):\n")
if (length(phase_results$ora_results) == 0) {
  cat("  [!] No bins reached significance. Try increasing bin_width or lowering min_overlap.\n")
} else {
  for (bin in names(phase_results$ora_results)) {
    df <- phase_results$ora_results[[bin]]
    for (i in seq_len(nrow(df))) {
      cat(sprintf("  %s | p=%.3g | Rich=%.2f | %s\n", 
                  bin, df$pvalue[i], df$RichFactor[i], df$Pathway[i]))
    }
  }
}

# ── Save ORA results to Excel ────────────────────────────────────────────────
ora_xlsx <- file.path(dir_pwA, paste0(focus_safe, "_phase_bin_enrich_ORA.xlsx"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  
  # Styling
  hdr_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#2F4F8F",
    halign = "center", textDecoration = "bold")
  sig_style <- openxlsx::createStyle(fgFill = "#E8F0FE") # Light blue for hits
  
  # 1. Summary Sheet (All Bins combined)
  if (length(phase_results$ora_results) > 0) {
    all_rows <- do.call(rbind, phase_results$ora_results)
    
    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", all_rows, rowNames = FALSE)
    openxlsx::setColWidths(wb, "Summary", cols = seq_len(ncol(all_rows)), widths = "auto")
    openxlsx::addStyle(wb, "Summary", hdr_style, rows = 1, cols = seq_len(ncol(all_rows)))
    
    # 2. Per-Bin Individual Sheets
    for (bin in names(phase_results$ora_results)) {
      # Clean sheet name (Excel doesn't like ZT00-02 as a name sometimes)
      sname <- gsub("[^[:alnum:]_]", "_", bin)
      df    <- phase_results$ora_results[[bin]]
      
      openxlsx::addWorksheet(wb, sname)
      openxlsx::writeData(wb, sname, df, rowNames = FALSE)
      openxlsx::setColWidths(wb, sname, cols = seq_len(ncol(df)), widths = "auto")
      openxlsx::addStyle(wb, sname, hdr_style, rows = 1, cols = seq_len(ncol(df)))
      
      # Highlight rows with high RichFactor
      high_rich <- which(df$RichFactor > 0.5)
      if (length(high_rich) > 0) {
        openxlsx::addStyle(wb, sname, sig_style, rows = high_rich + 1, 
                           cols = seq_len(ncol(df)), gridExpand = TRUE)
      }
    }
    
    openxlsx::saveWorkbook(wb, ora_xlsx, overwrite = TRUE)
    cat(sprintf("\n  [Success] ORA results saved -> %s\n", ora_xlsx))
  } else {
    cat("\n  [Skip] No Excel created (zero enriched bins).\n")
  }
}

# -- 9. AUCELL SCORING ON PHASE-RESTRICTED GENE SETS --------------------------
# Each gene set contains only co-phased genes (pathway INTERSECT acrophase bin).
# All members oscillate in the same window -> AUCell score oscillates cleanly.
# min_gs_size = 3 allows small sets from tight 1-hr bins.

auc_cache_phase <- file.path(dir_pwA,
  "auc_matrix_phase_KEGG_Reactome_GOBP.rds")

if (file.exists(auc_cache_phase)) {
  cat("\nLoading cached phase-restricted AUC matrix...\n")
  auc_phase <- readRDS(auc_cache_phase)
  cat(sprintf("  %d phase gene sets x %d cells\n",
              nrow(auc_phase), ncol(auc_phase)))
} else {
  cat("\nScoring cells with AUCell on phase-restricted gene sets...\n")
  auc_phase <- auc_score_cells(
    obj          = data,
    genesets     = phase_results$phase_gs,
    use_norm     = TRUE,
    auc_max_rank = 0.05,
    n_cores      = 1L,
    min_gs_size  = 3L
  )
  saveRDS(auc_phase, auc_cache_phase)
  cat(sprintf("  Cached -> %s\n", auc_cache_phase))
}

# -- 9b. METHOD B: AUCell on FULL pathway gene sets (no phase restriction) -----
# Method A (above) uses phase-restricted sets: pathway genes that co-peak in one
# ZT bin.  Method B scores the *entire* pathway gene set and asks whether the
# resulting activity score itself oscillates circadianly.
#
# Difference in interpretation:
#   Method A -> which WINDOW of the pathway oscillates, and does it oscillate?
#   Method B -> does the WHOLE pathway activity oscillate? (cancellation risk if
#               genes in the set peak at opposite phases, but useful as a check)
#
# `genesets` is already built in Section 7 from pull_genesets().

auc_cache_full <- file.path(dir_pwB,
  "auc_matrix_full_KEGG_Reactome_GOBP.rds")

if (file.exists(auc_cache_full)) {
  cat("\nLoading cached full-pathway AUC matrix (Method B)...\n")
  auc_full <- readRDS(auc_cache_full)
  cat(sprintf("  %d pathways x %d cells\n", nrow(auc_full), ncol(auc_full)))
} else {
  cat("\nScoring cells with AUCell on full pathway gene sets (Method B)...\n")
  auc_full <- auc_score_cells(
    obj          = data,
    genesets     = genesets,      # full sets from Section 7
    use_norm     = TRUE,
    auc_max_rank = 0.05,
    n_cores      = 1L,
    min_gs_size  = 10L            # larger minimum — full sets should be bigger
  )
  saveRDS(auc_full, auc_cache_full)
  cat(sprintf("  Cached -> %s\n", auc_cache_full))
}

# Circadian cosinor on full pathway scores
cat(sprintf("\nRunning cosinor on full pathway scores (Method B) for '%s'...\n",
            focus_ct))
path_full <- pathway_cosinor(
  auc_mat      = auc_full,
  meta         = data@meta.data,
  celltype_col = celltype_col,
  zt_col       = zt_col,
  tmeta        = tmeta,
  target_ct    = focus_ct,
  period12     = period12,
  custom_zt    = custom_zt
)

conf_full <- path_full$stats[
  path_full$stats$pvalue < 0.05 & path_full$stats$pvalue_corr < 0.05, ]
cat(sprintf("  Confident full-pathway sets: %d / %d\n",
            nrow(conf_full), nrow(path_full$stats)))
cat("\nTop 10 (Method B):\n")
print(head(conf_full[, c("Pathway","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))

# Write Excel
xlsx_full <- file.path(dir_pwB,
  paste0(focus_safe, "_full_pathway_circadian.xlsx"))
write_pathway_results(path_full, xlsx_full, celltype = focus_ct)
cat(sprintf("  Results saved -> %s\n", xlsx_full))

# -- 10. CIRCADIAN COSINOR ON PHASE-RESTRICTED SCORES -------------------------
# Each row of auc_phase is "ZTa-b__PATHWAY" — the cosinor should recover the
# expected acrophase matching the bin label if the biology is real.

cat(sprintf("\nRunning cosinor on phase-restricted pathway scores for '%s'...\n",
            focus_ct))
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
cat("\nTop 10:\n")
print(head(conf_phase[, c("Pathway","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))

# Write Excel (All + Confident sheets)
xlsx_phase <- file.path(dir_pwA,
  paste0(focus_safe, "_phase_pathway_circadian.xlsx"))
write_pathway_results(path_phase, xlsx_phase, celltype = focus_ct)
cat(sprintf("  Results saved -> %s\n", xlsx_phase))

# 10a. Plot top 6 confident phase-pathway sets
top_n_paths   <- 6L

if (nrow(conf_phase) > 0) {
  top_pw   <- head(conf_phase$Pathway, top_n_paths)
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
    print(gridExtra::grid.arrange(grid_pw))
    ggplot2::ggsave(
      filename = file.path(dir_pwA,
                           paste0(focus_safe, "_top_phase_pathways_grid.png")),
      plot     = grid_pw,
      width    = 6*n_col, height = 5*n_row, dpi = 150, bg = "white"
    )
  }
}

# 10b. Batch-save all confident phase-pathway plots
cat("\nBatch saving phase-pathway plots...\n")
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
  outdir       = ct_dir
)

# -- 11. GRN TIME SERIES ------------------------------------------------------
# Approach:
#   1. Build gene pool: confident circadian genes + core clock genes in data
#   2. Compute N×N Pearson correlation on z-scored expression (ALL cells pooled)
#      -> pooling is correct: genes are confirmed oscillators, their global
#         correlation captures phase relationships + shared regulators
#   3. Find hub genes = top hub_pct% by degree (number of significant edges)
#   4a. Approach-A GRN: hub genes ∩ ORA-enriched oscillating pathway genes
#   4b. Approach-B GRN: hub genes ∩ any oscillating pathway genes (conf_full)
#   -> 5-10 focused nodes per network, time-varying edges via plot_grn_timeseries

# ── Parameters ────────────────────────────────────────────────────────────────
grn_cor_thresh     <- 0.30   # |r| threshold for an edge to count toward degree
grn_p_thresh       <- 0.05   # p-value threshold for edges
grn_hub_pct        <- 0.10   # top X% of genes by degree = hub genes
grn_min_hub        <- 5L     # minimum hubs regardless of percentile
grn_min_nodes      <- 3L     # minimum hub nodes in a pathway to build a GRN
grn_node_size      <- 4
grn_label_size     <- 4.5
grn_edge_width_max <- 3.0
grn_zt_title_size  <- 18
grn_zt_hjust       <- 0.5

# ── 1. Gene pool ──────────────────────────────────────────────────────────────
grn_pool <- unique(c(
  conf_genes,                                              # circadian identified
  clock_genes_ref[clock_genes_ref %in% rownames(data)]    # clock genes in data
))
cat(sprintf("\n  [GRN] Gene pool: %d circadian + clock genes\n", length(grn_pool)))

# ── 2-3. Hub selection via global co-expression degree ────────────────────────
cat("  [GRN] Computing co-expression hub genes (pooled across all ZT)...\n")
hub_result <- tryCatch(
  select_hub_genes(
    obj          = data,
    gene_pool    = grn_pool,
    target_ct    = focus_ct,
    celltype_col = celltype_col,
    use_norm     = TRUE,
    cor_thresh   = grn_cor_thresh,
    p_thresh     = grn_p_thresh,
    hub_pct      = grn_hub_pct,
    min_hub      = grn_min_hub
  ),
  error = function(e) { message("  [GRN] Hub selection failed: ", e$message); NULL }
)

if (!is.null(hub_result)) {
  hub_genes  <- hub_result$hub_genes
  hub_circ   <- hub_genes[hub_genes %in% conf_genes]
  hub_clock  <- hub_genes[hub_genes %in% clock_genes_ref]
  hub_unique <- setdiff(hub_genes, c(conf_genes, clock_genes_ref))

  cat(sprintf("  [GRN] Hub genes: %d total\n", length(hub_genes)))
  if (length(hub_circ)   > 0) cat(sprintf("    Circadian hubs : %s\n", paste(hub_circ,   collapse = ", ")))
  if (length(hub_clock)  > 0) cat(sprintf("    Clock hubs     : %s\n", paste(hub_clock,  collapse = ", ")))
  if (length(hub_unique) > 0) cat(sprintf("    Other hubs     : %s\n", paste(hub_unique, collapse = ", ")))

  # Helper: run one GRN per pathway list
  .run_grn_batch <- function(osc_df, gs_list, grn_dir, label) {
    if (nrow(osc_df) == 0) {
      cat(sprintf("  [GRN %s] No oscillating pathways — skipped.\n", label)); return(invisible(NULL))
    }
    if (!dir.exists(grn_dir)) dir.create(grn_dir, recursive = TRUE)
    pws     <- osc_df$Pathway[order(osc_df$pvalue_corr)]
    n_built <- 0L
    for (pw in pws) {
      pw_genes  <- gs_list[[pw]]
      grn_nodes <- intersect(hub_genes, pw_genes)
      if (length(grn_nodes) < grn_min_nodes) next
      pw_safe     <- substr(gsub("[^[:alnum:]_]", "_", pw), 1, 60)
      grn_outfile <- file.path(grn_dir,
        sprintf("%s_GRN_%s_%s.png", focus_safe, label, pw_safe))
      cat(sprintf("  [GRN %s] %s  (%d hub nodes)\n", label, pw, length(grn_nodes)))
      tryCatch(
        plot_grn_timeseries(
          obj            = data,
          circ_genes     = grn_nodes,
          pathway_genes  = character(0),
          meta           = data@meta.data,
          celltype_col   = celltype_col,
          zt_col         = zt_col,
          tmeta          = tmeta,
          target_ct      = focus_ct,
          cor_thresh     = grn_cor_thresh,
          p_thresh       = grn_p_thresh,
          use_norm       = TRUE,
          outfile        = grn_outfile,
          ncol           = NULL,
          node_size      = grn_node_size,
          label_size     = grn_label_size,
          edge_width_max = grn_edge_width_max,
          zt_title_size  = grn_zt_title_size,
          zt_title_hjust = grn_zt_hjust
        ),
        error = function(e) message("    Skip GRN: ", e$message)
      )
      n_built <- n_built + 1L
    }
    cat(sprintf("  [GRN %s] %d GRN(s) saved -> %s\n", label, n_built, grn_dir))
  }

  grn_base <- file.path(ct_dir, "04_GRN")

  # ── 4a. Approach-A GRN: hub genes ∩ ORA-enriched oscillating pathways ────────
  # Build genesets from the enriched pathway list (Section 7 genesets)
  cat("\n  [GRN] Approach A (phase-restricted pathways)...\n")
  # conf_phase comes from Section 10 (phase-restricted cosinor)
  .run_grn_batch(
    osc_df  = if (exists("conf_phase")) conf_phase else data.frame(),
    gs_list = phase_results$phase_gs,
    grn_dir = file.path(grn_base, "approach_A"),
    label   = "A"
  )

  # ── 4b. Approach-B GRN: hub genes ∩ full-pathway oscillating pathways ─────────
  cat("\n  [GRN] Approach B (full-pathway oscillating)...\n")
  # Build named list from genesets (same as Section 7)
  .run_grn_batch(
    osc_df  = if (exists("conf_full")) conf_full else data.frame(),
    gs_list = genesets,
    grn_dir = file.path(grn_base, "approach_B"),
    label   = "B"
  )
}


# =============================================================================
# POST-PROCESSING — Custom acrophase polar plots (run any time after Section 5)
# =============================================================================
# load_stage_results() reads all circadian_analysis_all.csv files from disk so
# you can re-draw the polar plot without re-running the analysis.
# Useful here for exploring gene sets with more than one cell type loaded.
#
# ── Step 1: load all cell types from the output directory ─────────────────────
# clock_results_post <- load_stage_results(out_dir, period12 = period12)
#
# ── Step 2: inspect which confident genes appear across cell types ─────────────
# all_genes <- sort(unique(unlist(lapply(clock_results_post,
#   function(df) df$Genes[df$pvalue < 0.05 & df$pvalue_corr < 0.05]))))
# head(all_genes, 30)
#
# ── Step 3: plot with custom gene list and cell-type grouping ─────────────────
# plot_clock_acrophase(
#   results_list     = clock_results_post,
#   stage            = "test",
#   outfile          = file.path(out_dir, "clock_acrophase_custom.png"),
#   gene_list        = c("Per2", "Nr1d1", "Cry1", "Dbp", "Arntl"),
#   cell_group_rules = list(
#     "T cells"  = c("CD8+ T cells", "CD4+ T cells", "T-Reg cells"),
#     "Myeloid"  = c("M1 macrophages", "M2 macrophages", "Monocytes", "DCs"),
#     "Stromal"  = c("Endothelial", "Fibroblasts"),
#     "Tumor"    = c("Tumor cells")
#   ),
#   strict = TRUE,
#   dpi    = 300
# )