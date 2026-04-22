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

base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm\\r_pre_process\\TimeSCape_R_tes"
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
sc_style   <- "dots"    # "dots" = jitter  |  "violin" = violin density

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
      filename = file.path(ct_dir, paste0(key, "_top", top_n_show, "_genes.png")),
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
explore_violin <- FALSE    # FALSE = jitter dots  |  TRUE = violin density
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

  # NOTE: outdir = ct_dir is required — the function reads its input CSVs from
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
      out_file <- file.path(ct_dir,
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

# -- 7. PULL GENE SETS (KEGG + Reactome + GO:BP) ------------------------------
# We pull THREE collections so that phyper and EnrichR test the SAME pathway
# universe — this is what allows consensus hits.
#
# phyper uses these local msigdbr sets (mouse gene symbols, custom background).
# EnrichR queries the matching databases on the server (human annotation).
# Consensus = significant in both → highest confidence hits.
#
#   KEGG     CP:KEGG_LEGACY  ~186  curated metabolic / signalling pathways
#   Reactome CP:REACTOME     ~1839 fine-grained pathway hierarchy (capped at 300 genes)
#   GO:BP    GO:BP           ~7538 biological process terms  (dedup removes redundancy)
#
# After size + dedup filters expect ~2000–4000 total sets.
# Total phyper runtime is still fast (pure R, no network).
#
# Required:  install.packages(c("msigdbr", "openxlsx"))
#            BiocManager::install("AUCell")

cat("\n-- Step 7: Pull gene sets (KEGG + Reactome + GO:BP) --\n")

gs_kegg <- pull_genesets(
  collection  = "C2",
  subcategory = "CP:KEGG_LEGACY",
  species     = "Mus musculus",
  min_size    = 10L,
  max_size    = 500L,        # KEGG is non-redundant; keep full size range
  deduplicate = FALSE
)
cat(sprintf("  KEGG     : %d gene sets\n", length(gs_kegg)))

gs_reactome <- pull_genesets(
  collection  = "C2",
  subcategory = "CP:REACTOME",
  species     = "Mus musculus",
  min_size    = 30L,         # 30 removes virus/disease-specific sets (ribosomal noise)
                             # while keeping genuine immune, metabolic, and RNA pathways
  max_size    = 500L,        # immune pathways are legitimately large
  deduplicate = FALSE        # Reactome hierarchy is meaningful; keep sub-pathways
)
cat(sprintf("  Reactome : %d gene sets\n", length(gs_reactome)))

# GO:BP is very large (~7500 terms).  For CD8+ T cells most terms are irrelevant
# (liver, neuronal, etc.), so we filter more aggressively:
#   min_size = 20  — removes the most specific low-power terms (< 20 genes)
#   max_size = 200 — removes terms so broad they match almost everything
#   deduplicate    — drops POSITIVE/NEGATIVE_REGULATION_OF_X when X exists
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

# -- 8. PHASE-BIN ENRICHMENT --------------------------------------------------
# KEY IDEA: instead of scoring whole pathways, we:
#   1. Bin confident circadian genes into narrow acrophase windows
#   2. Run ORA per bin — which pathways are enriched in genes peaking at THIS phase?
#      TWO methods run in parallel, now testing the SAME databases:
#        a. phyper (hypergeometric): local msigdbr sets, custom tested-gene background
#        b. EnrichR API: same collections (KEGG + Reactome + GO:BP), human annotation
#      Consensus = significant in both → highest confidence.
#      phyper_only / enrichr_only = found in one method only.
#   3. Build PHASE-RESTRICTED gene sets = pathway genes INTERSECT bin genes
#      -> every member co-peaks in the same window -> no cancellation in AUCell
#
# bin_width = 3 hr  -> one ZT interval per bin (8 bins for 8-ZT designs)
# bin_width = 1 hr  -> tighter co-regulation windows (more bins, fewer genes each)

cat(sprintf("\n-- Step 8: Phase-bin enrichment for '%s' --\n", focus_ct))
phase_results <- phase_bin_analysis(
  T1            = T1,
  conf_mask     = conf_mask,
  genesets      = genesets,
  bin_width     = 3,      # hours; 3 = one ZT interval, 1 = tightest co-regulation
  n_top         = 5L,     # top enriched pathways to keep per bin (consensus ranked first)
  min_overlap   = 3L,     # min gene overlap for phyper ORA
  min_bin_genes = 3L,     # min genes in a bin to attempt ORA
  p_thresh      = 0.05,   # p-value threshold (applied to both phyper and EnrichR)
  use_padj      = FALSE,  # FALSE = raw p (recommended: bins too small for BH correction)
                          # TRUE  = BH-adjusted (only reliable with many conf. genes)
  use_enrichr   = F,   # TRUE  = query EnrichR API for consensus validation
                          # FALSE = phyper only (use when offline)
  enrichr_dbs   = c(      # must match the collections pulled above
    "KEGG_2026",
    "Reactome_Pathways_2024",
    "GO_Biological_Process_2025"
  )
  # exclude_patterns = NULL  (default — keep all ORA hits)
  # Irrelevant pathways (viral, cardiac, neuronal) are filtered naturally by
  # the cosinor in section 10: a pathway only survives if its phase-restricted
  # AUCell score actually oscillates rhythmically in these cells.
)

# Print enrichment summary (source column: consensus / phyper_only / enrichr_only)
cat("\nEnriched pathways per phase bin:\n")
for (bin in names(phase_results$ora_results)) {
  df  <- phase_results$ora_results[[bin]]
  for (i in seq_len(nrow(df))) {
    cat(sprintf("  %s [%s]  %s\n", bin, df$source[i], df$Pathway[i]))
  }
}

# Save ORA results to Excel (one sheet per bin + a Summary sheet)
# source column: "consensus" = sig in both phyper AND EnrichR  (highest confidence)
#                "phyper_only" = sig in phyper only (correct custom background)
#                "enrichr_only" = sig in EnrichR only (genome-wide background)
ora_xlsx <- file.path(out_dir, paste0(focus_safe, "_phase_bin_ORA.xlsx"))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()

  # Color styles: consensus = green, phyper_only = yellow, enrichr_only = blue
  style_consensus  <- openxlsx::createStyle(fgFill = "#C6EFCE")
  style_phyper     <- openxlsx::createStyle(fgFill = "#FFEB9C")
  style_enrichr    <- openxlsx::createStyle(fgFill = "#BDD7EE")
  hdr_style        <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#2F4F8F",
    halign = "center", textDecoration = "bold")

  # Summary sheet
  all_rows <- do.call(rbind, lapply(names(phase_results$ora_results), function(b) {
    df <- phase_results$ora_results[[b]]
    cbind(Bin = b, df)
  }))
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", all_rows, rowNames = FALSE)
  openxlsx::setColWidths(wb, "Summary", cols = seq_len(ncol(all_rows)), widths = "auto")
  openxlsx::addStyle(wb, "Summary", hdr_style, rows = 1L,
                     cols = seq_len(ncol(all_rows)), gridExpand = TRUE)
  # Color summary rows by source
  for (ri in seq_len(nrow(all_rows))) {
    sty <- switch(all_rows$source[ri],
                  "consensus"    = style_consensus,
                  "phyper_only"  = style_phyper,
                  "enrichr_only" = style_enrichr,
                  NULL)
    if (!is.null(sty))
      openxlsx::addStyle(wb, "Summary", sty, rows = ri + 1L,
                         cols = seq_len(ncol(all_rows)),
                         gridExpand = TRUE, stack = TRUE)
  }

  # Per-bin sheets
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
                    "enrichr_only" = style_enrichr,
                    NULL)
      if (!is.null(sty))
        openxlsx::addStyle(wb, sname, sty, rows = ri + 1L,
                           cols = seq_len(ncol(df)),
                           gridExpand = TRUE, stack = TRUE)
    }
  }

  openxlsx::saveWorkbook(wb, ora_xlsx, overwrite = TRUE)
  cat(sprintf("  ORA results saved -> %s\n", ora_xlsx))
  cat("  Color key: GREEN=consensus  YELLOW=phyper-only  BLUE=EnrichR-only\n")
}

# -- 9. AUCELL SCORING ON PHASE-RESTRICTED GENE SETS --------------------------
# Each gene set contains only co-phased genes (pathway INTERSECT acrophase bin).
# All members oscillate in the same window -> AUCell score oscillates cleanly.
# min_gs_size = 3 allows small sets from tight 1-hr bins.

auc_cache_phase <- file.path(out_dir,
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
xlsx_phase <- file.path(out_dir,
  paste0(focus_safe, "_phase_pathway_circadian.xlsx"))
write_pathway_results(path_phase, xlsx_phase, celltype = focus_ct)
cat(sprintf("  Results saved -> %s\n", xlsx_phase))

# 10a. Plot top 6 confident phase-pathway sets
top_n_paths   <- 6L
path_plot_dir <- file.path(out_dir, paste0(focus_safe, "_phase_pathway_plots"))
if (!dir.exists(path_plot_dir)) dir.create(path_plot_dir, recursive = TRUE)

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
      filename = file.path(path_plot_dir,
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
  outdir       = path_plot_dir
)

# -- 11. GRN TIME SERIES ------------------------------------------------------
# Gene selection:
#   circ_genes    = top 1-3 confident genes globally by (Abs_Amp + Mesor)
#                   -> strongest oscillators with highest baseline expression
#                 + any core clock genes present in the confident list
#                   -> always anchor the network with clock machinery if detected
#   pathway_genes = all genes from confident phase-restricted sets (all bins)
#                   -> full pathway biology across the 24-hr cycle
#
# Nodes:
#   orange-red = top circadian / clock gene
#   blue       = phase-restricted pathway gene
#   purple     = both
#
# Required: install.packages(c("igraph","ggraph","cowplot"))

# -- Core clock gene list (mouse symbols) ------------------------------------
clock_genes_ref <- c(
  "Arntl", "Bmal1", "Clock", "Npas2",        # positive arm (BMAL1/CLOCK)
  "Per1",  "Per2",  "Per3",                    # period genes
  "Cry1",  "Cry2",                             # cryptochrome genes
  "Nr1d1", "Nr1d2",                            # REV-ERB alpha/beta
  "Rora",  "Rorb",  "Rorc",                   # ROR alpha/beta/gamma
  "Dbp",   "Tef",   "Hlf",   "Nfil3",        # PAR/D-box TFs
  "Ciart", "Bhlhe40", "Bhlhe41",              # DEC1/DEC2
  "Timeless", "Csnk1d", "Csnk1e"             # CK1 delta/epsilon
)

# -- Top N per phase bin by (Abs_Amp + Mesor) ---------------------------------
# Each active bin contributes its 1-3 strongest genes so every ZT window
# has a representative in the network.
top_per_bin   <- 3L      # 1-3 recommended
T1_conf_all   <- T1[conf_mask, , drop = FALSE]
T1_conf_grn   <- phase_results$bin_table   # same rows + phase_bin column

bin_rep_genes <- character(0)
for (bin in names(phase_results$ora_results)) {
  bin_mask  <- as.character(T1_conf_grn$phase_bin) == bin
  bin_genes <- T1_conf_grn$Genes[bin_mask]
  bin_stats <- T1_conf_all[T1_conf_all$Genes %in% bin_genes, , drop = FALSE]
  if (nrow(bin_stats) == 0) next
  score  <- bin_stats$Abs_Amp + bin_stats$Mesor
  top_g  <- bin_stats$Genes[order(-score)]
  top_g  <- head(top_g, top_per_bin)
  bin_rep_genes <- c(bin_rep_genes, top_g)
  cat(sprintf("  [%s] top genes: %s\n", bin, paste(top_g, collapse = ", ")))
}
top_circ <- unique(bin_rep_genes)
cat(sprintf("\n  Phase-bin representative genes: %d\n", length(top_circ)))

# -- Clock genes present in the confident circadian set -----------------------
clock_in_data <- conf_genes[conf_genes %in% clock_genes_ref]
if (length(clock_in_data) > 0)
  cat(sprintf("  Clock genes detected: %s\n",
              paste(clock_in_data, collapse = ", ")))

# Combined circ_genes: top oscillators + any clock genes
circ_genes_grn <- unique(c(top_circ, clock_in_data))
cat(sprintf("  Total circ_genes for GRN: %d\n", length(circ_genes_grn)))

# -- One GRN per confident phase-restricted pathway ---------------------------
# Each GRN contains:
#   circ_genes    = top-3-per-bin representatives + clock genes (already built)
#   pathway_genes = genes of THIS specific phase-restricted set
# One network per pathway keeps the signal clean and comparable.

# ── GRN APPEARANCE ────────────────────────────────────────────────────────────
grn_node_size      <- 4       # node dot size
grn_label_size     <- 4.5    # gene label font size (ggplot pts)
grn_edge_width_max <- 3.0    # max edge line width (scaled by |r|)
grn_zt_title_size  <- 18     # ZT panel title size (pts)
grn_zt_hjust       <- 0.5    # 0 = left | 0.5 = centred | 1 = right
grn_cor_thresh     <- 0.2    # min |r| to draw an edge
grn_p_thresh       <- 0.05   # max p-value to draw an edge
# ─────────────────────────────────────────────────────────────────────────────

if (nrow(conf_phase) > 0) {
  conf_gs_names <- conf_phase$Pathway
  grn_keys <- names(phase_results$phase_gs)[
    sapply(names(phase_results$phase_gs), function(k)
      any(sapply(conf_gs_names, function(p) grepl(p, k, fixed = TRUE))))
  ]
} else {
  grn_keys <- names(phase_results$phase_gs)
}

if (length(grn_keys) == 0) {
  cat("  No confident phase-restricted sets for GRN.\n")
} else {
  cat(sprintf("\nBuilding GRN time series -- %d pathway(s)...\n", length(grn_keys)))
  grn_dir <- file.path(ct_dir, paste0(focus_safe, "_GRN"))
  if (!dir.exists(grn_dir)) dir.create(grn_dir, recursive = TRUE)

  for (pw_key in grn_keys) {
    pw_genes    <- phase_results$phase_gs[[pw_key]]
    pw_safe     <- gsub("[^[:alnum:]_]", "_", pw_key)
    grn_outfile <- file.path(grn_dir,
      sprintf("%s_GRN_%s.png", focus_safe, pw_safe))

    cat(sprintf("  [%s]  %d pathway genes\n", pw_key, length(pw_genes)))

    grn_plot <- tryCatch(
      plot_grn_timeseries(
        obj            = data,
        circ_genes     = circ_genes_grn,
        pathway_genes  = pw_genes,
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
      error = function(e) {
        message("  Skip GRN for ", pw_key, ": ", e$message)
        NULL
      }
    )
  }
  cat(sprintf("  GRN plots saved -> %s\n", grn_dir))
}
