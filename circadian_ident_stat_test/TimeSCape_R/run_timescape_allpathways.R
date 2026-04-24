# ============================================================================
# TimeSCape v0.2 -- Approach B: Score ALL pathways, detect oscillating programs
#
# Logic:
#   1. Load Seurat object (no gene-level cosinor pre-filter)
#   2. Build TERM2GENE from msigdbr (KEGG + Reactome + GO:BP + GO:MF)
#   3. Score ALL cells for ALL pathways (AUCell or AddModuleScore)
#      -- no enrichment filter, no gene selection step
#   4. Run cosinor on pathway scores -> find which gene programs oscillate
#
# Rationale:
#   - Genes that are not individually circadian can collectively form a
#     circadian gene program (phase-coherent expression across cells).
#   - By scoring all pathway programs and running cosinor on the scores,
#     we detect emergent rhythmic programs missed by gene-level filtering.
#   - Complement to Approach A: A detects pathways enriched IN circadian
#     genes; B detects pathways whose COLLECTIVE expression is rhythmic.
#
# Key differences from Approach A (run_timescape_enricher.R):
#   - No gene-level cosinor step (Steps 5 is optional, just for reference)
#   - No enricher() ORA -- all pathways are scored unconditionally
#   - Much larger score matrix (thousands of pathways vs. tens)
#   - Cosinor FDR correction is the sole filter for pathway rhythm
#
# Required packages:
#   install.packages("msigdbr")
#   BiocManager::install("AUCell")  # only if use_aucell = TRUE
# ============================================================================

library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)
library(dplyr)

set.seed(123)

# =============================================================================
# 1. PATHS
# =============================================================================

base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm\\r_pre_process\\TimeSCape_R_tes"
src_path  <- "C:\\Users\\selim\\Documentos\\vscode_working_dir\\single_cell_projects\\circadian_ident_stat_test\\TimeSCape_R\\R"
out_dir   <- file.path(base_dir, "TimeSCape_output_allpathways")

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
  "circadian_mouse_breast_decontX_QC_2_23_2026_subann_early.rds"))
cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(data), nrow(data)))

# =============================================================================
# 4. PARAMETERS
# =============================================================================
# Define all variants to collapse
cd8_variants <- c("CD8+ T cells", "Cyc. CD8+ T cells")

# Create new column, copy existing labels, overwrite CD8 variants
data$cell_types_merged <- as.character(data$sub_cell_types2)
data$cell_types_merged[data$sub_cell_types2 %in% cd8_variants] <- "CD8+ T cells"

# Verify
table(data$cell_types_merged)[grep("CD8|T cell", names(table(data$cell_types_merged)))]

celltype_col <- "broad_cell_types"
zt_col       <- "ZT_time_str"
period12     <- FALSE
norm_str     <- "logcounts"
n_workers    <- 2L
focus_ct     <- "CD8+ T cells"
focus_ct     <- "Endothelial"
focus_ct     <- "Tumor"
focus_safe   <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))

# ── SCORING METHOD ────────────────────────────────────────────────────────────
# TRUE  = AUCell  (rank-based, robust, slower; needs BiocManager::install("AUCell"))
# FALSE = AddModuleScore (Seurat built-in, fast, good for prototyping)
use_aucell <- TRUE

# ── PATHWAY DATABASE SELECTION ────────────────────────────────────────────────
# Comment/uncomment collections to include. All four is comprehensive but slow.
# Reactome alone is a good balance of coverage and speed.
use_kegg     <- TRUE    # ~200 pathways
use_reactome <- TRUE    # ~1,500 pathways
use_gobp     <- TRUE    # ~7,500 terms -- large; set FALSE for a quick run
use_gomf     <- FALSE   # ~1,000 terms -- often overlaps GO:BP, optional

# ── GENE SET SIZE FILTER ──────────────────────────────────────────────────────
# Applied before scoring to skip trivially small or genome-wide sets.
min_gs_size <- 10L
max_gs_size <- 500L

# ── COSINOR THRESHOLDS ────────────────────────────────────────────────────────
# Applied AFTER pathway scoring to select rhythmic programs.
cosinor_pval      <- 0.05
cosinor_pval_corr <- 0.05

tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)
ct_dir <- file.path(out_dir, focus_safe)
if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

# =============================================================================
# 5. (OPTIONAL) GENE-LEVEL COSINOR -- for reference only, not used as filter
# =============================================================================
# Uncomment if you want to also see which individual genes oscillate.
# The results are NOT used to filter pathways in Approach B.
#
# cat(sprintf("\n=== Step 5 [optional]: Cosinor reference -- '%s' ===\n", focus_ct))
# results_ref <- run_timescape(
#   obj             = data,
#   celltype_col    = celltype_col,
#   zt_col          = zt_col,
#   tmeta           = tmeta,
#   rm_low_conf     = TRUE,
#   period12        = period12,
#   custom_celltype = focus_ct,
#   plot_heat       = TRUE,
#   norm_str        = norm_str,
#   outdir          = out_dir,
#   n_workers       = n_workers
# )
# key_ref   <- names(results_ref)[1]
# T1_ref    <- results_ref[[key_ref]]$T1
# conf_mask <- T1_ref$pvalue < 0.05 & T1_ref$pvalue_corr < 0.05
# cat(sprintf("  Confident circadian genes (reference): %d / %d\n",
#             sum(conf_mask), nrow(T1_ref)))

# =============================================================================
# 6. BUILD TERM2GENE FROM msigdbr (mouse gene symbols)
# =============================================================================

cat("\n=== Step 6: Building TERM2GENE from msigdbr ===\n")

if (!requireNamespace("msigdbr", quietly = TRUE))
  stop("Install msigdbr:  install.packages('msigdbr')")

t2g_list <- list()

if (use_kegg) {
  cat("  Pulling KEGG_LEGACY...\n")
  t2g_list$kegg <- msigdbr::msigdbr(species = "Mus musculus",
                                      collection = "C2",
                                      subcollection = "CP:KEGG_LEGACY") |>
    dplyr::select(gs_name, gene_symbol)
}

if (use_reactome) {
  cat("  Pulling Reactome...\n")
  t2g_list$reactome <- msigdbr::msigdbr(species = "Mus musculus",
                                          collection = "C2",
                                          subcollection = "CP:REACTOME") |>
    dplyr::select(gs_name, gene_symbol)
}

if (use_gobp) {
  cat("  Pulling GO:BP...\n")
  t2g_list$gobp <- msigdbr::msigdbr(species = "Mus musculus",
                                      collection = "C5",
                                      subcollection = "GO:BP") |>
    dplyr::select(gs_name, gene_symbol)
}

if (use_gomf) {
  cat("  Pulling GO:MF...\n")
  t2g_list$gomf <- msigdbr::msigdbr(species = "Mus musculus",
                                      collection = "C5",
                                      subcollection = "GO:MF") |>
    dplyr::select(gs_name, gene_symbol)
}

TERM2GENE <- dplyr::bind_rows(t2g_list)
cat(sprintf("  Raw TERM2GENE: %d entries, %d unique pathways\n",
            nrow(TERM2GENE), length(unique(TERM2GENE$gs_name))))

# =============================================================================
# 7. BUILD GENE SET LIST (all pathways, size-filtered)
# =============================================================================

cat(sprintf("\n=== Step 7: Building gene set list (size %d-%d) ===\n",
            min_gs_size, max_gs_size))

# Count genes per pathway intersected with genes present in the object
genes_in_obj <- rownames(data)

gs_sizes <- TERM2GENE |>
  dplyr::filter(gene_symbol %in% genes_in_obj) |>
  dplyr::group_by(gs_name) |>
  dplyr::summarise(n_in_obj = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_in_obj >= min_gs_size, n_in_obj <= max_gs_size)

cat(sprintf("  Pathways passing size filter: %d / %d\n",
            nrow(gs_sizes), length(unique(TERM2GENE$gs_name))))

# Build named list: pathway -> gene symbols (intersection with object)
all_gs <- lapply(gs_sizes$gs_name, function(pw) {
  genes <- unique(TERM2GENE$gene_symbol[TERM2GENE$gs_name == pw])
  intersect(genes, genes_in_obj)
})
names(all_gs) <- gs_sizes$gs_name

cat(sprintf("  Total gene sets to score: %d\n", length(all_gs)))

if (length(all_gs) == 0)
  stop("No gene sets passed the size filter. Check min_gs_size/max_gs_size or database selection.")

# =============================================================================
# 8. SCORE CELLS (AUCell or AddModuleScore) -- ALL pathways
# =============================================================================

cat(sprintf("\n=== Step 8: Scoring cells (%s) for %d pathways ===\n",
            if (use_aucell) "AUCell" else "AddModuleScore",
            length(all_gs)))

# CRITICAL: subset to focus cell type before scoring.
# AddModuleScore samples control genes from a background distribution computed
# across ALL cells in the object -- mixing cell types biases the correction
# for any individual cell type. AUCell is per-cell and less sensitive, but
# subsetting also speeds it up significantly (fewer cells to rank).
focus_cells <- colnames(data)[data@meta.data[[celltype_col]] == focus_ct]
data_sub    <- data[, focus_cells]
cat(sprintf("  Subsetting to '%s': %d / %d cells\n",
            focus_ct, ncol(data_sub), ncol(data)))

score_cache <- file.path(ct_dir,
  paste0(focus_safe, "_allpathways_scores_",
         if (use_aucell) "aucell" else "ams", ".rds"))

if (file.exists(score_cache)) {
  cat("  Loading cached scores...\n")
  auc_mat <- readRDS(score_cache)
  cat(sprintf("  %d pathways x %d cells\n", nrow(auc_mat), ncol(auc_mat)))

} else if (use_aucell) {
  # ── AUCell ────────────────────────────────────────────────────────────────
  # Rank-based enrichment score per cell, robust to library size variation.
  # auc_score_cells() is defined in pathway_circadian.R
  auc_mat <- auc_score_cells(
    obj          = data_sub,
    genesets     = all_gs,
    use_norm     = TRUE,
    auc_max_rank = 0.05,
    n_cores      = 1L,
    min_gs_size  = min_gs_size
  )
  saveRDS(auc_mat, score_cache)
  cat(sprintf("  Cached -> %s\n", score_cache))

} else {
  # ── AddModuleScore (Seurat built-in) ──────────────────────────────────────
  # Scores in batches of 50 to keep memory manageable.
  # ctrl=5 (down from default 100) for speed; raise to 20-50 for publication.
  batch_size <- 50L
  pw_names   <- names(all_gs)
  n_pw       <- length(pw_names)
  score_rows <- vector("list", n_pw)

  cat(sprintf("  Running in batches of %d...\n", batch_size))

  for (start in seq(1, n_pw, by = batch_size)) {
    end      <- min(start + batch_size - 1L, n_pw)
    batch_nm <- pw_names[start:end]
    batch_gs <- all_gs[batch_nm]

    tmp <- Seurat::AddModuleScore(
      object   = data_sub,
      features = batch_gs,
      name     = "score",
      seed     = 42,
      ctrl     = 5        # lower for speed; use 20-50 for publication quality
    )
    idx_cols <- paste0("score", seq_along(batch_gs))
    for (j in seq_along(batch_nm)) {
      score_rows[[start + j - 1L]] <- tmp@meta.data[[idx_cols[j]]]
    }

    if ((end %% 500) == 0 || end == n_pw)
      cat(sprintf("  Scored %d / %d pathways...\n", end, n_pw))
  }

  # Pathways x cells matrix
  auc_mat <- do.call(rbind, score_rows)
  rownames(auc_mat) <- pw_names
  colnames(auc_mat) <- colnames(data_sub)
  saveRDS(auc_mat, score_cache)
  cat(sprintf("  Cached -> %s\n", score_cache))
}

# =============================================================================
# 9. PATHWAY COSINOR ON ALL SCORED PATHWAYS
# =============================================================================

cat(sprintf("\n=== Step 9: Pathway cosinor -- '%s' ===\n", focus_ct))
cat(sprintf("  Fitting cosinor for %d pathways...\n", nrow(auc_mat)))

path_res <- pathway_cosinor(
  auc_mat      = auc_mat,
  meta         = data_sub@meta.data,   # must match auc_mat columns (subsetted cells)
  celltype_col = celltype_col,
  zt_col       = zt_col,
  tmeta        = tmeta,
  target_ct    = focus_ct,
  period12     = period12
)

conf_phase <- path_res$stats[
  path_res$stats$pvalue      < cosinor_pval &
  path_res$stats$pvalue_corr < cosinor_pval_corr, ]

cat(sprintf("  Oscillating pathways: %d / %d\n",
            nrow(conf_phase), nrow(path_res$stats)))

if (nrow(conf_phase) > 0) {
  cat("\n  Top oscillating pathways:\n")
  print(head(conf_phase[order(conf_phase$pvalue_corr),
                         c("Pathway","Abs_Amp","Acrophase_24",
                           "pvalue","pvalue_corr")], 20))
}

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("\n=== Step 10: Saving results ===\n")

# Helper: attach gene members to a stats data frame from all_gs
add_gene_column <- function(df, gs_list) {
  df$Genes <- sapply(df$Pathway, function(pw) {
    g <- gs_list[[pw]]
    if (is.null(g) || length(g) == 0) NA_character_
    else paste(sort(g), collapse = "; ")
  })
  df
}

# Full cosinor stats for all pathways -- with gene lists
all_stats <- add_gene_column(path_res$stats, all_gs)
all_stats  <- all_stats[order(all_stats$pvalue_corr), ]

if (requireNamespace("openxlsx", quietly = TRUE)) {
  xlsx_all <- file.path(ct_dir, paste0(focus_safe, "_allpathways_cosinor_all.xlsx"))
  wb_all   <- openxlsx::createWorkbook()
  hdr      <- openxlsx::createStyle(fontColour="#FFFFFF", fgFill="#2F4F8F",
                                     halign="center", textDecoration="bold")
  wrap_st  <- openxlsx::createStyle(wrapText=TRUE, valign="top")
  openxlsx::addWorksheet(wb_all, "All_Pathways")
  openxlsx::writeData(wb_all, "All_Pathways", all_stats, rowNames=FALSE)
  openxlsx::addStyle(wb_all, "All_Pathways", hdr, rows=1,
                     cols=seq_len(ncol(all_stats)), gridExpand=TRUE)
  # Wrap the Genes column
  gene_col_idx <- which(names(all_stats) == "Genes")
  openxlsx::addStyle(wb_all, "All_Pathways", wrap_st,
                     rows=seq(2, nrow(all_stats)+1),
                     cols=gene_col_idx, gridExpand=TRUE)
  openxlsx::setColWidths(wb_all, "All_Pathways",
                         cols=seq_len(ncol(all_stats) - 1L), widths="auto")
  openxlsx::setColWidths(wb_all, "All_Pathways",
                         cols=gene_col_idx, widths=60)
  openxlsx::saveWorkbook(wb_all, xlsx_all, overwrite=TRUE)
  cat(sprintf("  All pathways -> %s\n", xlsx_all))
} else {
  write_pathway_results(path_res, xlsx_all, celltype = focus_ct)
  cat(sprintf("  All pathways -> %s\n", xlsx_all))
}

# Significant pathways only -- with gene lists
if (nrow(conf_phase) > 0 &&
    requireNamespace("openxlsx", quietly = TRUE)) {
  conf_with_genes <- add_gene_column(
    conf_phase[order(conf_phase$pvalue_corr), ], all_gs)

  xlsx_sig <- file.path(ct_dir, paste0(focus_safe, "_allpathways_cosinor_significant.xlsx"))
  wb <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(fontColour="#FFFFFF", fgFill="#2F4F8F",
                                halign="center", textDecoration="bold")
  wrap_st <- openxlsx::createStyle(wrapText=TRUE, valign="top")
  openxlsx::addWorksheet(wb, "Oscillating_Pathways")
  openxlsx::writeData(wb, "Oscillating_Pathways", conf_with_genes, rowNames=FALSE)
  openxlsx::addStyle(wb, "Oscillating_Pathways", hdr, rows=1,
                     cols=seq_len(ncol(conf_with_genes)), gridExpand=TRUE)
  gene_col_idx <- which(names(conf_with_genes) == "Genes")
  openxlsx::addStyle(wb, "Oscillating_Pathways", wrap_st,
                     rows=seq(2, nrow(conf_with_genes)+1),
                     cols=gene_col_idx, gridExpand=TRUE)
  openxlsx::setColWidths(wb, "Oscillating_Pathways",
                         cols=seq_len(ncol(conf_with_genes) - 1L), widths="auto")
  openxlsx::setColWidths(wb, "Oscillating_Pathways",
                         cols=gene_col_idx, widths=60)
  openxlsx::saveWorkbook(wb, xlsx_sig, overwrite=TRUE)
  cat(sprintf("  Significant only -> %s\n", xlsx_sig))
}

# =============================================================================
# 11. PATHWAY PLOTS (oscillating programs)
# =============================================================================

if (nrow(conf_phase) == 0) {
  cat("\n  No significant oscillating pathways to plot.\n")
  cat("  Try relaxing cosinor_pval_corr (e.g. 0.10) or expanding databases.\n")
} else {
  path_plot_dir <- file.path(ct_dir, paste0(focus_safe, "_allpathways_plots"))
  if (!dir.exists(path_plot_dir)) dir.create(path_plot_dir, recursive = TRUE)

  # ── Top-6 grid ──────────────────────────────────────────────────────────────
  top_pw   <- head(conf_phase$Pathway[order(conf_phase$pvalue_corr)], 6L)
  pw_plots <- lapply(top_pw, function(pw) {
    tryCatch(
      plot_pathway_single(auc_mat, path_res, data_sub@meta.data,
                          celltype_col, zt_col, tmeta, focus_ct,
                          pw, period12, use_violin = TRUE),
      error = function(e) { message("  Skip ", pw, ": ", e$message); NULL })
  })
  pw_plots <- Filter(Negate(is.null), pw_plots)

  if (length(pw_plots) > 0) {
    n_col   <- min(3L, length(pw_plots))
    grid_pw <- gridExtra::arrangeGrob(grobs = pw_plots, ncol = n_col)
    top_grid_png <- file.path(path_plot_dir,
      paste0(focus_safe, "_top_oscillating_pathways_grid.png"))
    ggplot2::ggsave(top_grid_png, grid_pw,
                    width=6*n_col, height=5*ceiling(length(pw_plots)/n_col),
                    dpi=150, bg="white")
    cat(sprintf("\n  Top-6 grid -> %s\n", top_grid_png))
  }

  # ── Batch save all significant pathways ────────────────────────────────────
  cat(sprintf("  Saving individual plots for %d oscillating pathways...\n",
              nrow(conf_phase)))
  save_batch_pathway_plots(
    auc_mat      = auc_mat,
    path_results = path_res,
    meta         = data_sub@meta.data,
    celltype_col = celltype_col,
    zt_col       = zt_col,
    tmeta        = tmeta,
    target_ct    = focus_ct,
    n_top        = nrow(conf_phase),
    period12     = period12,
    use_violin   = TRUE,
    outdir       = path_plot_dir
  )
  cat(sprintf("  Individual plots -> %s\n", path_plot_dir))
}

cat("\n=== Approach B complete ===\n")
cat(sprintf("Outputs -> %s\n", ct_dir))
cat("\nSummary:\n")
cat(sprintf("  Pathways screened : %d\n", nrow(path_res$stats)))
cat(sprintf("  Oscillating (sig) : %d\n", nrow(conf_phase)))

if (nrow(conf_phase) > 0){
   cat(sprintf("  Adj-p range       : %.4f -- %.4f\n",min(conf_phase$pvalue_corr)
}
              