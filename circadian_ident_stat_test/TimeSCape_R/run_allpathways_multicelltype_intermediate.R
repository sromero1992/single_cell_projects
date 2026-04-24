# ============================================================================
# TimeSCape v0.2 -- Approach B (All-Pathway Screen) | MULTI-CELL-TYPE | INTERMEDIATE
#
# Logic:
#   For each cell type in the dataset:
#   1. Subset cells to that cell type
#   2. Score ALL msigdbr pathways (AUCell or AddModuleScore)
#   3. Run cosinor on pathway scores -> detect oscillating gene programs
#   4. Save Excel (stats + gene members) + top-6 plot grid
#
# Efficiency note:
#   TERM2GENE and the filtered gene set list (all_gs) are built ONCE before
#   the loop -- they depend only on the object's gene universe, not on cell
#   type. Only the subsetting, scoring, and cosinor run per cell type.
#
# Stage: INTERMEDIATE
# ============================================================================

library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)
library(dplyr)

set.seed(123)

# =============================================================================
# 1. PATHS  (change stage label here only)
# =============================================================================

stage    <- "intermediate"
base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm\\r_pre_process\\TimeSCape_R_tes"
src_path <- "C:\\Users\\selim\\Documentos\\vscode_working_dir\\single_cell_projects\\circadian_ident_stat_test\\TimeSCape_R\\R"
out_dir  <- file.path(base_dir, sprintf("TimeSCape_output_allpathways_%s", stage))

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
# 3. LOAD DATA  (change RDS filename here only)
# =============================================================================

cat("Loading Seurat object...\n")
data <- readRDS(file.path(base_dir,
  "circadian_mouse_breast_decontX_QC_2_23_2026_subann_intermediate.rds"))
cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(data), nrow(data)))

# =============================================================================
# 4. CELL TYPE MERGING  (collapse subtypes before analysis)
# =============================================================================

data$cell_types_merged <- as.character(data$sub_cell_types2)

# Merge CD8+ T cell subtypes
data$cell_types_merged[data$sub_cell_types2 %in%
  c("CD8+ T cells", "Cyc. CD8+ T cells", "TCF7+ CD8+ T cells")] <- "CD8+ T cells"

# Merge T-Reg subtypes
data$cell_types_merged[data$sub_cell_types2 %in%
  c("T-Reg cells", "TCF7+ T-Reg cells", "Cyc. T-Reg cells")] <- "T-Reg cells"

# Merge CD4+ T cell subtypes
data$cell_types_merged[data$sub_cell_types2 %in%
  c("CD4+ T cells", "TCF7+ CD4+ T cells")] <- "CD4+ T cells"

# Merge M2 macrophage subtypes
data$cell_types_merged[data$sub_cell_types2 %in%
  c("Pro-resolution M2", "IFNγ-act chkpt M2",
    "Ccr2+ Fibro-remodeling M2", "Immunoreg scavenger M2",
    "Cd163+/Mrc1+ M2")] <- "M2 macrophages"

# Merge tumor subtypes
data$cell_types_merged[data$sub_cell_types2 %in%
  c("Tumor cells", "Prol. Tumor cells")] <- "Tumor cells"

cat("\nMerged cell type counts:\n")
print(sort(table(data$cell_types_merged), decreasing = TRUE))

# =============================================================================
# 5. PARAMETERS
# =============================================================================

celltype_col <- "cell_types_merged"
zt_col       <- "ZT_time_str"
period12     <- FALSE
n_workers    <- 2L

# ── SCORING METHOD ────────────────────────────────────────────────────────────
# TRUE  = AUCell  (rank-based, robust; needs BiocManager::install("AUCell"))
# FALSE = AddModuleScore (fast prototyping; ctrl=5 for speed)
use_aucell <- TRUE

# ── PATHWAY DATABASE SELECTION ────────────────────────────────────────────────
use_kegg     <- TRUE    # ~200 pathways
use_reactome <- TRUE    # ~1,500 pathways
use_gobp     <- TRUE    # ~7,500 terms -- set FALSE for a quick run
use_gomf     <- FALSE   # ~1,000 terms -- optional

# ── GENE SET SIZE FILTER ──────────────────────────────────────────────────────
min_gs_size <- 10L
max_gs_size <- 500L

# ── COSINOR THRESHOLDS ────────────────────────────────────────────────────────
cosinor_pval      <- 0.05
cosinor_pval_corr <- 0.05

# ── MINIMUM CELLS PER ZT (skip cell types with too few cells) ─────────────────
min_cells_per_zt <- 5L

# ── CELL TYPE SUBSET (NULL = run all) ────────────────────────────────────────
# Set to a character vector to run only specific cell types, e.g.:
# ct_subset <- c("CD8+ T cells", "Endothelial")
ct_subset <- NULL

tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)

# Detect cell types to run
all_labels <- sort(unique(na.omit(as.character(data@meta.data[[celltype_col]]))))
ct_targets <- if (!is.null(ct_subset)) intersect(ct_subset, all_labels) else all_labels
cat(sprintf("\nCell types to process: %d\n", length(ct_targets)))
cat(paste0("  ", ct_targets, collapse="\n"), "\n")

# =============================================================================
# 6. BUILD TERM2GENE FROM msigdbr  -- done ONCE outside the loop
# =============================================================================

cat("\n=== Building TERM2GENE from msigdbr (once for all cell types) ===\n")

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
# 7. BUILD GENE SET LIST  -- done ONCE, based on genes in the full object
# =============================================================================

cat(sprintf("\n=== Building gene set list (size %d-%d) ===\n",
            min_gs_size, max_gs_size))

genes_in_obj <- rownames(data)

gs_sizes <- TERM2GENE |>
  dplyr::filter(gene_symbol %in% genes_in_obj) |>
  dplyr::group_by(gs_name) |>
  dplyr::summarise(n_in_obj = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_in_obj >= min_gs_size, n_in_obj <= max_gs_size)

cat(sprintf("  Pathways passing size filter: %d / %d\n",
            nrow(gs_sizes), length(unique(TERM2GENE$gs_name))))

all_gs <- lapply(gs_sizes$gs_name, function(pw) {
  genes <- unique(TERM2GENE$gene_symbol[TERM2GENE$gs_name == pw])
  intersect(genes, genes_in_obj)
})
names(all_gs) <- gs_sizes$gs_name
cat(sprintf("  Total gene sets to score per cell type: %d\n\n", length(all_gs)))

if (length(all_gs) == 0)
  stop("No gene sets passed the size filter.")

# Helper: attach gene members column to a stats data frame
add_gene_column <- function(df, gs_list) {
  df$Genes <- sapply(df$Pathway, function(pw) {
    g <- gs_list[[pw]]
    if (is.null(g) || length(g) == 0) NA_character_
    else paste(sort(g), collapse = "; ")
  })
  df
}

# =============================================================================
# 8. MAIN LOOP -- one cell type at a time
# =============================================================================

run_log <- data.frame(
  CellType          = character(),
  N_cells           = integer(),
  Pathways_scored   = integer(),
  Oscillating       = integer(),
  Skipped           = logical(),
  stringsAsFactors  = FALSE
)

for (focus_ct in ct_targets) {

  focus_safe <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))
  ct_dir     <- file.path(out_dir, focus_safe)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

  cat(sprintf("\n%s\n=== Cell type: %s ===\n%s\n",
              strrep("=", 60), focus_ct, strrep("=", 60)))

  # ── Subset to this cell type ────────────────────────────────────────────────
  focus_cells <- colnames(data)[data@meta.data[[celltype_col]] == focus_ct]
  data_sub    <- data[, focus_cells]
  cat(sprintf("  Cells: %d\n", ncol(data_sub)))

  # ── Check minimum cells per ZT ──────────────────────────────────────────────
  zt_counts <- table(data_sub@meta.data[[zt_col]])
  if (any(zt_counts < min_cells_per_zt)) {
    cat(sprintf("  SKIP: some ZT points have < %d cells (%s)\n",
                min_cells_per_zt,
                paste(names(zt_counts)[zt_counts < min_cells_per_zt],
                      collapse=", ")))
    run_log <- rbind(run_log, data.frame(
      CellType=focus_ct, N_cells=ncol(data_sub),
      Pathways_scored=0L, Oscillating=0L, Skipped=TRUE))
    next
  }

  # ── Score cells ─────────────────────────────────────────────────────────────
  score_cache <- file.path(ct_dir,
    paste0(focus_safe, "_allpathways_scores_",
           if (use_aucell) "aucell" else "ams", ".rds"))

  if (file.exists(score_cache)) {
    cat("  Loading cached scores...\n")
    auc_mat <- readRDS(score_cache)
    cat(sprintf("  %d pathways x %d cells\n", nrow(auc_mat), ncol(auc_mat)))

  } else if (use_aucell) {
    cat(sprintf("  AUCell: scoring %d pathways...\n", length(all_gs)))
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
    cat(sprintf("  AddModuleScore: %d pathways in batches of 50...\n", length(all_gs)))
    batch_size <- 50L
    pw_names   <- names(all_gs)
    n_pw       <- length(pw_names)
    score_rows <- vector("list", n_pw)

    for (start in seq(1, n_pw, by = batch_size)) {
      end      <- min(start + batch_size - 1L, n_pw)
      batch_nm <- pw_names[start:end]
      batch_gs <- all_gs[batch_nm]
      tmp <- Seurat::AddModuleScore(
        object = data_sub, features = batch_gs,
        name = "score", seed = 42, ctrl = 5
      )
      idx_cols <- paste0("score", seq_along(batch_gs))
      for (j in seq_along(batch_nm))
        score_rows[[start + j - 1L]] <- tmp@meta.data[[idx_cols[j]]]
      if ((end %% 500) == 0 || end == n_pw)
        cat(sprintf("  Scored %d / %d...\n", end, n_pw))
    }
    auc_mat <- do.call(rbind, score_rows)
    rownames(auc_mat) <- pw_names
    colnames(auc_mat) <- colnames(data_sub)
    saveRDS(auc_mat, score_cache)
    cat(sprintf("  Cached -> %s\n", score_cache))
  }

  # ── Pathway cosinor ─────────────────────────────────────────────────────────
  cat(sprintf("  Fitting cosinor for %d pathways...\n", nrow(auc_mat)))
  path_res <- tryCatch(
    pathway_cosinor(
      auc_mat      = auc_mat,
      meta         = data_sub@meta.data,
      celltype_col = celltype_col,
      zt_col       = zt_col,
      tmeta        = tmeta,
      target_ct    = focus_ct,
      period12     = period12
    ),
    error = function(e) {
      message("  ERROR in pathway_cosinor: ", e$message)
      NULL
    }
  )

  if (is.null(path_res)) {
    run_log <- rbind(run_log, data.frame(
      CellType=focus_ct, N_cells=ncol(data_sub),
      Pathways_scored=nrow(auc_mat), Oscillating=NA_integer_, Skipped=TRUE))
    next
  }

  conf_phase <- path_res$stats[
    path_res$stats$pvalue      < cosinor_pval &
    path_res$stats$pvalue_corr < cosinor_pval_corr, ]

  cat(sprintf("  Oscillating pathways: %d / %d\n",
              nrow(conf_phase), nrow(path_res$stats)))

  # ── Save Excel (all pathways + significant, both with gene lists) ────────────
  hdr     <- openxlsx::createStyle(fontColour="#FFFFFF", fgFill="#2F4F8F",
                                    halign="center", textDecoration="bold")
  wrap_st <- openxlsx::createStyle(wrapText=TRUE, valign="top")

  # All pathways
  all_stats      <- add_gene_column(path_res$stats[order(path_res$stats$pvalue_corr),], all_gs)
  gene_col_idx   <- which(names(all_stats) == "Genes")
  wb_all         <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb_all, "All_Pathways")
  openxlsx::writeData(wb_all, "All_Pathways", all_stats, rowNames=FALSE)
  openxlsx::addStyle(wb_all, "All_Pathways", hdr, rows=1,
                     cols=seq_len(ncol(all_stats)), gridExpand=TRUE)
  openxlsx::addStyle(wb_all, "All_Pathways", wrap_st,
                     rows=seq(2, nrow(all_stats)+1), cols=gene_col_idx, gridExpand=TRUE)
  openxlsx::setColWidths(wb_all, "All_Pathways",
                         cols=seq_len(ncol(all_stats)-1L), widths="auto")
  openxlsx::setColWidths(wb_all, "All_Pathways", cols=gene_col_idx, widths=60)
  openxlsx::saveWorkbook(wb_all,
    file.path(ct_dir, paste0(focus_safe, "_allpathways_cosinor_all.xlsx")),
    overwrite=TRUE)

  # Significant only
  if (nrow(conf_phase) > 0) {
    conf_genes_df  <- add_gene_column(conf_phase[order(conf_phase$pvalue_corr),], all_gs)
    gene_col_idx2  <- which(names(conf_genes_df) == "Genes")
    wb_sig         <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_sig, "Oscillating_Pathways")
    openxlsx::writeData(wb_sig, "Oscillating_Pathways", conf_genes_df, rowNames=FALSE)
    openxlsx::addStyle(wb_sig, "Oscillating_Pathways", hdr, rows=1,
                       cols=seq_len(ncol(conf_genes_df)), gridExpand=TRUE)
    openxlsx::addStyle(wb_sig, "Oscillating_Pathways", wrap_st,
                       rows=seq(2, nrow(conf_genes_df)+1), cols=gene_col_idx2, gridExpand=TRUE)
    openxlsx::setColWidths(wb_sig, "Oscillating_Pathways",
                           cols=seq_len(ncol(conf_genes_df)-1L), widths="auto")
    openxlsx::setColWidths(wb_sig, "Oscillating_Pathways", cols=gene_col_idx2, widths=60)
    openxlsx::saveWorkbook(wb_sig,
      file.path(ct_dir, paste0(focus_safe, "_allpathways_cosinor_significant.xlsx")),
      overwrite=TRUE)
    cat(sprintf("  Excel saved -> %s/\n", ct_dir))
  }

  # ── Top-6 pathway plot grid ──────────────────────────────────────────────────
  if (nrow(conf_phase) > 0) {
    top_pw   <- head(conf_phase$Pathway[order(conf_phase$pvalue_corr)], 6L)
    pw_plots <- lapply(top_pw, function(pw) {
      tryCatch(
        plot_pathway_single(auc_mat, path_res, data_sub@meta.data,
                            celltype_col, zt_col, tmeta, focus_ct,
                            pw, period12, use_violin = TRUE),
        error = function(e) { message("  Skip plot ", pw, ": ", e$message); NULL })
    })
    pw_plots <- Filter(Negate(is.null), pw_plots)
    if (length(pw_plots) > 0) {
      n_col   <- min(3L, length(pw_plots))
      grid_pw <- gridExtra::arrangeGrob(grobs=pw_plots, ncol=n_col)
      ggplot2::ggsave(
        file.path(ct_dir, paste0(focus_safe, "_top_oscillating_pathways_grid.png")),
        grid_pw, width=6*n_col, height=5*ceiling(length(pw_plots)/n_col),
        dpi=150, bg="white")
      cat(sprintf("  Top-6 grid saved.\n"))
    }
  }

  run_log <- rbind(run_log, data.frame(
    CellType        = focus_ct,
    N_cells         = ncol(data_sub),
    Pathways_scored = nrow(auc_mat),
    Oscillating     = nrow(conf_phase),
    Skipped         = FALSE))
}

# =============================================================================
# 9. RUN SUMMARY
# =============================================================================

cat(sprintf("\n%s\n=== Run complete: %s stage ===\n%s\n",
            strrep("=", 60), stage, strrep("=", 60)))
print(run_log)

# Save run log
write.csv(run_log,
  file.path(out_dir, sprintf("run_log_%s.csv", stage)),
  row.names = FALSE)
cat(sprintf("\nAll outputs -> %s\n", out_dir))
