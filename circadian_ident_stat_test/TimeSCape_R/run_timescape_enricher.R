# ============================================================================
# TimeSCape v0.2 -- Approach A: enricher() on confident circadian genes
#
# Logic:
#   1. Run cosinor to find confident circadian genes (same as always)
#   2. Run clusterProfiler::enricher() on those genes against msigdbr gene sets
#      -- one standard hypergeometric ORA, no phase binning, clean output
#   3. Score ALL cells for each enriched pathway (AUCell or AddModuleScore)
#   4. Run cosinor on pathway scores -> find which enriched pathways oscillate
#
# Key differences from the phase-bin approach:
#   - No phase binning, no phyper, no phase-restricted gene sets
#   - enricher() handles the ORA with proper BH correction
#   - Pathway genes are used in full (not intersected with a single time bin)
#   - AddModuleScore option for fast prototyping (no AUCell install needed)
#
# Required packages:
#   install.packages("msigdbr")
#   BiocManager::install(c("clusterProfiler", "AUCell"))  # AUCell optional
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
out_dir   <- file.path(base_dir, "TimeSCape_output_enricher")

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

celltype_col <- "sub_cell_types2"
zt_col       <- "ZT_time_str"
period12     <- FALSE
norm_str     <- "logcounts"
n_workers    <- 2L
focus_ct     <- "CD8+ T cells"
focus_safe   <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))

# ── SCORING METHOD ────────────────────────────────────────────────────────────
# TRUE  = AUCell  (rank-based, robust, slower, needs BiocManager::install("AUCell"))
# FALSE = AddModuleScore (Seurat built-in, fast, good for prototyping)
use_aucell <- FALSE

# ── ENRICHER PARAMETERS ───────────────────────────────────────────────────────
enrich_pval   <- 0.05     # p-value cutoff (enricher)
enrich_padj   <- 0.20     # adjusted p-value cutoff (BH) -- relaxed to capture more candidates
enrich_min_gs <- 10L      # minimum gene set size
enrich_max_gs <- 500L     # maximum gene set size

tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)
ct_dir <- file.path(out_dir, focus_safe)
if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

# =============================================================================
# 5. GENE-LEVEL CIRCADIAN ANALYSIS
# =============================================================================

cat(sprintf("\n=== Step 5: Cosinor -- '%s' ===\n", focus_ct))
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
  n_workers       = n_workers
)

key        <- names(results)[1]
T1         <- results[[key]]$T1
conf_mask  <- T1$pvalue < 0.05 & T1$pvalue_corr < 0.05
conf_genes <- T1$Genes[conf_mask]

cat(sprintf("\n  Confident circadian genes: %d / %d\n",
            length(conf_genes), nrow(T1)))
print(head(T1[conf_mask,
              c("Genes","Abs_Amp","Acrophase_24","pvalue","pvalue_corr")], 10))

# Top-6 gene grid
top_n_show    <- 6L
genes_to_show <- head(conf_genes, top_n_show)
gene_plots <- lapply(genes_to_show, function(g) {
  tryCatch(
    plot_gene_single(tmeta=tmeta, cust_cells=key, period12=period12,
                     cust_gene=g, print_scdata=TRUE, sce=data,
                     celltype_col=celltype_col, zt_col=zt_col,
                     use_violin=FALSE, outdir=ct_dir),
    error = function(e) { message("  Skip ", g, ": ", e$message); NULL })
})
gene_plots <- Filter(Negate(is.null), gene_plots)
if (length(gene_plots) > 0) {
  n_col  <- min(3L, length(gene_plots))
  grid_p <- gridExtra::arrangeGrob(grobs = gene_plots, ncol = n_col)
  ggplot2::ggsave(file.path(ct_dir, paste0(focus_safe, "_top_genes.png")),
                  grid_p, width=6*n_col, height=5*ceiling(top_n_show/n_col),
                  dpi=150, bg="white")
}

tryCatch(
  generate_heatmap(celltype=key, strict=TRUE, circ=FALSE,
                   period12=period12, outdir=ct_dir),
  error = function(e) message("  Heatmap skipped: ", e$message)
)

if (length(conf_genes) == 0)
  stop("No confident circadian genes found -- cannot proceed.")

# =============================================================================
# 6. BUILD TERM2GENE FROM msigdbr (mouse gene symbols)
# =============================================================================

cat("\n=== Step 6: Building TERM2GENE from msigdbr ===\n")

if (!requireNamespace("msigdbr", quietly = TRUE))
  stop("Install msigdbr:  install.packages('msigdbr')")

# Pull collections and combine into a single TERM2GENE data frame
# (gs_name | gene_symbol) as required by clusterProfiler::enricher()
t2g_kegg <- msigdbr::msigdbr(species = "Mus musculus",
                               collection = "C2", subcollection = "CP:KEGG_LEGACY") |>
  dplyr::select(gs_name, gene_symbol)

t2g_reactome <- msigdbr::msigdbr(species = "Mus musculus",
                                   collection = "C2", subcollection = "CP:REACTOME") |>
  dplyr::select(gs_name, gene_symbol)

t2g_gobp <- msigdbr::msigdbr(species = "Mus musculus",
                               collection = "C5", subcollection = "GO:BP") |>
  dplyr::select(gs_name, gene_symbol)

t2g_gomf <- msigdbr::msigdbr(species = "Mus musculus",
                               collection = "C5", subcollection = "GO:MF") |>
  dplyr::select(gs_name, gene_symbol)

TERM2GENE <- dplyr::bind_rows(t2g_kegg, t2g_reactome, t2g_gobp, t2g_gomf)
cat(sprintf("  TERM2GENE: %d gene-set entries across %d unique pathways\n",
            nrow(TERM2GENE), length(unique(TERM2GENE$gs_name))))

# =============================================================================
# 7. enricher() ORA ON CONFIDENT CIRCADIAN GENES
# =============================================================================

cat(sprintf("\n=== Step 7: enricher() ORA on %d confident genes ===\n",
            length(conf_genes)))

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  stop("Install:  BiocManager::install('clusterProfiler')")

enrich_res <- clusterProfiler::enricher(
  gene          = conf_genes,
  TERM2GENE     = TERM2GENE,
  universe      = T1$Genes,       # all tested genes = proper background
  pvalueCutoff  = enrich_pval,
  pAdjustMethod = "BH",
  qvalueCutoff  = enrich_padj,
  minGSSize     = enrich_min_gs,
  maxGSSize     = enrich_max_gs
)

if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0) {
  stop("enricher() returned no significant pathways. ",
       "Try relaxing enrich_padj (e.g. 0.5) or enrich_pval (e.g. 0.1).")
}

enrich_df <- as.data.frame(enrich_res)
cat(sprintf("  Significant pathways: %d\n", nrow(enrich_df)))
print(head(enrich_df[, c("ID","GeneRatio","BgRatio","pvalue","p.adjust","Count")], 15))

# Save enrichment table
enrich_xlsx <- file.path(ct_dir, paste0(focus_safe, "_enricher_ORA.xlsx"))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb  <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(fontColour="#FFFFFF", fgFill="#2F4F8F",
                                halign="center", textDecoration="bold")
  openxlsx::addWorksheet(wb, "Enricher_ORA")
  openxlsx::writeData(wb, "Enricher_ORA", enrich_df, rowNames=FALSE)
  openxlsx::addStyle(wb, "Enricher_ORA", hdr, rows=1,
                     cols=seq_len(ncol(enrich_df)), gridExpand=TRUE)
  openxlsx::setColWidths(wb, "Enricher_ORA",
                         cols=seq_len(ncol(enrich_df)), widths="auto")
  openxlsx::saveWorkbook(wb, enrich_xlsx, overwrite=TRUE)
  cat(sprintf("  Saved -> %s\n", enrich_xlsx))
}

# Build gene set list from enriched pathways
# (use the genes from TERM2GENE, not just the overlap genes)
enrich_gs <- lapply(enrich_df$ID, function(pw) {
  unique(TERM2GENE$gene_symbol[TERM2GENE$gs_name == pw])
})
names(enrich_gs) <- enrich_df$ID
cat(sprintf("  Gene sets for scoring: %d\n", length(enrich_gs)))

# =============================================================================
# 8. SCORE CELLS (AUCell or AddModuleScore) -- subset to focus cell type first
# =============================================================================

cat(sprintf("\n=== Step 8: Scoring cells (%s) ===\n",
            if (use_aucell) "AUCell" else "AddModuleScore"))

# CRITICAL: subset to focus cell type before scoring.
# AddModuleScore computes control gene backgrounds from ALL cells in the
# object -- mixing cell types biases the correction. AUCell is per-cell but
# subsetting still avoids wasting time on irrelevant cells.
focus_cells <- colnames(data)[data@meta.data[[celltype_col]] == focus_ct]
data_sub    <- data[, focus_cells]
cat(sprintf("  Subsetting to '%s': %d / %d cells\n",
            focus_ct, ncol(data_sub), ncol(data)))

score_cache <- file.path(ct_dir,
  paste0(focus_safe, "_enricher_scores_",
         if (use_aucell) "aucell" else "ams", ".rds"))

if (file.exists(score_cache)) {
  cat("  Loading cached scores...\n")
  auc_mat <- readRDS(score_cache)
  cat(sprintf("  %d pathways x %d cells\n", nrow(auc_mat), ncol(auc_mat)))

} else if (use_aucell) {
  # ── AUCell ────────────────────────────────────────────────────────────────
  auc_mat <- auc_score_cells(
    obj          = data_sub,
    genesets     = enrich_gs,
    use_norm     = TRUE,
    auc_max_rank = 0.05,
    n_cores      = 1L,
    min_gs_size  = 3L
  )
  saveRDS(auc_mat, score_cache)
  cat(sprintf("  Cached -> %s\n", score_cache))

} else {
  # ── AddModuleScore (Seurat built-in) ──────────────────────────────────────
  # Scores cells in batches of 50 pathways and collects results into a
  # pathways x cells matrix compatible with pathway_cosinor().
  # ctrl=5 (down from default 100) for speed; raise to 20-50 for publication.
  cat(sprintf("  Running AddModuleScore on %d pathways in batches...\n",
              length(enrich_gs)))

  batch_size <- 50L
  pw_names   <- names(enrich_gs)
  n_pw       <- length(pw_names)
  score_rows <- vector("list", n_pw)

  for (start in seq(1, n_pw, by = batch_size)) {
    end      <- min(start + batch_size - 1L, n_pw)
    batch_nm <- pw_names[start:end]
    batch_gs <- enrich_gs[batch_nm]

    # AddModuleScore adds columns "score1", "score2", ... to meta.data
    tmp <- Seurat::AddModuleScore(
      object   = data_sub,
      features = batch_gs,
      name     = "score",
      seed     = 42,
      ctrl     = 5        # lower for speed; use 20-50 for publication quality
    )
    # Extract the newly added columns (suffix 1..length(batch_gs))
    idx_cols <- paste0("score", seq_along(batch_gs))
    for (j in seq_along(batch_nm)) {
      score_rows[[start + j - 1L]] <- tmp@meta.data[[idx_cols[j]]]
    }
    cat(sprintf("  Scored %d / %d pathways...\r", end, n_pw))
  }
  cat("\n")

  # Assemble matrix: pathways x cells
  auc_mat <- do.call(rbind, score_rows)
  rownames(auc_mat) <- pw_names
  colnames(auc_mat) <- colnames(data_sub)
  saveRDS(auc_mat, score_cache)
  cat(sprintf("  Cached -> %s\n", score_cache))
}

# =============================================================================
# 9. PATHWAY COSINOR ON ENRICHED PATHWAY SCORES
# =============================================================================

cat(sprintf("\n=== Step 9: Pathway cosinor -- '%s' ===\n", focus_ct))

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
  path_res$stats$pvalue < 0.05 & path_res$stats$pvalue_corr < 0.05, ]
cat(sprintf("  Oscillating pathways: %d / %d\n",
            nrow(conf_phase), nrow(path_res$stats)))
if (nrow(conf_phase) > 0)
  print(head(conf_phase[, c("Pathway","Abs_Amp","Acrophase_24",
                             "pvalue","pvalue_corr")], 15))

# Save Excel
xlsx_out <- file.path(ct_dir, paste0(focus_safe, "_enricher_pathway_cosinor.xlsx"))
write_pathway_results(path_res, xlsx_out, celltype = focus_ct)
cat(sprintf("  Saved -> %s\n", xlsx_out))

# =============================================================================
# 10. PATHWAY PLOTS (confident oscillating pathways)
# =============================================================================

path_plot_dir <- file.path(ct_dir, paste0(focus_safe, "_enricher_pathway_plots"))
if (!dir.exists(path_plot_dir)) dir.create(path_plot_dir, recursive = TRUE)

if (nrow(conf_phase) > 0) {
  # Top-6 grid
  top_pw   <- head(conf_phase$Pathway, 6L)
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
    ggplot2::ggsave(
      file.path(path_plot_dir, paste0(focus_safe, "_top_pathways_grid.png")),
      grid_pw, width=6*n_col, height=5*ceiling(length(pw_plots)/n_col),
      dpi=150, bg="white")
  }

  # Batch save all confident
  save_batch_pathway_plots(
    auc_mat=auc_mat, path_results=path_res, meta=data_sub@meta.data,
    celltype_col=celltype_col, zt_col=zt_col, tmeta=tmeta,
    target_ct=focus_ct, n_top=nrow(conf_phase), period12=period12,
    use_violin=TRUE, outdir=path_plot_dir
  )
}

cat("\n=== Approach A complete ===\n")
cat(sprintf("Outputs -> %s\n", ct_dir))
