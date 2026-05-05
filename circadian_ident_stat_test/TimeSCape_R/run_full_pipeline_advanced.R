# ============================================================================
# TimeSCape v0.2 -- Full Circadian Pipeline | ALL CELL TYPES | ADVANCED STAGE
#
# Per cell type this script runs three sequential analyses:
#
#   STEP A  Gene-level cosinor
#           run_timescape() -> saves all CSVs + confident spreadsheet
#           generate_heatmap() -> saves heatmap PNG
#           (no individual gene plots -- use run_timescape_test.R for those)
#
#   STEP B  Approach A: enricher() ORA on confident circadian genes
#           clusterProfiler::enricher() on conf_genes vs TERM2GENE
#           -> Score subsetted cells (AUCell or AddModuleScore)
#           -> pathway_cosinor() -> save Excel (with gene members) + top-6 grid
#
#   STEP C  Approach B: All-pathway circadian screen
#           Score subsetted cells for ALL size-filtered pathways
#           -> pathway_cosinor() -> save Excel (with gene members) + top-6 grid
#
# Efficiency:
#   TERM2GENE and all_gs (size-filtered gene set list) are built ONCE before
#   the loop and shared across all cell types.
#
# Stage: ADVANCED
# ============================================================================

library(Seurat)
library(future)
library(future.apply)
library(minpack.lm)
library(dplyr)

set.seed(123)

# =============================================================================
# 1. PATHS  -- only these three lines change between stage copies
# =============================================================================

stage    <- "advanced"
base_dir <- "Z:\\selim_working_dir\\2025_sato_anestacia_circadian_rhythm\\r_pre_process\\TimeSCape_R_tes"
src_path <- "C:\\Users\\selim\\Documentos\\vscode_working_dir\\single_cell_projects\\circadian_ident_stat_test\\TimeSCape_R\\R"
rds_file <- "circadian_mouse_breast_decontX_QC_2_23_2026_subann_advanced.rds"
out_dir  <- file.path(base_dir, sprintf("TimeSCape_output_%s", stage))

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
data <- readRDS(file.path(base_dir, rds_file))
cat(sprintf("  Loaded: %d cells x %d genes\n", ncol(data), nrow(data)))

# =============================================================================
# 4. CELL TYPE MERGING  -- collapse subtypes into clean groups
# =============================================================================

data$cell_types_merged <- as.character(data$sub_cell_types2)

data$cell_types_merged[data$sub_cell_types2 %in%
  c("CD8+ T cells", "Cyc. CD8+ T cells", "TCF7+ CD8+ T cells")] <- "CD8+ T cells"

data$cell_types_merged[data$sub_cell_types2 %in%
  c("T-Reg cells", "TCF7+ T-Reg cells", "Cyc. T-Reg cells")] <- "T-Reg cells"

data$cell_types_merged[data$sub_cell_types2 %in%
  c("CD4+ T cells", "TCF7+ CD4+ T cells")] <- "CD4+ T cells"

data$cell_types_merged[data$sub_cell_types2 %in%
  c("Pro-resolution M2", "IFNγ-act chkpt M2",
    "Ccr2+ Fibro-remodeling M2", "Immunoreg scavenger M2",
    "Cd163+/Mrc1+ M2")] <- "M2 macrophages"

data$cell_types_merged[data$sub_cell_types2 %in%
  c("Tumor cells", "Prol. Tumor cells")] <- "Tumor cells"

cat("\nMerged cell type counts:\n")
print(sort(table(data$cell_types_merged), decreasing = TRUE))

# =============================================================================
# 5. GLOBAL PARAMETERS
# =============================================================================

celltype_col <- "cell_types_merged"
zt_col       <- "ZT_time_str"
period12     <- FALSE
norm_str     <- "logcounts"
n_workers    <- 2L

# ── Scoring method (Step B & C) ───────────────────────────────────────────────
# TRUE  = AUCell  (rank-based, robust; BiocManager::install("AUCell"))
# FALSE = AddModuleScore (fast; ctrl=5 for speed, raise to 20-50 for publication)
use_aucell <- TRUE

# ── Pathway databases ─────────────────────────────────────────────────────────
use_kegg     <- TRUE    # ~200 pathways
use_reactome <- TRUE    # ~1,500 pathways
use_gobp     <- TRUE    # ~7,500 terms (set FALSE for a quick test run)
use_gomf     <- FALSE   # ~1,000 terms (overlaps GO:BP; optional)

# ── Gene set size filter (applied after intersecting with object genes) ───────
min_gs_size <- 10L
max_gs_size <- 500L

# ── Step B: enricher() ORA parameters ────────────────────────────────────────
enrich_pval   <- 0.05
enrich_padj   <- 0.20   # relaxed BH cutoff -- tightened later by pathway cosinor
enrich_min_gs <- 10L
enrich_max_gs <- 500L

# ── Cosinor significance thresholds (Steps B & C) ────────────────────────────
cosinor_pval      <- 0.05
cosinor_pval_corr <- 0.05

# ── Minimum cells per ZT to accept a cell type ───────────────────────────────
min_cells_per_zt <- 5L

# ── Cell type subset (NULL = run all) ────────────────────────────────────────
# e.g. ct_subset <- c("CD8+ T cells", "Endothelial")
ct_subset <- NULL

tmeta <- build_tmeta(data, zt_col = zt_col)
cat("\nParsed tmeta:\n"); print(tmeta)

all_labels <- sort(unique(na.omit(as.character(data@meta.data[[celltype_col]]))))
ct_targets <- if (!is.null(ct_subset)) intersect(ct_subset, all_labels) else all_labels
cat(sprintf("\nCell types to process (%d):\n", length(ct_targets)))
cat(paste0("  ", ct_targets, collapse = "\n"), "\n")

# =============================================================================
# 6. BUILD TERM2GENE FROM msigdbr  -- once, shared across all cell types
# =============================================================================

cat("\n=== Building TERM2GENE (once for all cell types) ===\n")

if (!requireNamespace("msigdbr", quietly = TRUE))
  stop("Install msigdbr:  install.packages('msigdbr')")

t2g_list <- list()
if (use_kegg) {
  cat("  Pulling KEGG_LEGACY...\n")
  t2g_list$kegg <- msigdbr::msigdbr(species = "Mus musculus",
    collection = "C2", subcollection = "CP:KEGG_LEGACY") |>
    dplyr::select(gs_name, gene_symbol)
}
if (use_reactome) {
  cat("  Pulling Reactome...\n")
  t2g_list$reactome <- msigdbr::msigdbr(species = "Mus musculus",
    collection = "C2", subcollection = "CP:REACTOME") |>
    dplyr::select(gs_name, gene_symbol)
}
if (use_gobp) {
  cat("  Pulling GO:BP...\n")
  t2g_list$gobp <- msigdbr::msigdbr(species = "Mus musculus",
    collection = "C5", subcollection = "GO:BP") |>
    dplyr::select(gs_name, gene_symbol)
}
if (use_gomf) {
  cat("  Pulling GO:MF...\n")
  t2g_list$gomf <- msigdbr::msigdbr(species = "Mus musculus",
    collection = "C5", subcollection = "GO:MF") |>
    dplyr::select(gs_name, gene_symbol)
}

TERM2GENE <- dplyr::bind_rows(t2g_list)
cat(sprintf("  %d entries, %d unique pathways\n",
            nrow(TERM2GENE), length(unique(TERM2GENE$gs_name))))

# =============================================================================
# 7. BUILD ALL_GS (size-filtered gene set list)  -- once, for Step C
# =============================================================================

cat(sprintf("\n=== Building gene set list (size %d-%d) ===\n",
            min_gs_size, max_gs_size))

genes_in_obj <- rownames(data)

gs_sizes <- TERM2GENE |>
  dplyr::filter(gene_symbol %in% genes_in_obj) |>
  dplyr::group_by(gs_name) |>
  dplyr::summarise(n_in_obj = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_in_obj >= min_gs_size, n_in_obj <= max_gs_size)

all_gs <- lapply(gs_sizes$gs_name, function(pw) {
  intersect(unique(TERM2GENE$gene_symbol[TERM2GENE$gs_name == pw]), genes_in_obj)
})
names(all_gs) <- gs_sizes$gs_name

cat(sprintf("  %d / %d pathways pass size filter\n",
            length(all_gs), length(unique(TERM2GENE$gs_name))))

# Helper: attach gene members as a semicolon-separated column
add_gene_column <- function(df, gs_list) {
  df$Genes <- sapply(df$Pathway, function(pw) {
    g <- gs_list[[pw]]
    if (is.null(g) || length(g) == 0) NA_character_
    else paste(sort(g), collapse = "; ")
  })
  df
}

# Helper: write a stats data frame to a formatted Excel sheet
write_excel_sheet <- function(wb, sheet_name, df, gene_col = "Genes") {
  hdr     <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#2F4F8F",
                                    halign = "center", textDecoration = "bold")
  wrap_st <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, rowNames = FALSE)
  openxlsx::addStyle(wb, sheet_name, hdr, rows = 1,
                     cols = seq_len(ncol(df)), gridExpand = TRUE)
  if (gene_col %in% names(df)) {
    gi <- which(names(df) == gene_col)
    openxlsx::addStyle(wb, sheet_name, wrap_st,
                       rows = seq(2, nrow(df) + 1), cols = gi, gridExpand = TRUE)
    openxlsx::setColWidths(wb, sheet_name,
                           cols = seq_len(ncol(df))[-gi], widths = "auto")
    openxlsx::setColWidths(wb, sheet_name, cols = gi, widths = 60)
  } else {
    openxlsx::setColWidths(wb, sheet_name,
                           cols = seq_len(ncol(df)), widths = "auto")
  }
}

# Helper: score cells, using cache
score_cells_cached <- function(data_sub, gs_list, cache_path,
                                use_aucell, min_gs_size, ctrl = 5L) {
  if (file.exists(cache_path)) {
    cat("    Loading cached scores...\n")
    mat <- readRDS(cache_path)
    cat(sprintf("    %d pathways x %d cells\n", nrow(mat), ncol(mat)))
    return(mat)
  }

  if (use_aucell) {
    cat(sprintf("    AUCell: scoring %d pathways...\n", length(gs_list)))
    mat <- auc_score_cells(obj = data_sub, genesets = gs_list,
                           use_norm = TRUE, auc_max_rank = 0.05,
                           n_cores = 1L, min_gs_size = min_gs_size)
  } else {
    cat(sprintf("    AddModuleScore: %d pathways in batches of 50...\n", length(gs_list)))
    pw_names   <- names(gs_list)
    n_pw       <- length(pw_names)
    score_rows <- vector("list", n_pw)
    for (start in seq(1, n_pw, by = 50L)) {
      end      <- min(start + 49L, n_pw)
      batch_nm <- pw_names[start:end]
      batch_gs <- gs_list[batch_nm]
      tmp <- Seurat::AddModuleScore(object = data_sub, features = batch_gs,
                                    name = "score", seed = 42L, ctrl = ctrl)
      for (j in seq_along(batch_nm))
        score_rows[[start + j - 1L]] <- tmp@meta.data[[paste0("score", j)]]
      if ((end %% 500) == 0 || end == n_pw)
        cat(sprintf("    Scored %d / %d pathways...\n", end, n_pw))
    }
    mat <- do.call(rbind, score_rows)
    rownames(mat) <- pw_names
    colnames(mat) <- colnames(data_sub)
  }

  saveRDS(mat, cache_path)
  cat(sprintf("    Cached -> %s\n", cache_path))
  mat
}

# =============================================================================
# 8. MAIN LOOP
# =============================================================================

# Core clock gene reference (mouse symbols) -- used by GRN anchor selection
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

run_log <- data.frame(
  CellType        = character(),
  N_cells         = integer(),
  Conf_genes      = integer(),
  EnrichA_pathways = integer(),
  OscA_pathways   = integer(),
  OscB_pathways   = integer(),
  Skipped         = logical(),
  stringsAsFactors = FALSE
)

# Accumulate T1 tables for the post-loop clock acrophase polar plot
clock_results <- list()

for (focus_ct in ct_targets) {

  focus_safe <- gsub("[^[:alnum:]_]", "_", trimws(focus_ct))
  ct_dir     <- file.path(out_dir, focus_safe)
  dir_gene   <- file.path(ct_dir, "01_gene_circadian")   # Step A: gene cosinor
  dir_pwA    <- file.path(ct_dir, "02_pathway_A")         # Step B: enricher ORA pathway
  dir_pwB    <- file.path(ct_dir, "03_pathway_B")         # Step C: all-pathway screen
  # 04_GRN created on demand inside the GRN block
  for (d in c(ct_dir, dir_gene, dir_pwA, dir_pwB))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)

  cat(sprintf("\n%s\n### %s ###\n%s\n",
              strrep("=", 64), focus_ct, strrep("=", 64)))

  # ── Subset + minimum cells check ────────────────────────────────────────────
  focus_cells <- colnames(data)[data@meta.data[[celltype_col]] == focus_ct]
  data_sub    <- data[, focus_cells]
  cat(sprintf("  Cells: %d\n", ncol(data_sub)))

  zt_counts <- table(data_sub@meta.data[[zt_col]])
  if (any(zt_counts < min_cells_per_zt)) {
    cat(sprintf("  SKIP: ZT points with < %d cells: %s\n",
                min_cells_per_zt,
                paste(names(zt_counts)[zt_counts < min_cells_per_zt], collapse = ", ")))
    run_log <- rbind(run_log, data.frame(CellType = focus_ct, N_cells = ncol(data_sub),
      Conf_genes = 0L, EnrichA_pathways = 0L, OscA_pathways = 0L,
      OscB_pathways = 0L, Skipped = TRUE))
    next
  }

  conf_genes       <- character(0)
  T1               <- NULL
  key              <- NULL
  n_enrichA        <- 0L
  n_oscA           <- 0L
  n_oscB           <- 0L
  # Pre-initialise Step-B objects so the GRN block can reference them safely
  enrich_df        <- data.frame(stringsAsFactors = FALSE)
  enrich_gs        <- list()
  conf_A           <- data.frame(stringsAsFactors = FALSE)
  auc_A            <- NULL

  # ============================================================
  #  STEP A — Gene-level cosinor
  # ============================================================

  cat(sprintf("\n  -- Step A: Gene-level cosinor --\n"))

  step_A_ok <- tryCatch({
    results <- run_timescape(
      obj             = data,
      celltype_col    = celltype_col,
      zt_col          = zt_col,
      tmeta           = tmeta,
      rm_low_conf     = TRUE,
      period12        = period12,
      custom_celltype = focus_ct,
      plot_heat       = FALSE,       # heatmap generated separately below
      norm_str        = norm_str,
      outdir          = out_dir,
      n_workers       = n_workers
    )

    key        <- names(results)[1]
    T1         <- results[[key]]$T1
    conf_mask  <- T1$pvalue < 0.05 & T1$pvalue_corr < 0.05
    conf_genes <- T1$Genes[conf_mask]

    cat(sprintf("    Confident circadian genes: %d / %d\n",
                length(conf_genes), nrow(T1)))
    if (length(conf_genes) > 0)
      print(head(T1[conf_mask, c("Genes","Abs_Amp","Acrophase_24",
                                  "pvalue","pvalue_corr")], 8))

    # Save confident genes as a clean standalone Excel
    if (length(conf_genes) > 0 &&
        requireNamespace("openxlsx", quietly = TRUE)) {
      wb_genes <- openxlsx::createWorkbook()
      write_excel_sheet(wb_genes, "Confident_Genes",
                        T1[conf_mask, ][order(T1[conf_mask, ]$Acrophase_24), ],
                        gene_col = "Genes")
      openxlsx::saveWorkbook(wb_genes,
        file.path(dir_gene, paste0(focus_safe, "_confident_genes.xlsx")),
        overwrite = TRUE)
    }

    # Heatmap (confident genes, sorted by acrophase)
    tryCatch(
      generate_heatmap(celltype = key, strict = TRUE, circ = FALSE,
                       period12 = period12, outdir = dir_gene),
      error = function(e) message("    Heatmap skipped: ", e$message)
    )

    # Accumulate T1 for post-loop clock acrophase polar plot
    clock_results[[focus_ct]] <- T1

    TRUE
  }, error = function(e) {
    message("  ERROR in Step A: ", e$message)
    FALSE
  })

  if (!step_A_ok) {
    run_log <- rbind(run_log, data.frame(CellType = focus_ct, N_cells = ncol(data_sub),
      Conf_genes = 0L, EnrichA_pathways = 0L, OscA_pathways = 0L,
      OscB_pathways = 0L, Skipped = TRUE))
    next
  }

  # ============================================================
  #  STEP B — Approach A: enricher() ORA -> score -> cosinor
  # ============================================================

  cat(sprintf("\n  -- Step B: Approach A (enricher ORA -> pathway cosinor) --\n"))

  if (length(conf_genes) < 5L) {
    cat(sprintf("    SKIP: only %d confident genes (need >= 5 for ORA)\n",
                length(conf_genes)))
  } else {
    tryCatch({
      if (!requireNamespace("clusterProfiler", quietly = TRUE))
        stop("Install:  BiocManager::install('clusterProfiler')")

      enrich_res <- clusterProfiler::enricher(
        gene          = conf_genes,
        TERM2GENE     = TERM2GENE,
        universe      = T1$Genes,
        pvalueCutoff  = enrich_pval,
        pAdjustMethod = "BH",
        qvalueCutoff  = enrich_padj,
        minGSSize     = enrich_min_gs,
        maxGSSize     = enrich_max_gs
      )

      if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0) {
        cat("    enricher() returned no significant pathways -- Step B skipped.\n")
        cat("    Tip: relax enrich_padj (e.g. 0.50) at the top of the script.\n")
      } else {
        enrich_df <- as.data.frame(enrich_res)
        n_enrichA <- nrow(enrich_df)
        cat(sprintf("    Enriched pathways: %d\n", n_enrichA))

        # Save enrichment ORA table
        if (requireNamespace("openxlsx", quietly = TRUE)) {
          wb_ora <- openxlsx::createWorkbook()
          write_excel_sheet(wb_ora, "Enricher_ORA", enrich_df, gene_col = "geneID")
          openxlsx::saveWorkbook(wb_ora,
            file.path(dir_pwA, paste0(focus_safe, "_enricher_ORA.xlsx")),
            overwrite = TRUE)
        }

        # Gene sets for enriched pathways (full pathway membership, not just overlap)
        enrich_gs <- lapply(enrich_df$ID, function(pw)
          intersect(unique(TERM2GENE$gene_symbol[TERM2GENE$gs_name == pw]),
                    genes_in_obj))
        names(enrich_gs) <- enrich_df$ID

        # Score subsetted cells
        cache_A <- file.path(dir_pwA, paste0(focus_safe, "_enricherA_scores_",
                             if (use_aucell) "aucell" else "ams", ".rds"))
        auc_A <- score_cells_cached(data_sub, enrich_gs, cache_A,
                                     use_aucell, min_gs_size)

        # Pathway cosinor
        path_A <- pathway_cosinor(
          auc_mat      = auc_A,
          meta         = data_sub@meta.data,
          celltype_col = celltype_col,
          zt_col       = zt_col,
          tmeta        = tmeta,
          target_ct    = focus_ct,
          period12     = period12
        )

        conf_A  <- path_A$stats[path_A$stats$pvalue < cosinor_pval &
                                 path_A$stats$pvalue_corr < cosinor_pval_corr, ]
        n_oscA  <- nrow(conf_A)
        cat(sprintf("    Oscillating pathways (Approach A): %d / %d\n",
                    n_oscA, nrow(path_A$stats)))

        # Save Excel: all + significant (with gene members)
        if (requireNamespace("openxlsx", quietly = TRUE)) {
          all_A   <- add_gene_column(path_A$stats[order(path_A$stats$pvalue_corr),], enrich_gs)
          wb_A    <- openxlsx::createWorkbook()
          write_excel_sheet(wb_A, "All_Pathways", all_A)
          if (n_oscA > 0) {
            sig_A <- add_gene_column(conf_A[order(conf_A$pvalue_corr),], enrich_gs)
            write_excel_sheet(wb_A, "Oscillating_Pathways", sig_A)
          }
          openxlsx::saveWorkbook(wb_A,
            file.path(dir_pwA, paste0(focus_safe, "_approachA_pathway_cosinor.xlsx")),
            overwrite = TRUE)
        }

        # Top-6 pathway grid
        if (n_oscA > 0) {
          top_pw <- head(conf_A$Pathway[order(conf_A$pvalue_corr)], 6L)
          pw_plots <- Filter(Negate(is.null), lapply(top_pw, function(pw) {
            tryCatch(plot_pathway_single(auc_A, path_A, data_sub@meta.data,
                       celltype_col, zt_col, tmeta, focus_ct, pw, period12,
                       use_violin = TRUE),
                     error = function(e) { message("  Skip ", pw); NULL })
          }))
          if (length(pw_plots) > 0) {
            nc  <- min(3L, length(pw_plots))
            ggplot2::ggsave(
              file.path(dir_pwA, paste0(focus_safe, "_approachA_top_pathways.png")),
              gridExtra::arrangeGrob(grobs = pw_plots, ncol = nc),
              width = 6*nc, height = 5*ceiling(length(pw_plots)/nc),
              dpi = 150, bg = "white")
          }
        }
      }
    }, error = function(e) message("  ERROR in Step B: ", e$message))
  }

  # ============================================================
  #  STEP C — Approach B: all pathways -> score -> cosinor
  # ============================================================

  cat(sprintf("\n  -- Step C: Approach B (all pathways -> pathway cosinor) --\n"))
  cat(sprintf("    Scoring %d pathways on %d cells...\n",
              length(all_gs), ncol(data_sub)))

  tryCatch({
    cache_B <- file.path(dir_pwB, paste0(focus_safe, "_allpathwaysB_scores_",
                         if (use_aucell) "aucell" else "ams", ".rds"))
    auc_B <- score_cells_cached(data_sub, all_gs, cache_B, use_aucell, min_gs_size)

    path_B <- pathway_cosinor(
      auc_mat      = auc_B,
      meta         = data_sub@meta.data,
      celltype_col = celltype_col,
      zt_col       = zt_col,
      tmeta        = tmeta,
      target_ct    = focus_ct,
      period12     = period12
    )

    conf_B <- path_B$stats[path_B$stats$pvalue < cosinor_pval &
                            path_B$stats$pvalue_corr < cosinor_pval_corr, ]
    n_oscB <- nrow(conf_B)
    cat(sprintf("    Oscillating pathways (Approach B): %d / %d\n",
                n_oscB, nrow(path_B$stats)))

    if (nrow(conf_B) > 0)
      print(head(conf_B[order(conf_B$pvalue_corr),
                         c("Pathway","Abs_Amp","Acrophase_24",
                           "pvalue","pvalue_corr")], 10))

    # Save Excel: all + significant (with gene members)
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      all_B   <- add_gene_column(path_B$stats[order(path_B$stats$pvalue_corr),], all_gs)
      wb_B    <- openxlsx::createWorkbook()
      write_excel_sheet(wb_B, "All_Pathways", all_B)
      if (n_oscB > 0) {
        sig_B <- add_gene_column(conf_B[order(conf_B$pvalue_corr),], all_gs)
        write_excel_sheet(wb_B, "Oscillating_Pathways", sig_B)
      }
      openxlsx::saveWorkbook(wb_B,
        file.path(dir_pwB, paste0(focus_safe, "_approachB_allpathways_cosinor.xlsx")),
        overwrite = TRUE)
    }

    # Top-6 pathway grid
    if (n_oscB > 0) {
      top_pw <- head(conf_B$Pathway[order(conf_B$pvalue_corr)], 6L)
      pw_plots <- Filter(Negate(is.null), lapply(top_pw, function(pw) {
        tryCatch(plot_pathway_single(auc_B, path_B, data_sub@meta.data,
                   celltype_col, zt_col, tmeta, focus_ct, pw, period12,
                   use_violin = TRUE),
                 error = function(e) { message("  Skip ", pw); NULL })
      }))
      if (length(pw_plots) > 0) {
        nc  <- min(3L, length(pw_plots))
        ggplot2::ggsave(
          file.path(dir_pwB, paste0(focus_safe, "_approachB_top_pathways.png")),
          gridExtra::arrangeGrob(grobs = pw_plots, ncol = nc),
          width = 6*nc, height = 5*ceiling(length(pw_plots)/nc),
          dpi = 150, bg = "white")
      }
    }
  }, error = function(e) message("  ERROR in Step C: ", e$message))


  # ============================================================
  #  STEP D — GRN time series (hub co-expression approach)
  # ============================================================
  # Gene pool: conf_genes (Step A) + clock genes present in data
  # Hub selection: global Pearson correlation across ALL cells pooled,
  #   top hub_pct% by degree (number of significant edges)
  # GRN nodes: hub_genes ∩ pathway_genes (per oscillating pathway)
  #   -> Approach-A pathways (ORA-enriched, conf_A)  -> 04_GRN/approach_A/
  #   -> Approach-B pathways (all-pathway screen, conf_B) -> 04_GRN/approach_B/

  tryCatch({
    grn_pool <- unique(c(
      conf_genes,
      clock_genes_ref[clock_genes_ref %in% rownames(data)]
    ))
    cat(sprintf("  [GRN] Gene pool: %d genes\n", length(grn_pool)))

    hub_result <- select_hub_genes(
      obj          = data_sub,
      gene_pool    = grn_pool,
      target_ct    = focus_ct,
      celltype_col = celltype_col,
      use_norm     = TRUE,
      cor_thresh   = 0.30,
      p_thresh     = 0.05,
      hub_pct      = 0.10,
      min_hub      = 5L
    )
    hub_genes  <- hub_result$hub_genes
    hub_circ   <- hub_genes[hub_genes %in% conf_genes]
    hub_clock  <- hub_genes[hub_genes %in% clock_genes_ref]
    hub_unique <- setdiff(hub_genes, c(conf_genes, clock_genes_ref))

    cat(sprintf("  [GRN] Hubs: %d total\n", length(hub_genes)))
    if (length(hub_circ)   > 0) cat(sprintf("    Circadian: %s\n", paste(hub_circ,   collapse = ", ")))
    if (length(hub_clock)  > 0) cat(sprintf("    Clock    : %s\n", paste(hub_clock,  collapse = ", ")))
    if (length(hub_unique) > 0) cat(sprintf("    Other    : %s\n", paste(hub_unique, collapse = ", ")))

    .run_grn_batch <- function(osc_df, gs_list, grn_dir, label) {
      if (nrow(osc_df) == 0) {
        cat(sprintf("  [GRN %s] No oscillating pathways.\n", label)); return(invisible(NULL))
      }
      if (!dir.exists(grn_dir)) dir.create(grn_dir, recursive = TRUE)
      pws     <- osc_df$Pathway[order(osc_df$pvalue_corr)]
      n_built <- 0L
      for (pw in pws) {
        pw_genes  <- gs_list[[pw]]
        grn_nodes <- intersect(hub_genes, pw_genes)
        if (length(grn_nodes) < 3L) next
        pw_safe     <- substr(gsub("[^[:alnum:]_]", "_", pw), 1, 60)
        grn_outfile <- file.path(grn_dir,
          sprintf("%s_GRN_%s_%s.png", focus_safe, label, pw_safe))
        cat(sprintf("  [GRN %s] %s  (%d nodes)\n", label, pw, length(grn_nodes)))
        tryCatch(
          plot_grn_timeseries(
            obj            = data_sub,
            circ_genes     = grn_nodes,
            pathway_genes  = character(0),
            meta           = data_sub@meta.data,
            celltype_col   = celltype_col,
            zt_col         = zt_col,
            tmeta          = tmeta,
            target_ct      = focus_ct,
            cor_thresh     = 0.20,
            p_thresh       = 0.05,
            use_norm       = TRUE,
            outfile        = grn_outfile,
            ncol           = NULL,
            node_size      = 4,
            label_size     = 4.5,
            edge_width_max = 3.0,
            zt_title_size  = 18,
            zt_title_hjust = 0.5
          ),
          error = function(e) message("    Skip GRN [", pw, "]: ", e$message)
        )
        n_built <- n_built + 1L
      }
      cat(sprintf("  [GRN %s] %d GRN(s) -> %s\n", label, n_built, grn_dir))
    }

    grn_base <- file.path(ct_dir, "04_GRN")

    # Approach A: hub genes ∩ ORA-enriched oscillating pathways (Step B)
    .run_grn_batch(conf_A, enrich_gs, file.path(grn_base, "approach_A"), "A")

    # Approach B: hub genes ∩ all-pathway oscillating pathways (Step C)
    .run_grn_batch(conf_B, all_gs,    file.path(grn_base, "approach_B"), "B")

  }, error = function(e) message("  ERROR in Step D (GRN): ", e$message))


  # ── Log this cell type ───────────────────────────────────────────────────────
  run_log <- rbind(run_log, data.frame(
    CellType         = focus_ct,
    N_cells          = ncol(data_sub),
    Conf_genes       = length(conf_genes),
    EnrichA_pathways = n_enrichA,
    OscA_pathways    = n_oscA,
    OscB_pathways    = n_oscB,
    Skipped          = FALSE
  ))

  cat(sprintf("\n  Done: %d conf genes | %d enriched (A) | %d osc-A | %d osc-B\n",
              length(conf_genes), n_enrichA, n_oscA, n_oscB))
}

# =============================================================================
# 9. CLOCK ACROPHASE POLAR PLOT  -- all cell types combined
# =============================================================================

if (length(clock_results) > 0) {
  cat("\n=== Generating clock gene acrophase polar plot ===\n")
  # ── Customisation options ──────────────────────────────────────────────────
  # gene_list (default NULL = built-in core clock genes):
  #   Pass a character vector to plot any genes instead, e.g.:
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
  #   Cell types not matched are placed in "Other" automatically.
  # ──────────────────────────────────────────────────────────────────────────
  tryCatch(
    plot_clock_acrophase(
      results_list = clock_results,
      stage        = stage,
      outfile      = file.path(out_dir, sprintf("clock_acrophase_%s.png", stage)),
      dpi          = 300,
      strict       = TRUE
      # gene_list        = NULL   # NULL = core clock genes; or c("Per2", "Nr1d1", ...)
      # cell_group_rules = NULL   # NULL = built-in bins;    or list("T cells" = c(...), ...)
    ),
    error = function(e) message("  Clock acrophase plot skipped: ", e$message)
  )
} else {
  cat("\n  No Step-A results accumulated -- clock acrophase plot skipped.\n")
}

# =============================================================================
# 10. RUN SUMMARY
# =============================================================================

cat(sprintf("\n%s\n=== Pipeline complete: %s stage ===\n%s\n",
            strrep("=", 64), stage, strrep("=", 64)))
print(run_log)

write.csv(run_log,
  file.path(out_dir, sprintf("run_log_%s.csv", stage)),
  row.names = FALSE)
cat(sprintf("\nAll outputs -> %s\n", out_dir))

# =============================================================================
# POST-PROCESSING — Custom acrophase polar plots (run any time after pipeline)
# =============================================================================
# load_stage_results() reads all circadian_analysis_all.csv files from disk so
# you can re-draw the polar plot without re-running the analysis.
# It recovers original cell type names from run_log so cell_group_rules can
# use the same names (with spaces / special characters) as in your data.
#
# Workflow:
#   1. Inspect which genes are present across all cell types
#   2. Choose a custom gene_list  (any confident circadian genes you like)
#   3. Define cell_group_rules    (group cell types into biological categories)
#   4. Call plot_clock_acrophase()
#
# ── Step 1: load ─────────────────────────────────────────────────────────────
# clock_results_post <- load_stage_results(out_dir, period12 = period12)
#
# ── Step 2: inspect available genes ──────────────────────────────────────────
# all_genes <- sort(unique(unlist(lapply(clock_results_post,
#   function(df) df$Genes[df$pvalue < 0.05 & df$pvalue_corr < 0.05]))))
# head(all_genes, 30)   # look for genes present across multiple cell types
#
# ── Step 3 & 4: custom plot ───────────────────────────────────────────────────
# plot_clock_acrophase(
#   results_list     = clock_results_post,
#   stage            = stage,
#   outfile          = file.path(out_dir,
#                        sprintf("clock_acrophase_%s_custom.png", stage)),
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
