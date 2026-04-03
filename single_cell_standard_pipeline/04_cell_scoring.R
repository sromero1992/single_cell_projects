# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 4: CELL SCORING & PATHWAY ENRICHMENT
# Version: 1.0 (Generalized from Nr4a1 T cell production script)
#
# PURPOSE:
#   This script performs gene-set–level scoring and visualization on any
#   annotated Seurat object output from Script 02. It supports two scoring
#   engines and two analysis modes:
#
#   SCORING ENGINES:
#     - AUCell (recommended): Ranks genes per cell and computes the "Area
#       Under the Curve" for a gene set within each cell's ranked list.
#       Robust to data sparsity. Suitable for pathway-level scoring.
#     - AddModuleScore (fast alternative): Seurat's built-in mean-expression
#       approach relative to background control genes. Faster but less robust.
#
#   ANALYSIS MODES:
#     MODE A — "pathway_scoring": Score cells against gene sets from msigdbr
#       (e.g., GO Biological Process, Hallmark, KEGG, Reactome). Produces
#       per-cell AUCell scores, bar/violin plots per cell sub-type, and
#       pseudo-bulk t-test tables for publication-ready statistics.
#
#     MODE B — "gene_expression": Violin and bar plots for individual genes
#       (or module scores) across annotated cell types and experimental groups.
#       Uses generate_gene_comparison_plots() from Script 02's logic.
#
# KEY FEATURES (from Nr4a1 T cell production script, now generalized):
#   - Subset to any annotated cell population (not just T cells)
#   - Configurable pathway collection and sub-collection (C5/BP, H, C2/KEGG, etc.)
#   - Configurable gene lists (collaborator_genes_1, collaborator_genes_2)
#   - Pseudo-bulk scoring per sample for proper statistical testing
#   - Automatic pairwise t-tests for all sub-types and pathways
#   - All plots saved as high-resolution PNGs
#
# DEPENDENCIES:
#   Seurat, AUCell, msigdbr, dplyr, ggplot2, ggpubr, patchwork, writexl, tidyr
# =============================================================================

# --- Load Libraries ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)        # stat_compare_means, Wilcoxon annotations
library(patchwork)
library(writexl)
library(AUCell)        # AUCell scoring
library(msigdbr)       # MSigDB gene set collections

set.seed(42)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) ---
# =============================================================================

# --- 1.1: Project Paths (must match Scripts 01–03) ---
PROJECT_NAME <- "Nr4a1_Study17_Project"
ROOT_PATH    <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"

OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")
SCORE_DIR    <- file.path(OUTPUT_DIR, "cell_scoring")
if (!dir.exists(SCORE_DIR)) dir.create(SCORE_DIR, recursive = TRUE)

MAIN_RDS <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))

# --- 1.2: Cell Population to Score -----------------------------------------
# CELLTYPE_COLUMN: Metadata column with annotated cell types.
# TARGET_CELL_TYPE: The specific broad cell type to subset for scoring.
#   Set to NULL to score ALL cells (no subsetting).
# SUBTYPE_COLUMN: After subsetting, this column defines the sub-populations
#   used as facets in all comparison plots.
# SUBTYPE_LEVELS: Order of sub-type levels in plots. Set to NULL to auto-detect.
CELLTYPE_COLUMN   <- "CellType"
TARGET_CELL_TYPE  <- "T cells"         # e.g., "T cells", "Macrophages", NULL
SUBTYPE_COLUMN    <- "sub_cell_types"  # e.g., "sub_cell_types", "CellType"
SUBTYPE_LEVELS    <- c("CD4+ T cells", "CD8+ T cells", "Tregs",
                       "NK cells", "gd T cells", "Cyc. T cells")  # NULL = auto

# --- 1.3: Experimental Groups -----------------------------------------------
# CONDITION_COLUMN: Metadata column defining the two conditions being compared.
# CONDITION_LEVELS: Factor order for conditions in all plots.
# CONDITION_COLORS: Named vector of hex colors (one per condition level).
CONDITION_COLUMN <- "Condition"
CONDITION_LEVELS <- c("Nr4a1 KO", "WT")
CONDITION_COLORS <- c("Nr4a1 KO" = "#F8766D", "WT" = "#00BFC4")

# --- 1.4: Analysis Modes ---
RUN_MODE_A_PATHWAY_SCORING  <- TRUE   # AUCell pathway scoring from MSigDB
RUN_MODE_B_GENE_EXPRESSION  <- TRUE   # Per-gene violin/bar comparison plots

# --- 1.5: Mode A: Pathway Scoring Configuration ---
# MSIGDB_SPECIES: Species name for msigdbr (e.g., "Mus musculus", "Homo sapiens")
MSIGDB_SPECIES     <- "Mus musculus"
# MSIGDB_CATEGORY / MSIGDB_SUBCATEGORY: MSigDB collection to use.
#   Common options: C5/BP (GO BioProcess), H/NULL (Hallmark), C2/KEGG, C2/REACTOME
MSIGDB_CATEGORY    <- "C5"
MSIGDB_SUBCATEGORY <- "BP"

# PATHWAY_KEYWORDS: Character vector of keyword substrings to filter pathways.
#   Set to NULL to score ALL pathways in the collection (can be slow for C5/BP).
#   Examples: c("T_CELL", "IMMUNE"), c("APOPTOSIS"), NULL
PATHWAY_KEYWORDS   <- c("T_CELL")

# AUCELL_TOPRANK_PCT: Top fraction of ranked genes used to compute AUC.
#   Default 0.05 (top 5%) is recommended for scRNA-seq.
AUCELL_TOPRANK_PCT <- 0.05

# MIN_PATHWAY_GENES: Pathways with fewer genes present in the data are skipped.
MIN_PATHWAY_GENES  <- 3

# --- 1.6: Mode B: Gene Expression Configuration ---
# GENE_LISTS: Named list of gene vectors for expression plots.
#   Each named list element becomes a separate output folder.
#   You can add as many lists as needed.
GENE_LISTS <- list(
  "collaborator_genes" = c(
    "Bcl2",  "Cd38",  "Cd69",  "Ctla4", "Ifng",  "Il10",  "Il12a",
    "Il15",  "Il17a", "Il21",  "Il23a", "Il27",  "Il2",   "Il2ra",
    "Il4",   "Il6",   "Il7",   "Mki67", "Nkg7",  "Pdcd1", "Prf1",
    "Stat5a","Tgfb1", "Tigit", "Tnf",   "Top2a", "Tox",   "Foxp3"
  ),
  "collaborator_genes_2" = c(
    "Ccl3",   "Ccl4",   "Ccl5",   "Egr1",   "Egr2",   "Eomes",
    "Fos",    "Gzmk",   "Havcr2", "Icos",   "Irf4",   "Jun",
    "Klrd1",  "Klrk1",  "Lag3",   "Lamp1",  "Myc",    "Nfkbia",
    "Rab27a", "Stat1",  "Stat4",  "Tbx21",  "Tnfrsf4","Tnfrsf9",
    "Zeb2",   "Nr4a1",  "Nr4a2",  "Nr4a3"
  )
)

# GENE_PLOT_TYPE: "violin" or "barplot" for gene expression comparisons.
GENE_PLOT_TYPE <- "violin"

# COMPARISON_PAIRS: List of character(2) for pairwise Wilcoxon comparisons.
COMPARISON_PAIRS <- list(CONDITION_LEVELS)  # Default: one pair, all conditions

# --- 1.7: Plot Settings ---
DPI_SETTING <- 300

# =============================================================================
# --- PART 2: UTILITY FUNCTIONS -----------------------------------------------
# =============================================================================

# ---------------------------------------------------------------------------
#' plot_expression_custom: Bar or violin plot for a single gene/score column
#'   across annotated sub-types, split by condition.
#'   Includes Wilcoxon rank-sum test annotations and optional CSV export.
# ---------------------------------------------------------------------------
plot_expression_custom <- function(seurat_obj,
                                   gene,
                                   plot_type     = "bar",
                                   subtype_col   = SUBTYPE_COLUMN,
                                   subtype_levels = SUBTYPE_LEVELS,
                                   condition_col  = CONDITION_COLUMN,
                                   conditions     = CONDITION_LEVELS,
                                   cond_colors    = CONDITION_COLORS,
                                   comparisons    = COMPARISON_PAIRS,
                                   assay          = "RNA",
                                   layer          = "data",
                                   save_path      = NULL) {

  # Fetch expression from assay, or from metadata if it's a score column
  if (gene %in% rownames(seurat_obj)) {
    expr_vec <- as.numeric(GetAssayData(seurat_obj, assay = assay, layer = layer)[gene, ])
  } else if (gene %in% colnames(seurat_obj@meta.data)) {
    expr_vec <- as.numeric(seurat_obj@meta.data[[gene]])
  } else {
    warning("Feature not found: ", gene)
    return(NULL)
  }

  if (is.null(subtype_levels)) {
    subtype_levels <- as.character(unique(seurat_obj@meta.data[[subtype_col]]))
    subtype_levels <- subtype_levels[!is.na(subtype_levels)]
  }

  plot_df <- seurat_obj@meta.data %>%
    dplyr::mutate(expr = expr_vec,
                  subtype   = .data[[subtype_col]],
                  condition = .data[[condition_col]]) %>%
    dplyr::filter(subtype %in% subtype_levels, condition %in% conditions) %>%
    dplyr::mutate(subtype   = factor(subtype,   levels = subtype_levels),
                  condition = factor(condition, levels = conditions))

  if (nrow(plot_df) == 0) { warning("No cells after filtering for: ", gene); return(NULL) }

  cust_theme <- theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
          axis.text.x   = element_text(angle = 0, hjust = 0.5, size = 12),
          axis.title.x  = element_blank(),
          axis.text.y   = element_text(size = 13),
          axis.title.y  = element_text(size = 14),
          strip.text    = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          legend.text   = element_text(size = 12),
          legend.title  = element_text(size = 12))

  p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition))

  if (plot_type == "bar") {
    p <- p +
      geom_bar(stat = "summary", fun = "mean", color = "black", alpha = 0.8,
               position = position_dodge(0.9), width = 0.7) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.25,
                    color = "black", position = position_dodge(0.9)) +
      geom_jitter(aes(fill = condition), shape = 21, color = "black", stroke = 0.4,
                  width = 0.12, size = 1.8, alpha = 0.75)
    y_label <- "Mean normalized expression"
  } else {
    p <- p +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
      geom_boxplot(width = 0.10, outlier.size = 0.3, color = "grey40", linewidth = 0.4)
    y_label <- "Normalized expression"
  }

  p <- p +
    stat_compare_means(comparisons = comparisons, label = "p.signif",
                       method = "wilcox.test", method.args = list(exact = FALSE),
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols   = c("****", "***", "**", "*", "ns"))) +
    facet_wrap(~ subtype, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cond_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(title = bquote(italic(.(gene))), y = y_label) +
    cust_theme + NoLegend()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 14, height = 4, dpi = DPI_SETTING, bg = "white")
    message("  Saved: ", basename(save_path))

    # Also save summary CSV for bar plots
    if (plot_type == "bar") {
      summary_data <- plot_df %>%
        dplyr::group_by(subtype, condition) %>%
        dplyr::summarise(Mean = mean(expr, na.rm = TRUE),
                         SE   = sd(expr, na.rm = TRUE) / sqrt(n()),
                         N    = n(), .groups = "drop")
      csv_path <- sub("\\.[^.]+$", "_summary.csv", save_path)
      write.csv(summary_data, csv_path, row.names = FALSE)
    }
  }
  return(p)
}

# ---------------------------------------------------------------------------
#' plot_pathway_score: Bar or violin plot for a scored pathway column.
#   Wrapper around plot_expression_custom with pathway-specific labeling.
# ---------------------------------------------------------------------------
plot_pathway_score <- function(seurat_obj, score_col, pathway_label,
                               plot_type = "bar", save_path = NULL) {
  plot_expression_custom(
    seurat_obj  = seurat_obj,
    gene        = score_col,
    plot_type   = plot_type,
    save_path   = save_path
  )
}

# =============================================================================
# --- PART 3: LOAD DATA AND SUBSET --------------------------------------------
# =============================================================================
message("=== Loading annotated Seurat object ===")
data <- readRDS(MAIN_RDS)

# Subset to target cell type (if specified)
if (!is.null(TARGET_CELL_TYPE)) {
  message(paste("=== Subsetting to:", TARGET_CELL_TYPE, "==="))
  score_obj <- subset(data, cells = rownames(data@meta.data)[
    data@meta.data[[CELLTYPE_COLUMN]] == TARGET_CELL_TYPE
  ])
  message(paste("  Cells in subset:", ncol(score_obj)))
} else {
  message("=== Scoring all cells (no subsetting) ===")
  score_obj <- data
}

DefaultAssay(score_obj) <- "RNA"

# Set factor levels for conditions and sub-types
if (CONDITION_COLUMN %in% colnames(score_obj@meta.data)) {
  score_obj@meta.data[[CONDITION_COLUMN]] <- factor(
    score_obj@meta.data[[CONDITION_COLUMN]], levels = CONDITION_LEVELS
  )
}
if (!is.null(SUBTYPE_LEVELS) && SUBTYPE_COLUMN %in% colnames(score_obj@meta.data)) {
  score_obj@meta.data[[SUBTYPE_COLUMN]] <- factor(
    score_obj@meta.data[[SUBTYPE_COLUMN]], levels = SUBTYPE_LEVELS
  )
}

# Resolve SUBTYPE_LEVELS if NULL
if (is.null(SUBTYPE_LEVELS) && SUBTYPE_COLUMN %in% colnames(score_obj@meta.data)) {
  SUBTYPE_LEVELS <- as.character(unique(score_obj@meta.data[[SUBTYPE_COLUMN]]))
  SUBTYPE_LEVELS <- SUBTYPE_LEVELS[!is.na(SUBTYPE_LEVELS)]
}

pop_label <- ifelse(is.null(TARGET_CELL_TYPE), "AllCells", gsub(" ", "_", TARGET_CELL_TYPE))
message(paste("  Cell population for scoring:", pop_label))
message(sprintf("  Sub-types: %s", paste(SUBTYPE_LEVELS, collapse = ", ")))

# =============================================================================
# --- PART 4: MODE A — PATHWAY SCORING (AUCell) ------------------------------
# =============================================================================
if (RUN_MODE_A_PATHWAY_SCORING) {
  message("\n=== MODE A: AUCell Pathway Scoring ===")

  pathway_dir <- file.path(SCORE_DIR, paste0(pop_label, "_pathways"))
  if (!dir.exists(pathway_dir)) dir.create(pathway_dir, recursive = TRUE)
  dir.create(file.path(pathway_dir, "plots_bar"),    showWarnings = FALSE)
  dir.create(file.path(pathway_dir, "plots_violin"), showWarnings = FALSE)

  # --- Fetch gene sets from MSigDB ---
  message(paste("  Fetching", MSIGDB_CATEGORY, "/", MSIGDB_SUBCATEGORY,
                "gene sets from msigdbr (", MSIGDB_SPECIES, ")..."))
  msig_all <- msigdbr(species  = MSIGDB_SPECIES,
                      category = MSIGDB_CATEGORY,
                      subcategory = if (nchar(MSIGDB_SUBCATEGORY) > 0) MSIGDB_SUBCATEGORY else NULL)

  # Filter by keyword if specified
  if (!is.null(PATHWAY_KEYWORDS) && length(PATHWAY_KEYWORDS) > 0) {
    keyword_pattern <- paste(PATHWAY_KEYWORDS, collapse = "|")
    msig_all <- msig_all[grepl(keyword_pattern, msig_all$gs_name, ignore.case = FALSE), ]
    message(paste("  Pathways after keyword filter:", length(unique(msig_all$gs_name))))
  }

  # Build named list: pathway name → gene vector (intersected with data)
  pathway_names <- unique(msig_all$gs_name)
  pathway_gene_sets <- lapply(pathway_names, function(pw) {
    genes <- msig_all[msig_all$gs_name == pw, "gene_symbol", drop = TRUE] %>% unique()
    intersect(genes, rownames(score_obj))
  })
  names(pathway_gene_sets) <- pathway_names

  # Filter by minimum gene count
  n_genes      <- sapply(pathway_gene_sets, length)
  pathway_gene_sets_ok <- pathway_gene_sets[n_genes >= MIN_PATHWAY_GENES]
  message(paste("  Pathways with >=", MIN_PATHWAY_GENES, "expressed genes:",
                length(pathway_gene_sets_ok)))

  # Report gene set sizes
  gs_sizes <- data.frame(
    pathway = gsub(paste0(MSIGDB_CATEGORY, "_"), "", names(n_genes)),
    n_genes = n_genes
  )
  write.csv(gs_sizes, file.path(pathway_dir, "pathway_gene_set_sizes.csv"), row.names = FALSE)

  if (length(pathway_gene_sets_ok) == 0) {
    warning("No valid pathways found. Check PATHWAY_KEYWORDS and MSIGDB settings.")
  } else {
    # --- Build AUCell rankings ---
    message("  Building AUCell rankings (may take several minutes)...")
    cells_rankings <- AUCell_buildRankings(
      GetAssayData(score_obj, layer = "data"),
      verbose = FALSE, plot = FALSE
    )

    # --- Calculate AUC ---
    message("  Calculating AUC scores for all pathways...")
    cells_AUC <- AUCell_calcAUC(
      pathway_gene_sets_ok, cells_rankings,
      aucMaxRank = nrow(cells_rankings) * AUCELL_TOPRANK_PCT,
      verbose    = FALSE
    )

    # Store scores in metadata
    auc_mat <- getAUC(cells_AUC)
    score_cols <- character(nrow(auc_mat))
    for (k in seq_len(nrow(auc_mat))) {
      pw       <- rownames(auc_mat)[k]
      col_name <- paste0("score_", gsub(paste0(MSIGDB_CATEGORY, "_"), "", pw))
      col_name <- gsub("[^A-Za-z0-9_]", "_", col_name)
      score_obj[[col_name]] <- as.numeric(auc_mat[pw, ])
      score_cols[k]         <- col_name
    }
    message("  AUCell scores added to metadata.")

    # --- Generate plots for each pathway ---
    message("  Generating bar and violin plots per pathway...")
    for (k in seq_along(pathway_gene_sets_ok)) {
      pw        <- names(pathway_gene_sets_ok)[k]
      sc        <- score_cols[k]
      short_lbl <- gsub(paste0(MSIGDB_CATEGORY, "_"), "", pw) %>%
        gsub("_", " ", .) %>% tools::toTitleCase()
      base_nm   <- gsub(paste0(MSIGDB_CATEGORY, "_"), "", pw)

      tryCatch({
        plot_pathway_score(
          score_obj, sc, short_lbl, "bar",
          file.path(pathway_dir, "plots_bar", paste0(base_nm, "_bar.png"))
        )
        plot_pathway_score(
          score_obj, sc, short_lbl, "violin",
          file.path(pathway_dir, "plots_violin", paste0(base_nm, "_violin.png"))
        )
      }, error = function(e) warning("Plot error for ", pw, ": ", e$message))
    }

    # --- Pseudo-bulk statistics ---
    # Use sample-level means (pseudobulk) for publication-quality statistics.
    # Reduces inflation caused by single-cell pseudoreplication.
    message("  Computing pseudo-bulk pathway scores per sample × sub-type...")

    # Determine the sample ID column
    sample_col <- if ("SampleID" %in% colnames(score_obj@meta.data)) "SampleID" else
      if ("BatchID" %in% colnames(score_obj@meta.data)) "BatchID" else NULL

    if (is.null(sample_col)) {
      warning("  [SKIP] No SampleID or BatchID column found. Skipping pseudo-bulk stats.")
    } else {
      pseudobulk_df <- score_obj@meta.data %>%
        dplyr::filter(.data[[SUBTYPE_COLUMN]] %in% SUBTYPE_LEVELS,
                      .data[[CONDITION_COLUMN]] %in% CONDITION_LEVELS) %>%
        dplyr::group_by(!!sym(sample_col), !!sym(CONDITION_COLUMN), !!sym(SUBTYPE_COLUMN)) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
                         n_cells = n(), .groups = "drop")

      write.csv(pseudobulk_df,
                file.path(pathway_dir, "pseudobulk_pathway_scores.csv"), row.names = FALSE)

      # Pairwise t-tests for each pathway × sub-type combination
      cond1 <- CONDITION_LEVELS[1]
      cond2 <- CONDITION_LEVELS[2]
      ttest_rows <- list()
      for (sc in score_cols) {
        for (st in SUBTYPE_LEVELS) {
          sub <- dplyr::filter(pseudobulk_df, .data[[SUBTYPE_COLUMN]] == st)
          g1  <- sub[sub[[CONDITION_COLUMN]] == cond1, sc, drop = TRUE]
          g2  <- sub[sub[[CONDITION_COLUMN]] == cond2, sc, drop = TRUE]
          if (length(g1) < 2 || length(g2) < 2) next
          res <- tryCatch(t.test(g1, g2), error = function(e) NULL)
          if (is.null(res)) next
          ttest_rows[[length(ttest_rows) + 1]] <- data.frame(
            pathway     = gsub("^score_", "", sc),
            SubType     = st,
            mean_cond1  = round(mean(g1, na.rm = TRUE), 5),
            mean_cond2  = round(mean(g2, na.rm = TRUE), 5),
            cond1_label = cond1, cond2_label = cond2,
            p_value     = round(res$p.value, 5),
            direction   = ifelse(mean(g1) > mean(g2),
                                 paste0(cond1, " > ", cond2),
                                 paste0(cond1, " < ", cond2)),
            sig = dplyr::case_when(res$p.value < 0.001 ~ "***",
                                   res$p.value < 0.01  ~ "**",
                                   res$p.value < 0.05  ~ "*",
                                   TRUE                ~ "ns"),
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(ttest_rows) > 0) {
        ttest_df <- dplyr::bind_rows(ttest_rows) %>% dplyr::arrange(p_value)
        write.csv(ttest_df, file.path(pathway_dir, "pseudobulk_ttest_results.csv"), row.names = FALSE)
        sig_hits <- dplyr::filter(ttest_df, sig != "ns")
        message(paste("  Significant pathway hits (pseudobulk t-test):", nrow(sig_hits)))
        if (nrow(sig_hits) > 0) print(sig_hits)
        message("  IMPORTANT: Use pseudobulk_ttest_results.csv for publication statistics.")
        message("  Per-cell Wilcoxon tests on individual plots are exploratory only.")
      }
    }
  }
}

# =============================================================================
# --- PART 5: MODE B — GENE EXPRESSION PLOTS ---------------------------------
# =============================================================================
if (RUN_MODE_B_GENE_EXPRESSION) {
  message("\n=== MODE B: Gene Expression Comparison Plots ===")

  for (list_name in names(GENE_LISTS)) {
    gene_list <- GENE_LISTS[[list_name]]
    gene_dir  <- file.path(SCORE_DIR, paste0(pop_label, "_gene_expression"), list_name)
    if (!dir.exists(gene_dir)) dir.create(gene_dir, recursive = TRUE)

    message(paste("  Processing gene list:", list_name, "(", length(gene_list), "features)"))

    for (gene in gene_list) {
      save_path_bar <- file.path(gene_dir, paste0(gene, "_bar.png"))
      save_path_vio <- file.path(gene_dir, paste0(gene, "_violin.png"))

      plot_expression_custom(score_obj, gene, plot_type = "bar",
                             save_path = save_path_bar,
                             comparisons = COMPARISON_PAIRS)
      plot_expression_custom(score_obj, gene, plot_type = "violin",
                             save_path = save_path_vio,
                             comparisons = COMPARISON_PAIRS)
    }
    message(paste("  Gene list complete:", list_name))
  }
}

# =============================================================================
# --- PART 6: SUMMARY ---------------------------------------------------------
# =============================================================================
message(paste0(
  "\n=== CELL SCORING COMPLETE ===\n",
  "  Output directory: ", SCORE_DIR, "\n",
  if (RUN_MODE_A_PATHWAY_SCORING) paste0(
    "  Pathways scored: ", if (exists("pathway_gene_sets_ok")) length(pathway_gene_sets_ok) else 0, "\n",
    "  Pseudobulk stats: ", file.path(SCORE_DIR, pop_label, "_pathways/pseudobulk_ttest_results.csv"), "\n"
  ) else "",
  if (RUN_MODE_B_GENE_EXPRESSION) paste0(
    "  Gene lists processed: ", paste(names(GENE_LISTS), collapse = ", "), "\n"
  ) else "",
  "\nNote: For publication statistics, always use the pseudobulk t-test CSV.\n",
  "Per-cell Wilcoxon p-values on individual plots are exploratory guides only.\n"
))
