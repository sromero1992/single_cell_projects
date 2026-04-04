# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 2: INTERACTIVE ANNOTATION STUDIO
# Version: 3.0 (Generalized, CSV-Driven, Weighted Pre-Scoring)
#
# PURPOSE:
#   This script is the interactive biological annotation hub. It loads the
#   processed Seurat object from Script 01 and provides:
#     1. WEIGHTED PRE-SCORING: Automatically scores every cluster against all
#        broad cell type markers using a specificity-weighted gene expression
#        metric (both z-scored and raw). This gives a quantitative "best guess"
#        annotation BEFORE manual review — the visual and numeric output should
#        be used as a guide, not a final answer.
#     2. MANUAL BROAD ANNOTATION: The biologist fills in BROAD_ANNOTATION_MAP
#        after reviewing the UMAP + DotPlot + weighted score results.
#     3. COMPOSITIONAL ANALYSIS: Cell type proportion plots by SampleID and
#        any additional metadata groups (e.g., Genotype_Diet, Condition).
#     4. SUB-CLUSTERING: High-resolution re-clustering of any user-specified
#        broad cell type, with weighted sub-scoring from the CSV.
#     5. GENE EXPRESSION COMPARISON: Violin/barplot comparisons of any gene
#        across groups with statistical tests (Wilcoxon).
#
# MARKER SYSTEM:
#   All marker genes are loaded from a CSV file (cell_type_markers.csv).
#   This CSV contains broad and sub-type markers for all expected cell types.
#   Sub-type markers are automatically selected based on CELL_TYPE_TO_SUBCLUSTER.
#   Tregs (Foxp3, Il2ra, Ikzf2, Ctla4) are included in the T cell sub-panel.
#
# HOW TO USE:
#   1. Fill in Part 1 (paths and options).
#   2. Run the script block by block in RStudio.
#   3. After each weighted scoring UMAP, review the top5 report table.
#   4. Fill in the annotation maps (ACTION 4 and ACTION 7).
#   5. Run the rest of the script for final plots and saves.
# =============================================================================

# --- Load Required Libraries -------------------------------------------------
library(Seurat)
library(harmony)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(writexl)
library(Matrix)

set.seed(123)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) ---
# =============================================================================

# --- 1.1: Project Identity ---
# These MUST match the values set in Script 01.
PROJECT_NAME <- "Nr4a1_Study17_Project"
ROOT_PATH    <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"  # Windows

# --- 1.2: Path Configuration (Auto-generated; do not edit) ---
OUTPUT_DIR       <- file.path(ROOT_PATH, "seurat_output")
INPUT_RDS        <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
MARKERS_CSV_FILE <- file.path(ROOT_PATH, "cell_type_markers.csv")
# Alternative: if CSV is in the same folder as this script:
# MARKERS_CSV_FILE <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "cell_type_markers.csv")

# --- 1.3: Annotation View Mode ---
# "harmony"    = Use Harmony-corrected clusters (recommended for multi-sample data).
# "no_harmony" = Use standard PCA clusters (useful for single-sample or diagnostic runs).
VIEW_MODE <- "harmony"

# --- 1.4: Sub-Clustering Configuration ---
# RUN_SUBCLUSTERING: If TRUE, runs the full sub-clustering workflow (Step 4).
# CELL_TYPE_TO_SUBCLUSTER: ANY broad cell type name (must match the CSV parent_cell_type).
#   Examples: "Colonocytes", "T cells", "Macrophages", "Fibroblasts"
#   Sub-marker genes and the default resolution are read automatically from the CSV.
RUN_SUBCLUSTERING        <- TRUE
CELL_TYPE_TO_SUBCLUSTER  <- "Colonocytes"
SUBCLUSTER_N_HVG         <- 2000    # HVGs for sub-clustering PCA
SUBCLUSTER_N_PCS         <- 50      # PCs used for sub-cluster graph
SUBCLUSTER_K_NEIGHBORS   <- 15      # k for kNN graph in sub-clustering
SUBCLUSTER_MIN_DIST      <- 0.3     # UMAP min.dist for sub-cluster embedding
# SUBCLUSTER_RESOLUTION: Set to NULL to auto-read from the CSV; override here if needed.
SUBCLUSTER_RESOLUTION    <- NULL

# --- 1.5: Compositional Analysis Groups ---
#
# =============================================================================
# USING METADATA COLUMNS FOR GROUPING (read this section carefully)
# =============================================================================
# All metadata columns from your Excel file (Condition, Sex, Diet, Genotype,
# Timepoint, etc.) are automatically available in the Seurat object after
# Script 01 runs. Any of them can be used here for compositional analysis,
# UMAP coloring, or pre-scoring grouping in this script, and for DE / scoring
# comparisons in Scripts 03 and 04.
#
# ADDITIONAL_GROUPS_TO_PLOT controls which metadata columns produce proportion
# stacked-bar plots (and proportion stats Excel files). Compositional plots are
# ALWAYS generated by SampleID — this list adds extra breakdowns on top.
#
# Simple columns (direct from Excel):
#   ADDITIONAL_GROUPS_TO_PLOT <- c("Condition")          # KO vs WT
#   ADDITIONAL_GROUPS_TO_PLOT <- c("Diet")               # inulin vs cellulose
#   ADDITIONAL_GROUPS_TO_PLOT <- c("Sex")                # M vs F
#
# Concatenated columns (auto-created by underscore pattern):
#   "Genotype_Diet"  → auto-built from Genotype + Diet columns (must exist)
#   "Condition_Sex"  → auto-built from Condition + Sex columns (must exist)
#   Any "ColA_ColB" pair works as long as ColA and ColB exist in the metadata.
#
# Multiple groups at once:
#   ADDITIONAL_GROUPS_TO_PLOT <- c("Condition", "Genotype_Diet", "Sex")
#   This produces one set of proportion plots + Excel per listed column.
#
# Downstream (Scripts 03 and 04): use the same column name in CONDITION_COLUMN.
#   If you use "Genotype_Diet" here, set CONDITION_COLUMN <- "Genotype_Diet"
#   in Scripts 03/04 and GROUP_1/GROUP_2 to the concatenated values
#   (e.g., GROUP_1 = "Nr4a1_KO_inulin", GROUP_2 = "WT_cellulose").
# =============================================================================
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_Diet")

# --- 1.6: Gene Expression Comparison (generate_gene_comparison_plots) ---
# COMPARISON_X_AXIS: Metadata column for x-axis grouping (e.g., "Genotype_Diet").
# COMPARISON_GROUPS: Character vector of the specific levels to compare.
#   CONVENTION: list CONTROL groups before TREATMENT groups within each factor
#   so that control (WT) always appears on the left of the x-axis.
# COMPARISON_PAIRS: List of character(2) vectors for pairwise stat comparisons.
#   CONVENTION: always c("control", "treatment") — control first in each pair.
COMPARISON_X_AXIS  <- "Genotype_Diet"
COMPARISON_GROUPS  <- c("WT_cellulose", "KO_cellulose", "WT_inulin", "KO_inulin")
COMPARISON_PAIRS   <- list(
  c("WT_inulin",    "KO_inulin"),
  c("WT_cellulose", "KO_cellulose"),
  c("WT_cellulose", "WT_inulin")
)

# --- 1.7: Output Resolution ---
DPI_SETTING <- 300

# =============================================================================
# --- PART 2: LOAD DATA & MARKERS --------------------------------------------
# =============================================================================
message("=== Loading processed Seurat object ===")
data <- readRDS(INPUT_RDS)
message(paste("  Loaded:", ncol(data), "cells,", nrow(data), "genes"))

# --- Load and parse markers from CSV ----------------------------------------
message("=== Loading markers from CSV ===")
if (!file.exists(MARKERS_CSV_FILE)) {
  stop(paste("[ERROR] Marker CSV not found:", MARKERS_CSV_FILE,
             "\nPlease place cell_type_markers.csv in your ROOT_PATH."))
}
markers_df <- read.csv(MARKERS_CSV_FILE, stringsAsFactors = FALSE)

# Helper: parse pipe-separated marker string into a character vector
parse_markers <- function(marker_string) {
  strsplit(trimws(marker_string), "\\|")[[1]]
}

# Build BROAD_MARKERS_LIST from tier == "broad"
broad_df          <- markers_df[markers_df$tier == "broad", ]
BROAD_MARKERS_LIST <- setNames(
  lapply(broad_df$markers, parse_markers),
  broad_df$cell_type
)

# Build SUB_MARKERS_LIST for the chosen cell type from tier == "sub"
sub_df <- markers_df[markers_df$tier == "sub" &
                       markers_df$parent_cell_type == CELL_TYPE_TO_SUBCLUSTER, ]
SUB_MARKERS_LIST <- if (nrow(sub_df) > 0) {
  setNames(lapply(sub_df$markers, parse_markers), sub_df$cell_type)
} else {
  list()
}

# Read suggested sub-clustering resolution from CSV (use first row's value)
if (is.null(SUBCLUSTER_RESOLUTION)) {
  SUBCLUSTER_RESOLUTION <- if (nrow(sub_df) > 0 && !is.na(sub_df$subcluster_resolution[1])) {
    as.numeric(sub_df$subcluster_resolution[1])
  } else {
    3.0  # default fallback
  }
}

# Build flat dotplot marker vectors from CSV dotplot_markers column
broad_dotplot_markers <- unique(unlist(lapply(broad_df$dotplot_markers, parse_markers)))
sub_dotplot_markers   <- if (nrow(sub_df) > 0) {
  unique(unlist(lapply(sub_df$dotplot_markers, parse_markers)))
} else {
  unique(unlist(SUB_MARKERS_LIST))
}

message(paste("  Broad cell types loaded:", length(BROAD_MARKERS_LIST)))
message(paste("  Sub-types for", CELL_TYPE_TO_SUBCLUSTER, ":", length(SUB_MARKERS_LIST)))
message(paste("  Sub-clustering resolution:", SUBCLUSTER_RESOLUTION))

# --- Set up clustering/UMAP keys based on VIEW_MODE -------------------------
if (VIEW_MODE == "harmony") {
  cluster_col        <- "clusters_harmony"
  umap_col           <- "umap_harmony"
  plot_title_suffix  <- "(Harmony Corrected)"
} else {
  cluster_col        <- "clusters_none"
  umap_col           <- "umap_none"
  plot_title_suffix  <- "(No Harmony)"
}

# =============================================================================
# --- PART 3: UTILITY FUNCTIONS -----------------------------------------------
# =============================================================================

# ---------------------------------------------------------------------------
#' get_weighted_annotation
#'
#' Scores every cell in a Seurat object against a named list of marker gene
#' sets using a specificity-weighted, mean cluster-level expression approach.
#'
#' ALGORITHM:
#'   1. Build a binary gene × cell-type presence matrix.
#'   2. Weight each gene by 1/occurrence (cross-cell-type specificity weight).
#'   3. Normalize each cell-type column so weights sum to 1.
#'   4. Optionally z-score gene expression across cells (standardize = TRUE).
#'   5. Compute per-cluster mean score for each cell type via matrix multiply.
#'   6. Assign the top-scoring cell type to all cells in each cluster.
#'
#' @param seurat_obj           Seurat object with normalized RNA data.
#' @param marker_genes         Named list: cell_type → character vector of genes.
#' @param cluster_key          Metadata column with cluster IDs.
#' @param standardize_expression If TRUE, z-score expression per gene across cells.
#'                              TRUE emphasizes relative expression; FALSE uses raw values.
#' @return List with: annotation_vector (per-cell character), top5_report (data.frame)
# ---------------------------------------------------------------------------
get_weighted_annotation <- function(seurat_obj,
                                    marker_genes,
                                    cluster_key,
                                    standardize_expression = TRUE) {
  # --- Sync markers with dataset ---
  all_obj_genes  <- rownames(seurat_obj)
  marker_genes   <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
  marker_genes   <- marker_genes[sapply(marker_genes, length) > 0]
  if (length(marker_genes) == 0) stop("No valid marker genes found in dataset.")

  all_marker_genes <- sort(unique(unlist(marker_genes)))
  cell_types       <- names(marker_genes)

  # --- Build specificity-weighted, column-normalized weight matrix ---
  W_binary <- matrix(0, nrow = length(all_marker_genes), ncol = length(cell_types),
                     dimnames = list(all_marker_genes, cell_types))
  for (ct in cell_types) { W_binary[marker_genes[[ct]], ct] <- 1.0 }

  gene_occurrence <- rowSums(W_binary)
  W_specificity   <- W_binary / gene_occurrence              # 1/occurrence weight
  column_sums     <- colSums(W_specificity)
  W_final         <- sweep(W_specificity, 2, column_sums, "/")  # column normalize → sums to 1
  W_final[is.na(W_final)] <- 0

  # --- Pull normalized expression for marker genes ---
  X_subset <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[all_marker_genes, ]

  if (standardize_expression) {
    # Z-score per gene across all cells (makes high-expression cells relative to others)
    X_subset <- t(scale(t(as.matrix(X_subset))))
    X_subset[is.na(X_subset)] <- 0
  } else {
    X_subset <- as.matrix(X_subset)
  }

  # --- Score per cluster ---
  W          <- t(W_final)   # (CellType × Gene)
  clusters   <- seurat_obj[[cluster_key, drop = TRUE]]
  unique_cls <- unique(clusters)

  annotation_vector  <- character(length = ncol(seurat_obj))
  names(annotation_vector) <- colnames(seurat_obj)
  cluster_score_data <- list()

  for (cluster_id in unique_cls) {
    cell_idx <- which(clusters == cluster_id)
    if (length(cell_idx) == 0) next

    # (CellType × Gene) @ (Gene × Cells_in_cluster) → mean per cell type
    raw_scores  <- W %*% X_subset[, cell_idx, drop = FALSE]
    mean_scores <- rowMeans(raw_scores)

    # Sort descending and keep top 5
    top_idx  <- order(mean_scores, decreasing = TRUE)
    top_ct   <- names(mean_scores)[top_idx]
    top_sc   <- mean_scores[top_idx]

    annotation_vector[cell_idx] <- top_ct[1]

    n_top <- min(5, length(top_ct))
    cluster_score_data[[as.character(cluster_id)]] <- data.frame(
      cluster   = cluster_id,
      Rank      = 1:n_top,
      Cell_Type = top_ct[1:n_top],
      Score     = top_sc[1:n_top]
    )
  }

  top5_report <- do.call(rbind, cluster_score_data)
  return(list(annotation_vector = annotation_vector, top5_report = top5_report))
}

# ---------------------------------------------------------------------------
#' summarize_annotations: Print and return a frequency table for an annotation column.
# ---------------------------------------------------------------------------
summarize_annotations <- function(seurat_obj, annotation_column, print_summary = TRUE) {
  if (!annotation_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", annotation_column, "not found in metadata."))
  }
  counts  <- table(seurat_obj[[annotation_column]])
  summary_df <- data.frame(
    Count      = as.numeric(counts),
    Percentage = as.numeric(prop.table(counts) * 100),
    row.names  = names(counts)
  ) %>% arrange(desc(Percentage))
  if (print_summary) {
    print_df <- summary_df
    print_df$Percentage <- sprintf("%.2f%%", print_df$Percentage)
    message(paste("\n--- Annotation Summary:", annotation_column, "---"))
    print(print_df)
  }
  return(summary_df)
}

# ---------------------------------------------------------------------------
#' plot_cell_proportions: Stacked bar charts of cell type composition.
#'   Generates and saves three plots: by SampleID, by a group column, and global.
#'   Also exports a multi-sheet Excel file with the underlying statistics.
#'
#' @param seurat_obj  Seurat object
#' @param cluster_col Metadata column containing cell type labels
#' @param group_col   Metadata column for the group-level breakdown (e.g., "Genotype_Diet")
#' @param output_prefix  Prefix for saved files
#' @param output_dir  Output directory path
# ---------------------------------------------------------------------------
plot_cell_proportions <- function(seurat_obj, cluster_col, group_col, output_prefix, output_dir) {

  meta <- seurat_obj@meta.data
  make_label <- function(p) { ifelse(p > 2.5, paste0(round(p, 1), "%"), "") }

  # Attempt to auto-create combined group column if not present
  if (!group_col %in% colnames(meta)) {
    parts <- strsplit(group_col, "_")[[1]]
    if (all(parts %in% colnames(meta))) {
      message(paste("  Auto-creating combined column:", group_col))
      seurat_obj[[group_col]] <- apply(meta[, parts, drop = FALSE], 1, paste, collapse = "_")
      meta <- seurat_obj@meta.data
    } else {
      warning(paste("  Skipping proportion plot: column '", group_col, "' not found."))
      return(NULL)
    }
  }

  df_sample <- meta %>%
    group_by(SampleID, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(SampleID) %>%
    mutate(percentage = (n / sum(n)) * 100)

  df_group <- meta %>%
    group_by(!!sym(group_col), !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!sym(group_col)) %>%
    mutate(percentage = (n / sum(n)) * 100)

  df_global <- meta %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(percentage = (n / sum(n)) * 100)

  # Export to Excel
  excel_path <- file.path(output_dir, paste0(output_prefix, "_Stats.xlsx"))
  write_xlsx(list(By_Sample = df_sample, By_Group = df_group, Global = df_global),
             path = excel_path)

  # Custom theme
  prop_theme <- theme_classic() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          axis.title.x = element_blank(),
          strip.text   = element_text(size = 12, face = "bold"),
          legend.position = "bottom")

  p1 <- ggplot(df_sample, aes(x = SampleID, y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = paste("Cell Proportions by Sample:", cluster_col), y = "Percentage (%)") + prop_theme

  p2 <- ggplot(df_group, aes(x = !!sym(group_col), y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = paste("Cell Proportions by", group_col, ":", cluster_col), y = "Percentage (%)") + prop_theme

  p3 <- ggplot(df_global, aes(x = "Global", y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", width = 0.6) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_stack(vjust = 0.5), size = 4.5) +
    labs(title = paste("Global Distribution:", cluster_col), y = "Percentage (%)") + prop_theme

  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_by_Sample.png")), p1,
         width = max(8, 1.5 * n_distinct(meta$SampleID)), height = 7, bg = "white", dpi = DPI_SETTING)
  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_by_", group_col, ".png")), p2,
         width = max(8, 1.5 * n_distinct(meta[[group_col]])), height = 7, bg = "white", dpi = DPI_SETTING)
  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_Global.png")), p3,
         width = 7, height = 7, bg = "white", dpi = DPI_SETTING)

  message(paste("  Proportion plots saved for:", cluster_col, "by", group_col))
  return(list(plots = list(p1, p2, p3), stats = list(By_Sample = df_sample, By_Group = df_group, Global = df_global)))
}

# ---------------------------------------------------------------------------
#' generate_gene_comparison_plots: Violin or bar plot for a gene/score column
#'   across cell types (facets) with group-level x-axis and Wilcoxon tests.
#'
#' @param seurat_obj    Seurat object
#' @param score_col     Gene name or metadata score column
#' @param group_by      Metadata column for faceting (e.g., "broad_cell_types")
#' @param x_axis        Metadata column for x-axis grouping (e.g., "Genotype_Diet")
#' @param comparisons   List of character(2) for pairwise stat comparisons
#' @param plot_type     "violin" or "barplot"
#' @param output_prefix File prefix for saved PNG
#' @param plot_title    Main plot title
#' @param y_label       Y-axis label
#' @param fig_width     Plot width in inches
#' @param fig_height    Plot height in inches
#' @param output_dir    Directory where plot is saved (default: OUTPUT_DIR)
# ---------------------------------------------------------------------------
generate_gene_comparison_plots <- function(seurat_obj,
                                           score_col,
                                           group_by,
                                           x_axis,
                                           comparisons,
                                           plot_type    = "violin",
                                           output_prefix = "",
                                           plot_title   = score_col,
                                           y_label      = "Expression",
                                           fig_width    = 16,
                                           fig_height   = 7,
                                           output_dir   = OUTPUT_DIR) {
  message(paste("  Generating comparison plot for:", score_col))

  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis)) %>%
    dplyr::rename(Expression = 1) %>%
    drop_na()

  cust_theme <- theme_classic() +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 18, face = "bold"),
      strip.text    = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      axis.title.y  = element_text(size = 16, face = "bold"),
      axis.title.x  = element_blank(),
      axis.text.y   = element_text(size = 13),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 13, face = "bold"),
      legend.position = "bottom",
      panel.spacing = unit(1.5, "lines")
    )

  p <- ggplot(df_plot, aes(x = !!sym(x_axis), y = Expression, fill = !!sym(x_axis)))

  if (plot_type == "barplot") {
    p <- p +
      stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.8, linewidth = 0.8) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8)
  } else {
    p <- p +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5, linewidth = 0.6)
  }

  p <- p +
    stat_compare_means(
      comparisons = comparisons, label = "p.signif", method = "wilcox.test",
      method.args = list(exact = FALSE),
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                         symbols   = c("****", "***", "**", "*", "ns")),
      step.increase = 0.1, size = 6, bracket.size = 0.8
    ) +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(title = plot_title, y = y_label) +
    scale_fill_brewer(palette = "Set1") +
    cust_theme

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plot_file <- file.path(output_dir, paste0(output_prefix, score_col, ".png"))
  ggsave(plot_file, p, width = fig_width, height = fig_height, dpi = DPI_SETTING, bg = "white")
  message(paste("  Saved:", basename(plot_file)))
  return(p)
}

# ---------------------------------------------------------------------------
#' process_and_extract_cell_type: Subset, re-cluster, and apply Harmony to a
#'   single broad cell type for high-resolution sub-annotation.
#'
#' @param data           Main annotated Seurat object
#' @param cell_type_name Cell type label to extract (matched against CellType column)
#' @param num_hvg        Number of HVGs for sub-clustering PCA
#' @param dims_pca       Number of PCs for graph construction
#' @param resolution     Leiden clustering resolution
#' @param min_dist       UMAP min.dist
#' @param kneigh         k-nearest neighbors for graph
# ---------------------------------------------------------------------------
process_and_extract_cell_type <- function(data,
                                          cell_type_name,
                                          num_hvg    = 2000,
                                          dims_pca   = 50,
                                          resolution = 3.0,
                                          min_dist   = 0.3,
                                          kneigh     = 15) {
  message(paste("  Subsetting and re-clustering:", cell_type_name))
  data_sub <- subset(data, subset = CellType == cell_type_name)

  # Remove old reductions / graphs — they are stale in the subset context
  data_sub@reductions <- list()
  data_sub@graphs     <- list()

  # Find sub-population-specific HVGs and re-run PCA
  data_sub <- FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = num_hvg) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = dims_pca, reduction.name = "pca", verbose = FALSE)
  gc()

  # --- Track A: Standard PCA ---
  data_sub <- FindNeighbors(data_sub, dims = 1:dims_pca, reduction = "pca",
                             k.param = kneigh, graph.name = "pca_nn", verbose = FALSE) %>%
    FindClusters(resolution = resolution, graph.name = "pca_nn",
                 cluster.name = "clusters_none", verbose = FALSE) %>%
    RunUMAP(dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh,
            min.dist = min_dist, n.epochs = 500,
            reduction.name = "umap_none", verbose = FALSE)
  gc()

  # --- Track B: Harmony ---
  message("  Running Harmony for sub-clustering...")
  data_sub <- RunHarmony(data_sub, group.by.vars = "SampleID",
                          reduction = "pca", reduction.save = "harmony", verbose = FALSE)
  gc()
  data_sub <- FindNeighbors(data_sub, dims = 1:dims_pca, reduction = "harmony",
                             k.param = kneigh, graph.name = "harmony_nn", verbose = FALSE) %>%
    FindClusters(resolution = resolution, graph.name = "harmony_nn",
                 cluster.name = "clusters_harmony", verbose = FALSE) %>%
    RunUMAP(dims = 1:dims_pca, reduction = "harmony", n.neighbors = kneigh,
            min.dist = min_dist, n.epochs = 500,
            reduction.name = "umap_harmony", verbose = FALSE)
  gc()
  return(data_sub)
}

# =============================================================================
# --- PART 3: BROAD ANNOTATION WORKFLOW ----------------------------------------
# =============================================================================
message("\n=== STEP 3: Broad Annotation Workflow ===")

# --- 3.1: Initial UMAP and DotPlot ---
p_umap_clusters <- DimPlot(data, reduction = umap_col, group.by = cluster_col, label = TRUE) +
  ggtitle(paste("Clusters", plot_title_suffix)) + NoLegend()
print(p_umap_clusters)
ggsave(file.path(OUTPUT_DIR, paste0("ANNOTATE_01_UMAP_clusters_", VIEW_MODE, ".png")),
       p_umap_clusters, width = 8, height = 7, dpi = DPI_SETTING)

p_dotplot_broad <- DotPlot(data, features = broad_dotplot_markers,
                           group.by = cluster_col, dot.min = 0.05,
                           cols = "RdBu", scale = TRUE) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Broad Markers", plot_title_suffix))
print(p_dotplot_broad)
ggsave(file.path(OUTPUT_DIR, paste0("ANNOTATE_01_DotPlot_broad_markers_", VIEW_MODE, ".png")),
       p_dotplot_broad, width = 14, height = 12, dpi = DPI_SETTING, bg = "white")

# --- 3.2: WEIGHTED PRE-SCORING (standardized = z-scored) ---------------------
# This is the "badass pre-scoring" — a quantitative first guess before manual annotation.
# standardize_expression = TRUE : z-scores expression per gene, emphasizing relative patterns.
# standardize_expression = FALSE: raw normalized counts, proportional to absolute expression.
message("  Running weighted pre-scoring (standardized)...")
results_std <- get_weighted_annotation(
  seurat_obj             = data,
  marker_genes           = BROAD_MARKERS_LIST,
  cluster_key            = cluster_col,
  standardize_expression = TRUE
)
data$broad_weighted_std <- results_std$annotation_vector

message("  Running weighted pre-scoring (non-standardized / raw)...")
results_raw <- get_weighted_annotation(
  seurat_obj             = data,
  marker_genes           = BROAD_MARKERS_LIST,
  cluster_key            = cluster_col,
  standardize_expression = FALSE
)
data$broad_weighted_raw <- results_raw$annotation_vector

# --- 3.3: Display pre-scoring results ----------------------------------------
# Review these to guide your BROAD_ANNOTATION_MAP below.
# The "Top 5 Report" shows score rankings per cluster — this is your key reference.
message("\n--- Pre-scoring Top-5 Report (standardized) ---")
print(results_std$top5_report)
write_xlsx(results_std$top5_report, file.path(OUTPUT_DIR, "PRESCORE_top5_report_standardized.xlsx"))

message("\n--- Pre-scoring Top-5 Report (non-standardized) ---")
print(results_raw$top5_report)
write_xlsx(results_raw$top5_report, file.path(OUTPUT_DIR, "PRESCORE_top5_report_raw.xlsx"))

# Summarize overall cell type distribution from pre-scoring
summarize_annotations(data, "broad_weighted_std")
summarize_annotations(data, "broad_weighted_raw")

# UMAP comparison: raw clusters vs. both weighted annotations
p_prescore_std <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_std",
                          label = TRUE, repel = TRUE) + ggtitle("Pre-Score: Standardized (z-scored)")
p_prescore_raw <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_raw",
                          label = TRUE, repel = TRUE) + ggtitle("Pre-Score: Non-standardized (raw)")
p_prescore_comp <- p_prescore_std | p_prescore_raw
print(p_prescore_comp)
ggsave(file.path(OUTPUT_DIR, paste0("ANNOTATE_02_UMAP_prescore_comparison_", VIEW_MODE, ".png")),
       p_prescore_comp, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

# =============================================================================
# --- ACTION 4: FILL IN YOUR FINAL BROAD ANNOTATIONS HERE ---
# =============================================================================
# Review: (1) ANNOTATE_01_UMAP_clusters, (2) ANNOTATE_01_DotPlot_broad_markers,
#         (3) ANNOTATE_02_UMAP_prescore_comparison, (4) PRESCORE_top5_report xlsx
#
# Map each cluster number to a cell type string.
# Example (fill in your own mappings):
BROAD_ANNOTATION_MAP <- c(
  '0'  = 'Colonocytes',  '1'  = 'Colonocytes',  '2'  = 'Colonocytes',
  '3'  = 'Colonocytes',  '4'  = 'Colonocytes',  '5'  = 'Colonocytes',
  '6'  = 'Colonocytes',  '7'  = 'Colonocytes',  '8'  = 'Plasma B cells',
  '9'  = 'Colonocytes',  '10' = 'Colonocytes',  '11' = 'B cells',
  '12' = 'Colonocytes',  '13' = 'Colonocytes',  '14' = 'T cells',
  '15' = 'Macrophages',  '16' = 'Fibroblasts',  '17' = 'SMCs',
  '18' = 'Colonocytes',  '19' = 'T cells',      '20' = 'SMCs',
  '21' = 'Fibroblasts',  '22' = 'VECs',         '23' = 'cDCs',
  '24' = 'Colonocytes',  '25' = 'SMCs',         '26' = 'LECs',
  '27' = 'Plasma B cells','28' = 'Colonocytes', '29' = 'SMCs',
  '30' = 'Plasma B cells','31' = 'SMCs',        '32' = 'SMCs',
  '33' = 'Colonocytes',  '34' = 'Colonocytes',  '35' = 'Colonocytes',
  '36' = 'Colonocytes',  '37' = 'Colonocytes',  '38' = 'EGCs',
  '39' = 'B cells',      '40' = 'B cells',      '41' = 'Fibroblasts',
  '42' = 'ENs',          '43' = 'Colonocytes',  '44' = 'Adipocytes',
  '45' = 'LECs',         '46' = 'Colonocytes',  '47' = 'Colonocytes',
  '48' = 'Colonocytes',  '49' = 'Colonocytes',  '50' = 'SMCs',
  '51' = 'LECs',         '52' = 'T cells',      '53' = 'Colonocytes'
)

# =============================================================================
# --- PART 4: APPLY BROAD ANNOTATIONS AND SAVE RESULTS -----------------------
# =============================================================================
message("\n=== STEP 4: Applying Broad Annotations ===")

# Apply the manual annotation map
data$broad_cell_types <- recode_factor(data[[cluster_col, drop = TRUE]], !!!BROAD_ANNOTATION_MAP)
data$CellType         <- data$broad_cell_types  # Alias used by Script 03 and 04

# Set factor level order for consistent color assignment in all plots
broad_type_levels <- unique(c(broad_df$cell_type, as.character(data$broad_cell_types)))
broad_type_levels <- intersect(broad_type_levels, as.character(unique(data$broad_cell_types)))
data$broad_cell_types         <- factor(data$broad_cell_types, levels = broad_type_levels)
data$broad_weighted_std       <- factor(data$broad_weighted_std, levels = broad_type_levels)
data$broad_weighted_raw       <- factor(data$broad_weighted_raw, levels = broad_type_levels)

# --- Final annotation UMAPs: manual + weighted side-by-side ---
p_final_manual <- DimPlot(data, reduction = umap_col, group.by = "broad_cell_types",
                          label = TRUE, repel = TRUE) +
  ggtitle(paste("Final Manual Annotation", plot_title_suffix))
p_final_std <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_std",
                       label = TRUE, repel = TRUE) + ggtitle("Weighted Pre-Score (z-scored)")
p_final_raw <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_raw",
                       label = TRUE, repel = TRUE) + ggtitle("Weighted Pre-Score (raw)")

p_triple <- p_final_manual | p_final_std | p_final_raw
ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotation_comparison.png"),
       p_triple, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")
ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotated.png"),
       p_final_manual, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# --- Final DotPlot ---
p_final_dot <- DotPlot(data, features = broad_dotplot_markers,
                       group.by = "broad_cell_types", dot.min = 0.05,
                       cols = "RdBu", scale = TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11))
ggsave(file.path(OUTPUT_DIR, "FINAL_DotPlot_broad_annotated.png"),
       p_final_dot, width = 8, height = 14, dpi = DPI_SETTING, bg = "white")

message("--- Broad annotation applied and plots saved ---")

# =============================================================================
# --- PART 5: COMPOSITIONAL ANALYSIS ------------------------------------------
# =============================================================================
message("\n=== STEP 5: Compositional Analysis ===")

# Always generate by SampleID
plot_cell_proportions(data, "broad_cell_types", "SampleID",
                      "FINAL_broad_proportions", OUTPUT_DIR)

# Generate for each additional user-specified group
for (group in ADDITIONAL_GROUPS_TO_PLOT) {
  plot_cell_proportions(data, "broad_cell_types", group,
                        paste0("FINAL_broad_proportions_by_", group), OUTPUT_DIR)
}

# --- Faceted UMAP split by primary group (if ADDITIONAL_GROUPS_TO_PLOT is set) ---
if (length(ADDITIONAL_GROUPS_TO_PLOT) > 0) {
  primary_group <- ADDITIONAL_GROUPS_TO_PLOT[1]
  if (primary_group %in% colnames(data@meta.data)) {
    n_levels <- n_distinct(data@meta.data[[primary_group]])
    p_facet <- DimPlot(data, reduction = umap_col, group.by = "broad_cell_types",
                       split.by = primary_group, label = FALSE, repel = TRUE) +
      facet_wrap(as.formula(paste("~", primary_group)), ncol = 2) +
      theme(strip.text = element_text(size = 13, face = "bold"),
            legend.text = element_text(size = 12)) +
      guides(color = guide_legend(override.aes = list(size = 4)))
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_UMAP_faceted_by_", primary_group, ".png")),
           p_facet, width = 12, height = 5 * ceiling(n_levels / 2), dpi = DPI_SETTING, bg = "white")
  }
}
message("--- Compositional analysis complete ---")

# --- Save intermediate object ---
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_broad_annotated.rds")))

# =============================================================================
# --- PART 6: SUB-CLUSTERING WORKFLOW -----------------------------------------
# =============================================================================
if (RUN_SUBCLUSTERING) {
  message(paste("\n=== STEP 6: Sub-Clustering of:", CELL_TYPE_TO_SUBCLUSTER, "==="))

  if (length(SUB_MARKERS_LIST) == 0) {
    warning(paste("[SKIP] No sub-markers found in CSV for:", CELL_TYPE_TO_SUBCLUSTER,
                  "\nAdd sub-tier rows with parent_cell_type =", CELL_TYPE_TO_SUBCLUSTER,
                  "to cell_type_markers.csv"))
  } else {

    # --- 6.1: Extract and re-cluster the target cell type ---
    data_sub <- process_and_extract_cell_type(
      data           = data,
      cell_type_name = CELL_TYPE_TO_SUBCLUSTER,
      num_hvg        = SUBCLUSTER_N_HVG,
      dims_pca       = SUBCLUSTER_N_PCS,
      resolution     = SUBCLUSTER_RESOLUTION,
      min_dist       = SUBCLUSTER_MIN_DIST,
      kneigh         = SUBCLUSTER_K_NEIGHBORS
    )

    # --- 6.2: Diagnostic plots ---
    p_sub_h_vs_pca <- (
      DimPlot(data_sub, reduction = "umap_none",    group.by = "SampleID") +
        ggtitle("Sub PCA (No Harmony)")
    ) + (
      DimPlot(data_sub, reduction = "umap_harmony", group.by = "SampleID") +
        ggtitle("Sub Harmony")
    )
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_01_Harmony_comparison_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_h_vs_pca, width = 16, height = 7, dpi = DPI_SETTING)

    p_sub_clust <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony",
                           label = TRUE) + NoLegend() +
      ggtitle(paste("Sub-Clusters:", CELL_TYPE_TO_SUBCLUSTER))
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_02_ANNOTATE_THIS_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_clust, width = 9, height = 8, dpi = DPI_SETTING)

    p_sub_dot <- DotPlot(data_sub, features = sub_dotplot_markers,
                         group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu") +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
      ggtitle(paste("Sub-Cluster Markers:", CELL_TYPE_TO_SUBCLUSTER))
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_dot, width = 14, height = 9, dpi = DPI_SETTING, bg = "white")

    # --- 6.3: Weighted pre-scoring for sub-types ---
    message("  Running sub-type weighted pre-scoring (standardized)...")
    sub_results_std <- get_weighted_annotation(
      seurat_obj             = data_sub,
      marker_genes           = SUB_MARKERS_LIST,
      cluster_key            = "clusters_harmony",
      standardize_expression = TRUE
    )
    data_sub$sub_weighted_std <- sub_results_std$annotation_vector

    message("  Running sub-type weighted pre-scoring (non-standardized)...")
    sub_results_raw <- get_weighted_annotation(
      seurat_obj             = data_sub,
      marker_genes           = SUB_MARKERS_LIST,
      cluster_key            = "clusters_harmony",
      standardize_expression = FALSE
    )
    data_sub$sub_weighted_raw <- sub_results_raw$annotation_vector

    message("\n--- Sub-cluster Pre-scoring Top-5 Report (standardized) ---")
    print(sub_results_std$top5_report)
    write_xlsx(sub_results_std$top5_report,
               file.path(OUTPUT_DIR, paste0("SUBCLUSTER_PRESCORE_top5_", CELL_TYPE_TO_SUBCLUSTER, "_std.xlsx")))
    write_xlsx(sub_results_raw$top5_report,
               file.path(OUTPUT_DIR, paste0("SUBCLUSTER_PRESCORE_top5_", CELL_TYPE_TO_SUBCLUSTER, "_raw.xlsx")))

    p_sub_prescore <- (
      DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std",
              label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Standardized")
    ) | (
      DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw",
              label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Raw")
    )
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_03_prescore_comparison_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_prescore, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

    # =========================================================================
    # --- ACTION 7: FILL IN SUB-CLUSTER ANNOTATIONS HERE ---
    # =========================================================================
    # Review: SUBCLUSTER_02_ANNOTATE_THIS_UMAP, SUBCLUSTER_02_DOTPLOT,
    #         SUBCLUSTER_03_prescore_comparison, and the top5 Excel report.
    SUB_ANNOTATION_MAP <- c(
      '0'  = 'Goblet cells', '1'  = 'Abs. colonocytes', '2'  = 'Abs. colonocytes',
      '3'  = 'Goblet cells', '4'  = 'Stem cells',       '5'  = 'Goblet cells',
      '6'  = 'Abs. colonocytes', '7' = 'Goblet cells',  '8'  = 'Abs. colonocytes',
      '9'  = 'Abs. colonocytes', '10' = 'TA cells',     '11' = 'Abs. colonocytes',
      '12' = 'Abs. colonocytes', '13' = 'TA cells',     '14' = 'TA cells',
      '15' = 'Abs. colonocytes', '16' = 'Abs. colonocytes', '17' = 'TA cells',
      '18' = 'Stem cells', '19' = 'Abs. colonocytes',   '20' = 'Goblet cells',
      '21' = 'Abs. colonocytes', '22' = 'Goblet cells', '23' = 'Abs. colonocytes',
      '24' = 'Goblet cells', '25' = 'Goblet cells',     '26' = 'Goblet cells',
      '27' = 'Abs. colonocytes', '28' = 'Abs. colonocytes', '29' = 'TA cells',
      '30' = 'Abs. colonocytes', '31' = 'Abs. colonocytes', '32' = 'Stem cells',
      '33' = 'Abs. colonocytes', '34' = 'Goblet cells', '35' = 'Abs. colonocytes',
      '36' = 'TA cells', '37' = 'Abs. colonocytes',     '38' = 'TA cells',
      '39' = 'Abs. colonocytes', '40' = 'EECs',         '41' = 'Goblet cells',
      '42' = 'TA cells', '43' = 'Goblet cells',         '44' = 'Goblet cells',
      '45' = 'Tuft cells', '46' = 'Tuft cells',         '47' = 'EECs',
      '48' = 'Goblet cells', '49' = 'Abs. colonocytes', '50' = 'Abs. colonocytes',
      '51' = 'Goblet cells', '52' = 'Abs. colonocytes', '53' = 'Goblet cells',
      '54' = 'Goblet cells', '55' = 'Goblet cells',     '56' = 'Abs. colonocytes',
      '57' = 'Goblet cells', '58' = 'Stem cells',       '59' = 'Goblet cells',
      '60' = 'Abs. colonocytes', '61' = 'TA cells',     '62' = 'Abs. colonocytes'
    )

    # Apply sub annotations
    data_sub$sub_cell_types <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
    data_sub$CellType        <- data_sub$sub_cell_types
    data_sub$seurat_clusters <- data_sub$clusters_harmony

    # Set factor levels for consistent coloring
    sub_type_levels <- names(SUB_MARKERS_LIST)
    sub_type_levels <- intersect(sub_type_levels, as.character(unique(data_sub$sub_cell_types)))
    data_sub$sub_cell_types   <- factor(data_sub$sub_cell_types, levels = sub_type_levels)
    data_sub$sub_weighted_std <- factor(data_sub$sub_weighted_std, levels = sub_type_levels)
    data_sub$sub_weighted_raw <- factor(data_sub$sub_weighted_raw, levels = sub_type_levels)

    # --- Final sub-type plots ---
    p_sub_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                           label = TRUE, repel = TRUE) +
      ggtitle(paste("Sub-Cell Types:", CELL_TYPE_TO_SUBCLUSTER))
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_SUBCLUSTER_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_final, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

    p_sub_compare <- DimPlot(data_sub, reduction = "umap_harmony",
                             group.by = c("sub_cell_types", "sub_weighted_std", "sub_weighted_raw"),
                             label = FALSE, repel = TRUE) +
      ggtitle("Sub-type: Manual vs. Weighted (std vs raw)")
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_SUBCLUSTER_comparison_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_compare, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")

    p_sub_dot_final <- DotPlot(data_sub, features = sub_dotplot_markers,
                               group.by = "sub_cell_types", dot.min = 0.05,
                               cols = "RdBu", scale = TRUE) + coord_flip() +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11))
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_DotPlot_sub_", CELL_TYPE_TO_SUBCLUSTER, ".png")),
           p_sub_dot_final, width = 8, height = 12, dpi = DPI_SETTING, bg = "white")

    # --- Sub-type proportional analysis ---
    plot_cell_proportions(data_sub, "sub_cell_types", "SampleID",
                          paste0("FINAL_sub_proportions_", CELL_TYPE_TO_SUBCLUSTER), OUTPUT_DIR)
    for (group in ADDITIONAL_GROUPS_TO_PLOT) {
      plot_cell_proportions(data_sub, "sub_cell_types", group,
                            paste0("FINAL_sub_proportions_", CELL_TYPE_TO_SUBCLUSTER, "_by_", group),
                            OUTPUT_DIR)
    }

    # --- Faceted sub-type UMAP ---
    if (length(ADDITIONAL_GROUPS_TO_PLOT) > 0 &&
        ADDITIONAL_GROUPS_TO_PLOT[1] %in% colnames(data_sub@meta.data)) {
      pgrp <- ADDITIONAL_GROUPS_TO_PLOT[1]
      n_levels <- n_distinct(data_sub@meta.data[[pgrp]])
      p_sub_facet <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                             split.by = pgrp, label = FALSE, repel = TRUE) +
        facet_wrap(as.formula(paste("~", pgrp)), ncol = 2) +
        theme(strip.text = element_text(size = 13, face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 4)))
      ggsave(file.path(OUTPUT_DIR, paste0("FINAL_sub_UMAP_faceted_", CELL_TYPE_TO_SUBCLUSTER, "_by_", pgrp, ".png")),
             p_sub_facet, width = 12, height = 5 * ceiling(n_levels / 2), dpi = DPI_SETTING, bg = "white")
    }

    # --- Save sub-cluster object ---
    saveRDS(data_sub,
            file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_", CELL_TYPE_TO_SUBCLUSTER, "_subclustered.rds")))
    message(paste("  Sub-clustering complete. Object saved for:", CELL_TYPE_TO_SUBCLUSTER))
  }
}

# =============================================================================
# --- PART 7: FINAL SAVE ------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving Final Annotated Object ===")
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds")))
message(paste0(
  "\n=== ANNOTATION COMPLETE ===\n",
  "  Final object: ", PROJECT_NAME, "_final_annotated.rds\n",
  "  Total cells:  ", ncol(data), "\n",
  "  Cell types:   ", paste(levels(data$broad_cell_types), collapse = ", "), "\n",
  "\nNext steps:\n",
  "  - Script 03: Differential expression analysis\n",
  "  - Script 04: Pathway scoring (AUCell) and gene expression comparisons\n"
))
