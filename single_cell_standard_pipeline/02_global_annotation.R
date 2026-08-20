# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 2: GLOBAL ANNOTATION
# Version: 3.1 (CSV-Driven, Weighted Pre-Scoring, Broad Cell Types Only)
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
#
# PURPOSE:
#   This script is the broad/global annotation hub. It loads the processed
#   Seurat object from Script 01 and provides:
#     1. WEIGHTED PRE-SCORING: Automatically scores every cluster against all
#        broad cell type markers using a specificity-weighted gene expression
#        metric (both z-scored and raw). This gives a quantitative "best guess"
#        annotation BEFORE manual review — the visual and numeric output should
#        be used as a guide, not a final answer.
#     2. MANUAL BROAD ANNOTATION: The biologist fills in BROAD_ANNOTATION_MAP
#        after reviewing the UMAP + DotPlot + weighted score results.
#     3. COMPOSITIONAL ANALYSIS: Cell type proportion plots by SampleID and
#        any additional metadata groups (e.g., Genotype_Diet, Condition).
#     4. GENE EXPRESSION COMPARISON: Violin/barplot comparisons of any gene
#        across groups with statistical tests (Wilcoxon).
#
#   NOTE: Sub-population identification has been moved to dedicated scripts:
#     → 03_tcell_subannotation.R      (T cell sub-types)
#     → 04_colonocyte_subannotation.R (Colonocyte sub-types)
#   Add additional sub-annotation scripts following the same pattern for any
#   other broad cell type (Macrophages, Fibroblasts, etc.).
#
# MARKER SYSTEM:
#   All marker genes are loaded from a CSV file (cell_type_markers.csv).
#   Required CSV columns:
#     cell_type            — name of the cell type
#     tier                 — "broad" (used here) or "sub" (used in 03/04)
#     parent_cell_type     — for sub-tier rows, the broad type they belong to
#     markers              — pipe-separated gene list (e.g. "Cd3e|Cd4|Cd8a")
#     dotplot_markers      — pipe-separated genes for the diagnostic DotPlot
#     subcluster_resolution— suggested Leiden resolution for sub-clustering
#
# HOW TO USE:
#   1. Fill in Part 1 (paths and options).
#   2. Run the script block by block in RStudio.
#   3. After each weighted scoring UMAP, review the top5 report table.
#   4. Fill in BROAD_ANNOTATION_MAP (ACTION 4).
#   5. Run the rest of the script for final plots and save.
#   6. Proceed to 03_tcell_subannotation.R / 04_colonocyte_subannotation.R.
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
PROJECT_NAME <- "Wu_Diet_project2"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"  # Windows
ROOT_PATH <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"

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

# --- 1.4: Compositional Analysis Groups ---
# See Script 01/02 comments for full explanation of concatenated columns.
# "Genotype_Diet" is auto-built from Genotype + Diet if both exist in metadata.

ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_Diet")

# --- 1.5: Gene Expression Comparison (generate_gene_comparison_plots) ---
data$Genotype_sex <- data$group
COMPARISON_X_AXIS  <- "Genotype_sex"
COMPARISON_GROUPS  <- c("WT_Female", "Polyp_Female",  "Polyp_NR4a1_KO_Female", 
                        "WT_Male", "Polyp_Male", "Polyp_NR4a1_KO_Male")
                        
# Define the custom order
custom_levels <- c(
  "WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
  "WT_Male", "Polyp_Male",   "Polyp_NR4a1_KO_Male"
)

# Apply it to the Seurat object
data$Genotype_sex <- data$group
data$Genotype_sex <- factor(
  data$Genotype_sex, 
  levels = custom_levels
)

COMPARISON_PAIRS   <- list(
  c("Polyp_Female",    "Polyp_NR4a1_KO_Female"),
  c("Polyp_Male",      "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_NR4a1_KO_Female"),
  c("WT_Male",         "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_Female"),
  c("WT_Male",         "Polyp_Male")
)

# --- 1.6: Output Resolution ---
DPI_SETTING <- 300

# =============================================================================
# --- PART 2: LOAD DATA & MARKERS --------------------------------------------
# =============================================================================
message("=== Loading processed Seurat object ===")
data <- readRDS(INPUT_RDS)
message(paste("  Loaded:", ncol(data), "cells,", nrow(data), "genes"))

# --- Load and parse markers via {TamuScDSC} ------------------------------------
# get_markers() reads the same CSV and explodes the pipe-separated gene lists,
# so BROAD_MARKERS_LIST is identical to the old hand-parsed version but the
# parsing lives in one tested place. Point it at MARKERS_CSV_FILE to keep your
# project's file as the source, or drop the argument to use the copy shipped in
# the package (single source of truth across all datasets).
message("=== Loading markers via TamuScDSC::get_markers() ===")
markers_df <- TamuScDSC::get_markers(MARKERS_CSV_FILE)

broad_df           <- markers_df[markers_df$tier == "broad", ]
BROAD_MARKERS_LIST <- TamuScDSC::markers_as_list(broad_df)   # cell_type -> genes, CSV order

# Flatten (first occurrence wins) for the dotplot gene axis.
broad_dotplot_markers <- unique(broad_df$gene)
print(broad_dotplot_markers)

message(paste("  Broad cell types loaded:", length(BROAD_MARKERS_LIST)))

# --- Set up clustering/UMAP keys based on VIEW_MODE -------------------------
if (VIEW_MODE == "harmony") {
  cluster_col       <- "clusters_harmony"
  umap_col          <- "umap_harmony"
  plot_title_suffix <- "(Harmony Corrected)"
} else {
  cluster_col       <- "clusters_none"
  umap_col          <- "umap_none"
  plot_title_suffix <- "(No Harmony)"
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
#'   1. Build a binary gene x cell-type presence matrix.
#'   2. Weight each gene by 1/occurrence (cross-cell-type specificity weight).
#'   3. Normalize each cell-type column so weights sum to 1.
#'   4. Optionally z-score gene expression across cells (standardize = TRUE).
#'   5. Compute per-cluster mean score for each cell type via matrix multiply.
#'   6. Assign the top-scoring cell type to all cells in each cluster.
#'
#' @param seurat_obj           Seurat object with normalized RNA data.
#' @param marker_genes         Named list: cell_type -> character vector of genes.
#' @param cluster_key          Metadata column with cluster IDs.
#' @param standardize_expression If TRUE, z-score expression per gene.
#' @return List with: annotation_vector (per-cell character), top5_report (data.frame)
# ---------------------------------------------------------------------------
get_weighted_annotation <- function(seurat_obj,
                                    marker_genes,
                                    cluster_key,
                                    standardize_expression = TRUE) {
  all_obj_genes <- rownames(seurat_obj)
  marker_genes  <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
  marker_genes  <- marker_genes[sapply(marker_genes, length) > 0]
  if (length(marker_genes) == 0) stop("No valid marker genes found in dataset.")

  all_marker_genes <- sort(unique(unlist(marker_genes)))
  cell_types       <- names(marker_genes)

  W_binary <- matrix(0, nrow = length(all_marker_genes), ncol = length(cell_types),
                     dimnames = list(all_marker_genes, cell_types))
  for (ct in cell_types) { W_binary[marker_genes[[ct]], ct] <- 1.0 }

  gene_occurrence <- rowSums(W_binary)
  W_specificity   <- W_binary / gene_occurrence
  column_sums     <- colSums(W_specificity)
  W_final         <- sweep(W_specificity, 2, column_sums, "/")
  W_final[is.na(W_final)] <- 0

  X_subset <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[all_marker_genes, ]

  if (standardize_expression) {
    X_subset <- t(scale(t(as.matrix(X_subset))))
    X_subset[is.na(X_subset)] <- 0
  } else {
    X_subset <- as.matrix(X_subset)
  }

  W          <- t(W_final)
  clusters   <- seurat_obj[[cluster_key, drop = TRUE]]
  unique_cls <- unique(clusters)

  annotation_vector  <- character(length = ncol(seurat_obj))
  names(annotation_vector) <- colnames(seurat_obj)
  cluster_score_data <- list()

  for (cluster_id in unique_cls) {
    cell_idx <- which(clusters == cluster_id)
    if (length(cell_idx) == 0) next
    raw_scores  <- W %*% X_subset[, cell_idx, drop = FALSE]
    mean_scores <- rowMeans(raw_scores)
    top_idx     <- order(mean_scores, decreasing = TRUE)
    top_ct      <- names(mean_scores)[top_idx]
    top_sc      <- mean_scores[top_idx]
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
#' generate_gene_comparison_plots: Violin or bar plot for a gene/score column.
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
                                           ncols  = 3,
                                           output_dir   = OUTPUT_DIR) {
  message(paste("  Generating comparison plot for:", score_col))
  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis)) %>%
    dplyr::rename(Expression = 1) %>% drop_na()

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
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y", ncol = ncols) +
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
#' process_and_extract_cell_type
#'
#' Subsets a broad cell type, re-runs HVG selection, PCA, and produces BOTH
#' a standard (non-Harmony) and a Harmony-corrected embedding + clustering.
#' Used by sub-annotation scripts (03_tcell_subannotation.R, etc.).
#'
#' @param data           Seurat object with CellType metadata column
#' @param cell_type_name Cell type label to extract
#' @param num_hvg        HVGs for sub-cluster PCA
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
                                          kneigh     = 30) {
  message(paste("  Subsetting and re-clustering:", cell_type_name))
  data_sub <- subset(data, subset = CellType == cell_type_name)

  # Clear stale reductions / graphs from the parent object
  data_sub@reductions <- list()
  data_sub@graphs     <- list()

  # Re-find HVGs specific to this sub-population and run PCA
  data_sub <- FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = num_hvg) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = dims_pca, reduction.name = "pca", verbose = FALSE)
  gc()

  # --- Track A: Standard PCA (no batch correction) ---
  data_sub <- FindNeighbors(data_sub, dims = 1:dims_pca, reduction = "pca",
                             k.param = kneigh, graph.name = "pca_nn", verbose = FALSE) %>%
    FindClusters(resolution = resolution, graph.name = "pca_nn",
                 cluster.name = "clusters_none", verbose = FALSE) %>%
    RunUMAP(dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh,
            min.dist = min_dist, n.epochs = 500,
            reduction.name = "umap_none", verbose = FALSE)
  gc()

  # --- Track B: Harmony batch-corrected ---
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
# --- PART 4: BROAD ANNOTATION WORKFLOW ----------------------------------------
# =============================================================================
message("\n=== STEP 4: Broad Annotation Workflow ===")

# --- 4.1: Initial UMAP and DotPlot ---
p_umap_clusters <- DimPlot(data, reduction = umap_col, group.by = cluster_col, label = TRUE) +
  ggtitle(paste("Clusters", plot_title_suffix)) + NoLegend()
print(p_umap_clusters)
ggsave(file.path(OUTPUT_DIR, paste0("ANNOTATE_01_UMAP_clusters_", VIEW_MODE, ".png")),
       p_umap_clusters, width = 8, height = 7, dpi = DPI_SETTING)

p_dotplot_broad <- DotPlot(data, features = broad_dotplot_markers,
                            group.by = cluster_col, dot.min = 0.05,
                            cols = "RdBu", scale = TRUE, cluster=TRUE) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Broad Markers", plot_title_suffix))
print(p_dotplot_broad)
ggsave(file.path(OUTPUT_DIR, paste0("ANNOTATE_01_DotPlot_broad_markers_", VIEW_MODE, ".png")),
       p_dotplot_broad, width = 14, height = 12, dpi = DPI_SETTING, bg = "white")

# --- 4.2: WEIGHTED PRE-SCORING -----------------------------------------------
message("  Running weighted pre-scoring (standardized / z-scored)...")
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

# --- 4.3: Review pre-scoring results -----------------------------------------
message("\n--- Pre-scoring Top-5 Report (standardized) ---")
print(results_std$top5_report)
write_xlsx(results_std$top5_report, file.path(OUTPUT_DIR, "PRESCORE_top5_report_standardized.xlsx"))
write_xlsx(results_raw$top5_report, file.path(OUTPUT_DIR, "PRESCORE_top5_report_raw.xlsx"))

summarize_annotations(data, "broad_weighted_std")
summarize_annotations(data, "broad_weighted_raw")

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
# Map each cluster NUMBER (as a string) to a broad cell type name.
# Add/remove entries to match your actual cluster count.
BROAD_ANNOTATION_MAP <- c(
  '1'  = 'Colonocytes',  
  '2'  = 'Colonocytes',
  '3'  = 'Colonocytes', 
  '4'  = 'Colonocytes',  
  '5'  = 'Colonocytes',
  '6'  = 'Colonocytes',  
  '7'  = 'Colonocytes',  
  '8'  = 'Fibroblasts',
  '9'  = 'Colonocytes',  
  '10' = 'Colonocytes',  
  '11' = 'Colonocytes',
  '12' = 'Macrophages',  
  '13' = 'Colonocytes',  
  '14' = 'Colonocytes',
  '15' = 'Fibroblasts',  
  '16' = 'Colonocytes',  
  '17' = 'SMCs',
  '18' = 'T cells',  
  '19' = 'SMCs',      
  '20' = 'Plasma B cells',
  '21' = 'SMCs',  
  '22' = 'Colonocytes',         
  '23' = 'cDCs',
  '24' = 'Colonocytes',  
  '25' = 'VECs',         
  '26' = 'Fibroblasts',
  '27' = 'B cells',
  '28' = 'Colonocytes', 
  '29' = 'LECs',
  '30' = 'Colonocytes',
  '31' = 'Colonocytes',        
  '32' = 'SMCs',
  '33' = 'Colonocytes',  
  '34' = 'Pericytes', # not sure
  '35' = 'vSMCs',
  '36' = 'ENs',
  '37' = 'EGCs',
  '38' = 'Colonocytes',
  '39' = 'Colonocytes',
  '40' = 'Adipocytes'
)
# Set the identity to your custom cluster names first
#Idents(data) <-  "clusters_harmony"
# Find markers for Unknown1 compared to everything else
# unk1_markers <- FindMarkers(data, 
#                             ident.1 = "39", 
#                             min.pct = 0.25, 
#                             logfc.threshold = 0.25)
# 
# top30_unk1 <- unk1_markers %>%
#   mutate(gene = rownames(.)) %>%
#   arrange(desc(avg_log2FC)) %>%
#   head(30)
#cat(top30_unk1$gene, sep = "\n")

# =============================================================================
# --- PART 5: APPLY BROAD ANNOTATIONS AND VISUALIZE --------------------------
# =============================================================================
message("\n=== STEP 5: Applying Broad Annotations ===")

data$broad_cell_types <- recode_factor(data[[cluster_col, drop = TRUE]], !!!BROAD_ANNOTATION_MAP)
data$CellType         <- data$broad_cell_types   # Alias used by sub-annotation scripts

broad_type_levels <- unique(c(broad_df$cell_type, as.character(data$broad_cell_types)))
broad_type_levels <- intersect(broad_type_levels, as.character(unique(data$broad_cell_types)))
data$broad_cell_types   <- factor(data$broad_cell_types, levels = broad_type_levels)
data$broad_weighted_std <- factor(data$broad_weighted_std, levels = broad_type_levels)
data$broad_weighted_raw <- factor(data$broad_weighted_raw, levels = broad_type_levels)

# Final annotation UMAPs: manual + both weighted side-by-side
p_final_manual <- DimPlot(data, reduction = umap_col, group.by = "broad_cell_types",
                           label = TRUE, repel = TRUE) +
  ggtitle(paste("Final Manual Annotation", plot_title_suffix))
p_final_std <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_std",
                        label = TRUE, repel = TRUE) + ggtitle("Weighted Pre-Score (z-scored)")
p_final_raw <- DimPlot(data, reduction = umap_col, group.by = "broad_weighted_raw",
                        label = TRUE, repel = TRUE) + ggtitle("Weighted Pre-Score (raw)")

p_triple <- p_final_manual | p_final_std | p_final_raw
p_triple
ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotation_comparison.png"),
       p_triple, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")
ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotated.png"),
       p_final_manual, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# Filter not found markers
# Force the plot to follow the CSV order for the Y-axis
found_cell_types <- unique(data$broad_cell_types)

filtered_markers_list <- BROAD_MARKERS_LIST[names(BROAD_MARKERS_LIST) %in% found_cell_types]
final_dotplot_markers <- unique(unlist(filtered_markers_list, use.names = FALSE))

# data$broad_cell_types <- factor(data$broad_cell_types, 
#                                 levels = names(filtered_markers_list))

p_final_dot <- DotPlot(data, features = final_dotplot_markers,
                        group.by = "broad_cell_types", dot.min = 0.05,
                        cols = "RdBu", scale = TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11))
p_final_dot
ggsave(file.path(OUTPUT_DIR, "FINAL_DotPlot_broad_annotated.png"),
       p_final_dot, width = 8, height = 14, dpi = DPI_SETTING, bg = "white")

message("--- Broad annotation applied and plots saved ---")

# =============================================================================
# --- PART 6: COMPOSITIONAL ANALYSIS ------------------------------------------
# =============================================================================
message("\n=== STEP 6: Compositional Analysis ===")

# --- 1.4: Compositional Analysis Groups --------------------------------------
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_sex")

# --- 1.5: Gene Expression Comparison -----------------------------------------
library(dplyr); library(tidyr)
library(ggplot2); library(ggpubr); library(writexl)

my_comparisons <- COMPARISON_PAIRS

# --- Step 1: Calculate % per sample ---
meta <- data@meta.data

df_pct <- meta %>%
  count(SampleID, Genotype_sex, broad_cell_types) %>%
  group_by(SampleID) %>%
  mutate(Percentage = (n / sum(n)) * 100) %>%
  ungroup()

# --- Step 2: Plot ---
p_stats <- ggplot(df_pct, aes(x = Genotype_sex, y = Percentage, fill = Genotype_sex)) +
  stat_summary(fun = "mean", geom = "bar", width = 0.7, color = "black", alpha = 0.8) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_jitter(shape = 21, color = "black", size = 2, width = 0.15) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",
                     label = "p.signif", size = 4) +
  facet_wrap(~ broad_cell_types, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Cell Type Proportions by Genotype/Sex",
       y = "% of Total Cells", x = NULL, fill = "Genotype_sex") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "top")
# --- Add these specific adjustments to your theme ---
p_stats <- p_stats +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # Keeps Females on top, Males on bottom
  theme(
    legend.position = "top",
    legend.justification = "left",   # This snaps the legend to the left side
    legend.box.just = "left"
  )

p_stats

# --- Step 3: Save ---
#dir.create("proportion_analysis", showWarnings = FALSE)
savef <- paste0(OUTPUT_DIR, "/broad_cell_types_stats_cell_props.png")
ggsave(savef,p_stats, width = 10, height = 10, dpi=DPI_SETTING)
savef <- paste0(OUTPUT_DIR, "/broad_cell_types_cell_props_long_data.xlsx")
write_xlsx(df_pct, savef)



# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 16), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text   = element_text(size = 16, face = "bold"),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 14)
)

# =============================================================================
# plot_expression_custom() - Final Clean Version with Bottom Legend
# =============================================================================
plot_expression_custom <- function(seurat_obj,
                                   gene,
                                   plot_type      = "bar",
                                   group_by       = "SubCellType",
                                   condition_col  = "Genotype",
                                   comparisons    = NULL,     
                                   facet_ncol     = NULL,     
                                   p_width        = NULL,     
                                   p_height       = NULL,
                                   x_angle        = 0,
                                   show_legend    = TRUE,     # New Argument
                                   hide_x_text    = TRUE,     # New Argument
                                   conditions     = NULL,     
                                   groups_to_keep = NULL,    
                                   assay          = "RNA",
                                   layer          = "data",
                                   save_path      = NULL,
                                   colors         = NULL,
                                   dpi            = 300) {
  
  # 1. Gene Check
  if (!gene %in% rownames(seurat_obj)) {
    warning("Gene not found: ", gene)
    return(NULL)
  }
  
  # 2. Extract Data
  expr_vec <- as.numeric(GetAssayData(seurat_obj, assay = assay, layer = layer)[gene, ])
  
  # 3. Handle Defaults
  if (is.null(conditions)) {
    conditions <- levels(factor(seurat_obj@meta.data[[condition_col]]))
  }
  if (is.null(groups_to_keep)) {
    groups_to_keep <- levels(factor(seurat_obj@meta.data[[group_by]]))
  }
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(conditions))
    names(colors) <- conditions
  }
  
  # 4. Filter and Prep Data
  plot_df <- seurat_obj@meta.data %>%
    dplyr::mutate(
      expr      = expr_vec,
      group     = .data[[group_by]],
      condition = .data[[condition_col]]
    ) %>%
    dplyr::filter(group %in% groups_to_keep, condition %in% conditions) %>%
    dplyr::mutate(
      group     = factor(group,     levels = groups_to_keep),
      condition = factor(condition, levels = conditions)
    )
  
  # 5. Comparisons and Theme
  if (is.null(comparisons) && length(conditions) == 2) {
    comparisons <- list(conditions)
  }
  
  plot_theme <- if(exists("cust_theme")) cust_theme else theme_minimal()
  
  # Helper to handle saving logic
  save_plot_and_data <- function(plot_obj, path, summary_df = NULL, w_factor) {
    if (!is.null(path)) {
      ext <- tolower(tools::file_ext(path))
      dev <- if(ext == "svg") "svg" else NULL
      final_w <- if(!is.null(p_width)) p_width else (w_factor * (if(!is.null(facet_ncol)) facet_ncol else length(groups_to_keep)))
      final_h <- if(!is.null(p_height)) p_height else 5
      
      ggplot2::ggsave(filename = path, plot = plot_obj, width = final_w, height = final_h, dpi = dpi, device = dev)
      if (!is.null(summary_df)) {
        write.csv(summary_df, sub("\\.[^.]+$", "_summary.csv", path))
      }
    }
  }
  
  # --- PLOT GENERATION ---
  p <- ggplot(plot_df, aes(x = condition, y = expr, fill = condition))
  
  if (plot_type == "bar") {
    p <- p +
      geom_bar(stat = "summary", fun = "mean", color = "black", width = 0.7) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.25) +
      geom_jitter(shape = 21, color = "black", stroke = 0.4, width = 0.12, size = 1.2, alpha = 0.5) +
      labs(y = "Mean normalized expression")
  } else {
    p <- p +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
      geom_boxplot(width = 0.1, outlier.size = 0.1, fill = "white") +
      labs(y = "Expression")
  }
  
  p <- p + facet_wrap(~ group, ncol = facet_ncol, scales = "free_y") +
    scale_fill_manual(values = colors, name = "Condition") + # Named for the legend
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) + 
    labs(title = bquote(italic(.(gene))), x = NULL) +
    theme_classic() + plot_theme + 
    theme(
      legend.position = if(show_legend) "bottom" else "none",
      # Logic to hide x-axis text
      axis.text.x = if(hide_x_text) element_blank() else element_text(angle = x_angle, hjust = if(x_angle != 0) 1 else 0.5),
      axis.ticks.x = if(hide_x_text) element_blank() else element_line()
    )
  
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(
      comparisons = comparisons, 
      method      = "wilcox.test", 
      label       = "p.signif",
      step.increase = 0.12
    )
  }

  p <- p +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # Keeps Females on top, Males on bottom
    theme(
      legend.position = "top",
      legend.justification = "left",   # This snaps the legend to the left side
      legend.box.just = "left"
    )  
  
  summary_data <- plot_df %>%
    dplyr::group_by(group, condition) %>%
    dplyr::summarise(Mean = mean(expr), SE = sd(expr)/sqrt(n()), n = n(), .groups = "drop")
  
  save_plot_and_data(p, save_path, summary_data, w_factor = 2.5)
  
  return(list(plot = p, summary = summary_data))
}

  
# =============================================================================
# --- PROBE / EXON KO-SUPPRESSION VERIFICATION PLOTS (conditional) ------------
# =============================================================================
# These plots verify knockdown/knockout at the probe level and are meaningful
# ONLY when the run included a probe-capture assay (ADD_PROBE_DATA = TRUE in
# Script 01). They reference Nr4a1 exon probe columns / the 'Nr4a1-cust'
# feature, which do not exist for studies without probes (e.g. the Wu diet
# study). Detect the probe data and skip the whole block cleanly if absent.
PROBE_DATA_PRESENT <- any(grepl("probe_count", colnames(data@meta.data))) ||
  ("Nr4a1-cust" %in% rownames(data)) || ("Nr4a1" %in% rownames(data))

if (!PROBE_DATA_PRESENT) {
  message("  [SKIP] No probe/exon data detected in this object - skipping the ",
          "KO-suppression verification plots (Nr4a1 exon/probe panels).")
}

if (PROBE_DATA_PRESENT) {

# Correct call for Nr4a1
# Recommended call for Nr4a1 using your new function
plot_results <- plot_expression_custom(
  seurat_obj    = data,
  gene          = "Nr4a1",
  plot_type     = "bar",
  group_by      = "broad_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 12,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_stats_plot.png"),
  dpi           = DPI_SETTING
)

# Correct call for Nr4a1-cust
# Recommended call for Nr4a1 using your new function
plot_results <- plot_expression_custom(
  seurat_obj    = data,
  gene          = "Nr4a1-cust",
  plot_type     = "bar", 
  group_by      = "broad_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 12,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_custom_stats_plot.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data,
  gene          = "Nr4a2",
  plot_type     = "bar", 
  group_by      = "broad_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 12,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a2_stats_plot.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data,
  gene          = "Nr4a3",
  plot_type     = "bar", 
  group_by      = "broad_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 12,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a3_stats_plot.png"),
  dpi           = DPI_SETTING
)




# Subset the data to only include the target epithelial cells
Idents(data) <- "broad_cell_types"
colonocytes_only <- subset(data, idents = "Colonocytes")
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Prepare data from the subsetted object
plot_df_functional <- colonocytes_only@meta.data %>%
  dplyr::select(
    Genotype_sex,
    Exon3 = `Nr4a1_Exon.3._KO.Target_probe_count`,
    Exon5 = `Nr4a1_Exon.5._KO.Target_probe_count`
  ) %>%
  # Apply "AND" logic: Both must be > 0 to be "Functional"
  mutate(Functional_Detection = ifelse(Exon3 > 0 & Exon5 > 0, "Functional", "Non-Functional (KO)")) %>%
  group_by(Genotype_sex, Functional_Detection) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(Percentage = (n / sum(n)) * 100)

# 2. Create the side-by-side comparison plot
p_functional <- ggplot(plot_df_functional, aes(x = Genotype_sex, y = Percentage, fill = Functional_Detection)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 4, fontface = "bold") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Functional Nr4a1 Status in Colonocytes",
       subtitle = "Percentage of cells expressing BOTH Exon 3 and Exon 5",
       x = "Condition & Sex",
       y = "% of Colonocytes")

# 3. Save the high-resolution image
save_path_functional <- file.path(OUTPUT_DIR, "Nr4a1_Colonocyte_Functional_Efficiency.png")
ggsave(filename = save_path_functional, plot = p_functional, width = 12, height = 8, dpi = 300, bg = "white")

p_functional

# 2. Reshape and plot
plot_df_single <- colonocytes_only@meta.data %>%
  dplyr::select(
    Genotype_sex,
    `Exon 3` = `Nr4a1_Exon.3._KO.Target_probe_count`,
    `Exon 5` = `Nr4a1_Exon.5._KO.Target_probe_count`,
    `Exon 7` = `Nr4a1_Exon.7_probe_count`
  ) %>%
  tidyr::pivot_longer(cols = starts_with("Exon"), names_to = "Exon_Region", values_to = "Counts") %>%
  mutate(Status = factor(ifelse(Counts == 0, "Zero", "Positive"), levels = c("Zero", "Positive"))) %>%
  group_by(Exon_Region, Genotype_sex, Status) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(Percentage = (n / sum(n)) * 100)

p_single <- ggplot(plot_df_single, aes(x = Status, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.5, size = 3, fontface = "bold") +
  facet_grid(Exon_Region ~ Genotype_sex) +
  labs(title = "Nr4a1 Individual Exon Detection (Colonocytes Only)", y = "% of Colonocytes") +
  theme_bw()

save_path_functional <- file.path(OUTPUT_DIR, "Nr4a1_Exon_Individual_Colonocytes.png")
ggsave(filename = save_path_functional, plot = p_single, width = 12, height = 8, dpi = 300, bg = "white")

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Prepare long-format data with raw/normalized counts per cell
# We do NOT summarise here because the violin needs every individual cell's value
plot_df_violin <- colonocytes_only@meta.data %>%
  dplyr::select(
    Genotype_sex,
    `Exon 3` = `Nr4a1_Exon.3._KO.Target_probe_count`,
    `Exon 5` = `Nr4a1_Exon.5._KO.Target_probe_count`,
    `Exon 7` = `Nr4a1_Exon.7_probe_count`
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("Exon"), 
    names_to = "Exon_Region", 
    values_to = "Expression"
  )

# 2. Generate the Violin Plot
p_violin <- ggplot(plot_df_violin, aes(x = Genotype_sex, y = Expression, fill = Genotype_sex)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  # Adding a boxplot helps show the median/quantiles since most cells are at 0
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~Exon_Region, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Nr4a1 Exon Expression Intensity in Colonocytes",
    subtitle = "Distribution of counts per cell across targeted exons",
    x = "Condition",
    y = "Expression Level (Counts)"
  )

p_violin
# 3. Save the High-Resolution Image
save_path_vln <- file.path(OUTPUT_DIR, "Nr4a1_Exon_Violin_Colonocytes.png")
ggsave(filename = save_path_vln, plot = p_violin, width = 12, height = 15, dpi = 300, bg = "white")

p_violin

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Prepare the mean expression data
plot_df_mean <- colonocytes_only@meta.data %>%
  dplyr::select(
    Genotype_sex,
    `Exon 3` = `Nr4a1_Exon.3._KO.Target_probe_count`,
    `Exon 5` = `Nr4a1_Exon.5._KO.Target_probe_count`,
    `Exon 7` = `Nr4a1_Exon.7_probe_count`
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("Exon"), 
    names_to = "Exon_Region", 
    values_to = "Expression"
  ) %>%
  # Calculate mean and standard error for error bars
  group_by(Exon_Region, Genotype_sex) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE = sd(Expression) / sqrt(n()),
    .groups = "drop"
  )

# 2. Generate the Mean Barplot
p_mean_bar <- ggplot(plot_df_mean, aes(x = Genotype_sex, y = Mean_Expression, fill = Genotype_sex)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE), width = 0.2) +
  facet_wrap(~Exon_Region, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Mean Nr4a1 Exon Expression (Colonocytes)",
    subtitle = "Aggregated signal intensity with SEM error bars",
    x = "Condition",
    y = "Mean Counts per Cell"
  )

# 3. Save the Plot
save_path_mean <- file.path(OUTPUT_DIR, "Nr4a1_Exon_Mean_Barplots_Colonocytes.png")
ggsave(filename = save_path_mean, plot = p_mean_bar, width = 10, height = 12, dpi = 300, bg = "white")

p_mean_bar

}  # end if (PROBE_DATA_PRESENT) - probe/exon KO-suppression verification plots

message("--- Compositional analysis complete ---")

# =============================================================================
# --- PART 7: FINAL SAVE ------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving Final Annotated Object ===")
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds")))
message(paste0(
  "\n=== GLOBAL ANNOTATION COMPLETE ===\n",
  "  Final object: ", PROJECT_NAME, "_final_annotated.rds\n",
  "  Total cells:  ", ncol(data), "\n",
  "  Cell types:   ", paste(levels(data$broad_cell_types), collapse = ", "), "\n",
  "\nNext steps:\n",
  "  - Script 03: T cell sub-annotation   (03_tcell_subannotation.R)\n",
  "  - Script 04: Colonocyte sub-annotation (04_colonocyte_subannotation.R)\n",
  "  - Script 05: Differential expression + ANOVA (05_DE_and_two_way_ANOVA.R)\n"
))

