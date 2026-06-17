# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 4: COLONOCYTE SUB-ANNOTATION
# Version: 1.0 (CSV-Driven, Seurat Wrappers, Harmony + Standard Clustering)
#
# PURPOSE:
#   Loads the globally annotated Seurat object from Script 02 and performs
#   high-resolution sub-clustering of the Colonocyte compartment. Provides:
#     1. Subset extraction with fresh HVG selection, PCA, Harmony, and UMAP.
#        BOTH a Harmony-corrected and a standard (non-Harmony) embedding are
#        produced for side-by-side comparison.
#     2. Weighted pre-scoring against colonocyte sub-type markers (CSV-driven).
#        Both z-scored (standardized) and raw (non-standardized) variants run
#        automatically — review the Top-5 report before manual annotation.
#     3. Manual sub-annotation via SUB_ANNOTATION_MAP (Action 5).
#     4. Compositional analysis (proportions by SampleID and group).
#     5. Gene expression violin/bar plots.
#
# MARKER SYSTEM (cell_type_markers.csv):
#   This script reads rows where tier == "sub" AND parent_cell_type == "Colonocytes".
#   Expected sub-types:
#     Abs. colonocytes, Goblet cells, EECs, Tuft cells, TA cells, Stem cells
#   If no sub-rows exist in the CSV the built-in DEFAULT_SUB_MARKERS are used.
#
# HOW TO USE:
#   1. Set paths/parameters in Part 1.
#   2. Run through Part 4 — review SUBCLUSTER_01 and SUBCLUSTER_02 PNGs.
#   3. Fill in SUB_ANNOTATION_MAP in Part 5 (Action 5).
#   4. Run the rest for final plots and save.
#
# INPUT:  {PROJECT_NAME}_final_annotated.rds   (output of 02_global_annotation.R)
# OUTPUT: {PROJECT_NAME}_Colonocytes_subclustered.rds
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

# --- 1.1: Project Paths (must match Scripts 01 and 02) -----------------------
PROJECT_NAME <- "Nr4a1_s17_ack"
ROOT_PATH <- "/home/ssromerogon/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"  # Windows

OUTPUT_DIR       <- file.path(ROOT_PATH, "seurat_output")
MAIN_RDS         <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))
MARKERS_CSV_FILE <- file.path(ROOT_PATH, "cell_type_markers.csv")

# --- 1.2: Target Cell Type ---------------------------------------------------
# Must exactly match a label in the broad_cell_types / CellType column from 02.
PARENT_CELL_TYPE <- "Colonocytes"

# --- 1.3: Sub-Clustering Parameters ------------------------------------------
SUBCLUSTER_N_HVG       <- 2000   # HVGs for sub-clustering PCA
SUBCLUSTER_N_PCS       <- 50     # PCs used for kNN graph
SUBCLUSTER_K_NEIGHBORS <- 30     # k for kNN
SUBCLUSTER_MIN_DIST    <- 0.2    # UMAP min.dist
# Resolution: set to NULL to read from CSV (subcluster_resolution column),
# or override with a number here (e.g., 3.0 for fine-grained colonocyte clusters).
SUBCLUSTER_RESOLUTION  <- NULL


# --- 1.4: Compositional Analysis Groups --------------------------------------
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_sex")

# --- 1.5: Gene Expression Comparison -----------------------------------------
COMPARISON_X_AXIS  <- "Genotype_sex"
COMPARISON_GROUPS  <- c("WT_Female", "Polyp_Female",  "Polyp_NR4a1_KO_Female", 
                        "WT_Male", "Polyp_Male", "Polyp_NR4a1_KO_Male")

# Define the custom order
custom_levels <- c(
  "WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
  "WT_Male", "Polyp_Male",   "Polyp_NR4a1_KO_Male"
)

COMPARISON_PAIRS   <- list(
  c("Polyp_Female",    "Polyp_NR4a1_KO_Female"),
  c("Polyp_Male",      "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_NR4a1_KO_Female"),
  c("WT_Male",         "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_Female"),
  c("WT_Male",         "Polyp_Male")
)

# --- 1.6: Output Resolution --------------------------------------------------
DPI_SETTING <- 300

# =============================================================================
# --- PART 2: LOAD DATA & MARKERS --------------------------------------------
# =============================================================================
message("=== Loading globally annotated Seurat object ===")
data <- readRDS(MAIN_RDS)
message(paste("  Loaded:", ncol(data), "cells"))

# --- Parse markers from CSV --------------------------------------------------
parse_markers <- function(marker_string) {
  strsplit(trimws(marker_string), "\\|")[[1]]
}

# Default fallback sub-markers for Colonocytes
DEFAULT_SUB_MARKERS <- list(
  "Abs. colonocytes" = c("Alpi", "Ces2c", "Slc26a2", "Ceacam1", "Aqp8"),
  "Goblet cells"     = c("Muc2", "Spink4", "Agr2", "Tff3"),
  "EECs"             = c("Chga", "Chgb", "Tph1", "Cck", "Scgn", "Scg5"),
  "Tuft cells"       = c("Dclk1", "Trpm5", "Avil", "Sh2d6", "Plcg2"),
  "TA cells"         = c("Mki67", "Top2a", "Birc5", "Pcna", "Stmn1"),
  "Stem cells"       = c("Lgr5", "Lrig1", "Ascl2", "Slc12a2", "Smoc2",
                          "Kcnq1", "Gpx2", "Ephb2", "Bmpr1a", "Hopx", "Sox9")
)

SUB_MARKERS_LIST    <- DEFAULT_SUB_MARKERS
sub_dotplot_markers <- unique(unlist(SUB_MARKERS_LIST))
SUBCLUSTER_RESOLUTION_FINAL <- if (!is.null(SUBCLUSTER_RESOLUTION)) SUBCLUSTER_RESOLUTION else 3.0

if (file.exists(MARKERS_CSV_FILE)) {
  markers_df <- read.csv(MARKERS_CSV_FILE, stringsAsFactors = FALSE)
  sub_df     <- markers_df[markers_df$tier == "sub" &
                             markers_df$parent_cell_type == PARENT_CELL_TYPE, ]
  if (nrow(sub_df) > 0) {
    SUB_MARKERS_LIST    <- setNames(lapply(sub_df$markers, parse_markers), sub_df$cell_type)
    sub_dotplot_markers <- unique(unlist(lapply(sub_df$dotplot_markers, parse_markers)))
    if (length(sub_dotplot_markers) == 0) sub_dotplot_markers <- unique(unlist(SUB_MARKERS_LIST))
    if (is.null(SUBCLUSTER_RESOLUTION) && !is.na(sub_df$subcluster_resolution[1])) {
      SUBCLUSTER_RESOLUTION_FINAL <- as.numeric(sub_df$subcluster_resolution[1])
    }
    message(paste("  CSV sub-markers loaded for:", PARENT_CELL_TYPE,
                  "| sub-types:", length(SUB_MARKERS_LIST),
                  "| resolution:", SUBCLUSTER_RESOLUTION_FINAL))
  } else {
    warning(paste0("[WARN] No sub-tier rows found for '", PARENT_CELL_TYPE,
                   "' in CSV. Using built-in defaults. Add rows with:\n",
                   "  tier=sub, parent_cell_type=Colonocytes, cell_type=<subtype>, markers=<pipe-sep>"))
  }
} else {
  warning(paste("[WARN] Marker CSV not found at:", MARKERS_CSV_FILE,
                "\nUsing built-in default colonocyte sub-markers."))
}

message(paste("  Sub-clustering resolution:", SUBCLUSTER_RESOLUTION_FINAL))

# =============================================================================
# --- PART 3: UTILITY FUNCTIONS -----------------------------------------------
# =============================================================================

get_weighted_annotation <- function(seurat_obj, marker_genes, cluster_key,
                                    standardize_expression = TRUE) {
  all_obj_genes <- rownames(seurat_obj)
  marker_genes  <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
  marker_genes  <- marker_genes[sapply(marker_genes, length) > 0]
  if (length(marker_genes) == 0) stop("No valid marker genes found in dataset.")

  all_marker_genes <- sort(unique(unlist(marker_genes)))
  cell_types       <- names(marker_genes)

  W_binary <- matrix(0, nrow = length(all_marker_genes), ncol = length(cell_types),
                     dimnames = list(all_marker_genes, cell_types))
  for (ct in cell_types) W_binary[marker_genes[[ct]], ct] <- 1.0

  gene_occurrence <- rowSums(W_binary)
  W_specificity   <- W_binary / gene_occurrence
  W_final         <- sweep(W_specificity, 2, colSums(W_specificity), "/")
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
  annotation_vector <- character(length = ncol(seurat_obj))
  names(annotation_vector) <- colnames(seurat_obj)
  cluster_score_data <- list()

  for (cluster_id in unique_cls) {
    cell_idx <- which(clusters == cluster_id)
    if (length(cell_idx) == 0) next
    raw_scores  <- W %*% X_subset[, cell_idx, drop = FALSE]
    mean_scores <- rowMeans(raw_scores)
    top_idx     <- order(mean_scores, decreasing = TRUE)
    annotation_vector[cell_idx] <- names(mean_scores)[top_idx[1]]
    n_top <- min(5, length(top_idx))
    cluster_score_data[[as.character(cluster_id)]] <- data.frame(
      cluster   = cluster_id,
      Rank      = 1:n_top,
      Cell_Type = names(mean_scores)[top_idx[1:n_top]],
      Score     = mean_scores[top_idx[1:n_top]]
    )
  }
  top5_report <- do.call(rbind, cluster_score_data)
  return(list(annotation_vector = annotation_vector, top5_report = top5_report))
}


process_and_extract_cell_type <- function(data, cell_type_name,
                                          num_hvg    = 2000,
                                          dims_pca   = 50,
                                          resolution = 3.0,
                                          min_dist   = 0.3,
                                          kneigh     = 15) {
  message(paste("  Subsetting and re-clustering:", cell_type_name))
  data_sub <- subset(data, subset = CellType == cell_type_name)
  data_sub@reductions <- list(); data_sub@graphs <- list()

  data_sub <- FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = num_hvg) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = dims_pca, reduction.name = "pca", verbose = FALSE)
  gc()

  # Track A: Standard PCA
  data_sub <- FindNeighbors(data_sub, dims = 1:dims_pca, reduction = "pca",
                             k.param = kneigh, graph.name = "pca_nn", verbose = FALSE) %>%
    FindClusters(resolution = resolution, graph.name = "pca_nn",
                 cluster.name = "clusters_none", verbose = FALSE) %>%
    RunUMAP(dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh,
            min.dist = min_dist, n.epochs = 500,
            reduction.name = "umap_none", verbose = FALSE)
  gc()

  # Track B: Harmony
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

generate_gene_comparison_plots <- function(seurat_obj, score_col, group_by, x_axis,
                                           comparisons, plot_type = "violin",
                                           output_prefix = "", plot_title = score_col,
                                           y_label = "Expression",
                                           fig_width = 16, fig_height = 7,
                                           output_dir = OUTPUT_DIR) {
  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis)) %>%
    dplyr::rename(Expression = 1) %>% drop_na()

  cust_theme <- theme_classic() + theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, face = "bold"),
    legend.position = "bottom", panel.spacing = unit(1.5, "lines")
  )
  p <- ggplot(df_plot, aes(!!sym(x_axis), Expression, fill = !!sym(x_axis)))
  if (plot_type == "barplot") {
    p <- p + stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.8) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else {
    p <- p + geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5)
  }
  p <- p +
    stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox.test",
                       method.args = list(exact = FALSE),
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols = c("****", "***", "**", "*", "ns")),
                       step.increase = 0.1, size = 6, bracket.size = 0.8) +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(title = plot_title, y = y_label) + scale_fill_brewer(palette = "Set1") + cust_theme

  plot_file <- file.path(output_dir, paste0(output_prefix, score_col, ".png"))
  ggsave(plot_file, p, width = fig_width, height = fig_height, dpi = DPI_SETTING, bg = "white")
  return(p)
}

# =============================================================================
# --- PART 4: SUB-CLUSTERING --------------------------------------------------
# =============================================================================
message(paste("\n=== STEP 4: Sub-Clustering of:", PARENT_CELL_TYPE, "==="))

if (!"CellType" %in% colnames(data@meta.data)) {
  if ("broad_cell_types" %in% colnames(data@meta.data)) {
    data$CellType <- data$broad_cell_types
  } else {
    stop("[ERROR] Neither 'CellType' nor 'broad_cell_types' found. Run Script 02 first.")
  }
}

n_parent <- sum(data@meta.data$CellType == PARENT_CELL_TYPE, na.rm = TRUE)
if (n_parent < 50) stop(paste("[ERROR] Only", n_parent, PARENT_CELL_TYPE,
                               "cells found. Check PARENT_CELL_TYPE matches exactly."))
message(paste("  Found", n_parent, PARENT_CELL_TYPE, "cells for sub-clustering."))

data_sub <- process_and_extract_cell_type(
  data           = data,
  cell_type_name = PARENT_CELL_TYPE,
  num_hvg        = SUBCLUSTER_N_HVG,
  dims_pca       = SUBCLUSTER_N_PCS,
  resolution     = SUBCLUSTER_RESOLUTION_FINAL,
  min_dist       = SUBCLUSTER_MIN_DIST,
  kneigh         = SUBCLUSTER_K_NEIGHBORS
)

# --- 4.1: Harmony vs Standard Comparison ---
p_harm_comp <- (
  DimPlot(data_sub, reduction = "umap_none",    group.by = "SampleID") + ggtitle("Standard PCA (No Harmony)")
) + (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "SampleID") + ggtitle("Harmony Corrected")
) & theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_01_Harmony_comparison_Colonocytes.png"),
       p_harm_comp, width = 16, height = 8, dpi = DPI_SETTING)

# --- 4.2: Cluster UMAPs ---
p_clust_none <- DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_none",
                         label = TRUE) + NoLegend() + ggtitle("Clusters: Standard PCA")
p_clust_harm <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony",
                         label = TRUE) + NoLegend() + ggtitle("Clusters: Harmony")
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_COMPARE_UMAP_Colonocytes.png"),
       p_clust_none + p_clust_harm, width = 16, height = 8, dpi = DPI_SETTING)
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_THIS_UMAP_Colonocytes.png"),
       p_clust_harm, width = 12, height = 8, dpi = DPI_SETTING)

# --- 4.3: DotPlot ---
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers, scale = T,
                      group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu", cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_Colonocytes.png"),
       p_dot_sub, width = 14, height = 20, dpi = DPI_SETTING, bg = "white")

# not scaled
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers, scale = F, 
                     cols = c("Dark Violet", "Magenta"),
                     group.by = "clusters_harmony", dot.min = 0.05, cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_NOT_SCALED_Colonocytes.png"),
       p_dot_sub, width = 14, height = 20, dpi = DPI_SETTING, bg = "white")

# =============================================================================
# --- PART 5: WEIGHTED PRE-SCORING -------------------------------------------
# =============================================================================
message("  Running sub-type weighted pre-scoring (standardized)...")
sub_results_std <- get_weighted_annotation(
  seurat_obj = data_sub, marker_genes = SUB_MARKERS_LIST,
  cluster_key = "clusters_harmony", standardize_expression = TRUE
)
data_sub$sub_weighted_std <- sub_results_std$annotation_vector

message("  Running sub-type weighted pre-scoring (non-standardized)...")
sub_results_raw <- get_weighted_annotation(
  seurat_obj = data_sub, marker_genes = SUB_MARKERS_LIST,
  cluster_key = "clusters_harmony", standardize_expression = FALSE
)
data_sub$sub_weighted_raw <- sub_results_raw$annotation_vector

message("\n--- Sub-cluster Pre-scoring Top-5 Report (standardized) ---")
print(sub_results_std$top5_report)
write_xlsx(sub_results_std$top5_report,
           file.path(OUTPUT_DIR, "SUBCLUSTER_PRESCORE_top5_Colonocytes_std.xlsx"))
write_xlsx(sub_results_raw$top5_report,
           file.path(OUTPUT_DIR, "SUBCLUSTER_PRESCORE_top5_Colonocytes_raw.xlsx"))

p_prescore <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Raw")
)
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_03_prescore_comparison_Colonocytes.png"),
       p_prescore, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

# =============================================================================
# --- ACTION 5: FILL IN SUB-CLUSTER ANNOTATIONS HERE ---
# =============================================================================
# Review:
#   SUBCLUSTER_02_ANNOTATE_THIS_UMAP_Colonocytes.png
#   SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_Colonocytes.png
#   SUBCLUSTER_03_prescore_comparison_Colonocytes.png
#   SUBCLUSTER_PRESCORE_top5_Colonocytes_std.xlsx
#
# Available sub-types (defaults or from CSV):
#   "Abs. colonocytes", "Goblet cells", "EECs", "Tuft cells", "TA cells", "Stem cells"
SUB_ANNOTATION_MAP <- c(
  '0'  = 'Stem cells',  
  '1'  = 'Stem cells',
  '2'  = 'Abs. colonocytes',
  '3'  = 'Goblet cells',
  '4'  = 'Goblet cells',
  '5'  = 'Abs. colonocytes',
  '6'  = 'Abs. colonocytes',
  '7'  = 'Abs. colonocytes',
  '8'  = 'Abs. colonocytes',
  '9'  = 'Abs. colonocytes',
  '10' = 'Abs. colonocytes',
  '11' = 'Abs. colonocytes',
  '12' = 'TA cells',
  '13' = 'Stem cells', # Great stem signature
  '14' = 'Abs. colonocytes',
  '15' = 'Prol. Goblet cells', # but looks like TA/Stem signature
  '16' = 'Stem cells',
  '17' = 'Abs. colonocytes',
  '18' = 'Abs. colonocytes',
  '19' = 'Goblet cells',
  '20' = 'Abs. colonocytes',
  '21' = 'Prol. Goblet cells',# but looks like TA/Stem signature
  '22' = 'Abs. colonocytes',
  '23' = 'Abs. colonocytes',
  '24' = 'Abs. colonocytes',
  '25' = 'Abs. colonocytes',
  '26' = 'Abs. colonocytes',
  '27' = 'TA cells',
  '28' = 'Goblet cells',
  '29' = 'Goblet cells',
  '30' = 'TA cells',
  '31' = 'Abs. colonocytes',
  '32' = 'Abs. colonocytes',
  '33' = 'Abs. colonocytes',
  '34' = 'Abs. colonocytes',
  '35' = 'Abs. colonocytes',
  '36' = 'Abs. colonocytes',
  '37' = 'Goblet cells',
  '38' = 'Abs. colonocytes',
  '39' = 'Abs. colonocytes', # or Stem cells
  '40' = 'Prol. Goblet cells', # but looks like TA/Stem signature
  '41' = 'Goblet cells',
  '42' = 'TA cells', # or Stem cells
  '43' = 'TA cells', # or Stem cells
  '44' = 'TA cells',
  '45' = 'Abs. colonocytes',
  '46' = 'Abs. colonocytes',
  '47' = 'EECs',
  '48' = 'Goblet cells',
  '49' = 'Goblet cells',
  '50' = 'Tuft cells',
  '51' = 'EECs',
  '52' = 'TA cells',
  '53' = 'Prol. Goblet cells', # This may be cancer derived
  '54' = 'Abs. colonocytes',
  '55' = 'Prol. Goblet cells', # but looks like TA/Stem signature
  '56' = 'Stem cells',
  '57' = 'Abs. colonocytes'
  # Add/remove entries matching your actual cluster count
)

# =============================================================================
# --- PART 6: APPLY SUB-ANNOTATIONS & FINAL VISUALIZATIONS -------------------
# =============================================================================
message("\n=== STEP 6: Applying Sub-Annotations ===")

data_sub$sub_cell_types  <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony

DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
        label = TRUE, repel = TRUE)

DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
        label = TRUE, repel = TRUE, split.by = "Genotype_sex")

cell_table <- data_sub@meta.data %>%
  group_by(Genotype_sex, sub_cell_types) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Genotype_sex) %>%
  mutate(pct = n / sum(n) * 100) %>%
  select(-n) %>%
  pivot_wider(names_from = Genotype_sex, values_from = pct, values_fill = 0)

print(cell_table)

# sub_type_levels <- names(SUB_MARKERS_LIST)
# sub_type_levels <- intersect(sub_type_levels,
#                               as.character(unique(data_sub$sub_cell_types)))
# data_sub$sub_cell_types   <- factor(data_sub$sub_cell_types,   levels = sub_type_levels)
# data_sub$sub_weighted_std <- factor(data_sub$sub_weighted_std, levels = sub_type_levels)
# data_sub$sub_weighted_raw <- factor(data_sub$sub_weighted_raw, levels = sub_type_levels)


sub_type_levels2 <- c("EECs", "Tuft cells", "Goblet cells", "Prol. Goblet cells",
                      "TA cells", "Stem cells", "Abs. colonocytes")
data_sub$sub_cell_types   <- factor(data_sub$sub_cell_types,   levels = sub_type_levels2)


# Final UMAP
p_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = FALSE, repel = TRUE) +
  ggtitle(paste("Sub-Cell Types:", PARENT_CELL_TYPE)) +
  theme(legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(OUTPUT_DIR, "FINAL_SUBCLUSTER_UMAP_Colonocytes.png"),
       p_final, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# Comparison: manual vs weighted
p_compare <- DimPlot(data_sub, reduction = "umap_harmony",
                     group.by = c("sub_cell_types", "sub_weighted_std", "sub_weighted_raw"),
                     label = FALSE, repel = TRUE)
ggsave(file.path(OUTPUT_DIR, "FINAL_SUBCLUSTER_comparison_Colonocytes.png"),
       p_compare, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")





# =============================================================================
# 1. DEFINE LEVELS AND MAPPING OF DOTPLOT
# =============================================================================
# This is the order you want for your plot axes/legend
sub_type_levels2 <- c("EECs", "Tuft cells", "Goblet cells", "Prol. Goblet cells",
                      "TA cells", "Stem cells", "Abs. colonocytes")
# This maps your Factor Levels (left) to the SUB_MARKERS_LIST names (right)
# This handles the "Cyc. CD4+ T cells" -> "Cyc. T cells" mismatch
name_mapping <- c(
  "EECs"   =  "EECs", 
  "Tuft cells" = "Tuft cells", 
  "Goblet cells" = "Goblet cells", 
  "Prol. Goblet cells" = "Goblet cells",
  "TA cells" = "TA cells", 
  "Stem cells" = "Stem cells", 
  "Abs. colonocytes" = "Abs. colonocytes"
  )


# Apply levels to the Seurat object
data_sub$sub_cell_types <- factor(
  data_sub$sub_cell_types, 
  levels = sub_type_levels2
)

# =============================================================================
# 2. EXTRACT AND ORCHESTRATE GENE LIST
# =============================================================================

# A. Set the mandatory lead genes
lead_genes <- c("Epcam")

# B. Extract markers based on the ordered factor levels and the mapping
# We use lapply to ensure we follow 'sub_type_levels2' order exactly
ordered_markers <- unlist(lapply(sub_type_levels2, function(lvl) {
  marker_key <- name_mapping[lvl]
  
  # Check if the key exists in your marker list to avoid errors
  if (!is.na(marker_key) && marker_key %in% names(SUB_MARKERS_LIST)) {
    return(SUB_MARKERS_LIST[[marker_key]])
  } else {
    return(NULL)
  }
}))

# C. Combine, remove duplicates, and verify existence in the dataset
# unique() ensures Cd3e/g stay at the top and don't repeat later
final_gene_list <- c(lead_genes, ordered_markers)
final_gene_list <- unique(final_gene_list)

# D. Final check against the Seurat object's actual genes
final_gene_list <- intersect(final_gene_list, rownames(data_sub))

# --- Final DotPlot ---
p_final_dot <- DotPlot(data_sub, features = final_gene_list,
                        group.by = "sub_cell_types", dot.min = 0.05,
                        cols = "RdBu", scale = TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11))
ggsave(file.path(OUTPUT_DIR, "FINAL_DotPlot_sub_Colonocytes.png"),
       p_final_dot, width = 8, height = 12, dpi = DPI_SETTING, bg = "white")


# Faceted UMAP
if (length(ADDITIONAL_GROUPS_TO_PLOT) > 0 &&
    ADDITIONAL_GROUPS_TO_PLOT[1] %in% colnames(data_sub@meta.data)) {
  pgrp     <- ADDITIONAL_GROUPS_TO_PLOT[1]
  n_levels <- n_distinct(data_sub@meta.data[[pgrp]])
  p_facet  <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                      split.by = pgrp, label = FALSE, repel = TRUE) +
    facet_wrap(as.formula(paste("~", pgrp)), ncol = 3) +
    theme(strip.text = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 4))) 
  ggsave(file.path(OUTPUT_DIR, paste0("FINAL_sub_UMAP_faceted_Colonocytes_by_", pgrp, ".png")),
         p_facet, width = 5 * ceiling(n_levels / 2), height = 10, dpi = DPI_SETTING, bg = "white")
}



# =============================================================================
# --- PART 7: SAVE ------------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving Colonocyte Sub-Cluster Object ===")
saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_Colonocytes_subclustered.rds")))
message(paste0(
  "\n=== COLONOCYTE SUB-ANNOTATION COMPLETE ===\n",
  "  Saved: ", PROJECT_NAME, "_Colonocytes_subclustered.rds\n",
  "  Sub-types found: ", paste(levels(data_sub$sub_cell_types), collapse = ", "), "\n",
  "\nNext step:\n",
  "  - Script 05: DE + Two-Way ANOVA (05_DE_and_two_way_ANOVA.R)\n"
))

# =============================================================================
# --- PART 8: COMPOSITIONAL ANALYSIS ------------------------------------------
# =============================================================================
message("\n=== STEP 8: Compositional Analysis ===")

# --- 1.4: Compositional Analysis Groups --------------------------------------
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_sex")

# --- 1.5: Gene Expression Comparison -----------------------------------------
library(dplyr); library(tidyr)
library(ggplot2); library(ggpubr); library(writexl)

my_comparisons <- COMPARISON_PAIRS

# --- Step 1: Calculate % per sample ---
meta <- data_sub@meta.data

df_pct <- meta %>%
  count(SampleID, Genotype_sex, sub_cell_types) %>%
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
  facet_wrap(~ sub_cell_types, scales = "free_y") +
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
savef <- paste0(OUTPUT_DIR, "/colonocytes_stats_cell_props.png")
ggsave(savef,p_stats, width = 10, height = 10, dpi=DPI_SETTING)
savef <- paste0(OUTPUT_DIR, "/colonocytes_cell_props_long_data.xlsx")
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


# Correct call for Nr4a1
# Recommended call for Nr4a1 using your new function
plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Nr4a1",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_stats_plot_colonocytes.png"),
  dpi           = DPI_SETTING
)

# Correct call for Nr4a1-cust
# Recommended call for Nr4a1 using your new function
plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Nr4a1-cust",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_custom_stats_plot_colonocytes.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Nr4a2",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a2_stats_plot_colonocytes.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Nr4a3",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(OUTPUT_DIR, "Nr4a3_stats_plot_colonocytes.png"),
  dpi           = DPI_SETTING
)


