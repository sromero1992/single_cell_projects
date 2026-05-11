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
PROJECT_NAME <- "Nr4a1_Study17_Project"
ROOT_PATH    <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"
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
SUBCLUSTER_K_NEIGHBORS <- 15     # k for kNN
SUBCLUSTER_MIN_DIST    <- 0.3    # UMAP min.dist
# Resolution: set to NULL to read from CSV (subcluster_resolution column),
# or override with a number here (e.g., 3.0 for fine-grained colonocyte clusters).
SUBCLUSTER_RESOLUTION  <- NULL

# --- 1.4: Compositional Analysis Groups --------------------------------------
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_Diet")

# --- 1.5: Gene Expression Comparison -----------------------------------------
COMPARISON_X_AXIS <- "Genotype_Diet"
COMPARISON_GROUPS <- c("WT_cellulose", "KO_cellulose", "WT_inulin", "KO_inulin")
COMPARISON_PAIRS  <- list(
  c("WT_inulin",    "KO_inulin"),
  c("WT_cellulose", "KO_cellulose"),
  c("WT_cellulose", "WT_inulin")
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

summarize_annotations <- function(seurat_obj, annotation_column, print_summary = TRUE) {
  if (!annotation_column %in% colnames(seurat_obj@meta.data))
    stop(paste("Column", annotation_column, "not found in metadata."))
  counts <- table(seurat_obj[[annotation_column]])
  summary_df <- data.frame(
    Count      = as.numeric(counts),
    Percentage = as.numeric(prop.table(counts) * 100),
    row.names  = names(counts)
  ) %>% arrange(desc(Percentage))
  if (print_summary) {
    tmp <- summary_df; tmp$Percentage <- sprintf("%.2f%%", tmp$Percentage)
    message(paste("\n--- Annotation Summary:", annotation_column, "---")); print(tmp)
  }
  return(summary_df)
}

plot_cell_proportions <- function(seurat_obj, cluster_col, group_col,
                                  output_prefix, output_dir) {
  meta <- seurat_obj@meta.data
  make_label <- function(p) ifelse(p > 2.5, paste0(round(p, 1), "%"), "")

  if (!group_col %in% colnames(meta)) {
    parts <- strsplit(group_col, "_")[[1]]
    if (all(parts %in% colnames(meta))) {
      seurat_obj[[group_col]] <- apply(meta[, parts, drop = FALSE], 1, paste, collapse = "_")
      meta <- seurat_obj@meta.data
    } else { warning(paste("Skipping:", group_col, "not found.")); return(NULL) }
  }

  df_sample <- meta %>% group_by(SampleID, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>% group_by(SampleID) %>%
    mutate(percentage = n / sum(n) * 100)
  df_group  <- meta %>% group_by(!!sym(group_col), !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>% group_by(!!sym(group_col)) %>%
    mutate(percentage = n / sum(n) * 100)
  df_global <- meta %>% group_by(!!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>% mutate(percentage = n / sum(n) * 100)

  write_xlsx(list(By_Sample = df_sample, By_Group = df_group, Global = df_global),
             path = file.path(output_dir, paste0(output_prefix, "_Stats.xlsx")))

  pt <- theme_classic() + theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
                                 axis.title.x = element_blank(), legend.position = "bottom")
  p1 <- ggplot(df_sample, aes(SampleID, percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = paste("Proportions by Sample:", cluster_col), y = "Percentage (%)") + pt
  p2 <- ggplot(df_group, aes(!!sym(group_col), percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = paste("Proportions by", group_col), y = "Percentage (%)") + pt
  p3 <- ggplot(df_global, aes("Global", percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", width = 0.6) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_stack(vjust = 0.5), size = 4.5) +
    labs(title = paste("Global Distribution:", cluster_col), y = "Percentage (%)") + pt

  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_by_Sample.png")),    p1,
         width = max(8, 1.5 * n_distinct(meta$SampleID)), height = 7, bg = "white", dpi = DPI_SETTING)
  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_by_", group_col, ".png")), p2,
         width = max(8, 1.5 * n_distinct(meta[[group_col]])), height = 7, bg = "white", dpi = DPI_SETTING)
  ggsave(file.path(output_dir, paste0(output_prefix, "_pct_Global.png")), p3,
         width = 7, height = 7, bg = "white", dpi = DPI_SETTING)
  return(invisible(NULL))
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
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_THIS_UMAP_Colonocytes.png"),
       p_clust_none + p_clust_harm, width = 16, height = 8, dpi = DPI_SETTING)

# --- 4.3: DotPlot ---
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers,
                      group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_Colonocytes.png"),
       p_dot_sub, width = 14, height = 9, dpi = DPI_SETTING, bg = "white")

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
  '0'  = 'Goblet cells',      # REPLACE WITH YOUR ACTUAL ANNOTATIONS
  '1'  = 'Abs. colonocytes',
  '2'  = 'Abs. colonocytes',
  '3'  = 'Goblet cells',
  '4'  = 'Stem cells',
  '5'  = 'Goblet cells',
  '6'  = 'Abs. colonocytes',
  '7'  = 'Goblet cells',
  '8'  = 'Abs. colonocytes',
  '9'  = 'Abs. colonocytes',
  '10' = 'TA cells',
  '11' = 'Abs. colonocytes',
  '12' = 'Abs. colonocytes',
  '13' = 'TA cells',
  '14' = 'TA cells',
  '15' = 'Abs. colonocytes',
  '16' = 'Abs. colonocytes',
  '17' = 'TA cells',
  '18' = 'Stem cells',
  '19' = 'Abs. colonocytes',
  '20' = 'Goblet cells',
  '21' = 'Abs. colonocytes',
  '22' = 'Goblet cells',
  '23' = 'Abs. colonocytes',
  '24' = 'Goblet cells',
  '25' = 'Goblet cells',
  '26' = 'Goblet cells',
  '27' = 'Abs. colonocytes',
  '28' = 'Abs. colonocytes',
  '29' = 'TA cells',
  '30' = 'Abs. colonocytes',
  '31' = 'Abs. colonocytes',
  '32' = 'Stem cells',
  '33' = 'Abs. colonocytes',
  '34' = 'Goblet cells',
  '35' = 'Abs. colonocytes',
  '36' = 'TA cells',
  '37' = 'Abs. colonocytes',
  '38' = 'TA cells',
  '39' = 'Abs. colonocytes',
  '40' = 'EECs',
  '41' = 'Goblet cells',
  '42' = 'TA cells',
  '43' = 'Goblet cells',
  '44' = 'Goblet cells',
  '45' = 'Tuft cells',
  '46' = 'Tuft cells',
  '47' = 'EECs',
  '48' = 'Goblet cells',
  '49' = 'Abs. colonocytes',
  '50' = 'Abs. colonocytes',
  '51' = 'Goblet cells',
  '52' = 'Abs. colonocytes',
  '53' = 'Goblet cells',
  '54' = 'Goblet cells',
  '55' = 'Goblet cells',
  '56' = 'Abs. colonocytes',
  '57' = 'Goblet cells',
  '58' = 'Stem cells',
  '59' = 'Goblet cells',
  '60' = 'Abs. colonocytes',
  '61' = 'TA cells',
  '62' = 'Abs. colonocytes'
  # Add/remove entries matching your actual cluster count
)

# =============================================================================
# --- PART 6: APPLY SUB-ANNOTATIONS & FINAL VISUALIZATIONS -------------------
# =============================================================================
message("\n=== STEP 6: Applying Sub-Annotations ===")

data_sub$sub_cell_types  <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony

sub_type_levels <- names(SUB_MARKERS_LIST)
sub_type_levels <- intersect(sub_type_levels,
                              as.character(unique(data_sub$sub_cell_types)))
data_sub$sub_cell_types   <- factor(data_sub$sub_cell_types,   levels = sub_type_levels)
data_sub$sub_weighted_std <- factor(data_sub$sub_weighted_std, levels = sub_type_levels)
data_sub$sub_weighted_raw <- factor(data_sub$sub_weighted_raw, levels = sub_type_levels)

summarize_annotations(data_sub, "sub_cell_types")

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

# Final DotPlot
p_final_dot <- DotPlot(data_sub, features = sub_dotplot_markers,
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
    facet_wrap(as.formula(paste("~", pgrp)), ncol = 2) +
    theme(strip.text = element_text(size = 13, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  ggsave(file.path(OUTPUT_DIR, paste0("FINAL_sub_UMAP_faceted_Colonocytes_by_", pgrp, ".png")),
         p_facet, width = 12, height = 5 * ceiling(n_levels / 2), dpi = DPI_SETTING, bg = "white")
}

# Proportional analysis
plot_cell_proportions(data_sub, "sub_cell_types", "SampleID",
                      "FINAL_sub_proportions_Colonocytes", OUTPUT_DIR)
for (group in ADDITIONAL_GROUPS_TO_PLOT) {
  plot_cell_proportions(data_sub, "sub_cell_types", group,
                        paste0("FINAL_sub_proportions_Colonocytes_by_", group), OUTPUT_DIR)
}

# Example gene expression comparison
p_pparg <- generate_gene_comparison_plots(
  seurat_obj    = data_sub,
  score_col     = "Pparg",
  group_by      = "sub_cell_types",
  x_axis        = COMPARISON_X_AXIS,
  comparisons   = COMPARISON_PAIRS,
  plot_type     = "violin",
  output_prefix = "Pparg_expression_Colonocytes_",
  plot_title    = "Pparg Expression — Colonocyte Sub-types",
  y_label       = "Log-Normalized Expression",
  fig_width     = 14, fig_height = 8
)

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
