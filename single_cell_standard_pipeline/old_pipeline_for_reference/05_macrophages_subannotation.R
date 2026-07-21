# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 5: Macrophages CELL SUB-ANNOTATION
# Version: 1.0 (CSV-Driven, Seurat Wrappers, Harmony + Standard Clustering)
#
# PURPOSE:
#   Loads the globally annotated Seurat object from Script 02 and performs
#   high-resolution sub-clustering of the T cell compartment. Provides:
#     1. Subset extraction with fresh HVG selection, PCA, Harmony, and UMAP.
#        BOTH a Harmony-corrected and a standard (non-Harmony) embedding are
#        produced for side-by-side comparison.
#     2. Weighted pre-scoring against T cell sub-type markers (CSV-driven).
#        Both z-scored (standardized) and raw (non-standardized) variants run
#        automatically — review the Top-5 report before manual annotation.
#     3. Manual sub-annotation via SUB_ANNOTATION_MAP (Action 5).
#     4. Compositional analysis (proportions by SampleID and group).
#     5. Gene expression violin/bar plots.
#
# MARKER SYSTEM (cell_type_markers.csv):
#   This script reads rows where tier == "sub" AND parent_cell_type == "T cells".
#   The CSV should contain entries for:
#     CD4+ T cells, CD8+ T cells, Tregs, NK cells, NKT cells,
#     gamma-delta T cells, Exhausted T cells, etc.
#   If no sub-rows exist for "T cells" in the CSV the script will warn and
#   fall back to the DEFAULT_SUB_MARKERS defined below.
#
# HOW TO USE:
#   1. Set paths/parameters in Part 1.
#   2. Run through Part 4 — review SUBCLUSTER_01 and SUBCLUSTER_02 PNGs.
#   3. Fill in SUB_ANNOTATION_MAP in Part 5 (Action 5).
#   4. Run the rest for final plots and save.
#
# INPUT:  {PROJECT_NAME}_final_annotated.rds   (output of 02_global_annotation.R)
# OUTPUT: {PROJECT_NAME}_T_cells_subclustered.rds
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
PARENT_CELL_TYPE <- "Macrophages"

# --- 1.3: Sub-Clustering Parameters ------------------------------------------
SUBCLUSTER_N_HVG       <- 2000   # HVGs for sub-clustering PCA
SUBCLUSTER_N_PCS       <- 50     # PCs used for kNN graph
SUBCLUSTER_K_NEIGHBORS <- 30     # k for kNN
SUBCLUSTER_MIN_DIST    <- 0.2    # UMAP min.dist
# Resolution: set to NULL to read from CSV (subcluster_resolution column),
# or override with a number here.
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

# Default fallback sub-markers for T cells (used if CSV has no "T cells" sub rows)
DEFAULT_SUB_MARKERS <-NULL
SUB_MARKERS_LIST    <- DEFAULT_SUB_MARKERS
sub_dotplot_markers <- unique(unlist(SUB_MARKERS_LIST))
SUBCLUSTER_RESOLUTION_FINAL <- if (!is.null(SUBCLUSTER_RESOLUTION)) SUBCLUSTER_RESOLUTION else 1.5

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
                   "  tier=sub, parent_cell_type=T cells, cell_type=<subtype>, markers=<pipe-sep>"))
  }
} else {
  warning(paste("[WARN] Marker CSV not found at:", MARKERS_CSV_FILE,
                "\nUsing built-in default T cell sub-markers."))
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

  W <- t(W_final)
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

# =============================================================================
# --- PART 4: SUB-CLUSTERING --------------------------------------------------
# =============================================================================
message(paste("\n=== STEP 4: Sub-Clustering of:", PARENT_CELL_TYPE, "==="))

# Ensure CellType column is set (alias from broad_cell_types in Script 02)
if (!"CellType" %in% colnames(data@meta.data)) {
  if ("broad_cell_types" %in% colnames(data@meta.data)) {
    data$CellType <- data$broad_cell_types
  } else {
    stop("[ERROR] Neither 'CellType' nor 'broad_cell_types' found. Run Script 02 first.")
  }
}

# Rename cDCs to Macrophages in both CellType layers so extraction works
data$CellType <- as.character(data$CellType)
data$broad_cell_types <- as.character(data$broad_cell_types)

data$CellType[data$CellType == "cDCs"] <- "Macrophages"
data$broad_cell_types[data$broad_cell_types == "cDCs"] <- "Macrophages"

# Verify
table(data$CellType == "Macrophages")

PARENT_CELL_TYPE <- "Macrophages"

n_parent <- sum(data@meta.data$CellType == PARENT_CELL_TYPE, na.rm = TRUE)
if (n_parent < 50) stop(paste("[ERROR] Only", n_parent, PARENT_CELL_TYPE,
                              "cells found. Check that PARENT_CELL_TYPE matches exactly."))
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
Idents(data_sub) <- "clusters_harmony"

# --- 4.1: Harmony vs Standard Comparison ---
p_harm_comp <- (
  DimPlot(data_sub, reduction = "umap_none",    group.by = "SampleID") + ggtitle("Standard PCA (No Harmony)")
) + (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "SampleID") + ggtitle("Harmony Corrected")
) & theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_01_Harmony_comparison_Macrophages.png")),
       p_harm_comp, width = 16, height = 8, dpi = DPI_SETTING)

# --- 4.2: Cluster UMAPs (annotate from these) ---
p_clust_none <- DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_none",
                         label = TRUE) + NoLegend() + ggtitle("Clusters: Standard PCA")
p_clust_harm <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony",
                         label = TRUE) + NoLegend() + ggtitle("Clusters: Harmony")
p_clust_comp <- p_clust_none + p_clust_harm
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_COMPARE_UMAP_Macrophages.png"),
       p_clust_comp, width = 16, height = 8, dpi = DPI_SETTING)
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_THIS_UMAP_Macrophages.png"),
       p_clust_harm, width = 10, height = 8, dpi = DPI_SETTING)

# --- 4.3: DotPlot for cluster characterization ---
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers,
                      group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11),  cluster.idents = T) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_Macrophages.png"),
       p_dot_sub, width = 14, height = 9, dpi = DPI_SETTING, bg = "white")

# =============================================================================
# --- PART 5: WEIGHTED PRE-SCORING -------------------------------------------
# =============================================================================
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
           file.path(OUTPUT_DIR, "SUBCLUSTER_PRESCORE_top5_Macrophages_std.xlsx"))
write_xlsx(sub_results_raw$top5_report,
           file.path(OUTPUT_DIR, "SUBCLUSTER_PRESCORE_top5_Macrophages_raw.xlsx"))

p_prescore <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Raw")
)
ggsave(file.path(OUTPUT_DIR, "SUBCLUSTER_03_prescore_comparison_Macrophages.png"),
       p_prescore, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

# =============================================================================
# --- ACTION 5: FILL IN SUB-CLUSTER ANNOTATIONS HERE ---
# =============================================================================
# Review:
#   SUBCLUSTER_02_ANNOTATE_THIS_UMAP_T_cells.png    — cluster layout
#   SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_T_cells.png — marker expression
#   SUBCLUSTER_03_prescore_comparison_T_cells.png   — weighted pre-scores
#   SUBCLUSTER_PRESCORE_top5_T_cells_std.xlsx       — top-5 ranked cell types per cluster
#
# Map each cluster number (as a string) to a T cell sub-type name.
# Available sub-types from CSV (or defaults if CSV had no entries):
#   "CD4+ T cells", "CD8+ T cells", "Tregs", "NK cells",
#   "NKT cells", "gdT cells", "Exhausted T", "ILC2s"
#FeaturePlot(data_sub, reduction = "umap_harmony", features = "scDblFinder_score")

# 'cDCs', # Itgax+ Clec9a+ → cDCs (no Cd68)
# 'TAM resident',#  Folr2+ Timd4+  C1q1/b+ and Retnla- Arg1-
# 'Monocytes', # Cd14+ Ccr2+ Ly6c2+ → Monocytes (no Timd4 no Adgre1)
# 'pDCs', # Siglech+ Bst2+  → pDCs (ignore everything else)
SUB_ANNOTATION_MAP <- c(
  '0'  = 'TAM resident',
  '1'  = 'M0 Macrophages', # ?
  '2'  = 'TAM resident',
  '3'  = 'M0 Macrophages', # ?
  '4'  = 'M0 Macrophages', # M1 macrophages
  '5'  = 'TAM resident', # or M2 macrophage
  '6'  = 'cDC1',
  '7'  = 'cDC1',
  '8'  = 'TAM resident',
  '9'  = 'Cyc. cDC1',
  '10' = 'M0 Macrophages',
  '11' = 'cDC2',
  '12' = 'Monocytes',
  '13' = 'cDC2',
  '14' = 'TAM SPP1', # or M1 macrophages or TAM inflammatory
  '15' = 'mRegDCs',
  '16' = 'Doublet', #?
  '17' = 'Monocytes',
  '18' = 'Doublet', #?
  '19' = 'M1 Macrophages', # or TAM IFN
  '20' = 'Cyc. Macrophages',
  '21' = 'Doublet', # or TAM IFN
  '22' = 'Monocytes',
  '23' = 'Neutrophils',
  '24' = 'Doublet',
  '25' = 'M2 Macrophages', # or M0 macrophages
  '26' = 'TAM resident',
  '27' = 'Doublet', # or unlikely M0 macrophages or TAM SPP1? because Gpnmb+?
  '28' = 'TAM IFN',
  '29' = 'Doublet', # or M2 macrophages
  '30' = 'pDCs',
  '31' = 'Doublet' # Everything is weak here...
)

# =============================================================================
# --- PART 6: APPLY SUB-ANNOTATIONS & FINAL VISUALIZATIONS -------------------
# =============================================================================
message("\n=== STEP 6: Applying Sub-Annotations ===")
#data_sub_bak <- data_sub
Idents(data_sub) <- "clusters_harmony"

all_markers <- FindAllMarkers(
  data_sub,
  min.pct        = 0.05,
  logfc.threshold = 0.25,
  only.pos       = TRUE
)

# Top 30 per cluster by log2FC
top30 <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100) %>%
  ungroup()

write_xlsx(top30, file.path(OUTPUT_DIR, "macrophage_top50_markers_per_cluster.xlsx"))
message("Saved: macrophage_top30_markers_per_cluster.xlsx")


# --- Contaminant marker FeaturePlots ---
# 16 = Fibroblasts, 18 = SMCs, 21 = T cells, 24 = Mast cells
# 25 = Plasma B cells, 27 = Neurons, 29 = VECs, 30 = T cells, 31 = B cells

contaminant_markers <- c(
  "Col1a1",  # 16 Fibroblasts
  "Myh11",   # 18 SMCs
  "Cd3e",    # 21 T cells
  "Cpa3",    # 24 Mast cells
  "Mzb1",    # 25 Plasma B cells
  "Slc1a2",  # 27 Neurons
  "Cldn5",   # 29 VECs
  "Cd8b1",   # 30 T cells/lymphoid
  "Cd19"     # 31 B cells
)

p2 <- FeaturePlot(data_sub, 
                  features = contaminant_markers,
                  reduction = "umap_harmony",
                  ncol = 3) &
  theme(plot.title = element_text(size = 10))
p2
ggsave(file.path(OUTPUT_DIR, "MAC_02_contaminant_markers_featureplot.png"),
       p2, width = 12, height = 10, dpi = DPI_SETTING, bg = "white")

message("Saved: MAC_01_annotation_umap.png")
message("Saved: MAC_02_contaminant_markers_featureplot.png")

p2<- FeaturePlot(data_sub,
            features = c(
              "Hexb",    # cluster 1  - TAM_resident/M2
              "Cd209a",  # cluster 11 - cDC2
              "Siglecg", # cluster 13 - pDC/regulatory
              "Spp1",    # cluster 14 - TAM_inflammatory
              "Ccr7",    # cluster 15 - Migratory cDCs
              "Ccr2",    # cluster 17 - Monocytes
              "Ccl3",    # cluster 19 - M1
              "Bub1",    # cluster 20 - Cyc. Macrophages
              "S100a8",  # cluster 23 - Neutrophils
              "Tmem119",   # cluster 26 - TAM_resident
              "Ifit1"    # cluster 28 - TAM_IFN
            ),
            reduction = "umap_harmony",
            ncol = 4) &
  theme(plot.title = element_text(size = 10))
p2
ggsave(file.path(OUTPUT_DIR, "MAC_02_new_markers_featureplot2.png"),
       p2, width = 12, height = 10, dpi = DPI_SETTING, bg = "white")



data_sub$sub_cell_types <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony
#data_sub_bak <- data_sub

# The find markers below and feature plots helped to determine the doublets of those cells
data_sub <- subset(data_sub, subset = sub_cell_types != "Doublet")

DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
        label = TRUE, repel = TRUE)


# 
# # Load the necessary libraries
# library(RColorBrewer)
# library(ggplot2)
# library(patchwork)
# 
# # 1. Plot nCount_RNA (Total molecules)
# p1 <- FeaturePlot(data_sub, features = "nCount_RNA", reduction = "umap_harmony",  pt.size = 0.5) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
#   ggtitle("Total RNA Counts (nCount)")
# 
# # 2. Plot nFeature_RNA (Unique genes)
# p2 <- FeaturePlot(data_sub, features = "nFeature_RNA", reduction = "umap_harmony", pt.size = 0.5) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
#   ggtitle("Unique Genes (nFeature)")
# 
# # 3. View them side-by-side
# p1 | p2

sub_type_levels_macrophages <- c(
  "pDCs",
  "cDC2",
  "mRegDCs",
  "cDC1", 
  "Cyc. cDC1",
  "Monocytes",
  "Cyc. Macrophages",
  "M0 Macrophages",
  "TAM SPP1",
  "M2 Macrophages",
  "TAM resident",
  "M1 Macrophages",
  "TAM IFN",
  "Neutrophils"
)

data_sub$sub_cell_types <- factor(data_sub$sub_cell_types, levels = sub_type_levels_macrophages)

# --- Final UMAP: manual annotation ---
p_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) +
  ggtitle(paste("Sub-Cell Types:", PARENT_CELL_TYPE)) +
  theme(legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(OUTPUT_DIR, "FINAL_SUBCLUSTER_UMAP_Macrophages.png"),
       p_final, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# --- Comparison: manual vs weighted ---
p_compare <- DimPlot(data_sub, reduction = "umap_harmony",
                     group.by = c("sub_cell_types", "sub_weighted_std", "sub_weighted_raw"),
                     label = FALSE, repel = TRUE)
ggsave(file.path(OUTPUT_DIR, "FINAL_SUBCLUSTER_comparison_Macrophages.png"),
       p_compare, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")


# =============================================================================
# 1. DEFINE LEVELS AND MAPPING OF DOTPLOT
# =============================================================================
# This is the order you want for your plot axes/legend
sub_type_levels2 <- sub_type_levels_macrophages
# This maps your Factor Levels (left) to the SUB_MARKERS_LIST names (right)
# This handles the "Cyc. CD4+ T cells" -> "Cyc. T cells" mismatch
name_mapping <- c(
  "pDCs"             = "pDCs",
  "cDC1"             = "cDC1",
  "Cyc. cDC1"        = "cDC1",
  "cDC2"             = "cDC2",
  "mRegDCs"          = "mRegDCs",
  "Monocytes"        = "Monocytes",
  "M0 Macrophages"   = "M0 macrophages",
  "M2 Macrophages"   = "M2 macrophages",
  "TAM resident"     = "TAM_resident",
  "M1 Macrophages"   = "M1 macrophages",
  "TAM IFN"          = "TAM_IFN",
  "TAM SPP1"         = "TAM_SPP1",
  "Cyc. Macrophages" = "Cyc. Macrophages",
  "Neutrophils"      = "Neutrophils"
)

# Apply levels to the Seurat object
data_sub$sub_cell_types <- factor(
  data_sub$sub_cell_types, 
  levels = sub_type_levels2
)

# p_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
#                    label = TRUE, repel = TRUE) +
#   ggtitle(paste("Sub-Cell Types:", PARENT_CELL_TYPE)) +
#   theme(legend.text  = element_text(size = 12),
#         legend.title = element_text(size = 14, face = "bold")) +
#   guides(color = guide_legend(override.aes = list(size = 4)))
# 
# p_final
# =============================================================================
# 2. EXTRACT AND ORCHESTRATE GENE LIST
# =============================================================================

# A. Set the mandatory lead genes
lead_genes <- NULL

ordered_markers <- unlist(lapply(sub_type_levels2, function(lvl) {
  marker_key <- name_mapping[lvl]
  message("marker key: ", marker_key)
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

p_final_dot
ggsave(file.path(OUTPUT_DIR, "FINAL_DotPlot_sub_macrophages.png"),
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
  ggsave(file.path(OUTPUT_DIR, paste0("FINAL_sub_UMAP_faceted_Macrophages_by_", pgrp, ".png")),
         p_facet, width = 5 * ceiling(n_levels / 2), height = 10, dpi = DPI_SETTING, bg = "white")
}

# =============================================================================
# --- PART 7: SAVE ------------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving T Cell Sub-Cluster Object ===")
saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_Macrophages_subclustered.rds")))
message(paste0(
  "\n=== T CELL SUB-ANNOTATION COMPLETE ===\n",
  "  Saved: ", PROJECT_NAME, "_Macrophages_subclustered.rds\n",
  "  Sub-types found: ", paste(levels(data_sub$sub_cell_types), collapse = ", "), "\n",
  "\nNext steps:\n",
  "  - Script 04: Colonocyte sub-annotation (04_colonocyte_subannotation.R)\n",
  "  - Script 05: DE + Two-Way ANOVA (05_DE_and_two_way_ANOVA.R)\n"
))


data_sub <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_Macrophages_subclustered.rds")))

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
  ungroup() %>%
  # Fill missing cell types within each SampleID only
  complete(nesting(SampleID, Genotype_sex), sub_cell_types,
           fill = list(n = 0, Percentage = 0)) %>%
  ungroup()

# --- Step 2: Plot ---
p_stats <- ggplot(df_pct, aes(x = Genotype_sex, y = Percentage, fill = Genotype_sex)) +
  stat_summary(fun = "mean", geom = "bar", width = 0.7, color = "black", alpha = 0.8) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_jitter(shape = 21, color = "black", size = 2, width = 0.15) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",
                     label = "p.signif", size = 4) +
  facet_wrap(~ sub_cell_types, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
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
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + # Keeps Females on top, Males on bottom
  theme(
    legend.position = "top",
    legend.justification = "left",   # This snaps the legend to the left side
    legend.box.just = "left"
  )

p_stats

# --- Step 3: Save ---
#dir.create("proportion_analysis", showWarnings = FALSE)
savef <- paste0(OUTPUT_DIR, "/macrophages_stats_cell_props.png")
ggsave(savef,p_stats, width = 12, height = 10, dpi=DPI_SETTING)
savef <- paste0(OUTPUT_DIR, "/macrophages_cell_props_long_data.xlsx")
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
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_stats_plot_macrophages.png"),
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
  save_path     = file.path(OUTPUT_DIR, "Nr4a1_custom_stats_plot_macrophages.png"),
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
  save_path     = file.path(OUTPUT_DIR, "Nr4a2_stats_plot_macrophages.png"),
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
  save_path     = file.path(OUTPUT_DIR, "Nr4a3_stats_plot_macrophages.png"),
  dpi           = DPI_SETTING
)



library(Seurat)
library(AUCell)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

message("\n##### Starting Targeted Macrophage Functional Scoring #####")
data_sub <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_Macrophages_subclustered.rds")))

# =============================================================================
# ENVIRONMENT PATHS
# =============================================================================
#PROJECT_NAME <- "Nr4a1_s17_ack"
#ROOT_PATH    <- "/home/ssromerogon/2026_nr4a1_ack/r_process"
#OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# 1. SUBSETTING (Isolating Macrophage Continuum from Myeloid Object)
# =============================================================================
mac_lineages <- c("Monocytes", "Cyc. Macrophages", "M0 Macrophages", 
                  "TAM SPP1", "M2 Macrophages", "TAM resident", 
                  "M1 Macrophages", "TAM IFN")

data_mac <- subset(data_sub, subset = sub_cell_types %in% mac_lineages)

# =============================================================================
# 2. DEFINE MACROPHAGE FUNCTIONAL GENE SETS
# =============================================================================
mac_gene_sets <- list(
  M1_ProInflammatory   = c("Stat1", "Irf1", "Nos2", "Irf5", "Il6", "Cxcl10", "Cxcl9", "Tnf", "Il1b"),
  M2_Immunosuppressive = c("Arg1", "Mrc1", "Cd163", "Chil3", "Retnla", "C1qa", "C1qb", "C1qc"),
  TAM_SPP1_Remodeling  = c("Spp1", "Fn1", "Vcan", "Gpnmb", "Lgals3", "Apoc1")
)

# =============================================================================
# 3. AUCell PIPELINE (Seurat v5 Layer & 10% Threshold)
# =============================================================================
expr_matrix_mac <- GetAssayData(object = data_mac, assay = "RNA", layer = "counts")
rankings_mac    <- AUCell_buildRankings(expr_matrix_mac, plotStats = FALSE)

valid_mac_sets <- list()
for (set_name in names(mac_gene_sets)) {
  genes_in_set  <- mac_gene_sets[[set_name]]
  present_genes <- genes_in_set[genes_in_set %in% rownames(expr_matrix_mac)]
  if (length(present_genes) > 0) valid_mac_sets[[set_name]] <- present_genes
}

auc_mac    <- AUCell_calcAUC(valid_mac_sets, rankings_mac, aucMaxRank = nrow(rankings_mac) * 0.10)
matrix_mac <- getAUC(auc_mac)

for (score_name in rownames(matrix_mac)) {
  data_mac[[score_name]] <- as.numeric(matrix_mac[score_name, ])
}

# =============================================================================
# 4. PLOTTING LOOP WITH REFINED TITLES & LARGE FONTS
# =============================================================================
custom_levels <- c(
  "WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
  "WT_Male", "Polyp_Male",   "Polyp_NR4a1_KO_Male"
)

COMPARISON_PAIRS <- list(
  c("Polyp_Female",    "Polyp_NR4a1_KO_Female"),
  c("Polyp_Male",      "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_NR4a1_KO_Female"),
  c("WT_Male",         "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_Female"),
  c("WT_Male",         "Polyp_Male")
)

scores_to_plot <- c("M1_ProInflammatory", "M2_Immunosuppressive", "TAM_SPP1_Remodeling")

for (score_col in scores_to_plot) {
  message("Generating plot for: ", score_col)
  
  plot_df <- data_mac@meta.data %>%
    select(all_of(c(score_col, "Genotype_sex", "sub_cell_types"))) %>%
    drop_na()
  
  plot_df$Genotype_sex <- factor(plot_df$Genotype_sex, levels = custom_levels)
  
  # Format specific unified title layout requested
  clean_title <- gsub("_", " ", score_col)
  display_title <- paste("Macrophage", clean_title, "Across Subpopulations")
  
  violin_p <- ggplot(plot_df, aes(x = Genotype_sex, y = .data[[score_col]], fill = Genotype_sex)) +
    geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    
    stat_compare_means(
      comparisons = COMPARISON_PAIRS,
      method = "wilcox.test",
      method.args = list(exact = FALSE, correct = TRUE),
      label = "p.signif", 
      step.increase = 0.05, 
      tip.length = 0.01
    ) +
    
    # Grid layout grouped by your specific macrophage UMAP clusters
    facet_wrap(~ sub_cell_types, scales = "free_y", ncol = 2) +
    
    labs(
      title = display_title,
      y = "Cell Score (AUCell)"
    ) +
    theme_classic() + 
    theme(
      plot.title   = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 14),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
      strip.text   = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )
  
  output_filename <- paste0("Macrophage_", score_col, "_Explicit_Groups.png")
  
  ggsave(
    filename = file.path(OUTPUT_DIR, output_filename), 
    plot = violin_p, 
    device = "png",
    dpi = 300,
    height = 14,    # Increased height slightly to accommodate the 8 faceted clusters cleanly
    width = 12, 
    units = "in"
  )
  
  print(violin_p)
}

# =============================================================================
# GLOBAL MACROPHAGE SCORING (NO SUBPOPULATION SPLITTING)
# =============================================================================
message("\n##### Generating Global Macrophage Overview Plots #####")

# Ensure your factor levels and pairs are set
custom_levels <- c(
  "WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
  "WT_Male", "Polyp_Male",   "Polyp_NR4a1_KO_Male"
)

COMPARISON_PAIRS <- list(
  c("Polyp_Female",    "Polyp_NR4a1_KO_Female"),
  c("Polyp_Male",      "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_NR4a1_KO_Female"),
  c("WT_Male",         "Polyp_NR4a1_KO_Male"),
  c("WT_Female",       "Polyp_Female"),
  c("WT_Male",         "Polyp_Male")
)

scores_to_plot <- c("M1_ProInflammatory", "M2_Immunosuppressive", "TAM_SPP1_Remodeling")

for (score_col in scores_to_plot) {
  message("Generating global plot for: ", score_col)
  
  # 1. Pull metadata for all macrophages together
  plot_df <- data_mac@meta.data %>%
    select(all_of(c(score_col, "Genotype_sex"))) %>%
    drop_na()
  
  plot_df$Genotype_sex <- factor(plot_df$Genotype_sex, levels = custom_levels)
  
  clean_title <- gsub("_", " ", score_col)
  display_title <- paste("Global Macrophage", clean_title, "Profile")
  
  # 2. Build the unified plot (No facet_wrap here)
  global_p <- ggplot(plot_df, aes(x = Genotype_sex, y = .data[[score_col]], fill = Genotype_sex)) +
    geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    
    stat_compare_means(
      comparisons = COMPARISON_PAIRS,
      method = "wilcox.test",
      method.args = list(exact = FALSE, correct = TRUE),
      label = "p.signif", 
      step.increase = 0.06, 
      tip.length = 0.01
    ) +
    
    labs(
      title = display_title,
      y = "Cell Score (AUCell)"
    ) +
    theme_classic() + 
    theme(
      plot.title   = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 14),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
      legend.position = "none"
    )
  
  # 3. Save as a cleaner, single-panel shape (smaller dimensions since there are no facets)
  output_filename <- paste0("Global_Macrophage_", score_col, "_Overview.png")
  
  ggsave(
    filename = file.path(OUTPUT_DIR, output_filename), 
    plot = global_p, 
    device = "png",
    dpi = 300,
    height = 7,    
    width = 8, 
    units = "in"
  )
  
  print(global_p)
}
