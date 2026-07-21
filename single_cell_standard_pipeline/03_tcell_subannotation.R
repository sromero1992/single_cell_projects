# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 3: T CELL SUB-ANNOTATION
# Version: 1.0 (CSV-Driven, Seurat Wrappers, Harmony + Standard Clustering)
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
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
#ROOT_PATH <- "/home/ssromerogon/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"  # Windows
ROOT_PATH <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"

OUTPUT_DIR       <- file.path(ROOT_PATH, "seurat_output")
TCELL_DIR        <- file.path(OUTPUT_DIR, "tcell_subannotation")   # <-- all T cell plots go here
MAIN_RDS         <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))
MARKERS_CSV_FILE <- file.path(ROOT_PATH, "cell_type_markers.csv")

# --- 1.2: Target Cell Type ---------------------------------------------------
# Must exactly match a label in the broad_cell_types / CellType column from 02.
PARENT_CELL_TYPE <- "T cells"

# --- 1.3: Sub-Clustering Parameters ------------------------------------------
SUBCLUSTER_N_HVG       <- 2000   # HVGs for sub-clustering PCA
SUBCLUSTER_N_PCS       <- 50     # PCs used for kNN graph
SUBCLUSTER_K_NEIGHBORS <- 15     # k for kNN
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
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(TCELL_DIR))  dir.create(TCELL_DIR,  recursive = TRUE)
message(paste("  T cell output folder:", TCELL_DIR))
data <- readRDS(MAIN_RDS)
message(paste("  Loaded:", ncol(data), "cells"))

# --- Parse markers from CSV --------------------------------------------------
parse_markers <- function(marker_string) {
  strsplit(trimws(marker_string), "\\|")[[1]]
}

# Default fallback sub-markers for T cells (used if CSV has no "T cells" sub rows)
DEFAULT_SUB_MARKERS <- list(
  "CD4+ T cells"     = c("Cd4"),
  "CD8+ T cells"     = c("Cd8a", "Cd8b1"),
  "Tregs"            = c("Foxp3", "Il2ra", "Ikzf2", "Tigit"),
  "NK cells"         = c("Ncr1", "Klrb1c", "Nkg7", "Klrd1", "Xcl1"),
  "NKT cells"        = c("Cd3e", "Klrb1c", "Cd8a"),
  "γδ T cells"        = c("Trdc", "Trgc1", "Trgc2"),
  "Exhausted T"      = c("Pdcd1", "Lag3", "Havcr2", "Tigit", "Entpd1"),
  "ILC2s"            = c("Il1rl1", "Gata3", "Il13", "Il5")
)

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

# Convenience negation of %in% (used in dirty-cluster removal)
`%nin%` <- Negate(`%in%`)

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


run_pca_umap <- function(data,
                         num_hvg  = 2000,
                         dims_pca = 50,
                         min_dist = 0.3,
                         kneigh   = 15) {
  
  message("  Running HVG selection, scaling, and PCA...")
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = num_hvg) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = dims_pca, reduction.name = "pca", verbose = FALSE)
  gc()
  
  # Track A: Standard PCA → UMAP
  message("  Running UMAP on PCA...")
  data <- RunUMAP(data, dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh,
                  min.dist = min_dist, n.epochs = 500,
                  reduction.name = "umap_none", verbose = FALSE)
  gc()
  
  # Track B: Harmony → UMAP
  message("  Running Harmony...")
  data <- RunHarmony(data, group.by.vars = "SampleID",
                     reduction = "pca", reduction.save = "harmony", verbose = FALSE)
  gc()
  
  message("  Running UMAP on Harmony...")
  data <- RunUMAP(data, dims = 1:dims_pca, reduction = "harmony", n.neighbors = kneigh,
                  min.dist = min_dist, n.epochs = 500,
                  reduction.name = "umap_harmony", verbose = FALSE)
  gc()
  
  return(data)
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
ggsave(file.path(TCELL_DIR, paste0("SUBCLUSTER_01_Harmony_comparison_T_cells.png")),
       p_harm_comp, width = 16, height = 8, dpi = DPI_SETTING)

# --- 4.2: Cluster UMAPs (annotate from these) ---
p_clust_none <- DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_none",
                        label = TRUE) + NoLegend() + ggtitle("Clusters: Standard PCA")
p_clust_harm <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony",
                        label = TRUE) + NoLegend() + ggtitle("Clusters: Harmony")
p_clust_comp <- p_clust_none + p_clust_harm
ggsave(file.path(TCELL_DIR, "SUBCLUSTER_02_COMPARE_UMAP_T_cells.png"),
       p_clust_comp, width = 16, height = 8, dpi = DPI_SETTING)
ggsave(file.path(TCELL_DIR, "SUBCLUSTER_02_ANNOTATE_THIS_UMAP_T_cells.png"),
       p_clust_harm, width = 10, height = 8, dpi = DPI_SETTING)

# --- 4.3: DotPlot for cluster characterization ---
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers,
                     group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11),  cluster.idents = T) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(TCELL_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_T_cells.png"),
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
           file.path(TCELL_DIR, "SUBCLUSTER_PRESCORE_top5_T_cells_std.xlsx"))
write_xlsx(sub_results_raw$top5_report,
           file.path(TCELL_DIR, "SUBCLUSTER_PRESCORE_top5_T_cells_raw.xlsx"))

p_prescore <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Raw")
)
ggsave(file.path(TCELL_DIR, "SUBCLUSTER_03_prescore_comparison_T_cells.png"),
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
# Doublet score column depends on DOUBLET_METHOD used in Script 01:
#   DoubletFinder -> "DF_score"   |   scDblFinder -> "scDblFinder_score"
#FeaturePlot(data_sub, reduction = "umap_harmony", features = "DF_score")

SUB_ANNOTATION_MAP <- c(
  '0'  = 'ILC2',
  '1'  = 'CD8+ T cells',
  '2'  = 'Tregs',
  '3'  = 'γδ T cells',
  '4'  = 'Tregs',
  '5'  = 'Cyc. CD4+ T cells',
  '6'  = 'γδ T cells', 
  '7'  = 'CD4+ T cells', 
  '8'  = 'CD8+ T cells',
  '9'  = 'CD4+ T cells', 
  '10' = 'NK cells', # Not ILC1, it is Prf1+
  '11' = 'CD8+ T cells',
  '12' = 'Tregs', # ??
  '13' = 'Tregs',
  '14' = 'Cyc. CD4+ T cells',
  '15' = 'ILC3', 
  '16' = 'Tregs', 
  '17' = 'CD4+ T cells',
  '18' = 'Tregs',
  '19' = 'CD4+ T cells', 
  '20' = 'CD8+ T cells', 
  '21' = 'CD8+ T cells', 
  '22' = 'Doublet', # ??
  '23' = 'γδ T cells',
  '24' = 'Tregs',
  '25' = 'Doublet', 
  '26' = 'Doublet' # ??
)

# =============================================================================
# --- PART 6: FIRST-PASS ANNOTATION & DIRTY-CLUSTER QC -----------------------
# =============================================================================
message("\n=== STEP 6a: First-Pass Sub-Annotation (pre-cleaning) ===")

# ACTION 5 annotation map is applied here (see above).
# Label anything you are unsure about as "Doublet" or "Unknown" — those get
# removed in the cleaning step below.
data_sub$sub_cell_types  <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony

# --- 6a.1: Pre-cleaning overview UMAP (dirty, includes flagged clusters) ---
p_dirty <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) +
  ggtitle("Pre-Cleaning: All Clusters Labelled (including Doublet/Unknown)") +
  theme(legend.text = element_text(size = 11))
ggsave(file.path(TCELL_DIR, "CLEANING_00_pre_removal_overview.png"),
       p_dirty, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# --- 6a.2: QC FeaturePlots for suspicious clusters -------------------------
# Edit the cluster numbers in the FindMarkers calls and the feature lists here
# as you discover new dirty clusters. Save them all to TCELL_DIR for records.

message("  Generating QC FeaturePlots for dirty cluster inspection...")

# Broad identity QC
p_qc1 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = c("Cd3e", "Cd3g", "Cd8a", "Cd8b1", "Cd4", "Foxp3", "Mki67", "nCount_RNA"),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(TCELL_DIR, "CLEANING_01_QC_broad_identity_markers.png"),
       p_qc1, width = 16, height = 8, dpi = DPI_SETTING, bg = "white")

# Contamination / non-T markers
p_qc2 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = c("Muc2", "Krt20", "Epcam", "Ptprc", "Cd19", "Ms4a1"),
                     ncol = 3) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(TCELL_DIR, "CLEANING_02_QC_contamination_markers.png"),
       p_qc2, width = 14, height = 10, dpi = DPI_SETTING, bg = "white")

# NK / ILC disambiguation
nk_markers   <- c("Eomes", "Gnly", "Nkg7", "Prf1", "Gzmb", "Fgfbp2", "Klrd1")
ilc1_markers <- c("Tbx21", "Ifng", "Il7r", "Cxcr6", "Itga1", "Cd69", "Tnfsf10", "Znf683")
nk_ilc_features <- intersect(c(nk_markers, ilc1_markers), rownames(data_sub))
p_qc3 <- FeaturePlot(data_sub, features = nk_ilc_features, reduction = "umap_harmony", ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(TCELL_DIR, "CLEANING_03_NK_ILC1_disambiguation.png"),
       p_qc3, width = 16, height = 12, dpi = DPI_SETTING, bg = "white")

# ILC / Th17 markers
p_qc4 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Rorc", "Rora", "Il17a", "Il22", "Prf1",
                                            "Ncr1", "Tbx21", "Rorc", "Stat3", "Batf",
                                            "Irf4", "Ncam1", "Klrb1", "Il1r1", "Klrd1"),
                                          rownames(data_sub)),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(TCELL_DIR, "CLEANING_04_ILC_Th17_markers.png"),
       p_qc4, width = 16, height = 12, dpi = DPI_SETTING, bg = "white")

# Exhaustion markers
p_qc5 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Eomes", "Prdm1", "Pdcd1", "Havcr2", "Entpd1", "Lag3"),
                                          rownames(data_sub)),
                     ncol = 3) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(TCELL_DIR, "CLEANING_05_exhaustion_markers.png"),
       p_qc5, width = 14, height = 10, dpi = DPI_SETTING, bg = "white")


# --- 6a.3: Per-suspect-cluster FindMarkers (edit cluster IDs as needed) -----
Idents(data_sub) <- "clusters_harmony"
message("  Running FindAllMarkers on all clusters...")
all_cluster_markers <- FindAllMarkers(
  data_sub,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

top30_by_cluster <- all_cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 30) %>%
  dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::ungroup()

write_xlsx(
  as.data.frame(top30_by_cluster),
  file.path(TCELL_DIR, "CLEANING_top30_per_cluster.xlsx")
)

message(paste("  Marker table saved —",
              length(unique(all_cluster_markers$cluster)), "clusters,",
              nrow(top30_by_cluster), "rows total."))

# =============================================================================
# --- PART 6b: REMOVE DIRTY CLUSTERS & RE-EMBED ----------------------------
# =============================================================================
message("\n=== STEP 6b: Removing dirty clusters and re-running UMAP ===")

# ACTION: set the label(s) used for dirty clusters in SUB_ANNOTATION_MAP above.
# Default is "Doublet" — add "Unknown" or others if you used them.
DIRTY_LABELS <- c("Doublet", "Unknown")   # <-- ACTION: edit if you used different labels

n_before <- ncol(data_sub)
data_sub  <- subset(data_sub, subset = sub_cell_types %nin% DIRTY_LABELS)
n_after   <- ncol(data_sub)
message(paste("  Removed", n_before - n_after, "dirty cells.",
              n_after, "clean cells remain."))

# Drop used reductions so run_pca_umap starts fresh
data_sub@reductions <- list()
data_sub@graphs     <- list()

message("  Re-running HVG / PCA / Harmony / UMAP on clean cells...")
data_sub <- run_pca_umap(
  data     = data_sub,
  num_hvg  = SUBCLUSTER_N_HVG,
  dims_pca = SUBCLUSTER_N_PCS,
  min_dist = SUBCLUSTER_MIN_DIST,
  kneigh   = SUBCLUSTER_K_NEIGHBORS
)
gc()

# Re-cluster on the clean embedding so cluster IDs are stable for re-annotation
SUBCLUSTER_N_PCS <- 20
message("  Re-clustering clean cells on harmony embedding...")
data_sub <- FindNeighbors(data_sub, dims = 1:SUBCLUSTER_N_PCS, reduction = "harmony",
                          k.param = SUBCLUSTER_K_NEIGHBORS, graph.name = "harmony_nn",
                          verbose = FALSE) %>%
  FindClusters(resolution = SUBCLUSTER_RESOLUTION_FINAL, graph.name = "harmony_nn",
               cluster.name = "clusters_harmony_clean", verbose = FALSE)
Idents(data_sub) <- "clusters_harmony_clean"
gc()

# --- 6b.1: Save clean cluster UMAPs for re-annotation review ---------------
p_clean_clust <- DimPlot(data_sub, reduction = "umap_harmony",
                         group.by = "clusters_harmony_clean",
                         label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("CLEAN Clusters: Harmony (re-embedded)")
ggsave(file.path(TCELL_DIR, "CLEANING_06_clean_clusters_ANNOTATE_THIS.png"),
       p_clean_clust, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

p_clean_both <- (
  DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_harmony_clean",
          label = TRUE) + NoLegend() + ggtitle("Clean — Standard PCA")
) + (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony_clean",
          label = TRUE) + NoLegend() + ggtitle("Clean — Harmony")
)
ggsave(file.path(TCELL_DIR, "CLEANING_07_clean_compare_both_embeddings.png"),
       p_clean_both, width = 16, height = 8, dpi = DPI_SETTING)

# --- 6b.2: Re-run weighted pre-scoring on clean clusters -------------------
message("  Re-running weighted pre-scoring on clean clusters...")
sub_results_clean_std <- get_weighted_annotation(
  seurat_obj             = data_sub,
  marker_genes           = SUB_MARKERS_LIST,
  cluster_key            = "clusters_harmony_clean",
  standardize_expression = TRUE
)
sub_results_clean_raw <- get_weighted_annotation(
  seurat_obj             = data_sub,
  marker_genes           = SUB_MARKERS_LIST,
  cluster_key            = "clusters_harmony_clean",
  standardize_expression = FALSE
)
data_sub$sub_weighted_std_clean <- sub_results_clean_std$annotation_vector
data_sub$sub_weighted_raw_clean <- sub_results_clean_raw$annotation_vector

message("\n--- Clean Sub-cluster Pre-scoring Top-5 Report (standardized) ---")
print(sub_results_clean_std$top5_report)
write_xlsx(sub_results_clean_std$top5_report,
           file.path(TCELL_DIR, "CLEANING_08_prescore_clean_std.xlsx"))
write_xlsx(sub_results_clean_raw$top5_report,
           file.path(TCELL_DIR, "CLEANING_08_prescore_clean_raw.xlsx"))

p_prescore_clean <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std_clean",
          label = TRUE, repel = TRUE) + ggtitle("Clean Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw_clean",
          label = TRUE, repel = TRUE) + ggtitle("Clean Pre-Score: Raw")
)
ggsave(file.path(TCELL_DIR, "CLEANING_09_prescore_clean_comparison.png"),
       p_prescore_clean, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

# --- 6b.3: DotPlot for clean cluster characterization ----------------------
p_dot_clean <- DotPlot(data_sub, features = sub_dotplot_markers,
                       group.by = "clusters_harmony_clean", dot.min = 0.05,
                       cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Clean Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(TCELL_DIR, "CLEANING_10_clean_ANNOTATE_USING_THIS_DOTPLOT.png"),
       p_dot_clean, width = 14, height = 9, dpi = DPI_SETTING, bg = "white")

# --- 6c: Per-suspect-cluster FindMarkers (edit cluster IDs as needed) -----
Idents(data_sub) <- "clusters_harmony_clean"
message("  Running FindAllMarkers on all clusters...")
all_cluster_markers <- FindAllMarkers(
  data_sub,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

top30_by_cluster <- all_cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 30) %>%
  dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::ungroup()

write_xlsx(
  as.data.frame(top30_by_cluster),
  file.path(TCELL_DIR, "CLEANING_top30_per_cluster_clean.xlsx")
)

message(paste("  Marker table saved —",
              length(unique(all_cluster_markers$cluster)), "clusters,",
              nrow(top30_by_cluster), "rows total."))

# =============================================================================
# --- ACTION 6: FILL IN CLEAN CLUSTER ANNOTATIONS HERE -----------------------
# =============================================================================
# Review (in TCELL_DIR):
#   CLEANING_06_clean_clusters_ANNOTATE_THIS.png     — new cluster layout
#   CLEANING_10_clean_ANNOTATE_USING_THIS_DOTPLOT.png — marker expression
#   CLEANING_08_prescore_clean_std.xlsx              — top-5 per clean cluster
#   CLEANING_09_prescore_clean_comparison.png        — automated suggestions
#
# Map each NEW cluster number (as a string) to a sub-type name.
# These will differ from the original SUB_ANNOTATION_MAP because cells were
# re-embedded and re-clustered after doublet removal.

SUB_ANNOTATION_MAP_CLEAN <- c(
  '0'  = 'ILC2', #      
  '1'  = 'CD8+ T cells', #
  '2'  = 'Tregs', #
  '3'  = 'Tregs', # 
  '4'  = 'γδ T cells', #
  '5'  = 'γδ T cells', #  
  '6'  = 'Cyc. T cells', # 
  '7'  = 'CD4+ T cells', #
  '8'  = 'CD8+ T cells', #
  '9'  = 'CD4+ T cells', #
  '10' = 'Tregs', # 
  '11' = 'NK cells', #
  '12' = 'CD8+ T cells', #
  '13' = 'CD8+ T cells', #
  '14' = 'Doublet',#  
  '15' = 'CD4+ T cells', #
  '16' = 'CD4+ T cells', #
  '17' = 'CD4+ T cells', # 
  '18' = 'CD4+ T cells', #
  '19' = 'CD4+ T cells',#
  '20' = 'Cyc. T cells', #
  '21' = 'γδ T cells', #
  '22' = 'Cyc. T cells', #
  '23' = 'CD4+ T cells', #
  '24' = 'NK cells', #
  '25' = 'Tregs'
)

# =============================================================================
# --- PART 6c: APPLY CLEAN ANNOTATIONS & FINAL VISUALIZATIONS ---------------
# =============================================================================
message("\n=== STEP 6c: Applying clean sub-annotations ===")

data_sub$sub_cell_types  <- recode_factor(data_sub$clusters_harmony_clean, !!!SUB_ANNOTATION_MAP_CLEAN)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony_clean

# Remove any remaining dirty clusters before finalizing levels
DIRTY_LABELS_CLEAN <- c("Doublet", "Unknown")   # <-- edit if you used different labels
n_before <- ncol(data_sub)
data_sub <- subset(data_sub, subset = sub_cell_types %nin% DIRTY_LABELS_CLEAN)
message(paste("  Removed", n_before - ncol(data_sub), "additional dirty cells.",
              ncol(data_sub), "final clean cells."))

sub_type_levels2 <- c(
  "CD8+ T cells", "γδ T cells", "Cyc. T cells",
  "CD4+ T cells", "Tregs",
  "NK cells", "ILC2"
)
data_sub$sub_cell_types <- factor(data_sub$sub_cell_types, levels = sub_type_levels2)


# Drop reductions and re-embed without the removed cells
data_sub@reductions <- list()
data_sub@graphs     <- list()

data_sub <- run_pca_umap(
  data     = data_sub,
  num_hvg  = SUBCLUSTER_N_HVG,
  dims_pca = SUBCLUSTER_N_PCS,
  min_dist = SUBCLUSTER_MIN_DIST,
  kneigh   = SUBCLUSTER_K_NEIGHBORS
)
gc()
# Re-cluster on the clean embedding so cluster IDs are stable for re-annotation
SUBCLUSTER_N_PCS <- 20
message("  Re-clustering clean cells on harmony embedding...")
data_sub <- FindNeighbors(data_sub, dims = 1:SUBCLUSTER_N_PCS, reduction = "harmony",
                          k.param = SUBCLUSTER_K_NEIGHBORS, graph.name = "harmony_nn",
                          verbose = FALSE) %>%
  FindClusters(resolution = SUBCLUSTER_RESOLUTION_FINAL, graph.name = "harmony_nn",
               cluster.name = "clusters_harmony_clean", verbose = FALSE)
Idents(data_sub) <- "clusters_harmony_clean"
gc()



# --- Final UMAP: clean manual annotation ---
p_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) +
  ggtitle(paste("Sub-Cell Types:", PARENT_CELL_TYPE, "(Clean)")) +
  theme(legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(TCELL_DIR, "FINAL_SUBCLUSTER_UMAP_T_cells.png"),
       p_final, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

# --- Comparison: manual vs automated pre-scores (on clean embedding) ---
p_compare <- DimPlot(data_sub, reduction = "umap_harmony",
                     group.by = c("sub_cell_types",
                                  "sub_weighted_std_clean",
                                  "sub_weighted_raw_clean"),
                     label = FALSE, repel = TRUE)
ggsave(file.path(TCELL_DIR, "FINAL_SUBCLUSTER_comparison_T_cells.png"),
       p_compare, width = 26, height = 8, dpi = DPI_SETTING, bg = "white")


# =============================================================================
# 1. DEFINE LEVELS AND MAPPING OF DOTPLOT
# =============================================================================
# sub_type_levels2 is already defined and applied above in Part 6c.
name_mapping <- c(
  "CD8+ T cells"      = "CD8+ T cells",
  "γδ T cells"       = "γδ T cells",
  "Cyc. T cells" = "Cyc. T cells",
  "CD4+ T cells"      = "CD4+ T cells",
  "Tregs"             = "Tregs",
  "NK cells"          = "NK cells",
  "ILC2"              = "ILC2"
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
lead_genes <- c("Cd3e", "Cd3g", "Cd3d")

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

p_final_dot
ggsave(file.path(TCELL_DIR, "FINAL_DotPlot_sub_T_cells.png"),
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
  ggsave(file.path(TCELL_DIR, paste0("FINAL_sub_UMAP_faceted_T_cells_by_", pgrp, ".png")),
         p_facet, width = 5 * ceiling(n_levels / 2), height = 10, dpi = DPI_SETTING, bg = "white")
}

# =============================================================================
# --- PART 7: SAVE ------------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving T Cell Sub-Cluster Object ===")
saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds")))
message(paste0(
  "\n=== T CELL SUB-ANNOTATION COMPLETE ===\n",
  "  Saved: ", PROJECT_NAME, "_T_cells_subclustered.rds\n",
  "  Sub-types found: ", paste(levels(data_sub$sub_cell_types), collapse = ", "), "\n",
  "\nNext steps:\n",
  "  - Script 04: Colonocyte sub-annotation (04_colonocyte_subannotation.R)\n",
  "  - Script 05: DE + Two-Way ANOVA (05_DE_and_two_way_ANOVA.R)\n"
))


# =============================================================================
# --- PART 8: COMPOSITIONAL ANALYSIS ------------------------------------------
# =============================================================================
message("\n=== STEP 8: Compositional Analysis ===")
data_sub <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_T_cells_subclustered.rds")))

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
savef <- paste0(TCELL_DIR, "/tcells_stats_cell_props.png")
ggsave(savef,p_stats, width = 10, height = 10, dpi=DPI_SETTING)
savef <- paste0(TCELL_DIR, "/tcells_cell_props_long_data.xlsx")
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
  save_path     = file.path(TCELL_DIR, "Nr4a1_stats_plot_tcells.png"),
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
  save_path     = file.path(TCELL_DIR, "Nr4a1_custom_stats_plot_tcells.png"),
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
  save_path     = file.path(TCELL_DIR, "Nr4a2_stats_plot_tcells.png"),
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
  save_path     = file.path(TCELL_DIR, "Nr4a3_stats_plot_tcells.png"),
  dpi           = DPI_SETTING
)


plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Pdcd1",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(TCELL_DIR, "Pdcd1_stats_plot_tcells.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Ctla4",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(TCELL_DIR, "Ctla4_stats_plot_tcells.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Lag3",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(TCELL_DIR, "Lag3_stats_plot_tcells.png"),
  dpi           = DPI_SETTING
)

plot_results <- plot_expression_custom(
  seurat_obj    = data_sub,
  gene          = "Havcr2",
  plot_type     = "bar", 
  group_by      = "sub_cell_types",
  condition_col = "Genotype_sex",
  hide_x_text   = TRUE,      # Removes the crowded labels from the bottom of every plot
  show_legend   = TRUE,      # Puts a clean color key at the very bottom
  comparisons   = COMPARISON_PAIRS, # Use your 6 pairs
  facet_ncol    = 4,               # Arrange in 4 columns
  p_width     = 12,  # Force 12 inches wide
  p_height    = 8,  # Force 10 inches tall for multiple rows,
  save_path     = file.path(TCELL_DIR, "Havcr2_stats_plot_tcells.png"),
  dpi           = DPI_SETTING
)


library(Seurat)
library(AUCell)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

message("\n##### Starting Multi-T Cell Exhaustion Analysis #####")

# =============================================================================
# ENVIRONMENT & OUTPUT PATHS SETUP
# =============================================================================
#PROJECT_NAME <- "Nr4a1_s17_ack"
#ROOT_PATH    <- "/home/ssromerogon/2026_nr4a1_ack/r_process"
#OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
# (TCELL_DIR is already created at load time above)

# =============================================================================
# 1. SUBSETTING
# =============================================================================
t_cells_of_interest <- c("CD8+ T cells", "CD4+ T cells", "Tregs", "γδ T cells", "Cyc. T cells")
data_t <- subset(data_sub, subset = sub_cell_types %in% t_cells_of_interest)

# --- Final UMAP: clean manual annotation ---
DimPlot(data_t, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) 

# =============================================================================
# 2. GENE SIGNATURE PANELS
# =============================================================================
cd8_gene_sets <- list(
  Exhaustion_Score = c("Lag3", "Pdcd1", "Havcr2", "Tigit", "Ctla4", "Tox"),
  Effector_Score = c("Gzmb", "Gzma", "Gzmk", "Ifng", "Tnf", "Il2", "Tbx21", "Prf1")
)

# =============================================================================
# 3. AUCell ENRICHMENT PIPELINE (10% Threshold)
# =============================================================================
# Using universal matrix extraction
expr_matrix <- GetAssayData(object = data_t, assay = "RNA", layer = "counts")

cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)

# Filter out features missing from the expression matrix
valid_gene_sets <- list()
for (set_name in names(cd8_gene_sets)) {
  genes_in_set <- cd8_gene_sets[[set_name]]
  present_genes <- genes_in_set[genes_in_set %in% rownames(expr_matrix)]
  
  if (length(present_genes) > 0) {
    valid_gene_sets[[set_name]] <- present_genes
  }
}

# Run execution (Cleaned up the duplicate chunk)
cells_auc <- AUCell_calcAUC(valid_gene_sets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.10)
auc_matrix <- getAUC(cells_auc)

for (score_name in rownames(auc_matrix)) {
  data_t[[score_name]] <- as.numeric(auc_matrix[score_name, ])
}

# =============================================================================
# 4. EXPLICIT PLOTTING LOOP (PNG / 300 DPI / Custom Font Sizes)
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

scores_to_plot <- c("Exhaustion_Score", "Effector_Score")

for (score_col in scores_to_plot) {
  message("Generating plot for: ", score_col)
  
  plot_df <- data_t@meta.data %>%
    select(all_of(c(score_col, "Genotype_sex", "sub_cell_types"))) %>%
    drop_na()
  
  plot_df$Genotype_sex <- factor(plot_df$Genotype_sex, levels = custom_levels)
  
  clean_title <- gsub("_", " ", score_col)
  display_title <- paste("T cell", clean_title, "Across Subpopulations")
  
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
    
    facet_wrap(~ sub_cell_types, scales = "free_y", ncol = 3) +
    
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
  
  output_filename <- paste0("Multi_T_", score_col, "_Explicit_Groups.png")
  
  ggsave(
    filename = file.path(TCELL_DIR, output_filename), 
    plot = violin_p, 
    device = "png",
    dpi = 300,
    height = 10,
    width = 12, 
    units = "in"
  )
  
  print(violin_p)
}
