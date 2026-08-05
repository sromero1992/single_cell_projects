# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 5: MACROPHAGE SUB-ANNOTATION
# Version: 1.0 (CSV-Driven, Seurat Wrappers, Harmony + Standard Clustering)
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
#
# PURPOSE:
#   Loads the globally annotated Seurat object from Script 02 and performs
#   high-resolution sub-clustering of the Macrophage/Myeloid compartment. Provides:
#     1. Subset extraction with fresh HVG selection, PCA, Harmony, and UMAP.
#        BOTH a Harmony-corrected and a standard (non-Harmony) embedding are
#        produced for side-by-side comparison.
#     2. Weighted pre-scoring against macrophage sub-type markers (CSV-driven).
#        Both z-scored (standardized) and raw (non-standardized) variants run
#        automatically — review the Top-5 report before manual annotation.
#     3. Manual sub-annotation via SUB_ANNOTATION_MAP (Action 5).
#     4. Doublet removal → re-UMAP → re-clustering → re-annotation pass.
#     5. Compositional analysis (proportions by SampleID and group).
#     6. Gene expression violin/bar plots + AUCell functional scoring.
#
# MARKER SYSTEM (cell_type_markers.csv):
#   This script reads rows where tier == "sub" AND parent_cell_type == "Macrophages".
#   If no sub-rows exist the script will warn and fall back to DEFAULT_SUB_MARKERS.
#
# HOW TO USE:
#   1. Set paths/parameters in Part 1.
#   2. Run through Part 4 — review SUBCLUSTER_01 and SUBCLUSTER_02 PNGs.
#   3. Fill in SUB_ANNOTATION_MAP in Action 5.
#   4. Run Part 6a QC block — review CLEANING_ PNGs and marker xlsx.
#   5. Fill in SUB_ANNOTATION_MAP_CLEAN in Action 6.
#   6. Run the rest for final plots and save.
#
# INPUT:  {PROJECT_NAME}_final_annotated.rds   (output of 02_global_annotation.R)
# OUTPUT: {PROJECT_NAME}_Macrophages_subclustered.rds
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
PROJECT_NAME <- "Wu_Diet_project2"
#ROOT_PATH <- "/home/ssromerogon/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"  # Windows
ROOT_PATH <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"

OUTPUT_DIR       <- file.path(ROOT_PATH, "seurat_output")
MAC_DIR          <- file.path(OUTPUT_DIR, "macrophages_subannotation")  # <-- all macrophage plots go here
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
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(MAC_DIR))    dir.create(MAC_DIR,    recursive = TRUE)
message(paste("  Macrophage output folder:", MAC_DIR))
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
  message("  Running UMAP on PCA...")
  data <- RunUMAP(data, dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh,
                  min.dist = min_dist, n.epochs = 500,
                  reduction.name = "umap_none", verbose = FALSE)
  gc()
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
ggsave(file.path(MAC_DIR, "SUBCLUSTER_01_Harmony_comparison_Macrophages.png"),
       p_harm_comp, width = 16, height = 8, dpi = DPI_SETTING)

# --- 4.2: Cluster UMAPs (annotate from these) ---
p_clust_none <- DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_none",
                        label = TRUE) + NoLegend() + ggtitle("Clusters: Standard PCA")
p_clust_harm <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony",
                        label = TRUE) + NoLegend() + ggtitle("Clusters: Harmony")
p_clust_comp <- p_clust_none + p_clust_harm
ggsave(file.path(MAC_DIR, "SUBCLUSTER_02_COMPARE_UMAP_Macrophages.png"),
       p_clust_comp, width = 16, height = 8, dpi = DPI_SETTING)
ggsave(file.path(MAC_DIR, "SUBCLUSTER_02_ANNOTATE_THIS_UMAP_Macrophages.png"),
       p_clust_harm, width = 10, height = 8, dpi = DPI_SETTING)

# --- 4.3: DotPlot for cluster characterization ---
p_dot_sub <- DotPlot(data_sub, features = sub_dotplot_markers,
                     group.by = "clusters_harmony", dot.min = 0.05, cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11),  cluster.idents = T) +
  ggtitle(paste("Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(MAC_DIR, "SUBCLUSTER_02_ANNOTATE_USING_THIS_DOTPLOT_Macrophages.png"),
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
           file.path(MAC_DIR, "SUBCLUSTER_PRESCORE_top5_Macrophages_std.xlsx"))
write_xlsx(sub_results_raw$top5_report,
           file.path(MAC_DIR, "SUBCLUSTER_PRESCORE_top5_Macrophages_raw.xlsx"))

p_prescore <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw",
          label = TRUE, repel = TRUE) + ggtitle("Sub Pre-Score: Raw")
)
ggsave(file.path(MAC_DIR, "SUBCLUSTER_03_prescore_comparison_Macrophages.png"),
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

# 'cDCs', # Itgax+ Clec9a+ → cDCs (no Cd68)
# 'TAM resident',#  Folr2+ Timd4+  C1q1/b+ and Retnla- Arg1-
# 'Monocytes', # Cd14+ Ccr2+ Ly6c2+ → Monocytes (no Timd4 no Adgre1)
# 'pDCs', # Siglech+ Bst2+  → pDCs (ignore everything else)
SUB_ANNOTATION_MAP <- c(
  '0'  = 'M2 Macrophages',#
  '1'  = 'M0 Macrophages',#
  '2'  = 'M2 Macrophages',#
  '3'  = 'M0 Macrophages',#
  '4'  = 'M0 Macrophages',#
  '5'  = 'M0 Macrophages',#
  '6'  = 'cDC1',#
  '7'  = 'cDC1',#
  '8'  = 'M2 Macrophages',#
  '9'  = 'Cyc. cDC1',#
  '10' = 'M0 Macrophages',#
  '11' = 'cDC2',#
  '12' = 'Monocytes',#
  '13' = 'cDC2',#
  '14' = 'M1 Macrophages',#
  '15' = 'mRegDCs',#
  '16' = 'M0 Macrophages',#
  '17' = 'Monocytes',#
  '18' = 'M2 Macrophages',# 
  '19' = 'M1 Macrophages',#
  '20' = 'Cyc. Macrophages',#
  '21' = 'Doublet', #
  '22' = 'Monocytes',#
  '23' = 'Neutrophils',#
  '24' = 'Doublet',#
  '25' = 'M2 Macrophages',#
  '26' = 'M2 Macrophages',#
  '27' = 'Doublet',#
  '28' = 'M1 Macrophages',#
  '29' = 'Doublet',#
  '30' = 'pDCs', #
  '31' = 'Doublet'#
)

# =============================================================================
# --- PART 6: FIRST-PASS ANNOTATION & DIRTY-CLUSTER QC -----------------------
# =============================================================================
message("\n=== STEP 6a: First-Pass Sub-Annotation (pre-cleaning) ===")

# Apply ACTION 5 annotation map.
# Label anything uncertain as "Doublet" or "Unknown" — removed in 6b.
data_sub$sub_cell_types  <- recode_factor(data_sub$clusters_harmony, !!!SUB_ANNOTATION_MAP)
data_sub$CellType        <- data_sub$sub_cell_types
data_sub$seurat_clusters <- data_sub$clusters_harmony

# --- 6a.1: Pre-cleaning overview UMAP ----------------------------------------
p_dirty <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) +
  ggtitle("Pre-Cleaning: All Clusters Labelled (including Doublet/Unknown)") +
  theme(legend.text = element_text(size = 11))
ggsave(file.path(MAC_DIR, "CLEANING_00_pre_removal_overview.png"),
       p_dirty, width = 12, height = 8, dpi = DPI_SETTING, bg = "white")

# --- 6a.2: QC FeaturePlots ---------------------------------------------------
message("  Generating QC FeaturePlots for dirty cluster inspection...")

# Broad myeloid identity markers
p_qc1 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Cd68", "Adgre1", "Csf1r", "Itgam",
                                            "Mki67", "Fcgr3", "nCount_RNA", "nFeature_RNA"),
                                          rownames(data_sub)),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(MAC_DIR, "CLEANING_01_QC_broad_myeloid_markers.png"),
       p_qc1, width = 16, height = 8, dpi = DPI_SETTING, bg = "white")

# Contamination / non-myeloid markers
# Non-myeloid contamination markers (defined here for use in p_qc2)
contaminant_markers <- c(
  "Col1a1",  # Fibroblasts
  "Myh11",   # Smooth muscle cells
  "Cd3e",    # T cells
  "Cpa3",    # Mast cells
  "Mzb1",    # Plasma B cells
  "Slc1a2",  # Neurons/enteric glia
  "Cldn5",   # Vascular endothelial cells
  "Cd8b1",   # T cells / lymphoid
  "Cd19"     # B cells
)
p_qc2 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(contaminant_markers, rownames(data_sub)),
                     ncol = 3) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(MAC_DIR, "CLEANING_02_QC_contamination_markers.png"),
       p_qc2, width = 14, height = 10, dpi = DPI_SETTING, bg = "white")

# DC subtype disambiguation
p_qc3 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Itgax", "Clec9a", "Xcr1", "Cd209a",
                                            "Siglech", "Bst2", "Ccr7", "Ccl17"),
                                          rownames(data_sub)),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(MAC_DIR, "CLEANING_03_DC_disambiguation.png"),
       p_qc3, width = 16, height = 8, dpi = DPI_SETTING, bg = "white")

# TAM / M1 / M2 state markers
p_qc4 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Spp1", "Fn1", "Gpnmb", "Arg1", "Mrc1",
                                            "Nos2", "Ifit1", "Isg15", "Folr2", "Timd4"),
                                          rownames(data_sub)),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(MAC_DIR, "CLEANING_04_TAM_M1_M2_markers.png"),
       p_qc4, width = 16, height = 10, dpi = DPI_SETTING, bg = "white")

# Monocyte / Neutrophil markers
p_qc5 <- FeaturePlot(data_sub, reduction = "umap_harmony",
                     features = intersect(c("Ly6c2", "Ccr2", "Cd14", "S100a8",
                                            "S100a9", "Csf3r", "Cxcr2"),
                                          rownames(data_sub)),
                     ncol = 4) &
  theme(plot.title = element_text(size = 13, face = "italic"), axis.title = element_blank())
ggsave(file.path(MAC_DIR, "CLEANING_05_monocyte_neutrophil_markers.png"),
       p_qc5, width = 14, height = 8, dpi = DPI_SETTING, bg = "white")

# --- 6a.3: FindAllMarkers for all first-pass clusters -------------------------
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
  file.path(MAC_DIR, "CLEANING_top30_per_cluster.xlsx")
)
message(paste("  Marker table saved —",
              length(unique(all_cluster_markers$cluster)), "clusters,",
              nrow(top30_by_cluster), "rows total."))

# =============================================================================
# --- PART 6b: REMOVE DIRTY CLUSTERS & RE-EMBED ------------------------------
# =============================================================================
message("\n=== STEP 6b: Removing dirty clusters and re-running UMAP ===")

DIRTY_LABELS <- c("Doublet", "Unknown")   # <-- ACTION: edit if you used different labels

n_before <- ncol(data_sub)
data_sub  <- subset(data_sub, subset = sub_cell_types %nin% DIRTY_LABELS)
n_after   <- ncol(data_sub)
message(paste("  Removed", n_before - n_after, "dirty cells.",
              n_after, "clean cells remain."))

# Drop reductions so run_pca_umap starts fresh
data_sub@reductions <- list()
data_sub@graphs     <- list()

message("  Re-running HVG / PCA / Harmony / UMAP on clean cells...")
SUBCLUSTER_K_NEIGHBORS = 15
SUBCLUSTER_N_PCS = 20
data_sub <- run_pca_umap(
  data     = data_sub,
  num_hvg  = SUBCLUSTER_N_HVG,
  dims_pca = SUBCLUSTER_N_PCS,
  min_dist = SUBCLUSTER_MIN_DIST,
  kneigh   = SUBCLUSTER_K_NEIGHBORS
)
gc()

# Re-cluster on clean embedding
message("  Re-clustering clean cells on harmony embedding...")
SUBCLUSTER_RESOLUTION_FINAL = 1.5
data_sub <- FindNeighbors(data_sub, dims = 1:SUBCLUSTER_N_PCS, reduction = "harmony",
                          k.param = SUBCLUSTER_K_NEIGHBORS, graph.name = "harmony_nn",
                          verbose = FALSE) %>%
  FindClusters(resolution = SUBCLUSTER_RESOLUTION_FINAL, graph.name = "harmony_nn",
               cluster.name = "clusters_harmony_clean", verbose = FALSE)
Idents(data_sub) <- "clusters_harmony_clean"
gc()

# --- 6b.1: Save clean cluster UMAPs for re-annotation review ----------------
p_clean_clust <- DimPlot(data_sub, reduction = "umap_harmony",
                         group.by = "clusters_harmony_clean",
                         label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("CLEAN Clusters: Harmony (re-embedded)")
ggsave(file.path(MAC_DIR, "CLEANING_06_clean_clusters_ANNOTATE_THIS.png"),
       p_clean_clust, width = 10, height = 8, dpi = DPI_SETTING, bg = "white")

p_clean_both <- (
  DimPlot(data_sub, reduction = "umap_none",    group.by = "clusters_harmony_clean",
          label = TRUE) + NoLegend() + ggtitle("Clean — Standard PCA")
) + (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "clusters_harmony_clean",
          label = TRUE) + NoLegend() + ggtitle("Clean — Harmony")
)
ggsave(file.path(MAC_DIR, "CLEANING_07_clean_compare_both_embeddings.png"),
       p_clean_both, width = 16, height = 8, dpi = DPI_SETTING)

# --- 6b.2: Re-run weighted pre-scoring on clean clusters --------------------
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
           file.path(MAC_DIR, "CLEANING_08_prescore_clean_std.xlsx"))
write_xlsx(sub_results_clean_raw$top5_report,
           file.path(MAC_DIR, "CLEANING_08_prescore_clean_raw.xlsx"))

p_prescore_clean <- (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_std_clean",
          label = TRUE, repel = TRUE) + ggtitle("Clean Pre-Score: Standardized")
) | (
  DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_weighted_raw_clean",
          label = TRUE, repel = TRUE) + ggtitle("Clean Pre-Score: Raw")
)
ggsave(file.path(MAC_DIR, "CLEANING_09_prescore_clean_comparison.png"),
       p_prescore_clean, width = 18, height = 8, dpi = DPI_SETTING, bg = "white")

# --- 6b.3: DotPlot for clean cluster characterization -----------------------
p_dot_clean <- DotPlot(data_sub, features = sub_dotplot_markers,
                       group.by = "clusters_harmony_clean", dot.min = 0.05,
                       cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11)) +
  ggtitle(paste("Clean Sub-Cluster Markers:", PARENT_CELL_TYPE))
ggsave(file.path(MAC_DIR, "CLEANING_10_clean_ANNOTATE_USING_THIS_DOTPLOT.png"),
       p_dot_clean, width = 14, height = 9, dpi = DPI_SETTING, bg = "white")

# --- 6b.4: FindAllMarkers on clean clusters for final verification ----------
message("  Running FindAllMarkers on clean clusters...")
all_cluster_markers_clean <- FindAllMarkers(
  data_sub,
  only.pos        = TRUE,
  min.pct         = 0.05,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

top30_clean <- all_cluster_markers_clean %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 30) %>%
  dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::ungroup()

write_xlsx(
  as.data.frame(top30_clean),
  file.path(MAC_DIR, "CLEANING_top30_per_cluster_clean.xlsx")
)
message(paste("  Clean marker table saved —",
              length(unique(all_cluster_markers_clean$cluster)), "clusters,",
              nrow(top30_clean), "rows total."))

# =============================================================================
# --- ACTION 6: FILL IN CLEAN CLUSTER ANNOTATIONS HERE -----------------------
# =============================================================================
# Review (in MAC_DIR):
#   CLEANING_06_clean_clusters_ANNOTATE_THIS.png      — new cluster layout
#   CLEANING_10_clean_ANNOTATE_USING_THIS_DOTPLOT.png — marker expression
#   CLEANING_08_prescore_clean_std.xlsx               — top-5 per clean cluster
#   CLEANING_09_prescore_clean_comparison.png         — automated suggestions
#   CLEANING_top30_per_cluster_clean.xlsx             — full top 30 DE genes
#
# Map each NEW cluster number (as a string) to a macrophage sub-type name.
# These will differ from SUB_ANNOTATION_MAP because cells were re-embedded
# and re-clustered after doublet removal.

  SUB_ANNOTATION_MAP_CLEAN <- c(
  '0'  = 'M2 Macrophages',#      
  '1'  = 'M0 Macrophages',#
  '2'  = 'M2 Macrophages',#
  '3'  = 'cDC1',#
  '4'  = 'M0 Macrophages',#
  '5'  = 'M2 Macrophages',#
  '6'  = 'cDC2',#
  '7'  = 'M1 Macrophages',#
  '8'  = 'M1 Macrophages',#
  '9'  = 'Monocytes',#
  '10' = 'M2 Macrophages',#
  '11' = 'cDC2',#
  '12' = 'Cyc. cDC1',#
  '13' = 'cDC1',#
  '14' = 'mRegDCs',#
  '15' = 'M1 Macrophages',#
  '16' = 'M2 Macrophages',#
  '17' = 'Monocytes',#
  '18' = 'Cyc. Macrophages',#
  '19' = 'M2 Macrophages',#
  '20' = 'Monocytes',#
  '21' = 'Cyc. cDC2',#
  '22' = 'Neutrophils',#
  '23' = 'M2 Macrophages',#
  '24' = 'cDC2',#
  '25' = 'M1 Macrophages',#
  '26' = 'pDCs' #
)

# =============================================================================
# --- PART 6c: APPLY CLEAN ANNOTATIONS & FINAL VISUALIZATIONS ----------------
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

sub_type_levels_macrophages <- c(
  "pDCs",
  "cDC1", "Cyc. cDC1", "cDC2", "Cyc. cDC2", "mRegDCs",
  "Monocytes",
  "M0 Macrophages", "M1 Macrophages", "M2 Macrophages", "Cyc. Macrophages",
  "Neutrophils"
)
data_sub$sub_cell_types <- factor(data_sub$sub_cell_types, levels = sub_type_levels_macrophages)

# Drop reductions and re-embed without the removed cells
data_sub@reductions <- list()
data_sub@graphs     <- list()

message("  Re-running UMAP after final clean subset...")
data_sub <- run_pca_umap(
  data     = data_sub,
  num_hvg  = SUBCLUSTER_N_HVG,
  dims_pca = SUBCLUSTER_N_PCS,
  min_dist = SUBCLUSTER_MIN_DIST,
  kneigh   = SUBCLUSTER_K_NEIGHBORS
)
gc()

# --- Final UMAP: clean manual annotation ------------------------------------
p_final <- DimPlot(data_sub, reduction = "umap_harmony", group.by = "sub_cell_types",
                   label = TRUE, repel = TRUE) +
  ggtitle(paste("Sub-Cell Types:", PARENT_CELL_TYPE, "(Clean)")) +
  theme(legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(file.path(MAC_DIR, "FINAL_SUBCLUSTER_UMAP_Macrophages.png"),
       p_final, width = 12, height = 8, dpi = DPI_SETTING, bg = "white")

# --- Comparison: manual vs automated pre-scores (on clean embedding) --------
p_compare <- DimPlot(data_sub, reduction = "umap_harmony",
                     group.by = c("sub_cell_types",
                                  "sub_weighted_std_clean",
                                  "sub_weighted_raw_clean"),
                     label = FALSE, repel = TRUE)
ggsave(file.path(MAC_DIR, "FINAL_SUBCLUSTER_comparison_Macrophages.png"),
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
  "Cyc. cDC1"        = "cDC1",       # maps to cDC1 markers
  "cDC2"             = "cDC2",
  "Cyc. cDC2"        = "cDC2",       # maps to cDC2 markers
  "mRegDCs"          = "mRegDCs",
  "Monocytes"        = "Monocytes",
  "M0 Macrophages"   = "M0 macrophages",
  "M1 Macrophages"   = "M1 macrophages",
  "M2 Macrophages"   = "M2 macrophages",
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
ggsave(file.path(MAC_DIR, "FINAL_DotPlot_sub_macrophages.png"),
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
  ggsave(file.path(MAC_DIR, paste0("FINAL_sub_UMAP_faceted_Macrophages_by_", pgrp, ".png")),
         p_facet, width = 5 * ceiling(n_levels / 2), height = 10, dpi = DPI_SETTING, bg = "white")
}

# =============================================================================
# --- PART 7: SAVE ------------------------------------------------------------
# =============================================================================
message("\n=== STEP 7: Saving Macrophage Sub-Cluster Object ===")
saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds")))
message(paste0(
  "\n=== MACROPHAGE SUB-ANNOTATION COMPLETE ===\n",
  "  Saved: ", PROJECT_NAME, "_macrophages_subclustered.rds\n",
  "  Sub-types found: ", paste(levels(data_sub$sub_cell_types), collapse = ", "), "\n",
  "\nNext steps:\n",
  "  - Script 06: DE + Two-Way ANOVA (06_DE_and_two_way_ANOVA.R)\n"
))


data_sub <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds")))

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
savef <- paste0(MAC_DIR, "/macrophages_stats_cell_props.png")
ggsave(savef,p_stats, width = 12, height = 10, dpi=DPI_SETTING)
savef <- paste0(MAC_DIR, "/macrophages_cell_props_long_data.xlsx")
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
  save_path     = file.path(MAC_DIR, "Nr4a1_stats_plot_macrophages.png"),
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
  save_path     = file.path(MAC_DIR, "Nr4a1_custom_stats_plot_macrophages.png"),
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
  save_path     = file.path(MAC_DIR, "Nr4a2_stats_plot_macrophages.png"),
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
  save_path     = file.path(MAC_DIR, "Nr4a3_stats_plot_macrophages.png"),
  dpi           = DPI_SETTING
)



library(Seurat)
library(AUCell)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

message("\n##### Starting Targeted Macrophage Functional Scoring #####")
data_sub <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds")))

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
mac_lineages <- c(
  "Monocytes", "Cyc. Macrophages",
  "M0 Macrophages", "M1 Macrophages", "M2 Macrophages"
)

data_mac <- subset(data_sub, subset = sub_cell_types %in% mac_lineages)

# =============================================================================
# 2. DEFINE MACROPHAGE FUNCTIONAL GENE SETS
# =============================================================================
mac_gene_sets <- list(
  M1_ProInflammatory   = c("Stat1", "Irf1", "Nos2", "Irf5", "Il6",
                           "Cxcl10", "Cxcl9", "Tnf", "Il1b"),
  M2_Immunosuppressive = c("Arg1", "Chil3", "Retnla", "Mrc1", "Il10",
                           "Tgfb1", "Ccl17", "Ccl22", "Cd200r1", "Socs1",
                           "C1qa", "C1qb", "C1qc"),
  IFN_Response         = c("Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1",
                           "Oas1a", "Irf7", "Stat1", "Cxcl10", "Isg20"),
  # --- NEW: Estrogen-NR4a1 axis ---
  Estrogen_Signaling        = c("Esr1", "Esr2", "Gper1",
                                "Foxa1", "Gata3",
                                "Tff1", "Cav1", "Cited1"),
  NR4a1_Activity            = c("Nr4a1", "Nr4a2", "Nr4a3",
                                "Sik1", "Dusp1", "Dusp6",
                                "Thbd", "Ccl2", "Il6"),
  Estrogen_AntiInflammatory = c("Il10", "Tgfb1", "Foxp3",
                                "Socs1", "Socs3",
                                "Ncoa1", "Ncoa2", "Nfkbia")
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

scores_to_plot <- c("M1_ProInflammatory", "M2_Immunosuppressive",
                    "IFN_Response",
                    "Estrogen_Signaling", "NR4a1_Activity",
                    "Estrogen_AntiInflammatory")

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
  
  output_filename <- paste0("Macrophage_", score_col, "_Explicit_Groups.png")
  
  ggsave(
    filename = file.path(MAC_DIR, output_filename), 
    plot = violin_p, 
    device = "png",
    dpi = 300,
    height = 10,    # Increased height slightly to accommodate the 8 faceted clusters cleanly
    width = 14, 
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

scores_to_plot <- c("M1_ProInflammatory", "M2_Immunosuppressive", 
                    "IFN_Response", "TAM_SPP1", "TAM_Resident")
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
    filename = file.path(MAC_DIR, output_filename), 
    plot = global_p, 
    device = "png",
    dpi = 300,
    height = 7,    
    width = 8, 
    units = "in"
  )
  
  print(global_p)
}

