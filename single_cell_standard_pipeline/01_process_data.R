# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING ENGINE (v7 - HIGH-RES OUTPUT)
#
# PURPOSE:
# This script performs all heavy computational work and generates all necessary
# diagnostic plots for quality assurance. It outputs a single, fully processed
# Seurat object ready for interactive annotation in Script 2.
#
# =============================================================================
# --- PART 1: USER CONFIGURATION ---
library(Seurat); library(harmony); library(openxlsx); library(dplyr)
library(ggplot2); library(patchwork); library(celda); library(DoubletFinder)
library(writexl); library(Matrix)
set.seed(123)

PROJECT_NAME <- "Wu_IEC_Project"
ROOT_PATH <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
METADATA_FILE <- file.path(ROOT_PATH, "Wu_metadata.xlsx")
H5_DIR <- file.path(ROOT_PATH, "h5_files")
OUTPUT_DIR <- file.path(ROOT_PATH, "seurat_output")

# --- QC Parameters ---
PRE_MIN_GENES_PER_CELL <- 500; PRE_MAX_MT_PERCENT <- 30.0
POST_MIN_GENES <- 500; POST_MAX_GENES <- 14000; POST_MIN_UMIS <- 1500
POST_MAX_UMIS <- 100000; POST_MAX_MT <- 20.0; POST_MIN_CELLS_PER_GENE <- 15

# --- DoubletFinder Parameters ---
DOUBLET_RATE <- 0.08; DF_PK_RANGE_MIN <- 0.01
DF_PK_RANGE_MAX <- 0.15; DF_PK_FALLBACK <- 0.09

# --- Core Processing Parameters ---
N_VARIABLE_FEATURES <- 5000; N_PCS_TO_USE <- 50; CLUSTER_RESOLUTION <- 1.0
UMAP_N_NEIGHBORS <- 15; UMAP_MIN_DIST <- 0.3

# --- Workflow & Output Parameters ---
RUN_DECONTX <- TRUE
DPI_SETTING <- 300 # DPI for all saved PNG images

# =============================================================================
# --- PART 2: PIPELINE EXECUTION (DO NOT EDIT) ---
# =============================================================================
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR, recursive = TRUE) }

# --- 2.1: Load and Merge Data ---
message("--- Step 2.1: Loading and Merging Data ---")
metadata <- read.xlsx(METADATA_FILE)
seurat_objects_list <- list()
for (i in 1:nrow(metadata)) {
  sample_info <- metadata[i, ]; sample_id <- as.character(sample_info$SampleID)
  h5_path <- file.path(H5_DIR, sample_id, "sample_filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) { warning(paste("H5 file not found:", sample_id)); next }
  counts_matrix <- Read10X_h5(h5_path)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = sample_id, min.cells = 5)
  for (col_name in colnames(sample_info)) { seurat_obj[[col_name]] <- sample_info[[col_name]] }
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= PRE_MIN_GENES_PER_CELL & percent.mt <= PRE_MAX_MT_PERCENT)
  seurat_objects_list[[sample_id]] <- seurat_obj
}
if (length(seurat_objects_list) == 0) stop("No data loaded.")
data <- merge(x = seurat_objects_list[[1]], y = seurat_objects_list[-1])
data <- JoinLayers(data)

# --- 2.2: Post-Merge QC and Diagnostic Plotting ---
message("\n--- Step 2.2: Post-Merge QC and Plotting ---")
p_before <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "01a_qc_violin_before_filtering.png"), plot = p_before, width = 10, height = 6, dpi = DPI_SETTING)

data <- subset(data, subset = nFeature_RNA >= POST_MIN_GENES & nFeature_RNA <= POST_MAX_GENES &
                 nCount_RNA >= POST_MIN_UMIS & nCount_RNA <= POST_MAX_UMIS & percent.mt <= POST_MAX_MT)
genes_to_keep <- rownames(data)[Matrix::rowSums(GetAssayData(data, layer = "counts") > 0) >= POST_MIN_CELLS_PER_GENE]
data <- subset(data, features = genes_to_keep)

p_after <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "01b_qc_violin_after_filtering.png"), plot = p_after, width = 10, height = 6, dpi = DPI_SETTING)
message(paste("Cells after QC:", ncol(data), "| Genes remaining:", nrow(data)));

# --- 2.3: Doublet Detection with Diagnostic Plotting ---
message("\n--- Step 2.3: Doublet Detection and Plotting ---")
pk_plot_dir <- file.path(OUTPUT_DIR, "doublet_finder_plots"); if (!dir.exists(pk_plot_dir)) dir.create(pk_plot_dir)
data_list <- SplitObject(data, split.by = "SampleID")
results_list <- list()
for (i in seq_along(data_list)) {
  sample_name <- names(data_list)[i]; message(paste("... finding doublets in", sample_name))
  seu_tmp <- data_list[[i]]
  seu_tmp <- NormalizeData(seu_tmp, verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(npcs=N_PCS_TO_USE, verbose=F)
  sweep.res.list <- paramSweep(seu_tmp, PCs = 1:N_PCS_TO_USE, sct = F); sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
  bcmvn <- find.pK(sweep.stats); bcmvn$pK <- as.numeric(as.character(bcmvn$pK));
  
  initial_optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  final_pk <- initial_optimal_pk
  if (is.na(final_pk) || final_pk < DF_PK_RANGE_MIN || final_pk > DF_PK_RANGE_MAX) {
    bcmvn_filtered <- bcmvn[bcmvn$pK > DF_PK_RANGE_MIN & bcmvn$pK < DF_PK_RANGE_MAX, ]
    if (nrow(bcmvn_filtered) > 0 && any(is.finite(bcmvn_filtered$BCmetric))) { final_pk <- bcmvn_filtered$pK[which.max(bcmvn_filtered$BCmetric)]
    } else { final_pk <- DF_PK_FALLBACK }
  }
  
  pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_line() + geom_point() + geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
    ggtitle(paste0("pK Finder for ", sample_name), subtitle = paste0("Final Selected pK = ", final_pk)) + theme_minimal()
  ggsave(file.path(pk_plot_dir, paste0(sample_name, "_pk_plot.png")), plot = pk_plot, width = 7, height = 5, dpi = DPI_SETTING)
  
  nExp_val <- round(ncol(seu_tmp) * DOUBLET_RATE)
  seu_tmp <- doubletFinder(seu_tmp, PCs = 1:N_PCS_TO_USE, pK = final_pk, nExp = nExp_val, sct = FALSE)
  res_col <- tail(grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE), 1)
  results_list[[sample_name]] <- seu_tmp@meta.data[, res_col, drop = FALSE]
}
all_res <- do.call(rbind, lapply(results_list, function(df) setNames(df, "Doublet_Status")))
data$Doublet_Status <- all_res[rownames(data@meta.data), "Doublet_Status"]

# Visualize doublets on a temporary UMAP (Memory Efficient)
message("...generating diagnostic plot for doublet visualization.")
data <- NormalizeData(data, verbose = F) %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(verbose = F) %>%
        RunPCA(npcs = N_PCS_TO_USE, verbose = F, reduction.name = "pca_temp_doublets") %>%
        RunUMAP(dims = 1:N_PCS_TO_USE, reduction = "pca_temp_doublets", reduction.name = "umap_temp_doublets", verbose = F)

p_doublets <- DimPlot(data, reduction = "umap_temp_doublets", group.by = "Doublet_Status",
                      cols = c("Singlet" = "grey", "Doublet" = "red")) +
              ggtitle("Doublet Visualization Before Removal")
ggsave(file.path(OUTPUT_DIR, "01c_DIAGNOSTIC_doublet_visualization.png"), plot = p_doublets, width = 8, height = 7, dpi = DPI_SETTING)

data@reductions$pca_temp_doublets <- NULL
data@reductions$umap_temp_doublets <- NULL
rm(p_doublets); gc()

# Filter out the identified doublets
message(paste("Total cells before doublet filtering:", ncol(data)))
data <- subset(data, subset = Doublet_Status != "Doublet" | is.na(Doublet_Status))
message(paste("Total cells after doublet filtering:", ncol(data)))

# --- 2.4: Ambient RNA Correction (Optional) ---
if (RUN_DECONTX) {
    message("\n--- Step 2.4: DecontX ---")
    counts_sparse <- GetAssayData(object = data, layer = "counts")
    decontx_results <- decontX(x = counts_sparse)
    data[["RNA"]]$counts <- decontx_results$decontXcounts; data[["RNA"]]$data <- NULL
    data$nCount_RNA <- colSums(data[["RNA"]]$counts); data$nFeature_RNA <- colSums(data[["RNA"]]$counts > 0)
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^mt-")
    data <- subset(data, subset = nFeature_RNA >= POST_MIN_GENES & nFeature_RNA <= POST_MAX_GENES &
                     nCount_RNA >= POST_MIN_UMIS & nCount_RNA <= POST_MAX_UMIS & percent.mt <= POST_MAX_MT)
}

# --- 2.5: Processing, Integration, and Clustering with Diagnostic Plotting ---
message("\n--- Step 2.5: Processing, Integration, and Clustering ---")
data <- NormalizeData(data, verbose=F) %>%
        FindVariableFeatures(nfeatures=N_VARIABLE_FEATURES, verbose=F) %>%
        ScaleData(verbose=F) %>%
        RunPCA(npcs=N_PCS_TO_USE, verbose=F)

# VERSION 1: Non-Harmony (Standard PCA)
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "pca", graph.name = "pca_nn", verbose = FALSE)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, graph.name = "pca_nn", cluster.name = "clusters_none", verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "pca", n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST, reduction.name = "umap_none", verbose = FALSE)

# VERSION 2: Harmony (Batch Corrected)
data <- RunHarmony(data, group.by.vars = "SampleID", reduction = "pca", reduction.save = "harmony", verbose = FALSE)
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "harmony", graph.name = "harmony_nn", verbose = FALSE)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, graph.name = "harmony_nn", cluster.name = "clusters_harmony", verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "harmony", n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST, reduction.name = "umap_harmony", verbose = FALSE)

# Generate Diagnostic Comparison Plots
p1 <- DimPlot(data, reduction = "umap_none", group.by = "SampleID") + ggtitle("Standard PCA (No Harmony)")
p2 <- DimPlot(data, reduction = "umap_harmony", group.by = "SampleID") + ggtitle("Harmony Corrected")
ggsave(file.path(OUTPUT_DIR, "02_DIAGNOSTIC_UMAP_Harmony_vs_NoHarmony.png"), plot = p1 + p2, width = 16, height = 7, dpi = DPI_SETTING)

p3 <- DimPlot(data, reduction = "umap_none", group.by = "clusters_none", label = TRUE) + ggtitle("Clusters: Standard PCA") + NoLegend()
p4 <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE) + ggtitle("Clusters: Harmony") + NoLegend()
ggsave(file.path(OUTPUT_DIR, "03_DIAGNOSTIC_UMAP_Cluster_Comparison.png"), plot = p3 + p4, width = 16, height = 7, dpi = DPI_SETTING)

# --- 2.6: Saving Processed Object for Annotation ---
message("\n--- Step 2.6: Saving Processed Object ---")
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds")))
message(paste0("\n\n--- PROCESSING COMPLETE ---",
               "\n\nNext Step: Open and run '02_annotate_data.R' to explore and annotate the data."))

