# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING
#
# PURPOSE:
# This script performs all the heavy computational work: loading data, QC,
# doublet removal, ambient RNA correction, normalization, and clustering.
# It outputs a single Seurat object ready for annotation, along with the
# essential plots needed for that task.
#
# INSTRUCTIONS:
# 1. Configure all parameters in the "USER CONFIGURATION" section.
# 2. Run the entire script.
# 3. After completion, proceed to '02_annotate_data.R' for cell type annotation.
#
# =============================================================================
# --- PART 1: USER CONFIGURATION ---
# --- A. Setup: Load Libraries ---
library(Seurat); library(harmony); library(openxlsx); library(dplyr)
library(ggplot2); library(patchwork); library(celda); library(DoubletFinder)
library(writexl); library(Matrix)
set.seed(123)

# --- B. File and Directory Paths ---
PROJECT_NAME <- "Wu_IEC_Project"
ROOT_PATH <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
METADATA_FILE <- file.path(ROOT_PATH, "Wu_metadata.xlsx")
H5_DIR <- file.path(ROOT_PATH, "h5_files")
OUTPUT_DIR <- file.path(ROOT_PATH, "seurat_output")

# --- C. Quality Control (QC) Parameters ---
PRE_MIN_GENES_PER_CELL <- 500; PRE_MAX_MT_PERCENT <- 30.0
POST_MIN_GENES <- 500; POST_MAX_GENES <- 14000; POST_MIN_UMIS <- 1500
POST_MAX_UMIS <- 100000; POST_MAX_MT <- 20.0; POST_MIN_CELLS_PER_GENE <- 15

# --- D. Doublet Detection Parameters ---
DOUBLET_RATE <- 0.08; DF_PK_RANGE_MIN <- 0.01
DF_PK_RANGE_MAX <- 0.15; DF_PK_FALLBACK <- 0.09

# --- E. Core Processing & Clustering Parameters ---
N_VARIABLE_FEATURES <- 5000; N_PCS_TO_USE <- 50; CLUSTER_RESOLUTION <- 1.0
UMAP_N_NEIGHBORS <- 15; UMAP_MIN_DIST <- 0.3; UMAP_N_EPOCHS <- 500

# --- F. Workflow Control ---
RUN_DECONTX <- TRUE
RUN_HARMONY_COMPARISON <- TRUE

# --- G. Marker Genes for Annotation Plots ---
# These markers are used to generate the dot plot for the biologist.
BROAD_MARKERS_LIST <- list(
  "Plasma B cells" = c("Ighg1", "Mzb1"), "B cells" = c("Ms4a1", "Cd19"),
  "T cells" = c("Cd3e", "Cd4", "Cd8a", "Trbc2"), "cDCs" = c("Clec9a", "Xcr1"),
  "Macrophages" = c("Cd163", "Csf1r", "Cd68", "Mrc1", "Aif1"),
  "Fibroblasts" = c("Col1a1", "Col3a1", "Col6a2"), "ENs" = c("Vip", "Tubb3", "Vat1l"),
  "EGCs" = c("Sox10", "Plp1"), "LECs" = c("Lyve1", "Flt4", "Pecam1"),
  "VECs" = c("Plvap", "Flt1"), "Colonocytes" = c("Epcam", "Muc2", "Krt19", "Krt20", "Vil1"),
  "SMCs" = c("Myh11", "Tagln", "Hhip", "Des"), "Adipocytes" = c("Adipoq", "Plin1", "Plin4")
)

# =============================================================================
# --- PART 2: PIPELINE EXECUTION (DO NOT EDIT) ---
# =============================================================================
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR, recursive = TRUE) }

# --- 2.1: Load, Merge, and Prepare Data ---
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
rm(seurat_objects_list); gc()
data <- JoinLayers(data)

# --- 2.2: Post-Merge QC ---
message("\n--- Step 2.2: Post-Merge QC ---")
p_before <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "01_qc_violin_before_filtering.png"), plot = p_before, width = 10, height = 6)
data <- subset(data, subset = nFeature_RNA >= POST_MIN_GENES & nFeature_RNA <= POST_MAX_GENES &
                 nCount_RNA >= POST_MIN_UMIS & nCount_RNA <= POST_MAX_UMIS & percent.mt <= POST_MAX_MT)
genes_to_keep <- rownames(data)[Matrix::rowSums(GetAssayData(data, layer = "counts") > 0) >= POST_MIN_CELLS_PER_GENE]
data <- subset(data, features = genes_to_keep)
message(paste("Cells after QC:", ncol(data), "| Genes remaining:", nrow(data))); gc()

# --- 2.3: Doublet Detection ---
message("\n--- Step 2.3: Doublet Detection ---")
pk_plot_dir <- file.path(OUTPUT_DIR, "doublet_finder_plots"); if (!dir.exists(pk_plot_dir)) dir.create(pk_plot_dir)
data_list <- SplitObject(data, split.by = "SampleID")
results_list <- list()
for (i in seq_along(data_list)) {
  sample_name <- names(data_list)[i]; message(paste("... finding doublets in", sample_name))
  seu_tmp <- data_list[[i]]
  seu_tmp <- NormalizeData(seu_tmp, verbose=F) %>% FindVariableFeatures(verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(npcs=N_PCS_TO_USE, verbose=F)
  sweep.res.list <- paramSweep(seu_tmp, PCs = 1:N_PCS_TO_USE, sct = F); sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
  bcmvn <- find.pK(sweep.stats); bcmvn$pK <- as.numeric(as.character(bcmvn$pK)); final_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  nExp_val <- round(ncol(seu_tmp) * DOUBLET_RATE)
  seu_tmp <- doubletFinder(seu_tmp, PCs = 1:N_PCS_TO_USE, pK = final_pk, nExp = nExp_val, sct = FALSE)
  res_col <- tail(grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE), 1)
  results_list[[sample_name]] <- seu_tmp@meta.data[, res_col, drop = FALSE]
}
all_res <- do.call(rbind, lapply(results_list, function(df) setNames(df, "Doublet_Status")))
data$Doublet_Status <- all_res[rownames(data@meta.data), "Doublet_Status"]
data <- subset(data, subset = Doublet_Status != "Doublet" | is.na(Doublet_Status))
rm(data_list, results_list, all_res); gc()

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
    rm(decontx_results, counts_sparse); gc()
}

# --- 2.5: Processing, Integration, and Clustering ---
message("\n--- Step 2.5: Processing, Integration, and Clustering ---")
data <- NormalizeData(data, verbose=F) %>% FindVariableFeatures(nfeatures=N_VARIABLE_FEATURES, verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(npcs=N_PCS_TO_USE, verbose=F)

# Always run Harmony and save both results for comparison and use
data <- RunHarmony(data, group.by.vars = "SampleID", reduction = "pca", reduction.save = "harmony", verbose = FALSE)
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "harmony", graph.name = "harmony_nn", verbose = FALSE)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, graph.name = "harmony_nn", cluster.name = "clusters_harmony", verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "harmony", n.neighbors = UMAP_N_NEIGHBORS, min.dist=UMAP_MIN_DIST, reduction.name = "umap_harmony", verbose = FALSE)

# --- 2.6: Generate Outputs for Annotation ---
message("\n--- Step 2.6: Generating outputs for annotation ---")
# Plot 1: UMAP of clusters
p_clust <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE) + NoLegend()
ggsave(file.path(OUTPUT_DIR, "ANNOTATE_THIS_UMAP.png"), plot = p_clust, width = 8, height = 7)

# Plot 2: DotPlot of markers vs clusters
dot_plot <- DotPlot(data, features = unique(unlist(BROAD_MARKERS_LIST)), group.by = "clusters_harmony", dot.min = 0.05, cols = 'RdBu') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "ANNOTATE_USING_THIS_DOTPLOT.png"), plot = dot_plot, width = 14, height = 10, bg = "white")

# Save the processed object for the next script
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_for_annotation.rds")))

message(paste0("\n\n--- PROCESSING COMPLETE ---",
               "\n\nNext Steps:",
               "\n1. Open the '", OUTPUT_DIR, "' folder.",
               "\n2. Examine 'ANNOTATE_THIS_UMAP.png' and 'ANNOTATE_USING_THIS_DOTPLOT.png'.",
               "\n3. Open the '02_annotate_data.R' script and fill in your cell type annotations.",
               "\n4. Run '02_annotate_data.R' to complete the analysis."))
# =============================================================================
