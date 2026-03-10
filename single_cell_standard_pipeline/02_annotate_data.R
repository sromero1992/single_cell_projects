# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 2: INTERACTIVE ANNOTATION STUDIO
#
# PURPOSE:
# This script is an interactive environment for biologists to explore and
# annotate the processed data from Script 1. You can dynamically generate
# plots, choose between different data integrations (Harmony vs. none),
# and then finalize your cell type labels.
#
# =============================================================================
# --- PART 1: SETUP & CONFIGURATION ---
library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(writexl)
set.seed(123)

# >>> ACTION 1: Define project name and ROOT path. These MUST match Script 1. <<<
PROJECT_NAME <- "Wu_IEC_Project"
ROOT_PATH <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"

# --- Automatic Path Generation (Do not edit) ---
OUTPUT_DIR <- file.path(ROOT_PATH, "seurat_output")
data <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds")))

# >>> ACTION 2: Define your marker gene lists here. <<<
# You can add, remove, or change genes and cell type names as needed.
BROAD_MARKERS_LIST <- list(
  "Plasma B cells" = c("Ighg1", "Mzb1"), "B cells" = c("Ms4a1", "Cd19"),
  "T cells" = c("Cd3e", "Cd4", "Cd8a", "Trbc2"), "cDCs" = c("Clec9a", "Xcr1"),
  "Macrophages" = c("Cd163", "Csf1r", "Cd68", "Mrc1", "Aif1"),
  "Fibroblasts" = c("Col1a1", "Col3a1", "Col6a2"), "ENs" = c("Vip", "Tubb3", "Vat1l"),
  "EGCs" = c("Sox10", "Plp1"), "LECs" = c("Lyve1", "Flt4", "Pecam1"),
  "VECs" = c("Plvap", "Flt1"), "Colonocytes" = c("Epcam", "Muc2", "Krt19", "Krt20", "Vil1"),
  "SMCs" = c("Myh11", "Tagln", "Hhip", "Des"), "Adipocytes" = c("Adipoq", "Plin1", "Plin4")
)

SUB_MARKERS_LIST <- list(
  "Abs. colonocytes" = c("Alpi", "Ces2c", "Slc26a2", "Ceacam1", "Aqp8"),
  "Goblet cells" = c("Muc2", "Spink4", "Agr2", "Tff3"),
  "EECs" = c("Chga", "Chgb", "Tph1", "Cck", "Scgn", "Scg5"),
  "Tuft cells" = c("Dclk1", "Trpm5", "Avil", "Sh2d6", "Plcg2"),
  "TA cells" = c("Mki67", "Top2a", "Birc5", "Pcna", "Stmn1"),
  "Stem cells" = c("Lgr5", "Lrig1", "Ascl2", "Slc12a2", "Smoc2", "Kcnq1", "Gpx2", "Ephb2", "Bmpr1a", "Hopx", "Sox9")
)

# =============================================================================
# --- PART 2: INTERACTIVE ANNOTATION WORKFLOW ---
#
# HOW TO USE THIS SECTION:
# 1. Set the 'VIEW_MODE' variable below to either "harmony" or "no_harmony".
#    - "harmony" (Recommended): Uses the batch-corrected data. Best for comparing across samples.
#    - "no_harmony": Uses the standard PCA data. Useful for diagnosing issues.
# 2. Run the code block under "Generate Plots". This will display the UMAP and DotPlot
#    in your RStudio 'Plots' panel.
# 3. Examine the plots.
# 4. Fill in the 'BROAD_ANNOTATION_MAP' with your cell type assignments.
# 5. When you are satisfied, run the rest of the script from Part 3 onwards.
# -----------------------------------------------------------------------------

# >>> ACTION 3: CHOOSE YOUR VIEW MODE <<<
VIEW_MODE <- "harmony"  # Options: "harmony" or "no_harmony"

# --- Set up variables based on view mode (No editing needed here) ---
if (VIEW_MODE == "harmony") {
  cluster_col <- "clusters_harmony"
  umap_col <- "umap_harmony"
  plot_title_suffix <- "(Harmony Corrected)"
} else {
  cluster_col <- "clusters_none"
  umap_col <- "umap_none"
  plot_title_suffix <- "(No Harmony)"
}

# --- Generate Plots (Run this block to see your UMAP and DotPlot) ---
# UMAP Plot
p_umap <- DimPlot(data, reduction = umap_col, group.by = cluster_col, label = TRUE) +
  ggtitle(paste("UMAP of Clusters", plot_title_suffix)) + NoLegend()

# DotPlot
p_dotplot <- DotPlot(data, features = unique(unlist(BROAD_MARKERS_LIST)), group.by = cluster_col, dot.min = 0.05, cols = 'RdBu', scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste("Marker Expression", plot_title_suffix))

# Display plots on screen
print(p_umap)
print(p_dotplot)

# =============================================================================
# --- PART 3: FINALIZE ANNOTATION & SAVE (ACTION REQUIRED) ---
#
# After examining the plots generated above, fill in the annotation map.
# Make sure the cluster IDs match the view mode you chose (e.g., 'clusters_harmony').
# -----------------------------------------------------------------------------

# >>> ACTION 4: FILL IN YOUR FINAL ANNOTATIONS HERE <<<
BROAD_ANNOTATION_MAP <- c(
  '0' = 'Colonocytes', '1' = 'Colonocytes', '2' = 'Colonocytes', '3' = 'Colonocytes', '4' = 'Colonocytes', '5' = 'Colonocytes', '6' = 'Colonocytes', '7' = 'Colonocytes', '8' = 'Plasma B cells', '9' = 'Colonocytes', '10' = 'Colonocytes', '11' = 'B cells', '12' = 'Colonocytes', '13' = 'Colonocytes', '14' = 'T cells', '15' = 'Macrophages', '16' = 'Fibroblasts', '17' = 'SMCs', '18' = 'Colonocytes', '19' = 'T cells', '20' = 'SMCs', '21' = 'Fibroblasts', '22' = 'VECs', '23' = 'cDCs', '24' = 'Colonocytes', '25' = 'SMCs', '26' = 'LECs', '27' = 'Plasma B cells', '28' = 'Colonocytes', '29' = 'SMCs', '30' = 'Plasma B cells', '31' = 'SMCs', '32' = 'SMCs', '33' = 'Colonocytes', '34' = 'Colonocytes', '35' = 'Colonocytes', '36' = 'Colonocytes', '37' = 'Colonocytes', '38' = 'EGCs', '39' = 'B cells', '40' = 'B cells', '41' = 'Fibroblasts', '42' = 'ENs', '43' = 'Colonocytes', '44' = 'Adipocytes', '45' = 'LECs', '46' = 'Colonocytes', '47' = 'Colonocytes', '48' = 'Colonocytes', '49' = 'Colonocytes', '50' = 'SMCs', '51' = 'LECs', '52' = 'T cells', '53' = 'Colonocytes'
)

# --- Apply annotations and save final plots (No editing needed below) ---
data$broad_cell_types <- recode_factor(data[[cluster_col, drop=TRUE]], !!!BROAD_ANNOTATION_MAP)

# Final UMAP and DotPlot grouped by your new cell type names
p_final_annot <- DimPlot(data, reduction = umap_col, group.by = "broad_cell_types", label = T, repel = T) + ggtitle(paste("Final Annotated Cell Types", plot_title_suffix))
p_final_dot <- DotPlot(data, features = unique(unlist(BROAD_MARKERS_LIST)), group.by = "broad_cell_types", dot.min=0.05, cols='RdBu', scale = T) + coord_flip()

ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotated.png"), plot = p_final_annot, width = 10, height = 8)
ggsave(file.path(OUTPUT_DIR, "FINAL_DotPlot_annotated.png"), plot = p_final_dot, width = 8, height = 12, bg="white")

message("--- Broad cell type annotation complete and saved! ---")


# ###############################################################
# ###            COMPLETE SECTION STARTS HERE                 ###
# ###############################################################

# =============================================================================
# --- PART 4: DOWNSTREAM ANALYSIS ---
#
# This section runs compositional analysis and an optional sub-clustering
# workflow on your newly annotated data.
# =============================================================================

# --- 4.1: Compositional Analysis ---

# >>> ACTION 5: Specify additional metadata columns for proportion plots <<<
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_Diet") # You can define your own classes for this 

message("\n--- Starting Compositional Analysis ---")

# Function to plot proportions
plot_cell_proportions <- function(seurat_obj, cluster_col, group_col, output_dir) {
    if (!group_col %in% colnames(seurat_obj@meta.data)) {
        # Auto-create combined columns if they don't exist (e.g., "Genotype_Diet")
        if (grepl("_", group_col)) {
            cols_to_combine <- strsplit(group_col, "_")[[1]]
            if (all(cols_to_combine %in% colnames(seurat_obj@meta.data))) {
                message(paste("Creating combined grouping column:", group_col))
                seurat_obj[[group_col]] <- apply(seurat_obj@meta.data[, cols_to_combine], 1, paste, collapse = "_")
            } else {
                warning(paste("Cannot create '", group_col, "' as its component columns are missing. Skipping."))
                return(NULL)
            }
        } else {
             warning(paste("Grouping column '", group_col, "' not found in metadata. Skipping."))
             return(NULL)
        }
    }
    
    df <- seurat_obj@meta.data %>%
        group_by(!!sym(group_col), !!sym(cluster_col)) %>% summarise(n = n(), .groups = "drop") %>%
        group_by(!!sym(group_col)) %>% mutate(percentage = (n / sum(n)) * 100) %>% ungroup()
    
    output_prefix <- paste0("by_", group_col)
    write_xlsx(df, path = file.path(output_dir, paste0("FINAL_Stats_", output_prefix, ".xlsx")))
    
    p <- ggplot(df, aes(x = !!sym(group_col), y = percentage, fill = !!sym(cluster_col))) +
        geom_bar(stat = "identity", position = "stack", color = "white") +
        labs(title = paste("Cell Proportions by", group_col), y = "Percentage (%)", x = group_col) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title=element_blank())
    
    ggsave(file.path(output_dir, paste0("FINAL_Proportions_", output_prefix, ".png")), p, width = max(8, 1.5 * n_distinct(seurat_obj@meta.data[[group_col]])), height = 7, bg = "white")
    message(paste("Generated proportion plot and stats for group:", group_col))
}

# 1. Always plot by SampleID
plot_cell_proportions(seurat_obj = data, cluster_col = "broad_cell_types", group_col = "SampleID", output_dir = OUTPUT_DIR)

# 2. Plot any additional user-defined groups
if (!is.null(ADDITIONAL_GROUPS_TO_PLOT)) {
  for (group in ADDITIONAL_GROUPS_TO_PLOT) {
    plot_cell_proportions(seurat_obj = data, cluster_col = "broad_cell_types", group_col = group, output_dir = OUTPUT_DIR)
  }
}
message("--- Compositional analysis complete. ---")


# --- 4.2: Sub-Clustering Workflow (Optional) ---
# >>> ACTION 6: Configure and enable sub-clustering <<<
RUN_SUBCLUSTERING <- TRUE
CELL_TYPE_TO_SUBCLUSTER <- "Colonocytes" 

if (RUN_SUBCLUSTERING) {
    message(paste("\n--- Starting Sub-clustering for:", CELL_TYPE_TO_SUBCLUSTER, "---"))
    
    # 1. Isolate the cell type and re-cluster at high resolution
    data_sub <- subset(data, subset = broad_cell_types == CELL_TYPE_TO_SUBCLUSTER)
    data_sub@reductions <- list(); data_sub@graphs <- list()
    data_sub <- FindVariableFeatures(data_sub, nfeatures=2000, verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(npcs=50, verbose=F)
    data_sub <- RunHarmony(data_sub, "SampleID", reduction.save="harmony_sub", verbose=F)
    data_sub <- FindNeighbors(data_sub, reduction="harmony_sub", dims=1:50, verbose=F) %>%
                FindClusters(resolution=3.0, cluster.name="sub_clusters", verbose=F) %>%
                RunUMAP(reduction="harmony_sub", dims=1:50, reduction.name="umap_sub", verbose=F)

    # 2. Generate new UMAP and DotPlot for sub-cluster annotation
    p_sub_clust <- DimPlot(data_sub, reduction = "umap_sub", group.by = "sub_clusters", label = TRUE) + NoLegend()
    dot_plot_sub <- DotPlot(data_sub, features = unique(unlist(SUB_MARKERS_LIST)), group.by = "sub_clusters", dot.min=0.05, cols='RdBu') +
                    theme(axis.text.x = element_text(angle=45, hjust=1))

    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_ANNOTATE_THIS_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")), p_sub_clust, width=8, height=7)
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_ANNOTATE_USING_THIS_DOTPLOT_", CELL_TYPE_TO_SUBCLUSTER, ".png")), dot_plot_sub, width=12, height=8, bg="white")
    
    # 3. User provides sub-cluster annotations
    # >>> ACTION 7: FILL IN SUB-CLUSTER ANNOTATIONS HERE <<<
    SUB_ANNOTATION_MAP <- c('0'='Goblet cells', '1'='Abs. colonocytes', '2'='Abs. colonocytes', '3'='Goblet cells', '4'='Stem cells', '5'='Goblet cells', '6'='Abs. colonocytes', '7'='Goblet cells', '8'='Abs. colonocytes', '9'='Abs. colonocytes', '10'='TA cells', '11'='Abs. colonocytes', '12'='Abs. colonocytes', '13'='TA cells', '14'='TA cells', '15'='Abs. colonocytes', '16'='Abs. colonocytes', '17'='TA cells', '18'='Stem cells', '19'='Abs. colonocytes', '20'='Goblet cells', '21'='Abs. colonocytes', '22'='Goblet cells', '23'='Abs. colonocytes', '24'='Goblet cells', '25'='Goblet cells', '26'='Goblet cells', '27'='Abs. colonocytes', '28'='Abs. colonocytes', '29'='TA cells', '30'='Abs. colonocytes', '31'='Abs. colonocytes', '32'='Stem cells', '33'='Abs. colonocytes', '34'='Goblet cells', '35'='Abs. colonocytes', '36'='TA cells', '37'='Abs. colonocytes', '38'='TA cells', '39'='Abs. colonocytes', '40'='EECs', '41'='Goblet cells', '42'='TA cells', '43'='Goblet cells', '44'='Goblet cells', '45'='Tuft cells', '46'='Tuft cells', '47'='EECs', '48'='Goblet cells', '49'='Abs. colonocytes', '50'='Abs. colonocytes', '51'='Goblet cells', '52'='Abs. colonocytes', '53'='Goblet cells', '54'='Goblet cells', '55'='Goblet cells', '56'='Abs. colonocytes', '57'='Goblet cells', '58'='Stem cells', '59'='Goblet cells', '60'='Abs. colonocytes', '61'='TA cells', '62'='Abs. colonocytes')
    
    # 4. Apply annotations and save final sub-cluster results
    data_sub$sub_cell_types <- recode_factor(data_sub$sub_clusters, !!!SUB_ANNOTATION_MAP)
    p_sub_annot <- DimPlot(data_sub, reduction = "umap_sub", group.by = "sub_cell_types", label = TRUE, repel = TRUE)
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_SUBCLUSTER_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")), p_sub_annot, width = 10, height = 8)
    saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_", CELL_TYPE_TO_SUBCLUSTER, "_subset_annotated.rds")))
    message("--- Sub-clustering analysis complete. ---")
}

# --- Final save of the main object with all annotations ---
saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds")))

message("\n\n--- ANALYSIS COMPLETE! All final data and plots are in the output directory. ---")
