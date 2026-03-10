# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 2: ANNOTATION & DOWNSTREAM ANALYSIS
#
# PURPOSE:
# This script is for the biologist. It loads the processed data from Script 1
# and allows you to perform the crucial step of manual cell type annotation.
# It then generates final plots, performs compositional analysis, and can
# run an optional sub-clustering workflow.
#
# =============================================================================
# --- PART 1: SETUP ---
# Load libraries and define project paths (MUST MATCH SCRIPT 1)
# -----------------------------------------------------------------------------
library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(writexl)
set.seed(123)

# Define project name and output directory (must be the same as in 01_process_data.R)
PROJECT_NAME <- "Wu_IEC_Project"
OUTPUT_DIR <- "Z:/selim_working_dir/2025_wu_iec_project/r_process/seurat_output"

# Load the processed Seurat object from Script 1
data <- readRDS(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_for_annotation.rds")))

# =============================================================================
# --- PART 2: BROAD CELL TYPE ANNOTATION (ACTION REQUIRED) ---
#
# HOW TO ANNOTATE YOUR DATA:
#
# 1. BROWSE THE PLOTS:
#    Go to your output folder ('seurat_output') and open these two files:
#    - 'ANNOTATE_THIS_UMAP.png': Shows where each cluster is located spatially.
#    - 'ANNOTATE_USING_THIS_DOTPLOT.png': This is your main tool!
#
# 2. READ THE DOT PLOT:
#    - Each ROW is a cell type defined by its marker genes (e.g., "T cells").
#    - Each COLUMN is a numbered cluster from the UMAP.
#    - The COLOR (e.g., red) shows the AVERAGE EXPRESSION of the markers in that cluster.
#    - The SIZE of the dot shows the PERCENTAGE of cells in that cluster expressing the markers.
#    - A large, bright red dot means that cluster strongly represents that cell type.
#
# 3. FILL IN THE MAP:
#    Based on the DotPlot, match each cluster number to a cell type name below.
# -----------------------------------------------------------------------------

# >>> ACTION: FILL IN YOUR ANNOTATIONS HERE <<<
BROAD_ANNOTATION_MAP <- c(
  '0' = 'Colonocytes', '1' = 'Colonocytes', '2' = 'Colonocytes', '3' = 'Colonocytes',
  '4' = 'Colonocytes', '5' = 'Colonocytes', '6' = 'Colonocytes', '7' = 'Colonocytes',
  '8' = 'Plasma B cells', '9' = 'Colonocytes', '10' = 'Colonocytes', '11' = 'B cells',
  '12' = 'Colonocytes', '13' = 'Colonocytes', '14' = 'T cells', '15' = 'Macrophages',
  '16' = 'Fibroblasts', '17' = 'SMCs', '18' = 'Colonocytes', '19' = 'T cells',
  '20' = 'SMCs', '21' = 'Fibroblasts', '22' = 'VECs', '23' = 'cDCs', '24' = 'Colonocytes',
  '25' = 'SMCs', '26' = 'LECs', '27' = 'Plasma B cells', '28' = 'Colonocytes', '29' = 'SMCs',
  '30' = 'Plasma B cells', '31' = 'SMCs', '32' = 'SMCs', '33' = 'Colonocytes', '34' = 'Colonocytes',
  '35' = 'Colonocytes', '36' = 'Colonocytes', '37' = 'Colonocytes', '38' = 'EGCs', '39' = 'B cells',
  '40' = 'B cells', '41' = 'Fibroblasts', '42' = 'ENs', '43' = 'Colonocytes', '44' = 'Adipocytes',
  '45' = 'LECs', '46' = 'Colonocytes', '47' = 'Colonocytes', '48' = 'Colonocytes',
  '49' = 'Colonocytes', '50' = 'SMCs', '51' = 'LECs', '52' = 'T cells', '53' = 'Colonocytes'
)

# --- Apply annotations and generate final plots (No editing needed below) ---
data$broad_cell_types <- recode_factor(data$clusters_harmony, !!!BROAD_ANNOTATION_MAP)
p_final_annot <- DimPlot(data, reduction = "umap_harmony", group.by = "broad_cell_types", label = TRUE, repel = TRUE)
ggsave(file.path(OUTPUT_DIR, "FINAL_UMAP_annotated.png"), plot = p_final_annot, width = 10, height = 8, dpi = 300)
message("--- Broad cell type annotation complete! ---")


# =============================================================================
# --- PART 3: COMPOSITIONAL ANALYSIS (REVISED & GENERALIZED) ---
#
# This section calculates the percentage of each cell type.
# - A plot grouped by 'SampleID' is always created.
# - You can specify additional metadata columns to create more plots.
# -----------------------------------------------------------------------------

# >>> ACTION: Specify additional metadata columns for proportion plots <<<
#
# List any other columns from your metadata file you want to group by.
# For example: c("Condition", "Timepoint")
#
# If you created a combined column like 'Genotype_Diet' in your Excel file,
# you can list it here.
#
# If you have no other groups, set this to NULL, like this:
# ADDITIONAL_GROUPS_TO_PLOT <- NULL
#
ADDITIONAL_GROUPS_TO_PLOT <- c("Genotype_Diet") # Example: c("Condition", "Genotype_Diet")


# --- Plotting function and execution (No editing needed below) ---
message("\n--- Starting Compositional Analysis ---")

# Function to plot proportions
plot_cell_proportions <- function(seurat_obj, cluster_col, group_col, output_dir) {
    if (!group_col %in% colnames(seurat_obj@meta.data)) {
        warning(paste("Grouping column '", group_col, "' not found in metadata. Skipping this plot."))
        return(NULL)
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
    # If a combined group like "Genotype_Diet" doesn't exist, try to create it
    if (!group %in% colnames(data@meta.data) && grepl("_", group)) {
        cols_to_combine <- strsplit(group, "_")[[1]]
        if (all(cols_to_combine %in% colnames(data@meta.data))) {
            message(paste("Creating combined grouping column:", group))
            data[[group]] <- apply(data@meta.data[, cols_to_combine], 1, paste, collapse = "_")
        }
    }
    plot_cell_proportions(seurat_obj = data, cluster_col = "broad_cell_types", group_col = group, output_dir = OUTPUT_DIR)
  }
}

message("--- Compositional analysis complete. ---")


# =============================================================================
# --- PART 4: SUB-CLUSTERING (OPTIONAL) ---
# (This part remains the same as the previous version)
# -----------------------------------------------------------------------------
# --- A. Configure Sub-clustering ---
RUN_SUBCLUSTERING <- TRUE
CELL_TYPE_TO_SUBCLUSTER <- "Colonocytes" 

# --- B. Execute Sub-clustering (No editing needed below if RUN_SUBCLUSTERING is TRUE) ---
if (RUN_SUBCLUSTERING) {
    message(paste("\n--- Starting Sub-clustering for:", CELL_TYPE_TO_SUBCLUSTER, "---"))
    data_sub <- subset(data, subset = broad_cell_types == CELL_TYPE_TO_SUBCLUSTER)
    data_sub@reductions <- list(); data_sub@graphs <- list()
    data_sub <- FindVariableFeatures(data_sub, nfeatures=2000, verbose=F) %>% ScaleData(verbose=F) %>% RunPCA(npcs=50, verbose=F)
    data_sub <- RunHarmony(data_sub, "SampleID", reduction.save="harmony_sub", verbose=F)
    data_sub <- FindNeighbors(data_sub, reduction="harmony_sub", dims=1:50, verbose=F) %>%
                FindClusters(resolution=3.0, cluster.name="sub_clusters", verbose=F) %>%
                RunUMAP(reduction="harmony_sub", dims=1:50, reduction.name="umap_sub", verbose=F)
    
    p_sub_clust <- DimPlot(data_sub, reduction = "umap_sub", group.by = "sub_clusters", label = TRUE)
    ggsave(file.path(OUTPUT_DIR, paste0("SUBCLUSTER_ANNOTATE_THIS_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")), p_sub_clust, width=8, height=7)
    
    # >>> ACTION: FILL IN SUB-CLUSTER ANNOTATIONS HERE <<<
    SUB_ANNOTATION_MAP <- c('0'='Goblet cells', '1'='Abs. colonocytes', '2'='Abs. colonocytes', '3'='Goblet cells', '4'='Stem cells', '5'='Goblet cells', '6'='Abs. colonocytes', '7'='Goblet cells', '8'='Abs. colonocytes', '9'='Abs. colonocytes', '10'='TA cells', '11'='Abs. colonocytes', '12'='Abs. colonocytes', '13'='TA cells', '14'='TA cells', '15'='Abs. colonocytes', '16'='Abs. colonocytes', '17'='TA cells', '18'='Stem cells', '19'='Abs. colonocytes', '20'='Goblet cells', '21'='Abs. colonocytes', '22'='Goblet cells', '23'='Abs. colonocytes', '24'='Goblet cells', '25'='Goblet cells', '26'='Goblet cells', '27'='Abs. colonocytes', '28'='Abs. colonocytes', '29'='TA cells', '30'='Abs. colonocytes', '31'='Abs. colonocytes', '32'='Stem cells', '33'='Abs. colonocytes', '34'='Goblet cells', '35'='Abs. colonocytes', '36'='TA cells', '37'='Abs. colonocytes', '38'='TA cells', '39'='Abs. colonocytes', '40'='EECs', '41'='Goblet cells', '42'='TA cells', '43'='Goblet cells', '44'='Goblet cells', '45'='Tuft cells', '46'='Tuft cells', '47'='EECs', '48'='Goblet cells', '49'='Abs. colonocytes', '50'='Abs. colonocytes', '51'='Goblet cells', '52'='Abs. colonocytes', '53'='Goblet cells', '54'='Goblet cells', '55'='Goblet cells', '56'='Abs. colonocytes', '57'='Goblet cells', '58'='Stem cells', '59'='Goblet cells', '60'='Abs. colonocytes', '61'='TA cells', '62'='Abs. colonocytes')
    
    data_sub$sub_cell_types <- recode_factor(data_sub$sub_clusters, !!!SUB_ANNOTATION_MAP)
    p_sub_annot <- DimPlot(data_sub, reduction = "umap_sub", group.by = "sub_cell_types", label = TRUE, repel = TRUE)
    ggsave(file.path(OUTPUT_DIR, paste0("FINAL_SUBCLUSTER_UMAP_", CELL_TYPE_TO_SUBCLUSTER, ".png")), p_sub_annot, width = 10, height = 8)
    saveRDS(data_sub, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_", CELL_TYPE_TO_SUBCLUSTER, "_subset_annotated.rds")))
    message("--- Sub-clustering analysis complete. ---")
}

saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds")))
message("\n\n--- ANALYSIS COMPLETE! ---")
