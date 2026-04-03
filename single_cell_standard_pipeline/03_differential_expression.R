# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 3: DIFFERENTIAL EXPRESSION ANALYSIS
# Version: 2.0 (SplineDV + Wilcoxon DE; no pseudo-bulk)
#
# PURPOSE:
#   Performs differential gene expression (DGE) and differential variability
#   (DV) analysis between user-defined groups across annotated cell populations.
#
# ANALYSIS MODES:
#   MODE A  — Single-cell Wilcoxon (FindMarkers): Fast exploratory DE across
#             all cell types. Use for DISCOVERY and supplementary plots.
#             NOTE: p-values are inflated by pseudoreplication (each cell is
#             treated as an independent observation, but cells from the same
#             animal are not). These are hypothesis-generating, not confirmatory.
#
#   MODE B  — Cell-type marker genes (FindAllMarkers): Identifies genes
#             enriched in each broad cell type vs. all others. Useful for
#             validating annotations and generating cell-type gene signatures.
#
#   MODE C  — Within-subtype Wilcoxon: Loads the sub-cluster RDS from
#             Script 02 and runs Mode A at the sub-type level.
#
#   MODE D  — SplineDV + DE Overlap + fgsea ORA Enrichment: For each cell
#             type, runs SplineDV to find differentially variable (DV) genes,
#             then overlaps DV genes with up- and down-regulated DE genes from
#             Mode A, and runs fgsea over-representation analysis (fora) on
#             each overlap gene set against MSigDB pathways.
#
#             Conceptual logic:
#               DV genes (SplineDV)  = genes where VARIANCE changes between groups
#               Up-DE genes (Mode A) = genes where MEAN goes UP in GROUP_1
#               Down-DE genes (Mode A) = genes where MEAN goes DOWN in GROUP_1
#
#               Overlap A = DV ∩ Up-DE  → pathway enrichment → "activated & more variable"
#               Overlap B = DV ∩ Down-DE → pathway enrichment → "repressed & more variable"
#
#             Biological interpretation:
#               Genes appearing in BOTH DV and DE have the strongest evidence
#               of regulation: their expression level AND their cell-to-cell
#               variability both change. These are priority candidates.
#
# METADATA NOTE:
#   All scripts use the metadata Excel file to group cells. The CONDITION_COLUMN
#   and SAMPLE_COLUMN values must exactly match column names in that Excel file.
#   Any column from the metadata (Condition, Genotype, Diet, Sex, Timepoint, etc.)
#   is usable here for subsetting or comparison. See Section 1.3 for details.
#
# DEPENDENCIES:
#   Core:    Seurat, dplyr, ggplot2, ggrepel, writexl, patchwork
#   Mode D:  SplineDV (GitHub: Xenon8778/SplineDV), fgsea, msigdbr, stringr
#
# NEXT STEP: Script 04 — Cell Scoring & Pathway Enrichment.
# =============================================================================

# --- Load Core Libraries -----------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(writexl)
library(patchwork)
library(Matrix)

set.seed(123)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) ---
# =============================================================================

# --- 1.1: Project Paths (must match Scripts 01 and 02) -----------------------
PROJECT_NAME <- "Nr4a1_Study17_Project"
ROOT_PATH    <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"
#ROOT_PATH   <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"

OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")
DE_DIR       <- file.path(OUTPUT_DIR, "differential_expression")
if (!dir.exists(DE_DIR)) dir.create(DE_DIR, recursive = TRUE)

MAIN_RDS <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))

# --- 1.2: Analysis Mode Toggles ----------------------------------------------
RUN_MODE_A_WILCOXON          <- TRUE   # Single-cell Wilcoxon DE (exploratory)
RUN_MODE_B_BETWEEN_CELLTYPES <- FALSE  # Cell-type marker genes (FindAllMarkers)
RUN_MODE_C_WITHIN_SUBTYPE    <- FALSE  # DE at sub-cluster level
RUN_MODE_D_SPLINEDV          <- TRUE   # Differential variability + DE overlap + enrichment

# --- 1.3: Comparison Definition & Metadata Columns --------------------------
#
# HOW METADATA COLUMNS WORK IN THIS PIPELINE:
#   The metadata Excel file (defined in Script 01, METADATA_FILE parameter) is
#   the single source of truth for all sample-level metadata. Any column in
#   that Excel file is automatically attached to every cell from that sample.
#
#   For example, if your Excel file has columns:
#     SampleID | Condition | Sex | Diet | Genotype | Timepoint
#   Then after Script 01 runs, your Seurat object will have all of those as
#   metadata columns accessible here.
#
#   IMPORTANT — Matching: Column names here must EXACTLY match the Excel header.
#   IMPORTANT — Concatenated columns: Script 02 auto-creates combined columns
#     like "Genotype_Diet" from "Genotype" + "Diet" if both exist. You can
#     use those combined columns in CONDITION_COLUMN below.
#
# CELLTYPE_COLUMN: Metadata column containing broad cell type labels.
#   Created by Script 02 (always named "CellType" by default).
CELLTYPE_COLUMN  <- "CellType"

# CONDITION_COLUMN: The metadata column that defines your two comparison groups.
#   Examples: "Condition" (KO vs WT), "Diet" (inulin vs cellulose),
#             "Genotype" (Nr4a1_KO vs WT), "Genotype_Diet" (combined).
CONDITION_COLUMN <- "Condition"

# GROUP_1: The "test" group. Positive log2FC means higher in GROUP_1.
# GROUP_2: The "reference" group. Must match values in CONDITION_COLUMN.
GROUP_1 <- "KO"
GROUP_2 <- "WT"

# SAMPLE_COLUMN: Column identifying biological replicates.
#   Used in Mode D to compute per-sample DV statistics.
#   Must be the same column used as SAMPLE_COLUMN in Script 01 (typically "SampleID").
SAMPLE_COLUMN <- "SampleID"

# CELL_TYPES_TO_TEST: Character vector of cell types to analyse.
#   Set to NULL to automatically run on ALL cell types found in CELLTYPE_COLUMN.
#   Example: c("T cells", "Macrophages", "Colonocytes")
CELL_TYPES_TO_TEST <- NULL

# MINIMUM_CELLS_PER_GROUP: Skip a cell type × group if either group has fewer cells.
MINIMUM_CELLS_PER_GROUP <- 20

# --- 1.4: DGE Parameters (Mode A) --------------------------------------------
DE_TEST         <- "wilcox"   # "wilcox", "t", or "MAST"
DE_MIN_PCT      <- 0.10       # Min fraction of cells expressing a gene
DE_LOGFC_THRESH <- 0.25       # Min |log2FC| to retain a gene in the results
DE_PADJ_THRESH  <- 0.05       # Adjusted p-value cutoff for "significant"
DE_LFC_THRESH   <- 0.5        # |log2FC| threshold for volcano plot coloring
DE_TOP_N_LABEL  <- 15         # Top N genes labeled on each volcano plot

# --- 1.5: Mode B Parameters --------------------------------------------------
MODE_B_MIN_PCT <- 0.25
MODE_B_LOGFC   <- 0.5
MODE_B_TOP_N   <- 20

# --- 1.6: Mode C Parameters --------------------------------------------------
# Sub-clustering object from Script 02. Must have been generated with
# RUN_SUBCLUSTERING = TRUE for the specified cell type.
SUBTYPE_CELL_TYPE <- "Colonocytes"
SUBTYPE_COLUMN    <- "sub_cell_types"
SUBTYPE_RDS       <- file.path(OUTPUT_DIR,
                                paste0(PROJECT_NAME, "_", SUBTYPE_CELL_TYPE, "_subclustered.rds"))

# --- 1.7: Mode D — SplineDV + Overlap + Enrichment Parameters ---------------

# DV_PADJ_THRESH: Adjusted p-value threshold for calling a gene "differentially
#   variable". SplineDV detects changes in variance (not mean expression).
DV_PADJ_THRESH <- 0.05

# Enrichment settings for fgsea::fora() over-representation analysis.
# The ORA tests whether DV ∩ DE overlap genes are enriched in MSigDB pathways,
# using all DE-tested genes as the background universe (more conservative
# than using the whole genome).
ENRICHMENT_SPECIES     <- "Mus musculus"  # or "Homo sapiens"
ENRICHMENT_CATEGORY    <- "C5"            # C5 = GO, H = Hallmark, C2 = KEGG/Reactome
ENRICHMENT_SUBCATEGORY <- "BP"            # BP = Biological Process
ENRICHMENT_MIN_SIZE    <- 5    # Min overlap size to test a pathway
ENRICHMENT_PADJ_THRESH <- 0.05 # Significance cutoff for ORA
ENRICHMENT_TOP_N_PLOT  <- 15   # Top N pathways shown in bar plots

# --- 1.8: Plot Settings -------------------------------------------------------
DPI_SETTING    <- 300
VOLCANO_WIDTH  <- 9
VOLCANO_HEIGHT <- 7

# =============================================================================
# --- PART 2: UTILITY FUNCTIONS -----------------------------------------------
# =============================================================================

# ---------------------------------------------------------------------------
#' run_de_wilcoxon: Per-cell-type Wilcoxon DE via Seurat FindMarkers.
#   EXPLORATORY ONLY — inflated p-values due to pseudoreplication.
# ---------------------------------------------------------------------------
run_de_wilcoxon <- function(seurat_obj, cell_type, celltype_col, condition_col,
                            group1, group2,
                            test.use = "wilcox", min.pct = 0.10,
                            logfc.threshold = 0.25) {
  cells_ct <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == cell_type]
  obj_ct   <- subset(seurat_obj, cells = cells_ct)

  n1 <- sum(obj_ct@meta.data[[condition_col]] == group1, na.rm = TRUE)
  n2 <- sum(obj_ct@meta.data[[condition_col]] == group2, na.rm = TRUE)
  if (n1 < MINIMUM_CELLS_PER_GROUP || n2 < MINIMUM_CELLS_PER_GROUP) {
    message(paste("  [SKIP]", cell_type, "| n(", group1, ")=", n1,
                  "n(", group2, ")=", n2, "< min (", MINIMUM_CELLS_PER_GROUP, ")"))
    return(NULL)
  }
  message(paste("  [Wilcoxon]", cell_type, "|", group1, "(n=", n1, ") vs", group2, "(n=", n2, ")"))

  Idents(obj_ct) <- condition_col
  result <- tryCatch(
    FindMarkers(obj_ct, ident.1 = group1, ident.2 = group2,
                test.use = test.use, min.pct = min.pct,
                logfc.threshold = logfc.threshold, verbose = FALSE),
    error = function(e) { message("    ERROR: ", e$message); return(NULL) }
  )
  if (is.null(result) || nrow(result) == 0) return(NULL)

  result$gene       <- rownames(result)
  result$cell_type  <- cell_type
  result$comparison <- paste0(group1, "_vs_", group2)
  result$method     <- "Wilcoxon_FindMarkers"
  result <- result %>%
    dplyr::rename(log2FC = avg_log2FC, padj = p_val_adj) %>%
    arrange(padj, desc(abs(log2FC)))
  return(result)
}

# ---------------------------------------------------------------------------
#' make_volcano_plot: Publication-ready volcano plot.
# ---------------------------------------------------------------------------
make_volcano_plot <- function(de_result, title,
                              padj_cut = DE_PADJ_THRESH,
                              lfc_cut  = DE_LFC_THRESH,
                              top_n    = DE_TOP_N_LABEL) {
  df <- de_result %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj + 1e-300),
      significance   = dplyr::case_when(
        padj < padj_cut & log2FC >  lfc_cut ~ "Up",
        padj < padj_cut & log2FC < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  n_up   <- sum(df$significance == "Up")
  n_down <- sum(df$significance == "Down")

  top_up   <- df %>% filter(significance == "Up")  %>% arrange(padj, desc(log2FC)) %>% head(ceiling(top_n / 2))
  top_down <- df %>% filter(significance == "Down") %>% arrange(padj, log2FC)       %>% head(floor(top_n / 2))

  ggplot(df, aes(x = log2FC, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.45, size = 1.4) +
    geom_text_repel(data = bind_rows(top_up, top_down), aes(label = gene),
                    size = 3.5, max.overlaps = 20, fontface = "italic",
                    box.padding = 0.4, point.padding = 0.3) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
               color = "black", linewidth = 0.6) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
               color = "black", linewidth = 0.6) +
    annotate("text", x =  max(abs(df$log2FC), na.rm = TRUE) * 0.75,
             y = max(df$neg_log10_padj, na.rm = TRUE) * 0.97,
             label = paste0("Up: ", n_up),   color = "#E41A1C", size = 4, fontface = "bold") +
    annotate("text", x = -max(abs(df$log2FC), na.rm = TRUE) * 0.75,
             y = max(df$neg_log10_padj, na.rm = TRUE) * 0.97,
             label = paste0("Down: ", n_down), color = "#377EB8", size = 4, fontface = "bold") +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
    labs(title = title,
         x     = expression(log[2] * "FC (" * Group[1] * "/" * Group[2] * ")"),
         y     = expression(-log[10] * "(adj. p-value)"),
         color = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          legend.position = "bottom")
}

# ---------------------------------------------------------------------------
#' save_de_results: Save DE result to Excel and volcano PNG.
# ---------------------------------------------------------------------------
save_de_results <- function(de_result, cell_type, label, out_dir) {
  if (is.null(de_result) || nrow(de_result) == 0) return(invisible(NULL))
  safe_ct  <- gsub("[^A-Za-z0-9_]", "_", cell_type)
  cmp      <- unique(de_result$comparison)[1]
  write_xlsx(de_result,
             file.path(out_dir, paste0(label, "_", safe_ct, "_", cmp, ".xlsx")))
  p_vol <- make_volcano_plot(de_result,
                             title = paste0(cell_type, "\n", cmp))
  ggsave(file.path(out_dir, paste0(label, "_Volcano_", safe_ct, "_", cmp, ".png")),
         p_vol, width = VOLCANO_WIDTH, height = VOLCANO_HEIGHT,
         dpi = DPI_SETTING, bg = "white")
  invisible(de_result)
}

# ---------------------------------------------------------------------------
#' run_splinedv: Wrapper around SplineDV to find differentially variable genes.
#
#   SplineDV (Xenon8778/SplineDV) fits a spline to the mean–variance
#   relationship per gene and tests whether variance differs between conditions.
#   This is COMPLEMENTARY to DE:
#     - DE asks: does the mean expression change?
#     - DV asks: does the cell-to-cell variability change?
#   A gene can be DV without being DE (same mean, different spread) or
#   both DV and DE (different mean AND spread — highest biological priority).
#
#   API: SplineDV::SplineDV(count_matrix_group1, count_matrix_group2)
#   Output: data.frame with gene, CV statistics, p-value, padj (FDR)
#   Docs: https://github.com/Xenon8778/SplineDV
# ---------------------------------------------------------------------------
run_splinedv <- function(seurat_obj, cell_type, celltype_col, condition_col,
                         group1, group2) {
  if (!requireNamespace("SplineDV", quietly = TRUE)) {
    stop("[ERROR] SplineDV not installed. Run 00_rlibs_installation.R to install it.")
  }
  library(SplineDV)

  cells_ct <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == cell_type]
  obj_ct   <- subset(seurat_obj, cells = cells_ct)
  obj_ct   <- subset(obj_ct, cells = rownames(obj_ct@meta.data)[
    obj_ct@meta.data[[condition_col]] %in% c(group1, group2)])

  n1 <- sum(obj_ct@meta.data[[condition_col]] == group1)
  n2 <- sum(obj_ct@meta.data[[condition_col]] == group2)
  if (n1 < MINIMUM_CELLS_PER_GROUP || n2 < MINIMUM_CELLS_PER_GROUP) {
    message(paste("  [SKIP SplineDV]", cell_type, "| insufficient cells"))
    return(NULL)
  }
  message(paste("  [SplineDV]", cell_type, "|", group1, "(n=", n1, ") vs", group2, "(n=", n2, ")"))

  cells_g1 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group1]
  cells_g2 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group2]

  # SplineDV expects raw count matrices (not normalized)
  counts_g1 <- GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g1, drop = FALSE]
  counts_g2 <- GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g2, drop = FALSE]

  dv_result <- tryCatch(
    SplineDV::SplineDV(counts_g1, counts_g2),
    error = function(e) {
      message("    [ERROR] SplineDV failed: ", e$message)
      message("    Verify API at: https://github.com/Xenon8778/SplineDV")
      return(NULL)
    }
  )
  if (is.null(dv_result) || nrow(dv_result) == 0) return(NULL)

  # Normalize column names across SplineDV versions
  if (!"gene" %in% colnames(dv_result) && "Gene" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, gene = Gene)
  if (!"gene" %in% colnames(dv_result))
    dv_result$gene <- rownames(dv_result)
  if (!"padj" %in% colnames(dv_result) && "adj.pvalue" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, padj = adj.pvalue)
  if (!"padj" %in% colnames(dv_result) && "FDR" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, padj = FDR)

  dv_result$cell_type  <- cell_type
  dv_result$comparison <- paste0(group1, "_vs_", group2)
  dv_result$n_cells_g1 <- n1
  dv_result$n_cells_g2 <- n2
  return(dv_result %>% arrange(padj))
}

# ---------------------------------------------------------------------------
#' run_fgsea_ora: Over-representation analysis via fgsea::fora().
#
#   Tests whether a gene list (the DV ∩ DE overlap) is enriched in MSigDB
#   pathways relative to the DE-tested gene universe (not the whole genome).
#   Using the DE universe as the background is CONSERVATIVE and appropriate
#   because genes not tested in DE were excluded due to low expression.
# ---------------------------------------------------------------------------
run_fgsea_ora <- function(query_genes, universe, pathways,
                          min_size    = ENRICHMENT_MIN_SIZE,
                          padj_thresh = ENRICHMENT_PADJ_THRESH) {
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("[ERROR] fgsea not installed. Run: BiocManager::install('fgsea')")
  library(fgsea)

  if (length(query_genes) < 2) {
    message("    [SKIP ORA] < 2 query genes — skipping.")
    return(NULL)
  }
  res <- tryCatch(
    fgsea::fora(pathways = pathways, genes = query_genes,
                universe = universe, minSize = min_size),
    error = function(e) { message("    [ERROR] fora: ", e$message); return(NULL) }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$overlapGenes <- sapply(res$overlapGenes, paste, collapse = "|")
  res %>% dplyr::filter(padj < padj_thresh) %>% dplyr::arrange(padj)
}

# ---------------------------------------------------------------------------
#' plot_ora_barplot: Horizontal bar plot of top enriched pathways.
# ---------------------------------------------------------------------------
plot_ora_barplot <- function(ora_result, title,
                             top_n      = ENRICHMENT_TOP_N_PLOT,
                             fill_color = "#E41A1C") {
  if (is.null(ora_result) || nrow(ora_result) == 0) return(NULL)

  df <- ora_result %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      pathway_label  = gsub("^GOBP_|^HALLMARK_|^KEGG_|^REACTOME_", "", pathway) %>%
        gsub("_", " ", .) %>% tools::toTitleCase() %>%
        stringr::str_wrap(width = 50),
      neg_log10_padj = -log10(padj + 1e-300)
    ) %>%
    dplyr::arrange(neg_log10_padj) %>%
    dplyr::mutate(pathway_label = factor(pathway_label, levels = pathway_label))

  ggplot(df, aes(x = neg_log10_padj, y = pathway_label)) +
    geom_col(fill = fill_color, alpha = 0.85, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = -log10(ENRICHMENT_PADJ_THRESH),
               linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_text(aes(label = paste0(overlap, "/", size)), hjust = -0.1, size = 3.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title    = title,
         x        = expression(-log[10] * "(adj. p-value)"),
         y        = NULL,
         subtitle = paste0(nrow(ora_result), " significant pathways | top ", top_n, " shown")) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
          axis.text.y   = element_text(size = 10),
          axis.text.x   = element_text(size = 10))
}

# =============================================================================
# --- PART 3: LOAD DATA -------------------------------------------------------
# =============================================================================
message("=== Loading annotated Seurat object ===")
data <- readRDS(MAIN_RDS)
message(paste("  Loaded:", ncol(data), "cells"))

if (is.null(CELL_TYPES_TO_TEST)) {
  CELL_TYPES_TO_TEST <- sort(as.character(unique(data@meta.data[[CELLTYPE_COLUMN]])))
  CELL_TYPES_TO_TEST <- CELL_TYPES_TO_TEST[!is.na(CELL_TYPES_TO_TEST)]
}
message(paste("  Cell types to analyse:", paste(CELL_TYPES_TO_TEST, collapse = ", ")))

# Pre-load MSigDB gene sets once (avoids repeated API calls)
msig_pathways <- NULL
if (RUN_MODE_D_SPLINEDV) {
  library(msigdbr)
  library(stringr)
  message("  Pre-loading MSigDB gene sets...")
  msig_df <- msigdbr(species = ENRICHMENT_SPECIES,
                     category = ENRICHMENT_CATEGORY,
                     subcategory = if (nchar(ENRICHMENT_SUBCATEGORY) > 0) ENRICHMENT_SUBCATEGORY else NULL)
  msig_pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
  message(paste("  Loaded", length(msig_pathways), "pathways from",
                ENRICHMENT_CATEGORY, "/", ENRICHMENT_SUBCATEGORY))
}

# =============================================================================
# --- PART 4: MODE A — WILCOXON DE (EXPLORATORY) ------------------------------
# =============================================================================
wilcox_results_all <- list()

if (RUN_MODE_A_WILCOXON) {
  message(paste("\n=== MODE A (Wilcoxon — Exploratory) |", GROUP_1, "vs", GROUP_2, "==="))

  mode_a_dir <- file.path(DE_DIR, paste0("ModeA_Wilcoxon_", GROUP_1, "_vs_", GROUP_2))
  if (!dir.exists(mode_a_dir)) dir.create(mode_a_dir, recursive = TRUE)

  for (ct in CELL_TYPES_TO_TEST) {
    result <- run_de_wilcoxon(data, ct, CELLTYPE_COLUMN, CONDITION_COLUMN,
                              GROUP_1, GROUP_2, test.use = DE_TEST,
                              min.pct = DE_MIN_PCT, logfc.threshold = DE_LOGFC_THRESH)
    if (is.null(result)) next
    wilcox_results_all[[ct]] <- result
    save_de_results(result, ct, "Wilcoxon", mode_a_dir)
  }

  if (length(wilcox_results_all) > 0) {
    write_xlsx(wilcox_results_all,
               file.path(mode_a_dir, paste0("Wilcoxon_ALL_", GROUP_1, "_vs_", GROUP_2, ".xlsx")))
    summary_w <- do.call(rbind, lapply(names(wilcox_results_all), function(ct) {
      df <- wilcox_results_all[[ct]]
      data.frame(CellType = ct,
                 N_Sig  = sum(df$padj < DE_PADJ_THRESH & abs(df$log2FC) > DE_LFC_THRESH, na.rm = TRUE),
                 N_Up   = sum(df$padj < DE_PADJ_THRESH & df$log2FC >  DE_LFC_THRESH, na.rm = TRUE),
                 N_Down = sum(df$padj < DE_PADJ_THRESH & df$log2FC < -DE_LFC_THRESH, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
    write_xlsx(summary_w, file.path(mode_a_dir, "Wilcoxon_DEG_summary.xlsx"))
    message("  Mode A complete.")
    print(summary_w)
  }
}

# =============================================================================
# --- PART 5: MODE B — CELL-TYPE MARKER GENES ---------------------------------
# =============================================================================
if (RUN_MODE_B_BETWEEN_CELLTYPES) {
  message("\n=== MODE B: Cell-Type Marker Genes ===")
  mode_b_dir <- file.path(DE_DIR, "ModeB_CellType_Markers")
  if (!dir.exists(mode_b_dir)) dir.create(mode_b_dir, recursive = TRUE)

  Idents(data) <- CELLTYPE_COLUMN
  all_markers  <- FindAllMarkers(data, only.pos = TRUE, min.pct = MODE_B_MIN_PCT,
                                 logfc.threshold = MODE_B_LOGFC, test.use = DE_TEST,
                                 verbose = FALSE)
  if (!is.null(all_markers) && nrow(all_markers) > 0) {
    all_markers <- all_markers %>%
      dplyr::rename(log2FC = avg_log2FC, padj = p_val_adj) %>%
      filter(padj < DE_PADJ_THRESH) %>%
      arrange(cluster, padj, desc(log2FC))
    top_markers <- all_markers %>% group_by(cluster) %>%
      slice_head(n = MODE_B_TOP_N) %>% ungroup()
    write_xlsx(split(all_markers, all_markers$cluster),
               file.path(mode_b_dir, "CellType_Markers_All.xlsx"))
    write_xlsx(top_markers,
               file.path(mode_b_dir, paste0("CellType_Top", MODE_B_TOP_N, "_Markers.xlsx")))
    message(paste("  Mode B complete.", nrow(all_markers), "significant marker genes."))
  }
}

# =============================================================================
# --- PART 6: MODE C — WITHIN-SUBTYPE WILCOXON --------------------------------
# =============================================================================
if (RUN_MODE_C_WITHIN_SUBTYPE) {
  message(paste("\n=== MODE C: Within-Subtype |", SUBTYPE_CELL_TYPE, "==="))
  if (!file.exists(SUBTYPE_RDS)) {
    warning(paste("[SKIP] Sub-cluster RDS not found:", SUBTYPE_RDS))
  } else {
    data_sub   <- readRDS(SUBTYPE_RDS)
    mode_c_dir <- file.path(DE_DIR, paste0("ModeC_", SUBTYPE_CELL_TYPE, "_", GROUP_1, "_vs_", GROUP_2))
    if (!dir.exists(mode_c_dir)) dir.create(mode_c_dir, recursive = TRUE)

    sub_types <- sort(as.character(unique(data_sub@meta.data[[SUBTYPE_COLUMN]])))
    all_sub   <- list()
    for (st in sub_types[!is.na(sub_types)]) {
      r <- run_de_wilcoxon(data_sub, st, SUBTYPE_COLUMN, CONDITION_COLUMN, GROUP_1, GROUP_2)
      if (is.null(r)) next
      all_sub[[st]] <- r
      save_de_results(r, st, "Wilcoxon", mode_c_dir)
    }
    if (length(all_sub) > 0)
      write_xlsx(all_sub, file.path(mode_c_dir, paste0("ModeC_ALL_SUBTYPES_", GROUP_1, "_vs_", GROUP_2, ".xlsx")))
    message("  Mode C complete.")
  }
}

# =============================================================================
# --- PART 7: MODE D — SplineDV + DE OVERLAP + fgsea ORA ---------------------
# =============================================================================
if (RUN_MODE_D_SPLINEDV) {
  message(paste("\n=== MODE D: SplineDV + DE Overlap + Enrichment |",
                GROUP_1, "vs", GROUP_2, "==="))

  mode_d_dir <- file.path(DE_DIR, paste0("ModeD_SplineDV_", GROUP_1, "_vs_", GROUP_2))
  if (!dir.exists(mode_d_dir)) dir.create(mode_d_dir, recursive = TRUE)

  if (!RUN_MODE_A_WILCOXON || length(wilcox_results_all) == 0) {
    stop("[ERROR] Mode D requires Mode A Wilcoxon results. Set RUN_MODE_A_WILCOXON = TRUE.")
  }

  dv_summary_rows <- list()

  for (ct in CELL_TYPES_TO_TEST) {
    safe_ct <- gsub("[^A-Za-z0-9_]", "_", ct)
    ct_dir  <- file.path(mode_d_dir, safe_ct)
    if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)
    message(paste("\n  --- Mode D:", ct, "---"))

    # --- Step D.1: Run SplineDV ---
    dv_result <- run_splinedv(data, ct, CELLTYPE_COLUMN, CONDITION_COLUMN, GROUP_1, GROUP_2)
    if (is.null(dv_result)) next

    write.csv(dv_result, file.path(ct_dir, paste0(safe_ct, "_SplineDV_full_results.csv")),
              row.names = FALSE)

    dv_sig_genes <- dv_result$gene[!is.na(dv_result$padj) & dv_result$padj < DV_PADJ_THRESH]
    message(paste("    DV significant genes (padj <", DV_PADJ_THRESH, "):", length(dv_sig_genes)))
    if (length(dv_sig_genes) == 0) {
      message("    No significant DV genes — skipping overlap.")
      next
    }

    # --- Step D.2: Get DE results for this cell type (from Mode A) ---
    if (!ct %in% names(wilcox_results_all)) {
      message(paste("    No Wilcoxon DE results for", ct, "— skipping overlap."))
      next
    }
    de_ct      <- wilcox_results_all[[ct]]
    de_universe <- de_ct$gene   # All DE-tested genes = ORA background

    up_genes   <- de_ct$gene[!is.na(de_ct$padj) & de_ct$padj < DE_PADJ_THRESH & de_ct$log2FC >  DE_LFC_THRESH]
    down_genes <- de_ct$gene[!is.na(de_ct$padj) & de_ct$padj < DE_PADJ_THRESH & de_ct$log2FC < -DE_LFC_THRESH]

    message(paste("    DE up:", length(up_genes), "| DE down:", length(down_genes),
                  "| DE universe:", length(de_universe)))

    # --- Step D.3: Compute overlaps ---
    overlap_up   <- intersect(dv_sig_genes, up_genes)
    overlap_down <- intersect(dv_sig_genes, down_genes)
    message(paste("    DV ∩ Up-DE:", length(overlap_up),
                  "| DV ∩ Down-DE:", length(overlap_down)))

    # Save overlap gene table with all statistics
    all_overlap_genes <- sort(unique(c(overlap_up, overlap_down)))
    if (length(all_overlap_genes) > 0) {
      overlap_df <- data.frame(
        gene         = all_overlap_genes,
        in_DV_sig    = all_overlap_genes %in% dv_sig_genes,
        in_Up_DE     = all_overlap_genes %in% up_genes,
        in_Down_DE   = all_overlap_genes %in% down_genes,
        DV_padj      = dv_result$padj[match(all_overlap_genes, dv_result$gene)],
        DE_log2FC    = de_ct$log2FC[match(all_overlap_genes, de_ct$gene)],
        DE_padj      = de_ct$padj[match(all_overlap_genes, de_ct$gene)],
        stringsAsFactors = FALSE
      )
      write.csv(overlap_df, file.path(ct_dir, paste0(safe_ct, "_DV_DE_overlap_genes.csv")),
                row.names = FALSE)
    }

    # Summary row for this cell type
    dv_summary_rows[[ct]] <- data.frame(
      CellType        = ct,
      N_DV_sig        = length(dv_sig_genes),
      N_DE_up         = length(up_genes),
      N_DE_down       = length(down_genes),
      N_Overlap_Up    = length(overlap_up),
      N_Overlap_Down  = length(overlap_down),
      stringsAsFactors = FALSE
    )

    # --- Step D.4: fgsea ORA on overlaps ---
    ora_sheets <- list()

    # Overlap Up: DV ∩ Up-DE
    if (length(overlap_up) >= ENRICHMENT_MIN_SIZE) {
      message(paste("    ORA: DV ∩ Up-DE (n =", length(overlap_up), ")..."))
      ora_up <- run_fgsea_ora(overlap_up, de_universe, msig_pathways)
      if (!is.null(ora_up) && nrow(ora_up) > 0) {
        ora_sheets[["DV_UpDE"]] <- ora_up
        write.csv(ora_up, file.path(ct_dir, paste0(safe_ct, "_ORA_DV_UpDE.csv")),
                  row.names = FALSE)
        p_up <- plot_ora_barplot(
          ora_up, fill_color = "#E41A1C",
          title = paste0(ct, " | DV ∩ Up-DE\n", GROUP_1, " vs ", GROUP_2)
        )
        if (!is.null(p_up))
          ggsave(file.path(ct_dir, paste0(safe_ct, "_ORA_DV_UpDE_barplot.png")),
                 p_up, width = 10, height = 7, dpi = DPI_SETTING, bg = "white")
        message(paste("    Sig. pathways (Up):", nrow(ora_up)))
      } else { message("    No significant pathways for DV ∩ Up-DE.") }
    } else {
      message(paste0("    [SKIP ORA Up] n=", length(overlap_up), " < min (", ENRICHMENT_MIN_SIZE, ")"))
    }

    # Overlap Down: DV ∩ Down-DE
    if (length(overlap_down) >= ENRICHMENT_MIN_SIZE) {
      message(paste("    ORA: DV ∩ Down-DE (n =", length(overlap_down), ")..."))
      ora_down <- run_fgsea_ora(overlap_down, de_universe, msig_pathways)
      if (!is.null(ora_down) && nrow(ora_down) > 0) {
        ora_sheets[["DV_DownDE"]] <- ora_down
        write.csv(ora_down, file.path(ct_dir, paste0(safe_ct, "_ORA_DV_DownDE.csv")),
                  row.names = FALSE)
        p_down <- plot_ora_barplot(
          ora_down, fill_color = "#377EB8",
          title = paste0(ct, " | DV ∩ Down-DE\n", GROUP_1, " vs ", GROUP_2)
        )
        if (!is.null(p_down))
          ggsave(file.path(ct_dir, paste0(safe_ct, "_ORA_DV_DownDE_barplot.png")),
                 p_down, width = 10, height = 7, dpi = DPI_SETTING, bg = "white")
        message(paste("    Sig. pathways (Down):", nrow(ora_down)))
      } else { message("    No significant pathways for DV ∩ Down-DE.") }
    } else {
      message(paste0("    [SKIP ORA Down] n=", length(overlap_down), " < min (", ENRICHMENT_MIN_SIZE, ")"))
    }

    if (length(ora_sheets) > 0)
      write_xlsx(ora_sheets, file.path(ct_dir, paste0(safe_ct, "_ORA_enrichment.xlsx")))
  }

  # Cross-cell-type summary
  if (length(dv_summary_rows) > 0) {
    dv_summary_df <- do.call(rbind, dv_summary_rows)
    write.csv(dv_summary_df,
              file.path(mode_d_dir, "SplineDV_DE_overlap_summary.csv"),
              row.names = FALSE)
    message("\n--- Mode D Summary ---")
    print(dv_summary_df)
  }
  message("  Mode D complete.")
}

# =============================================================================
# --- PART 8: FINAL MESSAGE ---------------------------------------------------
# =============================================================================
message(paste0(
  "\n=== DIFFERENTIAL EXPRESSION COMPLETE ===\n",
  "  Output directory: ", DE_DIR, "\n\n",
  if (RUN_MODE_A_WILCOXON) paste0(
    "  Mode A  (Wilcoxon):  ",
    file.path(DE_DIR, paste0("ModeA_Wilcoxon_", GROUP_1, "_vs_", GROUP_2)), "\n") else "",
  if (RUN_MODE_D_SPLINEDV) paste0(
    "  Mode D  (SplineDV):  ",
    file.path(DE_DIR, paste0("ModeD_SplineDV_", GROUP_1, "_vs_", GROUP_2)), "\n") else "",
  "\nNOTE: Wilcoxon p-values are exploratory (pseudoreplication). Use them\n",
  "for discovery; pair with SplineDV DV evidence for priority genes.\n",
  "\nNext step: Run '04_cell_scoring.R' for pathway-level AUCell scoring.\n"
))
