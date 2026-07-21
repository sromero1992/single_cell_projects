# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 5: DE + TWO-WAY ANOVA ANALYSIS
# Version: 3.0 (SplineDV + Wilcoxon DE + Two-Way ANOVA; no pseudo-bulk)
#
# PURPOSE:
#   Performs differential gene expression (DGE), differential variability
#   (DV), pathway enrichment, and two-way ANOVA analysis between user-defined
#   groups across annotated cell populations.
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
#             Scripts 03/04 and runs Mode A at the sub-type level.
#
#   MODE D  — SplineDV + DE Overlap + fgsea ORA Enrichment: For each cell
#             type, runs SplineDV to find differentially variable (DV) genes,
#             then overlaps SIGNIFICANT DV genes (pval < DV_PVAL_THRESH) with
#             padj-filtered DE genes from Mode A (all directions pooled), then
#             splits by log2FC sign to get directional subsets, and runs fgsea
#             over-representation analysis (fora) on each overlap gene set.
#
#             KEY CHANGE (v3.0): DV genes filtered by raw pval < DV_PVAL_THRESH.
#             DE genes filtered by padj < DE_PADJ_THRESH (all directions first).
#             Overlap is computed undirected, then split by log2FC sign.
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
#   MODE E  — Two-Way ANOVA (Genotype × Diet or any two factors): Runs a
#             parallel, memory-efficient two-way ANOVA on all genes for a
#             user-specified cell type. Genes are ranked by the interaction
#             term p-value (Genotype × Diet). Top 250 interaction genes are
#             sent to Enrichr. Produces a sorted results CSV + a violin plot
#             for the top interaction gene.
#             Requires: parallel, broom, enrichR
#
# DEPENDENCIES:
#   Core:    Seurat, dplyr, ggplot2, ggrepel, writexl, patchwork
#   Mode D:  SplineDV (GitHub: Xenon8778/SplineDV), fgsea, msigdbr, stringr
#   Mode E:  parallel, broom, enrichR
#
# NEXT STEP: Script 06 — Cell Scoring & Pathway Enrichment (06_cell_scoring.R).
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
RUN_MODE_E_TWOWAY_ANOVA      <- FALSE  # Two-way ANOVA (Genotype x Diet) per cell type

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
# GROUP_2: The "reference" / control group. Must match values in CONDITION_COLUMN.
#
# CONVENTION — DE direction vs. plot order:
#   GROUP_1 = KO  → positive log2FC = upregulated in KO (keep this order for DE).
#   For plot axis ordering (WT on the left), set CONDITION_LEVELS = c("WT", ...)
#   in Script 04. These two conventions are intentionally independent.
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

# --- 1.7: Mode D — SplineDV Parameters --------------------------------------

# DV_PVAL_THRESH: Raw (unadjusted) p-value threshold for calling a gene
#   "differentially variable". SplineDV detects changes in variance (not mean
#   expression). We use raw pval here because the overlap is then further
#   restricted to genes that also pass padj filtering on the DE side — the
#   combined requirement keeps the joint list stringent without double-correcting.
DV_PVAL_THRESH <- 0.05

# --- 1.8: GSEA Parameters (ranked gene-set enrichment — ALL tested DE genes) -
#
# GSEA uses the FULL ranked list of tested genes (not just significant ones).
# Genes are ranked by: sign(log2FC) × -log10(padj) — a signed significance
# score that emphasises genes that are both strongly directional AND significant.
# This is more statistically robust than ORA because it does not require an
# arbitrary significance cutoff to define the gene list.
#
# Gene sets come from MSigDB (msigdbr) — same species and category settings
# apply to both GSEA and are loaded once in Part 3.
GSEA_SPECIES     <- "Mus musculus"  # or "Homo sapiens"
GSEA_CATEGORY    <- "C5"            # C5 = GO, H = Hallmark, C2 = KEGG/Reactome
GSEA_SUBCATEGORY <- "BP"            # BP = Biological Process; "" for Hallmark
GSEA_MIN_SIZE    <- 15              # Min genes in pathway to test (recommend 15)
GSEA_MAX_SIZE    <- 500             # Max genes in pathway to test
GSEA_NPERM       <- 1000            # Permutations (increase to 10000 for publication)
GSEA_PADJ_THRESH <- 0.05            # FDR threshold for significant GSEA results
GSEA_TOP_N_PLOT  <- 15              # Top N pathways shown in the NES bar plot

# --- 1.9: Enrichr ORA Parameters (on specific ranked gene lists) -------------
#
# Enrichr ORA is run on five gene lists per cell type:
#   1. DE up    — top ENRICHR_TOP_N_GENES up-regulated genes (by padj then log2FC)
#   2. DE down  — top ENRICHR_TOP_N_GENES down-regulated genes
#   3. DV       — top ENRICHR_TOP_N_GENES differentially variable genes (by padj)
#   4. DV ∩ DE up   — intersection of lists 3 and 1 (no directionality on DV)
#   5. DV ∩ DE down — intersection of lists 3 and 2
#
# Three databases are queried for each list. Verify exact database names with:
#   enrichR::listEnrichrDbs()   (requires internet connection)
#
# NOTE on mouse data: Enrichr accepts mouse gene symbols (title-case, e.g. Foxp3)
#   alongside human. Most GO and KEGG databases map both conventions. If you see
#   very low overlap counts, check symbol capitalisation or use mouse-specific
#   databases (e.g. "GO_Biological_Process_2023_Mouse").
ENRICHR_DATABASES   <- c(
  "GO_Biological_Process_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2026"
)
ENRICHR_TOP_N_GENES  <- 250   # Top N genes taken from each ranked DE / DV list
ENRICHR_TOP_N_PLOT   <- 15    # Top N significant terms shown per bar plot
ENRICHR_PADJ_THRESH  <- 0.05  # Adjusted p-value cutoff for Enrichr terms

# --- 1.10: Mode E — Two-Way ANOVA Parameters ---------------------------------
#
# ANOVA_CELL_TYPE: The cell type to analyse (subset from CELLTYPE_COLUMN).
#   Must exactly match a value in CELLTYPE_COLUMN.
#   Examples: "Stem cells", "Goblet cells", "T cells"
ANOVA_CELL_TYPE    <- "Stem cells"

# ANOVA_FACTOR1 / ANOVA_FACTOR2: Metadata column names for the two ANOVA factors.
#   These columns must already exist in the Seurat metadata (Script 01 imports).
ANOVA_FACTOR1      <- "Genotype"   # e.g., WT vs KO
ANOVA_FACTOR2      <- "Diet"       # e.g., cellulose vs inulin

# ANOVA_INTERACTION: If TRUE, tests Genotype * Diet interaction term.
#   Set FALSE for an additive model (Genotype + Diet).
ANOVA_INTERACTION  <- TRUE

# ANOVA_N_CORES: Number of CPU cores for parallel ANOVA. Set to detectCores()-1
#   for maximum throughput on a multi-core server.
ANOVA_N_CORES      <- 4

# ANOVA_CHUNK_SIZE: Genes processed per parallel chunk. Larger = fewer chunks
#   but higher peak RAM. 500 is a safe default for most datasets.
ANOVA_CHUNK_SIZE   <- 500

# ANOVA_TOP_N_ENRICHR: Number of top interaction genes sent to Enrichr.
ANOVA_TOP_N_ENRICHR <- 250

# ANOVA_MIN_CELLS: Skip cell type if either factor level has fewer cells.
ANOVA_MIN_CELLS    <- 50

# --- 1.11: Plot Settings -----------------------------------------------------
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
  # Normalize raw p-value column name across SplineDV versions
  if (!"pval" %in% colnames(dv_result) && "pvalue" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, pval = pvalue)
  if (!"pval" %in% colnames(dv_result) && "p.value" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, pval = p.value)
  if (!"pval" %in% colnames(dv_result) && "PValue" %in% colnames(dv_result))
    dv_result <- dplyr::rename(dv_result, pval = PValue)

  dv_result$cell_type  <- cell_type
  dv_result$comparison <- paste0(group1, "_vs_", group2)
  dv_result$n_cells_g1 <- n1
  dv_result$n_cells_g2 <- n2
  return(dv_result %>% arrange(padj))
}

# ---------------------------------------------------------------------------
#' run_gsea_ranked
#'
#'   Standard GSEA via fgsea::fgsea() on ALL tested DE genes.
#'   Unlike ORA, no significance cutoff is applied — every gene tested in
#'   Mode A enters the ranked list. This makes GSEA more sensitive to modest
#'   but coordinated expression shifts that would be missed by ORA.
#'
#'   Ranking metric: sign(log2FC) × -log10(padj + 1e-300)
#'     → up-regulated significant genes score highest (large positive)
#'     → down-regulated significant genes score lowest (large negative)
#'     → non-significant genes cluster near zero
#'
#'   Gene sets come from MSigDB (pre-loaded via msigdbr in Part 3).
# ---------------------------------------------------------------------------
run_gsea_ranked <- function(de_result, pathways, label = "") {
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("[ERROR] fgsea not installed. Run 00_rlibs_installation.R.")
  library(fgsea)

  ranked_df <- de_result %>%
    dplyr::filter(!is.na(log2FC), !is.na(padj)) %>%
    dplyr::mutate(rank_score = sign(log2FC) * -log10(padj + 1e-300)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(desc(rank_score))

  ranked <- setNames(ranked_df$rank_score, ranked_df$gene)
  if (length(ranked) < GSEA_MIN_SIZE * 2) {
    message(paste("  [SKIP GSEA]", label, "— too few ranked genes:", length(ranked)))
    return(NULL)
  }
  message(paste("  [GSEA]", label, "| ranked genes:", length(ranked)))

  result <- tryCatch(
    fgsea::fgsea(pathways    = pathways,
                 stats       = ranked,
                 minSize     = GSEA_MIN_SIZE,
                 maxSize     = GSEA_MAX_SIZE,
                 nPermSimple = GSEA_NPERM),
    error = function(e) { message("  [ERROR GSEA]: ", e$message); return(NULL) }
  )
  if (is.null(result) || nrow(result) == 0) return(NULL)
  result$leadingEdge <- sapply(result$leadingEdge, paste, collapse = "|")
  result %>% dplyr::arrange(padj)
}

# ---------------------------------------------------------------------------
#' plot_gsea_barplot
#'   Horizontal NES bar chart. Positive NES = enriched in GROUP_1 up direction.
#'   Top N up- and top N down-enriched pathways are shown (balanced).
# ---------------------------------------------------------------------------
plot_gsea_barplot <- function(gsea_result, title,
                               top_n       = GSEA_TOP_N_PLOT,
                               padj_thresh = GSEA_PADJ_THRESH) {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) return(NULL)
  df_sig <- gsea_result %>% dplyr::filter(padj < padj_thresh)
  if (nrow(df_sig) == 0) return(NULL)

  df_up   <- df_sig %>% dplyr::filter(NES > 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = ceiling(top_n / 2))
  df_down <- df_sig %>% dplyr::filter(NES < 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = floor(top_n / 2))
  df_plot <- dplyr::bind_rows(df_up, df_down) %>%
    dplyr::mutate(
      pathway_label = gsub("^GOBP_|^GOMF_|^HALLMARK_|^KEGG_|^REACTOME_", "", pathway) %>%
        gsub("_", " ", .) %>% tools::toTitleCase() %>% stringr::str_wrap(width = 50),
      Direction = ifelse(NES > 0, paste0("Up in ", GROUP_1), paste0("Down in ", GROUP_1))
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway_label = factor(pathway_label, levels = pathway_label))

  if (nrow(df_plot) == 0) return(NULL)

  dir_colors <- setNames(c("#E41A1C", "#377EB8"),
                         c(paste0("Up in ", GROUP_1), paste0("Down in ", GROUP_1)))

  ggplot(df_plot, aes(x = NES, y = pathway_label, fill = Direction)) +
    geom_col(alpha = 0.85, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = dir_colors) +
    labs(title    = title,
         subtitle = paste0(nrow(df_sig), " sig. pathways (padj < ", padj_thresh, ") | top ", top_n, " shown"),
         x        = "Normalized Enrichment Score (NES)",
         y        = NULL, fill = NULL) +
    theme_classic() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
          axis.text.y   = element_text(size = 9),
          legend.position = "bottom")
}

# ---------------------------------------------------------------------------
#' save_enrichr_results
#'
#'   Run Enrichr ORA on a gene list against ENRICHR_DATABASES, save:
#'     - One multi-sheet Excel (one sheet per database, significant terms only)
#'     - One horizontal bar plot PNG per database
#'
#'   Returns the raw enrichR result list invisibly (one element per database).
#'
#' @param gene_list  Character vector of gene symbols to test.
#' @param label      Short identifier used in file names and plot titles.
#' @param out_dir    Directory to write outputs into.
#' @param fill_color Bar fill color for the bar plots.
# ---------------------------------------------------------------------------
save_enrichr_results <- function(gene_list, label, out_dir,
                                  fill_color = "#E41A1C") {
  if (!requireNamespace("enrichR", quietly = TRUE))
    stop("[ERROR] enrichR not installed. Add 'enrichR' to cran_pkgs in Script 00 and re-run.")
  library(enrichR)

  if (length(gene_list) == 0) {
    message(paste("  [SKIP Enrichr]", label, "— empty gene list"))
    return(invisible(NULL))
  }
  message(paste("  [Enrichr]", label, "| genes:", length(gene_list),
                "| databases:", paste(ENRICHR_DATABASES, collapse = ", ")))

  enrichr_res <- tryCatch({
    enrichR::setEnrichrSite("Enrichr")
    enrichR::enrichr(gene_list, ENRICHR_DATABASES)
  }, error = function(e) {
    message(paste("  [ERROR Enrichr]", label, ":", e$message))
    return(NULL)
  })
  if (is.null(enrichr_res)) return(invisible(NULL))

  # ---- Save multi-sheet Excel (one sheet per database, sig terms only) -----
  xl_sheets <- lapply(ENRICHR_DATABASES, function(db) {
    df <- enrichr_res[[db]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df_sig <- df %>% dplyr::filter(Adjusted.P.value < ENRICHR_PADJ_THRESH) %>%
      dplyr::arrange(Adjusted.P.value)
    if (nrow(df_sig) == 0) return(NULL)
    df_sig
  })
  names(xl_sheets) <- ENRICHR_DATABASES
  xl_sheets <- xl_sheets[!sapply(xl_sheets, is.null)]
  if (length(xl_sheets) > 0)
    write_xlsx(xl_sheets,
               file.path(out_dir, paste0(label, "_Enrichr.xlsx")))

  # ---- Bar plot per database -----------------------------------------------
  db_colors <- c("#E41A1C", "#377EB8", "#4DAF4A")  # red, blue, green per DB
  for (k in seq_along(ENRICHR_DATABASES)) {
    db     <- ENRICHR_DATABASES[k]
    df     <- enrichr_res[[db]]
    if (is.null(df) || nrow(df) == 0) next
    df_sig <- df %>% dplyr::filter(Adjusted.P.value < ENRICHR_PADJ_THRESH)
    if (nrow(df_sig) == 0) {
      message(paste("    No sig. Enrichr terms for", label, "/", db))
      next
    }
    df_plot <- df_sig %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::slice_head(n = ENRICHR_TOP_N_PLOT) %>%
      dplyr::mutate(
        term_label     = gsub("\\s*\\(GO:\\d+\\)$", "", Term) %>%
          gsub("_", " ", .) %>%
          stringr::str_wrap(width = 55),
        neg_log10_padj = -log10(Adjusted.P.value + 1e-300)
      ) %>%
      dplyr::arrange(neg_log10_padj) %>%
      dplyr::mutate(term_label = factor(term_label, levels = term_label))

    p <- ggplot(df_plot, aes(x = neg_log10_padj, y = term_label)) +
      geom_col(fill = db_colors[k], alpha = 0.85, color = "black", linewidth = 0.3) +
      geom_vline(xintercept = -log10(ENRICHR_PADJ_THRESH),
                 linetype = "dashed", color = "grey40", linewidth = 0.7) +
      geom_text(aes(label = Overlap), hjust = -0.1, size = 3.2) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(title    = paste0(label),
           subtitle = paste0(db, " | ", nrow(df_sig),
                             " sig. terms | top ", ENRICHR_TOP_N_PLOT, " shown"),
           x        = expression(-log[10] * "(adj. p-value)"),
           y        = NULL) +
      theme_classic() +
      theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey50"),
            axis.text.y   = element_text(size = 9),
            axis.text.x   = element_text(size = 10))

    db_safe <- gsub("[^A-Za-z0-9_]", "_", db)
    plot_h  <- max(4, nrow(df_plot) * 0.38 + 2)
    ggsave(file.path(out_dir, paste0(label, "_Enrichr_", db_safe, ".png")),
           p, width = 11, height = plot_h, dpi = DPI_SETTING, bg = "white")
    message(paste("    Sig. terms (", db, "):", nrow(df_sig)))
  }
  invisible(enrichr_res)
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
  msig_df <- msigdbr(species = GSEA_SPECIES,
                     category = GSEA_CATEGORY,
                     subcategory = if (nchar(GSEA_SUBCATEGORY) > 0) GSEA_SUBCATEGORY else NULL)
  msig_pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
  message(paste("  Loaded", length(msig_pathways), "pathways from",
                GSEA_CATEGORY, "/", GSEA_SUBCATEGORY))
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
# --- PART 7: MODE D — SplineDV + GSEA + Enrichr ORA -------------------------
#
#  Per-cell-type workflow:
#    D.1  SplineDV        → full DV results CSV + top-N DV gene list
#    D.2  DE gene lists   → top-N up and top-N down lists (no padj cutoff,
#                           just rank — so Enrichr always gets a meaningful list)
#    D.3  Gene-list CSVs  → saved individually for all 5 lists + overlap tables
#                           with both DE and DV statistics per gene
#    D.4  GSEA            → fgsea::fgsea on ALL tested DE genes (ranked),
#                           MSigDB gene sets (GSEA_SPECIES / GSEA_CATEGORY /
#                           GSEA_SUBCATEGORY configured in Section 1.8)
#    D.5  Enrichr ORA     → enrichR::enrichr on 5 gene lists against
#                           ENRICHR_DATABASES (configured in Section 1.9):
#                             1. DE up top ENRICHR_TOP_N_GENES
#                             2. DE down top ENRICHR_TOP_N_GENES
#                             3. DV top ENRICHR_TOP_N_GENES (no direction)
#                             4. DV ∩ DE up (overlap)
#                             5. DV ∩ DE down (overlap)
#    D.6  Summary CSV     → cross-cell-type table of list sizes
#
#  Key design decisions:
#    • "Top N" uses NO significance cutoff — genes are ranked by padj (DV) or
#      by padj + direction (DE). This guarantees a non-empty query list even for
#      cell types with few significant hits, making Enrichr results comparable
#      across cell types. Change ENRICHR_TOP_N_GENES (Section 1.9) to adjust.
#    • SplineDV has NO directionality — DV genes are simply "more variable in
#      one group". Overlaps with DE up / DE down capture genes whose variance
#      AND mean both shift.
#    • GSEA uses ALL genes tested in Mode A — no cutoff applied, ranked by
#      sign(log2FC) × -log10(padj + 1e-300). Pathways loaded from MSigDB via
#      msigdbr (pre-loaded in Part 3).
# =============================================================================
if (RUN_MODE_D_SPLINEDV) {
  message(paste("\n=== MODE D: SplineDV + GSEA + Enrichr ORA |",
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

    # -------------------------------------------------------------------------
    # Step D.1: SplineDV — differential variability
    # -------------------------------------------------------------------------
    dv_result <- run_splinedv(data, ct, CELLTYPE_COLUMN, CONDITION_COLUMN, GROUP_1, GROUP_2)
    if (is.null(dv_result)) next

    write.csv(dv_result,
              file.path(ct_dir, paste0(safe_ct, "_SplineDV_full_results.csv")),
              row.names = FALSE)

    # Top-N DV genes by padj — no significance cutoff (used for Enrichr DV list)
    dv_top_n <- dv_result %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = ENRICHR_TOP_N_GENES) %>%
      dplyr::pull(gene)

    # DV significant genes: raw pval < DV_PVAL_THRESH (unadjusted).
    # We intentionally use the raw p-value here. The overlap is further
    # restricted to genes that also pass padj filtering on the DE side,
    # so the joint gene set is stringent without double-correcting.
    dv_sig_genes <- dv_result %>%
      dplyr::filter(!is.na(pval), pval < DV_PVAL_THRESH) %>%
      dplyr::pull(gene)

    n_dv_sig <- length(dv_sig_genes)
    message(paste("    DV sig (pval <", DV_PVAL_THRESH, "):", n_dv_sig,
                  "| top", length(dv_top_n), "taken for Enrichr DV list"))

    # -------------------------------------------------------------------------
    # Step D.2: DE gene lists for this cell type (from Mode A)
    # -------------------------------------------------------------------------
    if (!ct %in% names(wilcox_results_all)) {
      message(paste("    No Wilcoxon DE results for", ct, "— skipping."))
      next
    }
    de_ct <- wilcox_results_all[[ct]]

    # All DE significant genes (padj < DE_PADJ_THRESH, any direction).
    # Used for the undirected DV/DE overlap; direction is assigned afterward.
    de_sig_all <- de_ct %>%
      dplyr::filter(!is.na(padj), !is.na(log2FC), padj < DE_PADJ_THRESH) %>%
      dplyr::pull(gene)

    # Top-N up-regulated: padj < DE_PADJ_THRESH AND log2FC > 0,
    #   sorted by padj asc then log2FC desc — used for Enrichr DE up list
    de_up_n <- de_ct %>%
      dplyr::filter(!is.na(padj), !is.na(log2FC),
                    padj < DE_PADJ_THRESH, log2FC > 0) %>%
      dplyr::arrange(padj, desc(log2FC)) %>%
      dplyr::slice_head(n = ENRICHR_TOP_N_GENES) %>%
      dplyr::pull(gene)

    # Top-N down-regulated: padj < DE_PADJ_THRESH AND log2FC < 0,
    #   sorted by padj asc then log2FC asc — used for Enrichr DE down list
    de_dn_n <- de_ct %>%
      dplyr::filter(!is.na(padj), !is.na(log2FC),
                    padj < DE_PADJ_THRESH, log2FC < 0) %>%
      dplyr::arrange(padj, log2FC) %>%
      dplyr::slice_head(n = ENRICHR_TOP_N_GENES) %>%
      dplyr::pull(gene)

    # Overlap: intersect DV sig (pval) with ALL DE sig (padj), direction-
    # agnostic first. Then split by log2FC sign to get directional subsets.
    overlap_all <- intersect(de_sig_all, dv_sig_genes)
    overlap_up  <- overlap_all[
      overlap_all %in% de_ct$gene[!is.na(de_ct$log2FC) & de_ct$log2FC > 0]
    ]
    overlap_dn  <- overlap_all[
      overlap_all %in% de_ct$gene[!is.na(de_ct$log2FC) & de_ct$log2FC < 0]
    ]

    message(paste("    DE sig (padj<", DE_PADJ_THRESH, "):", length(de_sig_all),
                  "| DV sig (pval<", DV_PVAL_THRESH, "):", n_dv_sig))
    message(paste("    Overlap (undirected):", length(overlap_all),
                  "| -> up:", length(overlap_up),
                  "| -> dn:", length(overlap_dn)))


    # -------------------------------------------------------------------------
    # Step D.3: Save all 5 gene lists as individual CSVs
    # -------------------------------------------------------------------------
    enrichr_pfx <- paste0(safe_ct, "_", GROUP_1, "_vs_", GROUP_2)

    write.csv(data.frame(gene = de_up_n),
              file.path(ct_dir, paste0(safe_ct, "_genes_DE_up_top", ENRICHR_TOP_N_GENES, ".csv")),
              row.names = FALSE)
    write.csv(data.frame(gene = de_dn_n),
              file.path(ct_dir, paste0(safe_ct, "_genes_DE_dn_top", ENRICHR_TOP_N_GENES, ".csv")),
              row.names = FALSE)
    write.csv(data.frame(gene = dv_top_n),
              file.path(ct_dir, paste0(safe_ct, "_genes_DV_top", ENRICHR_TOP_N_GENES, ".csv")),
              row.names = FALSE)

    # Overlap tables include DE + DV statistics for each gene
    if (length(overlap_up) > 0) {
      write.csv(
        data.frame(
          gene      = overlap_up,
          DE_log2FC = de_ct$log2FC[match(overlap_up, de_ct$gene)],
          DE_padj   = de_ct$padj[match(overlap_up, de_ct$gene)],
          DV_padj   = dv_result$padj[match(overlap_up, dv_result$gene)],
          stringsAsFactors = FALSE
        ),
        file.path(ct_dir, paste0(safe_ct, "_genes_overlap_DV_DEup.csv")),
        row.names = FALSE
      )
    }
    if (length(overlap_dn) > 0) {
      write.csv(
        data.frame(
          gene      = overlap_dn,
          DE_log2FC = de_ct$log2FC[match(overlap_dn, de_ct$gene)],
          DE_padj   = de_ct$padj[match(overlap_dn, de_ct$gene)],
          DV_padj   = dv_result$padj[match(overlap_dn, dv_result$gene)],
          stringsAsFactors = FALSE
        ),
        file.path(ct_dir, paste0(safe_ct, "_genes_overlap_DV_DEdn.csv")),
        row.names = FALSE
      )
    }

    # -------------------------------------------------------------------------
    # Step D.4: GSEA on ALL ranked DE genes (fgsea::fgsea)
    # -------------------------------------------------------------------------
    if (!is.null(msig_pathways)) {
      gsea_result <- run_gsea_ranked(de_ct, msig_pathways, label = ct)
      if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
        write.csv(gsea_result,
                  file.path(ct_dir, paste0(safe_ct, "_GSEA_results.csv")),
                  row.names = FALSE)
        p_gsea <- plot_gsea_barplot(
          gsea_result,
          title       = paste0(ct, " | GSEA (all DE genes)\n", GROUP_1, " vs ", GROUP_2),
          top_n       = GSEA_TOP_N_PLOT,
          padj_thresh = GSEA_PADJ_THRESH
        )
        if (!is.null(p_gsea))
          ggsave(file.path(ct_dir, paste0(safe_ct, "_GSEA_NES_barplot.png")),
                 p_gsea, width = 11, height = 8, dpi = DPI_SETTING, bg = "white")
        message(paste("    GSEA sig. pathways (padj <", GSEA_PADJ_THRESH, "):",
                      sum(gsea_result$padj < GSEA_PADJ_THRESH, na.rm = TRUE)))
      } else {
        message(paste("    No significant GSEA results for", ct))
      }
    }

    # -------------------------------------------------------------------------
    # Step D.5: Enrichr ORA on 5 gene lists
    #   Colors: up=red, dn=blue, DV=purple, overlap-up=orange, overlap-dn=brown
    # -------------------------------------------------------------------------

    # List 1 — DE up top-N
    if (length(de_up_n) > 0)
      save_enrichr_results(
        gene_list  = de_up_n,
        label      = paste0(enrichr_pfx, "_DE_up_top", ENRICHR_TOP_N_GENES),
        out_dir    = ct_dir,
        fill_color = "#E41A1C"   # red
      )

    # List 2 — DE down top-N
    if (length(de_dn_n) > 0)
      save_enrichr_results(
        gene_list  = de_dn_n,
        label      = paste0(enrichr_pfx, "_DE_dn_top", ENRICHR_TOP_N_GENES),
        out_dir    = ct_dir,
        fill_color = "#377EB8"   # blue
      )

    # List 3 — DV top-N (no directionality)
    if (length(dv_top_n) > 0)
      save_enrichr_results(
        gene_list  = dv_top_n,
        label      = paste0(enrichr_pfx, "_DV_top", ENRICHR_TOP_N_GENES),
        out_dir    = ct_dir,
        fill_color = "#984EA3"   # purple
      )

    # List 4 — DV (sig padj) ∩ DE up
    if (length(overlap_up) > 0)
      save_enrichr_results(
        gene_list  = overlap_up,
        label      = paste0(enrichr_pfx, "_overlap_DV_sig_DEup"),
        out_dir    = ct_dir,
        fill_color = "#FF7F00"   # orange
      )

    # List 5 — DV (sig padj) ∩ DE down
    if (length(overlap_dn) > 0)
      save_enrichr_results(
        gene_list  = overlap_dn,
        label      = paste0(enrichr_pfx, "_overlap_DV_sig_DEdn"),
        out_dir    = ct_dir,
        fill_color = "#A65628"   # brown
      )

    # -------------------------------------------------------------------------
    # Step D.6: Cross-cell-type summary row
    # -------------------------------------------------------------------------
    dv_summary_rows[[ct]] <- data.frame(
      CellType             = ct,
      N_DV_sig_padj05      = n_dv_sig,
      N_DV_top_N           = length(dv_top_n),
      N_DE_up_top_N        = length(de_up_n),
      N_DE_dn_top_N        = length(de_dn_n),
      N_Overlap_DV_DEup    = length(overlap_up),
      N_Overlap_DV_DEdn    = length(overlap_dn),
      stringsAsFactors     = FALSE
    )
  }

  # Cross-cell-type summary table
  if (length(dv_summary_rows) > 0) {
    dv_summary_df <- do.call(rbind, dv_summary_rows)
    write.csv(dv_summary_df,
              file.path(mode_d_dir, "SplineDV_DE_enrichment_summary.csv"),
              row.names = FALSE)
    message("\n--- Mode D Summary ---")
    print(dv_summary_df)
  }
  message("  Mode D complete.")
}

# =============================================================================
# --- PART 9: MODE E — TWO-WAY ANOVA (Genotype × Diet) -----------------------
#
#  Memory-efficient parallel two-way ANOVA over all genes for a single
#  user-specified cell type.
#
#  Workflow:
#    E.1  Identify cells for analysis WITHOUT subsetting the Seurat object.
#         This keeps RAM low — only the expression slice is materialised.
#    E.2  Build model formula: expression ~ Factor1 * Factor2 (or additive).
#    E.3  Run parallel chunked ANOVA using base R's parLapply.
#    E.4  Adjust p-values (BH), sort by interaction term, save CSV.
#    E.5  Enrichr ORA on top ANOVA_TOP_N_ENRICHR interaction genes.
#    E.6  Violin plot for top interaction gene.
#
#  NOTE on factors: the metadata columns ANOVA_FACTOR1 and ANOVA_FACTOR2 must
#  already exist in the Seurat object. Script 01 imports them from your Excel
#  metadata file. Use exactly the column names from that file.
# =============================================================================
if (RUN_MODE_E_TWOWAY_ANOVA) {
  message(paste("\n=== MODE E: Two-Way ANOVA |",
                ANOVA_FACTOR1, "x", ANOVA_FACTOR2,
                "| Cell type:", ANOVA_CELL_TYPE, "==="))
  library(parallel)
  library(broom)

  mode_e_dir <- file.path(DE_DIR, paste0("ModeE_ANOVA_", ANOVA_CELL_TYPE))
  if (!dir.exists(mode_e_dir)) dir.create(mode_e_dir, recursive = TRUE)

  # --- E.1: Identify cells ---------------------------------------------------
  if (!all(c(ANOVA_FACTOR1, ANOVA_FACTOR2) %in% colnames(data@meta.data))) {
    warning(paste("[SKIP Mode E] Factor columns not found in metadata:",
                  ANOVA_FACTOR1, "/", ANOVA_FACTOR2,
                  "\nVerify these columns exist in your metadata Excel file."))
  } else {

    # Ensure factors are properly typed
    data@meta.data[[ANOVA_FACTOR1]] <- as.factor(data@meta.data[[ANOVA_FACTOR1]])
    data@meta.data[[ANOVA_FACTOR2]] <- as.factor(data@meta.data[[ANOVA_FACTOR2]])

    Idents(data) <- CELLTYPE_COLUMN
    cells_anova  <- WhichCells(data, idents = ANOVA_CELL_TYPE)
    meta_anova   <- data@meta.data[cells_anova, c(ANOVA_FACTOR1, ANOVA_FACTOR2)]
    genes_anova  <- rownames(data)

    message(paste("  Cells for ANOVA:", length(cells_anova),
                  "| Genes:", length(genes_anova)))

    # Check minimum cell count per group
    group_table <- table(meta_anova[[ANOVA_FACTOR1]], meta_anova[[ANOVA_FACTOR2]])
    print(group_table)
    if (any(group_table < ANOVA_MIN_CELLS)) {
      warning(paste("[WARN Mode E] Some factor level combinations have <",
                    ANOVA_MIN_CELLS, "cells. Results for sparse groups may be unreliable."))
    }

    # --- E.2: Model formula --------------------------------------------------
    if (ANOVA_INTERACTION) {
      anova_formula <- as.formula(paste("expression ~",
                                        ANOVA_FACTOR1, "*", ANOVA_FACTOR2))
    } else {
      anova_formula <- as.formula(paste("expression ~",
                                        ANOVA_FACTOR1, "+", ANOVA_FACTOR2))
    }
    message(paste("  ANOVA model formula:", deparse(anova_formula)))

    # Interaction term name (used for sorting and column naming)
    interaction_term <- paste0(ANOVA_FACTOR1, ":", ANOVA_FACTOR2)

    # --- E.3: Parallel chunked ANOVA -----------------------------------------
    gene_chunks  <- split(genes_anova,
                          ceiling(seq_along(genes_anova) / ANOVA_CHUNK_SIZE))
    num_chunks   <- length(gene_chunks)
    message(paste("  Split", length(genes_anova), "genes into",
                  num_chunks, "chunks of", ANOVA_CHUNK_SIZE))

    cl <- makeCluster(ANOVA_N_CORES)
    clusterExport(cl, varlist = c("anova_formula", "ANOVA_FACTOR1", "ANOVA_FACTOR2"))
    clusterEvalQ(cl, { library(broom); library(dplyr) })

    all_chunk_results <- list()
    start_time <- Sys.time()
    message(paste("  Starting parallel ANOVA on", ANOVA_N_CORES, "cores..."))

    for (i in seq_len(num_chunks)) {
      if (i %% 10 == 0 || i == 1 || i == num_chunks)
        message(paste0("    Chunk ", i, " / ", num_chunks, " ..."))

      chunk_genes <- gene_chunks[[i]]
      chunk_data  <- FetchData(data,
                                vars = c(chunk_genes, ANOVA_FACTOR1, ANOVA_FACTOR2),
                                cells = cells_anova)

      chunk_res <- parLapply(cl, chunk_genes, function(gene, cdata) {
        gdf <- data.frame(
          expression = cdata[[gene]],
          Genotype   = cdata[[ANOVA_FACTOR1]],
          Diet       = cdata[[ANOVA_FACTOR2]]
        )
        # Use generic column names so the formula resolves inside the worker
        colnames(gdf)[2:3] <- c(ANOVA_FACTOR1, ANOVA_FACTOR2)
        aov_m <- tryCatch(aov(anova_formula, data = gdf), error = function(e) NULL)
        if (is.null(aov_m)) return(NULL)
        tidy_r <- broom::tidy(aov_m)
        tidy_r$gene <- gene
        tidy_r
      }, cdata = chunk_data)

      all_chunk_results[[i]] <- chunk_res
    }
    stopCluster(cl)

    end_time <- Sys.time()
    message(paste("  ANOVA completed in:",
                  round(difftime(end_time, start_time, units = "mins"), 2), "mins"))

    # --- E.4: Adjust p-values, format, sort ----------------------------------
    anova_df <- bind_rows(unlist(all_chunk_results, recursive = FALSE))
    interaction_p_adj_col <- paste0(interaction_term, ".p.adj")

    anova_wide <- anova_df %>%
      dplyr::filter(term != "Residuals") %>%
      dplyr::group_by(term) %>%
      dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, term, statistic, p.value, p.adj) %>%
      tidyr::pivot_wider(
        names_from  = term,
        values_from = c(statistic, p.value, p.adj),
        names_glue  = "{term}.{.value}"
      ) %>%
      dplyr::arrange(!!sym(interaction_p_adj_col))

    anova_csv <- file.path(mode_e_dir,
                            paste0("ANOVA_results_", gsub("[^A-Za-z0-9_]", "_", ANOVA_CELL_TYPE), ".csv"))
    write.csv(anova_wide, anova_csv, row.names = FALSE)
    message(paste("  ANOVA results saved:", anova_csv))
    message("  Top 10 genes by interaction term:")
    print(head(anova_wide, 10))

    # --- E.5: Enrichr on top ANOVA_TOP_N_ENRICHR interaction genes -----------
    top_anova_genes <- head(anova_wide$gene, ANOVA_TOP_N_ENRICHR)
    top_anova_genes <- top_anova_genes[!is.na(top_anova_genes)]

    if (length(top_anova_genes) > 0 && requireNamespace("enrichR", quietly = TRUE)) {
      safe_ct_anova <- gsub("[^A-Za-z0-9_]", "_", ANOVA_CELL_TYPE)
      save_enrichr_results(
        gene_list  = top_anova_genes,
        label      = paste0(safe_ct_anova, "_ANOVA_top", ANOVA_TOP_N_ENRICHR, "_interaction"),
        out_dir    = mode_e_dir,
        fill_color = "#1B7837"
      )
    } else {
      message("  [SKIP] Enrichr not available or gene list empty — skipping enrichment for Mode E.")
    }

    # --- E.6: Violin plot for top interaction gene ---------------------------
    if (nrow(anova_wide) > 0) {
      top_gene <- anova_wide$gene[1]
      message(paste("  Plotting top interaction gene:", top_gene))

      pval_row <- anova_wide %>% dplyr::filter(gene == top_gene)
      genotype_p_col     <- paste0(ANOVA_FACTOR1, ".p.adj")
      diet_p_col         <- paste0(ANOVA_FACTOR2, ".p.adj")

      plot_subtitle <- sprintf(
        "2-Way ANOVA Adj. P-values:\n  %s = %.2e\n  %s = %.2e\n  Interaction = %.2e",
        ANOVA_FACTOR1, pval_row[[genotype_p_col]],
        ANOVA_FACTOR2, pval_row[[diet_p_col]],
        pval_row[[interaction_p_adj_col]]
      )

      plotting_sub <- subset(data, cells = cells_anova)
      plot_df <- FetchData(plotting_sub,
                            vars = c(top_gene, ANOVA_FACTOR1, ANOVA_FACTOR2))
      colnames(plot_df)[1] <- "Expression"

      p_anova <- ggplot(plot_df, aes(x = !!sym(ANOVA_FACTOR2), y = Expression,
                                      fill = !!sym(ANOVA_FACTOR2))) +
        geom_violin(trim = TRUE, scale = "width") +
        geom_jitter(height = 0, width = 0.1, size = 0.1, alpha = 0.3) +
        facet_wrap(as.formula(paste("~", ANOVA_FACTOR1))) +
        labs(
          title    = paste("Expression of", top_gene, "in", ANOVA_CELL_TYPE),
          subtitle = plot_subtitle,
          x        = ANOVA_FACTOR2, y = "Log-Normalized Expression"
        ) +
        theme_classic() +
        theme(
          plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1),
          strip.text    = element_text(face = "bold", size = 12),
          legend.position = "none"
        )
      print(p_anova)
      ggsave(
        file.path(mode_e_dir, paste0(top_gene, "_", gsub("[^A-Za-z0-9_]", "_", ANOVA_CELL_TYPE),
                                      "_ANOVA_VlnPlot.png")),
        plot = p_anova, width = 10, height = 7, dpi = DPI_SETTING
      )
      rm(plotting_sub, plot_df)
    }
    message("  Mode E complete.")
  }
}

# =============================================================================
# --- PART 10: FINAL MESSAGE --------------------------------------------------
# =============================================================================
message(paste0(
  "\n=== DIFFERENTIAL EXPRESSION + ANOVA COMPLETE ===\n",
  "  Output directory: ", DE_DIR, "\n\n",
  if (RUN_MODE_A_WILCOXON) paste0(
    "  Mode A  (Wilcoxon):  ",
    file.path(DE_DIR, paste0("ModeA_Wilcoxon_", GROUP_1, "_vs_", GROUP_2)), "\n"