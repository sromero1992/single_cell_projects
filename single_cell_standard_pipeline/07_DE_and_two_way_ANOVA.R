# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 6: MAST DIFFERENTIAL EXPRESSION ANALYSIS
# Version: 1.1 (SplineDV fix)
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
#
# PURPOSE:
#   Runs MAST-based differential expression across all cell types in a
#   Seurat object, accounting for sample-level variability via latent.vars.
#   Designed to run on EITHER the broad annotated object OR any subtype RDS.
#
# DESIGN:
#   4 contrasts per cell type (sex-separated):
#     Female: Polyp_Female         vs WT_Female
#     Female: Polyp_NR4a1_KO_Female vs Polyp_Female
#     Male:   Polyp_Male           vs WT_Male
#     Male:   Polyp_NR4a1_KO_Male  vs Polyp_Male
#
#   Loops automatically over all cell types found in CELLTYPE_COLUMN.
#   Outputs one xlsx per cell type with one sheet per contrast.
#
# ANOVA:
#   Optional one-way or two-way ANOVA controlled by USE_TWO_WAY_ANOVA flag.
#   One-way: ~ Genotype_sex (recommended for this design)
#   Two-way: ~ Genotype * Sex (formal interaction term)
#
# HOW TO USE:
#   BROAD:   Set MODE = "broad", point RDS_PATH to final_annotated.rds
#            CELLTYPE_COLUMN = "CellType", OUT_DIR = DE_DIR/broad/
#   SUBTYPE: Set MODE = "subtype", point RDS_PATH to subclustered.rds
#            CELLTYPE_COLUMN = "sub_cell_types"
#            OUT_DIR = e.g. tcell_subannotation/DE/
#
# DEPENDENCIES:
#   Seurat, MAST, dplyr, tidyr, writexl, tibble, ggplot2, patchwork
#   Optional ANOVA: parallel, broom, enrichR
# =============================================================================
library(TamuScDSC)   # DE/DV/enrichment helpers now live in the package
library(Seurat)
library(MAST)
library(dplyr)
library(tidyr)
library(tibble)
library(writexl)
library(ggplot2)
library(patchwork)
set.seed(123)

# =============================================================================
# --- HELPERS: Ensembl mapping + sheet finalization ---------------------------
# =============================================================================
# symbol_to_ensembl(): mouse SYMBOL -> ENSMUSG, trying the official symbol first
# then ALIAS for anything unmapped. Degrades to NA (never errors) if the
# annotation packages are absent, so a missing org.Mm.eg.db cannot stop the run.
# Adapted from add_ensembl_for_pathvisio.R.
.ensembl_ok <- ADD_ENSEMBL &&
  requireNamespace("AnnotationDbi", quietly = TRUE) &&
  requireNamespace("org.Mm.eg.db", quietly = TRUE)
if (ADD_ENSEMBL && !.ensembl_ok) {
  message("  [NOTE] ADD_ENSEMBL=TRUE but org.Mm.eg.db/AnnotationDbi not installed; ",
          "Ensembl column will be skipped. Install via 00_rlibs_installation.R.")
}

# finalize_sheet(): guarantee 'gene' is the FIRST column and 'Ensembl' the
# SECOND in every exported table. Applied to every DE/DV/overlap sheet just
# before writing, so all downstream files are PathVisio-ready with one call.
# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
# --- 1.1: Mode ---------------------------------------------------------------
# "broad"   → runs on broad cell type labels (CellType column)
# "subtype" → runs on sub-annotated object (sub_cell_types column)
MODE <- "broad"   # <-- ACTION: change to "subtype" for subtype runs
#MODE <- "subtype"   # <-- ACTION: change to "subtype" for subtype runs

# --- 1.2: Paths --------------------------------------------------------------
PROJECT_NAME <- "Wu_Diet_project2"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# RDS to load — change per run
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))

# For subtypes, use one of:
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))

# Output directory — automatically set based on MODE.
# Everything from Script 07 goes into ONE master folder (no sub-subfolders).
if (MODE == "broad") {
  DE_OUT_DIR <- file.path(OUTPUT_DIR, "DE_results")
} else if (MODE == "subtype") {
  # Detect subtype from RDS filename automatically
  rds_base   <- tools::file_path_sans_ext(basename(RDS_PATH))
  DE_OUT_DIR <- file.path(OUTPUT_DIR,
                          gsub(paste0(PROJECT_NAME, "_"), "", rds_base) %>%
                            gsub("_subclustered", "_subannotation", .) %>%
                            paste0("/DE"))
}
if (!dir.exists(DE_OUT_DIR)) dir.create(DE_OUT_DIR, recursive = TRUE)
message(paste("  Output folder:", DE_OUT_DIR))

# --- 1.3: Cell type column ---------------------------------------------------
# Automatically set based on MODE; override here if needed
CELLTYPE_COLUMN <- if (MODE == "broad") "CellType" else "sub_cell_types"

# --- 1.4: MAST Parameters ----------------------------------------------------
CONDITION_COLUMN  <- "Diet"       # Wu diet study: Cellulose / Fiber / PUFA / Mix
SAMPLE_COLUMN     <- "SampleID"
DE_LOGFC_THRESH   <- 0.5     # discovery threshold — moderate+ effects in xlsx
DE_MIN_PCT        <- 0.10
DE_PADJ_THRESH    <- 0.05
MIN_CELLS_GROUP   <- 20      # skip group if fewer cells than this
DV_PVAL_THRESH    <- 0.05    # SplineDV raw p-value threshold for DV genes

# --- 1.5: Contrasts ----------------------------------------------------------
# Wu diet study. Cellulose is the control diet, so every contrast is
# <treatment> vs Cellulose: a positive log2FC means UP in the treatment diet.
# CONFIRM these names match the Diet column of Wu_p2_metadata.xlsx exactly
# (case-sensitive) before running.
CONTRASTS_LIST <- list(
  Fiber_vs_Cellulose = c("Fiber", "Cellulose"),
  PUFA_vs_Cellulose  = c("PUFA",  "Cellulose"),
  Mix_vs_Cellulose   = c("Mix",   "Cellulose"),
  PUFA_vs_Fiber      = c("PUFA",  "Fiber"),
  Mix_vs_PUFA        = c("Mix",   "PUFA")
)

# --- 1.6: ANOVA Parameters ---------------------------------------------------
RUN_ANOVA         <- TRUE
USE_TWO_WAY_ANOVA <- FALSE   # FALSE = one-way (~Genotype_sex)
# TRUE  = two-way (~Genotype * Sex)
ANOVA_FACTOR1     <- "Diet"      # used only if USE_TWO_WAY_ANOVA = TRUE
ANOVA_FACTOR2     <- "SetDay"    # used only if USE_TWO_WAY_ANOVA = TRUE
ANOVA_MIN_CELLS   <- 20
ANOVA_TOP_N       <- 250         # top N genes for Enrichr

# One-way ANOVA across all four diets in a single model.
# Cellulose is listed FIRST, which makes it the reference level.
ANOVA_GROUPS <- list(
  AllDiets = c("Cellulose", "Fiber", "PUFA", "Mix")
)

# --- 1.7: Enrichr Parameters -------------------------------------------------
RUN_ENRICHR       <- TRUE
# WikiPathways added for PathVisio interoperability (WikiPathways pathways are
# what PathVisio consumes). Use the current mouse release.
ENRICHR_DBS       <- c("GO_Biological_Process_2026",
                       "GO_Molecular_Function_2026",
                       "KEGG_2026",
                       "WikiPathways_2024_Mouse")
ENRICHR_TOP_N     <- 10     # top N pathways to plot per gene list

# --- 1.7b: Ensembl IDs for PathVisio -----------------------------------------
# ADD_ENSEMBL: map each gene symbol to its mouse Ensembl ID (ENSMUSG...) and
# write it as the SECOND column of every DE / DV / overlap sheet. PathVisio and
# WikiPathways key on Ensembl, and symbols are lossy (aliases/renames), so
# providing Ensembl directly maximizes the pathway match rate.
# Requires org.Mm.eg.db + AnnotationDbi. If missing, Ensembl is skipped with a
# warning and the gene-symbol column is still written first.
ADD_ENSEMBL <- TRUE

# --- 1.8: Output Resolution --------------------------------------------------
DPI_SETTING <- 300

# --- 1.9: Gene Hub Heatmap Parameters ----------------------------------------
GENE_HUB_TOP_N_GENES    <- 30    # top hub genes to show
GENE_HUB_TOP_N_PATHWAYS <- 20    # top pathways to show
GENE_HUB_MIN_PATHWAYS   <- 2     # gene must appear in >=N pathways
GENE_HUB_MIN_LOGFC      <- 0.5   # padj <= 0.05 AND |logFC| >= this
GENE_HUB_WRAP_WIDTH     <- 30    # pathway name wrap width

# =============================================================================
# --- PART 2: LOAD DATA -------------------------------------------------------
# =============================================================================
message("=== Loading Seurat object ===")
data <- readRDS(RDS_PATH)
message(paste("  Loaded:", ncol(data), "cells"))

# Ensure CellType column exists
if (!CELLTYPE_COLUMN %in% colnames(data@meta.data)) {
  stop(paste("[ERROR] Column", CELLTYPE_COLUMN, "not found in metadata."))
}

# --- Resolve the comparison column (generalised) -----------------------------
# CONDITION_COLUMN is used everywhere as data@meta.data[[CONDITION_COLUMN]], so
# it only has to EXIST. Three cases, handled in order:
#   1. It already exists in the metadata            -> use it as-is (most cases).
#   2. It is a compound "ColA_ColB" and both parts  -> build it by pasting the
#      exist as columns (e.g. "Genotype_Diet",         two columns with "_".
#      "Genotype_Sex", "Diet_Sex")                     This is the common
#                                                       genotype+diet / genotype+
#                                                       sex situation.
#   3. Neither                                       -> stop with the available
#                                                       columns listed.
data <- resolve_comparison_column(data, CONDITION_COLUMN)

# The compound column must also exist for the ANOVA factors, when two-way.
if (RUN_ANOVA && USE_TWO_WAY_ANOVA) {
  for (fac in c(ANOVA_FACTOR1, ANOVA_FACTOR2)) {
    if (!fac %in% colnames(data@meta.data)) {
      stop(paste0("[ERROR] Two-way ANOVA factor '", fac, "' not found in metadata."),
           call. = FALSE)
    }
  }
}

# Normalize if not already done
if (max(GetAssayData(data, layer = "data")) <= 0) {
  message("  Normalizing data...")
  data <- NormalizeData(data, verbose = FALSE)
}

# Get all cell types to loop over
all_cell_types <- sort(unique(as.character(data@meta.data[[CELLTYPE_COLUMN]])))
message(paste("  Cell types found:", paste(all_cell_types, collapse = ", ")))

# =============================================================================
# --- PART 3: MAST WRAPPER ----------------------------------------------------
# =============================================================================

# =============================================================================
# --- PART 3b: SPLINEDV WRAPPER -----------------------------------------------
# =============================================================================

# =============================================================================
# --- PART 4: ENRICHR WRAPPER -------------------------------------------------
# =============================================================================
if (RUN_ENRICHR) {
  library(enrichR)
  tryCatch({
    setEnrichrSite("Enrichr")
    message("  Enrichr connection established.")
  }, error = function(e) {
    message("  [WARN] Enrichr connection failed — enrichment will be skipped.")
    RUN_ENRICHR <<- FALSE
  })
}


# =============================================================================
# --- PART 5: MAIN MAST + SPLINEDV LOOP ---------------------------------------
# =============================================================================
de_summary <- data.frame(
  cell_type  = character(),
  contrast   = character(),
  direction  = character(),
  n_genes    = integer(),
  stringsAsFactors = FALSE
)

# Global collectors for pathway heatmaps
enr_de_all      <- list()   # key = "CellType|contrast|direction"
enr_overlap_all <- list()   # key = "CellType|contrast|direction"
enr_dv_all      <- list()   # key = "CellType|contrast"
de_data_all     <- list()   # key = "CellType|contrast|direction" — raw DE results

message("\n=== STEP 5: Running MAST DE + SplineDV across all cell types ===")
for (ct in all_cell_types) {
  safe_ct <- gsub("[^A-Za-z0-9_]", "_", ct)
  message(paste("\n--- Cell type:", ct, "---"))
  
  cells_ct <- rownames(data@meta.data)[data@meta.data[[CELLTYPE_COLUMN]] == ct]
  so_ct    <- subset(data, cells = cells_ct)
  
  ct_de_sheets       <- list()   # MAST DE — one entry per contrast
  ct_dv_sheets       <- list()   # SplineDV — one entry per contrast
  ct_overlap_sheets  <- list()   # DV ∩ DE — one entry per contrast_direction
  ct_enrichr_sheets  <- list()   # MAST Enrichr — one entry per contrast_direction
  ct_ov_enr_sheets   <- list()   # Overlap Enrichr — one entry per contrast_direction
  
  for (contrast_name in names(CONTRASTS_LIST)) {
    groups <- CONTRASTS_LIST[[contrast_name]]
    ident1 <- groups[1]
    ident2 <- groups[2]
    message(paste("  Contrast:", contrast_name, "|", ident1, "vs", ident2))
    
    # --- MAST DE -------------------------------------------------------------
    result <- run_mast(so_ct, ident.1 = ident1, ident.2 = ident2,
                       group_by = CONDITION_COLUMN, latent_vars = SAMPLE_COLUMN,
                       logfc_thresh = DE_LOGFC_THRESH, min_pct = DE_MIN_PCT,
                       min_cells = MIN_CELLS_GROUP)
    
    if (!is.null(result)) {
      n_sig <- sum(result$p_val_adj < DE_PADJ_THRESH, na.rm = TRUE)
      n_up  <- sum(result$p_val_adj < DE_PADJ_THRESH & result$avg_log2FC > 0, na.rm = TRUE)
      n_dn  <- sum(result$p_val_adj < DE_PADJ_THRESH & result$avg_log2FC < 0, na.rm = TRUE)
      message(paste("  → MAST:", n_sig, "sig genes | up:", n_up, "| down:", n_dn))
      
      de_summary <- rbind(de_summary,
                          data.frame(cell_type = ct, contrast = contrast_name,
                                     direction = "Up",   n_genes = n_up),
                          data.frame(cell_type = ct, contrast = contrast_name,
                                     direction = "Down", n_genes = n_dn)
      )
      
      # Save all significant genes to xlsx (logFC >= 0.5 + padj < 0.05)
      ct_de_sheets[[contrast_name]] <- result %>%
        dplyr::filter(p_val_adj < DE_PADJ_THRESH) %>%
        dplyr::arrange(desc(avg_log2FC))
      
      # Enrichr uses only high confidence genes (logFC >= 1.0) — cleaner pathways
      sig_up <- result$gene[result$p_val_adj < DE_PADJ_THRESH &
                              result$avg_log2FC > 0 &
                              result$confidence == "high"]
      sig_dn <- result$gene[result$p_val_adj < DE_PADJ_THRESH &
                              result$avg_log2FC < 0 &
                              result$confidence == "high"]
      
      # MAST Enrichr
      enr_up <- run_enrichr_ora(sig_up, paste0(contrast_name, "_UP"), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
      enr_dn <- run_enrichr_ora(sig_dn, paste0(contrast_name, "_DOWN"), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
      if (!is.null(enr_up)) {
        ct_enrichr_sheets[[paste0(contrast_name, "_UP")]] <- enr_up
        enr_de_all[[paste0(ct, "|", contrast_name, "|UP")]] <- enr_up %>%
          dplyr::mutate(cell_type = ct, contrast = contrast_name, direction = "UP")
      }
      if (!is.null(enr_dn)) {
        ct_enrichr_sheets[[paste0(contrast_name, "_DOWN")]] <- enr_dn
        enr_de_all[[paste0(ct, "|", contrast_name, "|DOWN")]] <- enr_dn %>%
          dplyr::mutate(cell_type = ct, contrast = contrast_name, direction = "DOWN")
      }
      
      # Collect raw DE data for gene hub heatmaps
      de_sig <- result %>%
        dplyr::filter(p_val_adj <= 0.05, abs(avg_log2FC) >= GENE_HUB_MIN_LOGFC) %>%
        dplyr::mutate(cell_type = ct, contrast = contrast_name)
      if (nrow(de_sig) > 0) {
        de_data_all[[paste0(ct, "|", contrast_name, "|UP")]] <- de_sig %>%
          dplyr::filter(avg_log2FC > 0)
        de_data_all[[paste0(ct, "|", contrast_name, "|DOWN")]] <- de_sig %>%
          dplyr::filter(avg_log2FC < 0)
      }
      
      
    } else {
      message(paste("  → MAST: No results for:", contrast_name))
      ct_de_sheets[[contrast_name]] <- data.frame(
        gene = NA_character_,
        note = "No results — skipped or too few cells",
        stringsAsFactors = FALSE
      )
      sig_up <- character(0)
      sig_dn <- character(0)
    }
    
    # --- SplineDV ------------------------------------------------------------
    dv_result <- run_splinedv(
      seurat_obj    = so_ct,
      cell_type     = ct,
      celltype_col  = CELLTYPE_COLUMN,
      condition_col = CONDITION_COLUMN,
      group1        = ident1,
      group2        = ident2,
      min_cells     = MIN_CELLS_GROUP
    )
    
    if (!is.null(dv_result)) {
      n_dv <- sum(dv_result$pval <= DV_PVAL_THRESH, na.rm = TRUE)
      message(paste("  → SplineDV:", n_dv, "DV genes (pval <=", DV_PVAL_THRESH, ")"))
      ct_dv_sheets[[contrast_name]] <- dv_result
      
      # SplineDV standalone Enrichr
      dv_sig_genes <- dv_result$gene[dv_result$pval <= DV_PVAL_THRESH & !is.na(dv_result$pval)]
      enr_dv <- run_enrichr_ora(dv_sig_genes, paste0(contrast_name, "_DV_only"), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
      if (!is.null(enr_dv)) {
        enr_dv_all[[paste0(ct, "|", contrast_name)]] <- enr_dv %>%
          dplyr::mutate(cell_type = ct, contrast = contrast_name)
      }
      
      # --- DV ∩ DE Overlap --------------------------------------------------
      if (!is.null(result)) {
        dv_genes <- dv_result$gene[dv_result$pval < DV_PVAL_THRESH & !is.na(dv_result$pval)]
        
        # Prepare DV columns to join
        dv_to_join <- dv_result %>%
          dplyr::select(-any_of(c("cell_type", "comparison", "n_cells_g1", "n_cells_g2")))
        dv_cols_to_rename <- setdiff(colnames(dv_to_join), "gene")
        names(dv_cols_to_rename) <- paste0("DV_", dv_cols_to_rename)
        dv_to_join <- dplyr::rename(dv_to_join, !!!dv_cols_to_rename)
        
        # Overlap UP — high confidence DE only
        ov_up <- result %>%
          dplyr::filter(p_val_adj < DE_PADJ_THRESH &
                          avg_log2FC > 0 &
                          confidence == "high" &
                          gene %in% dv_genes) %>%
          dplyr::left_join(dv_to_join, by = "gene") %>%
          dplyr::arrange(desc(avg_log2FC))
        if (nrow(ov_up) > 0) {
          ct_overlap_sheets[[paste0(contrast_name, "_UP")]] <- ov_up
          message(paste("  → Overlap UP:", nrow(ov_up), "genes"))
          enr_ov_up <- run_enrichr_ora(ov_up$gene, paste0(contrast_name, "_overlap_UP"), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
          if (!is.null(enr_ov_up)) {
            ct_ov_enr_sheets[[paste0(contrast_name, "_UP")]] <- enr_ov_up
            enr_overlap_all[[paste0(ct, "|", contrast_name, "|UP")]] <- enr_ov_up %>%
              dplyr::mutate(cell_type = ct, contrast = contrast_name, direction = "UP")
          }
        }
        
        # Overlap DOWN — high confidence DE only
        ov_dn <- result %>%
          dplyr::filter(p_val_adj < DE_PADJ_THRESH &
                          avg_log2FC < 0 &
                          confidence == "high" &
                          gene %in% dv_genes) %>%
          dplyr::left_join(dv_to_join, by = "gene") %>%
          dplyr::arrange(avg_log2FC)
        if (nrow(ov_dn) > 0) {
          ct_overlap_sheets[[paste0(contrast_name, "_DOWN")]] <- ov_dn
          message(paste("  → Overlap DOWN:", nrow(ov_dn), "genes"))
          enr_ov_dn <- run_enrichr_ora(ov_dn$gene, paste0(contrast_name, "_overlap_DOWN"), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
          if (!is.null(enr_ov_dn)) {
            ct_ov_enr_sheets[[paste0(contrast_name, "_DOWN")]] <- enr_ov_dn
            enr_overlap_all[[paste0(ct, "|", contrast_name, "|DOWN")]] <- enr_ov_dn %>%
              dplyr::mutate(cell_type = ct, contrast = contrast_name, direction = "DOWN")
          }
        }
      }
    } else {
      message(paste("  → SplineDV: No results for:", contrast_name))
      ct_dv_sheets[[contrast_name]] <- data.frame(
        gene = NA_character_,
        note = "No results — skipped or too few cells",
        stringsAsFactors = FALSE
      )
    }  # end if/else (!is.null(dv_result))
    
  }  # end for (contrast_name)
  
  # --- Save all files for this cell type -------------------------------------
  
  # DE file: one sheet per contrast, sorted log2FC desc
  de_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_MAST_DE.xlsx"))
  write_xlsx(finalize_all(ct_de_sheets, add_ensembl = ADD_ENSEMBL), de_file)
  message(paste("  Saved:", basename(de_file)))
  
  # SplineDV file: one sheet per contrast, sorted by DV pval
  if (any(sapply(ct_dv_sheets, function(x) !is.null(x) && nrow(x) > 1))) {
    dv_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_SplineDV.xlsx"))
    write_xlsx(finalize_all(ct_dv_sheets, add_ensembl = ADD_ENSEMBL), dv_file)
    message(paste("  Saved:", basename(dv_file)))
  }
  
  # DV ∩ DE Overlap file
  if (length(ct_overlap_sheets) > 0) {
    ov_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_DV_DE_Overlap.xlsx"))
    write_xlsx(finalize_all(ct_overlap_sheets, add_ensembl = ADD_ENSEMBL), ov_file)
    message(paste("  Saved:", basename(ov_file)))
  }
  
  # MAST Enrichr file
  if (length(ct_enrichr_sheets) > 0) {
    enr_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_MAST_Enrichr.xlsx"))
    write_xlsx(lapply(ct_enrichr_sheets, as.data.frame), enr_file)
    message(paste("  Saved:", basename(enr_file)))
  }
  
  # Overlap Enrichr file
  if (length(ct_ov_enr_sheets) > 0) {
    ov_enr_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_Overlap_Enrichr.xlsx"))
    write_xlsx(lapply(ct_ov_enr_sheets, as.data.frame), ov_enr_file)
    message(paste("  Saved:", basename(ov_enr_file)))
  }
  
  rm(so_ct, ct_de_sheets, ct_dv_sheets, ct_overlap_sheets,
     ct_enrichr_sheets, ct_ov_enr_sheets)
  gc()
  
}  # end for (ct)
# =============================================================================
# --- SUMMARY PLOT: DE gene counts per cell type per contrast -----------------
# =============================================================================
if (nrow(de_summary) > 0) {
  de_summary$n_genes_signed <- ifelse(
    de_summary$direction == "Down",
    -de_summary$n_genes,
    de_summary$n_genes
  )
  de_summary$cell_type <- factor(de_summary$cell_type,
                                 levels = rev(sort(unique(de_summary$cell_type))))
  de_summary$direction <- factor(de_summary$direction, levels = c("Up", "Down"))
  
  p_summary <- ggplot(de_summary,
                      aes(x = n_genes_signed, y = cell_type, fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
    scale_x_continuous(labels = abs) +
    facet_wrap(~ contrast, ncol = 2) +
    labs(title  = paste("DE Gene Counts per Cell Type —", MODE),
         x      = "Number of Significant Genes",
         y      = NULL,
         fill   = "Direction") +
    theme_bw(base_size = 13) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text   = element_text(size = 13, face = "bold"),
      axis.text.y  = element_text(size = 11),
      legend.position = "top"
    )
  
  n_celltypes <- length(unique(de_summary$cell_type))
  ggsave(file.path(DE_OUT_DIR, "SUMMARY_DE_gene_counts.png"),
         p_summary,
         width  = 14,
         height = max(6, n_celltypes * 0.4 + 3),
         dpi    = DPI_SETTING)
  message("  Summary plot saved: SUMMARY_DE_gene_counts.png")
  
  # Also save the raw counts as csv
  write.csv(
    de_summary %>% dplyr::select(-n_genes_signed) %>%
      tidyr::pivot_wider(names_from = direction, values_from = n_genes),
    file.path(DE_OUT_DIR, "SUMMARY_DE_gene_counts.csv"),
    row.names = FALSE
  )
}

# =============================================================================
# --- PATHWAY HEATMAPS --------------------------------------------------------
# =============================================================================
library(stringr)

HEATMAP_TOP_N_PATHWAYS   <- 20
HEATMAP_MIN_GENE_OVERLAP <- 3
HEATMAP_REDUNDANCY_RATIO <- 0.5
HEATMAP_MIN_CELL_TYPES   <- 2


# =============================================================================
# --- PATHWAY BARPLOTS (companion to the heatmaps) ----------------------------
# =============================================================================
# Top-N enriched GOBP pathways as a horizontal barplot:
#   bar length = Odds Ratio (ORA effect size; falls back to Combined Score,
#                then gene overlap if a column is absent)
#   bar colour = adjusted p-value (red = most significant)
#   number at bar end = gene count (overlap with the pathway)
# Pools cell types by keeping each pathway's most significant occurrence, then
# takes the top ENRICHR_TOP_N by significance. One per contrast/direction, saved
# next to the matching heatmap.

# --- 1: DE MAST — one per contrast per direction ---------------------------
message("\n=== Generating DE MAST pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_up   <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|UP$"),   names(enr_de_all))]
  enr_down <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"), names(enr_de_all))]
  make_pathway_heatmap(enr_up,
                       paste("DE MAST — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_DE_MAST_UP_", contrast_name, ".png")))
  make_pathway_barplot(enr_up,
                       paste("DE MAST — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("BARPLOT_DE_MAST_UP_", contrast_name, ".png")))
  make_pathway_heatmap(enr_down,
                       paste("DE MAST — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_DE_MAST_DOWN_", contrast_name, ".png")))
  make_pathway_barplot(enr_down,
                       paste("DE MAST — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("BARPLOT_DE_MAST_DOWN_", contrast_name, ".png")))
}

# --- 2: DV+DE Overlap — one per contrast per direction --------------------
message("\n=== Generating DV+DE Overlap pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_up   <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|UP$"),   names(enr_overlap_all))]
  enr_down <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"), names(enr_overlap_all))]
  make_pathway_heatmap(enr_up,
                       paste("DV+DE Overlap — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_Overlap_UP_", contrast_name, ".png")))
  make_pathway_barplot(enr_up,
                       paste("DV+DE Overlap — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("BARPLOT_Overlap_UP_", contrast_name, ".png")))
  make_pathway_heatmap(enr_down,
                       paste("DV+DE Overlap — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_Overlap_DOWN_", contrast_name, ".png")))
  make_pathway_barplot(enr_down,
                       paste("DV+DE Overlap — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("BARPLOT_Overlap_DOWN_", contrast_name, ".png")))
}

# --- 3: SplineDV — one per contrast (no direction) ------------------------
message("\n=== Generating SplineDV standalone pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_dv <- enr_dv_all[grep(paste0("\\|", contrast_name, "$"), names(enr_dv_all))]
  make_pathway_heatmap(enr_dv,
                       paste("SplineDV — DV Pathways —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_SplineDV_", contrast_name, ".png")))
  make_pathway_barplot(enr_dv,
                       paste("SplineDV — DV Pathways —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("BARPLOT_SplineDV_", contrast_name, ".png")))
}

# =============================================================================
# --- GENE HUB HEATMAPS -------------------------------------------------------
# =============================================================================
gene_hub_dir <- DE_OUT_DIR   # gene-hub heatmaps go in the master DE folder


# --- Run gene hub heatmaps per contrast -------------------------------------
message("\n=== Generating Gene Hub heatmaps ===")

for (contrast_name in names(CONTRASTS_LIST)) {
  
  # DE MAST UP
  enr_up  <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|UP$"),  names(enr_de_all))]
  de_up   <- de_data_all[grep(paste0("\\|", contrast_name, "\\|UP$"), names(de_data_all))]
  make_gene_hub_heatmap(
    enr_list = enr_up, de_list = de_up,
    title    = paste("Gene Hubs — DE MAST UP —", contrast_name, "—", MODE),
    out_file = file.path(gene_hub_dir,
                         paste0("GENE_HUB_DE_MAST_UP_", contrast_name, ".png"))
  )
  
  # DE MAST DOWN
  enr_down <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"),  names(enr_de_all))]
  de_down  <- de_data_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"), names(de_data_all))]
  make_gene_hub_heatmap(
    enr_list = enr_down, de_list = de_down,
    title    = paste("Gene Hubs — DE MAST DOWN —", contrast_name, "—", MODE),
    out_file = file.path(gene_hub_dir,
                         paste0("GENE_HUB_DE_MAST_DOWN_", contrast_name, ".png"))
  )
  
  # Overlap UP
  enr_ov_up <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|UP$"),
                                    names(enr_overlap_all))]
  make_gene_hub_heatmap(
    enr_list = enr_ov_up, de_list = de_up,
    title    = paste("Gene Hubs — Overlap UP —", contrast_name, "—", MODE),
    out_file = file.path(gene_hub_dir,
                         paste0("GENE_HUB_Overlap_UP_", contrast_name, ".png"))
  )
  
  # Overlap DOWN
  enr_ov_dn <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"),
                                    names(enr_overlap_all))]
  make_gene_hub_heatmap(
    enr_list = enr_ov_dn, de_list = de_down,
    title    = paste("Gene Hubs — Overlap DOWN —", contrast_name, "—", MODE),
    out_file = file.path(gene_hub_dir,
                         paste0("GENE_HUB_Overlap_DOWN_", contrast_name, ".png"))
  )
}
if (RUN_ANOVA) {
  library(broom)
  
  anova_dir <- DE_OUT_DIR   # ANOVA outputs go in the master DE folder
  
  # --- Auto-derive Genotype and Sex columns if needed (for two-way ANOVA) ---
  if (USE_TWO_WAY_ANOVA) {
    if (!"Genotype" %in% colnames(data@meta.data)) {
      message("  Deriving 'Genotype' column from Genotype_sex...")
      data$Genotype <- sub("_Female$|_Male$", "", data$Genotype_sex)
    }
    if (!"Sex" %in% colnames(data@meta.data)) {
      message("  Deriving 'Sex' column from Genotype_sex...")
      data$Sex <- ifelse(grepl("Female", data$Genotype_sex), "Female", "Male")
    }
    anova_formula <- as.formula(paste("expression ~", ANOVA_FACTOR1, "*", ANOVA_FACTOR2))
    sort_term     <- paste0(ANOVA_FACTOR1, ":", ANOVA_FACTOR2)
    message("\n=== STEP 6: Two-way ANOVA (~Genotype * Sex) across all cell types ===")
  } else {
    anova_formula <- as.formula(paste("expression ~", CONDITION_COLUMN))
    sort_term     <- CONDITION_COLUMN
    message("\n=== STEP 6: One-way ANOVA (~Genotype_sex, sex-separated) across all cell types ===")
  }
  
  # --- Reusable ANOVA runner ------------------------------------------------
  run_anova_for_object <- function(so_anova, genes_to_test, label) {
    message(paste("  Running ANOVA for:", label,
                  "| Genes:", length(genes_to_test)))
    
    all_results <- vector("list", length(genes_to_test))
    names(all_results) <- genes_to_test
    
    for (i in seq_along(genes_to_test)) {
      gene <- genes_to_test[i]
      
      fetch_vars <- if (USE_TWO_WAY_ANOVA) {
        c(gene, ANOVA_FACTOR1, ANOVA_FACTOR2)
      } else {
        c(gene, CONDITION_COLUMN)
      }
      
      gdf <- FetchData(so_anova, vars = fetch_vars)
      colnames(gdf)[1] <- "expression"
      
      aov_m <- tryCatch(aov(anova_formula, data = gdf), error = function(e) NULL)
      if (!is.null(aov_m)) {
        tidy_r       <- broom::tidy(aov_m)
        tidy_r$gene  <- gene
        all_results[[i]] <- tidy_r
      }
      
      if (i %% 1000 == 0)
        message(paste0("    ", i, " / ", length(genes_to_test), " genes done..."))
    }
    
    anova_df <- dplyr::bind_rows(all_results)
    
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
      )
    
    sort_col <- paste0(sort_term, ".p.adj")
    if (sort_col %in% colnames(anova_wide)) {
      anova_wide <- anova_wide %>% dplyr::arrange(!!sym(sort_col))
    }
    
    # Save CSV
    anova_file <- file.path(anova_dir, paste0(label, "_ANOVA.csv"))
    write.csv(anova_wide, anova_file, row.names = FALSE)
    message(paste("  Saved:", basename(anova_file)))
    
    # Enrichr on top N genes — flat table saved alongside ANOVA csv
    top_genes <- head(anova_wide$gene[!is.na(anova_wide[[sort_col]])], ANOVA_TOP_N)
    enr_anova <- run_enrichr_ora(top_genes, paste0(label, "_ANOVA_top", ANOVA_TOP_N), dbs = ENRICHR_DBS, enabled = RUN_ENRICHR)
    if (!is.null(enr_anova)) {
      enr_file <- file.path(anova_dir, paste0(label, "_ANOVA_Enrichr.xlsx"))
      write_xlsx(list(Enrichr = as.data.frame(enr_anova)), enr_file)
      message(paste("  ANOVA Enrichr saved:", basename(enr_file)))
    }
    
    return(invisible(NULL))
  }
  
  # --- Main ANOVA loop ------------------------------------------------------
  for (ct in all_cell_types) {
    safe_ct <- gsub("[^A-Za-z0-9_]", "_", ct)
    message(paste("\n  ANOVA for:", ct))
    
    cells_ct <- rownames(data@meta.data)[data@meta.data[[CELLTYPE_COLUMN]] == ct]
    so_ct    <- subset(data, cells = cells_ct)
    
    if (USE_TWO_WAY_ANOVA) {
      # --- Two-way: all cells together, ~ Genotype * Sex --------------------
      group_counts <- table(so_ct@meta.data[[ANOVA_FACTOR1]],
                            so_ct@meta.data[[ANOVA_FACTOR2]])
      if (any(group_counts < ANOVA_MIN_CELLS)) {
        message(paste("  [SKIP] Some Genotype × Sex combinations have <",
                      ANOVA_MIN_CELLS, "cells"))
        rm(so_ct); gc(); next
      }
      
      genes_ct <- rownames(so_ct)[
        Matrix::rowSums(GetAssayData(so_ct, layer = "counts") > 0) /
          ncol(so_ct) >= DE_MIN_PCT
      ]
      
      run_anova_for_object(so_ct, genes_ct, paste0(safe_ct, "_TwoWay"))
      
    } else {
      # --- One-way: sex-separated, ~ Genotype_sex ---------------------------
      for (sex_group in names(ANOVA_GROUPS)) {
        levels_to_keep <- ANOVA_GROUPS[[sex_group]]
        
        counts_sex <- table(so_ct@meta.data[[CONDITION_COLUMN]][
          so_ct@meta.data[[CONDITION_COLUMN]] %in% levels_to_keep
        ])
        
        if (length(counts_sex) < length(levels_to_keep) ||
            any(counts_sex < ANOVA_MIN_CELLS)) {
          message(paste("  [SKIP] Too few cells for", sex_group))
          next
        }
        
        cells_sex <- rownames(so_ct@meta.data)[
          so_ct@meta.data[[CONDITION_COLUMN]] %in% levels_to_keep
        ]
        so_sex <- subset(so_ct, cells = cells_sex)
        
        # Re-factor — WT first = reference
        so_sex@meta.data[[CONDITION_COLUMN]] <- factor(
          so_sex@meta.data[[CONDITION_COLUMN]],
          levels = levels_to_keep
        )
        
        genes_sex <- rownames(so_sex)[
          Matrix::rowSums(GetAssayData(so_sex, layer = "counts") > 0) /
            ncol(so_sex) >= DE_MIN_PCT
        ]
        
        run_anova_for_object(so_sex, genes_sex,
                             paste0(safe_ct, "_", sex_group))
        
        rm(so_sex); gc()
      }  # end sex_group loop
    }  # end one-way branch
    
    rm(so_ct); gc()
  }  # end cell type loop
}

# =============================================================================
# --- PART 7: FINAL MESSAGE ---------------------------------------------------
# =============================================================================
message(paste0(
  "\n=== MAST DE ANALYSIS COMPLETE ===\n",
  "  Mode:       ", MODE, "\n",
  "  RDS loaded: ", basename(RDS_PATH), "\n",
  "  Output:     ", DE_OUT_DIR, "\n",
  "  Cell types: ", length(all_cell_types), "\n",
  "  Contrasts:  ", paste(names(CONTRASTS_LIST), collapse = ", "), "\n",
  "  ANOVA:      ", ifelse(RUN_ANOVA,
                           ifelse(USE_TWO_WAY_ANOVA, "Two-way", "One-way"),
                           "Skipped"), "\n"
))

