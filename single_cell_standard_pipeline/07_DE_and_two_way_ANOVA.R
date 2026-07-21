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
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
# --- 1.1: Mode ---------------------------------------------------------------
# "broad"   → runs on broad cell type labels (CellType column)
# "subtype" → runs on sub-annotated object (sub_cell_types column)
MODE <- "broad"   # <-- ACTION: change to "subtype" for subtype runs
#MODE <- "subtype"   # <-- ACTION: change to "subtype" for subtype runs

# --- 1.2: Paths --------------------------------------------------------------
PROJECT_NAME <- "Nr4a1_s17_ack"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# RDS to load — change per run
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))

# For subtypes, use one of:
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))

# Output directory — automatically set based on MODE
if (MODE == "broad") {
  DE_OUT_DIR <- file.path(OUTPUT_DIR, "DE_results", "broad")
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
CONDITION_COLUMN  <- "Genotype_sex"
SAMPLE_COLUMN     <- "SampleID"
DE_LOGFC_THRESH   <- 0.5     # discovery threshold — moderate+ effects in xlsx
DE_MIN_PCT        <- 0.10
DE_PADJ_THRESH    <- 0.05
MIN_CELLS_GROUP   <- 20      # skip group if fewer cells than this
DV_PVAL_THRESH    <- 0.05    # SplineDV raw p-value threshold for DV genes

# --- 1.5: Contrasts ----------------------------------------------------------
# Sex-separated: Polyp vs WT, KO vs Polyp
CONTRASTS_LIST <- list(
  Female_Polyp_vs_WT         = c("Polyp_Female",         "WT_Female"),
  Female_KO_vs_Polyp         = c("Polyp_NR4a1_KO_Female", "Polyp_Female"),
  Male_Polyp_vs_WT           = c("Polyp_Male",            "WT_Male"),
  Male_KO_vs_Polyp           = c("Polyp_NR4a1_KO_Male",   "Polyp_Male"),
  Female_Polyp_vs_Male_Polyp = c("Polyp_Female",          "Polyp_Male")
)

# --- 1.6: ANOVA Parameters ---------------------------------------------------
RUN_ANOVA         <- TRUE
USE_TWO_WAY_ANOVA <- FALSE   # FALSE = one-way (~Genotype_sex)
# TRUE  = two-way (~Genotype * Sex)
ANOVA_FACTOR1     <- "Genotype"  # used only if USE_TWO_WAY_ANOVA = TRUE
ANOVA_FACTOR2     <- "Sex"       # used only if USE_TWO_WAY_ANOVA = TRUE
ANOVA_MIN_CELLS   <- 20
ANOVA_TOP_N       <- 250         # top N genes for Enrichr

# Sex-separated ANOVA groups — each run independently
# WT is first in each vector = reference level automatically
ANOVA_GROUPS <- list(
  Female = c("WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female"),
  Male   = c("WT_Male",   "Polyp_Male",   "Polyp_NR4a1_KO_Male")
)

# --- 1.7: Enrichr Parameters -------------------------------------------------
RUN_ENRICHR       <- TRUE
ENRICHR_DBS       <- c("GO_Biological_Process_2026",
                       "GO_Molecular_Function_2026",
                       "KEGG_2026")
ENRICHR_TOP_N     <- 10     # top N pathways to plot per gene list

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
run_mast <- function(so, ident.1, ident.2,
                     group.by     = CONDITION_COLUMN,
                     latent.vars  = SAMPLE_COLUMN,
                     logfc.thresh = DE_LOGFC_THRESH,
                     min.pct      = DE_MIN_PCT) {
  
  Idents(so) <- so[[group.by, drop = TRUE]]
  
  # Check minimum cells
  n1 <- sum(Idents(so) == ident.1)
  n2 <- sum(Idents(so) == ident.2)
  if (n1 < MIN_CELLS_GROUP || n2 < MIN_CELLS_GROUP) {
    message(paste("  [SKIP] Too few cells:", ident.1, "=", n1, "|", ident.2, "=", n2,
                  "(min =", MIN_CELLS_GROUP, ")"))
    return(NULL)
  }
  
  markers <- tryCatch({
    FindMarkers(
      object          = so,
      ident.1         = ident.1,
      ident.2         = ident.2,
      test.use        = "MAST",
      latent.vars     = latent.vars,
      logfc.threshold = logfc.thresh,
      min.pct         = min.pct,
      verbose         = FALSE
    )
  }, error = function(e) {
    message(paste("  [ERROR] MAST failed:", e$message))
    return(NULL)
  })
  
  if (is.null(markers) || nrow(markers) == 0) return(NULL)
  
  markers %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      contrast   = paste0(ident.1, "_vs_", ident.2),
      direction  = ifelse(avg_log2FC > 0,
                          paste0("up_in_", ident.1),
                          paste0("up_in_", ident.2)),
      confidence = dplyr::case_when(
        abs(avg_log2FC) >= 1.0 ~ "high",      # ≥2-fold: publication standard
        abs(avg_log2FC) >= 0.5 ~ "moderate",  # ≥1.4-fold: discovery
        TRUE                   ~ "low"
      )
    ) %>%
    dplyr::arrange(p_val_adj, desc(abs(avg_log2FC)))
}

# =============================================================================
# --- PART 3b: SPLINEDV WRAPPER -----------------------------------------------
# =============================================================================
run_splinedv <- function(seurat_obj, cell_type, celltype_col, condition_col,
                         group1, group2) {
  if (!requireNamespace("SplineDV", quietly = TRUE)) {
    stop("[ERROR] SplineDV not installed. Run 00_rlibs_installation.R to install it.")
  }
  library(SplineDV)
  
  # Subset to cell type
  cells_ct <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] == cell_type]
  obj_ct   <- subset(seurat_obj, cells = cells_ct)
  
  # Subset to the two groups being compared
  obj_ct <- subset(obj_ct, cells = rownames(obj_ct@meta.data)[
    obj_ct@meta.data[[condition_col]] %in% c(group1, group2)])
  
  # Check minimum cells per group
  n1 <- sum(obj_ct@meta.data[[condition_col]] == group1)
  n2 <- sum(obj_ct@meta.data[[condition_col]] == group2)
  if (n1 < MIN_CELLS_GROUP || n2 < MIN_CELLS_GROUP) {
    message(paste("  [SKIP SplineDV]", cell_type, "| insufficient cells:",
                  group1, "=", n1, "|", group2, "=", n2,
                  "(min =", MIN_CELLS_GROUP, ")"))
    return(NULL)
  }
  message(paste("  [SplineDV]", cell_type, "|", group1, "(n=", n1, ") vs",
                group2, "(n=", n2, ")"))
  
  # Get cells for each group
  cells_g1 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group1]
  cells_g2 <- rownames(obj_ct@meta.data)[obj_ct@meta.data[[condition_col]] == group2]
  
  # SplineDV expects raw count matrices (not normalized)
  counts_g1 <- GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g1, drop = FALSE]
  counts_g2 <- GetAssayData(obj_ct, assay = "RNA", layer = "counts")[, cells_g2, drop = FALSE]
  
  # Run splineDV
  dv_result <- tryCatch(
    splineDV(counts_g1, counts_g2),
    error = function(e) {
      message("    [ERROR] SplineDV failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dv_result) || nrow(dv_result) == 0) return(NULL)
  
  # Convert S4 DFrame to standard data.frame for dplyr compatibility
  dv_result <- as.data.frame(dv_result)
  
  # --- Standardize column names based on confirmed splineDV output ---
  # Columns are: genes, mu1, mu2, CV1, CV2, drop1, drop2, dist1, dist2,
  #   X_splinex, X_spliney, X_splinez, Y_splinex, Y_spliney, Y_splinez,
  #   vectorDist, Direction, Pval
  
  # Gene column: "genes" -> "gene"
  if ("genes" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, gene = genes)
  } else if ("Gene" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, gene = Gene)
  } else if (!"gene" %in% colnames(dv_result)) {
    dv_result$gene <- rownames(dv_result)
  }
  
  # P-value column: "Pval" -> "pval"
  if ("Pval" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = Pval)
  } else if ("PValue" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = PValue)
  } else if ("pvalue" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = pvalue)
  } else if ("p.value" %in% colnames(dv_result)) {
    dv_result <- dplyr::rename(dv_result, pval = p.value)
  }
  # If still no pval, assign NA
  
  if (!"pval" %in% colnames(dv_result)) {
    message("    [WARN] No p-value column found in SplineDV output. Columns: ",
            paste(colnames(dv_result), collapse = ", "))
    dv_result$pval <- NA_real_
  }
  
  # Compute adjusted p-value (BH) — splineDV does not provide one
  if (any(!is.na(dv_result$pval))) {
    dv_result$padj <- p.adjust(dv_result$pval, method = "BH")
  } else {
    dv_result$padj <- NA_real_
  }
  
  # Add metadata columns
  dv_result$cell_type  <- cell_type
  dv_result$comparison <- paste0(group1, "_vs_", group2)
  dv_result$n_cells_g1 <- n1
  dv_result$n_cells_g2 <- n2
  
  # Sort by pval
  if (any(!is.na(dv_result$pval))) {
    dv_result <- dv_result %>% dplyr::arrange(pval)
  }
  
  return(dv_result)
}

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

run_enrichr_ora <- function(gene_list, label) {
  # Returns a flat data frame with a 'database' column, or NULL if skipped
  if (!RUN_ENRICHR || length(gene_list) == 0) return(NULL)
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  if (length(gene_list) < 5) {
    message(paste("  [SKIP Enrichr] Too few genes for:", label))
    return(NULL)
  }
  Sys.sleep(2)   # respect Enrichr rate limit
  tryCatch({
    enrich_res <- enrichr(gene_list, ENRICHR_DBS)
    flat_rows  <- list()
    for (db in names(enrich_res)) {
      df <- enrich_res[[db]]
      if (!is.null(df) && nrow(df) > 0) {
        df <- df %>%
          dplyr::filter(Adjusted.P.value < 0.05) %>%
          dplyr::arrange(Adjusted.P.value) %>%
          head(100) %>%
          dplyr::mutate(database = db, .before = 1)
        flat_rows[[db]] <- df
      }
    }
    if (length(flat_rows) > 0) {
      message(paste("  Enrichr done:", label))
      return(dplyr::bind_rows(flat_rows))
    }
    return(NULL)
  }, error = function(e) {
    message(paste("  [ERROR Enrichr]", e$message))
    return(NULL)
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
    result <- run_mast(so_ct, ident.1 = ident1, ident.2 = ident2)
    
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
      enr_up <- run_enrichr_ora(sig_up, paste0(contrast_name, "_UP"))
      enr_dn <- run_enrichr_ora(sig_dn, paste0(contrast_name, "_DOWN"))
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
      group2        = ident2
    )
    
    if (!is.null(dv_result)) {
      n_dv <- sum(dv_result$pval < DV_PVAL_THRESH, na.rm = TRUE)
      message(paste("  → SplineDV:", n_dv, "DV genes (pval <", DV_PVAL_THRESH, ")"))
      ct_dv_sheets[[contrast_name]] <- dv_result
      
      # SplineDV standalone Enrichr
      dv_sig_genes <- dv_result$gene[dv_result$pval < DV_PVAL_THRESH & !is.na(dv_result$pval)]
      enr_dv <- run_enrichr_ora(dv_sig_genes, paste0(contrast_name, "_DV_only"))
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
          enr_ov_up <- run_enrichr_ora(ov_up$gene, paste0(contrast_name, "_overlap_UP"))
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
          enr_ov_dn <- run_enrichr_ora(ov_dn$gene, paste0(contrast_name, "_overlap_DOWN"))
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
  write_xlsx(lapply(ct_de_sheets, as.data.frame), de_file)
  message(paste("  Saved:", basename(de_file)))
  
  # SplineDV file: one sheet per contrast, sorted by DV pval
  if (any(sapply(ct_dv_sheets, function(x) !is.null(x) && nrow(x) > 1))) {
    dv_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_SplineDV.xlsx"))
    write_xlsx(lapply(ct_dv_sheets, as.data.frame), dv_file)
    message(paste("  Saved:", basename(dv_file)))
  }
  
  # DV ∩ DE Overlap file
  if (length(ct_overlap_sheets) > 0) {
    ov_file <- file.path(DE_OUT_DIR, paste0(safe_ct, "_DV_DE_Overlap.xlsx"))
    write_xlsx(lapply(ct_overlap_sheets, as.data.frame), ov_file)
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

make_pathway_heatmap <- function(enr_list, title, out_file) {
  if (length(enr_list) == 0) {
    message(paste("  [SKIP heatmap] No enrichr data for:", title))
    return(invisible(NULL))
  }
  all_enr <- dplyr::bind_rows(enr_list)
  if (nrow(all_enr) == 0) return(invisible(NULL))
  
  # Filter to GOBP only — most interpretable for biological pathway stories
  all_enr <- all_enr %>%
    dplyr::filter(database == "GO_Biological_Process_2026")
  
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP heatmap] No GOBP results for:", title))
    return(invisible(NULL))
  }
  
  all_enr <- all_enr %>%
    dplyr::mutate(
      gene_count     = as.integer(sub("/.*", "", Overlap)),
      neg_log10_padj = -log10(pmax(Adjusted.P.value, 1e-10))
    ) %>%
    dplyr::filter(Adjusted.P.value < 0.05, gene_count >= HEATMAP_MIN_GENE_OVERLAP)
  
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP heatmap] No significant pathways after filtering for:", title))
    return(invisible(NULL))
  }
  
  if ("direction" %in% colnames(all_enr)) {
    all_enr <- all_enr %>%
      dplyr::mutate(group_key = paste0(cell_type, "\n", contrast, "\n", direction))
  } else {
    all_enr <- all_enr %>%
      dplyr::mutate(group_key = paste0(cell_type, "\n", contrast))
  }
  
  deduplicate_pathways <- function(df) {
    if (nrow(df) <= 1) return(df)
    df <- df %>%
      dplyr::mutate(score = gene_count * neg_log10_padj) %>%
      dplyr::arrange(desc(score))
    gene_lists <- lapply(df$Genes, function(g) unlist(strsplit(g, ";\\s*|,\\s*")))
    keep <- rep(TRUE, nrow(df))
    n    <- nrow(df)
    for (i in seq_len(n)) {
      if (!keep[i]) next
      if (i == n) break
      for (j in seq(i + 1, n)) {
        if (!keep[j]) next
        gi <- gene_lists[[i]]; gj <- gene_lists[[j]]
        shared <- length(intersect(gi, gj))
        pct_i  <- if (length(gi) > 0) shared / length(gi) else 0
        pct_j  <- if (length(gj) > 0) shared / length(gj) else 0
        if (pct_j > HEATMAP_REDUNDANCY_RATIO || pct_i > HEATMAP_REDUNDANCY_RATIO)
          keep[j] <- FALSE
      }
    }
    df[keep, ] %>% dplyr::select(-score)
  }
  
  all_enr_dedup <- all_enr %>%
    dplyr::group_by(group_key) %>%
    dplyr::group_modify(~ deduplicate_pathways(.x) %>% head(HEATMAP_TOP_N_PATHWAYS)) %>%
    dplyr::ungroup()
  
  median_gene_count <- median(all_enr_dedup$gene_count, na.rm = TRUE)
  
  pathway_global <- all_enr_dedup %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(
      n_cell_types    = dplyr::n_distinct(cell_type),
      mean_gene_count = mean(gene_count, na.rm = TRUE),
      mean_score      = mean(gene_count * neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      n_cell_types >= HEATMAP_MIN_CELL_TYPES |
        mean_gene_count >= median_gene_count
    ) %>%
    dplyr::arrange(desc(n_cell_types), desc(mean_score)) %>%
    head(HEATMAP_TOP_N_PATHWAYS)
  
  selected_terms <- pathway_global$Term
  if (length(selected_terms) == 0) {
    message(paste("  [SKIP heatmap] No cross-representative pathways for:", title))
    return(invisible(NULL))
  }
  message(paste("  Pathways selected:", length(selected_terms)))
  
  heatmap_df <- all_enr_dedup %>%
    dplyr::filter(Term %in% selected_terms) %>%
    dplyr::group_by(Term, group_key) %>%
    dplyr::summarise(
      neg_log10_padj = max(neg_log10_padj, na.rm = TRUE),
      gene_count     = max(gene_count,     na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::complete(Term, group_key,
                    fill = list(neg_log10_padj = 0, gene_count = 0))
  
  term_order         <- pathway_global$Term
  term_order_wrapped <- stringr::str_wrap(term_order, width = 35)
  heatmap_df$Term    <- stringr::str_wrap(heatmap_df$Term, width = 35)
  heatmap_df$Term    <- factor(heatmap_df$Term, levels = rev(term_order_wrapped))
  heatmap_df$group_key <- factor(heatmap_df$group_key)
  heatmap_df$label   <- ifelse(heatmap_df$gene_count > 0,
                               as.character(heatmap_df$gene_count), "")
  
  n_terms  <- length(unique(heatmap_df$Term))
  n_groups <- length(unique(heatmap_df$group_key))
  
  p_heat <- ggplot(heatmap_df,
                   aes(x = group_key, y = Term, fill = neg_log10_padj)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.5, color = "black", fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#d73027",
                        name = "-log10\n(adj.p)", na.value = "grey95") +
    labs(title    = title,
         subtitle = paste0("Top ", n_terms, " pathways | Gene count in cells | ",
                           "logFC ≥ 1.0 genes only | GOBP only | ",
                           "Redundancy: ", HEATMAP_REDUNDANCY_RATIO * 100, "%"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y   = element_text(size = 13),
      legend.title  = element_text(size = 13),
      legend.text   = element_text(size = 12),
      panel.grid    = element_blank(),
      plot.margin   = margin(10, 10, 10, 20)
    )
  
  w <- max(10, n_groups * 1.6 + 4)
  h <- max(8,  n_terms  * 0.7 + 5)
  
  ggsave(out_file, p_heat, width = w, height = h,
         dpi = DPI_SETTING, bg = "white", limitsize = FALSE)
  message(paste("  Heatmap saved:", basename(out_file)))
  
  write.csv(
    heatmap_df %>%
      dplyr::left_join(
        pathway_global %>%
          dplyr::select(Term, n_cell_types, mean_gene_count) %>%
          dplyr::mutate(Term = stringr::str_wrap(Term, width = 35)),
        by = "Term"),
    sub("\\.png$", "_data.csv", out_file), row.names = FALSE
  )
  return(invisible(p_heat))
}

# --- 1: DE MAST — one per contrast per direction ---------------------------
message("\n=== Generating DE MAST pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_up   <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|UP$"),   names(enr_de_all))]
  enr_down <- enr_de_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"), names(enr_de_all))]
  make_pathway_heatmap(enr_up,
                       paste("DE MAST — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_DE_MAST_UP_", contrast_name, ".png")))
  make_pathway_heatmap(enr_down,
                       paste("DE MAST — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_DE_MAST_DOWN_", contrast_name, ".png")))
}

# --- 2: DV+DE Overlap — one per contrast per direction --------------------
message("\n=== Generating DV+DE Overlap pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_up   <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|UP$"),   names(enr_overlap_all))]
  enr_down <- enr_overlap_all[grep(paste0("\\|", contrast_name, "\\|DOWN$"), names(enr_overlap_all))]
  make_pathway_heatmap(enr_up,
                       paste("DV+DE Overlap — UP —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_Overlap_UP_", contrast_name, ".png")))
  make_pathway_heatmap(enr_down,
                       paste("DV+DE Overlap — DOWN —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_Overlap_DOWN_", contrast_name, ".png")))
}

# --- 3: SplineDV — one per contrast (no direction) ------------------------
message("\n=== Generating SplineDV standalone pathway heatmaps ===")
for (contrast_name in names(CONTRASTS_LIST)) {
  enr_dv <- enr_dv_all[grep(paste0("\\|", contrast_name, "$"), names(enr_dv_all))]
  make_pathway_heatmap(enr_dv,
                       paste("SplineDV — DV Pathways —", contrast_name, "—", MODE),
                       file.path(DE_OUT_DIR, paste0("HEATMAP_SplineDV_", contrast_name, ".png")))
}

# =============================================================================
# --- GENE HUB HEATMAPS -------------------------------------------------------
# =============================================================================
gene_hub_dir <- file.path(DE_OUT_DIR, "gene_heatmaps")
if (!dir.exists(gene_hub_dir)) dir.create(gene_hub_dir, recursive = TRUE)

make_gene_hub_heatmap <- function(enr_list, de_list, title, out_file) {
  # enr_list: named list of enrichr flat dfs (cell_type, contrast, direction cols)
  # de_list:  named list of DE result dfs (gene, avg_log2FC, p_val_adj, cell_type)
  
  if (length(enr_list) == 0 || length(de_list) == 0) {
    message(paste("  [SKIP gene hub] Insufficient data for:", title))
    return(invisible(NULL))
  }
  
  # --- 1. Build pathway universe from enrichr results -----------------------
  all_enr <- dplyr::bind_rows(enr_list) %>%
    dplyr::filter(database == "GO_Biological_Process_2026") %>%   # GOBP only
    dplyr::mutate(
      gene_count     = as.integer(sub("/.*", "", Overlap)),
      neg_log10_padj = -log10(pmax(Adjusted.P.value, 1e-10))
    ) %>%
    dplyr::filter(Adjusted.P.value < 0.05, gene_count >= HEATMAP_MIN_GENE_OVERLAP)
  
  if (nrow(all_enr) == 0) {
    message(paste("  [SKIP gene hub] No significant pathways for:", title))
    return(invisible(NULL))
  }
  
  # Select top pathways — same scoring as pathway heatmap
  top_pathways <- all_enr %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(
      mean_gene_count = mean(gene_count, na.rm = TRUE),
      mean_score      = mean(gene_count * neg_log10_padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_score)) %>%
    head(GENE_HUB_TOP_N_PATHWAYS)
  
  selected_pathways <- top_pathways$Term
  
  # --- 2. Build gene universe from DE results --------------------------------
  all_de <- dplyr::bind_rows(de_list) %>%
    dplyr::filter(p_val_adj <= 0.05, abs(avg_log2FC) >= GENE_HUB_MIN_LOGFC)
  
  if (nrow(all_de) == 0) {
    message(paste("  [SKIP gene hub] No qualifying DE genes for:", title))
    return(invisible(NULL))
  }
  
  # --- 3. Build gene × pathway membership matrix ----------------------------
  # Parse pathway gene lists
  pathway_gene_lists <- all_enr %>%
    dplyr::filter(Term %in% selected_pathways) %>%
    dplyr::select(Term, Genes) %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      gene_vec = list(trimws(unlist(strsplit(Genes, ";|,"))))
    ) %>%
    dplyr::ungroup()
  
  # Get unique qualifying DE genes
  de_genes <- unique(all_de$gene)
  
  # For each gene × pathway: is the gene in the pathway?
  membership <- tidyr::expand_grid(
    gene = de_genes,
    Term = selected_pathways
  ) %>%
    dplyr::left_join(
      pathway_gene_lists %>% dplyr::select(Term, gene_vec),
      by = "Term"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      member = gene %in% gene_vec[[1]]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, Term, member)
  
  # --- 4. Score genes by hub connectivity -----------------------------------
  gene_hub_scores <- membership %>%
    dplyr::filter(member) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_pathways = dplyr::n_distinct(Term),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_pathways >= GENE_HUB_MIN_PATHWAYS) %>%
    dplyr::left_join(
      all_de %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(
          avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "gene"
    ) %>%
    dplyr::mutate(hub_score = n_pathways * abs(avg_log2FC)) %>%
    dplyr::arrange(desc(hub_score)) %>%
    head(GENE_HUB_TOP_N_GENES)
  
  if (nrow(gene_hub_scores) == 0) {
    message(paste("  [SKIP gene hub] No hub genes found (need >=",
                  GENE_HUB_MIN_PATHWAYS, "pathway memberships) for:", title))
    return(invisible(NULL))
  }
  
  top_genes    <- gene_hub_scores$gene
  selected_pathways_final <- selected_pathways[
    selected_pathways %in% (membership %>%
                              dplyr::filter(gene %in% top_genes, member) %>%
                              dplyr::pull(Term) %>% unique())
  ]
  
  # --- 5. Build plot matrix -------------------------------------------------
  plot_df <- membership %>%
    dplyr::filter(gene %in% top_genes, Term %in% selected_pathways_final) %>%
    dplyr::left_join(
      all_de %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
                         .groups = "drop"),
      by = "gene"
    ) %>%
    dplyr::mutate(
      fill_val = ifelse(member, avg_log2FC, NA_real_),
      label    = ifelse(member, sprintf("%.1f", avg_log2FC), "")
    )
  
  # Order genes by hub score (top hub at top)
  gene_order    <- gene_hub_scores$gene
  pathway_order <- stringr::str_wrap(selected_pathways_final, width = GENE_HUB_WRAP_WIDTH)
  
  plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order))
  plot_df$Term <- factor(
    stringr::str_wrap(plot_df$Term, width = GENE_HUB_WRAP_WIDTH),
    levels = pathway_order
  )
  
  # Add n_pathways annotation to gene labels
  gene_labels <- gene_hub_scores %>%
    dplyr::mutate(label_full = paste0(gene, " (", n_pathways, ")")) %>%
    dplyr::select(gene, label_full)
  levels(plot_df$gene) <- gene_labels$label_full[
    match(levels(plot_df$gene), gene_labels$gene)
  ]
  
  n_genes    <- length(unique(plot_df$gene))
  n_pathways <- length(unique(plot_df$Term))
  
  # --- 6. Plot --------------------------------------------------------------
  p_hub <- ggplot(plot_df, aes(x = Term, y = gene)) +
    geom_tile(aes(fill = fill_val), color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.8, color = "black", fontface = "bold") +
    scale_fill_gradient2(
      low      = "#4575b4",
      mid      = "white",
      high     = "#d73027",
      midpoint = 0,
      name     = "avg\nlog2FC",
      na.value = "grey92"     # NA cells automatically get grey
    ) +
    labs(
      title    = title,
      subtitle = paste0(
        "Genes: padj ≤ 0.05 | logFC ≥ ", GENE_HUB_MIN_LOGFC,
        " | Red = up-regulated | Blue = down-regulated | ",
        "Number = logFC (1dp) | Gene label: name (n pathways) | ",
        "Hub score = n_pathways × |logFC| | GOBP only"
      ),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 17, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y   = element_text(size = 11),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 11),
      panel.grid    = element_blank(),
      plot.margin   = margin(10, 10, 10, 20)
    )
  
  w <- max(10, n_pathways * 1.4 + 4)
  h <- max(8,  n_genes    * 0.5 + 5)
  
  ggsave(out_file, p_hub, width = w, height = h,
         dpi = DPI_SETTING, bg = "white", limitsize = FALSE)
  message(paste("  Gene hub heatmap saved:", basename(out_file),
                "| Hub genes:", n_genes, "| Pathways:", n_pathways))
  
  # Save data
  write.csv(
    plot_df %>% dplyr::select(gene, Term, member, fill_val, label),
    sub("\\.png$", "_data.csv", out_file), row.names = FALSE
  )
  return(invisible(p_hub))
}

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
  
  anova_dir <- file.path(DE_OUT_DIR, "ANOVA")
  if (!dir.exists(anova_dir)) dir.create(anova_dir, recursive = TRUE)
  
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
    enr_anova <- run_enrichr_ora(top_genes, paste0(label, "_ANOVA_top", ANOVA_TOP_N))
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

