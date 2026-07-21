# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 07: ANNOTATION UNIFIER
# Version: 1.0
# UNIFIED BUILD: part of unified_pipeline/. Consumes the object produced by
#   01_process_data.R v11.0. Doublet calls arrive standardised in the
#   'Doublet_Status' column regardless of which caller ran, so this script
#   requires no changes when DOUBLET_METHOD is switched.
#
# PURPOSE:
#   Merges sub-annotations from T cells, Macrophages, and Colonocytes back
#   into the broad annotated object. Cells that were removed during sub-
#   annotation QC (doublets, dirty clusters) are dropped from the unified
#   object. Cells from other broad types retain their original CellType label.
#
#   The output is a single Seurat object with:
#     - data$CellType        = fine-grained label (sub_cell_types for T/Mac/Col,
#                              original CellType for all others)
#     - data$CellType_broad  = original broad label (always preserved)
#     - data$Genotype_sex    = preserved from original
#     - all other metadata   = preserved from original
#
# HOW TO USE:
#   1. Set paths in Part 1.
#   2. Run the full script.
#   3. Use the output RDS as input for CellChat (Script 08).
#
# INPUT:
#   - {PROJECT_NAME}_final_annotated.rds          (broad object)
#   - {PROJECT_NAME}_T_cells_subclustered.rds
#   - {PROJECT_NAME}_Macrophages_subclustered.rds
#   - {PROJECT_NAME}_Colonocytes_subclustered.rds
#
# OUTPUT:
#   - {PROJECT_NAME}_unified_annotated.rds
# =============================================================================

library(Seurat)
library(dplyr)

# =============================================================================
# --- PART 1: CONFIGURATION ---------------------------------------------------
# =============================================================================
PROJECT_NAME <- "Nr4a1_s17_ack"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# Input files
BROAD_RDS    <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))
TCELL_RDS    <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds"))
MAC_RDS      <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_macrophages_subclustered.rds"))
COLON_RDS    <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))

# Output file
UNIFIED_RDS  <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_unified_annotated.rds"))

# Broad cell type labels that have sub-annotations
# These must match exactly what's in data$CellType in the broad object
SUBANNOTATED_BROAD_TYPES <- c("T cells", "Macrophages", "Colonocytes")

# Column in sub-objects that holds the final sub-annotation
SUB_LABEL_COLUMN <- "sub_cell_types"

# =============================================================================
# --- PART 2: LOAD ALL OBJECTS ------------------------------------------------
# =============================================================================
message("=== Loading broad annotated object ===")
data_broad <- readRDS(BROAD_RDS)
message(paste("  Broad object:", ncol(data_broad), "cells,",
              length(unique(data_broad$CellType)), "cell types"))
message(paste("  Broad cell types:", paste(sort(unique(data_broad$CellType)), collapse = ", ")))

message("\n=== Loading sub-annotated objects ===")
data_tcell <- readRDS(TCELL_RDS)
message(paste("  T cells:", ncol(data_tcell), "cells,",
              length(unique(data_tcell[[SUB_LABEL_COLUMN, drop=TRUE]])), "sub-types"))

data_mac   <- readRDS(MAC_RDS)
message(paste("  Macrophages:", ncol(data_mac), "cells,",
              length(unique(data_mac[[SUB_LABEL_COLUMN, drop=TRUE]])), "sub-types"))

data_colon <- readRDS(COLON_RDS)
message(paste("  Colonocytes:", ncol(data_colon), "cells,",
              length(unique(data_colon[[SUB_LABEL_COLUMN, drop=TRUE]])), "sub-types"))

# =============================================================================
# --- PART 3: BUILD BARCODE → FINE LABEL MAPPING ------------------------------
# =============================================================================
message("\n=== Building barcode → fine label mapping ===")

# Extract barcodes and labels from each sub-object
tcell_map <- setNames(
  as.character(data_tcell[[SUB_LABEL_COLUMN, drop = TRUE]]),
  colnames(data_tcell)
)
mac_map <- setNames(
  as.character(data_mac[[SUB_LABEL_COLUMN, drop = TRUE]]),
  colnames(data_mac)
)
colon_map <- setNames(
  as.character(data_colon[[SUB_LABEL_COLUMN, drop = TRUE]]),
  colnames(data_colon)
)

# Combine all sub-annotation maps
all_sub_maps <- c(tcell_map, mac_map, colon_map)

# Check for barcode collisions across sub-objects (should not happen)
if (any(duplicated(names(all_sub_maps)))) {
  warning("[WARN] Duplicate barcodes found across sub-objects — check for overlapping cells.")
}

message(paste("  T cell barcodes in sub-object:", length(tcell_map)))
message(paste("  Macrophage barcodes in sub-object:", length(mac_map)))
message(paste("  Colonocyte barcodes in sub-object:", length(colon_map)))

# =============================================================================
# --- PART 4: IDENTIFY CELLS TO KEEP -----------------------------------------
# =============================================================================
message("\n=== Identifying cells to keep ===")

broad_barcodes  <- colnames(data_broad)
broad_celltypes <- as.character(data_broad$CellType)
names(broad_celltypes) <- broad_barcodes

total_broad <- length(broad_barcodes)

# Cells from non-subannotated broad types — keep all
non_sub_cells <- broad_barcodes[!broad_celltypes %in% SUBANNOTATED_BROAD_TYPES]

# Cells from subannotated broad types — keep only those that survived sub QC
sub_cells_in_broad <- broad_barcodes[broad_celltypes %in% SUBANNOTATED_BROAD_TYPES]
sub_cells_survived <- sub_cells_in_broad[sub_cells_in_broad %in% names(all_sub_maps)]
sub_cells_removed  <- sub_cells_in_broad[!sub_cells_in_broad %in% names(all_sub_maps)]

# Final cells to keep
cells_to_keep <- c(non_sub_cells, sub_cells_survived)

# --- Per-compartment breakdown -----------------------------------------------
message("\n  ┌─────────────────────────────────────────────────────────────────┐")
message(  "  │              CELL REMOVAL SUMMARY PER COMPARTMENT              │")
message(  "  ├──────────────────────┬───────────┬──────────┬──────────────────┤")
message(  "  │ Compartment          │  In broad │  Removed │  % of broad type │")
message(  "  ├──────────────────────┼───────────┼──────────┼──────────────────┤")

removal_summary <- data.frame(
  compartment     = character(),
  n_broad         = integer(),
  n_survived      = integer(),
  n_removed       = integer(),
  pct_of_type     = numeric(),
  pct_of_total    = numeric(),
  stringsAsFactors = FALSE
)

for (broad_type in SUBANNOTATED_BROAD_TYPES) {
  cells_broad <- broad_barcodes[broad_celltypes == broad_type]
  n_broad    <- length(cells_broad)
  n_survived <- sum(cells_broad %in% names(all_sub_maps))
  n_removed  <- n_broad - n_survived
  pct_type   <- n_removed / n_broad * 100
  pct_total  <- n_removed / total_broad * 100
  
  message(sprintf("  │ %-20s │ %9d │ %8d │ %14.1f%% │",
                  broad_type, n_broad, n_removed, pct_type))
  
  removal_summary <- rbind(removal_summary, data.frame(
    compartment  = broad_type,
    n_broad      = n_broad,
    n_survived   = n_survived,
    n_removed    = n_removed,
    pct_of_type  = pct_type,
    pct_of_total = pct_total,
    stringsAsFactors = FALSE
  ))
}

# Non-subannotated
n_non_sub <- length(non_sub_cells)
message(sprintf("  │ %-20s │ %9d │ %8d │ %14.1f%% │",
                "Other (not subann.)", n_non_sub, 0, 0.0))
message(  "  ├──────────────────────┼───────────┼──────────┼──────────────────┤")

# Totals
total_removed  <- length(sub_cells_removed)
total_kept     <- length(cells_to_keep)
pct_removed_total <- total_removed / total_broad * 100
pct_kept_total    <- total_kept    / total_broad * 100

message(sprintf("  │ %-20s │ %9d │ %8d │ %14.1f%% │",
                "TOTAL", total_broad, total_removed, pct_removed_total))
message(  "  └──────────────────────┴───────────┴──────────┴──────────────────┘")
message(sprintf("\n  Cells kept in unified object: %d / %d (%.1f%%)",
                total_kept, total_broad, pct_kept_total))
message(sprintf("  Cells dropped from broad object: %d (%.1f%% of total)",
                total_removed, pct_removed_total))

# Also print % of total dataset dropped per compartment
message("\n  Per-compartment removal as % of total broad dataset:")
for (i in seq_len(nrow(removal_summary))) {
  message(sprintf("    %-20s: %d removed = %.2f%% of total %d cells",
                  removal_summary$compartment[i],
                  removal_summary$n_removed[i],
                  removal_summary$pct_of_total[i],
                  total_broad))
}

# =============================================================================
# --- PART 5: SUBSET AND ASSIGN FINE LABELS -----------------------------------
# =============================================================================
message("\n=== Subsetting broad object to kept cells ===")
data_unified <- subset(data_broad, cells = cells_to_keep)

# Preserve broad label
data_unified$CellType_broad <- data_unified$CellType

# Assign fine labels
message("=== Assigning fine-grained sub-type labels ===")
fine_labels <- as.character(data_unified$CellType)
names(fine_labels) <- colnames(data_unified)

# Replace subannotated cells with their sub-type label
sub_mask <- names(fine_labels) %in% names(all_sub_maps)
fine_labels[sub_mask] <- all_sub_maps[names(fine_labels)[sub_mask]]

data_unified$CellType <- fine_labels
data_unified$CellType <- factor(data_unified$CellType)

message(paste("  Fine-grained cell types:", length(unique(data_unified$CellType))))
message(paste("  Labels:", paste(sort(unique(data_unified$CellType)), collapse = ", ")))

# =============================================================================
# --- PART 6: VALIDATION ------------------------------------------------------
# =============================================================================
message("\n=== Validation ===")

# Check no broad labels remain for subannotated compartments
broad_labels_remaining <- data_unified$CellType[
  data_unified$CellType_broad %in% SUBANNOTATED_BROAD_TYPES
]
if (any(broad_labels_remaining %in% SUBANNOTATED_BROAD_TYPES)) {
  warning("[WARN] Some cells still have broad labels — check mapping.")
} else {
  message("  ✓ All subannotated cells have fine-grained labels")
}

# Print composition table
message("\n  Final cell type composition:")
ct_table <- sort(table(data_unified$CellType), decreasing = TRUE)
for (i in seq_along(ct_table)) {
  pct <- ct_table[i] / ncol(data_unified) * 100
  message(sprintf("    %-35s %6d cells (%5.1f%%)",
                  names(ct_table)[i], ct_table[i], pct))
}

# =============================================================================
# --- PART 7: SAVE ------------------------------------------------------------
# =============================================================================
message("\n=== Saving unified object ===")
saveRDS(data_unified, UNIFIED_RDS)
message(paste("  Saved:", UNIFIED_RDS))
message(paste("  Final:", ncol(data_unified), "cells,",
              length(unique(data_unified$CellType)), "cell types"))
message("\n=== UNIFIER COMPLETE ===")
message("  Next step: Run 08_cellchat.R using the unified object.")