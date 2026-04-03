# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING ENGINE
# Version: 8.0 (Generalized Pipeline)
#
# PURPOSE:
#   This script is the computational backbone of the single-cell RNA-seq
#   pipeline. It sequentially performs:
#     1. Per-sample data loading from 10x Genomics H5 files
#     2. Optional probe/custom feature integration (e.g., KO-target probes)
#     3. Copy-number variation detection via SCEVAN
#     4. Pre-merge QC (loose filters) + post-merge QC (stringent filters)
#     5. Doublet detection via DoubletFinder with a ROLLBACK SAFEGUARD:
#        if the detected doublet fraction exceeds DOUBLET_ROLLBACK_THRESHOLD,
#        doublet removal is skipped for that sample to avoid over-filtering.
#     6. Ambient RNA correction via DecontX (optional)
#     7. Data normalization, dimensionality reduction, and batch-correction
#        via Harmony, producing parallel (PCA vs. Harmony) clusterings.
#     8. Saving a fully processed, un-annotated Seurat object for Script 02.
#
# CHECKPOINT SYSTEM:
#   Each sample is processed individually and saved to a checkpoint .rds file
#   inside SCEVAN_DIR/<SampleID>/. If a checkpoint exists on re-run, the
#   sample is loaded directly, saving hours of computation.
#
# MEMORY MANAGEMENT:
#   After processing each sample, large objects are explicitly removed and
#   garbage collection is forced (gc()) before the next iteration.
#
# OUTPUT:
#   - <PROJECT_NAME>_processed_for_annotation.rds  (main output for Script 02)
#   - Diagnostic QC violin plots (before/after filtering)
#   - Doublet pK-selection plots per sample
#   - Doublet UMAP visualization (before removal)
#   - Harmony vs. No-Harmony UMAP comparison plots
#
# NEXT STEP:
#   Open and run '02_annotate_data.R' for interactive annotation.
# =============================================================================

# --- Load Required Libraries -------------------------------------------------
# These must be installed via 00_rlibs_installation.R before running this script.
library(Seurat)       # Core single-cell analysis framework
library(harmony)      # Batch correction via Harmony integration
library(openxlsx)     # Reading .xlsx metadata files
library(dplyr)        # Data manipulation
library(ggplot2)      # Publication-quality plotting
library(patchwork)    # Combining multiple ggplot objects
library(celda)        # DecontX for ambient RNA correction
library(DoubletFinder)# Doublet detection
library(writexl)      # Writing .xlsx outputs
library(Matrix)       # Sparse matrix support
library(SCEVAN)       # Copy-number variation analysis

set.seed(123) # Set global random seed for full reproducibility

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) ---
# =============================================================================

# --- 1.1: Project Identity & Paths ---
# PROJECT_NAME: A short alphanumeric tag for this study. Used to name output files.
PROJECT_NAME <- "Nr4a1_Study17_Project"

# ROOT_PATH: The top-level working directory for this project.
# All other paths below are derived from this one.
#ROOT_PATH <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"   # Windows (RStudio local)
ROOT_PATH <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"  # Linux/HPC

# METADATA_FILE: Excel (.xlsx) file describing your samples.
#
# =============================================================================
# HOW TO CONFIGURE THE METADATA EXCEL FILE
# =============================================================================
# The metadata file is the single source of truth for all sample-level
# biological and technical information. Every column you add is automatically
# transferred to every cell from that sample — no extra code needed.
#
# REQUIRED COLUMN:
#   SampleID  — Must exactly match the subfolder name inside H5_DIR.
#               This is how the pipeline finds each sample's H5 files.
#               Example: if your H5 is at h5_files/S01_WT_M/, then SampleID = "S01_WT_M"
#
# COMMON OPTIONAL COLUMNS (add any combination you need):
#   Condition    — Primary experimental group (e.g., "KO" / "WT")
#   Genotype     — Genotype label (e.g., "Nr4a1_KO" / "WT")
#   Diet         — Dietary treatment (e.g., "inulin" / "cellulose")
#   Sex          — Biological sex ("M" / "F")
#   Timepoint    — Collection time ("Week4" / "Week8")
#   Treatment    — Drug/intervention label
#   BatchID      — Technical sequencing batch (use if you have batch effects
#                  beyond sample-to-sample variation; Harmony corrects for SampleID
#                  but BatchID can be added to group.by.vars if needed)
#   Age          — Animal age at collection
#
# AUTO-CREATED COLUMNS:
#   group        — Automatically created as Condition_Sex if both columns exist.
#                  You can reference this as a grouping variable in Scripts 02-04.
#
# CONCATENATED GROUPS IN SCRIPTS 02-04:
#   In Script 02, ADDITIONAL_GROUPS_TO_PLOT accepts ANY metadata column name,
#   including compound ones like "Genotype_Diet". If "Genotype_Diet" is NOT
#   already a column in your Excel file, Script 02 will auto-create it by
#   concatenating the "Genotype" and "Diet" columns with an underscore.
#   This works for any two-column combination using the pattern "ColA_ColB".
#   Example: "Condition_Sex" → auto-created from Condition + Sex columns.
#
# EXAMPLE METADATA TABLE:
#   SampleID     | Condition | Genotype   | Diet      | Sex | Timepoint
#   S01_WT_M     | WT        | WT         | cellulose | M   | Week4
#   S02_WT_F     | WT        | WT         | cellulose | F   | Week4
#   S03_KO_M     | KO        | Nr4a1_KO   | inulin    | M   | Week4
#   S04_KO_F     | KO        | Nr4a1_KO   | inulin    | F   | Week4
#
# After Script 01 runs, every cell from S03_KO_M will have:
#   $SampleID = "S03_KO_M",  $Condition = "KO",  $Genotype = "Nr4a1_KO",
#   $Diet = "inulin",  $Sex = "M",  $Timepoint = "Week4",  $group = "KO_M"
#
# In Script 02, you can then use:
#   ADDITIONAL_GROUPS_TO_PLOT <- c("Condition", "Genotype_Diet", "Sex")
# In Scripts 03-04:
#   CONDITION_COLUMN <- "Condition"   # or "Genotype_Diet", "Diet", etc.
#   GROUP_1 <- "KO"
#   GROUP_2 <- "WT"
# =============================================================================
METADATA_FILE <- file.path(ROOT_PATH, "Nr4a1_s17_metadata.xlsx")

# H5_DIR: Directory containing one subfolder per sample, each holding
#   the 10x H5 output files. Expected structure:
#   H5_DIR/<SampleID>/sample_filtered_feature_bc_matrix.h5
#   H5_DIR/<SampleID>/sample_raw_probe_bc_matrix.h5  (optional, for probe data)
H5_DIR        <- file.path(ROOT_PATH, "h5_files")

# OUTPUT_DIR: All processed objects and plots are saved here.
OUTPUT_DIR    <- file.path(ROOT_PATH, "seurat_output")

# SCEVAN_DIR: Dedicated directory for SCEVAN outputs and per-sample checkpoints.
SCEVAN_DIR    <- file.path(OUTPUT_DIR, "scevan_per_sample_results")

# --- 1.2: Pre-Merge QC Parameters (Loose / Per-Sample) ---
# Applied BEFORE merging samples. These are intentionally lenient to catch only
# clearly empty droplets and dying cells. Fine-tuning happens POST-merge.
PRE_MIN_GENES_PER_CELL <- 500   # Minimum genes detected per cell (discard empty droplets)
PRE_MAX_MT_PERCENT     <- 30.0  # Maximum mitochondrial % (discard dying cells, loose cutoff)

# --- 1.3: Post-Merge QC Parameters (Stringent) ---
# Applied AFTER merging. Examine the QC violin plots (01a/01b) to adjust these.
# The post-merge step also re-applies filtering after DecontX correction.
POST_MIN_GENES          <- 500    # Min genes per cell. Raise if low-quality cells persist.
POST_MAX_GENES          <- 14000  # Max genes per cell. High values may indicate multiplets.
POST_MIN_UMIS           <- 1500   # Min UMI count. Raise to remove low-depth cells.
POST_MAX_UMIS           <- 100000 # Max UMI count. Rarely needs changing.
POST_MAX_MT             <- 20.0   # Max mitochondrial %. Standard range: 15-25%.
POST_MIN_CELLS_PER_GENE <- 15     # Genes expressed in fewer cells than this are removed.

# --- 1.4: DoubletFinder Parameters ---
# DOUBLET_RATE: Expected fraction of doublets. Use the 10x rule of thumb:
#   ~0.8% per 1,000 cells loaded (e.g., 10,000 cells loaded → 0.08).
DOUBLET_RATE <- 0.08

# pK Selection: DoubletFinder sweeps a range of pK values and picks the one
# maximizing the BCmetric. The range below constrains the search to biologically
# sensible values. The fallback is used if no valid pK is found in range.
DF_PK_RANGE_MIN <- 0.01  # Minimum valid pK value
DF_PK_RANGE_MAX <- 0.15  # Maximum valid pK value
DF_PK_FALLBACK  <- 0.09  # Fallback pK if no valid value found in range

# DOUBLET_ROLLBACK_THRESHOLD: Safety mechanism.
#   If DoubletFinder classifies MORE than this fraction of cells as doublets
#   for a given sample, the result is suspicious (often a sign of bad pK
#   selection or unusual data). In that case, doublet removal is ROLLED BACK
#   for that sample (all cells kept), and a warning is issued.
#   Rationale: a >14% doublet rate is biologically implausible for standard
#   10x Genomics runs and most likely reflects a classification artifact.
DOUBLET_ROLLBACK_THRESHOLD <- 0.14  # 14% maximum acceptable doublet fraction

# --- 1.5: Core Dimensionality Reduction & Clustering Parameters ---
N_VARIABLE_FEATURES <- 5000  # Number of highly variable genes for PCA.
                              # Range: 2000-5000. Fewer HVGs can improve UMAP
                              # separation in heterogeneous datasets.
N_PCS_TO_USE        <- 50    # Number of PCs used for graph construction and UMAP.
                              # Increase for very complex datasets.
CLUSTER_RESOLUTION  <- 1.0   # Leiden/Louvain resolution. Higher = more clusters.
                              # Start at 1.0; tune after viewing the UMAP in Script 02.
UMAP_N_NEIGHBORS    <- 15    # UMAP: local neighborhood size. Higher = more global structure.
UMAP_MIN_DIST       <- 0.3   # UMAP: minimum distance between embedded points.
                              # Lower = tighter clusters; higher = more spread.

# --- 1.6: Workflow Toggles ---
RUN_DECONTX <- TRUE   # Run DecontX ambient RNA correction (recommended for all runs).
RUN_SCEVAN  <- TRUE   # Run SCEVAN copy-number variation analysis per sample.
DPI_SETTING <- 300    # DPI for all saved diagnostic plots.

# --- 1.7: Probe / KO Gene Integration Parameters ---
# Enable this if your 10x run included a probe-based capture assay (e.g., for
# tracking KO efficiency at the probe level). Set to FALSE to skip entirely.
ADD_PROBE_DATA <- TRUE

# PROBE_MAPPING: Named list mapping probe IDs (as they appear in the H5 file)
# to human-readable labels. These labels become metadata column names.
PROBE_MAPPING <- list(
  'Nr4a1|dbf6af3' = 'Exon 3 (KO Target)',
  'Nr4a1|99a0eaa' = 'Exon 5 (KO Target)',
  'Nr4a1|03c248a' = 'Exon 7'
)

# PROBES_FOR_CUSTOM_SUM: A subset of probe IDs whose counts are summed into
# a single custom feature (e.g., total KO-exon coverage). This custom feature
# is added as a new row ('Nr4a1_cust') to the count matrix.
PROBES_FOR_CUSTOM_SUM <- c('Nr4a1|dbf6af3', 'Nr4a1|99a0eaa')

# --- 1.8: SCEVAN Parameters ---
SCEVAN_ORGANISM  <- "mouse"  # "mouse" or "human"
SCEVAN_N_CORES   <- 4        # Parallel cores for SCEVAN. Increase on HPC.
SCEVAN_SUBCLONES <- TRUE     # Detect subclonal CNA populations.
SCEVAN_PLOTTREE  <- FALSE    # Plot phylogenetic tree of subclones (slow; set TRUE if needed).

# =============================================================================
# --- PART 2: PIPELINE EXECUTION (DO NOT EDIT BELOW THIS LINE) ---
# =============================================================================

# --- Create output directories if they do not exist -------------------------
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR, recursive = TRUE) }
if (!dir.exists(SCEVAN_DIR)) { dir.create(SCEVAN_DIR, recursive = TRUE) }

# =============================================================================
# --- STEP 2.1: Per-Sample Processing with Checkpoint System -----------------
# =============================================================================
# Each sample is processed individually. A checkpoint .rds file is saved after
# processing. On re-run, existing checkpoints are loaded directly to save time.
message("=== STEP 2.1: Processing or Loading Individual Samples ===")
metadata           <- read.xlsx(METADATA_FILE)
seurat_objects_list <- list()

for (i in 1:nrow(metadata)) {
  sample_info      <- metadata[i, ]
  sample_id        <- as.character(sample_info$SampleID)
  sample_scevan_dir <- file.path(SCEVAN_DIR, sample_id)
  checkpoint_file  <- file.path(sample_scevan_dir, paste0(sample_id, "_scevan_processed.rds"))

  # --- Load from checkpoint if available ---
  if (file.exists(checkpoint_file)) {
    message(paste("  [CHECKPOINT] Found checkpoint for", sample_id, "- loading from disk."))
    seurat_obj <- readRDS(checkpoint_file)

  } else {
    # --- Process sample from raw H5 ---
    message(paste("  [PROCESSING] No checkpoint for", sample_id, "- processing from scratch."))
    if (!dir.exists(sample_scevan_dir)) { dir.create(sample_scevan_dir, recursive = TRUE) }

    # -- Read RNA count matrix from 10x H5 file --
    rna_h5_path <- file.path(H5_DIR, sample_id, "sample_filtered_feature_bc_matrix.h5")
    if (!file.exists(rna_h5_path)) {
      warning(paste("  [SKIP] RNA H5 file not found for:", sample_id))
      next
    }
    counts_matrix <- Read10X_h5(rna_h5_path)

    # -- Optional: Read probe data and create custom summed feature --
    probe_matrix <- NULL
    if (ADD_PROBE_DATA) {
      probe_h5_path <- file.path(H5_DIR, sample_id, "sample_raw_probe_bc_matrix.h5")
      if (file.exists(probe_h5_path)) {
        message("    -> Reading probe data and creating custom 'Nr4a1_cust' feature...")
        probe_matrix <- Read10X_h5(probe_h5_path)
        # Handle nested list structure from probe H5 files
        if (is.list(probe_matrix) && "Probe Barcode" %in% names(probe_matrix)) {
          probe_matrix <- probe_matrix[["Probe Barcode"]]
        }
        # Sum counts across selected KO-target probes into one custom feature
        common_cells       <- intersect(colnames(counts_matrix), colnames(probe_matrix))
        custom_sum_counts  <- numeric(length(common_cells))
        names(custom_sum_counts) <- common_cells
        probes_found <- intersect(PROBES_FOR_CUSTOM_SUM, rownames(probe_matrix))
        if (length(probes_found) > 0) {
          summed_vector <- Matrix::colSums(probe_matrix[probes_found, common_cells, drop = FALSE])
          custom_sum_counts[names(summed_vector)] <- summed_vector
        }
        # Append custom feature row to the RNA count matrix
        custom_row <- Matrix::Matrix(0, 1, ncol(counts_matrix), sparse = TRUE)
        colnames(custom_row) <- colnames(counts_matrix)
        rownames(custom_row) <- "Nr4a1_cust"
        custom_row[1, names(custom_sum_counts)] <- custom_sum_counts
        counts_matrix <- rbind(counts_matrix, custom_row)
      } else {
        warning(paste("    -> [SKIP] Probe H5 not found for", sample_id))
      }
    }

    # -- Create Seurat object --
    seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = sample_id, min.cells = 5)

    # -- Attach all metadata columns from the Excel file --
    for (col_name in colnames(sample_info)) {
      seurat_obj[[col_name]] <- sample_info[[col_name]]
    }
    # Auto-create combined group label (Condition_Sex) if both columns are present
    if ("Condition" %in% colnames(seurat_obj@meta.data) &&
        "Sex" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj$group <- paste(seurat_obj$Condition, seurat_obj$Sex, sep = "_")
    }

    # -- Add individual probe counts as metadata columns --
    if (!is.null(probe_matrix)) {
      for (probe_id in names(PROBE_MAPPING)) {
        if (probe_id %in% rownames(probe_matrix)) {
          # Sanitize probe label to create a valid R column name
          col_name_clean <- gsub("[\\s\\(\\)\n]+", "_", PROBE_MAPPING[[probe_id]])
          col_name_clean <- gsub("_$", "", col_name_clean)
          col_name_final <- paste0("Nr4a1_", col_name_clean, "_probe_count")
          seurat_obj <- AddMetaData(seurat_obj,
                                    metadata = probe_matrix[probe_id, ],
                                    col.name = col_name_final)
        }
      }
    }

    # -- Optional: Run SCEVAN for copy-number variation detection --
    if (RUN_SCEVAN) {
      tryCatch({
        message("    -> Running SCEVAN for CNA detection...")
        scevan_res_df <- pipelineCNA(
          as.matrix(GetAssayData(seurat_obj, layer = "counts")),
          par_cores  = SCEVAN_N_CORES,
          SUBCLONES  = SCEVAN_SUBCLONES,
          plotTree   = SCEVAN_PLOTTREE,
          organism   = SCEVAN_ORGANISM,
          output_dir = sample_scevan_dir
        )
        if (!is.null(scevan_res_df)) {
          seurat_obj <- AddMetaData(seurat_obj, metadata = scevan_res_df)
          message("    -> SCEVAN metadata added to Seurat object.")
        }
      }, error = function(e) {
        message(paste("    -> [WARNING] SCEVAN failed for", sample_id, "| Error:", e$message))
      })
    }

    # -- Pre-merge QC: compute mitochondrial percentage and apply loose filters --
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
    seurat_obj <- subset(seurat_obj,
                         subset = nFeature_RNA >= PRE_MIN_GENES_PER_CELL &
                                  percent.mt   <= PRE_MAX_MT_PERCENT)

    # -- Save checkpoint for this sample --
    message("    -> Saving processed object to checkpoint file.")
    saveRDS(seurat_obj, file = checkpoint_file)

    # -- Explicit memory cleanup before next iteration --
    message("    -> Freeing memory before next sample...")
    rm(counts_matrix, seurat_obj)
    if (!is.null(probe_matrix)) { rm(probe_matrix) }
    gc()

    # Reload the object cleanly to add to the main list
    seurat_obj <- readRDS(checkpoint_file)
  }

  seurat_objects_list[[sample_id]] <- seurat_obj
}

# =============================================================================
# --- STEP 2.2: Merge All Processed Samples ----------------------------------
# =============================================================================
if (length(seurat_objects_list) == 0) stop("[ERROR] No samples were loaded. Check H5 file paths.")
message("\n=== STEP 2.2: Merging all processed samples ===")
data <- merge(
  x = seurat_objects_list[[1]],
  y = seurat_objects_list[-1],
  add.cell.ids = names(seurat_objects_list)
)
data <- JoinLayers(data)   # Required for Seurat v5: joins split count layers into one
rm(seurat_objects_list); gc()

# =============================================================================
# --- STEP 2.3: Post-Merge QC & Diagnostic Plots ----------------------------
# =============================================================================
message("\n=== STEP 2.3: Post-Merge QC and Diagnostic Plotting ===")

# Save violin plots BEFORE filtering so the full distribution is visible
p_before <- VlnPlot(data,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "01a_qc_violin_before_filtering.png"),
       plot = p_before, width = 10, height = 6, dpi = DPI_SETTING)

# Apply stringent post-merge filters (cell-level)
data <- subset(data,
               subset = nFeature_RNA >= POST_MIN_GENES  &
                        nFeature_RNA <= POST_MAX_GENES  &
                        nCount_RNA   >= POST_MIN_UMIS   &
                        nCount_RNA   <= POST_MAX_UMIS   &
                        percent.mt   <= POST_MAX_MT)

# Apply gene-level filter: remove genes expressed in too few cells
genes_to_keep <- rownames(data)[
  Matrix::rowSums(GetAssayData(data, layer = "counts") > 0) >= POST_MIN_CELLS_PER_GENE
]
data <- subset(data, features = genes_to_keep)

p_after <- VlnPlot(data,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "01b_qc_violin_after_filtering.png"),
       plot = p_after, width = 10, height = 6, dpi = DPI_SETTING)

message(paste("  Cells after QC:", ncol(data), "| Genes retained:", nrow(data)))

# =============================================================================
# --- STEP 2.4: Doublet Detection via DoubletFinder --------------------------
# =============================================================================
# Strategy: DoubletFinder is run per-sample (not on the merged object) to
# avoid confounding batch effects. The per-cell doublet labels are then
# concatenated and applied back to the merged object.
#
# ROLLBACK SAFEGUARD:
#   If a sample's detected doublet fraction exceeds DOUBLET_ROLLBACK_THRESHOLD,
#   it is likely a classification artifact. All cells in that sample are
#   retained (treated as singlets), and a warning is logged.
message("\n=== STEP 2.4: Doublet Detection (DoubletFinder) ===")
pk_plot_dir <- file.path(OUTPUT_DIR, "doublet_finder_plots")
if (!dir.exists(pk_plot_dir)) dir.create(pk_plot_dir)

data_list    <- SplitObject(data, split.by = "SampleID")
results_list <- list()
rollback_log <- list()  # Tracks any samples where rollback was triggered

for (i in seq_along(data_list)) {
  sample_name <- names(data_list)[i]
  message(paste("  Processing doublets for:", sample_name))
  seu_tmp <- data_list[[i]]

  # --- Preliminary processing required by DoubletFinder ---
  seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = N_PCS_TO_USE, verbose = FALSE)

  # --- pK selection: find the pK that maximizes the BCmetric ---
  sweep.res.list <- paramSweep(seu_tmp, PCs = 1:N_PCS_TO_USE, sct = FALSE)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)
  bcmvn$pK       <- as.numeric(as.character(bcmvn$pK))
  initial_pk     <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  final_pk       <- initial_pk

  # --- Constrain pK to biologically sensible range ---
  # If the global maximum falls outside [DF_PK_RANGE_MIN, DF_PK_RANGE_MAX],
  # find the best pK within range. If no valid value exists, use the fallback.
  if (is.na(final_pk) || final_pk < DF_PK_RANGE_MIN || final_pk > DF_PK_RANGE_MAX) {
    bcmvn_filtered <- bcmvn[bcmvn$pK > DF_PK_RANGE_MIN & bcmvn$pK < DF_PK_RANGE_MAX, ]
    if (nrow(bcmvn_filtered) > 0 && any(is.finite(bcmvn_filtered$BCmetric))) {
      final_pk <- bcmvn_filtered$pK[which.max(bcmvn_filtered$BCmetric)]
    } else {
      final_pk <- DF_PK_FALLBACK
      message(paste("    -> [WARNING] No valid pK in range for", sample_name,
                    "- using fallback pK =", DF_PK_FALLBACK))
    }
  }

  # --- Save pK selection diagnostic plot ---
  pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
    ggtitle(paste0("pK Finder: ", sample_name),
            subtitle = paste0("Selected pK = ", final_pk,
                              " | Initial optimal = ", initial_pk)) +
    theme_minimal()
  ggsave(file.path(pk_plot_dir, paste0(sample_name, "_pk_plot.png")),
         plot = pk_plot, width = 7, height = 5, dpi = DPI_SETTING)

  # --- Run DoubletFinder ---
  nExp_val <- round(ncol(seu_tmp) * DOUBLET_RATE)
  seu_tmp  <- doubletFinder(seu_tmp, PCs = 1:N_PCS_TO_USE,
                             pK = final_pk, nExp = nExp_val, sct = FALSE)
  res_col  <- tail(grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE), 1)

  # --- ROLLBACK CHECK ---
  # Count the fraction of cells classified as doublets for this sample
  n_doublets       <- sum(seu_tmp@meta.data[[res_col]] == "Doublet")
  doublet_fraction <- n_doublets / ncol(seu_tmp)

  if (doublet_fraction > DOUBLET_ROLLBACK_THRESHOLD) {
    # Rollback: retain all cells for this sample; log the event
    warning(paste0(
      "  [ROLLBACK] Sample '", sample_name, "': detected doublet fraction = ",
      round(doublet_fraction * 100, 1), "% exceeds threshold of ",
      round(DOUBLET_ROLLBACK_THRESHOLD * 100, 1), "%. ",
      "Doublet removal SKIPPED for this sample (all cells treated as Singlets)."
    ))
    rollback_result <- data.frame(
      Doublet_Status = rep("Singlet", ncol(seu_tmp)),
      row.names      = colnames(seu_tmp)
    )
    results_list[[sample_name]] <- rollback_result
    rollback_log[[sample_name]] <- list(
      fraction     = doublet_fraction,
      n_doublets   = n_doublets,
      n_cells      = ncol(seu_tmp),
      final_pk     = final_pk,
      action       = "ROLLBACK_APPLIED"
    )
  } else {
    # Normal path: use DoubletFinder classifications
    results_list[[sample_name]] <- seu_tmp@meta.data[, res_col, drop = FALSE]
    rollback_log[[sample_name]] <- list(
      fraction     = doublet_fraction,
      n_doublets   = n_doublets,
      n_cells      = ncol(seu_tmp),
      final_pk     = final_pk,
      action       = "DOUBLETS_REMOVED"
    )
  }
  message(paste0("    -> Doublet fraction: ",
                 round(doublet_fraction * 100, 1), "% | pK used: ", final_pk,
                 " | Action: ", rollback_log[[sample_name]]$action))
}

# --- Save rollback log as a summary Excel file ---
rollback_df <- do.call(rbind, lapply(names(rollback_log), function(s) {
  data.frame(
    SampleID         = s,
    N_Cells          = rollback_log[[s]]$n_cells,
    N_Doublets_Flagged = rollback_log[[s]]$n_doublets,
    Doublet_Fraction = round(rollback_log[[s]]$fraction * 100, 2),
    Threshold_Pct    = DOUBLET_ROLLBACK_THRESHOLD * 100,
    pK_Used          = rollback_log[[s]]$final_pk,
    Action           = rollback_log[[s]]$action,
    stringsAsFactors = FALSE
  )
}))
write_xlsx(rollback_df, file.path(OUTPUT_DIR, "doublet_finder_rollback_log.xlsx"))
message(paste("  Doublet summary saved to: doublet_finder_rollback_log.xlsx"))

# --- Merge per-sample doublet results back to the combined object ---
all_res       <- do.call(rbind, lapply(results_list, function(df) setNames(df, "Doublet_Status")))
data$Doublet_Status <- all_res[rownames(data@meta.data), "Doublet_Status"]

# --- Diagnostic UMAP of doublets BEFORE removal ---
message("  -> Generating doublet visualization UMAP (diagnostic only)...")
data <- NormalizeData(data, verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = N_PCS_TO_USE, verbose = FALSE, reduction.name = "pca_temp_doublets") %>%
  RunUMAP(dims = 1:N_PCS_TO_USE, reduction = "pca_temp_doublets",
          reduction.name = "umap_temp_doublets", verbose = FALSE)

p_doublets <- DimPlot(data, reduction = "umap_temp_doublets",
                      group.by = "Doublet_Status",
                      cols = c("Singlet" = "grey80", "Doublet" = "red")) +
  ggtitle("Doublet Visualization (Before Removal)")
ggsave(file.path(OUTPUT_DIR, "01c_DIAGNOSTIC_doublet_visualization.png"),
       plot = p_doublets, width = 8, height = 7, dpi = DPI_SETTING)

# Remove temporary reductions to keep the object clean
data@reductions$pca_temp_doublets  <- NULL
data@reductions$umap_temp_doublets <- NULL
rm(p_doublets); gc()

# --- Remove doublets ---
message(paste("  Total cells BEFORE doublet filtering:", ncol(data)))
data <- subset(data, subset = Doublet_Status != "Doublet" | is.na(Doublet_Status))
message(paste("  Total cells AFTER doublet filtering:", ncol(data)))

# =============================================================================
# --- STEP 2.5: Ambient RNA Correction (DecontX) ------------------------------
# =============================================================================
# DecontX models ambient RNA contamination per cell and produces a
# decontaminated count matrix. The corrected matrix replaces the raw counts,
# and QC metrics are recalculated and re-filtered.
if (RUN_DECONTX) {
  message("\n=== STEP 2.5: Ambient RNA Correction (DecontX) ===")
  counts_sparse    <- GetAssayData(object = data, layer = "counts")
  decontx_results  <- decontX(x = counts_sparse)

  # Replace raw counts with decontaminated counts
  data[["RNA"]]$counts <- decontx_results$decontXcounts
  data[["RNA"]]$data   <- NULL  # Clear any stale normalized layer

  # Recompute QC metrics on decontaminated matrix
  data$nCount_RNA   <- colSums(data[["RNA"]]$counts)
  data$nFeature_RNA <- colSums(data[["RNA"]]$counts > 0)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^mt-")

  # Re-apply stringent QC filters after decontamination
  data <- subset(data,
                 subset = nFeature_RNA >= POST_MIN_GENES &
                          nFeature_RNA <= POST_MAX_GENES &
                          nCount_RNA   >= POST_MIN_UMIS  &
                          nCount_RNA   <= POST_MAX_UMIS  &
                          percent.mt   <= POST_MAX_MT)
  message(paste("  Cells remaining after DecontX re-filtering:", ncol(data)))
}

# =============================================================================
# --- STEP 2.6: Normalization, Dimensionality Reduction & Integration --------
# =============================================================================
# Two parallel analytical tracks are produced:
#   (A) Standard PCA - no batch correction (useful to diagnose batch effects)
#   (B) Harmony-corrected PCA (recommended for downstream annotation)
# Both UMAPs and cluster assignments are saved in the Seurat object.
message("\n=== STEP 2.6: Normalization, Dimensionality Reduction, and Integration ===")

# --- Normalize, find HVGs, scale, and run PCA ---
data <- NormalizeData(data, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = N_VARIABLE_FEATURES, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = N_PCS_TO_USE, verbose = FALSE)

# --- Track A: Standard PCA (no batch correction) ---
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                      graph.name = "pca_nn", verbose = FALSE)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION,
                     graph.name = "pca_nn", cluster.name = "clusters_none",
                     verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                reduction.name = "umap_none", verbose = FALSE)

# --- Track B: Harmony integration (batch correction by SampleID) ---
data <- RunHarmony(data, group.by.vars = "SampleID",
                   reduction = "pca", reduction.save = "harmony",
                   verbose = FALSE)
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                      graph.name = "harmony_nn", verbose = FALSE)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION,
                     graph.name = "harmony_nn", cluster.name = "clusters_harmony",
                     verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                reduction.name = "umap_harmony", verbose = FALSE)

# --- Diagnostic plots: compare Harmony vs. no Harmony ---
p1 <- DimPlot(data, reduction = "umap_none",    group.by = "SampleID") +
  ggtitle("Standard PCA (No Harmony)")
p2 <- DimPlot(data, reduction = "umap_harmony", group.by = "SampleID") +
  ggtitle("Harmony Batch-Corrected")
ggsave(file.path(OUTPUT_DIR, "02_DIAGNOSTIC_UMAP_Harmony_vs_NoHarmony.png"),
       plot = p1 + p2, width = 16, height = 7, dpi = DPI_SETTING)

p3 <- DimPlot(data, reduction = "umap_none",    group.by = "clusters_none",
              label = TRUE) + ggtitle("Clusters: Standard PCA") + NoLegend()
p4 <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony",
              label = TRUE) + ggtitle("Clusters: Harmony") + NoLegend()
ggsave(file.path(OUTPUT_DIR, "03_DIAGNOSTIC_UMAP_Cluster_Comparison.png"),
       plot = p3 + p4, width = 16, height = 7, dpi = DPI_SETTING)

message("  Saved Harmony vs. PCA diagnostic UMAP plots.")

# =============================================================================
# --- STEP 2.7: Save Processed Object ----------------------------------------
# =============================================================================
message("\n=== STEP 2.7: Saving Processed Object for Annotation ===")
output_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
saveRDS(data, output_rds)
message(paste0(
  "\n",
  "=== PROCESSING COMPLETE ===\n",
  "  Output saved: ", output_rds, "\n",
  "  Total cells:  ", ncol(data), "\n",
  "  Total genes:  ", nrow(data), "\n",
  "\nNext step: Open and run '02_annotate_data.R' to annotate the data.\n"
))
