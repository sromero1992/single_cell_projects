# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING ENGINE
# Version: 9.0 (Generalized Pipeline)
#
# PURPOSE:
#   This script is the computational backbone of the single-cell RNA-seq
#   pipeline. It sequentially performs:
#     1. Per-sample data loading from 10x Genomics H5 files
#     2. Optional probe/custom feature integration (e.g., KO-target probes)
#     3. Copy-number variation detection via SCEVAN (optional, per-sample)
#     4. Pre-merge QC (loose filters, per-sample)
#     5. Ambient RNA correction via DecontX (optional, per-sample),
#        with optional rounding of corrected counts to integers.
#     6. Doublet detection + removal via DoubletFinder (per-sample),
#        with a ROLLBACK SAFEGUARD: if the detected doublet fraction exceeds
#        DOUBLET_ROLLBACK_THRESHOLD, removal is skipped for that sample.
#     7. Post-merge QC (stringent filters applied after merging all samples)
#     8. Data normalization, dimensionality reduction, and batch-correction
#        via Harmony, producing parallel (PCA vs. Harmony) clusterings.
#     9. Saving a fully processed, un-annotated Seurat object for Script 02.
#
# PER-SAMPLE PROCESSING ORDER (STEP 2.1a):
#   Load → [SCEVAN] → [DecontX] → DoubletFinder → Remove Doublets → Pre-QC
#   → Checkpoint
#   This order ensures that:
#     - SCEVAN sees unmodified raw counts (best for CNA detection)
#     - DoubletFinder operates on ambient-corrected counts
#     - Pre-merge QC metrics reflect the corrected, singlet-only matrix
#     - Each checkpoint contains clean, singlet-only cells
#
# CHECKPOINT SYSTEM (dual-checkpoint per sample):
#   Each sample produces two checkpoint .rds files inside SCEVAN_DIR/<SampleID>/:
#     _scevan_processed.rds       — saved right after SCEVAN (Checkpoint 1)
#     _decontx_dblt_processed.rds — saved after DecontX + DoubletFinder + Pre-QC
#                                   (Checkpoint 2 — the file used for merging)
#   On re-run, Checkpoint 2 is tried first (sample fully done → skip). If only
#   Checkpoint 1 exists, SCEVAN is skipped and the pipeline resumes from
#   DecontX. If neither exists, the full pipeline runs from scratch.
#
# MEMORY MANAGEMENT:
#   After processing each sample, large objects are explicitly removed and
#   garbage collection is forced (gc()) before the next iteration.
#   BPCells on-disk matrices are used to handle large datasets.
#
# OUTPUT:
#   - <PROJECT_NAME>_processed_for_annotation.rds  (main output for Script 02)
#   - Diagnostic QC violin plots (before/after filtering)
#   - Doublet pK-selection plots per sample (in SCEVAN_DIR/<SampleID>/)
#   - Doublet UMAP visualization per sample (in SCEVAN_DIR/<SampleID>/)
#   - Doublet rollback summary Excel file
#   - Harmony vs. No-Harmony UMAP comparison plots
#
# NEXT STEP:
#   Open and run '02_annotate_data.R' for interactive annotation.
# =============================================================================

# --- Load Required Libraries -------------------------------------------------
# These must be installed via 00_rlibs_installation.R before running this script.
library(Seurat)         # Core single-cell analysis framework
library(SeuratWrappers) # Batch correction via Harmony integration
library(openxlsx)       # Reading .xlsx metadata files
library(dplyr)          # Data manipulation
library(ggplot2)        # Publication-quality plotting
library(patchwork)      # Combining multiple ggplot objects
library(celda)          # DecontX for ambient RNA correction
library(DoubletFinder)  # Doublet detection
library(writexl)        # Writing .xlsx outputs
library(Matrix)         # Sparse matrix support
library(SCEVAN)         # Copy-number variation analysis
library(hdf5r)          # HDF5 file support
library(BPCells)        # On-disk matrix handling for large datasets

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
#   Example: "Condition_Sex" -> auto-created from Condition + Sex columns.
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
PRE_MAX_MT_PERCENT     <- 20.0  # Maximum mitochondrial % (discard dying cells, loose cutoff)

# --- 1.3: Post-Merge QC Parameters (Stringent) ---
# Applied AFTER merging. Examine the QC violin plots (01a/01b) to adjust these.
POST_MIN_GENES          <- 500    # Min genes per cell. Raise if low-quality cells persist.
POST_MAX_GENES          <- 15000  # Max genes per cell. High values may indicate multiplets.
POST_MIN_UMIS           <- 1500   # Min UMI count. Raise to remove low-depth cells.
POST_MAX_UMIS           <- 100000 # Max UMI count. Similarly look at QC violins to adjust.
POST_MAX_MT             <- 5.0   # Max mitochondrial %. Standard range: 5-20%.
POST_MIN_CELLS_PER_GENE <- 15     # Genes expressed in fewer cells than this are removed.

# --- 1.4: DoubletFinder Parameters ---
# DOUBLET_RATE: Expected fraction of doublets. Use the 10x rule of thumb:
#   ~0.8% per 1,000 cells loaded (e.g., 10,000 cells loaded -> 0.08).
DOUBLET_RATE <- 0.08

# pK Selection: DoubletFinder sweeps a range of pK values and picks the one
# maximizing the BCmetric. The range below constrains the search to biologically
# sensible values. The fallback is used if no valid pK is found in range.
DF_PK_RANGE_MIN <- 0.005  # Minimum valid pK value
DF_PK_RANGE_MAX <- 0.25   # Maximum valid pK value
DF_PK_FALLBACK  <- 0.09   # Fallback pK if no valid value found in range

# DOUBLET_ROLLBACK_THRESHOLD: Safety mechanism.
#   If DoubletFinder classifies MORE than this fraction of cells as doublets
#   for a given sample, the result is suspicious (often a sign of bad pK
#   selection or unusual data). In that case, doublet removal is ROLLED BACK
#   for that sample (all cells kept), and a warning is issued.
#   Rationale: a >14% doublet rate is biologically implausible for standard
#   10x Genomics runs and most likely reflects a classification artifact.
DOUBLET_ROLLBACK_THRESHOLD <- DF_PK_RANGE_MAX  # 14% maximum acceptable doublet fraction

# --- 1.5: Core Dimensionality Reduction & Clustering Parameters ---
N_VARIABLE_FEATURES <- 2000  # Number of highly variable genes for PCA.
                              # Range: 2000-5000. Fewer HVGs can improve UMAP
                              # separation in heterogeneous datasets.
N_PCS_TO_USE        <- 50    # Number of PCs used for graph construction and UMAP.
                              # Increase for very complex datasets.
CLUSTER_RESOLUTION  <- 1.0   # Leiden/Louvain resolution. Higher = more clusters.
                              # Start at 1.0; tune after viewing the UMAP in Script 02.
UMAP_N_NEIGHBORS    <- 30    # UMAP: local neighborhood size. Higher = more global structure.
UMAP_MIN_DIST       <- 0.3   # UMAP: minimum distance between embedded points.
                              # Lower = tighter clusters; higher = more spread.

# --- 1.6: Workflow Toggles ---
RUN_DECONTX            <- TRUE  # Run DecontX ambient RNA correction (recommended for all runs).

ROUND_DECONTX_COUNTS   <- FALSE # Round decontX corrected counts to nearest integer.
                                 # DecontX produces non-integer (fractional) values by default.
                                 # Set TRUE if downstream tools require integer count matrices.
                                 # See decontX vignette for details.

USE_BPCELLS            <- TRUE  # Use BPCells on-disk matrix handling (write_matrix_dir) during
                                 # the post-merge QC (Step 2.3) and integration (Step 2.4) steps.
                                 # Recommended for large cohorts (>8 samples) to prevent RAM
                                 # exhaustion. Set FALSE for small datasets or if the BPCells
                                 # package is not available; standard in-memory JoinLayers() is
                                 # used as fallback. Requires hdf5r and BPCells to be installed.

RUN_SCEVAN             <- TRUE  # Run SCEVAN copy-number variation analysis.
                                 # Enable if ANY sample in the cohort is cancer/tumor — run on ALL
                                 # samples including normals, as normals serve as the CNA reference
                                 # baseline. Set FALSE only if the entire study uses normal tissue.
DPI_SETTING            <- 300   # DPI for all saved diagnostic plots.

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
SCEVAN_N_CORES   <- 6        # Parallel cores for SCEVAN. Increase on HPC.
SCEVAN_SUBCLONES <- TRUE     # Detect subclonal CNA populations.
SCEVAN_PLOTTREE  <- FALSE    # Plot phylogenetic tree of subclones (slow; set TRUE if needed).

# =============================================================================
# --- PART 2: PIPELINE EXECUTION (DO NOT EDIT BELOW THIS LINE) ---
# =============================================================================

# --- Create output directories if they do not exist -------------------------
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR, recursive = TRUE) }
if (!dir.exists(SCEVAN_DIR)) { dir.create(SCEVAN_DIR, recursive = TRUE) }

# =============================================================================
# --- STEP 2.1a: Per-Sample Processing & Checkpoint Generation ---------------
# =============================================================================
# PURPOSE:
#   Process each sample independently and save checkpoint .rds files.
#   Each sample goes through the full per-sample pipeline in this order:
#     1. Load H5 data and optionally integrate probe data
#     2. Create Seurat object and attach metadata
#     3. SCEVAN (if enabled) — run on raw counts before any correction
#     4. DecontX (if enabled) — ambient RNA correction, with optional rounding
#     5. DoubletFinder — per-sample doublet detection and removal
#     6. Pre-merge QC (loose filters) — applied on corrected, singlet-only counts
#     7. Save checkpoint and free memory
#
# DUAL-CHECKPOINT SYSTEM (three-tier resume logic):
#   Each sample produces TWO checkpoint files inside SCEVAN_DIR/<SampleID>/:
#
#     Checkpoint 1 (_scevan_processed.rds):
#       Saved after SCEVAN runs. Contains raw Seurat object + SCEVAN metadata.
#       Used as a mid-pipeline recovery point so that expensive SCEVAN
#       computation is never repeated.
#
#     Checkpoint 2 (_decontx_dblt_processed.rds):
#       Saved after DecontX + DoubletFinder + Pre-QC. This is the FINAL
#       per-sample checkpoint and the file loaded for the merge step (2.1b).
#       It contains clean, singlet-only, ambient-corrected cells.
#
#   At the start of each iteration the loop checks:
#     (a) If Checkpoint 2 exists  → sample fully processed, SKIP entirely.
#     (b) If only Checkpoint 1 exists → load it, SKIP SCEVAN, run DecontX
#         + DoubletFinder + Pre-QC, then save Checkpoint 2.
#     (c) If neither exists → run the entire pipeline from scratch.
#
#   This means that if a run crashes mid-way through DecontX or DoubletFinder,
#   the next run will resume from the SCEVAN result rather than starting over.
#
# MEMORY STRATEGY:
#   Each sample object is removed from RAM and garbage-collected immediately
#   after its checkpoint is written. Only one sample occupies memory at a time.
# =============================================================================
message("=== STEP 2.1a: Per-Sample Processing (SCEVAN → DecontX → DoubletFinder → Pre-QC → Checkpoint) ===")
metadata <- read.xlsx(METADATA_FILE)

# Initialize rollback log (filled during loop, saved after)
rollback_log <- list()

for (i in 1:nrow(metadata)) {
  sample_info       <- metadata[i, ]
  sample_id         <- as.character(sample_info$SampleID)
  sample_scevan_dir <- file.path(SCEVAN_DIR, sample_id)
  checkpoint_file_1 <- file.path(sample_scevan_dir, paste0(sample_id, "_scevan_processed.rds"))
  checkpoint_file_2 <- file.path(sample_scevan_dir, paste0(sample_id, "_decontx_dblt_processed.rds"))

  # --- Three-tier checkpoint logic -------------------------------------------
  # (a) Full checkpoint exists: sample is completely done, skip.
  if (file.exists(checkpoint_file_2)) {
    message(paste("  [CHECKPOINT-2] Skipping", sample_id,
                  "- full checkpoint (decontX + doublets + QC) found."))
    next
  }

  # (b) SCEVAN-only checkpoint exists: load it, skip SCEVAN, run remaining steps.
  scevan_checkpoint_loaded <- FALSE
  if (file.exists(checkpoint_file_1)) {
    message(paste("  [CHECKPOINT-1] Loading SCEVAN checkpoint for", sample_id,
                  "- will run DecontX + DoubletFinder + Pre-QC only."))
    seurat_obj <- readRDS(checkpoint_file_1)
    scevan_checkpoint_loaded <- TRUE
  }
  # (c) No checkpoint: fall through to full processing below.

  message(paste("  [PROCESSING] Starting", sample_id,
                if (scevan_checkpoint_loaded) "(resuming from SCEVAN checkpoint)" else "(full run)",
                "..."))
  if (!dir.exists(sample_scevan_dir)) { dir.create(sample_scevan_dir, recursive = TRUE) }

  # ---- Steps 1–4: Load data, build Seurat object, run SCEVAN ---------------
  # These steps are SKIPPED when resuming from a SCEVAN checkpoint (case b above).
  # The flag `scevan_checkpoint_loaded` controls this — if TRUE, `seurat_obj`
  # was already loaded from checkpoint_file_1 above and already contains SCEVAN
  # metadata. `counts_matrix` and `probe_matrix` are left unset intentionally
  # (they are guarded by the same flag in the cleanup block at the bottom).
  if (!scevan_checkpoint_loaded) {

    # ---- 1. Read RNA count matrix from 10x H5 file -------------------------
    rna_h5_path <- file.path(H5_DIR, sample_id, "sample_filtered_feature_bc_matrix.h5")
    if (!file.exists(rna_h5_path)) {
      warning(paste("  [SKIP] RNA H5 file not found for:", sample_id))
      next
    }
    counts_matrix <- Read10X_h5(rna_h5_path)

    # ---- 2. Optional: Read probe data and create custom summed feature ------
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
        common_cells      <- intersect(colnames(counts_matrix), colnames(probe_matrix))
        custom_sum_counts <- numeric(length(common_cells))
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

    # ---- 3. Create Seurat object and attach metadata -----------------------
    seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = sample_id, min.cells = 5)

    # Attach all metadata columns from the Excel file
    for (col_name in colnames(sample_info)) {
      seurat_obj[[col_name]] <- sample_info[[col_name]]
    }
    # Auto-create combined group label (Condition_Sex) if both columns are present
    if ("Condition" %in% colnames(seurat_obj@meta.data) &&
        "Sex" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj$group <- paste(seurat_obj$Condition, seurat_obj$Sex, sep = "_")
    }

    # Add individual probe counts as metadata columns
    if (!is.null(probe_matrix)) {
      for (probe_id in names(PROBE_MAPPING)) {
        if (probe_id %in% rownames(probe_matrix)) {
          col_name_clean <- gsub("[\\s\\(\\)\n]+", "_", PROBE_MAPPING[[probe_id]])
          col_name_clean <- gsub("_$", "", col_name_clean)
          col_name_final <- paste0("Nr4a1_", col_name_clean, "_probe_count")
          seurat_obj <- AddMetaData(seurat_obj,
                                    metadata = probe_matrix[probe_id, ],
                                    col.name = col_name_final)
        }
      }
    }

    # Free raw count matrices — now absorbed into the Seurat object
    rm(counts_matrix)
    if (!is.null(probe_matrix)) {
      rm(probe_matrix)
      rm(list = intersect(c("common_cells", "custom_sum_counts", "probes_found",
                             "summed_vector", "custom_row"), ls()))
    }
    gc()

    # ---- 4. Optional: SCEVAN copy-number variation detection ---------------
    # SCEVAN is run on raw counts before any correction to ensure the most
    # accurate CNA signal. Running per-sample prevents RAM exhaustion in large
    # cohorts.
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
          rm(scevan_res_df)
          message("    -> SCEVAN metadata added to Seurat object.")
        }
      }, error = function(e) {
        message(paste("    -> [WARNING] SCEVAN failed for", sample_id, "| Error:", e$message))
      })
      gc()  # SCEVAN allocates heavily internally; reclaim memory immediately
    }

    # ---- Save Checkpoint 1 (SCEVAN result) ---------------------------------
    # Saved immediately after SCEVAN so that if DecontX or DoubletFinder crash
    # later, the next run can resume from here and skip the expensive SCEVAN step.
    message("    -> Saving SCEVAN checkpoint (checkpoint 1)...")
    saveRDS(seurat_obj, file = checkpoint_file_1)
    gc()
    message("    -> Checkpoint 1 saved.")

  } # end if (!scevan_checkpoint_loaded)

  # ---- 5. Optional: DecontX ambient RNA correction -------------------------
  # DecontX is run per-sample BEFORE doublet detection so that DoubletFinder
  # operates on counts free of ambient RNA contamination.
  if (RUN_DECONTX) {
    tryCatch({
      message(paste("    -> Running DecontX for ambient RNA correction on", sample_id, "..."))
      counts_sparse   <- GetAssayData(object = seurat_obj, layer = "counts")
      decontx_results <- decontX(x = counts_sparse)

      # Replace raw counts with decontaminated counts.
      # Optionally round to integer (ROUND_DECONTX_COUNTS = TRUE) — useful for
      # tools that require integer count matrices. See decontX vignette for details.
      if (ROUND_DECONTX_COUNTS) {
        seurat_obj[["RNA"]]$counts <- round(decontx_results$decontXcounts)
        message("    -> DecontX counts rounded to nearest integer.")
      } else {
        seurat_obj[["RNA"]]$counts <- decontx_results$decontXcounts
      }

      # Recompute QC metrics on the decontaminated matrix
      seurat_obj$nCount_RNA   <- colSums(seurat_obj[["RNA"]]$counts)
      seurat_obj$nFeature_RNA <- colSums(seurat_obj[["RNA"]]$counts > 0)
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")

      rm(decontx_results, counts_sparse)
      gc()
      message(paste("    -> DecontX complete for", sample_id, "."))
    }, error = function(e) {
      message(paste("    -> [WARNING] DecontX failed for", sample_id, "| Error:", e$message,
                    "| Proceeding with uncorrected counts."))
    })
  }

  # ---- 6. DoubletFinder: per-sample doublet detection and removal -----------
  # DoubletFinder requires a temporary normalized + PCA-reduced object.
  # We use a copy (seu_tmp) so the checkpoint contains only clean raw counts
  # without transient normalized/scaled layers or reductions.
  message(paste("    -> Running DoubletFinder for", sample_id, "..."))
  seu_tmp <- seurat_obj

  # Preliminary processing required by DoubletFinder
  seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = N_PCS_TO_USE, verbose = FALSE)

  # pK selection: find the pK that maximizes the BCmetric
  sweep.res.list <- paramSweep(seu_tmp, PCs = 1:N_PCS_TO_USE, sct = FALSE)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)
  bcmvn$pK       <- as.numeric(as.character(bcmvn$pK))

  # Identify the global BCmetric maximum
  initial_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]

  # Constrained pK selection: keep within biologically sensible range
  if (initial_pk < DF_PK_RANGE_MIN || initial_pk > DF_PK_RANGE_MAX) {
    in_range <- bcmvn[bcmvn$pK >= DF_PK_RANGE_MIN & bcmvn$pK <= DF_PK_RANGE_MAX, ]
    if (nrow(in_range) > 0) {
      final_pk <- in_range$pK[which.max(in_range$BCmetric)]
      message(paste("      -> Global peak (", initial_pk,
                    ") out of bounds. Using best in-range peak:", final_pk))
    } else {
      final_pk <- DF_PK_FALLBACK
      message(paste("      -> No peaks in range. Using fallback:", final_pk))
    }
  } else {
    final_pk <- initial_pk
  }

  # Save pK selection diagnostic plot to sample directory
  pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
    ggtitle(paste0("pK Finder: ", sample_id),
            subtitle = paste0("Selected pK = ", final_pk,
                              " | Initial optimal = ", initial_pk)) +
    theme_minimal()
  ggsave(file.path(sample_scevan_dir, paste0(sample_id, "_pk_plot.png")),
         plot = pk_plot, width = 7, height = 5, dpi = DPI_SETTING)
  rm(pk_plot)

  # Run DoubletFinder
  nExp_val <- round(ncol(seu_tmp) * DOUBLET_RATE)
  seu_tmp  <- doubletFinder(seu_tmp, PCs = 1:N_PCS_TO_USE,
                             pK = final_pk, nExp = nExp_val, sct = FALSE)
  res_col  <- tail(grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE), 1)

  # Rollback check: guard against biologically implausible doublet fractions
  n_doublets       <- sum(seu_tmp@meta.data[[res_col]] == "Doublet")
  doublet_fraction <- n_doublets / ncol(seu_tmp)

  if (doublet_fraction > DOUBLET_ROLLBACK_THRESHOLD) {
    warning(paste0(
      "  [ROLLBACK] Sample '", sample_id, "': detected doublet fraction = ",
      round(doublet_fraction * 100, 1), "% exceeds threshold of ",
      round(DOUBLET_ROLLBACK_THRESHOLD * 100, 1), "%. ",
      "Doublet removal SKIPPED for this sample (all cells treated as Singlets)."
    ))
    seurat_obj$Doublet_Status <- "Singlet"
    rollback_log[[sample_id]] <- list(
      fraction   = doublet_fraction,
      n_doublets = n_doublets,
      n_cells    = ncol(seu_tmp),
      final_pk   = final_pk,
      action     = "ROLLBACK_APPLIED"
    )
  } else {
    # Transfer doublet labels from the temp object to the main object
    doublet_labels <- seu_tmp@meta.data[[res_col]]
    names(doublet_labels) <- colnames(seu_tmp)
    seurat_obj$Doublet_Status <- doublet_labels[colnames(seurat_obj)]
    rollback_log[[sample_id]] <- list(
      fraction   = doublet_fraction,
      n_doublets = n_doublets,
      n_cells    = ncol(seu_tmp),
      final_pk   = final_pk,
      action     = "DOUBLETS_REMOVED"
    )
  }

  message(paste0("      -> Doublet fraction: ",
                 round(doublet_fraction * 100, 1), "% | pK used: ", final_pk,
                 " | Action: ", rollback_log[[sample_id]]$action))

  # Diagnostic UMAP: visualize doublets BEFORE removal (reuses PCA from DoubletFinder)
  tryCatch({
    seu_tmp <- RunUMAP(seu_tmp, dims = 1:N_PCS_TO_USE, verbose = FALSE,
                       reduction.name = "umap_doublets")
    p_dblt <- DimPlot(seu_tmp, reduction = "umap_doublets",
                      group.by = res_col) +
      coord_fixed() +
      ggtitle(paste0("Doublets: ", sample_id),
              subtitle = paste0(round(doublet_fraction * 100, 1), "% doublets | Action: ",
                                rollback_log[[sample_id]]$action))
    ggsave(file.path(sample_scevan_dir, paste0(sample_id, "_doublet_visualization.png")),
           plot = p_dblt, width = 7, height = 6, dpi = DPI_SETTING, bg = "white")
    rm(p_dblt)
  }, error = function(e) {
    message(paste("      -> [WARNING] Doublet UMAP failed for", sample_id, ":", e$message))
  })

  # Remove doublets from the main Seurat object
  cells_before_dblt <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj, subset = Doublet_Status == "Singlet")
  message(paste0("      -> Cells: ", cells_before_dblt, " -> ", ncol(seurat_obj),
                 " (removed ", cells_before_dblt - ncol(seurat_obj), " doublets)"))

  # Clean up DoubletFinder temporary objects
  rm(seu_tmp, sweep.res.list, sweep.stats, bcmvn)
  rm(list = intersect(c("doublet_labels", "in_range", "nExp_val",
                         "n_doublets", "doublet_fraction", "initial_pk",
                         "final_pk", "res_col", "cells_before_dblt"), ls()))
  gc()

  # ---- 7. Pre-merge QC: loose per-sample filters ---------------------------
  # Applied after DecontX and doublet removal so QC metrics reflect the
  # corrected, singlet-only count matrix. Removes empty droplets and
  # clearly dying cells before merging.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA >= PRE_MIN_GENES_PER_CELL &
                                percent.mt   <= PRE_MAX_MT_PERCENT)
  message(paste0("      -> Pre-merge QC: ", ncol(seurat_obj), " cells retained."))

  # ---- 8. Save Checkpoint 2 (decontX + doublets + QC) and free memory ------
  # This is the FINAL per-sample checkpoint. It contains clean, singlet-only,
  # ambient-corrected cells and is the file that Step 2.1b will load for merging.
  # Saving here means a crash AFTER this point (e.g., during merge) will not
  # require re-running DecontX or DoubletFinder on the next attempt.
  message("    -> Saving full checkpoint (decontX + doublets + QC) and releasing memory...")
  saveRDS(seurat_obj, file = checkpoint_file_2)
  rm(seurat_obj)
  gc()
  message(paste("  [DONE]", sample_id, "- checkpoint 2 saved, memory freed."))
}

# --- Save DoubletFinder rollback summary (covers all processed samples) ------
if (length(rollback_log) > 0) {
  rollback_df <- do.call(rbind, lapply(names(rollback_log), function(s) {
    data.frame(
      SampleID           = s,
      N_Cells            = rollback_log[[s]]$n_cells,
      N_Doublets_Flagged = rollback_log[[s]]$n_doublets,
      Doublet_Fraction   = round(rollback_log[[s]]$fraction * 100, 2),
      Threshold_Pct      = DOUBLET_ROLLBACK_THRESHOLD * 100,
      pK_Used            = rollback_log[[s]]$final_pk,
      Action             = rollback_log[[s]]$action,
      stringsAsFactors   = FALSE
    )
  }))
  write_xlsx(rollback_df, file.path(OUTPUT_DIR, "doublet_finder_rollback_log.xlsx"))
  message("  Doublet summary saved to: doublet_finder_rollback_log.xlsx")
}

# =============================================================================
# --- STEP 2.1b: Load All Checkpoints into List for Merging ------------------
# =============================================================================
# PURPOSE:
#   Now that all samples have been processed and checkpointed individually
#   (SCEVAN metadata, DecontX correction, and doublet removal already applied),
#   load them sequentially into a named list for the merge step below.
#
#   This is a separate loop so that the heavy per-sample computation (Step 2.1a)
#   never co-exists in memory with a growing list of Seurat objects.
# =============================================================================
message("=== STEP 2.1b: Loading Checkpoints for Merging ===")
seurat_objects_list <- list()

for (i in 1:nrow(metadata)) {
  sample_id         <- as.character(metadata$SampleID[i])
  sample_scevan_dir <- file.path(SCEVAN_DIR, sample_id)
  checkpoint_file_2 <- file.path(sample_scevan_dir, paste0(sample_id, "_decontx_dblt_processed.rds"))
  checkpoint_file_1 <- file.path(sample_scevan_dir, paste0(sample_id, "_scevan_processed.rds"))

  if (file.exists(checkpoint_file_2)) {
    # Normal case: full checkpoint (decontX + doublets + QC) is present.
    message(paste("  [LOAD]", sample_id, "(checkpoint 2 — decontX + doublets + QC)"))
    seurat_objects_list[[sample_id]] <- readRDS(checkpoint_file_2)
  } else if (file.exists(checkpoint_file_1)) {
    # Fallback: only the SCEVAN checkpoint exists (e.g., pipeline from a
    # previous version, or Step 2.1a crashed before saving checkpoint 2).
    # NOTE: this object has NOT had DecontX or DoubletFinder applied — re-run
    # Step 2.1a first to produce the complete checkpoint 2 before merging.
    warning(paste0(
      "  [FALLBACK] Sample '", sample_id, "': only SCEVAN checkpoint found. ",
      "DecontX and DoubletFinder may NOT have been applied. ",
      "Re-run Step 2.1a to generate the full checkpoint before proceeding."
    ))
    seurat_objects_list[[sample_id]] <- readRDS(checkpoint_file_1)
  } else {
    warning(paste("  [MISSING] No checkpoint found for", sample_id, "- skipping from merge."))
  }
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
rm(seurat_objects_list); gc()

output_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_merged_post_metadata_gen.rds"))
saveRDS(data, output_rds, compress = TRUE)

# =============================================================================
# --- STEP 2.3: Post-Merge QC & Diagnostic Plots -----------------------------
# =============================================================================
# Stringent QC filters applied after merging. The BPCells on-disk layer
# approach is used here to prevent out-of-memory errors on large datasets.
# =============================================================================
output_rds_layered <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_merged_post_metadata_gen.rds"))
data <- readRDS(output_rds_layered)

message("\n=== STEP 2.3: Post-Merge QC and Diagnostic Plotting ===")

# --- A) Calculate QC metrics and plot BEFORE filtering ----------------------
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^mt-")
data$all_cells <- PROJECT_NAME

p_before <- tryCatch({
  message("Attempting standard VlnPlot...")
  VlnPlot(
    object   = data,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "all_cells",
    ncol     = 3,
    pt.size  = 0,
    raster   = FALSE
  )
}, error = function(e) {
  message("Standard VlnPlot failed. Falling back to manual ggplot...")
  plot_df <- data@meta.data[, c("all_cells", "nFeature_RNA", "nCount_RNA", "percent.mt")]
  p1 <- ggplot(plot_df, aes(x = all_cells, y = nFeature_RNA, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "nFeature_RNA (Manual)", x = NULL)
  p2 <- ggplot(plot_df, aes(x = all_cells, y = nCount_RNA, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "nCount_RNA (Manual)", x = NULL)
  p3 <- ggplot(plot_df, aes(x = all_cells, y = percent.mt, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "percent.mt (Manual)", x = NULL)
  return(p1 + p2 + p3 + plot_layout(ncol = 3))
})
print(p_before)
ggsave(file.path(OUTPUT_DIR, "01a_qc_violin_before_filtering.png"),
       plot = p_before, width = 10, height = 6, dpi = DPI_SETTING)

# --- B) Apply cell-level QC filters -----------------------------------------
message("  -> Applying cell-level QC filters...")
cells_before <- ncol(data)
data <- subset(data,
               subset = nFeature_RNA >= POST_MIN_GENES  &
                        nFeature_RNA <= POST_MAX_GENES  &
                        nCount_RNA   >= POST_MIN_UMIS   &
                        nCount_RNA   <= POST_MAX_UMIS   &
                        percent.mt   <= POST_MAX_MT)
message(paste("  -> Removed", cells_before - ncol(data), "low-quality cells."))
gc()

# --- C) Layer joining: BPCells on-disk (recommended) or in-memory fallback ---
message("  -> Joining layers into a single matrix...")
if (USE_BPCELLS) {
  # Write each split layer to disk via BPCells before joining.
  # Prevents RAM exhaustion when merging many samples.
  for (i in Layers(data)) {
    write_matrix_dir(
      mat = data[["RNA"]][i],
      dir = file.path(OUTPUT_DIR, paste0("tmpdir/", PROJECT_NAME, "_merged_post_metadata_gen_", i))
    )
    data[["RNA"]][i] <- open_matrix_dir(
      dir = file.path(OUTPUT_DIR, paste0("tmpdir/", PROJECT_NAME, "_merged_post_metadata_gen_", i))
    )
    print(i)
  }
  gc()
  data <- JoinLayers(data)
  data[["RNA"]]$counts <- as(LayerData(data, assay = "RNA", layer = "counts"), "dgCMatrix")
} else {
  # Fallback: standard in-memory JoinLayers (fine for small/medium datasets)
  data <- JoinLayers(data)
  data[["RNA"]]$counts <- as(LayerData(data, assay = "RNA", layer = "counts"), "dgCMatrix")
}
gc()
message("  -> Layers successfully joined.")

# --- D) Apply gene-level filter ---------------------------------------------
message("  -> Applying gene-level filter...")
genes_to_keep <- rownames(data)[
  Matrix::rowSums(GetAssayData(data, layer = "counts") > 0) >= POST_MIN_CELLS_PER_GENE
]
data <- subset(data, features = genes_to_keep)
gc()

# --- E) Plot AFTER filtering ------------------------------------------------
p_after <- tryCatch({
  message("Attempting standard VlnPlot (After)...")
  VlnPlot(
    object   = data,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "all_cells",
    ncol     = 3,
    pt.size  = 0,
    raster   = FALSE
  )
}, error = function(e) {
  message("VlnPlot failed again. Using manual ggplot backup for 'p_after'...")
  plot_df_after <- data@meta.data[, c("all_cells", "nFeature_RNA", "nCount_RNA", "percent.mt")]
  pa1 <- ggplot(plot_df_after, aes(x = all_cells, y = nFeature_RNA, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "nFeature_RNA (Filtered)", x = NULL)
  pa2 <- ggplot(plot_df_after, aes(x = all_cells, y = nCount_RNA, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "nCount_RNA (Filtered)", x = NULL)
  pa3 <- ggplot(plot_df_after, aes(x = all_cells, y = percent.mt, fill = all_cells)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    theme_classic() + NoLegend() + labs(title = "percent.mt (Filtered)", x = NULL)
  return(pa1 + pa2 + pa3 + plot_layout(ncol = 3))
})
print(p_after)
ggsave(file.path(OUTPUT_DIR, "01b_qc_violin_after_filtering.png"),
       plot = p_after, width = 10, height = 6, dpi = DPI_SETTING)
message(paste("  Cells after all QC:", ncol(data), "| Genes retained:", nrow(data)))

# --- F) Save QC-filtered object ---------------------------------------------
output_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_merged_post_qc.rds"))
saveRDS(data, output_rds)

# =============================================================================
# --- STEP 2.4: Normalization, Dimensionality Reduction & Integration --------
# =============================================================================
# This modern workflow uses the IntegrateLayers function with HarmonyIntegration.
# https://satijalab.org/seurat/articles/seurat5_integration
#   1. The RNA assay is split into temporary layers, one per sample.
#   2. Normalization, HVG finding, Scaling, and PCA run on each layer.
#   3. An 'unintegrated' UMAP is created for diagnostics (Track A).
#   4. IntegrateLayers with HarmonyIntegration creates a corrected 'harmony'
#      reduction (Track B).
#   5. Layers are re-joined for downstream analysis.
# =============================================================================
output_rds_layered <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_merged_post_qc.rds"))
data <- readRDS(output_rds_layered)
data[["RNA"]]$data       <- NULL
data[["RNA"]]$scale.data <- NULL
gc()

message("\n=== STEP 2.4: Normalization, Dimensionality Reduction, and Integration (Seurat v5 Method) ===")

# --- Split the object into layers by SampleID for integration ---------------
data[["RNA"]] <- split(data[["RNA"]], f = data$SampleID)
message("  -> RNA assay split into layers by SampleID.")
gc()

# BPCells on-disk storage for split layers (skipped if USE_BPCELLS = FALSE)
if (USE_BPCELLS) {
  for (i in Layers(data)) {
    write_matrix_dir(
      mat = data[["RNA"]][i],
      dir = file.path(OUTPUT_DIR, paste0("tmpdir/", PROJECT_NAME, "_merged_qc_", i))
    )
    data[["RNA"]][i] <- open_matrix_dir(
      dir = file.path(OUTPUT_DIR, paste0("tmpdir/", PROJECT_NAME, "_merged_qc_", i))
    )
    print(i)
  }
  gc()
}

# --- Normalize, find HVGs, scale, and run PCA on each layer -----------------
data <- NormalizeData(data, verbose = TRUE)
data <- FindVariableFeatures(data, nfeatures = N_VARIABLE_FEATURES, verbose = TRUE)
data <- ScaleData(data, verbose = TRUE)
data <- RunPCA(data, npcs = N_PCS_TO_USE, verbose = TRUE)
message("  -> Per-layer normalization, scaling, and PCA complete.")
gc()

# --- Track A: Unintegrated PCA (diagnostic) ---------------------------------
message("  -> Generating unintegrated UMAP (Track A)...")
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                      graph.name = "pca_nn", verbose = FALSE, k.param = UMAP_N_NEIGHBORS)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, algorithm = "leiden",
                     graph.name = "pca_nn", cluster.name = "clusters_none", verbose = FALSE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                reduction.name = "umap_none", verbose = FALSE, n.epochs = 500)
gc()

# --- Track B: Harmony integration -------------------------------------------
message("  -> Running Harmony integration via IntegrateLayers (Track B)...")
data <- IntegrateLayers(
  object         = data,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by       = "SampleID",
  verbose        = TRUE
)
gc()

data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                      graph.name = "harmony_nn", verbose = TRUE, k.param = UMAP_N_NEIGHBORS)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, algorithm = "leiden",
                     graph.name = "harmony_nn", cluster.name = "clusters_harmony", verbose = TRUE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                reduction.name = "umap_harmony", verbose = TRUE, n.epochs = 500)
gc()

# --- Re-join layers for downstream analysis ---------------------------------
message("  -> Re-joining layers into a single matrix for downstream analysis...")
data <- JoinLayers(data)
if (USE_BPCELLS) {
  # Convert BPCells on-disk matrices back to sparse dgCMatrix in memory
  data[["RNA"]]$counts     <- as(LayerData(data, assay = "RNA", layer = "counts"),     "dgCMatrix")
  gc()
  data[["RNA"]]$data       <- as(LayerData(data, assay = "RNA", layer = "data"),       "dgCMatrix")
  gc()
  data[["RNA"]]$scale.data <- as(LayerData(data, assay = "RNA", layer = "scale.data"), "dgCMatrix")
  gc()
}
message("  -> Layers successfully joined.")

# --- Diagnostic UMAP plots --------------------------------------------------
p1 <- DimPlot(data, reduction = "umap_none",    group.by = "SampleID") +
  ggtitle("Unintegrated PCA") + coord_fixed()
p2 <- DimPlot(data, reduction = "umap_harmony", group.by = "SampleID") +
  ggtitle("Harmony Integrated") + coord_fixed()
p1 + p2
ggsave(file.path(OUTPUT_DIR, "02_DIAGNOSTIC_UMAP_Harmony_vs_NoHarmony.png"),
       plot = p1 + p2, width = 16, height = 7, dpi = DPI_SETTING, bg = "white")

p3 <- DimPlot(data, reduction = "umap_none",    group.by = "clusters_none",
              label = TRUE) + ggtitle("Clusters: Unintegrated PCA") + NoLegend() + coord_fixed()
p4 <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony",
              label = TRUE) + ggtitle("Clusters: Harmony") + NoLegend() + coord_fixed()
p3 + p4
ggsave(file.path(OUTPUT_DIR, "03_DIAGNOSTIC_UMAP_Cluster_Comparison.png"),
       plot = p3 + p4, width = 16, height = 7, dpi = DPI_SETTING)
message("  Saved Harmony vs. PCA diagnostic UMAP plots.")

# =============================================================================
# --- STEP 2.5: Save Processed Object ----------------------------------------
# =============================================================================
message("\n=== STEP 2.5: Saving Processed Object for Annotation ===")
output_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
saveRDS(data, output_rds, compress = TRUE)
message(paste0(
  "\n",
  "=== PROCESSING COMPLETE ===\n",
  "  Output saved: ", output_rds, "\n",
  "  Total cells:  ", ncol(data), "\n",
  "  Total genes:  ", nrow(data), "\n",
  "\nNext step: Open and run '02_annotate_data.R' to annotate the data.\n"
))

# Load if needed:
# output_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
# data <- readRDS(output_rds)
