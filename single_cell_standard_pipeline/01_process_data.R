# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING (Modular Pipeline)
# Version: 2.2 - EXPLICIT MODULAR STEPS
#
# Explicit step-by-step driver over {TamuScDSC}. The order is FIXED and not
# interchangeable: ingest -> (optional SCEVAN) -> DecontX -> per-sample cell QC
# -> doublet detection -> merge + dataset QC -> integrate -> cluster review.
# =============================================================================
library(TamuScDSC)
library(Seurat)

# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
PROJECT_NAME  <- "Nr4a1_ack17"
ROOT_PATH     <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process/debug_pipeline_pkg"
ROOT_PATH     <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process/debug_pipeline_pkg"
METADATA_FILE <- file.path(ROOT_PATH, "Nr4a1_s17_metadata.xlsx")
H5_DIR        <- file.path(ROOT_PATH, "h5_files")
OUTPUT_DIR    <- file.path(ROOT_PATH, "seurat_output")
DIAG_DIR      <- file.path(OUTPUT_DIR, "doublet_diagnostics")
CHECKPOINT_DIR<- file.path(OUTPUT_DIR, "per_sample_checkpoints")
SAMPLE_COL    <- "SampleID"

# Thresholds & Parameters
SPECIES <- "mouse"
mt_pat  <- if (tolower(SPECIES) %in% c("human", "hs", "homo_sapiens")) "^MT-" else "^mt-"

qc <- list(
  # --- Stage 1: Per-Sample Stringent CELL-ONLY QC ---
  pre = list(
    min_features       = 500,     # Stringent cell features cutoff
    max_features       = 16000,
    min_counts         = 1000,    # Stringent cell counts cutoff
    max_counts         = 100000,
    max_mt             = 5.0,     # Stringent mitochondrial cutoff
    mt_pattern         = mt_pat
    # NO min_cells_per_gene here to avoid dropping KO genes per sample!
  ),
  # --- Stage 2: Post-Merge Global DATASET QC ---
  post = list(
    min_cells_per_gene = 15,      # Evaluated across ALL samples combined
    min_features       = 500,     # Secondary sanity check on merged cells
    max_features       = 16000,
    min_counts         = 1000,
    max_counts         = 100000,
    max_mt             = 5.0,
    mt_pattern         = mt_pat
  )
)

dbl <- doublet_params(
  method                = "both",         # "DoubletFinder" | "scDblFinder" | "both"
  consensus_rule        = "intersect",
  rate                  = "auto",         # 10x per-1k-cell auto rate + paramSweep
  multiplet_rate_per_1k = 0.008,
  n_pcs                 = 50,
  n_variable_features   = 2000,
  df_pn                 = 0.25,
  scdbl_score_threshold = 0.5
)

INTEGRATION_METHOD <- "RunHarmony"
CLUSTER_RESOLUTION <- 1.0          # first-pass clustering resolution (1.0 or 1.5)
RUN_DECONTX        <- TRUE
# Use the raw/droplet h5 as the ambient background (DecontX-recommended: the
# ambient profile is estimated from ALL barcodes, not just the filtered cells).
# Expects <H5_DIR>/<SampleID>/<DECONTX_RAW_H5> next to the filtered h5. Set FALSE
# (or if the raw files are absent) to fall back to heuristic within-cell clusters.
DECONTX_USE_BACKGROUND <- TRUE
DECONTX_RAW_H5     <- "sample_raw_feature_bc_matrix.h5"
RUN_SCEVAN         <- FALSE        # optional CNA / malignant-cell inference
SCEVAN_CORES       <- 20          # processors for SCEVAN (par_cores)
PROBE              <- NULL

# Cluster-level doublet review is a TWO-PASS workflow:
#   1st run: REVIEW_ACTION = "flag"  -> re-clusters at 4.0, writes plots, removes
#            nothing. Inspect DIAG_DIR, confirm the flagged clusters are real
#            artefacts.
#   2nd run: REVIEW_ACTION = "remove" -> drops the flagged clusters and re-clusters
#            the clean object at REVIEW_RECLUSTER_RES (1.0 or 1.5).
REVIEW_ACTION        <- "flag"
REVIEW_HIRES         <- 4.0        # high-resolution clustering used for scoring
REVIEW_RECLUSTER_RES <- CLUSTER_RESOLUTION   # resolution after removal (1.0 / 1.5)

# Create working directories
dirs <- c(OUTPUT_DIR, DIAG_DIR, CHECKPOINT_DIR)
invisible(lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# =============================================================================
# --- PART 2: MODULAR STEP-BY-STEP EXECUTION ----------------------------------
# =============================================================================
metadata <- as.data.frame(readxl::read_excel(METADATA_FILE))

# -----------------------------------------------------------------------------
# STEP 1: INGESTION (Raw)
# -----------------------------------------------------------------------------
message("\n[STEP 1] Ingesting 10x Raw Matrices...")
sample_list <- read_10x_samples(metadata, H5_DIR, sample_col = SAMPLE_COL, probe = PROBE)

# -----------------------------------------------------------------------------
# STEP 1b (OPTIONAL): SCEVAN CNA inference on RAW counts, per-sample checkpoints
# -----------------------------------------------------------------------------
# SCEVAN is slow, so each sample's result is cached under CHECKPOINT_DIR and
# reused on re-runs. Runs on raw counts, BEFORE DecontX alters them.
if (RUN_SCEVAN) {
  message("\n[STEP 1b] Running SCEVAN (optional) with per-sample checkpoints...")
  scevan_dir <- file.path(CHECKPOINT_DIR, "scevan")
  if (!dir.exists(scevan_dir)) dir.create(scevan_dir, recursive = TRUE)
  sample_list <- setNames(lapply(names(sample_list), function(sid) {
    ckpt <- file.path(CHECKPOINT_DIR, paste0(sid, "_scevan.rds"))
    if (file.exists(ckpt)) {
      message("   [checkpoint] loading SCEVAN result for ", sid)
      return(readRDS(ckpt))
    }
    obj <- run_scevan(sample_list[[sid]], sample_id = sid,
                      par_cores = SCEVAN_CORES, organism = SPECIES,
                      out_dir = file.path(scevan_dir, sid))
    saveRDS(obj, ckpt)
    obj
  }), names(sample_list))
}

# -----------------------------------------------------------------------------
# STEP 2: AMBIENT RNA DECONTAMINATION (DecontX on Raw Data)
# -----------------------------------------------------------------------------
if (RUN_DECONTX) {
  message("\n[STEP 2] Running DecontX Ambient RNA Removal...")
  sample_list <- run_decontx(
    sample_list, assay = "RNA", round_counts = TRUE, mt_pattern = mt_pat,
    background_dir = if (DECONTX_USE_BACKGROUND) H5_DIR else NULL,
    background_h5  = DECONTX_RAW_H5, sample_col = SAMPLE_COL
  )
}

# -----------------------------------------------------------------------------
# STEP 3: STRINGENT PER-SAMPLE CELL QC (No Gene Filtering Yet)
# -----------------------------------------------------------------------------
message("\n[STEP 3] Running Stringent Per-Sample Cell Filtering...")
sample_list <- apply_qc(sample_list, p = qc$pre, apply = TRUE)

# -----------------------------------------------------------------------------
# STEP 4: DOUBLET DETECTION (paramSweep + auto per-1k rate, get scores)
# -----------------------------------------------------------------------------
message("\n[STEP 4] Detecting Doublets (DoubletFinder + scDblFinder)...")
sample_list <- detect_doublets(sample_list, p = dbl, diag_dir = DIAG_DIR)

# -----------------------------------------------------------------------------
# STEP 5: MERGE SAMPLES (Union) & GLOBAL DATASET QC (Gene + Secondary Cell QC)
# -----------------------------------------------------------------------------
message("\n[STEP 5] Merging Samples and Applying Dataset-Wide QC...")
merged_obj <- merge_samples(sample_list, project = PROJECT_NAME)
merged_obj <- apply_qc(merged_obj, p = qc$post, apply = TRUE)   # global gene threshold

# -----------------------------------------------------------------------------
# STEP 6: INTEGRATION, CLUSTERING & EMBEDDING (first-pass resolution)
# -----------------------------------------------------------------------------
message("\n[STEP 6] Performing Batch Integration (Harmony) & UMAP...")
merged_obj <- integrate_data(
  merged_obj,
  method     = INTEGRATION_METHOD,
  group_by   = SAMPLE_COL,
  n_pcs      = dbl$n_pcs,
  n_features = dbl$n_variable_features,
  resolution = CLUSTER_RESOLUTION
)

# -----------------------------------------------------------------------------
# STEP 7: CLUSTER-LEVEL DOUBLET REVIEW (high-res score -> flag -> remove)
# -----------------------------------------------------------------------------
# Scores the two doublet metrics over a high-resolution (REVIEW_HIRES) clustering,
# flags artefact clusters, and (when REVIEW_ACTION = "remove") drops them and
# re-clusters the clean object at REVIEW_RECLUSTER_RES. Needs method = "both".
if (dbl$method == "both") {
  message("\n[STEP 7] Cluster-Level Doublet Review (action = ", REVIEW_ACTION, ")...")
  rev_params <- cluster_review_params(
    recluster_resolution  = REVIEW_HIRES,
    recluster_reduction   = "harmony",
    recluster_dims        = dbl$n_pcs,
    annotation_resolution = REVIEW_RECLUSTER_RES
  )
  merged_obj <- cluster_doublet_review(
    merged_obj, p = rev_params, action = REVIEW_ACTION, out_dir = DIAG_DIR
  )
}


feats <- c("Epcam", "Cd3e", "Pecam1", "Cd19", "Cd68", "Sox10", "Tubb3", "Col1a1", "Myh11", "Sox4", "Sox9")
FeaturePlot(merged_obj, reduction = "umap_harmony", features = feats)
DimPlot(merged_obj, reduction = "umap_harmony", group.by = "clusters_harmony")

# =============================================================================
# --- PART 3: SAVE & AUDIT -----------------------------------------------------
# =============================================================================
out_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
saveRDS(merged_obj, out_rds)

message("\n=================================================================")
message(" Execution Summary & Provenance Log:")
message("=================================================================")
provenance(merged_obj)
message("\nCompleted successfully -> Output saved to: ", basename(out_rds))
message("  Cluster review action was '", REVIEW_ACTION,
        "'. Inspect ", DIAG_DIR, " then set REVIEW_ACTION <- \"remove\" and re-run STEP 7.")
