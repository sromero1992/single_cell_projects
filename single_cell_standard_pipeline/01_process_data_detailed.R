# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING (Modular Pipeline)
# Version: 2.1 - EXPLICIT MODULAR STEPS
# =============================================================================
library(TamuScDSC)
library(Seurat)

# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
PROJECT_NAME  <- "Nr4a1_ack17"
ROOT_PATH     <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process/debug_pipeline_pkg"
METADATA_FILE <- file.path(ROOT_PATH, "Nr4a1_s17_metadata.xlsx")
H5_DIR        <- file.path(ROOT_PATH, "h5_files")
OUTPUT_DIR    <- file.path(ROOT_PATH, "seurat_output")
DIAG_DIR      <- file.path(OUTPUT_DIR, "doublet_diagnostics")
CHECKPOINT_DIR<- file.path(OUTPUT_DIR, "per_sample_checkpoints")
SAMPLE_COL    <- "SampleID"

# Thresholds & Parameters
SPECIES <- "mouse"  

mt_pat <- if (tolower(SPECIES) %in% c("human", "hs", "homo_sapiens")) "^MT-" else "^mt-"

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
  method                = "both",
  consensus_rule        = "intersect",
  rate                  = "auto",
  multiplet_rate_per_1k = 0.008,
  n_pcs                 = 50,
  n_variable_features   = 2000,
  df_pn                 = 0.25,
  scdbl_score_threshold = 0.5
)


INTEGRATION_METHOD <- "RunHarmony"
CLUSTER_RESOLUTION <- 1.0
RUN_DECONTX        <- TRUE
DOUBLET_ACTION     <- "remove" # "label" or "remove"
PROBE              <- NULL

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
message("\n[STEP 1/7] Ingesting 10x Raw Matrices...")
sample_list <- read_10x_samples(metadata, H5_DIR, sample_col = SAMPLE_COL, probe = PROBE)

# -----------------------------------------------------------------------------
# STEP 2: AMBIENT RNA DECONTAMINATION (DecontX on Raw Data)
# -----------------------------------------------------------------------------
if (RUN_DECONTX) {
  message("\n[STEP 2/7] Running DecontX Ambient RNA Removal...")
  sample_list <- run_decontx(sample_list, assay = "RNA", round_counts = TRUE, mt_pattern = mt_pat)
}

# -----------------------------------------------------------------------------
# STEP 3: STRINGENT PER-SAMPLE CELL QC (No Gene Filtering Yet)
# -----------------------------------------------------------------------------
message("\n[STEP 3/7] Running Stringent Per-Sample Cell Filtering...")
sample_list <- apply_qc(sample_list, p = qc$pre, apply = TRUE)

# -----------------------------------------------------------------------------
# STEP 4: DOUBLET DETECTION
# -----------------------------------------------------------------------------
message("\n[STEP 4/7] Detecting Doublets (DoubletFinder + scDblFinder)...")
sample_list <- detect_doublets(sample_list, p = dbl, diag_dir = DIAG_DIR)

# -----------------------------------------------------------------------------
# STEP 5: MERGE SAMPLES (Union) & GLOBAL DATASET QC (Gene + Secondary Cell QC)
# -----------------------------------------------------------------------------
message("\n[STEP 5/7] Merging Samples and Applying Dataset-Wide QC...")
merged_obj <- merge_samples(sample_list, project = PROJECT_NAME) # Union merge

# Global gene threshold (e.g., min_cells_per_gene = 15 across the merged dataset)
merged_obj <- apply_qc(merged_obj, p = qc$post, apply = TRUE)


# -----------------------------------------------------------------------------
# STEP 6: INTEGRATION, CLUSTERING & EMBEDDING
# -----------------------------------------------------------------------------
message("\n[STEP 6/7] Performing Batch Integration (Harmony) & UMAP...")
merged_obj <- integrate_data(
  merged_obj,
  method     = INTEGRATION_METHOD,
  group_by   = SAMPLE_COL,
  n_pcs      = dbl$n_pcs,
  n_features = dbl$n_variable_features,
  resolution = CLUSTER_RESOLUTION
)

# -----------------------------------------------------------------------------
# STEP 7: OPTIONAL CLUSTER-LEVEL DOUBLET REVIEW (01b Logic)
# -----------------------------------------------------------------------------
if (dbl$method == "both") {
  message("\n[STEP 7/7] Executing Cluster-Level Doublet Review...")
  rev_params <- cluster_review_params(
    recluster_resolution = 4.0,
    recluster_reduction  = "harmony",
    recluster_dims       = dbl$n_pcs,
    annotation_resolution = CLUSTER_RESOLUTION
  )
  
  merged_obj <- cluster_doublet_review(
    merged_obj,
    p       = rev_params,
    action  = "flag", # Labels clusters without altering graph embeddings
    out_dir = DIAG_DIR
  )
}

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