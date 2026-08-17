# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING  (thin driver over {scprep})
# Version: 2.0 - PACKAGED
#
# This is the small replacement for the old 2,000-line 01_process_data.R. All
# of the heavy logic (QC, DecontX, the DoubletFinder + scDblFinder engine,
# merging, Harmony integration) now lives in the installable {scprep} package,
# so this file is just configuration plus a single recipe call.
#
#   - The full original script is preserved as 01_process_data_legacy.R.
#   - Install the package first:  devtools::install("scprep")
#   - Cluster-level doublet review is the separate 01b_cluster_doublet_review.R.
#
# THREE WAYS TO SUPPLY INPUT (pick one for `input` below):
#   1. Raw Cell Ranger H5s   -> read_10x_samples(metadata, H5_DIR)
#   2. Per-sample .rds files  -> read_checkpoints("<dir of per-sample rds>")
#   3. A big merged .rds       -> readRDS("<path>")   (split internally by SampleID)
# =============================================================================

library(scprep)
library(Seurat)

# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
PROJECT_NAME  <- "Wu_Diet_project2"
ROOT_PATH     <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
METADATA_FILE <- file.path(ROOT_PATH, "Wu_p2_metadata.xlsx")
H5_DIR        <- file.path(ROOT_PATH, "h5_files")
OUTPUT_DIR    <- file.path(ROOT_PATH, "seurat_output")
DIAG_DIR      <- file.path(OUTPUT_DIR, "doublet_diagnostics")
SAMPLE_COL    <- "SampleID"

# --- QC thresholds (loose per-sample -> stringent merged) --------------------
qc <- qc_params(
  min_features           = 500,    # pre-merge (loose)
  max_mt                 = 20,     # pre-merge (loose)
  mt_pattern             = "^MT-|^mt-",
  min_features_stringent = 500,    # post-merge (stringent)
  max_mt_stringent       = 5
)

# --- Doublet detection (your two-score engine) -------------------------------
# method: "DoubletFinder" | "scDblFinder" | "both"  (use "both" to enable 01b).
dbl <- doublet_params(
  method                = "DoubletFinder",
  consensus_rule        = "intersect",
  rate                  = "auto",
  expected_cells_per_sample = NULL,
  multiplet_rate_per_1k = 0.008,
  rate_floor            = 0.005,
  rate_ceiling          = 0.20,
  n_pcs                 = 50,
  n_variable_features   = 2000,
  df_pn                 = 0.25,
  scdbl_score_threshold = 0.5
)

INTEGRATION_METHOD <- "RunHarmony"   # "RunHarmony" | "IntegrateLayers"
CLUSTER_RESOLUTION <- 1.0
RUN_DECONTX        <- TRUE
# Doublets are LABELLED, not removed, so 01b can review them first. To drop them
# per-sample instead, set this to "remove" (acts on the consensus call).
DOUBLET_ACTION     <- "label"

# --- Optional probe summing (Flex / KO designs). NULL = off. -----------------
# PROBE <- list(h5_name = "sample_raw_probe_bc_matrix.h5", feature = "Nr4a1_cust",
#               probes_for_sum = c("<probe_id_1>", "<probe_id_2>"))
PROBE <- NULL

# =============================================================================
# --- PART 2: RUN -------------------------------------------------------------
# =============================================================================
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(DIAG_DIR))   dir.create(DIAG_DIR,   recursive = TRUE)

metadata <- as.data.frame(readxl::read_excel(METADATA_FILE))

# --- Input: choose ONE of the three scenarios ------------------------------
input <- read_10x_samples(metadata, H5_DIR, sample_col = SAMPLE_COL, probe = PROBE)
# input <- read_checkpoints(file.path(OUTPUT_DIR, "processed_samples"))   # scenario 2
# input <- readRDS(file.path(OUTPUT_DIR, "big_object.rds"))               # scenario 3

# --- One recipe call does QC -> DecontX -> doublets -> merge -> integrate ---
data <- recipe_per_sample(
  input, metadata = metadata, sample_col = SAMPLE_COL,
  qc = qc, dbl = dbl, decontx = RUN_DECONTX,
  integration = INTEGRATION_METHOD, n_pcs = dbl$n_pcs,
  n_features = dbl$n_variable_features, resolution = CLUSTER_RESOLUTION,
  diag_dir = DIAG_DIR
)

if (DOUBLET_ACTION == "remove")
  data <- apply_doublet_action(data, action = "remove", target = "consensus")

# --- Save for annotation (Script 02) ----------------------------------------
out_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
saveRDS(data, out_rds)

provenance(data)
message("\n=== Script 01 complete -> ", basename(out_rds), " ===")
message("  Next: 01b_cluster_doublet_review.R (if DOUBLET_METHOD = 'both'), then 02.")
