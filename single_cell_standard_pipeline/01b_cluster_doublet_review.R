# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1b: CLUSTER-LEVEL DOUBLET REVIEW
#   thin driver over {scprep}
# Version: 2.0 - PACKAGED
#
# Runs on the processed object BEFORE annotation. Re-clusters at high resolution,
# averages the two doublet scores per cluster, flags clusters with MAD-adaptive
# cutoffs, and (optionally) removes them and re-clusters the clean object.
# All logic lives in scprep::cluster_doublet_review(); the full original script
# is preserved as 01b_cluster_doublet_review_legacy.R.
#
# Requires Script 01 to have been run with DOUBLET_METHOD = "both" so that both
# DF_score and scDblFinder_score are present.
#
# REVIEW, DON'T AUTOPILOT: default action = "flag". Open the plots in REVIEW_DIR,
# confirm the flagged clusters are real artefacts, THEN set ACTION <- "remove".
# =============================================================================

library(scprep)
library(Seurat)

# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
PROJECT_NAME <- "Wu_Diet_project2"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")
RDS_PATH     <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_processed_for_annotation.rds"))
REVIEW_DIR   <- file.path(OUTPUT_DIR, "cluster_doublet_review")

ACTION <- "flag"   # "flag" (label + plots only) or "remove" (drop + re-cluster)

rev <- cluster_review_params(
  recluster            = TRUE,
  recluster_resolution = 4.0,     # 2.0-3.0 also useful; try a few
  recluster_reduction  = "harmony",
  recluster_dims       = 50,      # match n_pcs from Script 01
  flag_rule            = "hybrid",# "hybrid" | "either" | "both"
  weak_k               = 1.5,     # moderate-outlier line (median + k*MAD)
  strong_k             = 3.5,     # extreme-outlier line
  df_mean_floor        = 0.30,
  scdbl_mean_floor     = 0.30,
  fraction_threshold   = 0.60,
  annotation_resolution = 1.0     # re-cluster resolution after removal
)

# =============================================================================
# --- PART 2: RUN -------------------------------------------------------------
# =============================================================================
data <- readRDS(RDS_PATH)

data <- cluster_doublet_review(data, p = rev, action = ACTION, out_dir = REVIEW_DIR)

if (ACTION == "remove") {
  clean_path <- file.path(OUTPUT_DIR,
                          paste0(PROJECT_NAME, "_processed_for_annotation_dbl_removed.rds"))
  saveRDS(data, clean_path)
  message("  Cleaned object saved: ", basename(clean_path), "  -> point Script 02 here.")
}

message("\n=== Script 1b complete. Review plots in: ", REVIEW_DIR, " ===")
