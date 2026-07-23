# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 1: DATA PROCESSING ENGINE
# Version: 11.0 - UNIFIED (dual-method doublet detection)
#
# UNIFIED BUILD PROVENANCE:
#   This script merges the best components of three prior project pipelines:
#     - 2026_nr4a1 / github_pipeline (v10.0) : overall architecture, SCEVAN +
#       DecontX dual-checkpoint system, probe integration, Harmony integration,
#       BPCells memory handling, rollback safeguard, scDblFinder path.
#     - eithan_coffee (01_process_data_SG.R) : DoubletFinder paramSweep loop
#       with constrained pK selection ("safe zone" + fallback) and per-sample
#       rollback logging.
#     - wu_project1 (01_processing_samples.R) : enhanced pK diagnostic plot
#       (grey trace, red dashed line at selection, red diamond marking the
#       selected peak) and the explicit reasonable-range/fallback messaging.
#   See CHANGELOG.md in this folder for a full file-by-file provenance table.
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
#     6. Doublet detection + removal (per-sample) via a SELECTABLE METHOD
#        (DoubletFinder | scDblFinder | both), with a ROLLBACK SAFEGUARD: if
#        the detected doublet fraction exceeds DOUBLET_ROLLBACK_THRESHOLD,
#        removal is skipped for that sample.
#     7. Post-merge QC (stringent filters applied after merging all samples)
#     8. Data normalization, dimensionality reduction, and batch-correction
#        via Harmony, producing parallel (PCA vs. Harmony) clusterings.
#     9. Saving a fully processed, un-annotated Seurat object for Script 02.
#
# PER-SAMPLE PROCESSING ORDER (STEP 2.1a):
#   Load → [SCEVAN] → Pre-QC → [DecontX] → Doublet Detection →
#   Remove Doublets → Checkpoint
#   This order ensures that:
#     - SCEVAN sees unmodified raw counts (best for CNA detection)
#     - Doublet detection operates on ambient-corrected counts
#     - Empty drops are gone before doublet calling, but filtering stays light
#       (stringent filtering biases the doublet-rate estimate)
#     - Each checkpoint contains clean, singlet-only cells
#
# ============================================================================
# DOUBLET DETECTION: TWO METHODS, ONE SWITCH  (DOUBLET_METHOD, Section 1.4)
# ============================================================================
#   "DoubletFinder"  (DEFAULT) — McGinnis et al. 2019, Cell Systems.
#       Simulates artificial doublets by averaging random cell pairs, embeds
#       them with the real cells, and scores each real cell by the proportion
#       of artificial nearest neighbours (pANN) in a neighbourhood of size pK.
#       Requires a paramSweep to choose pK; this pipeline runs the sweep and
#       then CONSTRAINS the choice to a biologically plausible pK window
#       (DF_PK_RANGE_MIN..MAX) with a hardcoded fallback, because the raw
#       BCmetric global maximum frequently lands at implausible extremes.
#       Chosen as the default here on the basis of observed real-world
#       behaviour on this lab's colon/immune datasets, where scDblFinder was
#       over-aggressive and stripped biologically real transitional and
#       proliferating populations.
#
#   "scDblFinder"    — Germain et al. 2021, F1000Research.
#       Gradient-boosted classifier on kNN features from artificial doublets.
#       Fully automatic (no pK), faster, and better in published benchmarks.
#       Retained as a switchable alternative and as a cross-check.
#
#   "both"           — Runs BOTH methods, writes a CONCORDANCE table and a
#       confusion matrix per sample, and removes cells according to
#       DOUBLET_CONSENSUS_RULE ("union", "intersect", or "DoubletFinder" /
#       "scDblFinder" to let one method decide while still logging the other).
#       Use this once on a new dataset to decide which method to trust, then
#       switch to the single method for production runs.
#
#   Regardless of method, the resulting per-cell call is standardised into a
#   single metadata column, Doublet_Status ("Singlet" / "Doublet"), so every
#   downstream script (02-08) is method-agnostic.
# ============================================================================
#
# CHECKPOINT SYSTEM (dual-checkpoint per sample):
#   Each sample produces two checkpoint .rds files inside PROCESSED_DIR/<SampleID>/:
#     _scevan_processed.rds       — saved right after SCEVAN (Checkpoint 1)
#     _decontx_dblt_processed.rds — saved after DecontX + doublet removal + Pre-QC
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
#   - pK sweep plots, doublet score histograms and doublet UMAPs per sample
#     (all in OUTPUT_DIR/doublet_diagnostics/)
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
library(DoubletFinder)  # Doublet detection - paramSweep/pK method (DEFAULT)
library(scDblFinder)    # Doublet detection - ML method (alternative / cross-check)
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
PROJECT_NAME <- "Wu_Diet_project2"

# ROOT_PATH: The top-level working directory for this project.
# All other paths below are derived from this one.
#ROOT_PATH <- "Z:/selim_working_dir/2026_nr4a1_ack/r_process"   # Windows (RStudio local)
#ROOT_PATH <- "/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process"  # Linux/HPC
ROOT_PATH <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"

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
METADATA_FILE <- file.path(ROOT_PATH, "Wu_p2_metadata.xlsx")

# H5_DIR: Directory containing one subfolder per sample, each holding
#   the 10x H5 output files. Expected structure:
#   H5_DIR/<SampleID>/sample_filtered_feature_bc_matrix.h5
#   H5_DIR/<SampleID>/sample_raw_probe_bc_matrix.h5  (optional, for probe data)
H5_DIR        <- file.path(ROOT_PATH, "h5_files")

# OUTPUT_DIR: All processed objects and plots are saved here.
OUTPUT_DIR    <- file.path(ROOT_PATH, "seurat_output")

# PROCESSED_DIR: Per-sample working directory. Holds the resume checkpoints
#   for every sample, and (when RUN_SCEVAN = TRUE) that sample's SCEVAN output
#   in its own scevan/ subfolder. Named generically because it is used whether
#   or not SCEVAN runs - the checkpoints are always written here.
PROCESSED_DIR <- file.path(OUTPUT_DIR, "processed_samples")

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

# =============================================================================
# --- 1.4: DOUBLET DETECTION PARAMETERS ---------------------------------------
# =============================================================================

# --- 1.4.0: METHOD SELECTION -------------------------------------------------
# DOUBLET_METHOD: Which doublet caller to use. One of:
#   "DoubletFinder" — (DEFAULT) paramSweep + constrained pK. Slower but, on
#                     this lab's datasets, more conservative and less prone to
#                     deleting real transitional/proliferating populations.
#   "scDblFinder"   — automatic ML classifier. Faster, no pK tuning.
#   "both"          — run both, log concordance, and resolve the final call
#                     using DOUBLET_CONSENSUS_RULE below.
DOUBLET_METHOD <- "DoubletFinder"

# DOUBLET_CONSENSUS_RULE: Only used when DOUBLET_METHOD == "both".
#   "union"         — a cell is a doublet if EITHER method calls it (aggressive)
#   "intersect"     — a cell is a doublet only if BOTH agree (conservative;
#                     recommended, highest precision)
#   "DoubletFinder" — DoubletFinder decides; scDblFinder is logged only
#   "scDblFinder"   — scDblFinder decides; DoubletFinder is logged only
DOUBLET_CONSENSUS_RULE <- "intersect"

# --- 1.4.1: Expected doublet rate (applies to BOTH methods) ------------------
#
# *** THIS IS THE MOST CONSEQUENTIAL SETTING IN THE DOUBLET STEP. ***
# Under DoubletFinder the rate sets nExp directly (nExp = n_cells * rate),
# which is a HARD QUOTA, not a ceiling: it will call approximately that many
# doublets whether or not they exist. Too high silently deletes real cells.
#
# You only ever need to answer ONE question here:
#   "How many cells did each lane recover?"
# The 10x multiplet rate follows from that (~0.8% per 1,000 cells recovered):
#
#        5,000 cells ->  4%        16,000 cells -> 12.8%
#       10,000 cells ->  8%        20,000 cells -> 16.0%
#
# EXPECTED_CELLS_PER_SAMPLE — how to answer that question:
#   NULL      (RECOMMENDED) use each sample's OWN cell count, measured after
#             the loose pre-QC. Adapts automatically when lanes differ, which
#             they almost always do.
#   <number>  assume every lane targeted this many cells, e.g. 20000.
#             Use when you know the loading target and would rather state it
#             than have it inferred. NOTE this gives every sample the SAME
#             rate, so a lane that under-recovered will be over-corrected.
EXPECTED_CELLS_PER_SAMPLE <- NULL

# -----------------------------------------------------------------------------
# MULTIPLEXED RUNS (10x Flex / CellPlex / hashing) - IMPORTANT
# -----------------------------------------------------------------------------
# If several samples shared one GEM well, do NOT feed the GEM-well total here.
# Leave EXPECTED_CELLS_PER_SAMPLE = NULL and let each sample use its own count.
#
# Why that is the correct answer rather than an approximation:
#   Cell Ranger demultiplexes by probe barcode and ALREADY removes multiplets
#   whose cells came from different samples. What survives into each per-sample
#   H5 is only the UNDETECTED multiplets - those where both cells happened to
#   carry the same barcode. For k equally sized barcodes that is a fraction
#   1/k of all multiplets. So:
#
#       undetected rate  =  (rate on the whole GEM well) / k
#                        =  (c * n_total) / k          [c = slope per cell]
#                        =  c * (n_total / k)
#                        =  the plain rule applied to ONE sample's own count
#
#   The model is linear, so dividing the cells by k and dividing the rate by k
#   land on the same number. Using each sample's measured count is exact here,
#   not a shortcut.
#
#   Worked example (this project): a 4-plex Flex GEM well holding ~80,000 cells
#   yields ~20,000 cells per probe barcode.
#       whole well : 0.008 * 80  = 64%,  /4 barcodes = 16%
#       per sample : 0.008 * 20  = 16%          <- identical
#
# CAVEAT: 10x caps Flex at 20,000 cells per probe barcode and recommends 4,000
#   as a starting point, noting the undetected multiplet rate rises as you
#   approach the cap. At ~20k per barcode a ~16% rate is expected and real, not
#   an artifact - but the doublet step is then removing a large number of cells,
#   so read the Action and Target_Rate_Pct columns of the doublet log carefully.
# -----------------------------------------------------------------------------

# DOUBLET_RATE — normally leave as "auto".
#   "auto"    derive the rate from EXPECTED_CELLS_PER_SAMPLE above.
#   <number>  bypass the model entirely and force a fixed fraction (e.g. 0.16)
#             for every sample. Only if you have a specific reason.
DOUBLET_RATE <- "auto"

# --- Internals of the 10x model (rarely touched) -----------------------------
# Slope of the published linear multiplet curve, and sanity clamps on the
# result. The clamps only stop absurd values; they are not tuning knobs.
MULTIPLET_RATE_PER_1K <- 0.008   # 0.8% per 1,000 cells recovered
DOUBLET_RATE_FLOOR    <- 0.005   # never below 0.5%
DOUBLET_RATE_CEILING  <- 0.20    # never above 20% (must stay < rollback below)

# DOUBLET_ROLLBACK_THRESHOLD: Safety mechanism (applies to BOTH methods).
#   If the chosen method classifies MORE than this fraction of cells as
#   doublets for a given sample, the result is suspicious. In that case,
#   doublet removal is ROLLED BACK for that sample (all cells kept) and a
#   warning is issued. A >25% doublet rate is biologically implausible for
#   standard 10x runs and most likely reflects a classification artifact.
#   NOTE: if DOUBLET_RATE is set at/above this value, rollback will trigger on
#   every sample under DoubletFinder (since nExp is a hard target). Keep
#   DOUBLET_ROLLBACK_THRESHOLD comfortably above DOUBLET_RATE. With
#   DOUBLET_RATE = "auto" this is enforced automatically, since the automatic
#   estimate is capped at DOUBLET_RATE_CEILING.
DOUBLET_ROLLBACK_THRESHOLD <- 0.25

# --- 1.4.2: DoubletFinder-Specific Parameters --------------------------------
# DF_PN: Proportion of artificial doublets to generate, as a fraction of the
#   merged real+artificial cell pool. McGinnis et al. show results are largely
#   invariant to pN; 0.25 is the published default and should rarely change.
DF_PN <- 0.25

# DF_PK_RANGE_MIN / DF_PK_RANGE_MAX: The "safe zone" for pK selection.
#   The paramSweep BCmetric curve often peaks at implausible extremes (e.g.
#   pK = 0.005 or 0.30), which produces nonsense neighbourhoods. If the global
#   BCmetric maximum falls outside this window, the pipeline instead selects
#   the highest peak WITHIN the window. Range [0.01, 0.15] covers the values
#   reported as optimal across the published DoubletFinder validation datasets.
DF_PK_RANGE_MIN <- 0.01
DF_PK_RANGE_MAX <- 0.25

# DF_PK_FALLBACK: Used only if NO peak exists inside the safe zone at all.
#   0.09 is the median optimal pK across the DoubletFinder validation sets and
#   is a safe, well-behaved default.
DF_PK_FALLBACK <- 0.09

# DF_HOMOTYPIC_ADJUST: Adjust nExp downward for HOMOTYPIC doublets.
#   DoubletFinder can only detect HETEROtypic doublets (two different cell
#   types). Homotypic doublets (two cells of the SAME type) are statistically
#   invisible to it. modelHomotypic() estimates their proportion from the
#   cluster-annotation frequencies and shrinks nExp accordingly, so the
#   algorithm is not forced to over-call heterotypic doublets to hit its quota.
#   This is the McGinnis-recommended behaviour. Requires a clustering to exist
#   on the per-sample object, which this pipeline computes automatically.
#   Set FALSE to use the raw, unadjusted nExp.
DF_HOMOTYPIC_ADJUST <- TRUE

# DF_CLUSTER_RES: Clustering resolution used ONLY for the homotypic adjustment
#   above (not used for any downstream biology). Coarse is fine.
DF_CLUSTER_RES <- 0.5

# DF_SCT: Set TRUE only if the per-sample object was normalised with SCTransform.
#   This pipeline uses standard LogNormalize per sample, so keep FALSE.
DF_SCT <- FALSE

# --- 1.4.4: What to DO with the doublets once they are called ----------------
# REMOVE_DOUBLETS: whether Script 01 deletes them, or only labels them.
#
#   FALSE  (DEFAULT) LABEL ONLY. Every cell is kept and carries a
#          Doublet_Status of "Singlet" or "Doublet". Nothing is deleted at any
#          point in Script 01. You filter when YOU choose - typically after
#          annotation in Script 02, once you can see whether the flagged cells
#          form a sensible scatter or land on a real biological cluster.
#          Filter with:
#              data <- subset(data, subset = Doublet_Status == "Singlet")
#
#   TRUE   Delete them immediately, per sample, before the merge. Every
#          downstream object is singlet-only and the decision is irreversible
#          without re-running from the checkpoints.
#
# Keeping them (FALSE) is the safer default: a doublet call is a prediction,
# not a measurement, and this way a bad pK or an over-aggressive rate costs you
# nothing permanent. The trade is that doublets are present during clustering
# and annotation, where they can form their own small hybrid clusters - which
# is often informative in itself.
REMOVE_DOUBLETS <- FALSE

# --- 1.4.3: scDblFinder-Specific Parameters ----------------------------------
# SCDBLFINDER_SCORE_THRESHOLD: Minimum doublet score to call a cell a doublet.
#   scDblFinder outputs scores in [0,1] where 1 = most likely doublet.
#   The default threshold is 0.5 (cells with score >= 0.5 are called doublets).
#   Lower to be more aggressive (catch more doublets, more false positives).
#   Raise to be more conservative (fewer false positives, may miss real doublets).
SCDBLFINDER_SCORE_THRESHOLD <- 0.5

# SCDBLFINDER_DBR_NULL: If TRUE, pass dbr = NULL so scDblFinder estimates the
#   doublet rate from cell number itself instead of using DOUBLET_RATE.
#   Recommended when the cell loading density is unknown.
SCDBLFINDER_DBR_NULL <- FALSE

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

ROUND_DECONTX_COUNTS   <- TRUE # Round decontX corrected counts to nearest integer.
                                 # DecontX produces non-integer (fractional) values by default.
                                 # Set TRUE if downstream tools require integer count matrices.
                                 # See decontX vignette for details.

# INTEGRATION_METHOD: how Harmony batch correction is run. Same maths either
# way (Harmony grouping on SampleID, producing a 'harmony' reduction), but the
# object handling differs sharply in cost:
#
#   "RunHarmony"      (DEFAULT) The classic harmony::RunHarmony() call on the
#                     JOINED object. Normalization, HVGs, scaling and PCA run
#                     once on the whole matrix, then RunHarmony() corrects the
#                     PCA. No layer splitting at any point. Cheaper in both time
#                     and memory, and makes USE_BPCELLS largely unnecessary.
#                     This is the original wu_project1 behaviour.
#
#   "IntegrateLayers" The Seurat v5 layered workflow: split the RNA assay into
#                     one layer per sample, process per layer, then
#                     IntegrateLayers(HarmonyIntegration), then re-join. This is
#                     the split/join cycle that is expensive on large cohorts
#                     and the reason BPCells was added. Kept for reproducibility
#                     with runs made under that workflow.
#
# NOTE on the harmony package: both the CRAN release ("harmony") and the GitHub
# dev version expose RunHarmony() with args reduction.use / reduction.save /
# dims.use. The very old argument name 'reduction=' no longer works; this
# pipeline uses the current names, so either install is fine.
INTEGRATION_METHOD     <- "RunHarmony"  # "RunHarmony" | "IntegrateLayers"

USE_BPCELLS            <- FALSE  # Use BPCells on-disk matrix handling (write_matrix_dir) during
                                 # the post-merge QC (Step 2.3) and integration (Step 2.4) steps.
                                 # Recommended for large cohorts (>8 samples) to prevent RAM
                                 # exhaustion. Set FALSE for small datasets or if the BPCells
                                 # package is not available; standard in-memory JoinLayers() is
                                 # used as fallback. Requires hdf5r and BPCells to be installed.

RUN_SCEVAN             <- FALSE # Run SCEVAN copy-number variation analysis.
                                 # Enable if ANY sample in the cohort is cancer/tumor — run on ALL
                                 # samples including normals, as normals serve as the CNA reference
                                 # baseline. Set FALSE only if the entire study uses normal tissue.
DPI_SETTING            <- 300   # DPI for all saved diagnostic plots.

# --- 1.7: Probe / KO Gene Integration Parameters ---
# Enable this if your 10x run included a probe-based capture assay (e.g., for
# tracking KO efficiency at the probe level). Set to FALSE to skip entirely.
ADD_PROBE_DATA <- FALSE

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
if (!dir.exists(PROCESSED_DIR)) { dir.create(PROCESSED_DIR, recursive = TRUE) }

# DOUBLET_DIAG_DIR: all doublet diagnostics (pK plots, score histograms,
# concordance matrices) are collected here in one place for easy review.
DOUBLET_DIAG_DIR <- file.path(OUTPUT_DIR, "doublet_diagnostics")
if (!dir.exists(DOUBLET_DIAG_DIR)) { dir.create(DOUBLET_DIAG_DIR, recursive = TRUE) }

# =============================================================================
# --- CONFIGURATION VALIDATION ------------------------------------------------
# =============================================================================
# Fail fast and loudly on impossible parameter combinations, rather than
# discovering them 6 hours into a run.
validate_config <- function() {
  errs <- character(0)

  if (!is.logical(REMOVE_DOUBLETS) || length(REMOVE_DOUBLETS) != 1 ||
      is.na(REMOVE_DOUBLETS)) {
    errs <- c(errs, "REMOVE_DOUBLETS must be TRUE or FALSE.")
  }
  if (!INTEGRATION_METHOD %in% c("RunHarmony", "IntegrateLayers")) {
    errs <- c(errs, paste0("INTEGRATION_METHOD must be 'RunHarmony' or ",
                           "'IntegrateLayers'. Got: '", INTEGRATION_METHOD, "'"))
  }
  if (!DOUBLET_METHOD %in% c("DoubletFinder", "scDblFinder", "both")) {
    errs <- c(errs, paste0("DOUBLET_METHOD must be one of 'DoubletFinder', ",
                           "'scDblFinder', or 'both'. Got: '", DOUBLET_METHOD, "'"))
  }
  if (DOUBLET_METHOD == "both" &&
      !DOUBLET_CONSENSUS_RULE %in% c("union", "intersect", "DoubletFinder", "scDblFinder")) {
    errs <- c(errs, paste0("DOUBLET_CONSENSUS_RULE must be one of 'union', ",
                           "'intersect', 'DoubletFinder', 'scDblFinder'. Got: '",
                           DOUBLET_CONSENSUS_RULE, "'"))
  }
  # DOUBLET_RATE must be either the string "auto" or a proper fraction.
  rate_is_auto <- is.character(DOUBLET_RATE) && identical(DOUBLET_RATE, "auto")
  if (!rate_is_auto) {
    if (is.null(DOUBLET_RATE) || !is.numeric(DOUBLET_RATE) ||
        DOUBLET_RATE <= 0 || DOUBLET_RATE >= 1) {
      errs <- c(errs, paste0(
        "DOUBLET_RATE must be either \"auto\" or a number strictly between ",
        "0 and 1. Got: ", paste(deparse(DOUBLET_RATE), collapse = "")))
    }
  }
  if (DOUBLET_METHOD %in% c("DoubletFinder", "both")) {
    if (DF_PK_RANGE_MIN >= DF_PK_RANGE_MAX) {
      errs <- c(errs, "DF_PK_RANGE_MIN must be less than DF_PK_RANGE_MAX.")
    }
    if (DF_PN <= 0 || DF_PN >= 1) {
      errs <- c(errs, "DF_PN must be strictly between 0 and 1 (published default: 0.25).")
    }
  }
  # Auto-mode sanity.
  if (rate_is_auto) {
    if (!is.null(EXPECTED_CELLS_PER_SAMPLE)) {
      if (!is.numeric(EXPECTED_CELLS_PER_SAMPLE) || EXPECTED_CELLS_PER_SAMPLE <= 0) {
        errs <- c(errs, paste0(
          "EXPECTED_CELLS_PER_SAMPLE must be NULL or a positive number of ",
          "cells (e.g. 20000). Got: ",
          paste(deparse(EXPECTED_CELLS_PER_SAMPLE), collapse = "")))
      } else if (EXPECTED_CELLS_PER_SAMPLE < 100) {
        # Catches the classic mistake of entering a RATE (0.16) where a CELL
        # COUNT (20000) is expected.
        errs <- c(errs, paste0(
          "EXPECTED_CELLS_PER_SAMPLE is ", EXPECTED_CELLS_PER_SAMPLE,
          ", which looks like a rate rather than a cell count. It expects a ",
          "NUMBER OF CELLS, e.g. 20000 for a 20k-cell lane. To set a fixed ",
          "rate instead, use DOUBLET_RATE <- 0.16."))
      }
    }
    if (DOUBLET_RATE_FLOOR >= DOUBLET_RATE_CEILING) {
      errs <- c(errs, "DOUBLET_RATE_FLOOR must be less than DOUBLET_RATE_CEILING.")
    }
    if (DOUBLET_RATE_CEILING >= DOUBLET_ROLLBACK_THRESHOLD) {
      errs <- c(errs, paste0(
        "DOUBLET_RATE_CEILING (", DOUBLET_RATE_CEILING,
        ") must be LESS than DOUBLET_ROLLBACK_THRESHOLD (",
        DOUBLET_ROLLBACK_THRESHOLD, "), otherwise a densely loaded sample ",
        "would hit its own rollback and lose all doublet removal."))
    }
  } else if (is.numeric(DOUBLET_RATE) && DOUBLET_ROLLBACK_THRESHOLD <= DOUBLET_RATE) {
    # Rollback must sit above the target rate or DoubletFinder rolls back always.
    errs <- c(errs, paste0(
      "DOUBLET_ROLLBACK_THRESHOLD (", DOUBLET_ROLLBACK_THRESHOLD,
      ") must be GREATER than DOUBLET_RATE (", DOUBLET_RATE,
      "). DoubletFinder always calls ~DOUBLET_RATE doublets, so rollback ",
      "would trigger on every sample."))
  }
  if (length(errs) > 0) {
    stop(paste0("\n\n=== CONFIGURATION ERRORS (Section 1.4) ===\n  - ",
                paste(errs, collapse = "\n  - "), "\n"), call. = FALSE)
  }
  invisible(TRUE)
}
validate_config()

message(paste0("=== Doublet detection method: ", DOUBLET_METHOD,
               if (DOUBLET_METHOD == "both")
                 paste0(" (consensus rule: ", DOUBLET_CONSENSUS_RULE, ")") else "",
               " ==="))
if (REMOVE_DOUBLETS) {
  message("    Doublets will be REMOVED per sample. Downstream objects are singlet-only.")
} else {
  message("    Doublets will be LABELLED ONLY - no cells are deleted in Script 01.")
  message("    Filter later with: subset(data, subset = Doublet_Status == \"Singlet\")")
}

# =============================================================================
# --- DOUBLET RATE RESOLVER# --- DOUBLET RATE RESOLVER ---------------------------------------------------
# =============================================================================
# Returns the expected doublet rate for ONE sample.
#
#   DOUBLET_RATE == "auto"  -> estimate from this sample's own recovered cell
#                              count using the 10x linear multiplet model,
#                              clamped to [DOUBLET_RATE_FLOOR, ..._CEILING].
#   DOUBLET_RATE == number  -> that fixed number, for every sample.
#
# Called once per sample, immediately before doublet detection, using the cell
# count at that moment (post pre-QC, post-DecontX) - which is the best
# available proxy for "cells recovered" from the 10x lane.
resolve_doublet_rate <- function(n_cells, sample_id = "") {

  # Path 1: user forced a fixed fraction - use it verbatim.
  if (!identical(DOUBLET_RATE, "auto")) {
    message(paste0("    -> Doublet rate (fixed): ", round(DOUBLET_RATE * 100, 2),
                   "%  (expected ~", round(n_cells * DOUBLET_RATE), " doublets)"))
    return(list(rate = DOUBLET_RATE, source = "FIXED", raw = DOUBLET_RATE,
                basis_cells = n_cells))
  }

  # Path 2: "auto". Decide WHICH cell count feeds the 10x model.
  #   EXPECTED_CELLS_PER_SAMPLE = NULL   -> this sample's measured count
  #   EXPECTED_CELLS_PER_SAMPLE = number -> the stated loading target
  if (is.null(EXPECTED_CELLS_PER_SAMPLE)) {
    basis_cells <- n_cells
    basis_label <- "measured"
  } else {
    basis_cells <- EXPECTED_CELLS_PER_SAMPLE
    basis_label <- "user-stated target"
  }

  raw_rate <- (basis_cells / 1000) * MULTIPLET_RATE_PER_1K
  rate     <- min(max(raw_rate, DOUBLET_RATE_FLOOR), DOUBLET_RATE_CEILING)

  source <- if (raw_rate < DOUBLET_RATE_FLOOR) {
    "AUTO_FLOORED"
  } else if (raw_rate > DOUBLET_RATE_CEILING) {
    "AUTO_CAPPED"
  } else {
    "AUTO"
  }

  message(paste0("    -> Doublet rate (auto): ", basis_cells, " cells (",
                 basis_label, ") -> ", round(rate * 100, 2), "%",
                 if (source != "AUTO")
                   paste0("  [", source, "; unclamped estimate was ",
                          round(raw_rate * 100, 2), "%]") else "",
                 "  (expected ~", round(n_cells * rate), " doublets)"))

  if (source == "AUTO_CAPPED") {
    warning(paste0("  [RATE CAP] Sample '", sample_id, "': ", basis_cells,
                   " cells implies a ", round(raw_rate * 100, 1),
                   "% multiplet rate. Clamped to ",
                   round(DOUBLET_RATE_CEILING * 100, 1),
                   "%. Verify this lane was not overloaded."))
  }

  list(rate = rate, source = source, raw = raw_rate, basis_cells = basis_cells)
}

# =============================================================================
# --- DOUBLET DETECTION ENGINE ------------------------------------------------
# =============================================================================
# Two interchangeable detectors, each with an identical contract:
#
#   INPUT :  seu        — a per-sample Seurat object (post-DecontX, post-pre-QC)
#            sample_id  — character label, used for plot filenames
#            diag_dir   — directory to write diagnostic plots into
#            dbl_rate   — the expected doublet rate FOR THIS SAMPLE, already
#                         resolved by resolve_doublet_rate(). Passed in rather
#                         than read from the global, so that each sample can
#                         carry its own rate under DOUBLET_RATE = "auto".
#
#   OUTPUT:  a list with
#              $calls  — character vector ("Doublet"/"Singlet"), names = cell barcodes
#              $score  — numeric vector of doublet scores (pANN or scDblFinder score)
#              $info   — named list of method-specific diagnostics for the log
#
# Because both return the same shape, the caller does not care which ran.
# -----------------------------------------------------------------------------

# --- DETECTOR 1: DoubletFinder (paramSweep + constrained pK) -----------------
# Ported from eithan_coffee/01_process_data_SG.R, with the enhanced diagnostic
# plot from wu_project1/01_processing_samples.R and homotypic adjustment added.
run_doubletfinder <- function(seu, sample_id, diag_dir, dbl_rate) {

  # --- A) Preliminary processing required by DoubletFinder ------------------
  # DoubletFinder operates in PCA space, so the sample must be normalised,
  # scaled and PCA-reduced first. This is done on a THROWAWAY copy: none of
  # these reductions are kept, because the real clustering happens post-merge
  # on the Harmony-integrated object.
  n_cells <- ncol(seu)

  # Guard: PCA cannot request more PCs than cells-1. Small samples get fewer.
  npcs_use <- min(N_PCS_TO_USE, n_cells - 1)
  if (npcs_use < 5) {
    warning(paste0("  [SKIP] Sample '", sample_id, "' has only ", n_cells,
                   " cells - too few for DoubletFinder. All cells kept as Singlets."))
    return(list(
      calls = setNames(rep("Singlet", n_cells), colnames(seu)),
      score = setNames(rep(NA_real_, n_cells), colnames(seu)),
      info  = list(method = "DoubletFinder", pk = NA, pk_initial = NA,
                   nExp = 0, nExp_adj = 0, rate_used = dbl_rate,
                   homotypic_prop = NA, note = "TOO_FEW_CELLS")
    ))
  }

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = npcs_use, verbose = FALSE)

  # --- B) paramSweep: build the BCmetric curve over candidate pK values -----
  # paramSweep generates artificial doublets at many pN/pK combinations and
  # computes pANN distributions. summarizeSweep + find.pK then reduce this to
  # a single BCmetric (bimodality coefficient) per pK. The pK that maximises
  # BCmetric is, in principle, the one that best separates doublets from
  # singlets. In practice the curve is noisy at the extremes, hence step C.
  message(paste("    -> [DoubletFinder] Running paramSweep for", sample_id,
                "(this is the slow step)..."))
  sweep.res.list <- paramSweep(seu, PCs = 1:npcs_use, sct = DF_SCT)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)
  bcmvn$pK       <- as.numeric(as.character(bcmvn$pK))

  # --- C) Constrained pK selection ------------------------------------------
  # 1. Identify the global maximum first.
  initial_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  final_pk   <- initial_pk
  pk_note    <- "GLOBAL_PEAK_IN_RANGE"

  # 2. If it is outside the plausible window, look for the best in-window peak.
  if (final_pk < DF_PK_RANGE_MIN || final_pk > DF_PK_RANGE_MAX) {
    message(paste0("    -> [WARNING] Global optimal pK (", initial_pk,
                   ") is outside the plausible range [", DF_PK_RANGE_MIN,
                   " - ", DF_PK_RANGE_MAX, "]. Attempting constrained selection."))
    in_range <- bcmvn[bcmvn$pK >= DF_PK_RANGE_MIN & bcmvn$pK <= DF_PK_RANGE_MAX, ]

    if (nrow(in_range) > 0 && any(is.finite(in_range$BCmetric))) {
      final_pk <- in_range$pK[which.max(in_range$BCmetric)]
      pk_note  <- "CONSTRAINED_TO_RANGE"
      message(paste("      --> Using best in-range peak:", final_pk))
    } else {
      final_pk <- DF_PK_FALLBACK
      pk_note  <- "HARDCODED_FALLBACK"
      message(paste("      --> No usable peak in range. Reverting to fallback pK:", final_pk))
    }
  } else {
    message(paste0("    -> Optimal pK (", final_pk, ") is within the plausible range. Proceeding."))
  }

  # --- D) Diagnostic pK plot (enhanced version from wu_project1) ------------
  # Grey trace = full BCmetric curve; red dashed line = selected pK;
  # red diamond = the selected point itself (omitted if a fallback pK was used
  # that does not correspond to any swept value).
  tryCatch({
    selected_point_data <- bcmvn[bcmvn$pK == final_pk, ]
    pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
      geom_line(color = "grey60") +
      geom_point(color = "grey60") +
      annotate("rect", xmin = DF_PK_RANGE_MIN, xmax = DF_PK_RANGE_MAX,
               ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#4393C3") +
      geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
      { if (nrow(selected_point_data) > 0)
          geom_point(data = selected_point_data, aes(x = pK, y = BCmetric),
                     color = "red", size = 4, shape = 18)
      } +
      ggtitle(paste0("pK Finder: ", sample_id),
              subtitle = paste0("Selected pK = ", final_pk,
                                "  |  Global peak = ", initial_pk,
                                "  |  ", pk_note,
                                "\nShaded band = plausible range [",
                                DF_PK_RANGE_MIN, ", ", DF_PK_RANGE_MAX, "]")) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    ggsave(file.path(diag_dir, paste0(sample_id, "_DF_pk_plot.png")),
           plot = pk_plot, width = 7, height = 5, dpi = DPI_SETTING, bg = "white")
    rm(pk_plot)
  }, error = function(e) {
    message(paste("      -> [WARNING] pK plot failed:", e$message))
  })

  # --- E) nExp estimation, with optional homotypic adjustment ---------------
  nExp_raw       <- round(n_cells * dbl_rate)
  nExp_use       <- nExp_raw
  homotypic_prop <- NA_real_

  if (DF_HOMOTYPIC_ADJUST) {
    # modelHomotypic() needs cluster labels. A coarse clustering is enough:
    # it is only used to estimate what fraction of doublets would be
    # same-type (and therefore undetectable) given the cluster proportions.
    tryCatch({
      seu <- FindNeighbors(seu, dims = 1:npcs_use, verbose = FALSE)
      seu <- FindClusters(seu, resolution = DF_CLUSTER_RES, verbose = FALSE)
      homotypic_prop <- modelHomotypic(seu$seurat_clusters)
      nExp_use <- round(nExp_raw * (1 - homotypic_prop))
      message(paste0("    -> Homotypic adjustment: proportion = ",
                     round(homotypic_prop, 3), " | nExp ", nExp_raw,
                     " -> ", nExp_use))
    }, error = function(e) {
      message(paste("      -> [WARNING] Homotypic adjustment failed:", e$message,
                    "| Using unadjusted nExp =", nExp_raw))
    })
  }
  # Never let nExp collapse to 0 or exceed the cell count.
  nExp_use <- max(1, min(nExp_use, n_cells - 1))

  # --- F) Run DoubletFinder with the final, validated parameters -----------
  seu <- doubletFinder(seu,
                       PCs   = 1:npcs_use,
                       pN    = DF_PN,
                       pK    = final_pk,
                       nExp  = nExp_use,
                       sct   = DF_SCT)

  # DoubletFinder appends columns whose names embed the parameters used, e.g.
  # "DF.classifications_0.25_0.09_412" and "pANN_0.25_0.09_412". Grab the last
  # of each (the most recently added) rather than hardcoding the name.
  res_col  <- tail(grep("^DF.classifications", colnames(seu@meta.data), value = TRUE), 1)
  pann_col <- tail(grep("^pANN", colnames(seu@meta.data), value = TRUE), 1)

  calls <- setNames(as.character(seu@meta.data[[res_col]]), rownames(seu@meta.data))
  score <- if (length(pann_col) == 1 && !is.na(pann_col)) {
    setNames(as.numeric(seu@meta.data[[pann_col]]), rownames(seu@meta.data))
  } else {
    setNames(rep(NA_real_, n_cells), rownames(seu@meta.data))
  }

  list(
    calls = calls,
    score = score,
    info  = list(method = "DoubletFinder", pk = final_pk, pk_initial = initial_pk,
                 nExp = nExp_raw, nExp_adj = nExp_use, rate_used = dbl_rate,
                 homotypic_prop = homotypic_prop, note = pk_note)
  )
}

# --- DETECTOR 2: scDblFinder (automatic ML classifier) -----------------------
# Retained unchanged in behaviour from the 2026_nr4a1 v10.0 pipeline.
run_scdblfinder <- function(seu, sample_id, diag_dir, dbl_rate) {
  n_cells <- ncol(seu)

  sce_tmp <- as.SingleCellExperiment(seu)
  sce_tmp <- scDblFinder(
    sce_tmp,
    dbr     = if (SCDBLFINDER_DBR_NULL) NULL else dbl_rate,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 123),
    verbose = FALSE
  )

  score <- setNames(as.numeric(sce_tmp$scDblFinder.score), colnames(seu))
  calls <- setNames(
    ifelse(score >= SCDBLFINDER_SCORE_THRESHOLD, "Doublet", "Singlet"),
    colnames(seu)
  )

  # Diagnostic score-distribution histogram
  tryCatch({
    dblt_frac <- sum(calls == "Doublet") / n_cells
    score_df  <- data.frame(score = score, call = calls)
    p_score <- ggplot(score_df, aes(x = score, fill = call)) +
      geom_histogram(bins = 60, color = "black", linewidth = 0.2) +
      geom_vline(xintercept = SCDBLFINDER_SCORE_THRESHOLD,
                 linetype = "dashed", color = "red", linewidth = 0.8) +
      scale_fill_manual(values = c("Singlet" = "#4393C3", "Doublet" = "#D6604D")) +
      labs(title    = paste0("scDblFinder Score Distribution: ", sample_id),
           subtitle = paste0(round(dblt_frac * 100, 1),
                             "% doublets | threshold = ", SCDBLFINDER_SCORE_THRESHOLD),
           x = "Doublet Score", y = "Cell Count", fill = "Call") +
      theme_minimal()
    ggsave(file.path(diag_dir, paste0(sample_id, "_scDbl_score_distribution.png")),
           plot = p_score, width = 7, height = 5, dpi = DPI_SETTING, bg = "white")
    rm(p_score, score_df)
  }, error = function(e) {
    message(paste("      -> [WARNING] Score plot failed:", e$message))
  })

  rm(sce_tmp); gc()

  list(
    calls = calls,
    score = score,
    info  = list(method = "scDblFinder", pk = NA, pk_initial = NA,
                 nExp = NA, nExp_adj = NA, rate_used = dbl_rate,
                 homotypic_prop = NA,
                 note = paste0("threshold=", SCDBLFINDER_SCORE_THRESHOLD))
  )
}

# --- CONCORDANCE REPORTER (used only when DOUBLET_METHOD == "both") ----------
# Writes a 2x2 confusion matrix comparing the two callers and returns summary
# statistics (agreement rate, Cohen's kappa) for the master log.
report_concordance <- function(df_calls, sc_calls, sample_id, diag_dir) {
  common <- intersect(names(df_calls), names(sc_calls))
  a <- factor(df_calls[common], levels = c("Singlet", "Doublet"))
  b <- factor(sc_calls[common], levels = c("Singlet", "Doublet"))
  tab <- table(DoubletFinder = a, scDblFinder = b)

  n         <- length(common)
  agree     <- sum(diag(tab)) / n
  # Cohen's kappa: agreement corrected for chance.
  p_exp     <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  kappa     <- if (p_exp < 1) (agree - p_exp) / (1 - p_exp) else NA_real_

  tryCatch({
    tab_df <- as.data.frame(tab)
    p_conc <- ggplot(tab_df, aes(x = scDblFinder, y = DoubletFinder, fill = Freq)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = Freq), size = 6, fontface = "bold") +
      scale_fill_gradient(low = "#F7F7F7", high = "#4393C3") +
      labs(title    = paste0("Doublet Call Concordance: ", sample_id),
           subtitle = paste0("Agreement = ", round(agree * 100, 1),
                             "%  |  Cohen's kappa = ", round(kappa, 3),
                             "  |  n = ", n, " cells"),
           x = "scDblFinder call", y = "DoubletFinder call") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    ggsave(file.path(diag_dir, paste0(sample_id, "_doublet_concordance.png")),
           plot = p_conc, width = 6, height = 5, dpi = DPI_SETTING, bg = "white")
    rm(p_conc)
  }, error = function(e) {
    message(paste("      -> [WARNING] Concordance plot failed:", e$message))
  })

  message(paste0("    -> Concordance: ", round(agree * 100, 1),
                 "% agreement | kappa = ", round(kappa, 3)))

  list(agreement = agree, kappa = kappa,
       n_df_only = tab["Doublet", "Singlet"],
       n_sc_only = tab["Singlet", "Doublet"],
       n_both    = tab["Doublet", "Doublet"])
}

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
#     5. Doublet detection — per-sample, via DOUBLET_METHOD, then removal
#     6. Pre-merge QC (loose filters) — applied on corrected, singlet-only counts
#     7. Save checkpoint and free memory
#
# DUAL-CHECKPOINT SYSTEM (three-tier resume logic):
#   Each sample produces TWO checkpoint files inside PROCESSED_DIR/<SampleID>/:
#
#     Checkpoint 1 (_scevan_processed.rds):
#       Saved after SCEVAN runs. Contains raw Seurat object + SCEVAN metadata.
#       Used as a mid-pipeline recovery point so that expensive SCEVAN
#       computation is never repeated.
#
#     Checkpoint 2 (_decontx_dblt_processed.rds):
#       Saved after DecontX + doublet removal + Pre-QC. This is the FINAL
#       per-sample checkpoint and the file loaded for the merge step (2.1b).
#       It contains clean, singlet-only, ambient-corrected cells.
#
#   At the start of each iteration the loop checks:
#     (a) If Checkpoint 2 exists  → sample fully processed, SKIP entirely.
#     (b) If only Checkpoint 1 exists → load it, SKIP SCEVAN, run DecontX
#         + doublet detection + Pre-QC, then save Checkpoint 2.
#     (c) If neither exists → run the entire pipeline from scratch.
#
#   This means that if a run crashes mid-way through DecontX or doublet detection,
#   the next run will resume from the SCEVAN result rather than starting over.
#
# MEMORY STRATEGY:
#   Each sample object is removed from RAM and garbage-collected immediately
#   after its checkpoint is written. Only one sample occupies memory at a time.
# =============================================================================
message(paste0("=== STEP 2.1a: Per-Sample Processing (SCEVAN → Pre-QC → DecontX → ",
               DOUBLET_METHOD, " → Checkpoint) ==="))
metadata <- read.xlsx(METADATA_FILE)

# --- Clean up header artifacts -----------------------------------------------
# Excel/openxlsx headers frequently arrive dirty:
#   - read.xlsx() converts spaces to dots  ("Sample ID" -> "Sample.ID")
#   - headers pasted from web sources carry XML/HTML fragments
#     ('xml:space="preserve">Probes', 'Tissue.Type.&amp;.species')
# Strip those artifacts so column matching is not defeated by cosmetics.
clean_header <- function(x) {
  x <- gsub('xml:space="preserve">', "", x, fixed = TRUE)  # leaked XML attr
  x <- gsub("&amp;", "and", x, fixed = TRUE)               # HTML ampersand
  x <- gsub("&[a-z]+;", "", x)                             # any other entity
  x <- gsub("[.]+", " ", x)                                # dots back to spaces
  x <- trimws(gsub("\\s+", " ", x))                        # collapse whitespace
  x
}
colnames(metadata) <- clean_header(colnames(metadata))

# --- Resolve the SampleID column, tolerantly ---------------------------------
# Accept common spellings/casings (SampleID, Sample ID, Sample_ID, sampleid ...)
# by comparing on a normalised key: lowercase, alphanumerics only.
norm_key <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
sid_hits <- which(norm_key(colnames(metadata)) == "sampleid")

if (length(sid_hits) == 0) {
  stop(paste0(
    "\n\nNo SampleID column found in the metadata file.\n",
    "  File   : ", METADATA_FILE, "\n",
    "  Columns: ", paste(colnames(metadata), collapse = " | "), "\n",
    "Name the sample-identifier column 'SampleID' (or 'Sample ID' / ",
    "'Sample_ID' - the pipeline normalises those) and re-run."), call. = FALSE)
}
if (length(sid_hits) > 1) {
  stop(paste0("\n\nMultiple columns look like SampleID: ",
              paste(colnames(metadata)[sid_hits], collapse = ", "),
              ".\nKeep only one and re-run."), call. = FALSE)
}
if (colnames(metadata)[sid_hits] != "SampleID") {
  message(paste0("  [METADATA] Using column '", colnames(metadata)[sid_hits],
                 "' as SampleID."))
  colnames(metadata)[sid_hits] <- "SampleID"
}

metadata$SampleID <- trimws(as.character(metadata$SampleID))

bad_ids <- is.na(metadata$SampleID) | metadata$SampleID == "" | metadata$SampleID == "NA"
if (any(bad_ids)) {
  stop(paste0("\n\nThe SampleID column has ", sum(bad_ids),
              " empty/NA value(s) in row(s): ",
              paste(which(bad_ids), collapse = ", "),
              ".\nEvery row needs a SampleID that matches its H5 subfolder."),
       call. = FALSE)
}

# Confirm each SampleID has a matching H5 folder before doing any work.
missing_dirs <- metadata$SampleID[!dir.exists(file.path(H5_DIR, metadata$SampleID))]
if (length(missing_dirs) > 0) {
  message(paste0("  [WARNING] ", length(missing_dirs),
                 " SampleID(s) have no folder inside H5_DIR:"))
  message(paste0("    ", paste(missing_dirs, collapse = ", ")))
  message(paste0("  H5_DIR contains: ",
                 paste(head(list.dirs(H5_DIR, recursive = FALSE, full.names = FALSE), 20),
                       collapse = ", ")))
  message("  These samples will fail to load. Fix SampleID <-> folder-name mismatches.")
}
message(paste0("  Metadata OK: ", nrow(metadata), " samples, SampleID column present."))

# Initialize rollback log (filled during loop, saved after)
rollback_log <- list()

for (i in 1:nrow(metadata)) {
  sample_info       <- metadata[i, ]
  sample_id         <- as.character(sample_info$SampleID)
  sample_dir <- file.path(PROCESSED_DIR, sample_id)
  checkpoint_file_1 <- file.path(sample_dir, paste0(sample_id, "_scevan_processed.rds"))
  checkpoint_file_2 <- file.path(sample_dir, paste0(sample_id, "_decontx_dblt_processed.rds"))

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
                  "- will run DecontX + doublet detection only."))
    seurat_obj <- readRDS(checkpoint_file_1)
    scevan_checkpoint_loaded <- TRUE
  }
  # (c) No checkpoint: fall through to full processing below.

  message(paste("  [PROCESSING] Starting", sample_id,
                if (scevan_checkpoint_loaded) "(resuming from SCEVAN checkpoint)" else "(full run)",
                "..."))
  if (!dir.exists(sample_dir)) { dir.create(sample_dir, recursive = TRUE) }

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

    # ---- Attach ALL metadata columns from the Excel file --------------------
    # Every column of this sample's metadata row is broadcast to every cell
    # from this sample. sample_info[[col_name]] is a length-1 value, so R
    # recycles it across all cells. This is deliberate: any column you add to
    # the .xlsx becomes available in Scripts 02-10 with no code changes.
    #
    # Two guards are applied before assignment:
    #   1. Reserved names — a metadata column called e.g. "nCount_RNA" would
    #      silently overwrite the QC metric Seurat computed, corrupting every
    #      downstream filter. Such columns are skipped with a warning.
    #   2. Non-syntactic names — a header like "Pool ID" becomes a metadata
    #      column that needs backticks in R (data$`Pool ID`). Allowed, but
    #      flagged once so it is a deliberate choice rather than a surprise.
    reserved_names <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                        "Doublet_Status", "DF_score", "DF_class",
                        "scDblFinder_score", "scDblFinder_class",
                        "seurat_clusters", "CellType", "sub_cell_types")

    for (col_name in colnames(sample_info)) {
      if (col_name %in% reserved_names) {
        warning(paste0(
          "  [METADATA] Column '", col_name, "' in the metadata file collides ",
          "with a reserved pipeline/Seurat column and was SKIPPED for sample '",
          sample_id, "'. Rename it in the .xlsx to keep its values."))
        next
      }
      if (col_name != make.names(col_name)) {
        message(paste0("      -> [NOTE] Metadata column '", col_name,
                       "' is not a syntactic R name; refer to it with backticks",
                       " downstream (e.g. data$`", col_name, "`)."))
      }
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
      # SCEVAN writes several files under a FIXED "output" folder relative to
      # the current working directory, and does not fully honour output_dir.
      # Left alone it scatters those files into whatever directory R happens to
      # be in, and successive samples overwrite each other. We therefore run it
      # with the working directory temporarily set to this sample's own scevan/
      # subfolder, and restore the original afterwards no matter what happens.
      scevan_dir <- file.path(sample_dir, "scevan")
      if (!dir.exists(scevan_dir)) dir.create(scevan_dir, recursive = TRUE)
      .old_wd <- getwd()
      tryCatch({
        setwd(scevan_dir)
        message("    -> Running SCEVAN for CNA detection...")
        scevan_res_df <- pipelineCNA(
          as.matrix(GetAssayData(seurat_obj, layer = "counts")),
          par_cores  = SCEVAN_N_CORES,
          SUBCLONES  = SCEVAN_SUBCLONES,
          plotTree   = SCEVAN_PLOTTREE,
          organism   = SCEVAN_ORGANISM,
          output_dir = scevan_dir
        )
        if (!is.null(scevan_res_df)) {
          seurat_obj <- AddMetaData(seurat_obj, metadata = scevan_res_df)
          rm(scevan_res_df)
          message("    -> SCEVAN metadata added to Seurat object.")
        }
      }, error = function(e) {
        message(paste("    -> [WARNING] SCEVAN failed for", sample_id, "| Error:", e$message))
      }, finally = {
        # ALWAYS restore the working directory - a stranded setwd() would send
        # every later output of this run into the sample's scevan folder.
        setwd(.old_wd)
      })
      gc()  # SCEVAN allocates heavily internally; reclaim memory immediately
    }

    # ---- Save Checkpoint 1 (SCEVAN result) ---------------------------------
    # Saved immediately after SCEVAN so that if DecontX or doublet detection crash
    # later, the next run can resume from here and skip the expensive SCEVAN step.
    message("    -> Saving SCEVAN checkpoint (checkpoint 1)...")
    saveRDS(seurat_obj, file = checkpoint_file_1)
    gc()
    message("    -> Checkpoint 1 saved.")

  } # end if (!scevan_checkpoint_loaded)

  # ---- 5. Pre-merge QC: loose per-sample filters ---------------------------
  # Applied BEFORE DecontX and doublet detection so that empty droplets and
  # clearly dying cells are removed first, consistent with the recommendation
  # in Germain et al. (scDblFinder) and McGinnis et al. (DoubletFinder): the
  # object should not contain empty drops
  # but should not otherwise have undergone very stringent filtering.
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA >= PRE_MIN_GENES_PER_CELL &
                                percent.mt   <= PRE_MAX_MT_PERCENT)
  message(paste0("      -> Pre-merge QC: ", ncol(seurat_obj), " cells retained."))

  # ---- 6. Optional: DecontX ambient RNA correction -------------------------
  # DecontX is run per-sample BEFORE doublet detection so that the doublet
  # caller operates on counts free of ambient RNA contamination.
  if (RUN_DECONTX) {
    tryCatch({
      message(paste("    -> Running DecontX for ambient RNA correction on", sample_id, "..."))
      counts_sparse   <- GetAssayData(object = seurat_obj, layer = "counts")
      decontx_results <- decontX(x = counts_sparse)

      # Replace raw counts with decontaminated counts.
      # Optionally round to integer (ROUND_DECONTX_COUNTS = TRUE) — useful for
      # tools that require integer count matrices. See decontX vignette for details.
      if (ROUND_DECONTX_COUNTS) {
        # round() on a sparse matrix leaves explicit 0 entries for values that
        # rounded down to 0 (e.g. 0.3 → 0). drop0() removes those structural
        # zeros so the matrix stays truly sparse.
        seurat_obj[["RNA"]]$counts <- Matrix::drop0(round(decontx_results$decontXcounts))
        message("    -> DecontX counts rounded to nearest integer and structural zeros dropped.")
      } else {
        # drop0() also applied here: decontX can emit near-zero fractional values
        # that are structurally non-zero, wasting memory.
        seurat_obj[["RNA"]]$counts <- Matrix::drop0(decontx_results$decontXcounts)
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

  # ---- 7. Doublet detection and removal (method set by DOUBLET_METHOD) ------
  # Dispatches to run_doubletfinder() and/or run_scdblfinder(), both defined
  # near the top of PART 2. Whichever runs, the result is normalised into a
  # single Doublet_Status column so downstream scripts stay method-agnostic.
  #
  # Runs AFTER DecontX so that doublet calling sees ambient-corrected counts,
  # and after the loose pre-QC so that empty drops are gone - but no stringent
  # filtering has occurred, which would bias the doublet-rate estimate.
  message(paste0("    -> Running doublet detection (", DOUBLET_METHOD, ") for ", sample_id, "..."))

  # Resolve THIS sample's expected doublet rate before anything else. Under
  # DOUBLET_RATE = "auto" this is derived from the sample's own recovered cell
  # count, so a shallow sample and a deep sample are no longer forced to share
  # a single (and for one of them, wrong) quota.
  rate_info <- resolve_doublet_rate(ncol(seurat_obj), sample_id)

  tryCatch({

    df_res <- NULL   # DoubletFinder result
    sc_res <- NULL   # scDblFinder result
    conc   <- NULL   # concordance stats (only when method == "both")

    # ---- 7a. Run the requested detector(s) ---------------------------------
    if (DOUBLET_METHOD %in% c("DoubletFinder", "both")) {
      df_res <- run_doubletfinder(seurat_obj, sample_id, DOUBLET_DIAG_DIR, rate_info$rate)
      seurat_obj$DF_score  <- df_res$score[colnames(seurat_obj)]
      seurat_obj$DF_class  <- df_res$calls[colnames(seurat_obj)]
    }
    if (DOUBLET_METHOD %in% c("scDblFinder", "both")) {
      sc_res <- run_scdblfinder(seurat_obj, sample_id, DOUBLET_DIAG_DIR, rate_info$rate)
      seurat_obj$scDblFinder_score <- sc_res$score[colnames(seurat_obj)]
      seurat_obj$scDblFinder_class <- sc_res$calls[colnames(seurat_obj)]
    }

    # ---- 7b. Resolve the final call ----------------------------------------
    if (DOUBLET_METHOD == "DoubletFinder") {
      final_calls  <- df_res$calls[colnames(seurat_obj)]
      primary_info <- df_res$info

    } else if (DOUBLET_METHOD == "scDblFinder") {
      final_calls  <- sc_res$calls[colnames(seurat_obj)]
      primary_info <- sc_res$info

    } else {  # "both" - compare, log, then apply the consensus rule
      conc <- report_concordance(df_res$calls, sc_res$calls, sample_id, DOUBLET_DIAG_DIR)
      a <- df_res$calls[colnames(seurat_obj)]
      b <- sc_res$calls[colnames(seurat_obj)]
      final_calls <- switch(
        DOUBLET_CONSENSUS_RULE,
        # union: flagged by either method (aggressive, highest recall)
        "union"         = ifelse(a == "Doublet" | b == "Doublet", "Doublet", "Singlet"),
        # intersect: flagged by both methods (conservative, highest precision)
        "intersect"     = ifelse(a == "Doublet" & b == "Doublet", "Doublet", "Singlet"),
        # one method decides; the other is still recorded in metadata
        "DoubletFinder" = a,
        "scDblFinder"   = b
      )
      names(final_calls) <- colnames(seurat_obj)
      primary_info  <- df_res$info
    }

    seurat_obj$Doublet_Status <- unname(final_calls)

    n_doublets       <- sum(seurat_obj$Doublet_Status == "Doublet")
    doublet_fraction <- n_doublets / ncol(seurat_obj)

    # ---- 7c. Rollback check -------------------------------------------------
    # Guard against biologically implausible doublet fractions, which almost
    # always indicate a classification artifact rather than a real result.
    if (doublet_fraction > DOUBLET_ROLLBACK_THRESHOLD) {
      warning(paste0(
        "  [ROLLBACK] Sample '", sample_id, "': detected doublet fraction = ",
        round(doublet_fraction * 100, 1), "% exceeds threshold of ",
        round(DOUBLET_ROLLBACK_THRESHOLD * 100, 1), "%. ",
        "Doublet removal SKIPPED for this sample (all cells treated as Singlets)."
      ))
      seurat_obj$Doublet_Status <- "Singlet"
      action_taken <- "ROLLBACK_APPLIED"
    } else {
      action_taken <- "DOUBLETS_REMOVED"
    }

    rollback_log[[sample_id]] <- list(
      method         = DOUBLET_METHOD,
      consensus      = if (DOUBLET_METHOD == "both") DOUBLET_CONSENSUS_RULE else NA_character_,
      disposition    = if (REMOVE_DOUBLETS) "REMOVED" else "LABELLED_ONLY",
      rate_used      = rate_info$rate,
      rate_source    = rate_info$source,
      fraction       = doublet_fraction,
      n_doublets     = n_doublets,
      n_cells        = ncol(seurat_obj),
      pk             = primary_info$pk,
      pk_initial     = primary_info$pk_initial,
      nExp           = primary_info$nExp,
      nExp_adj       = primary_info$nExp_adj,
      homotypic_prop = primary_info$homotypic_prop,
      note           = primary_info$note,
      agreement      = if (!is.null(conc)) conc$agreement else NA_real_,
      kappa          = if (!is.null(conc)) conc$kappa     else NA_real_,
      action         = action_taken
    )

    message(paste0("      -> Doublet fraction: ",
                   round(doublet_fraction * 100, 1), "%",
                   if (!is.na(primary_info$pk)) paste0(" | pK: ", primary_info$pk) else "",
                   " | Action: ", action_taken))

    # ---- 7d. Diagnostic UMAPs: doublet scores before removal ---------------
    # One UMAP is computed and reused for every panel below.
    #
    # IMPORTANT: each method gets its OWN plot. The two scores are NOT on the
    # same scale and are never merged:
    #   DF_score          = pANN, the proportion of artificial nearest
    #                       neighbours. Range depends on pK and pN.
    #   scDblFinder_score = a calibrated classifier probability in [0, 1].
    # Averaging them would be arithmetic on incompatible units, so the
    # pipeline plots them side by side and lets you judge.
    tryCatch({
      seu_tmp_umap <- seurat_obj
      n_pc_umap    <- min(30, ncol(seu_tmp_umap) - 1)
      seu_tmp_umap <- NormalizeData(seu_tmp_umap, verbose = FALSE) %>%
        FindVariableFeatures(verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = n_pc_umap, verbose = FALSE) %>%
        RunUMAP(dims = 1:n_pc_umap, verbose = FALSE,
                reduction.name = "umap_doublets")

      sub_txt <- paste0(round(doublet_fraction * 100, 1),
                        "% doublets | Action: ", action_taken)

      save_score_umap <- function(col, label, fname) {
        if (!col %in% colnames(seu_tmp_umap@meta.data)) return(invisible(NULL))
        if (all(is.na(seu_tmp_umap@meta.data[[col]]))) return(invisible(NULL))
        p <- FeaturePlot(seu_tmp_umap, features = col,
                         reduction = "umap_doublets") +
          coord_fixed() +
          ggtitle(paste0(label, ": ", sample_id), subtitle = sub_txt)
        ggsave(file.path(DOUBLET_DIAG_DIR, paste0(sample_id, "_", fname, ".png")),
               plot = p, width = 7, height = 6, dpi = DPI_SETTING, bg = "white")
        rm(p)
      }

      # Per-method score maps - only for methods that actually ran.
      save_score_umap("DF_score",
                      "DoubletFinder pANN score", "doublet_umap_DoubletFinder")
      save_score_umap("scDblFinder_score",
                      "scDblFinder probability",  "doublet_umap_scDblFinder")

      # The final consensus call, whatever rule produced it.
      p_call <- DimPlot(seu_tmp_umap, group.by = "Doublet_Status",
                        reduction = "umap_doublets",
                        cols = c("Singlet" = "grey80", "Doublet" = "#D6604D")) +
        coord_fixed() +
        ggtitle(paste0("Final call (", DOUBLET_METHOD,
                       if (DOUBLET_METHOD == "both")
                         paste0("/", DOUBLET_CONSENSUS_RULE) else "",
                       "): ", sample_id),
                subtitle = sub_txt)
      ggsave(file.path(DOUBLET_DIAG_DIR, paste0(sample_id, "_doublet_umap_call.png")),
             plot = p_call, width = 7, height = 6, dpi = DPI_SETTING, bg = "white")
      rm(p_call)

      # METHOD-CLASSIFICATION MAP (both-mode only). The most informative panel:
      # every cell falls into one of four states - agreed doublet, agreed
      # singlet, or flagged by exactly one method.
      # cells the two methods disagree about. If they scatter randomly the
      # callers are merely noisy at the margin; if they concentrate on one
      # cluster, one method is systematically wrong about a real population -
      # which the agreement percentage and kappa alone will never reveal.
      if (DOUBLET_METHOD == "both" &&
          all(c("DF_class", "scDblFinder_class") %in% colnames(seu_tmp_umap@meta.data))) {
        seu_tmp_umap$.dblt_agreement <- ifelse(
          seu_tmp_umap$DF_class == "Doublet" & seu_tmp_umap$scDblFinder_class == "Doublet",
          "Both: doublet",
          ifelse(seu_tmp_umap$DF_class == "Doublet", "DoubletFinder only",
          ifelse(seu_tmp_umap$scDblFinder_class == "Doublet", "scDblFinder only",
                 "Both: singlet")))
        p_dis <- DimPlot(seu_tmp_umap, group.by = ".dblt_agreement",
                         reduction = "umap_doublets",
                         cols = c("Both: singlet"      = "grey85",
                                  "Both: doublet"      = "#B2182B",
                                  "DoubletFinder only" = "#2166AC",
                                  "scDblFinder only"   = "#F4A582")) +
          coord_fixed() +
          ggtitle(paste0("Doublet classification by method: ", sample_id),
                  subtitle = paste0("agreement ",
                                    if (!is.null(conc)) round(conc$agreement * 100, 1) else NA,
                                    "% | kappa ",
                                    if (!is.null(conc)) round(conc$kappa, 3) else NA))
        ggsave(file.path(DOUBLET_DIAG_DIR,
                         paste0(sample_id, "_doublet_umap_method_classification.png")),
               plot = p_dis, width = 7.5, height = 6, dpi = DPI_SETTING, bg = "white")
        rm(p_dis)
      }

      rm(seu_tmp_umap)
    }, error = function(e) {
      message(paste("      -> [WARNING] Doublet UMAP failed for", sample_id, ":", e$message))
    })

    # ---- 7e. Remove doublets, or keep them labelled ------------------------
    # Controlled by REMOVE_DOUBLETS (Section 1.4.4).
    #   TRUE  - drop them here; every downstream object is singlet-only
    #   FALSE - keep every cell and only LABEL it. Nothing is deleted; you
    #           decide when (and whether) to filter, e.g. after annotation so
    #           you can first see where the flagged cells actually sit.
    cells_before_dblt <- ncol(seurat_obj)
    if (REMOVE_DOUBLETS) {
      seurat_obj <- subset(seurat_obj, subset = Doublet_Status == "Singlet")
      message(paste0("      -> Cells: ", cells_before_dblt, " -> ", ncol(seurat_obj),
                     " (REMOVED ", cells_before_dblt - ncol(seurat_obj), " doublets)"))
    } else {
      n_flagged <- sum(seurat_obj$Doublet_Status == "Doublet")
      message(paste0("      -> Cells: ", cells_before_dblt, " retained (",
                     n_flagged, " LABELLED as doublets, none removed)."))
      message("         REMOVE_DOUBLETS = FALSE - filter later with:")
      message("           data <- subset(data, subset = Doublet_Status == \"Singlet\")")
      rm(n_flagged)
    }

    rm(df_res, sc_res, conc, n_doublets, doublet_fraction, cells_before_dblt)
    gc()

  }, error = function(e) {
    message(paste("    -> [WARNING] Doublet detection failed for", sample_id,
                  "| Error:", e$message,
                  "| Proceeding with NO doublet removal for this sample."))
    # NOTE: `<<-` is REQUIRED here, not `<-`.
    # An error handler passed to tryCatch is a *function*, so it gets its own
    # evaluation frame. A plain `<-` would create a local copy that is discarded
    # the instant the handler returns, meaning the sample would be checkpointed
    # with NO Doublet_Status column and NO entry in the log - a silent data
    # integrity failure. `<<-` assigns to the enclosing (script) environment.
    # (The success branch above does not need this: tryCatch evaluates its
    # main expression in the caller's frame, not a new one.)
    seurat_obj$Doublet_Status <<- "Singlet"
    rollback_log[[sample_id]] <<- list(
      method = DOUBLET_METHOD, consensus = NA_character_,
      rate_used = rate_info$rate, rate_source = rate_info$source,
      fraction = NA_real_, n_doublets = NA_integer_, n_cells = ncol(seurat_obj),
      pk = NA_real_, pk_initial = NA_real_, nExp = NA_integer_,
      nExp_adj = NA_integer_, homotypic_prop = NA_real_,
      note = paste("ERROR:", e$message),
      agreement = NA_real_, kappa = NA_real_,
      action = "ERROR_SKIPPED"
    )
  })

  # ---- 8. Save Checkpoint 2 (decontX + doublets + QC) and free memory ------
  # This is the FINAL per-sample checkpoint and the file Step 2.1b loads for
  # merging. Its contents depend on REMOVE_DOUBLETS: singlet-only when TRUE,
  # or all cells with a Doublet_Status label when FALSE.
  # Saving here means a crash AFTER this point (e.g., during merge) will not
  # require re-running DecontX or doublet detection on the next attempt.
  message("    -> Saving full checkpoint (decontX + doublets + QC) and releasing memory...")
  saveRDS(seurat_obj, file = checkpoint_file_2)
  rm(seurat_obj)
  gc()
  message(paste("  [DONE]", sample_id, "- checkpoint 2 saved, memory freed."))
}

# --- Save the master doublet log (covers all processed samples) --------------
# One row per sample recording: which method ran, the pK actually used and
# whether it had to be constrained, the homotypic adjustment, the resulting
# doublet fraction, whether rollback fired, and (for method "both") the
# agreement between the two callers. This file is the audit trail for the
# doublet step and should be reviewed before trusting any downstream result.
if (length(rollback_log) > 0) {
  # Helper: pull a field, returning NA if it is missing or NULL, so that a
  # partially-populated log entry (e.g. from an error path) cannot break rbind.
  .fld <- function(s, f, default = NA) {
    v <- rollback_log[[s]][[f]]
    if (is.null(v) || length(v) == 0) default else v[1]
  }

  rollback_df <- do.call(rbind, lapply(names(rollback_log), function(s) {
    data.frame(
      SampleID           = s,
      Method             = .fld(s, "method", NA_character_),
      Disposition        = .fld(s, "disposition", NA_character_),
      Consensus_Rule     = .fld(s, "consensus", NA_character_),
      N_Cells            = .fld(s, "n_cells"),
      N_Doublets_Flagged = .fld(s, "n_doublets"),
      Doublet_Fraction   = round(as.numeric(.fld(s, "fraction")) * 100, 2),
      Target_Rate_Pct    = round(as.numeric(.fld(s, "rate_used")) * 100, 2),
      Rate_Source        = .fld(s, "rate_source", NA_character_),
      Rollback_Thresh    = DOUBLET_ROLLBACK_THRESHOLD * 100,
      pK_Used            = .fld(s, "pk"),
      pK_Global_Peak     = .fld(s, "pk_initial"),
      nExp_Raw           = .fld(s, "nExp"),
      nExp_Homotypic_Adj = .fld(s, "nExp_adj"),
      Homotypic_Prop     = round(as.numeric(.fld(s, "homotypic_prop")), 4),
      Method_Agreement   = round(as.numeric(.fld(s, "agreement")) * 100, 2),
      Cohens_Kappa       = round(as.numeric(.fld(s, "kappa")), 4),
      Note               = as.character(.fld(s, "note", NA_character_)),
      Action             = .fld(s, "action", NA_character_),
      stringsAsFactors   = FALSE
    )
  }))
  log_name <- paste0("doublet_detection_log_", DOUBLET_METHOD, ".xlsx")
  write_xlsx(rollback_df, file.path(OUTPUT_DIR, log_name))
  message(paste("  Doublet summary saved to:", log_name))

  # Console summary so problems are visible without opening the Excel file.
  n_rollback <- sum(rollback_df$Action == "ROLLBACK_APPLIED", na.rm = TRUE)
  n_error    <- sum(rollback_df$Action == "ERROR_SKIPPED",    na.rm = TRUE)
  n_fallback <- sum(rollback_df$Note %in% c("HARDCODED_FALLBACK"), na.rm = TRUE)
  message(paste0("  Doublet step summary: ", nrow(rollback_df), " samples | ",
                 n_rollback, " rollback(s) | ", n_error, " error(s) | ",
                 n_fallback, " fallback pK"))
  if (identical(DOUBLET_RATE, "auto")) {
    rr <- suppressWarnings(as.numeric(rollback_df$Target_Rate_Pct))
    rr <- rr[is.finite(rr)]
    if (length(rr) > 0) {
      message(paste0("  Auto doublet rate across samples: ",
                     round(min(rr), 2), "% - ", round(max(rr), 2),
                     "% (median ", round(stats::median(rr), 2), "%)"))
    }
    n_capped <- sum(rollback_df$Rate_Source == "AUTO_CAPPED", na.rm = TRUE)
    if (n_capped > 0) {
      message(paste0("  [ATTENTION] ", n_capped,
                     " sample(s) hit the rate cap (DOUBLET_RATE_CEILING). ",
                     "Check for overloaded lanes."))
    }
  }
  if (n_rollback > 0 || n_error > 0) {
    message("  [ATTENTION] Review the doublet log before trusting downstream results.")
  }
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
  sample_dir <- file.path(PROCESSED_DIR, sample_id)
  checkpoint_file_2 <- file.path(sample_dir, paste0(sample_id, "_decontx_dblt_processed.rds"))
  checkpoint_file_1 <- file.path(sample_dir, paste0(sample_id, "_scevan_processed.rds"))

  if (file.exists(checkpoint_file_2)) {
    # Normal case: full checkpoint (decontX + doublets + QC) is present.
    message(paste("  [LOAD]", sample_id, "(checkpoint 2 — decontX + doublets + QC)"))
    seurat_objects_list[[sample_id]] <- readRDS(checkpoint_file_2)
  } else if (file.exists(checkpoint_file_1)) {
    # Fallback: only the SCEVAN checkpoint exists (e.g., pipeline from a
    # previous version, or Step 2.1a crashed before saving checkpoint 2).
    # NOTE: this object has NOT had DecontX or doublet removal applied — re-run
    # Step 2.1a first to produce the complete checkpoint 2 before merging.
    warning(paste0(
      "  [FALLBACK] Sample '", sample_id, "': only SCEVAN checkpoint found. ",
      "DecontX and doublet removal may NOT have been applied. ",
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
# Two interchangeable routes, chosen by INTEGRATION_METHOD (Section 1.6). Both
# end with the SAME reductions - "pca", "umap_none"/clusters_none (Track A,
# unintegrated diagnostic) and "harmony", "umap_harmony"/clusters_harmony
# (Track B, batch-corrected) - and a joined object, so Script 02 does not care
# which route ran.
#
#   "RunHarmony"      Work on the JOINED matrix throughout. Normalize -> HVG ->
#                     Scale -> PCA once, then harmony::RunHarmony() on the PCA.
#                     No split, no re-join. Cheaper; the recommended default.
#
#   "IntegrateLayers" Seurat v5 layered path: split by SampleID, process per
#                     layer, IntegrateLayers(HarmonyIntegration), then JoinLayers.
#                     Retained for reproducibility with older runs.
# =============================================================================
output_rds_layered <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_merged_post_qc.rds"))
data <- readRDS(output_rds_layered)
data[["RNA"]]$data       <- NULL
data[["RNA"]]$scale.data <- NULL
gc()

message(paste0("\n=== STEP 2.4: Normalization, Dimensionality Reduction, and ",
               "Integration (", INTEGRATION_METHOD, ") ==="))

if (INTEGRATION_METHOD == "RunHarmony") {
  # ---------------------------------------------------------------------------
  # ROUTE A: harmony::RunHarmony() on the joined object. No layer splitting.
  # The object reloaded above was already joined in Step 2.3, so there is
  # nothing to split or re-join here.
  # ---------------------------------------------------------------------------
  suppressPackageStartupMessages(library(harmony))

  data <- NormalizeData(data, verbose = TRUE)
  data <- FindVariableFeatures(data, nfeatures = N_VARIABLE_FEATURES, verbose = TRUE)
  data <- ScaleData(data, verbose = TRUE)
  data <- RunPCA(data, npcs = N_PCS_TO_USE, verbose = TRUE)
  message("  -> Normalization, scaling and PCA complete (joined matrix).")
  gc()

  # --- Track A: unintegrated PCA (diagnostic) ---
  message("  -> Generating unintegrated UMAP (Track A)...")
  data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                        graph.name = "pca_nn", verbose = FALSE, k.param = UMAP_N_NEIGHBORS)
  data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, algorithm = "leiden",
                       graph.name = "pca_nn", cluster.name = "clusters_none", verbose = FALSE)
  data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                  n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                  reduction.name = "umap_none", verbose = FALSE, n.epochs = 500)
  gc()

  # --- Track B: Harmony correction of the PCA ---
  # Current harmony (CRAN and GitHub) uses reduction.use / reduction.save /
  # dims.use. The legacy 'reduction=' argument is gone.
  message("  -> Running Harmony via harmony::RunHarmony (Track B)...")
  data <- RunHarmony(data,
                     group.by.vars  = "SampleID",
                     reduction.use  = "pca",
                     dims.use       = 1:N_PCS_TO_USE,
                     reduction.save = "harmony",
                     verbose        = TRUE)
  gc()

} else {
  # ---------------------------------------------------------------------------
  # ROUTE B: Seurat v5 layered workflow (split -> per-layer -> IntegrateLayers).
  # ---------------------------------------------------------------------------
  data[["RNA"]] <- split(data[["RNA"]], f = data$SampleID)
  message("  -> RNA assay split into layers by SampleID.")
  gc()

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

  data <- NormalizeData(data, verbose = TRUE)
  data <- FindVariableFeatures(data, nfeatures = N_VARIABLE_FEATURES, verbose = TRUE)
  data <- ScaleData(data, verbose = TRUE)
  data <- RunPCA(data, npcs = N_PCS_TO_USE, verbose = TRUE)
  message("  -> Per-layer normalization, scaling, and PCA complete.")
  gc()

  # --- Track A: unintegrated PCA (diagnostic) ---
  message("  -> Generating unintegrated UMAP (Track A)...")
  data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                        graph.name = "pca_nn", verbose = FALSE, k.param = UMAP_N_NEIGHBORS)
  data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, algorithm = "leiden",
                       graph.name = "pca_nn", cluster.name = "clusters_none", verbose = FALSE)
  data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "pca",
                  n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                  reduction.name = "umap_none", verbose = FALSE, n.epochs = 500)
  gc()

  # --- Track B: Harmony via IntegrateLayers ---
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

  # --- Re-join layers for downstream analysis ---
  message("  -> Re-joining layers into a single matrix for downstream analysis...")
  data <- JoinLayers(data)
  if (USE_BPCELLS) {
    data[["RNA"]]$counts     <- as(LayerData(data, assay = "RNA", layer = "counts"),     "dgCMatrix")
    gc()
    data[["RNA"]]$data       <- as(LayerData(data, assay = "RNA", layer = "data"),       "dgCMatrix")
    gc()
    data[["RNA"]]$scale.data <- as(LayerData(data, assay = "RNA", layer = "scale.data"), "dgCMatrix")
    gc()
  }
  message("  -> Layers successfully joined.")
}

# --- Track B graph, clustering and UMAP (shared by both routes) -------------
# RunHarmony and IntegrateLayers both leave a 'harmony' reduction; the
# neighbour graph, clusters and UMAP built on it are identical downstream.
data <- FindNeighbors(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                      graph.name = "harmony_nn", verbose = TRUE, k.param = UMAP_N_NEIGHBORS)
data <- FindClusters(data, resolution = CLUSTER_RESOLUTION, algorithm = "leiden",
                     graph.name = "harmony_nn", cluster.name = "clusters_harmony", verbose = TRUE)
data <- RunUMAP(data, dims = 1:N_PCS_TO_USE, reduction = "harmony",
                n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
                reduction.name = "umap_harmony", verbose = TRUE, n.epochs = 500)
gc()

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

