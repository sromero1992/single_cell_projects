# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 9: CELL SCORES
#   Developmental potency, stemness, and transcriptional entropy
# Version: 1.0 (UNIFIED build)
#
# UNIFIED BUILD: part of unified_pipeline/. Consumes the annotated object from
#   06_annotation_unifier.R. Adds per-cell potency/entropy scores and writes an
#   enriched .rds that scripts 07/08/10 can also consume.
#
# PURPOSE:
#   Quantify "how differentiated is this cell?" from three independent angles,
#   so that conclusions do not rest on any single method:
#
#     1. CytoTRACE 2  — supervised deep-learning prediction of ABSOLUTE
#        developmental potential. Returns a calibrated potency score in [0,1]
#        plus a discrete potency category (Differentiated ... Totipotent).
#        Calibrated across datasets, so scores are comparable between runs.
#        Kang et al., Nature Methods 2025. doi:10.1038/s41592-025-02857-2
#
#     2. CytoTRACE (v1) — the original unsupervised method. Its core insight is
#        that the NUMBER OF DETECTABLY EXPRESSED GENES per cell tracks
#        differentiation state. Returns a RELATIVE ordering within the dataset
#        (0 = most differentiated, 1 = least), not an absolute scale.
#        Gulati et al., Science 2020. doi:10.1126/science.aax0249
#
#     3. Transcriptional entropy — computed natively in this script, no
#        external package required, so it ALWAYS runs even if 1 and 2 fail to
#        install. Rooted in the observation that stem/progenitor cells have
#        more uniform (higher-entropy) transcriptomes, while differentiated
#        cells concentrate expression into a focused program.
#        Teschendorff & Enver, Nat Commun 2017. doi:10.1038/ncomms15599
#
#   Methods 1 and 2 are orthogonal in construction (supervised vs unsupervised)
#   and method 3 is model-free. Where all three agree, the ordering is solid.
#   Where they disagree, that disagreement is itself the finding — this script
#   quantifies it explicitly (Spearman correlation matrix + concordance plots).
#
# ============================================================================
# CRITICAL INPUT REQUIREMENT — READ BEFORE RUNNING
# ============================================================================
#   CytoTRACE 2 requires RAW COUNTS or CPM/TPM. It must NOT receive
#   log-transformed data. This script therefore pulls from the "counts" layer,
#   never "data". If your counts layer holds DecontX-corrected values, that is
#   fine and expected (Script 01 writes corrected counts back into "counts"),
#   but confirm they are non-negative and not log-scaled.
#
#   A minimum of ~500-1000 detected genes per cell is recommended for reliable
#   CytoTRACE 2 predictions. The pipeline's POST_MIN_GENES default of 500 sits
#   right at that boundary; this script reports the fraction of cells below
#   1000 genes so you can judge how much weight to put on the result.
#
# ============================================================================
# BATCHING STRATEGY
# ============================================================================
#   CytoTRACE is run PER SAMPLE by default (CT2_RUN_PER_SAMPLE = TRUE): only one
#   sample's cells are materialised at a time, so peak memory is bounded by the
#   biggest single sample rather than the whole dataset, and each sample is scored
#   intact (better than CytoTRACE 2's internal random subsampling of the mixed
#   object). The absolute outputs (CytoTRACE2_Score, _Potency) are calibrated and
#   stay comparable across samples; CytoTRACE2_Relative is rescaled within each
#   sample. Set CT2_RUN_PER_SAMPLE = FALSE to score the whole object in one pass.
#
# INSTALLATION:
#   CytoTRACE 2 and CytoTRACE v1 are NOT on CRAN/Bioconductor and are the two
#   most failure-prone installs in this pipeline. See INSTALL_NOTES.md in this
#   folder for step-by-step instructions and fixes for known errors.
#   This script DEGRADES GRACEFULLY: any method whose package is missing is
#   skipped with a clear message, and the remaining methods still run.
#
# OUTPUT:
#   - <PROJECT>_with_cell_scores.rds        (object + all score columns)
#   - cell_scores/ plots (UMAPs, boxplots, correlation heatmap, ridge plots)
#   - cell_scores/cell_scores_summary.xlsx  (per cell type / per group stats)
#   - cell_scores/cell_scores_per_cell.csv.gz (full per-cell table)
#
# NEXT STEP:
#   10_trajectory_cellrank.R — uses these scores to direct trajectory inference.
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)   # stat_compare_means() for the significance brackets/stars
library(patchwork)
library(writexl)
library(Matrix)

set.seed(123)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) --------------------------
# =============================================================================

# --- 1.1: Project Identity & Paths -------------------------------------------
PROJECT_NAME <- "Nr4a1_s17_ack"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# Input object — the final annotated object from Script 06.
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_unified_annotated.rds"))

# Or point at a subtype object to score within one lineage:
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))

SCORES_DIR <- file.path(OUTPUT_DIR, "cell_scores")

# Per-sample CytoTRACE results are checkpointed here as CSV. A crash at sample
# 10/16 keeps the first 9 - re-running reloads them and scores only what's left.
CELL_POTENCY_SCRATCH <- file.path(OUTPUT_DIR, "cell_potency_scratch")

# --- 1.2: Metadata Columns ---------------------------------------------------
MODE <- "subtypes"
CELLTYPE_COLUMN  <- if (MODE == "broad") "CellType_broad" else "CellType"
SAMPLE_COLUMN    <- "SampleID"
CONDITION_COLUMN <- "Genotype_sex"

# Order for plot axes. Put the control/reference level FIRST.
# Set to NULL to use whatever order the factor already has.
CONDITION_LEVELS <- c(
  "WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female",
  "WT_Male",   "Polyp_Male",   "Polyp_NR4a1_KO_Male"
)

# Pairwise contrasts: each is c(group1, group2). Used for the per-contrast score
# comparisons (Wilcoxon per cell type) on both the potency/entropy and the
# AUCell pathway scores.
CONTRASTS_LIST <- list(
  Nr4a1_KO_polyp_vs_polyp_Female = c("Polyp_NR4a1_KO_Female", "Polyp_Female"),
  Nr4a1_KO_polyp_vs_polyp_Male   = c("Polyp_NR4a1_KO_Male",   "Polyp_Male"),
  Polyp_vs_WT_Female             = c("Polyp_Female",          "WT_Female"),
  Polyp_vs_WT_Male               = c("Polyp_Male",            "WT_Male")
)

# Reduction used for score UMAPs. Script 01 produces "umap_harmony".
UMAP_REDUCTION <- "umap_harmony"

# --- 1.3: Which Methods To Run -----------------------------------------------
RUN_CYTOTRACE2 <- TRUE   # absolute potency (recommended primary method)
RUN_CYTOTRACE1 <- TRUE   # original relative CytoTRACE; skipped if not installed
RUN_ENTROPY    <- TRUE   # native entropy metrics; no dependencies, always works
RUN_PATHWAY_SCORES <- TRUE  # AUCell per-cell scoring of the GOBP gene sets below
RUN_CCAT  <- TRUE    # CCAT connectome-correlation potency (SCENT PPI). Fully guarded:
                     #   any failure in the fragile mouse->human->PPI chain is skipped.
RUN_SCENT <- FALSE   # SCENT signalling entropy (CompSR) - slow (min); fully guarded.

# RESUME_SCORES: if TRUE and a previous run already wrote
# cell_scores/potency_scores_per_cell.csv, load those columns onto the object and
# SKIP all (expensive) score computation - just redo the plots/tables/save. Use
# this after an OOM in plotting, or to re-plot without re-scoring.
RESUME_SCORES <- FALSE
CCAT_SPECIES <- "mouse"  # "mouse" -> homologene 10090->9606 to the human PPI; "human" = direct

# --- 1.3b: Pathway (AUCell) scoring ------------------------------------------
# Per-cell activity of curated GO Biological Process gene sets, scored with
# AUCell. Gene sets are pulled from org.Mm.eg.db by GO ID (GOALL = term +
# descendants); if that is unavailable, the curated fallback list is used.
# Column names become AUCell_<name>; these join the potency/entropy scores in
# every by-condition plot and in the per-contrast comparisons.
PATHWAY_SPECIES <- "mouse"    # "mouse" (org.Mm.eg.db) or "human" (org.Hs.eg.db)
PATHWAY_GO_TERMS <- list(
  canonical_WNT         = "GO:0060070",  # canonical Wnt signaling pathway
  noncanonical_WNT      = "GO:0035567",  # non-canonical Wnt signaling pathway
  apoptosis             = "GO:0006915",  # apoptotic process
  programmed_cell_death = "GO:0012501",  # programmed cell death
  ferroptosis           = "GO:0097707"   # ferroptosis
)
# Curated fallback gene sets (mouse symbols) used only if the GO lookup fails.
PATHWAY_FALLBACK <- list(
  canonical_WNT = c("Ctnnb1","Wnt3","Wnt3a","Lef1","Tcf7","Axin2","Ccnd1","Myc",
                    "Apc","Gsk3b","Dvl1","Fzd1","Lrp6","Sox9","Ascl2"),
  noncanonical_WNT = c("Wnt5a","Wnt11","Ror2","Vangl2","Prickle1","Daam1","Rhoa",
                       "Rac1","Jun","Nfatc1","Ryk","Camk2a","Ptk7"),
  apoptosis = c("Casp3","Casp8","Casp9","Bax","Bak1","Bcl2","Bcl2l1","Bid","Apaf1",
                "Trp53","Fas","Fadd","Cycs","Diablo"),
  programmed_cell_death = c("Casp3","Casp1","Ripk1","Ripk3","Mlkl","Tnf","Fas","Bax",
                            "Bcl2","Gsdmd","Casp8","Tnfrsf1a","Zbp1"),
  ferroptosis = c("Gpx4","Slc7a11","Acsl4","Ncoa4","Fth1","Ftl1","Slc3a2","Nfe2l2",
                  "Tfrc","Alox15","Sat1","Aifm2","Vdac2")
)
PATHWAY_MIN_GENES <- 5   # skip a gene set with fewer than this many genes present

# --- 1.4: CytoTRACE 2 Parameters ---------------------------------------------
# SPECIES: "mouse" or "human". CytoTRACE 2's feature set is mouse-based; for
#   human input it performs orthology mapping internally.
CT2_SPECIES <- "mouse"

# CT2_RUN_PER_SAMPLE: score each SampleID separately. Recommended for large data:
#   each run materialises only ONE sample's cells (memory bounded by the biggest
#   sample, not the whole dataset), keeps samples biologically intact, and is
#   checkpointed to CSV so a crash mid-loop is resumable. The absolute
#   CytoTRACE2_Score/_Potency stay calibrated across samples; only the RELATIVE
#   score is rescaled within each sample. Set FALSE to run the whole object once.
CT2_RUN_PER_SAMPLE <- TRUE

# CT2_SLOT: which layer to feed. MUST be raw/CPM counts, never log-transformed.
CT2_SLOT <- "counts"

# CT2_PARALLELIZE: run the model + smoothing on multiple threads. Default FALSE
#   (serial) = LOWEST memory, which matters on big objects: TRUE spawns workers
#   that each hold a copy of the expression block and can OOM a large run. Turn it
#   TRUE only if you have RAM headroom and want speed.
CT2_PARALLELIZE <- FALSE

# CT2_NCORES: cores when parallelizing. NULL = the package default (auto-detects
#   and uses half the cores; Windows forced to 1). The authors advise 1-2 on
#   machines with < 16 GB RAM, since each worker holds a copy of the block.
CT2_NCORES <- NULL

# Batch sizes (package defaults). To reproduce the manuscript on a big machine,
# the authors use batch_size = 100000, smooth_batch_size = 10000.
#   batch_size        — cells processed at once (recommended for >10K cells)
#   smooth_batch_size — subsample within a batch for the diffusion smoothing
CT2_BATCH_SIZE        <- 10000
CT2_SMOOTH_BATCH_SIZE <- 1000

# CT2_USE_PREKNN: use the pre-KNN-smoothing score as the primary value.
#   The final CytoTRACE2_Score is smoothed over neighbours, which stabilises
#   it but can drag rare populations (<= ~5 cells) toward their more abundant
#   neighbours. Set TRUE if rare-population potency is the point of the study.
CT2_USE_PREKNN <- FALSE

# CT2_MIN_CELLS: skip a sample with fewer cells than this (predictions on very
#   small inputs are unstable, and the smoothing step needs neighbours).
CT2_MIN_CELLS <- 50

# --- 1.5: Entropy Parameters -------------------------------------------------
# Entropy is computed on the counts layer, per cell, over detected genes only.
#
# ENTROPY_NORMALIZE: divide Shannon entropy by log(n_detected_genes) to give a
#   value in [0,1] that is independent of library complexity. Strongly
#   recommended: raw Shannon entropy is heavily confounded by sequencing depth,
#   and without this you will mostly be measuring depth, not biology.
ENTROPY_NORMALIZE <- TRUE

# ENTROPY_MIN_GENES: cells below this many detected genes get NA rather than a
#   misleadingly precise entropy value.
ENTROPY_MIN_GENES <- 200

# --- 1.6: Comparison & Statistics --------------------------------------------
RUN_GROUP_STATS <- TRUE   # Wilcoxon / Kruskal-Wallis across CONDITION_COLUMN
STATS_MIN_CELLS <- 20     # skip a cell type x group cell with fewer cells
STATS_PADJ_METHOD <- "BH" # multiple testing correction across cell types

# --- 1.7: Plotting -----------------------------------------------------------
DPI_SETTING  <- 300
PLOT_WIDTH   <- 9
PLOT_HEIGHT  <- 6
POINT_SIZE   <- 0.3

# Colour scale for potency UMAPs (low potency -> high potency)
POTENCY_COLORS <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")

# Canonical CytoTRACE 2 category order, low to high potency.
POTENCY_LEVELS <- c("Differentiated", "Unipotent", "Oligopotent",
                    "Multipotent", "Pluripotent", "Totipotent")

# =============================================================================
# --- PART 2: EXECUTION (DO NOT EDIT BELOW THIS LINE) -------------------------
# =============================================================================

if (!dir.exists(SCORES_DIR)) dir.create(SCORES_DIR, recursive = TRUE)

# --- Availability checks -----------------------------------------------------
# Detect once, up front, so the run either proceeds knowingly or stops early
# rather than failing after an hour of computation.
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

CT2_AVAILABLE <- has_pkg("CytoTRACE2")
CT1_AVAILABLE <- has_pkg("CytoTRACE")

message("=== Script 09: Cell Scores ===")
message(paste0("  CytoTRACE2 package : ", if (CT2_AVAILABLE) "FOUND" else "NOT INSTALLED"))
message(paste0("  CytoTRACE  package : ", if (CT1_AVAILABLE) "FOUND" else "NOT INSTALLED"))

if (RUN_CYTOTRACE2 && !CT2_AVAILABLE) {
  message("  [SKIP] RUN_CYTOTRACE2 is TRUE but the package is missing.")
  message("         See INSTALL_NOTES.md section 1. Continuing without it.")
  RUN_CYTOTRACE2 <- FALSE
}
if (RUN_CYTOTRACE1 && !CT1_AVAILABLE) {
  message("  [SKIP] RUN_CYTOTRACE1 is TRUE but the package is missing.")
  message("         See INSTALL_NOTES.md section 2. Continuing without it.")
  RUN_CYTOTRACE1 <- FALSE
}
if (!RUN_CYTOTRACE2 && !RUN_CYTOTRACE1 && !RUN_ENTROPY) {
  stop("All scoring methods are disabled or unavailable - nothing to do.", call. = FALSE)
}

# --- Load object -------------------------------------------------------------
message(paste("  Loading:", RDS_PATH))
if (!file.exists(RDS_PATH)) {
  stop(paste0("Input .rds not found:\n  ", RDS_PATH,
              "\nRun Script 06 first, or correct RDS_PATH."), call. = FALSE)
}
data <- readRDS(RDS_PATH)
DefaultAssay(data) <- "RNA"

# Seurat v5 keeps per-sample layers separate; scoring needs one matrix.
if (inherits(data[["RNA"]], "Assay5")) {
  lyrs <- Layers(data[["RNA"]], search = "counts")
  if (length(lyrs) > 1) {
    message(paste0("  Joining ", length(lyrs), " count layers into one matrix..."))
    data <- JoinLayers(data)
  }
}

message(paste0("  Object: ", ncol(data), " cells x ", nrow(data), " genes"))

# --- Validate configured columns exist --------------------------------------
missing_cols <- setdiff(
  c(CELLTYPE_COLUMN, SAMPLE_COLUMN, CONDITION_COLUMN),
  colnames(data@meta.data)
)
if (length(missing_cols) > 0) {
  stop(paste0("Metadata column(s) not found in the object: ",
              paste(missing_cols, collapse = ", "),
              "\nAvailable columns:\n  ",
              paste(colnames(data@meta.data), collapse = ", ")), call. = FALSE)
}

if (!is.null(CONDITION_LEVELS)) {
  present <- unique(as.character(data@meta.data[[CONDITION_COLUMN]]))
  unknown <- setdiff(CONDITION_LEVELS, present)
  if (length(unknown) > 0) {
    warning(paste0("CONDITION_LEVELS contains values absent from the data: ",
                   paste(unknown, collapse = ", ")))
  }
  data@meta.data[[CONDITION_COLUMN]] <- factor(
    data@meta.data[[CONDITION_COLUMN]],
    levels = intersect(CONDITION_LEVELS, present)
  )
}

# --- Input quality report ----------------------------------------------------
# CytoTRACE 2 is unreliable below ~500-1000 detected genes per cell. Report the
# exposure so the result can be interpreted honestly rather than assumed clean.
n_genes_per_cell <- data$nFeature_RNA
frac_low <- mean(n_genes_per_cell < 1000, na.rm = TRUE)
message(paste0("  Detected genes per cell: median = ",
               round(stats::median(n_genes_per_cell, na.rm = TRUE)),
               " | ", round(frac_low * 100, 1), "% of cells below 1000"))
if (frac_low > 0.5) {
  warning(paste0("Over half of cells have <1000 detected genes (",
                 round(frac_low * 100, 1), "%). CytoTRACE 2 predictions on ",
                 "such cells are less reliable - interpret with caution and ",
                 "lean on the entropy metrics for cross-checking."))
}

# =============================================================================
# --- RESUME: reload saved scores and skip computation ------------------------
# =============================================================================
# Prevents recomputing potency after an OOM/crash in the plotting below: load the
# saved per-cell scores, attach them by barcode, and turn every compute step off.
pathway_cols <- character(0)
ccat_cols    <- character(0)
.scores_csv  <- file.path(SCORES_DIR, "potency_scores_per_cell.csv")
if (RESUME_SCORES && file.exists(.scores_csv)) {
  message("\n=== RESUME: loading saved scores, skipping computation ===")
  sc <- utils::read.csv(.scores_csv, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(sc) <- sc$Barcode
  sc <- sc[match(colnames(data), rownames(sc)), , drop = FALSE]
  add_cols <- setdiff(colnames(sc),
                      c("Barcode", SAMPLE_COLUMN, CELLTYPE_COLUMN, CONDITION_COLUMN))
  for (cn in add_cols) data@meta.data[[cn]] <- sc[[cn]]
  # Restore the potency-category ordering lost through CSV round-trip.
  for (pc in intersect(c("CytoTRACE2_Potency", "preKNN_CytoTRACE2_Potency"), add_cols))
    data@meta.data[[pc]] <- factor(as.character(data@meta.data[[pc]]), levels = POTENCY_LEVELS)
  pathway_cols <- grep("^AUCell_", add_cols, value = TRUE)
  ccat_cols    <- intersect(c("CCAT_score", "SCENT_score"), add_cols)
  RUN_CYTOTRACE2 <- RUN_CYTOTRACE1 <- RUN_ENTROPY <- FALSE
  RUN_PATHWAY_SCORES <- RUN_CCAT <- RUN_SCENT <- FALSE
  message(paste0("  Loaded ", length(add_cols), " score columns; computation skipped."))
}

# =============================================================================
# --- STEP 1: TRANSCRIPTIONAL ENTROPY (native, no dependencies) ---------------
# =============================================================================
# Rationale (Teschendorff & Enver 2017): a stem or progenitor cell keeps many
# transcriptional programs simultaneously accessible, so its expression is
# spread comparatively evenly across genes - high entropy. As a cell commits,
# expression concentrates into a focused lineage program - low entropy.
#
# For each cell we treat the count vector as a probability distribution over
# detected genes, p_i = count_i / sum(counts), and compute Shannon entropy:
#
#     H = - sum_i p_i * log(p_i)
#
# CONFOUND AND CORRECTION:
#   Raw H scales with the number of detected genes, which scales with
#   sequencing depth. Comparing raw H between cells largely compares depth.
#   Normalising by the theoretical maximum, log(n_detected), gives
#
#     H_norm = H / log(n_detected)     in [0, 1]
#
#   which measures how EVENLY expression is distributed independent of how many
#   genes were captured. H_norm is the value to compare across cells; keep
#   ENTROPY_NORMALIZE = TRUE unless you have a specific reason not to.
# =============================================================================
if (RUN_ENTROPY) {
  message("\n=== STEP 1: Transcriptional entropy ===")

  counts_mat <- GetAssayData(data, assay = "RNA", layer = "counts")

  # Work in sparse form and column-by-column blocks to bound peak memory.
  # For a dgCMatrix the non-zero entries of column j are x[(p[j]+1):p[j+1]],
  # so entropy per cell can be computed without ever densifying the matrix.
  message("  Computing per-cell Shannon entropy on the counts layer...")

  compute_entropy_sparse <- function(m) {
    m <- as(m, "CsparseMatrix")
    n_cells_local <- ncol(m)
    H        <- numeric(n_cells_local)
    n_det    <- integer(n_cells_local)
    tot      <- numeric(n_cells_local)

    p_idx <- m@p
    xvals <- m@x

    for (j in seq_len(n_cells_local)) {
      lo <- p_idx[j] + 1L
      hi <- p_idx[j + 1L]
      if (hi < lo) {                 # completely empty cell
        H[j] <- NA_real_; n_det[j] <- 0L; tot[j] <- 0
        next
      }
      v <- xvals[lo:hi]
      v <- v[v > 0]                  # guard against explicit zeros
      s <- sum(v)
      if (s <= 0 || length(v) == 0) {
        H[j] <- NA_real_; n_det[j] <- 0L; tot[j] <- 0
        next
      }
      p <- v / s
      H[j]     <- -sum(p * log(p))
      n_det[j] <- length(v)
      tot[j]   <- s
    }
    list(H = H, n_det = n_det, total = tot)
  }

  ent <- compute_entropy_sparse(counts_mat)

  # Cells with too few genes get NA rather than a falsely precise number.
  too_few <- ent$n_det < ENTROPY_MIN_GENES
  ent$H[too_few] <- NA_real_

  data$entropy_shannon    <- ent$H
  data$entropy_n_detected <- ent$n_det

  # Normalised entropy (evenness): H / log(n_detected), in [0,1].
  # log(1) = 0 would divide by zero, so single-gene cells are excluded.
  denom <- log(ent$n_det)
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  data$entropy_normalized <- ent$H / denom

  # Primary entropy column used for all downstream comparisons.
  data$entropy_score <- if (ENTROPY_NORMALIZE) {
    data$entropy_normalized
  } else {
    data$entropy_shannon
  }

  # Gene counts: the raw determinant underlying original CytoTRACE. Kept
  # explicitly so the "more genes = less differentiated" signal can be
  # inspected directly, independent of any package.
  data$gene_counts_score <- ent$n_det

  n_na <- sum(is.na(data$entropy_score))
  message(paste0("  Entropy computed. ", ncol(data), " cells | ",
                 n_na, " NA (below ENTROPY_MIN_GENES = ", ENTROPY_MIN_GENES, ")"))
  message(paste0("  ", if (ENTROPY_NORMALIZE) "Normalized" else "Raw",
                 " entropy range: ",
                 paste(round(range(data$entropy_score, na.rm = TRUE), 4),
                       collapse = " - ")))

  rm(counts_mat, ent, denom, too_few); gc()
}

# =============================================================================
# --- STEP 2: CytoTRACE 2 (absolute developmental potential) ------------------
# =============================================================================
if (RUN_CYTOTRACE2) {
  message("\n=== STEP 2: CytoTRACE 2 ===")
  suppressPackageStartupMessages(library(CytoTRACE2))

  # Container for results (whole-dataset pass).
  ct2_cols <- c("CytoTRACE2_Score", "CytoTRACE2_Potency", "CytoTRACE2_Relative",
                "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency")
  ct2_all  <- NULL

  # ---- Helper: build the SPARSE counts CytoTRACE 2 will subsample ------------
  # We hand CytoTRACE 2 a lean Seurat object (is_seurat = TRUE) and let IT densify
  # one ~batch_size chunk at a time. So this returns a SPARSE matrix - it must
  # NEVER densify the whole thing (that ~25 GB copy was the OOM). Two things still
  # matter for correct scoring:
  #   1. Gene SYMBOL rownames. If the object carries Ensembl IDs they are mapped
  #      to symbols (mouse: org.Mm.eg.db, human: org.Hs.eg.db); otherwise
  #      CytoTRACE 2's internal gene panel matches almost nothing.
  #   2. Deduplicated rownames (keep the highest mean-expression row per symbol).
  prep_ct2_input <- function(obj) {
    m <- GetAssayData(obj, assay = "RNA", layer = CT2_SLOT)   # RAW counts, never log

    # Ensembl -> symbol, but only if the rownames actually look like Ensembl IDs.
    if (length(rownames(m)) && all(grepl("^ENS", rownames(m)))) {
      org_db <- if (CT2_SPECIES == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
      if (has_pkg(org_db) && has_pkg("AnnotationDbi")) {
        message("    Mapping Ensembl IDs -> gene symbols via ", org_db, " ...")
        sym <- suppressMessages(AnnotationDbi::mapIds(
          getExportedValue(org_db, org_db),
          keys = rownames(m), keytype = "ENSEMBL",
          column = "SYMBOL", multiVals = "first"))
        keep <- !is.na(sym) & sym != ""
        m <- m[keep, , drop = FALSE]
        rownames(m) <- sym[keep]
      } else {
        message("    [NOTE] Rownames look like Ensembl IDs but ", org_db,
                " is not installed; CytoTRACE 2 may score poorly.")
      }
    }

    # Deduplicate symbols: keep the highest mean-expression row for each.
    if (any(duplicated(rownames(m)))) {
      o <- order(Matrix::rowMeans(m), decreasing = TRUE)
      m <- m[o, , drop = FALSE]
      m <- m[!duplicated(rownames(m)), , drop = FALSE]
    }

    m   # SPARSE counts (symbol rownames, deduplicated); CT2 densifies per chunk
  }

  # ---- Helper: run CytoTRACE 2 (memory-lean) --------------------------------
  # Builds a COUNTS-ONLY Seurat object (no data/scale.data/reductions/other
  # assays) and passes it with is_seurat = TRUE, so CytoTRACE 2 pulls the sparse
  # matrix and densifies only one batch_size chunk at a time - never the whole
  # 200k x 16k matrix. parallelize_* follow CT2_PARALLELIZE (default FALSE =
  # serial = lowest RAM; TRUE spawns workers that each copy the block).
  # Do NOT add verbose = TRUE (invalid argument).
  run_ct2_block <- function(obj, label) {
    if (ncol(obj) < CT2_MIN_CELLS) {
      message(paste0("    [SKIP] ", label, ": only ", ncol(obj),
                     " cells (< CT2_MIN_CELLS = ", CT2_MIN_CELLS, ")"))
      return(NULL)
    }
    message(paste0("    -> ", label, ": ", ncol(obj), " cells..."))
    ct2_par <- isTRUE(CT2_PARALLELIZE) && .Platform$OS.type != "windows"

    res <- tryCatch({
      cnts <- prep_ct2_input(obj)                              # sparse, symbol rownames
      lean <- Seurat::CreateSeuratObject(counts = cnts)        # counts-only, minimal RAM
      rm(cnts); gc()
      out <- cytotrace2(
        lean,
        is_seurat             = TRUE,
        slot_type             = "counts",
        species               = CT2_SPECIES,
        batch_size            = CT2_BATCH_SIZE,
        smooth_batch_size     = CT2_SMOOTH_BATCH_SIZE,
        parallelize_models    = ct2_par,
        parallelize_smoothing = ct2_par,
        ncores                = CT2_NCORES,
        seed                  = 14
      )
      # is_seurat = TRUE returns a Seurat object; scores live in its metadata.
      md <- out@meta.data
      rm(out, lean); gc()
      keep <- intersect(ct2_cols, colnames(md))
      if (length(keep) == 0) stop("cytotrace2() returned no recognised prediction columns.")
      res_df <- md[, keep, drop = FALSE]                       # rownames = cell barcodes
      for (sc_col in intersect(c("CytoTRACE2_Score", "preKNN_CytoTRACE2_Score"),
                               colnames(res_df))) {
        res_df[[sc_col]] <- pmin(pmax(res_df[[sc_col]], 0), 1)
      }
      res_df
    }, error = function(e) {
      message(paste0("    [WARNING] CytoTRACE 2 failed on ", label, ": ", e$message))
      NULL
    })
    res
  }

  if (CT2_RUN_PER_SAMPLE) {
    # Per sample: materialise ONE sample at a time (memory bounded by the biggest
    # sample), checkpoint each to CSV, and free it before the next.
    samples <- unique(as.character(data@meta.data[[SAMPLE_COLUMN]]))
    message(paste0("  Running CytoTRACE 2 per sample (", length(samples), " samples)..."))
    if (!dir.exists(CELL_POTENCY_SCRATCH)) dir.create(CELL_POTENCY_SCRATCH, recursive = TRUE)

    res_list <- list()
    for (s in samples) {
      ckpt <- file.path(CELL_POTENCY_SCRATCH,
                        paste0("ct2_", gsub("[^A-Za-z0-9_.-]", "_", s), ".csv"))
      if (file.exists(ckpt)) {                                  # resume
        message(paste0("    [checkpoint] ", s, ": loading ", basename(ckpt)))
        res_list[[s]] <- utils::read.csv(ckpt, row.names = 1, check.names = FALSE,
                                         stringsAsFactors = FALSE)
        next
      }
      cells_s <- colnames(data)[as.character(data@meta.data[[SAMPLE_COLUMN]]) == s]
      obj_s   <- subset(data, cells = cells_s)
      r <- run_ct2_block(obj_s, s)
      if (!is.null(r)) {
        utils::write.csv(r, ckpt, row.names = TRUE)            # save immediately
        res_list[[s]] <- r
        message(paste0("    [saved] ", s, " -> ", basename(ckpt)))
      }
      rm(obj_s, r); gc()                                        # free before next sample
    }
    if (length(res_list) > 0) {
      common   <- Reduce(intersect, lapply(res_list, colnames))
      res_list <- lapply(res_list, function(d) d[, common, drop = FALSE])
      ct2_all  <- do.call(rbind, res_list)
    }
    rm(res_list); gc()
  } else {
    message("  Running CytoTRACE 2 on the full object (one pass)...")
    ct2_all <- run_ct2_block(data, "ALL")
    gc()
  }

  # ---- Attach results ------------------------------------------------------
  if (!is.null(ct2_all) && nrow(ct2_all) > 0) {
    # Align strictly by barcode; never assume row order matches.
    ct2_aligned <- ct2_all[match(colnames(data), rownames(ct2_all)), , drop = FALSE]
    rownames(ct2_aligned) <- colnames(data)

    for (cn in colnames(ct2_aligned)) {
      data@meta.data[[cn]] <- ct2_aligned[[cn]]
    }

    # Enforce the canonical low-to-high potency ordering on the category
    # columns, so every boxplot and table reads in biological order.
    for (pc in intersect(c("CytoTRACE2_Potency", "preKNN_CytoTRACE2_Potency"),
                         colnames(data@meta.data))) {
      data@meta.data[[pc]] <- factor(as.character(data@meta.data[[pc]]),
                                     levels = POTENCY_LEVELS)
    }

    # Primary potency column: smoothed by default, pre-KNN if the study is
    # about rare populations (see CT2_USE_PREKNN).
    if (CT2_USE_PREKNN && "preKNN_CytoTRACE2_Score" %in% colnames(data@meta.data)) {
      data$potency_score <- data$preKNN_CytoTRACE2_Score
      message("  Primary potency = preKNN_CytoTRACE2_Score (unsmoothed).")
    } else if ("CytoTRACE2_Score" %in% colnames(data@meta.data)) {
      data$potency_score <- data$CytoTRACE2_Score
      message("  Primary potency = CytoTRACE2_Score (KNN-smoothed).")
    }

    n_scored <- sum(!is.na(data$potency_score))
    message(paste0("  CytoTRACE 2 complete: ", n_scored, " / ", ncol(data),
                   " cells scored (",
                   round(n_scored / ncol(data) * 100, 1), "%)"))
    if (n_scored < ncol(data)) {
      message("  [NOTE] Unscored cells belong to samples that were skipped or failed.")
    }
  } else {
    message("  [WARNING] CytoTRACE 2 produced no results. Continuing without it.")
    RUN_CYTOTRACE2 <- FALSE
  }
  rm(ct2_all); gc()
}

# =============================================================================
# --- STEP 3: CytoTRACE v1 (relative differentiation order) -------------------
# =============================================================================
# The original method. Reported value is RELATIVE within the input dataset:
# 0 = most differentiated, 1 = least differentiated. It is not calibrated
# across datasets, so it is used here as a cross-check on the ordering, not as
# an absolute claim.
#
# NOTE: CytoTRACE v1 expects a plain counts matrix (genes x cells), not a
# Seurat object, and can be memory-hungry. It is run once over the whole dataset
# so the relative ordering is computed across all cells together.
# =============================================================================
if (RUN_CYTOTRACE1) {
  message("\n=== STEP 3: CytoTRACE (v1) ===")
  suppressPackageStartupMessages(library(CytoTRACE))

  run_ct1_block <- function(mat, label) {
    if (ncol(mat) < CT2_MIN_CELLS) {
      message(paste0("    [SKIP] ", label, ": only ", ncol(mat), " cells"))
      return(NULL)
    }
    message(paste0("    -> ", label, ": ", ncol(mat), " cells..."))
    tryCatch({
      # CytoTRACE() wants a dense matrix of counts with gene rownames.
      dense <- as.matrix(mat)
      res   <- CytoTRACE(dense, ncores = 1)
      out   <- data.frame(
        CytoTRACE1_Score   = as.numeric(res$CytoTRACE),
        CytoTRACE1_Rank    = as.numeric(res$CytoTRACErank),
        CytoTRACE1_GCS     = as.numeric(res$GCS),
        row.names          = names(res$CytoTRACE),
        stringsAsFactors   = FALSE
      )
      rm(dense); gc()
      out
    }, error = function(e) {
      message(paste0("    [WARNING] CytoTRACE v1 failed on ", label, ": ", e$message))
      NULL
    })
  }

  counts_all <- GetAssayData(data, assay = "RNA", layer = "counts")
  ct1_list   <- list()

  if (CT2_RUN_PER_SAMPLE) {
    # v1 densifies the ENTIRE matrix (it cannot chunk), so per sample is the only
    # memory-safe way on a large object - one sample's dense block at a time.
    samples <- unique(as.character(data@meta.data[[SAMPLE_COLUMN]]))
    for (s in samples) {
      idx <- which(as.character(data@meta.data[[SAMPLE_COLUMN]]) == s)
      r <- run_ct1_block(counts_all[, idx, drop = FALSE], s)
      if (!is.null(r)) ct1_list[[s]] <- r
      gc()
    }
  } else {
    r <- run_ct1_block(counts_all, "ALL")
    if (!is.null(r)) ct1_list[["ALL"]] <- r
  }

  if (length(ct1_list) > 0) {
    ct1_all <- do.call(rbind, ct1_list)
    ct1_all <- ct1_all[match(colnames(data), rownames(ct1_all)), , drop = FALSE]
    data$CytoTRACE1_Score <- ct1_all$CytoTRACE1_Score
    data$CytoTRACE1_Rank  <- ct1_all$CytoTRACE1_Rank
    data$CytoTRACE1_GCS   <- ct1_all$CytoTRACE1_GCS
    message(paste0("  CytoTRACE v1 complete: ",
                   sum(!is.na(data$CytoTRACE1_Score)), " cells scored."))
  } else {
    message("  [WARNING] CytoTRACE v1 produced no results.")
    RUN_CYTOTRACE1 <- FALSE
  }
  rm(counts_all, ct1_list); gc()
}

# =============================================================================
# --- STEP 3b: PATHWAY SCORES (AUCell over GOBP gene sets) --------------------
# =============================================================================
# Per-cell AUCell activity for the configured GO Biological Process gene sets.
# Adds AUCell_<name> columns and collects them in `pathway_cols`, which are
# folded into the by-condition plots and the per-contrast comparisons below.
# (pathway_cols was initialised near the top / RESUME block.)
if (RUN_PATHWAY_SCORES) {
  message("\n=== STEP 3b: AUCell pathway scores ===")
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    message("  [SKIP] AUCell not installed (BiocManager::install('AUCell')).")
  } else {
    present <- rownames(data)
    org_db  <- if (tolower(PATHWAY_SPECIES) %in% c("human", "hs")) "org.Hs.eg.db" else "org.Mm.eg.db"
    go_ok   <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
               requireNamespace(org_db, quietly = TRUE)

    build_set <- function(nm) {
      go_id <- PATHWAY_GO_TERMS[[nm]]
      genes <- character(0)
      if (go_ok && !is.null(go_id)) {
        genes <- tryCatch(
          unique(AnnotationDbi::select(getExportedValue(org_db, org_db), keys = go_id,
                                       keytype = "GOALL", columns = "SYMBOL")$SYMBOL),
          error = function(e) character(0))
        genes <- genes[!is.na(genes)]
      }
      if (length(intersect(genes, present)) < PATHWAY_MIN_GENES)
        genes <- PATHWAY_FALLBACK[[nm]]                       # curated fallback
      intersect(unique(genes), present)
    }

    gene_sets <- lapply(names(PATHWAY_GO_TERMS), build_set)
    names(gene_sets) <- names(PATHWAY_GO_TERMS)
    gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= PATHWAY_MIN_GENES]
    for (nm in names(gene_sets))
      message(sprintf("  %-22s %d genes", nm, length(gene_sets[[nm]])))

    if (length(gene_sets) > 0) {
      expr <- GetAssayData(data, assay = "RNA", layer = "counts")
      rk   <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
      auc  <- AUCell::AUCell_calcAUC(gene_sets, rk, verbose = FALSE)
      am   <- AUCell::getAUC(auc)                             # gene sets x cells
      for (nm in rownames(am)) {
        col <- paste0("AUCell_", nm)
        data@meta.data[[col]] <- as.numeric(am[nm, colnames(data)])
        pathway_cols <- c(pathway_cols, col)
      }
      message(paste0("  AUCell complete: ", length(pathway_cols), " pathway score(s) added."))
      rm(expr, rk, auc, am); gc()
    } else {
      message("  [SKIP] No gene sets had enough genes present.")
    }
  }
}

# =============================================================================
# --- STEP 3c: CCAT / SCENT (connectome correlation + signalling entropy) -----
# =============================================================================
# CCAT (Teschendorff) correlates each cell's expression with PPI hub degree; the
# PPI (net17Jan16.m from SCENT) is indexed by HUMAN Entrez, so mouse data must be
# mapped mouse-symbol -> human-symbol (homologene 10090->9606) -> human Entrez.
# This chain is fragile, so the WHOLE block is wrapped: any failure just skips
# CCAT/SCENT and the rest of Script 09 continues normally.
# (ccat_cols was initialised near the top / RESUME block.)
if (RUN_CCAT || RUN_SCENT) {
  message("\n=== STEP 3c: CCAT / SCENT ===")
  deps_ok <- requireNamespace("SCENT", quietly = TRUE) &&
             requireNamespace("AnnotationDbi", quietly = TRUE) &&
             (CCAT_SPECIES != "mouse" || requireNamespace("homologene", quietly = TRUE)) &&
             requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!deps_ok) {
    message("  [SKIP] Need SCENT + org.Hs.eg.db (+ homologene for mouse). ",
            "See INSTALL_NOTES.md (Potency benchmark methods). Continuing without CCAT/SCENT.")
  } else tryCatch({
    expr <- as.matrix(GetAssayData(data, assay = "RNA", layer = "counts"))
    mouse_syms <- rownames(expr)

    # 1) rownames -> HUMAN symbol
    if (CCAT_SPECIES == "mouse") {
      hom <- homologene::homologene(mouse_syms, inTax = 10090, outTax = 9606)
      hom <- hom[!is.na(hom[[2]]) & hom[[2]] != "", , drop = FALSE]
      hom <- hom[!duplicated(hom[[1]]), , drop = FALSE]        # unique source symbol
      sym2human <- stats::setNames(as.character(hom[[2]]), as.character(hom[[1]]))
      human_of <- sym2human[mouse_syms]
    } else {
      human_of <- stats::setNames(mouse_syms, mouse_syms)      # already human
    }

    # 2) HUMAN symbol -> Entrez
    hsym <- unique(human_of[!is.na(human_of)])
    ent  <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = hsym, keytype = "SYMBOL",
      column = "ENTREZID", multiVals = "first"))
    row_entrez <- ent[human_of]                                # per original row
    keep <- !is.na(row_entrez)
    expr_e <- expr[keep, , drop = FALSE]
    rownames(expr_e) <- as.character(row_entrez[keep])
    expr_e <- expr_e[!duplicated(rownames(expr_e)), , drop = FALSE]

    # 3) load the SCENT PPI (human Entrez) and take the shared genes
    utils::data("net17Jan16", package = "SCENT")
    net <- get("net17Jan16.m")
    common <- intersect(rownames(expr_e), rownames(net))
    if (length(common) < 100)
      stop("only ", length(common), " genes map to the PPI - aborting CCAT")
    exp_sub <- expr_e[common, , drop = FALSE]
    net_sub <- net[common, common, drop = FALSE]
    message(sprintf("  Mapped %d genes into the PPI (of %d).", length(common), nrow(expr)))

    # 4) integrate + CCAT (+ optional SCENT entropy)
    integ <- SCENT::DoIntegPPI(exp.m = exp_sub, ppiA.m = net_sub)
    if (RUN_CCAT) {
      ccat <- as.numeric(SCENT::CompCCAT(exp = integ$expMC, ppiA = integ$adjMC))
      if (length(ccat) == ncol(data)) {
        data$CCAT_score <- ccat
        ccat_cols <- c(ccat_cols, "CCAT_score")
        message("  CCAT complete (", ncol(data), " cells).")
      } else message("  [WARN] CCAT length mismatch; skipped.")
    }
    if (RUN_SCENT) {
      tryCatch({
        sr <- as.numeric(SCENT::CompSR(integ$expMC, integ$adjMC))
        if (length(sr) == ncol(data)) {
          data$SCENT_score <- sr
          ccat_cols <- c(ccat_cols, "SCENT_score")
          message("  SCENT signalling entropy complete.")
        }
      }, error = function(e)
        message("  [SKIP SCENT] CompSR failed: ", conditionMessage(e)))
    }
    rm(expr, expr_e); gc()
  }, error = function(e) {
    message("  [SKIP] CCAT/SCENT failed: ", conditionMessage(e),
            " -- continuing without it.")
  })
}

# =============================================================================
# --- CHECKPOINT: save scores + object BEFORE the heavy plotting ---------------
# =============================================================================
# Every per-cell score is computed by now. Save them immediately - a compact CSV
# of just the score columns AND the enriched .rds - so any failure (OOM, etc.) in
# the plotting/stats below never loses the expensive potency computation. To
# re-run only the plots later, set RESUME_SCORES <- TRUE at the top and this
# object already carries the columns.
message("\n=== Checkpoint: saving scores + enriched object ===")
.score_keep <- intersect(c(
  "entropy_shannon", "entropy_normalized", "entropy_score", "entropy_n_detected",
  "gene_counts_score", "CytoTRACE2_Score", "CytoTRACE2_Potency", "CytoTRACE2_Relative",
  "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency", "potency_score",
  "CytoTRACE1_Score", "CytoTRACE1_Rank", "CytoTRACE1_GCS",
  pathway_cols, ccat_cols), colnames(data@meta.data))
.id_keep <- intersect(c(SAMPLE_COLUMN, CELLTYPE_COLUMN, CONDITION_COLUMN),
                      colnames(data@meta.data))
tryCatch({
  utils::write.csv(
    data.frame(Barcode = colnames(data),
               data@meta.data[, c(.id_keep, .score_keep), drop = FALSE],
               check.names = FALSE),
    file.path(SCORES_DIR, "potency_scores_per_cell.csv"), row.names = FALSE)
  saveRDS(data, file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds")))
  message("  Saved potency_scores_per_cell.csv + ",
          PROJECT_NAME, "_with_cell_scores.rds (plotting next; safe to interrupt).")
}, error = function(e) message("  [WARNING] checkpoint save failed: ", e$message))

# =============================================================================
# --- STEP 4: METHOD CONCORDANCE ----------------------------------------------
# =============================================================================
# Every method above claims to order cells by differentiation state. If they
# genuinely measure the same thing they should correlate strongly. Spearman
# (rank) correlation is the right test: the methods share an ordering but not
# a common scale.
#
# Expected sign conventions:
#   CytoTRACE2_Score  HIGH = less differentiated (more potent)
#   CytoTRACE1_Score  HIGH = less differentiated
#   entropy_score     HIGH = less differentiated
#   gene_counts       HIGH = less differentiated
# So all pairwise correlations should be POSITIVE. A negative correlation is a
# red flag worth investigating before interpreting anything downstream.
# =============================================================================
message("\n=== STEP 4: Method concordance ===")

# Include any score column that is PRESENT on the object (whether freshly
# computed or reloaded via RESUME_SCORES) - a column exists only if it was scored.
score_cols <- c()
if ("potency_score"     %in% colnames(data@meta.data)) score_cols <- c(score_cols, "potency_score")
if ("CytoTRACE1_Score"  %in% colnames(data@meta.data)) score_cols <- c(score_cols, "CytoTRACE1_Score")
if ("entropy_score"     %in% colnames(data@meta.data)) score_cols <- c(score_cols, "entropy_score")
if ("gene_counts_score" %in% colnames(data@meta.data)) score_cols <- c(score_cols, "gene_counts_score")
if ("CCAT_score"        %in% colnames(data@meta.data)) score_cols <- c(score_cols, "CCAT_score")
if ("SCENT_score"       %in% colnames(data@meta.data)) score_cols <- c(score_cols, "SCENT_score")

cor_df <- NULL
if (length(score_cols) >= 2) {
  score_mat <- as.matrix(data@meta.data[, score_cols, drop = FALSE])
  cor_mat   <- suppressWarnings(
    stats::cor(score_mat, method = "spearman", use = "pairwise.complete.obs")
  )
  cor_df <- as.data.frame(cor_mat) %>% tibble::rownames_to_column("Method")

  message("  Spearman correlation between methods:")
  print(round(cor_mat, 3))

  neg <- which(cor_mat < 0 & upper.tri(cor_mat), arr.ind = TRUE)
  if (nrow(neg) > 0) {
    for (k in seq_len(nrow(neg))) {
      warning(paste0("  [CONCORDANCE] ", rownames(cor_mat)[neg[k, 1]], " and ",
                     colnames(cor_mat)[neg[k, 2]],
                     " are NEGATIVELY correlated (rho = ",
                     round(cor_mat[neg[k, 1], neg[k, 2]], 3),
                     "). These methods should agree in direction - investigate ",
                     "before interpreting potency results."))
    }
  }

  # Correlation heatmap
  tryCatch({
    cor_long <- as.data.frame(as.table(cor_mat))
    colnames(cor_long) <- c("Method1", "Method2", "rho")
    p_cor <- ggplot(cor_long, aes(x = Method1, y = Method2, fill = rho)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = round(rho, 2)), size = 4, fontface = "bold") +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(title    = "Concordance between potency / entropy methods",
           subtitle = "Spearman rho; all pairs should be POSITIVE",
           x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold"))
    ggsave(file.path(SCORES_DIR, "method_concordance_heatmap.png"),
           p_cor, width = 7, height = 6, dpi = DPI_SETTING, bg = "white")
  }, error = function(e) message(paste("  [WARNING] Concordance plot failed:", e$message)))

  rm(score_mat); gc()
} else {
  message("  Fewer than two score columns available - skipping concordance.")
}

# Pathway (AUCell) scores are activity levels, not differentiation ordering, so
# they were excluded from the concordance above but ARE included in every plot
# and comparison from here on.
if (length(pathway_cols) > 0) score_cols <- unique(c(score_cols, pathway_cols))

# =============================================================================
# --- STEP 5: VISUALISATION ---------------------------------------------------
# =============================================================================
message("\n=== STEP 5: Plots ===")

reduction_use <- if (UMAP_REDUCTION %in% names(data@reductions)) {
  UMAP_REDUCTION
} else if ("umap" %in% names(data@reductions)) {
  message(paste0("  [NOTE] '", UMAP_REDUCTION, "' not found; using 'umap'."))
  "umap"
} else {
  message("  [NOTE] No UMAP reduction found - UMAP plots will be skipped.")
  NA_character_
}

# --- 5a: Continuous score UMAPs ---------------------------------------------
if (!is.na(reduction_use)) {
  for (sc in score_cols) {
    tryCatch({
      p <- FeaturePlot(data, features = sc, reduction = reduction_use,
                       pt.size = POINT_SIZE) +
        scale_color_gradientn(colors = POTENCY_COLORS) +
        coord_fixed() +
        labs(title = sc, subtitle = "higher = less differentiated") +
        theme(plot.title = element_text(face = "bold"))
      ggsave(file.path(SCORES_DIR, paste0("umap_", sc, ".png")),
             p, width = PLOT_WIDTH, height = PLOT_HEIGHT,
             dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) {
      message(paste("  [WARNING] UMAP failed for", sc, ":", e$message))
    })
  }

  # --- 5b: Discrete potency category UMAP -----------------------------------
  if ("CytoTRACE2_Potency" %in% colnames(data@meta.data)) {
    tryCatch({
      p <- DimPlot(data, group.by = "CytoTRACE2_Potency",
                   reduction = reduction_use, pt.size = POINT_SIZE) +
        coord_fixed() +
        labs(title = "CytoTRACE 2 potency category") +
        theme(plot.title = element_text(face = "bold"))
      ggsave(file.path(SCORES_DIR, "umap_CytoTRACE2_Potency_category.png"),
             p, width = PLOT_WIDTH, height = PLOT_HEIGHT,
             dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) {
      message(paste("  [WARNING] Potency category UMAP failed:", e$message))
    })
  }
}

# --- 5c: Score by cell type --------------------------------------------------
# The key sanity check: known stem/progenitor compartments should sit at the
# top of these boxplots and terminally differentiated types at the bottom.
plot_by_group <- function(df, score, group, title, fname, angle = 45) {
  tryCatch({
    d <- df[!is.na(df[[score]]) & !is.na(df[[group]]), , drop = FALSE]
    if (nrow(d) == 0) return(invisible(NULL))
    # Order categories by median score so the gradient is readable at a glance.
    ord <- d %>% group_by(.data[[group]]) %>%
      summarise(m = stats::median(.data[[score]], na.rm = TRUE), .groups = "drop") %>%
      arrange(m)
    d[[group]] <- factor(as.character(d[[group]]), levels = as.character(ord[[group]]))

    p <- ggplot(d, aes(x = .data[[group]], y = .data[[score]], fill = .data[[group]])) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.6, linewidth = 0.3) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9, linewidth = 0.3) +
      labs(title = title, x = NULL, y = score) +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = angle, hjust = 1),
            plot.title  = element_text(face = "bold"))
    ggsave(file.path(SCORES_DIR, fname), p,
           width = PLOT_WIDTH + 2, height = PLOT_HEIGHT,
           dpi = DPI_SETTING, bg = "white")
    rm(p)
  }, error = function(e) {
    message(paste("  [WARNING] Plot", fname, "failed:", e$message))
  })
}

md <- data@meta.data
for (sc in score_cols) {
  plot_by_group(md, sc, CELLTYPE_COLUMN,
                paste0(sc, " by cell type"),
                paste0("box_", sc, "_by_celltype.png"))
}

# --- 5d: Score by condition, per cell type -----------------------------------
# Two views per score, ordered by CONDITION_LEVELS, with Wilcoxon significance
# brackets over each CONTRASTS_LIST pair:
#   barplot_<score>.png  - group MEAN bar + SE
#   violin_<score>.png   - violin + boxplot
generate_gene_comparison_plots <- function(seurat_obj, score_col, group_by, x_axis,
                                           comparisons, plot_type = "violin",
                                           output_prefix = "", plot_title = score_col,
                                           y_label = "Score",
                                           fig_width = 16, fig_height = 7,
                                           output_dir = SCORES_DIR) {
  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis)) %>%
    dplyr::rename(Expression = 1) %>% tidyr::drop_na()
  cust_theme <- theme_classic() + theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, face = "bold"),
    legend.position = "bottom", panel.spacing = unit(1.5, "lines")
  )
  p <- ggplot(df_plot, aes(!!sym(x_axis), Expression, fill = !!sym(x_axis)))
  if (plot_type == "barplot") {
    p <- p + stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.8) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  } else {
    p <- p + geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5)
  }
  p <- p +
    ggpubr::stat_compare_means(comparisons = comparisons, label = "p.signif",
                       method = "wilcox.test", method.args = list(exact = FALSE),
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols = c("****", "***", "**", "*", "ns")),
                       step.increase = 0.1, size = 6, bracket.size = 0.8) +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(title = plot_title, y = y_label) + scale_fill_brewer(palette = "Set1") + cust_theme
  ggsave(file.path(output_dir, paste0(output_prefix, score_col, ".png")),
         p, width = fig_width, height = fig_height, dpi = DPI_SETTING, bg = "white")
  invisible(p)
}

if (exists("CONTRASTS_LIST") && length(CONTRASTS_LIST) > 0) {
  ct_comparisons <- unname(CONTRASTS_LIST)   # list of c(group1, group2) pairs
  for (sc in score_cols) {
    for (pt in c("violin", "barplot")) {
      tryCatch(
        generate_gene_comparison_plots(
          data, score_col = sc, group_by = CELLTYPE_COLUMN, x_axis = CONDITION_COLUMN,
          comparisons = ct_comparisons, plot_type = pt,
          output_prefix = paste0(pt, "_"), plot_title = sc, y_label = sc),
        error = function(e)
          message("  [WARNING] ", pt, " plot failed for ", sc, ": ", e$message))
    }
  }
} else {
  message("  [NOTE] CONTRASTS_LIST empty - skipping the bracketed comparison plots.")
}

# --- 5e: Potency category composition stacked bars ---------------------------
if ("CytoTRACE2_Potency" %in% colnames(md)) {
  for (grp in unique(c(CELLTYPE_COLUMN, CONDITION_COLUMN))) {
    tryCatch({
      d <- md[!is.na(md$CytoTRACE2_Potency) & !is.na(md[[grp]]), , drop = FALSE]
      if (nrow(d) == 0) next
      comp <- d %>%
        group_by(.data[[grp]], CytoTRACE2_Potency) %>%
        summarise(n = dplyr::n(), .groups = "drop") %>%
        group_by(.data[[grp]]) %>%
        mutate(pct = n / sum(n) * 100) %>%
        ungroup()
      p <- ggplot(comp, aes(x = .data[[grp]], y = pct, fill = CytoTRACE2_Potency)) +
        geom_col(color = "white", linewidth = 0.2) +
        labs(title = paste0("Potency category composition by ", grp),
             x = NULL, y = "% of cells", fill = "Potency") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title  = element_text(face = "bold"))
      ggsave(file.path(SCORES_DIR, paste0("composition_potency_by_", grp, ".png")),
             p, width = PLOT_WIDTH + 2, height = PLOT_HEIGHT,
             dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) {
      message(paste("  [WARNING] Composition plot failed for", grp, ":", e$message))
    })
  }
}

# =============================================================================
# --- STEP 6: SUMMARY TABLES AND STATISTICS -----------------------------------
# =============================================================================
message("\n=== STEP 6: Summary tables and statistics ===")

sheets <- list()

# --- 6a: Per cell type summary ----------------------------------------------
if (length(score_cols) > 0) {
  summ_ct <- md %>%
    group_by(.data[[CELLTYPE_COLUMN]]) %>%
    summarise(
      N_Cells = dplyr::n(),
      dplyr::across(dplyr::all_of(score_cols),
                    list(mean   = ~mean(.x, na.rm = TRUE),
                         median = ~stats::median(.x, na.rm = TRUE),
                         sd     = ~stats::sd(.x, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(dplyr::desc(.data[[paste0(score_cols[1], "_median")]]))
  sheets[["By_CellType"]] <- as.data.frame(summ_ct)

  # --- 6b: Per cell type x condition -----------------------------------------
  summ_cc <- md %>%
    group_by(.data[[CELLTYPE_COLUMN]], .data[[CONDITION_COLUMN]]) %>%
    summarise(
      N_Cells = dplyr::n(),
      dplyr::across(dplyr::all_of(score_cols),
                    list(mean   = ~mean(.x, na.rm = TRUE),
                         median = ~stats::median(.x, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  sheets[["By_CellType_Condition"]] <- as.data.frame(summ_cc)

  # --- 6c: Per sample (QC view) ----------------------------------------------
  summ_s <- md %>%
    group_by(.data[[SAMPLE_COLUMN]]) %>%
    summarise(
      N_Cells = dplyr::n(),
      dplyr::across(dplyr::all_of(score_cols),
                    list(median = ~stats::median(.x, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  sheets[["By_Sample"]] <- as.data.frame(summ_s)
}

if (!is.null(cor_df)) sheets[["Method_Concordance"]] <- cor_df

# --- 6d: Group comparisons per cell type ------------------------------------
# Two-level CONDITION -> Wilcoxon rank-sum. More levels -> Kruskal-Wallis.
# Non-parametric throughout: potency scores are bounded and rarely normal.
#
# IMPORTANT CAVEAT recorded in the output: these tests treat individual cells
# as independent replicates. Cells from the same animal are not independent,
# so p-values are anti-conservative. Treat them as descriptive effect-size
# ranking, and confirm anything important with a sample-level test
# (e.g. per-sample medians compared across animals, n = number of animals).
if (RUN_GROUP_STATS && length(score_cols) > 0) {
  message("  Running group comparisons per cell type...")

  cond_vals <- unique(as.character(md[[CONDITION_COLUMN]]))
  cond_vals <- cond_vals[!is.na(cond_vals)]
  stat_rows <- list()

  for (sc in score_cols) {
    for (ct in unique(as.character(md[[CELLTYPE_COLUMN]]))) {
      if (is.na(ct)) next
      sub <- md[as.character(md[[CELLTYPE_COLUMN]]) == ct & !is.na(md[[sc]]), , drop = FALSE]
      if (nrow(sub) < STATS_MIN_CELLS) next

      grp_sizes <- table(as.character(sub[[CONDITION_COLUMN]]))
      grp_keep  <- names(grp_sizes)[grp_sizes >= STATS_MIN_CELLS]
      if (length(grp_keep) < 2) next
      sub <- sub[as.character(sub[[CONDITION_COLUMN]]) %in% grp_keep, , drop = FALSE]

      res <- tryCatch({
        if (length(grp_keep) == 2) {
          a <- sub[[sc]][as.character(sub[[CONDITION_COLUMN]]) == grp_keep[1]]
          b <- sub[[sc]][as.character(sub[[CONDITION_COLUMN]]) == grp_keep[2]]
          tt <- stats::wilcox.test(a, b)
          data.frame(
            Score       = sc,
            CellType    = ct,
            Test        = "Wilcoxon",
            Groups      = paste(grp_keep, collapse = " vs "),
            N_Total     = nrow(sub),
            Median_1    = stats::median(a, na.rm = TRUE),
            Median_2    = stats::median(b, na.rm = TRUE),
            Delta       = stats::median(a, na.rm = TRUE) - stats::median(b, na.rm = TRUE),
            Statistic   = unname(tt$statistic),
            P_Value     = tt$p.value,
            stringsAsFactors = FALSE
          )
        } else {
          tt <- stats::kruskal.test(
            stats::as.formula(paste0("`", sc, "` ~ `", CONDITION_COLUMN, "`")),
            data = sub
          )
          data.frame(
            Score       = sc,
            CellType    = ct,
            Test        = "Kruskal-Wallis",
            Groups      = paste(grp_keep, collapse = ", "),
            N_Total     = nrow(sub),
            Median_1    = NA_real_,
            Median_2    = NA_real_,
            Delta       = NA_real_,
            Statistic   = unname(tt$statistic),
            P_Value     = tt$p.value,
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)

      if (!is.null(res)) stat_rows[[paste(sc, ct, sep = "|")]] <- res
    }
  }

  if (length(stat_rows) > 0) {
    stats_df <- do.call(rbind, stat_rows)
    # Correct across cell types WITHIN each score, not across everything at
    # once - the scores are different questions, not one family of tests.
    stats_df <- stats_df %>%
      group_by(Score) %>%
      mutate(P_Adj = stats::p.adjust(P_Value, method = STATS_PADJ_METHOD)) %>%
      ungroup() %>%
      arrange(Score, P_Adj) %>%
      as.data.frame()
    stats_df$Significant <- stats_df$P_Adj < 0.05
    stats_df$CAVEAT <- "Cell-level test; cells within a sample are not independent. Confirm with sample-level statistics."
    sheets[["Group_Comparisons"]] <- stats_df

    n_sig <- sum(stats_df$Significant, na.rm = TRUE)
    message(paste0("  ", nrow(stats_df), " comparisons | ", n_sig,
                   " significant after ", STATS_PADJ_METHOD, " correction"))
  } else {
    message("  No comparisons met the minimum cell requirements.")
  }
}

# --- 6d-2: Per-contrast comparisons (CONTRASTS_LIST) -------------------------
# For each defined contrast (a specific pair of CONDITION groups) and each cell
# type, a Wilcoxon test on every score (potency/entropy + AUCell pathways). This
# is the focused, pairwise version of the omnibus test above, and is the one to
# read for the KO-vs-polyp / polyp-vs-WT questions.
if (RUN_GROUP_STATS && exists("CONTRASTS_LIST") && length(CONTRASTS_LIST) > 0 &&
    length(score_cols) > 0) {
  message("  Running per-contrast comparisons (", length(CONTRASTS_LIST), " contrasts)...")
  con_rows <- list()
  for (cn in names(CONTRASTS_LIST)) {
    g1 <- CONTRASTS_LIST[[cn]][1]; g2 <- CONTRASTS_LIST[[cn]][2]
    for (sc in score_cols) {
      for (ct in unique(as.character(md[[CELLTYPE_COLUMN]]))) {
        if (is.na(ct)) next
        sub <- md[as.character(md[[CELLTYPE_COLUMN]]) == ct & !is.na(md[[sc]]) &
                  as.character(md[[CONDITION_COLUMN]]) %in% c(g1, g2), , drop = FALSE]
        a <- sub[[sc]][as.character(sub[[CONDITION_COLUMN]]) == g1]
        b <- sub[[sc]][as.character(sub[[CONDITION_COLUMN]]) == g2]
        if (length(a) < STATS_MIN_CELLS || length(b) < STATS_MIN_CELLS) next
        tt <- tryCatch(stats::wilcox.test(a, b), error = function(e) NULL)
        if (is.null(tt)) next
        con_rows[[paste(cn, sc, ct, sep = "|")]] <- data.frame(
          Contrast = cn, Score = sc, CellType = ct,
          Group1 = g1, Group2 = g2, N1 = length(a), N2 = length(b),
          Median_1 = stats::median(a, na.rm = TRUE),
          Median_2 = stats::median(b, na.rm = TRUE),
          Delta    = stats::median(a, na.rm = TRUE) - stats::median(b, na.rm = TRUE),
          Statistic = unname(tt$statistic), P_Value = tt$p.value,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(con_rows) > 0) {
    con_df <- do.call(rbind, con_rows)
    # BH-correct within each Contrast x Score family (across cell types).
    con_df <- con_df %>%
      group_by(Contrast, Score) %>%
      mutate(P_Adj = stats::p.adjust(P_Value, method = STATS_PADJ_METHOD)) %>%
      ungroup() %>%
      arrange(Contrast, Score, P_Adj) %>%
      as.data.frame()
    con_df$Significant <- con_df$P_Adj < 0.05
    con_df$CAVEAT <- "Cell-level Wilcoxon; cells within a sample are not independent. Confirm at the sample level."
    sheets[["Contrast_Comparisons"]] <- con_df
    message(paste0("  ", nrow(con_df), " per-contrast comparisons | ",
                   sum(con_df$Significant, na.rm = TRUE), " significant after ",
                   STATS_PADJ_METHOD, "."))
  } else {
    message("  No per-contrast comparisons met the minimum cell requirements.")
  }
}

# --- 6e: Run metadata sheet --------------------------------------------------
sheets[["Run_Info"]] <- data.frame(
  Parameter = c("PROJECT_NAME", "RDS_PATH", "N_Cells", "N_Genes",
                "CELLTYPE_COLUMN", "CONDITION_COLUMN", "SAMPLE_COLUMN",
                "RUN_CYTOTRACE2", "RUN_CYTOTRACE1", "RUN_ENTROPY",
                "CT2_SPECIES", "CT2_SLOT",
                "CT2_USE_PREKNN", "ENTROPY_NORMALIZE", "ENTROPY_MIN_GENES",
                "Pct_Cells_Under_1000_Genes", "Date"),
  Value = c(PROJECT_NAME, RDS_PATH, ncol(data), nrow(data),
            CELLTYPE_COLUMN, CONDITION_COLUMN, SAMPLE_COLUMN,
            RUN_CYTOTRACE2, RUN_CYTOTRACE1, RUN_ENTROPY,
            CT2_SPECIES, CT2_SLOT,
            CT2_USE_PREKNN, ENTROPY_NORMALIZE, ENTROPY_MIN_GENES,
            round(frac_low * 100, 2), as.character(Sys.Date())),
  stringsAsFactors = FALSE
)

if (length(sheets) > 0) {
  xlsx_path <- file.path(SCORES_DIR, "cell_scores_summary.xlsx")
  write_xlsx(sheets, xlsx_path)
  message(paste("  Summary written to:", basename(xlsx_path)))
}

# --- 6f: Full per-cell table -------------------------------------------------
# Compressed CSV so the raw values are available for custom analysis without
# reloading the (large) Seurat object.
tryCatch({
  keep_cols <- unique(c(SAMPLE_COLUMN, CELLTYPE_COLUMN, CONDITION_COLUMN,
                        score_cols,
                        intersect(c("CytoTRACE2_Score", "CytoTRACE2_Potency",
                                    "CytoTRACE2_Relative",
                                    "preKNN_CytoTRACE2_Score",
                                    "preKNN_CytoTRACE2_Potency",
                                    "CytoTRACE1_Score", "CytoTRACE1_Rank",
                                    "CytoTRACE1_GCS",
                                    "entropy_shannon", "entropy_normalized",
                                    "entropy_n_detected"),
                                  colnames(md))))
  per_cell <- data.frame(Barcode = rownames(md), md[, keep_cols, drop = FALSE],
                         stringsAsFactors = FALSE)
  gz <- gzfile(file.path(SCORES_DIR, "cell_scores_per_cell.csv.gz"), "w")
  utils::write.csv(per_cell, gz, row.names = FALSE)
  close(gz)
  message("  Per-cell table written to: cell_scores_per_cell.csv.gz")
}, error = function(e) {
  message(paste("  [WARNING] Per-cell table failed:", e$message))
})

# =============================================================================
# --- STEP 7: SAVE ENRICHED OBJECT --------------------------------------------
# =============================================================================
out_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))
message(paste("\n=== STEP 7: Saving object to", basename(out_rds), "==="))
saveRDS(data, out_rds)

message("\n=== Script 09 complete ===")
message(paste0("  Scores added: ", paste(score_cols, collapse = ", ")))
message(paste0("  Plots and tables: ", SCORES_DIR))
message(paste0("  Enriched object:  ", out_rds))
message("\n  INTERPRETATION REMINDERS:")
message("   - CytoTRACE2_Score is ABSOLUTE and comparable across datasets.")
message("   - CytoTRACE2_Relative and CytoTRACE1_Score are RELATIVE to this run only.")
message("   - Check method_concordance_heatmap.png: all correlations should be positive.")
message("   - Group p-values are cell-level; confirm key findings at the sample level.")
message("\n  NEXT: 10_trajectory_cellrank.R to infer trajectories using these scores.")
