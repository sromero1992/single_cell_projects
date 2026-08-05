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
#   The CytoTRACE 2 authors recommend running SEPARATELY over each dataset
#   rather than integrating first, because the post-processing step smooths
#   predictions using neighbouring cells and can therefore be distorted by
#   batch effects. Absolute outputs (CytoTRACE2_Score, _Potency) are calibrated
#   and remain comparable across runs; CytoTRACE2_Relative is NOT comparable
#   across runs, as it is rescaled within each input.
#
#   RUN_PER_SAMPLE = TRUE (default) honours this: each SampleID is scored
#   independently and results are stitched back together. Set FALSE only if
#   you have few cells per sample and accept the caveat.
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
library(patchwork)
library(writexl)
library(Matrix)

set.seed(123)

# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) --------------------------
# =============================================================================

# --- 1.1: Project Identity & Paths -------------------------------------------
PROJECT_NAME <- "Wu_Diet_project2"
ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_wu_project2/r_process"
OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")

# Input object — the final annotated object from Script 06.
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))

# Or point at a subtype object to score within one lineage:
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_tcells_subclustered.rds"))
#RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))

SCORES_DIR <- file.path(OUTPUT_DIR, "cell_scores")

# --- 1.2: Metadata Columns ---------------------------------------------------
CELLTYPE_COLUMN  <- "CellType"        # "sub_cell_types" for subtype objects
SAMPLE_COLUMN    <- "SampleID"
CONDITION_COLUMN <- "Diet"            # Wu diet study grouping variable

# Order for plot axes. Put the control/reference level FIRST.
# Set to NULL to use whatever order the factor already has.
CONDITION_LEVELS <- NULL
# e.g. CONDITION_LEVELS <- c("WT_Female", "Polyp_Female", "Polyp_NR4a1_KO_Female")

# Reduction used for score UMAPs. Script 01 produces "umap_harmony".
UMAP_REDUCTION <- "umap_harmony"

# --- 1.3: Which Methods To Run -----------------------------------------------
RUN_CYTOTRACE2 <- TRUE   # absolute potency (recommended primary method)
RUN_CYTOTRACE1 <- TRUE   # original relative CytoTRACE; skipped if not installed
RUN_ENTROPY    <- TRUE   # native entropy metrics; no dependencies, always works

# --- 1.4: CytoTRACE 2 Parameters ---------------------------------------------
# SPECIES: "mouse" or "human". CytoTRACE 2's feature set is mouse-based; for
#   human input it performs orthology mapping internally.
CT2_SPECIES <- "mouse"

# RUN_PER_SAMPLE: score each SampleID separately (recommended — see header).
CT2_RUN_PER_SAMPLE <- TRUE

# CT2_SLOT: which layer to feed. MUST be raw/CPM counts, never log-transformed.
CT2_SLOT <- "counts"

# CT2_NCORES: parallel cores. The authors advise 1-2 on machines with <16 GB
#   RAM, since each worker holds a copy of the expression block.
CT2_NCORES <- 4

# Batch sizes. Defaults follow the package recommendations:
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

  # Container for results, filled either per-sample or in one pass.
  ct2_cols <- c("CytoTRACE2_Score", "CytoTRACE2_Potency", "CytoTRACE2_Relative",
                "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency")
  ct2_all  <- NULL

  # ---- Helper: build the exact input CytoTRACE 2 expects --------------------
  # Validated on the scDVEP benchmark (see run_cytotrace2.R, NOTES_debugging.txt
  # and cytotrace2_linux_setup.docx). Three things are enforced here because
  # each one, when wrong, silently collapses every score:
  #   1. DENSE data.frame (genes x cells). Sparse matrices cause SILENT failures.
  #   2. Gene SYMBOL rownames. If the object carries Ensembl IDs they are mapped
  #      to symbols (mouse: org.Mm.eg.db, human: org.Hs.eg.db); otherwise
  #      CytoTRACE 2's internal gene panel matches almost nothing.
  #   3. Deduplicated rownames (keep the highest mean-expression row per symbol).
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

    as.data.frame(as.matrix(m))   # DENSE data.frame - required
  }

  # ---- Helper: run CytoTRACE 2 on one Seurat object ------------------------
  # Uses the VALIDATED argument set. Do NOT re-enable the parallelize_* flags or
  # add verbose = TRUE: parallelize_smoothing = TRUE reintroduces socket-cluster
  # hangs and 'verbose' is not a valid argument ("unused argument" error). Both
  # were the source of repeated failed runs (see NOTES_debugging.txt section 1).
  run_ct2_block <- function(obj, label) {
    if (ncol(obj) < CT2_MIN_CELLS) {
      message(paste0("    [SKIP] ", label, ": only ", ncol(obj),
                     " cells (< CT2_MIN_CELLS = ", CT2_MIN_CELLS, ")"))
      return(NULL)
    }
    message(paste0("    -> ", label, ": ", ncol(obj), " cells..."))

    res <- tryCatch({
      # is_seurat = FALSE: feed a plain DENSE data.frame, not a Seurat object.
      ct2_input <- prep_ct2_input(obj)
      out <- cytotrace2(
        ct2_input,
        is_seurat             = FALSE,
        species               = CT2_SPECIES,
        batch_size            = CT2_BATCH_SIZE,
        smooth_batch_size     = CT2_SMOOTH_BATCH_SIZE,
        parallelize_models    = FALSE,   # validated: avoids parallel backend issues
        parallelize_smoothing = FALSE,   # validated: REQUIRED, avoids socket hangs
        ncores                = CT2_NCORES,
        seed                  = 14
      )
      # is_seurat = FALSE returns a data.frame: rows = cells, cols = score types.
      keep <- intersect(ct2_cols, colnames(out))
      if (length(keep) == 0) {
        stop("cytotrace2() returned no recognised prediction columns.")
      }
      res_df <- out[, keep, drop = FALSE]
      # Clip the continuous scores to [0,1] (validated guard against tiny spill).
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
    # Per the CytoTRACE 2 FAQ: run separately per dataset/batch, because the
    # KNN post-processing borrows information across neighbouring cells and is
    # therefore sensitive to batch structure.
    samples <- unique(as.character(data@meta.data[[SAMPLE_COLUMN]]))
    message(paste0("  Running per sample (", length(samples), " samples)..."))

    res_list <- list()
    for (s in samples) {
      cells_s <- colnames(data)[as.character(data@meta.data[[SAMPLE_COLUMN]]) == s]
      obj_s   <- subset(data, cells = cells_s)
      r <- run_ct2_block(obj_s, s)
      if (!is.null(r)) res_list[[s]] <- r
      rm(obj_s); gc()
    }

    if (length(res_list) > 0) {
      # Harmonise columns across samples before stacking (a failed optional
      # column in one sample must not silently drop it everywhere).
      common <- Reduce(intersect, lapply(res_list, colnames))
      res_list <- lapply(res_list, function(d) d[, common, drop = FALSE])
      ct2_all  <- do.call(rbind, res_list)
    }
    rm(res_list); gc()

  } else {
    message("  Running on the full object in one pass...")
    ct2_all <- run_ct2_block(data, "ALL")
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
# Seurat object, and can be memory-hungry. It is run per sample for the same
# batch-effect reason as CytoTRACE 2.
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

score_cols <- c()
if (RUN_CYTOTRACE2 && "potency_score"     %in% colnames(data@meta.data)) score_cols <- c(score_cols, "potency_score")
if (RUN_CYTOTRACE1 && "CytoTRACE1_Score"  %in% colnames(data@meta.data)) score_cols <- c(score_cols, "CytoTRACE1_Score")
if (RUN_ENTROPY    && "entropy_score"     %in% colnames(data@meta.data)) score_cols <- c(score_cols, "entropy_score")
if (RUN_ENTROPY    && "gene_counts_score" %in% colnames(data@meta.data)) score_cols <- c(score_cols, "gene_counts_score")

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
  plot_by_group(md, sc, CONDITION_COLUMN,
                paste0(sc, " by ", CONDITION_COLUMN),
                paste0("box_", sc, "_by_condition.png"))
}

# --- 5d: Score by condition, split by cell type ------------------------------
# The comparison that usually matters: within each cell type, does potency
# shift between experimental groups?
for (sc in score_cols) {
  tryCatch({
    d <- md[!is.na(md[[sc]]), , drop = FALSE]
    if (nrow(d) == 0) next
    p <- ggplot(d, aes(x = .data[[CONDITION_COLUMN]], y = .data[[sc]],
                       fill = .data[[CONDITION_COLUMN]])) +
      geom_boxplot(outlier.size = 0.15, linewidth = 0.3) +
      facet_wrap(stats::as.formula(paste0("~ `", CELLTYPE_COLUMN, "`")),
                 scales = "free_y") +
      labs(title = paste0(sc, " by condition, per cell type"),
           x = NULL, y = sc) +
      theme_bw() +
      theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "bottom",
            strip.text    = element_text(size = 7, face = "bold"),
            plot.title    = element_text(face = "bold"))
    ggsave(file.path(SCORES_DIR, paste0("facet_", sc, "_condition_by_celltype.png")),
           p, width = 14, height = 10, dpi = DPI_SETTING, bg = "white")
    rm(p)
  }, error = function(e) {
    message(paste("  [WARNING] Facet plot failed for", sc, ":", e$message))
  })
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

# --- 6e: Run metadata sheet --------------------------------------------------
sheets[["Run_Info"]] <- data.frame(
  Parameter = c("PROJECT_NAME", "RDS_PATH", "N_Cells", "N_Genes",
                "CELLTYPE_COLUMN", "CONDITION_COLUMN", "SAMPLE_COLUMN",
                "RUN_CYTOTRACE2", "RUN_CYTOTRACE1", "RUN_ENTROPY",
                "CT2_SPECIES", "CT2_RUN_PER_SAMPLE", "CT2_SLOT",
                "CT2_USE_PREKNN", "ENTROPY_NORMALIZE", "ENTROPY_MIN_GENES",
                "Pct_Cells_Under_1000_Genes", "Date"),
  Value = c(PROJECT_NAME, RDS_PATH, ncol(data), nrow(data),
            CELLTYPE_COLUMN, CONDITION_COLUMN, SAMPLE_COLUMN,
            RUN_CYTOTRACE2, RUN_CYTOTRACE1, RUN_ENTROPY,
            CT2_SPECIES, CT2_RUN_PER_SAMPLE, CT2_SLOT,
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
