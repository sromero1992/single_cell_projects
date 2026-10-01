# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 10: TRAJECTORY INFERENCE (CellRank CytoTRACEKernel)
# Version: 1.0 (UNIFIED build)
#
# UNIFIED BUILD: part of unified_pipeline/. Consumes the scored object from
#   09_cell_scores.R and infers directed trajectories from it.
#
# PURPOSE:
#   Script 09 answers "how potent is each cell?". This script answers the next
#   question: "given that potency gradient, where is each cell GOING?"
#
#   CellRank models differentiation as a Markov chain over cells. Each kernel
#   supplies a different notion of directionality; the CytoTRACEKernel derives
#   it from a CytoTRACE-style score (number of expressed genes per cell),
#   which means it works WITHOUT RNA velocity and therefore without spliced /
#   unspliced counts. That matters here: standard Cell Ranger output has no
#   velocity layers, so velocity-based kernels are not an option for this data
#   without re-running alignment through velocyto or kallisto.
#
#   From the transition matrix CellRank computes:
#     - terminal states (endpoints of differentiation)
#     - initial states (progenitor compartments)
#     - fate probabilities: for each cell, the probability of ending in each
#       terminal state
#     - a pseudotime ordering
#
# ============================================================================
# WHY THIS SCRIPT USES PYTHON THROUGH reticulate
# ============================================================================
#   CellRank and scanpy are Python-only; there is no R port. reticulate embeds
#   a Python interpreter inside the R session so the whole workflow stays in
#   one script and one object, rather than splitting into a manual R -> h5ad ->
#   Jupyter -> R round trip.
#
#   The Python environment is NOT created automatically here, because silently
#   building conda environments from inside an analysis script is a good way to
#   corrupt a working setup. Create it ONCE by following INSTALL_NOTES.md
#   section 5, then point PYTHON_ENV below at it.
#
#   Note the split: CytoTRACE 2 (Script 09) runs natively in R and needs no
#   Python. Only this trajectory script requires the Python stack.
#
#   LAUNCH NOTE (shared C++ runtime): system R may ship an older libstdc++ than
#   the conda env's compiled extensions (llvmlite/numba) need, giving a
#   "GLIBCXX_3.4.30 not found" error on import. Launch R with the env's libstdc++
#   preloaded so it is used from process start:
#     conda activate scanpy_env_311
#     LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6 Rscript 10_trajectory_cellrank.R
# ============================================================================
#
# INPUT:
#   <PROJECT>_with_cell_scores.rds   (from Script 09)
#
# OUTPUT:
#   - <PROJECT>_with_trajectory.rds  (object + pseudotime + fate probabilities)
#   - trajectory/ plots and tables
#   - trajectory/adata_cellrank.h5ad (the AnnData object, for Python follow-up)
#
# =============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)
library(Matrix)
library(reticulate)

# --- Serial / single-threaded execution --------------------------------------
# Force BLAS / OpenMP / numba to one thread so Python never spins up ~20 workers
# that duplicate large arrays (a memory-safety measure, not a speed one). These
# MUST be set BEFORE numpy/numba are imported, so keep this block at the very top.
# Set to a small number (e.g. "4") instead of "1" if you want some parallelism.
.threads <- "1"
Sys.setenv(OMP_NUM_THREADS        = .threads,
           OPENBLAS_NUM_THREADS   = .threads,
           MKL_NUM_THREADS        = .threads,
           NUMEXPR_NUM_THREADS    = .threads,
           NUMBA_NUM_THREADS      = .threads,
           VECLIB_MAXIMUM_THREADS = .threads)

set.seed(123)
# =============================================================================
# --- PART 1: USER CONFIGURATION (EDIT THIS SECTION) --------------------------
# =============================================================================
# --- 1.1: Python environment -------------------------------------------------
# PYTHON_ENV_TYPE: how to locate Python.
#   "conda"      — a named conda/mamba environment (recommended)
#   "virtualenv" — a virtualenv name or path
#   "python"     — a direct path to a python binary
#   "auto"       — let reticulate decide (least reproducible; use as fallback)
PYTHON_ENV_TYPE <- "conda"
# Name of the environment. Point this at an existing conda env that has
# scanpy + cellrank installed. See INSTALL_NOTES.md section 5. Defaults to the
# config.R value (PY_ENV_CELLRANK, overridable via NR4A1_PY_CELLRANK).
PYTHON_ENV <- if (exists("PY_ENV_CELLRANK")) PY_ENV_CELLRANK else "scanpy_env_311"
# For PYTHON_ENV_TYPE = "python", set the full binary path instead, e.g.:
# PYTHON_BIN <- "/home/ssromerogon/miniconda3/envs/cellrank_env/bin/python"
PYTHON_BIN <- NULL
# --- 1.2: Project Identity & Paths -------------------------------------------
# ---- Shared portable config: nr4a1 defaults, override via env vars (see config.R).
# Run from the pipeline directory, or set NR4A1_CONFIG=/full/path/to/config.R. ----
.NR4A1_CFG <- Sys.getenv("NR4A1_CONFIG", "config.R")
if (!file.exists(.NR4A1_CFG)) stop("config.R not found at '", .NR4A1_CFG,
  "' - cd to the pipeline directory or set NR4A1_CONFIG.", call. = FALSE)
source(.NR4A1_CFG)

# PROJECT_NAME <- "Nr4a1_s17_ack"   # [portable] now set in config.R
# ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"   # [portable] now set in config.R
# OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")   # [portable] now set in config.R
# Input object — the colonocyte sub-clustered/sub-annotated object.
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_colonocytes_subclustered.rds"))
TRAJ_DIR <- file.path(OUTPUT_DIR, "trajectory")
# --- 1.3: Metadata columns ---------------------------------------------------
CELLTYPE_COLUMN  <- "CellType"
SAMPLE_COLUMN    <- "SampleID"
CONDITION_COLUMN <- "Genotype_sex"
# --- 1.4: Subsetting ---------------------------------------------------------
# Trajectory inference is only meaningful within a connected differentiation
# hierarchy. Running it across unrelated lineages (T cells AND colonocytes AND
# macrophages at once) produces transitions between cell types that cannot
# actually interconvert, and the result is meaningless.
#
# SUBSET_TO: character vector of values in CELLTYPE_COLUMN to keep, or NULL
#   for all cells. STRONGLY recommended to set this to one lineage.
SUBSET_TO <- NULL
# e.g. epithelial hierarchy:
# SUBSET_TO <- c("Stem cells", "TA cells", "Colonocytes", "Goblet cells")
# Optionally restrict to one experimental group, so trajectories are not
# distorted by mixing genotypes. NULL = keep all.
SUBSET_CONDITION <- NULL
# e.g. SUBSET_CONDITION <- "WT_Female"
# MAX_CELLS: downsample above this many cells. CellRank's dense operations grow
#   quickly; 30-50K is comfortable on a typical workstation. NULL = no limit.
MAX_CELLS <- NULL
# --- 1.5: Preprocessing for the kernel ---------------------------------------
N_HVG        <- 2000   # highly variable genes for the PCA/neighbour graph
N_PCS        <- 20     # principal components
N_NEIGHBORS  <- 30     # kNN graph size; also used for the imputation step
# USE_HARMONY_EMBEDDING: reuse the batch-corrected Harmony embedding from
#   Script 01 instead of recomputing PCA inside Python. Recommended when the
#   trajectory spans multiple samples, since uncorrected PCA will place cells
#   from different animals on separate branches and CellRank will interpret
#   that batch structure as biology.
USE_HARMONY_EMBEDDING <- TRUE
HARMONY_REDUCTION     <- "harmony"
# --- 1.6: CytoTRACEKernel parameters -----------------------------------------
# CellRank's CytoTRACEKernel needs a smoothed ("imputed") expression layer,
# because raw counts are too sparse for a stable gene-count correlation.
#
# CT_IMPUTE_METHOD:
#   "moments"    — kNN averaging over the neighbour graph (built directly here;
#                  scvelo's moments needs spliced/unspliced, which plain scRNA
#                  lacks). Well established, simple.
#   "cellmapper" — the newer CellMapper-based imputation shown in current
#                  CellRank docs. Requires the cellmapper package.
#   "none"       — use the log-normalised matrix directly. Fastest, noisiest;
#                  acceptable for a quick look only.
CT_IMPUTE_METHOD <- "moments"
# THRESHOLD for calling terminal states.
N_STATES <- 5          # number of macrostates to compute
N_TERMINAL_STATES <- 3 # how many of those to designate terminal
# CLUSTER_KEY used by CellRank when naming macrostates.
CLUSTER_KEY <- CELLTYPE_COLUMN
# --- 1.7: Compute ------------------------------------------------------------
N_JOBS <- 1            # threads for CellRank/scanpy (1 = serial, memory-safe)
RANDOM_SEED <- 123
# --- 1.8: Plotting -----------------------------------------------------------
DPI_SETTING <- 300
SAVE_H5AD   <- TRUE    # keep the AnnData for follow-up work in Python
# =============================================================================
# --- PART 2: EXECUTION (DO NOT EDIT BELOW THIS LINE) -------------------------
# =============================================================================
if (!dir.exists(TRAJ_DIR)) dir.create(TRAJ_DIR, recursive = TRUE)
# =============================================================================
# --- STEP 1: BIND THE PYTHON ENVIRONMENT -------------------------------------
# =============================================================================
# This must happen BEFORE any Python module is touched. reticulate binds an
# interpreter on first use and cannot switch afterwards without restarting R.
message("=== STEP 1: Binding Python environment ===")
bind_python <- function() {
  ok <- tryCatch({
    if (PYTHON_ENV_TYPE == "conda") {
      reticulate::use_condaenv(PYTHON_ENV, required = TRUE)
    } else if (PYTHON_ENV_TYPE == "virtualenv") {
      reticulate::use_virtualenv(PYTHON_ENV, required = TRUE)
    } else if (PYTHON_ENV_TYPE == "python") {
      if (is.null(PYTHON_BIN)) stop("PYTHON_ENV_TYPE is 'python' but PYTHON_BIN is NULL.")
      reticulate::use_python(PYTHON_BIN, required = TRUE)
    } else if (PYTHON_ENV_TYPE == "auto") {
      message("  [NOTE] PYTHON_ENV_TYPE = 'auto'. Reproducibility is not guaranteed.")
    } else {
      stop(paste0("Unknown PYTHON_ENV_TYPE: '", PYTHON_ENV_TYPE, "'"))
    }
    TRUE
  }, error = function(e) {
    message(paste("  [ERROR] Could not bind Python environment:", e$message))
    FALSE
  })
  ok
}
if (!bind_python()) {
  stop(paste0(
    "\n\nFailed to bind the Python environment '", PYTHON_ENV, "'.\n",
    "Create it once by following INSTALL_NOTES.md section 5, then re-run.\n",
    "Quick check from a terminal:\n",
    "    conda env list\n",
    "    conda activate ", PYTHON_ENV, " && python -c 'import cellrank; print(cellrank.__version__)'\n"
  ), call. = FALSE)
}
py_cfg <- reticulate::py_config()
message(paste("  Python:", py_cfg$python))
message(paste("  Version:", py_cfg$version))
# --- Verify every required module is importable BEFORE doing any work -------
required_modules <- c("scanpy", "anndata", "cellrank", "numpy", "scipy", "pandas")
if (CT_IMPUTE_METHOD == "cellmapper") required_modules <- c(required_modules, "cellmapper")
missing_modules <- required_modules[!vapply(required_modules,
                                            reticulate::py_module_available,
                                            logical(1))]
if (length(missing_modules) > 0) {
  stop(paste0(
    "\n\nMissing Python modules in '", PYTHON_ENV, "':\n  - ",
    paste(missing_modules, collapse = "\n  - "),
    "\n\nInstall them into that environment, e.g.:\n",
    "    conda activate ", PYTHON_ENV, "\n",
    "    pip install ", paste(missing_modules, collapse = " "), "\n",
    "\nSee INSTALL_NOTES.md section 5.\n"
  ), call. = FALSE)
}
sc  <- reticulate::import("scanpy",   convert = FALSE)
ad  <- reticulate::import("anndata",  convert = FALSE)
cr  <- reticulate::import("cellrank", convert = FALSE)
np  <- reticulate::import("numpy",    convert = FALSE)
sp  <- reticulate::import("scipy.sparse", convert = FALSE)
sc$settings$n_jobs <- as.integer(N_JOBS)   # keep scanpy/CellRank serial too
message(paste("  scanpy  :", reticulate::py_to_r(sc$`__version__`)))
message(paste("  cellrank:", reticulate::py_to_r(cr$`__version__`)))
# =============================================================================
# --- STEP 2: LOAD AND SUBSET THE SEURAT OBJECT -------------------------------
# =============================================================================
message("\n=== STEP 2: Loading Seurat object ===")
if (!file.exists(RDS_PATH)) {
  stop(paste0("Input .rds not found:\n  ", RDS_PATH,
              "\nRun Script 09 first, or correct RDS_PATH."), call. = FALSE)
}
data <- readRDS(RDS_PATH)
DefaultAssay(data) <- "RNA"
if (inherits(data[["RNA"]], "Assay5") &&
    length(Layers(data[["RNA"]], search = "counts")) > 1) {
  data <- JoinLayers(data)
}
message(paste0("  Loaded: ", ncol(data), " cells"))
# --- Subset to a coherent lineage -------------------------------------------
if (!is.null(SUBSET_TO)) {
  present <- unique(as.character(data@meta.data[[CELLTYPE_COLUMN]]))
  unknown <- setdiff(SUBSET_TO, present)
  if (length(unknown) > 0) {
    warning(paste0("SUBSET_TO values not present in the data: ",
                   paste(unknown, collapse = ", ")))
  }
  keep <- colnames(data)[as.character(data@meta.data[[CELLTYPE_COLUMN]]) %in% SUBSET_TO]
  if (length(keep) == 0) stop("SUBSET_TO matched zero cells.", call. = FALSE)
  data <- subset(data, cells = keep)
  message(paste0("  Subset to ", length(SUBSET_TO), " cell type(s): ",
                 ncol(data), " cells"))
} else {
  message("  [WARNING] SUBSET_TO is NULL - running across ALL cell types.")
  message("            Trajectories between unrelated lineages are not meaningful.")
  message("            Consider setting SUBSET_TO to one differentiation hierarchy.")
}
if (!is.null(SUBSET_CONDITION)) {
  keep <- colnames(data)[as.character(data@meta.data[[CONDITION_COLUMN]]) %in% SUBSET_CONDITION]
  if (length(keep) == 0) stop("SUBSET_CONDITION matched zero cells.", call. = FALSE)
  data <- subset(data, cells = keep)
  message(paste0("  Subset to condition(s) ",
                 paste(SUBSET_CONDITION, collapse = ", "), ": ", ncol(data), " cells"))
}
# --- Downsample if needed ----------------------------------------------------
if (!is.null(MAX_CELLS) && ncol(data) > MAX_CELLS) {
  message(paste0("  Downsampling ", ncol(data), " -> ", MAX_CELLS, " cells..."))
  set.seed(RANDOM_SEED)
  data <- subset(data, cells = sample(colnames(data), MAX_CELLS))
}
if (ncol(data) < 200) {
  stop(paste0("Only ", ncol(data), " cells remain after subsetting. ",
              "Trajectory inference needs a reasonably dense manifold."), call. = FALSE)
}
# =============================================================================
# --- STEP 3: CONVERT SEURAT -> ANNDATA ---------------------------------------
# =============================================================================
# Built manually rather than via SeuratDisk, which is unmaintained and breaks
# regularly on Seurat v5 objects. A direct construction is more code but has no
# extra dependency and no version surprises.
#
# AnnData convention is CELLS x GENES, the transpose of Seurat.
message("\n=== STEP 3: Converting to AnnData ===")
counts <- GetAssayData(data, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
# Transpose to cells x genes, then hand the CSC triplet form to scipy.
message(paste0("  Building sparse matrix (", ncol(counts), " cells x ",
               nrow(counts), " genes)..."))
cx <- as(Matrix::t(counts), "CsparseMatrix")
X_py <- sp$csc_matrix(
  reticulate::tuple(
    np$array(cx@x),
    np$array(as.integer(cx@i)),
    np$array(as.integer(cx@p))
  ),
  shape = reticulate::tuple(as.integer(nrow(cx)), as.integer(ncol(cx)))
)
X_py <- X_py$tocsr()
# obs: cell metadata. Factors are converted to character so pandas does not
# invent its own categorical ordering.
md <- data@meta.data
obs_cols <- intersect(
  c(CELLTYPE_COLUMN, SAMPLE_COLUMN, CONDITION_COLUMN,
    "potency_score", "CytoTRACE2_Score", "CytoTRACE2_Potency",
    "CytoTRACE1_Score", "entropy_score", "gene_counts_score",
    "nFeature_RNA", "nCount_RNA", "percent.mt"),
  colnames(md)
)
obs_df <- md[, obs_cols, drop = FALSE]
for (cn in colnames(obs_df)) {
  if (is.factor(obs_df[[cn]])) obs_df[[cn]] <- as.character(obs_df[[cn]])
}
obs_df$barcode <- rownames(md)
pd     <- reticulate::import("pandas", convert = FALSE)
obs_py <- pd$DataFrame(reticulate::r_to_py(obs_df))
obs_py$index <- reticulate::r_to_py(as.character(rownames(md)))
var_py <- pd$DataFrame(reticulate::r_to_py(
  data.frame(gene_name = rownames(counts), stringsAsFactors = FALSE)
))
var_py$index <- reticulate::r_to_py(as.character(rownames(counts)))
adata <- ad$AnnData(X = X_py, obs = obs_py, var = var_py)
# --- Carry over the Harmony embedding ---------------------------------------
# Reusing Harmony means the neighbour graph is batch-corrected, so CellRank
# does not mistake sample-of-origin for a differentiation axis.
harmony_ok <- FALSE
if (USE_HARMONY_EMBEDDING && HARMONY_REDUCTION %in% names(data@reductions)) {
  emb <- Embeddings(data, reduction = HARMONY_REDUCTION)
  n_use <- min(N_PCS, ncol(emb))
  reticulate::py_set_item(adata$obsm, "X_harmony",
                          np$array(emb[, seq_len(n_use), drop = FALSE]))
  harmony_ok <- TRUE
  message(paste0("  Harmony embedding carried over (", n_use, " dims)."))
} else if (USE_HARMONY_EMBEDDING) {
  message(paste0("  [NOTE] Reduction '", HARMONY_REDUCTION,
                 "' not found; PCA will be recomputed in Python."))
}
# Carry over the existing UMAP so plots match the rest of the pipeline.
umap_name <- if ("umap_harmony" %in% names(data@reductions)) "umap_harmony" else
             if ("umap" %in% names(data@reductions)) "umap" else NA_character_
if (!is.na(umap_name)) {
  reticulate::py_set_item(adata$obsm, "X_umap",
                          np$array(Embeddings(data, reduction = umap_name)[, 1:2]))
  message(paste0("  UMAP embedding carried over from '", umap_name, "'."))
}
message(paste0("  AnnData built: ", reticulate::py_to_r(adata$n_obs), " cells x ",
               reticulate::py_to_r(adata$n_vars), " genes"))
rm(counts, cx, X_py); gc()
# =============================================================================
# --- STEP 4: PREPROCESSING IN SCANPY -----------------------------------------
# =============================================================================
message("\n=== STEP 4: Preprocessing (scanpy) ===")
# Keep a copy of raw counts: CytoTRACE's gene-count signal must be computed on
# counts, not on the log-normalised matrix.
reticulate::py_set_item(adata$layers, "counts", adata$X$copy())
sc$pp$filter_genes(adata, min_cells = 10L)
sc$pp$normalize_total(adata)
sc$pp$log1p(adata)
sc$pp$highly_variable_genes(adata, n_top_genes = as.integer(N_HVG))
if (harmony_ok) {
  # Use the precomputed batch-corrected space.
  sc$pp$neighbors(adata,
                  n_neighbors = as.integer(N_NEIGHBORS),
                  use_rep     = "X_harmony")
} else {
  sc$pp$pca(adata, n_comps = as.integer(N_PCS))
  sc$pp$neighbors(adata,
                  n_neighbors = as.integer(N_NEIGHBORS),
                  n_pcs       = as.integer(N_PCS))
}
message("  Neighbour graph built.")
# --- Imputation / smoothing --------------------------------------------------
# CytoTRACEKernel correlates each gene with the per-cell gene count. On raw
# sparse data those correlations are dominated by dropout, so a smoothed layer
# is required. The kernel then reads from that layer.
impute_layer <- NULL
if (CT_IMPUTE_METHOD == "moments") {
  message("  Building smoothed 'Ms' layer (kNN averaging over the neighbour graph)...")
  # scvelo's scv.pp.moments needs spliced/unspliced (velocity) layers, which plain
  # scRNA lacks -> it silently skips and never writes 'Ms'. The CytoTRACE 'first
  # moment' is just the kNN-smoothed expression, so build it directly from scanpy's
  # connectivity graph:  Ms = rownormalize(connectivities + I) @ X.
  .pymain <- reticulate::import_main()
  .pymain$adata <- adata               # expose adata to the __main__ namespace
  reticulate::py_run_string(paste(
    "import scipy.sparse as _sp, numpy as _np",
    "_conn = adata.obsp['connectivities']",
    "_A = _conn + _sp.eye(_conn.shape[0], format='csr')",
    "_d = _np.asarray(_A.sum(1)).ravel(); _d[_d == 0] = 1.0",
    "adata.layers['Ms'] = (_sp.diags(1.0 / _d) @ _A) @ adata.X",
    sep = "\n"))
  impute_layer <- "Ms"
} else if (CT_IMPUTE_METHOD == "cellmapper") {
  message("  Computing imputed layer (CellMapper)...")
  cm <- reticulate::import("cellmapper", convert = FALSE)
  mapper <- cm$CellMapper(adata)$map(layer_key = "X", use_rep = if (harmony_ok) "X_harmony" else "X_pca")
  reticulate::py_set_item(adata$layers, "imputed", mapper$query_imputed$X$copy())
  impute_layer <- "imputed"
} else {
  message("  [NOTE] CT_IMPUTE_METHOD = 'none'; using the log-normalised matrix.")
  impute_layer <- NULL
}
# =============================================================================
# --- STEP 5: CYTOTRACE KERNEL AND TRANSITION MATRIX --------------------------
# =============================================================================
message("\n=== STEP 5: CytoTRACEKernel ===")
ctk <- cr$kernels$CytoTRACEKernel(adata)
# compute_cytotrace() writes into the AnnData:
#   obs['ct_score']         — CytoTRACE-like score per cell
#   obs['ct_pseudotime']    — 1 - ct_score, i.e. a differentiation pseudotime
#   obs['ct_num_exp_genes'] — detected genes per cell
#   var['ct_gene_corr']     — per-gene correlation with the gene count
ctk <- if (is.null(impute_layer)) {
  ctk$compute_cytotrace()
} else {
  ctk$compute_cytotrace(layer = impute_layer)
}
message("  CytoTRACE score computed.")
ctk <- ctk$compute_transition_matrix()
message("  Transition matrix computed.")
# =============================================================================
# --- STEP 6: MACROSTATES, TERMINAL STATES, FATE PROBABILITIES ----------------
# =============================================================================
message("\n=== STEP 6: Estimating states and fates ===")
g <- cr$estimators$GPCCA(ctk)
# --- 6a: Schur decomposition + spectrum -------------------------------------
# Compute the Schur spectrum with a few extra components so the eigenvalue GAP is
# visible, and save the plot. Choose N_STATES = the number of eigenvalues BEFORE
# the largest gap (that is the data-driven macrostate count), then re-run. On a
# single lineage (e.g. Colonocytes only) this is usually small (~2-4), not 5.
n_schur <- as.integer(max(N_STATES + 3L, 10L))
tryCatch({
  g$compute_schur(n_components = n_schur)
  old_wd <- getwd(); setwd(TRAJ_DIR)
  sc$settings$figdir <- TRAJ_DIR
  sc$settings$set_figure_params(dpi_save = as.integer(DPI_SETTING), format = "png")
  g$plot_spectrum(real_only = TRUE, save = "schur_spectrum.png")
  setwd(old_wd)
  message("  Saved schur_spectrum.png -> set N_STATES to the #eigenvalues before the gap.")
}, error = function(e) message("  [WARNING] Schur/spectrum plot failed: ", e$message))
fate_ok <- tryCatch({
  g$compute_macrostates(n_states = as.integer(N_STATES),
                        cluster_key = CLUSTER_KEY)
  message(paste0("  ", N_STATES, " macrostates computed."))
  g$predict_terminal_states(n_states = as.integer(N_TERMINAL_STATES))
  message(paste0("  ", N_TERMINAL_STATES, " terminal states predicted."))
  g$compute_fate_probabilities()
  message("  Fate probabilities computed.")
  TRUE
}, error = function(e) {
  message(paste0("  [WARNING] State estimation failed: ", e$message))
  message("            Pseudotime is still available; fate probabilities are not.")
  message("            Try lowering N_STATES (see schur_spectrum.png), or a cleaner lineage.")
  FALSE
})
# =============================================================================
# --- STEP 7: PULL RESULTS BACK INTO SEURAT -----------------------------------
# =============================================================================
message("\n=== STEP 7: Importing results back into Seurat ===")
obs_back <- reticulate::py_to_r(adata$obs)
bc       <- rownames(obs_back)
add_col <- function(seu, values, name) {
  v <- values[match(colnames(seu), bc)]
  seu@meta.data[[name]] <- v
  seu
}
for (cn in intersect(c("ct_score", "ct_pseudotime", "ct_num_exp_genes"),
                     colnames(obs_back))) {
  data <- add_col(data, obs_back[[cn]], paste0("cellrank_", cn))
  message(paste0("  Added: cellrank_", cn))
}
fate_df <- NULL
if (fate_ok) {
  # NOTE: the tryCatch main expression evaluates in THIS frame, so the
  # `data <- add_col(...)` assignments below correctly update the object.
  # (Only an error = function(e) handler would need `<<-`.)
  invisible(tryCatch({
    fp     <- adata$obsm$get("lineages_fwd")
    fp_r   <- reticulate::py_to_r(np$array(fp))
    lin_names <- tryCatch(
      as.character(reticulate::py_to_r(g$terminal_states$cat$categories$tolist())),
      error = function(e) paste0("Fate_", seq_len(ncol(fp_r)))
    )
    if (length(lin_names) != ncol(fp_r)) {
      lin_names <- paste0("Fate_", seq_len(ncol(fp_r)))
    }
    colnames(fp_r) <- paste0("fate_", make.names(lin_names))
    rownames(fp_r) <- bc
    for (cn in colnames(fp_r)) {
      data <- add_col(data, fp_r[, cn], cn)
    }
    # Dominant fate + its confidence: usually more interpretable than the
    # full probability vector when plotting.
    dom_idx <- max.col(fp_r, ties.method = "first")
    data <- add_col(data, colnames(fp_r)[dom_idx], "fate_dominant")
    data <- add_col(data, apply(fp_r, 1, max),     "fate_confidence")
    fate_df <- as.data.frame(fp_r)
    message(paste0("  Added ", ncol(fp_r), " fate probability columns."))
    TRUE
  }, error = function(e) {
    message(paste("  [WARNING] Could not import fate probabilities:", e$message))
    FALSE
  }))
}
# =============================================================================
# --- STEP 7b: TERMINAL-STATE POTENCY SANITY CHECK ----------------------------
# =============================================================================
# A terminal state should be the MOST differentiated population, so its cells
# should carry LOW potency (low CytoTRACE2 / low CellRank ct_score). If a terminal
# state's median potency is ABOVE the global median, it looks stem-like and is
# probably a wrong call - lower N_STATES / N_TERMINAL_STATES or revisit the kernel.
if (exists("fate_ok") && isTRUE(fate_ok)) {
  tryCatch({
    ts <- reticulate::py_to_r(g$terminal_states)          # factor over cells, NA = not terminal
    ts <- stats::setNames(as.character(ts), names(ts))
    pot_col <- intersect(c("CytoTRACE2_Score", "potency_score", "cellrank_ct_score"),
                         colnames(data@meta.data))
    pot_col <- if (length(pot_col) > 0) pot_col[1] else NA_character_
    if (!is.na(pot_col) && any(!is.na(ts))) {
      pv <- data@meta.data[[pot_col]]; names(pv) <- colnames(data)
      glob_med <- stats::median(pv, na.rm = TRUE)
      message(sprintf("  Terminal-state potency check (%s; global median = %.3f):",
                      pot_col, glob_med))
      for (st in unique(ts[!is.na(ts)])) {
        cells_st <- names(ts)[which(ts == st)]
        m <- stats::median(pv[cells_st], na.rm = TRUE)
        flag <- if (isTRUE(m > glob_med))
          "  <-- ABOVE average: looks stem-like, suspicious for a terminal state" else ""
        message(sprintf("    '%s': median %s = %.3f  (n = %d)%s",
                        st, pot_col, m, length(cells_st), flag))
      }
      message("  Expectation: terminal states sit BELOW the global median (most differentiated).")
    } else {
      message("  [NOTE] No potency column (CytoTRACE2_Score/potency_score/cellrank_ct_score) ",
              "on the object - skipping the terminal-state potency check.")
    }
  }, error = function(e) message("  [WARNING] Terminal-state potency check failed: ", e$message))
}
# =============================================================================
# --- STEP 8: PLOTS AND TABLES ------------------------------------------------
# =============================================================================
message("\n=== STEP 8: Plots and tables ===")
reduction_use <- if ("umap_harmony" %in% names(data@reductions)) "umap_harmony" else
                 if ("umap" %in% names(data@reductions)) "umap" else NA_character_
traj_cols <- grep("^cellrank_|^fate_", colnames(data@meta.data), value = TRUE)
num_cols  <- traj_cols[vapply(traj_cols,
                              function(c) is.numeric(data@meta.data[[c]]),
                              logical(1))]
# --- Cell-type colour palette (shared R <-> Python) --------------------------
# Seurat's DimPlot colours identities with scales::hue_pal() over the sorted
# levels. Reproduce that exact palette, use it in the R DimPlots, AND push it into
# AnnData (categorical order + uns colours) so the scanpy/CellRank plots draw the
# SAME colours per cell type.
ct_cols <- NULL
if (CELLTYPE_COLUMN %in% colnames(data@meta.data)) {
  ct_levels <- sort(unique(as.character(data@meta.data[[CELLTYPE_COLUMN]])))
  ct_cols   <- scales::hue_pal()(length(ct_levels))
  names(ct_cols) <- ct_levels
  tryCatch({
    .pm <- reticulate::import_main()
    .pm$adata        <- adata
    .pm$ct_levels_r  <- as.character(ct_levels)
    .pm$ct_cols_r    <- as.character(unname(ct_cols))
    .pm$celltype_key <- CELLTYPE_COLUMN
    reticulate::py_run_string(paste(
      "import pandas as _pd",
      "adata.obs[celltype_key] = _pd.Categorical(adata.obs[celltype_key].astype(str), categories=list(ct_levels_r))",
      "adata.uns[celltype_key + '_colors'] = list(ct_cols_r)",
      sep = "\n"))
    message("  Cell-type palette pushed to AnnData (matches R DimPlot colours).")
  }, error = function(e)
    message("  [WARNING] Could not push cell-type colours to AnnData: ", e$message))
}
if (!is.na(reduction_use)) {
  for (cn in num_cols) {
    tryCatch({
      p <- FeaturePlot(data, features = cn, reduction = reduction_use, pt.size = 0.3) +
        scale_color_viridis_c() + coord_fixed() +
        labs(title = cn) + theme(plot.title = element_text(face = "bold"))
      ggsave(file.path(TRAJ_DIR, paste0("umap_", cn, ".png")),
             p, width = 8, height = 6.5, dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) {
      message(paste("  [WARNING] Plot failed for", cn, ":", e$message))
    })
  }
  if ("fate_dominant" %in% colnames(data@meta.data)) {
    tryCatch({
      p <- DimPlot(data, group.by = "fate_dominant", reduction = reduction_use,
                   pt.size = 0.3) +
        coord_fixed() + labs(title = "Dominant predicted fate") +
        theme(plot.title = element_text(face = "bold"))
      ggsave(file.path(TRAJ_DIR, "umap_fate_dominant.png"),
             p, width = 8, height = 6.5, dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) message(paste("  [WARNING]", e$message)))
  }
  # Cell-type reference on the SAME UMAP, so the flow/pseudotime maps can be read
  # against the annotation.
  if (CELLTYPE_COLUMN %in% colnames(data@meta.data)) {
    tryCatch({
      p <- DimPlot(data, group.by = CELLTYPE_COLUMN, reduction = reduction_use,
                   pt.size = 0.3, label = TRUE, repel = TRUE, cols = ct_cols) +
        coord_fixed() + labs(title = paste0("Cell type (", CELLTYPE_COLUMN, ")")) +
        theme(plot.title = element_text(face = "bold"))
      ggsave(file.path(TRAJ_DIR, "umap_celltype.png"),
             p, width = 9, height = 6.5, dpi = DPI_SETTING, bg = "white")
      rm(p)
    }, error = function(e) message(paste("  [WARNING] Cell-type UMAP failed:", e$message)))
  }
}
# --- Pseudotime by cell type -------------------------------------------------
if ("cellrank_ct_pseudotime" %in% colnames(data@meta.data)) {
  tryCatch({
    d <- data@meta.data[!is.na(data@meta.data$cellrank_ct_pseudotime), , drop = FALSE]
    ord <- d %>% dplyr::group_by(.data[[CELLTYPE_COLUMN]]) %>%
      dplyr::summarise(m = stats::median(cellrank_ct_pseudotime), .groups = "drop") %>%
      dplyr::arrange(m)
    d[[CELLTYPE_COLUMN]] <- factor(as.character(d[[CELLTYPE_COLUMN]]),
                                   levels = as.character(ord[[CELLTYPE_COLUMN]]))
    p <- ggplot(d, aes(x = .data[[CELLTYPE_COLUMN]], y = cellrank_ct_pseudotime,
                       fill = .data[[CELLTYPE_COLUMN]])) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.6, linewidth = 0.3) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, linewidth = 0.3) +
      labs(title = "CytoTRACE pseudotime by cell type",
           subtitle = "low = less differentiated (earlier)",
           x = NULL, y = "ct_pseudotime") +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold"))
    ggsave(file.path(TRAJ_DIR, "pseudotime_by_celltype.png"),
           p, width = 11, height = 6, dpi = DPI_SETTING, bg = "white")
    rm(p)
  }, error = function(e) message(paste("  [WARNING] Pseudotime plot failed:", e$message)))
}
# --- CellRank's own plots ----------------------------------------------------
# Rendered through scanpy's figure machinery, which writes into a 'figures/'
# subdirectory of the working directory.
tryCatch({
  old_wd <- getwd()
  setwd(TRAJ_DIR)
  sc$settings$figdir <- TRAJ_DIR
  sc$settings$set_figure_params(dpi_save = as.integer(DPI_SETTING), format = "png")
  if (fate_ok) {
    g$plot_macrostates(which = "all", save = "macrostates.png")
    g$plot_macrostates(which = "terminal", save = "terminal_states.png")
    g$plot_fate_probabilities(save = "fate_probabilities.png")
  }
  setwd(old_wd)
  message("  CellRank plots saved.")
}, error = function(e) {
  message(paste("  [WARNING] CellRank plotting failed:", e$message))
})
# --- Projected differentiation flow (streamlines + arrows) -------------------
# Projects the CytoTRACEKernel transition matrix onto the UMAP to draw the
# expected displacement field (velocity-style arrows), and colours it by cell
# type so the flow can be read against the annotation. Arrows point toward
# increasing differentiation (i.e. toward the terminal states).
tryCatch({
  old_wd <- getwd()
  setwd(TRAJ_DIR)
  sc$settings$figdir <- TRAJ_DIR
  sc$settings$set_figure_params(dpi_save = as.integer(DPI_SETTING), format = "png")
  # alpha = 1 keeps the cell scatter fully opaque (scvelo defaults to a pale ~0.3);
  # colours come from adata.uns[CellType_colors] set above, matching the R DimPlot.
  ctk$plot_projection(basis = "umap", color = CLUSTER_KEY, recompute = TRUE,
                      stream = TRUE,  alpha = 1, save = "flow_streamlines.png")
  ctk$plot_projection(basis = "umap", color = CLUSTER_KEY, recompute = FALSE,
                      stream = FALSE, alpha = 1, save = "flow_arrows.png")
  setwd(old_wd)
  message("  Saved flow_streamlines.png + flow_arrows.png (coloured by ", CLUSTER_KEY, ").")
}, error = function(e) {
  message(paste("  [WARNING] Flow projection plot failed:", e$message))
})
# --- Export UMAP coordinates (+ scalars) for native R re-plotting ------------
# Same embedding the Python plots use, so anything drawn from this CSV in R lines
# up cell-for-cell with the scanpy/CellRank figures.
if (!is.na(reduction_use)) {
  tryCatch({
    um <- Embeddings(data, reduction = reduction_use)[, 1:2, drop = FALSE]
    ex <- data.frame(Barcode = rownames(um), UMAP_1 = um[, 1], UMAP_2 = um[, 2],
                     stringsAsFactors = FALSE)
    for (cn in intersect(c(CELLTYPE_COLUMN, CONDITION_COLUMN, SAMPLE_COLUMN,
                           "cellrank_ct_pseudotime", "cellrank_ct_score",
                           "fate_dominant", "fate_confidence"),
                         colnames(data@meta.data)))
      ex[[cn]] <- data@meta.data[[cn]]
    utils::write.csv(ex, file.path(TRAJ_DIR, "umap_coordinates.csv"), row.names = FALSE)
    message("  Exported umap_coordinates.csv (", reduction_use, ").")
  }, error = function(e) message("  [WARNING] UMAP coordinate export failed: ", e$message))
}
# --- Summary tables ----------------------------------------------------------
sheets <- list()
if (length(num_cols) > 0) {
  sheets[["By_CellType"]] <- data@meta.data %>%
    dplyr::group_by(.data[[CELLTYPE_COLUMN]]) %>%
    dplyr::summarise(N_Cells = dplyr::n(),
              dplyr::across(dplyr::all_of(num_cols),
                            ~stats::median(.x, na.rm = TRUE),
                            .names = "{.col}_median"),
              .groups = "drop") %>%
    as.data.frame()
  sheets[["By_Condition"]] <- data@meta.data %>%
    dplyr::group_by(.data[[CELLTYPE_COLUMN]], .data[[CONDITION_COLUMN]]) %>%
    dplyr::summarise(N_Cells = dplyr::n(),
              dplyr::across(dplyr::all_of(num_cols),
                            ~stats::median(.x, na.rm = TRUE),
                            .names = "{.col}_median"),
              .groups = "drop") %>%
    as.data.frame()
}
sheets[["Run_Info"]] <- data.frame(
  Parameter = c("PROJECT_NAME", "N_Cells", "SUBSET_TO", "SUBSET_CONDITION",
                "CT_IMPUTE_METHOD", "USE_HARMONY_EMBEDDING", "N_HVG", "N_PCS",
                "N_NEIGHBORS", "N_STATES", "N_TERMINAL_STATES",
                "Fate_Probabilities_Computed", "Python_Env", "Date"),
  Value = c(PROJECT_NAME, ncol(data),
            if (is.null(SUBSET_TO)) "ALL" else paste(SUBSET_TO, collapse = "; "),
            if (is.null(SUBSET_CONDITION)) "ALL" else paste(SUBSET_CONDITION, collapse = "; "),
            CT_IMPUTE_METHOD, USE_HARMONY_EMBEDDING, N_HVG, N_PCS,
            N_NEIGHBORS, N_STATES, N_TERMINAL_STATES,
            fate_ok, PYTHON_ENV, as.character(Sys.Date())),
  stringsAsFactors = FALSE
)
write_xlsx(sheets, file.path(TRAJ_DIR, "trajectory_summary.xlsx"))
message("  Summary written to: trajectory_summary.xlsx")
# =============================================================================
# --- STEP 9: SAVE ------------------------------------------------------------
# =============================================================================
if (SAVE_H5AD) {
  tryCatch({
    h5_path <- file.path(TRAJ_DIR, "adata_cellrank.h5ad")
    # Categorical/enum columns can break h5ad writing; coerce to plain strings.
    adata$obs <- adata$obs$astype(reticulate::dict(list()))
    adata$write_h5ad(h5_path)
    message(paste("  AnnData written to:", basename(h5_path)))
  }, error = function(e) {
    message(paste("  [WARNING] Could not write h5ad:", e$message))
    message("            (Known CellRank issue: enum columns in .obs/.uns.)")
  })
}
out_rds <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_trajectory.rds"))
saveRDS(data, out_rds)
message("\n=== Script 10 complete ===")
message(paste0("  Trajectory columns added: ", paste(traj_cols, collapse = ", ")))
message(paste0("  Plots and tables: ", TRAJ_DIR))
message(paste0("  Object: ", out_rds))
message("\n  INTERPRETATION REMINDERS:")
message("   - ct_pseudotime is 1 - ct_score: LOW = less differentiated.")
message("   - Terminal states are statistical endpoints, not proven biology.")
message("     Cross-check them against known markers before drawing conclusions.")
message("   - If SUBSET_TO was NULL, transitions between unrelated lineages are")
message("     present in the model and the result should not be trusted.")
