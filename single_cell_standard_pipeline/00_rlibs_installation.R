# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 0: R ENVIRONMENT SETUP (FORCE MODE)
# Version: 4.0 - UNIFIED (Red Hat 9 Optimized)
#
# CHANGES IN THIS UNIFIED VERSION:
#   - BOTH doublet callers installed: DoubletFinder (GitHub) and scDblFinder
#     (Bioconductor), since Script 01 can now run either or both.
#   - Added BPCells: Script 01 calls library(BPCells) when USE_BPCELLS = TRUE,
#     but it was never in the install list in previous versions.
#   - Fixed 'seuratWrappers' -> 'SeuratWrappers'. The namespace is capitalised,
#     so the old lowercase spelling never matched an installed package and the
#     package was force-reinstalled on every single run.
#   - Added a hard verification gate at the end: the script now fails loudly if
#     a package Script 01 depends on is missing, instead of printing a table
#     that is easy to skim past.
# =============================================================================

setwd("/mnt/SCDC/Optimus/selim_working_dir/2026_nr4a1_ack/r_process/debug_pipeline_pkg")
# 0. THE "FORCE" PROTOCOL - Locking the environment to official system paths
#Sys.setenv(LD_LIBRARY_PATH = "/usr/lib64")
#Sys.setenv(TMPDIR = paste0(Sys.getenv("HOME"), "/Rtemp"))
options(timeout = 1200) # Prevents timeout on large GitHub builds like Monocle3

if (!dir.exists(Sys.getenv("TMPDIR"))) dir.create(Sys.getenv("TMPDIR"))

# 1. DEFINE PACKAGE LISTS
cran_pkgs <- c(
  "Seurat", "devtools", "dplyr", "ggplot2", "Matrix", "ggpubr", "tidyr", "patchwork",
  "stringr", "tibble", "cowplot", "openxlsx","writexl", "readxl", "parallelly", "hdf5r",
  "enrichR", "remotes", "R.utils", "sf", "harmony", "xgboost"
)

bioc_pkgs <- c(
  "BiocManager", "ComplexHeatmap", "Biobase", "BiocNeighbors", "BiocGenerics",
  "celda", "dittoSeq", "AUCell", "Gviz", "GenomicRanges", "rtracklayer", "MAST",
  "BiocSingular", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder",
  # Symbol <-> Ensembl mapping for the PathVisio/WikiPathways exports in Script 07
  # and for the Ensembl->symbol step in Script 09's CytoTRACE 2 preprocessing.
  # org.Hs.eg.db is needed only when CT2_SPECIES = "human"; cheap to keep here.
  "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db"
)

github_pkgs <- c(
  # DoubletFinder: DEFAULT doublet caller in Script 01 (paramSweep + pK).
  "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
  # CytoTRACE2: absolute developmental potency, Script 09.
  # NOTE the subdir - handled separately below, since install_github() needs
  # the `subdir` argument and this named-vector loop cannot express it.
  # CytoTRACE (v1): optional cross-check in Script 09. Old package with stale
  # dependencies; see INSTALL_NOTES.md section 2 if it fails.
  "CytoTRACE" = "gunsagargulati/CytoTRACE",
  # BPCells: on-disk matrix backend, required when USE_BPCELLS = TRUE.
  "BPCells" = "bnprks/BPCells/r",
  # NOTE: namespace is capital-S 'SeuratWrappers'. The previous lowercase
  # spelling never matched, causing a forced reinstall on every run.
  "SeuratWrappers" = "satijalab/seurat-wrappers",
  "presto" = "immunogenomics/presto",
  "yaGST" = "miccec/yaGST",
  "ggtree" = "YuLab-SMU/ggtree",
  "fgsea" = "alserglab/fgsea",
  "SCEVAN" = "AntonioDeFalco/SCEVAN",
  "SplineDV" = "Xenon8778/SplineDV",
  "CellChat" = "jinworks/CellChat",
  "scSGS" = "Xenon8778/scSGS",
  "leidenbase" = "cole-trapnell-lab/leidenbase",
  "monocle3" = "cole-trapnell-lab/monocle3"
  )

# 2. AUTOMATED INSTALLATION LOGIC
is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# --- Install CRAN Packages with Manual Overrides ---
cat("\n--- Installing CRAN packages with Force Overrides ---\n")
for (pkg in cran_pkgs) {
  if (!is_installed(pkg)) {
    cat("Force Installing:", pkg, "\n")
    
    # Specific overrides for the "problematic" packages
    if (pkg == "sf") {
      install.packages(pkg, type = "source", configure.args = "--with-gdal-config=/usr/bin/gdal-config --with-proj-lib=/usr/lib64")
    } else if (pkg == "hdf5r") {
      install.packages(pkg, type = "source", configure.args = "--with-hdf5=/usr/bin/h5cc")
    } else {
      install.packages(pkg, type = "source")
    }
  }
}

# --- Install Bioconductor Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com")
}

cat("\n--- Checking Bioconductor packages ---\n")
missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, is_installed)]

if (length(missing_bioc) > 0) {
  # 1. Use the 'Ncpus' option to speed it up
  # 2. Use 'checkBuilt = TRUE' to ensure they match your RHEL 9 build
  # 3. Use 'INSTALL_opts = "--no-test-load"' to prevent R from loading them mid-install
  BiocManager::install(missing_bioc, 
                       update = FALSE, 
                       ask = FALSE, 
                       checkBuilt = TRUE,
                       INSTALL_opts = "--no-test-load")
}


# --- Install GitHub Packages ---
cat("\n--- Checking GitHub packages ---\n")
for (pkg_name in names(github_pkgs)) {
  if (!is_installed(pkg_name)) {
    repo <- github_pkgs[pkg_name]
    cat("Installing", pkg_name, "...\n")
    remotes::install_github(repo, upgrade = "never", force = TRUE, quiet = FALSE)
  }
}

# --- Install CytoTRACE 2 (special case: lives in a repo subdirectory) -------
# devtools::install_github() needs subdir = "cytotrace2_r". Omitting it gives a
# confusing "does not appear to be a package" error. Handled outside the loop
# above because the named-vector form cannot carry extra arguments.
cat("\n--- Checking CytoTRACE 2 (Script 09) ---\n")
if (!is_installed("CytoTRACE2")) {
  cat("Installing CytoTRACE2 from digitalcytometry/cytotrace2 ...\n")
  tryCatch({
    remotes::install_github("digitalcytometry/cytotrace2",
                            subdir  = "cytotrace2_r",
                            upgrade = "never")
  }, error = function(e) {
    cat("[WARNING] CytoTRACE2 install failed:", conditionMessage(e), "\n")
    cat("          Script 09 will skip it and still run entropy metrics.\n")
    cat("          See INSTALL_NOTES.md section 1 for fixes.\n")
  })
} else {
  cat("CytoTRACE2 already installed.\n")
}

# --- reticulate: only needed for Script 10 (CellRank trajectories) ----------
# The Python environment itself is NOT created here. Building conda
# environments from inside an analysis script is a good way to corrupt a
# working setup. See INSTALL_NOTES.md section 5 for the one-time setup.
if (!is_installed("reticulate")) {
  cat("\nInstalling reticulate (needed only by Script 10)...\n")
  install.packages("reticulate")
}

# 3. FINAL VERIFICATION REPORT
cat("\n--- FINAL STATUS REPORT ---\n")
all_pkgs_to_check <- unique(c(cran_pkgs, bioc_pkgs, names(github_pkgs),
                              "CytoTRACE2", "reticulate"))
final_check <- data.frame(
  Package = all_pkgs_to_check,
  Status  = sapply(all_pkgs_to_check, is_installed)
)
print(final_check)

# 4. HARD GATE ON SCRIPT 01 DEPENDENCIES
# Script 01 will fail at library() time if any of these are missing. Better to
# find out here, in 2 seconds, than 40 minutes into a processing run.
cat("\n--- VERIFYING SCRIPT 01 CRITICAL DEPENDENCIES ---\n")
critical_pkgs <- c(
  "Seurat", "SeuratWrappers", "openxlsx", "dplyr", "ggplot2", "patchwork",
  "celda",          # DecontX ambient RNA correction
  "DoubletFinder",  # default doublet caller
  "scDblFinder",    # alternative doublet caller
  "writexl", "Matrix", "SCEVAN", "hdf5r"
)
missing_critical <- critical_pkgs[!sapply(critical_pkgs, is_installed)]

if (length(missing_critical) > 0) {
  cat("\n!!! MISSING CRITICAL PACKAGES !!!\n")
  cat(paste0("  - ", missing_critical, collapse = "\n"), "\n")
  cat("\nScript 01 CANNOT run until these are installed.\n")
  cat("Re-run this script, or install the listed packages manually.\n")
  cat("If DoubletFinder fails to build, confirm 'remotes' and a working\n")
  cat("compiler toolchain are available, then retry install_github().\n\n")
  stop("Environment setup incomplete - see the list above.", call. = FALSE)
} else {
  cat("All critical dependencies present. Environment is ready for Script 01.\n")
}

# 5. OPTIONAL PACKAGES - report status, but never block the run
# These gate specific features. Missing ones only disable that feature, so
# they are reported as notes rather than errors.
cat("\n--- OPTIONAL PACKAGES (feature-gating, not blocking) ---\n")
optional_notes <- list(
  BPCells    = "Only needed if USE_BPCELLS <- TRUE in Script 01. Leave FALSE to proceed.",
  CytoTRACE2 = "Needed for CytoTRACE 2 potency in Script 09. See INSTALL_NOTES.md section 1.",
  CytoTRACE  = "Optional cross-check in Script 09. Set RUN_CYTOTRACE1 <- FALSE to skip.",
  reticulate = "Needed only for Script 10 (CellRank trajectories).",
  CellChat   = "Needed for Script 08 cell-cell communication.",
  MAST       = "Needed for Script 07 differential expression."
)
for (p in names(optional_notes)) {
  if (is_installed(p)) {
    cat(sprintf("  [OK]      %-12s\n", p))
  } else {
    cat(sprintf("  [MISSING] %-12s %s\n", p, optional_notes[[p]]))
  }
}

# Script 10 additionally needs a Python environment, which this script does not
# and should not build. Point the user at the instructions instead.
if (is_installed("reticulate")) {
  cat("\n[NOTE] Script 10 also requires a Python env with scanpy + cellrank.\n")
  cat("       reticulate alone is not sufficient. One-time setup:\n")
  cat("         conda create -n cellrank_env python=3.11 -y\n")
  cat("         conda activate cellrank_env\n")
  cat("         pip install \"cellrank>=2.0\" scanpy scvelo anndata igraph leidenalg\n")
  cat("       Full details and troubleshooting: INSTALL_NOTES.md section 5.\n")
}

# =============================================================================
# 6. LOCAL PACKAGE: TamuScDSC (preprocessing engine for Scripts 01 / 01b / 02)
# =============================================================================
# TamuScDSC is NOT on CRAN/Bioc/GitHub - it ships in this repo as a folder. Install
# it here, LAST, so its Imports (Seurat, dplyr, ...) are already present. It is a
# fast source build. If the folder cannot be found, this only warns - install it
# by hand with devtools::install("<path>/TamuScDSC").
cat("\n--- Installing local package: TamuScDSC ---\n")
# Set TAMUSCDSC_PATH explicitly to skip the search (e.g. "/path/to/repo/TamuScDSC").
if (!exists("TAMUSCDSC_PATH")) TAMUSCDSC_PATH <- NULL
TamuScDSC_candidates <- c(
  TAMUSCDSC_PATH,
  "TamuScDSC",            # cwd is the repo root
  "./TamuScDSC",         # cwd is unified_pipeline/, package one level up
  file.path(Sys.getenv("HOME"), "TamuScDSC")
)
TamuScDSC_dir <- Filter(function(p) !is.null(p) &&
                       file.exists(file.path(p, "DESCRIPTION")), TamuScDSC_candidates)

if (length(TamuScDSC_dir) > 0) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  tryCatch({
    devtools::install(TamuScDSC_dir[[1]], upgrade = "never", quiet = FALSE)
    cat("TamuScDSC installed from: ", TamuScDSC_dir[[1]], "\n")
  }, error = function(e) {
    cat("[WARNING] TamuScDSC install failed:", conditionMessage(e), "\n")
    cat("          Install by hand: devtools::install('", TamuScDSC_dir[[1]], "')\n", sep = "")
  })
} else {
  cat("[NOTE] TamuScDSC/ folder not found from the current working directory.\n")
  cat("       Set TAMUSCDSC_PATH <- \"/path/to/repo/TamuScDSC\" and re-run, or install\n")
  cat("       by hand: devtools::install(\"/path/to/repo/TamuScDSC\").\n")
}

cat("\nSetup complete.\n")

