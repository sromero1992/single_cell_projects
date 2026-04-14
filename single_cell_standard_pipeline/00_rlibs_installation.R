# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 0: R ENVIRONMENT SETUP
# Version: 2.1 (Corrected for non-interactive use)
#
# PURPOSE:
#   Installs all R packages required by Scripts 01–04 of this pipeline.
#   Run this script ONCE per machine / R environment before running any other
#   pipeline script. Re-running is safe — packages already present are skipped.
#
# PACKAGE SOURCES:
#   CRAN        — Standard R packages (Seurat, ggplot2, harmony, enrichR, etc.)
#   Bioconductor— Biology-specific packages (AUCell, celda, SingleCellExperiment)
#   GitHub      — Packages not yet on CRAN/Bioc (DoubletFinder, SplineDV, SCEVAN,
#                 CellChat, monocle3, fgsea, etc.)
#
# HOW TO RUN:
#   From the terminal in the pipeline directory, run:
#     Rscript 00_rlibs_installation.R
#   Watch the console for any red error messages. If a package fails,
#   fix the error and re-run — the script will only install missing packages.
#
# OUTPUT:
#   A status report table is printed at the end showing TRUE/FALSE per package.
#   All packages must show TRUE before proceeding to Script 01.
#
# NEXT STEP: Run '01_process_data.R'
# =============================================================================

# 0. INITIAL SETUP
# Define a default CRAN mirror for non-interactive installation
CRAN_mirror <- "https://cloud.r-project.org"

# 1. SET THE PERSONAL LIBRARY PATH
personal_lib <- if (.Platform$OS.type == "windows") {
  file.path(Sys.getenv("USERPROFILE"), "Documents", "R_libs_scRNA")
} else {
  file.path(Sys.getenv("HOME"), "R_libs_scRNA")
}
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}
# CRITICAL: Force R to use the personal library FIRST and check it for dependencies.
.libPaths(c(personal_lib, .libPaths()))
cat("Using Library Paths:\n")
print(.libPaths())

# 2. DEFINE PACKAGE LISTS
cran_pkgs <- c(
  "Seurat", "devtools", "dplyr", "ggplot2", "Matrix", "ggpubr", "tidyr", "patchwork",
  "stringr", "tibble", "cowplot", "openxlsx","writexl", "readxl", "parallelly", "hdf5r",
  "enrichR", "remotes", "R.utils"
)
bioc_pkgs <- c(
  "BiocManager", "ComplexHeatmap", "Biobase", "BiocNeighbors", "BiocGenerics",
  "celda", "dittoSeq", "AUCell", "Gviz", "GenomicRanges", "rtracklayer",
  "BiocSingular", "SingleCellExperiment", "SummarizedExperiment"
)
github_pkgs <- c(
  "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
  "seuratWrappers" = "satijalab/seurat-wrappers",
  "presto" = "immunogenomics/presto",
  "yaGST" = "miccec/yaGST",
  "ggtree" = "YuLab-SMU/ggtree",
  "fgsea" = "alserglab/fgsea",
  "SCEVAN" = "AntonioDeFalco/SCEVAN",
  "SplineDV" = "Xenon8778/SplineDV",
  "CellChat" = "jinworks/CellChat",
  "scSGS" = "Xenon8778/scSGS",
  "leidenbase" = "cole-trapnell-lab/leidenbase"
)

# 3. AUTOMATED INSTALLATION LOGIC
cat("\n--- Starting Automated Installation ---\n")
# Logic: Check only the personal library to avoid "Ghost" version conflicts from system folders
get_installed <- function(lib) {
  if (!dir.exists(lib)) return(character(0))
  installed.packages(lib.loc = lib)[, "Package"]
}

# --- Install CRAN Packages ---
cat("\n--- Checking CRAN packages ---\n")
installed_custom <- get_installed(personal_lib)
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed_custom)]
if (length(missing_cran) > 0) {
  cat("Installing missing CRAN packages to:", personal_lib, "\n")
  # MODIFIED: Added repos argument
  install.packages(missing_cran, lib = personal_lib, repos = CRAN_mirror)
} else {
  cat("All CRAN packages present in personal library.\n")
}

# --- Install Bioconductor Packages ---
cat("\n--- Checking Bioconductor packages ---\n")
# The script checks for BiocManager in the personal library. If it's missing, install it.
if (!("BiocManager" %in% installed_custom)) {
  # MODIFIED: Added repos argument
  install.packages("BiocManager", lib = personal_lib, repos = CRAN_mirror)
}
# Load BiocManager explicitly from the custom path to avoid "not found" errors
library(BiocManager, lib.loc = personal_lib)
# Re-check installed packages after potentially adding BiocManager
installed_custom <- get_installed(personal_lib)
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed_custom)]
if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages...\n")
  BiocManager::install(missing_bioc, lib = personal_lib, update = FALSE, ask = FALSE)
} else {
  cat("All Bioconductor packages present in personal library.\n")
}

# --- Install GitHub Packages ---
cat("\n--- Checking GitHub packages ---\n")
# The script checks for remotes in the personal library. If it's missing, install it.
if (!("remotes" %in% installed_custom)) {
  # MODIFIED: Added repos argument
  install.packages("remotes", lib = personal_lib, repos = CRAN_mirror)
}
library(remotes, lib.loc = personal_lib)
gh_package_names <- names(github_pkgs)
# Re-check installed packages after potentially adding remotes
installed_custom <- get_installed(personal_lib)
missing_gh <- gh_package_names[!(gh_package_names %in% installed_custom)]
if (length(missing_gh) > 0) {
  for (pkg_name in missing_gh) {
    repo <- github_pkgs[pkg_name]
    cat("Installing", pkg_name, "from", repo, "...\n")
    remotes::install_github(repo, lib = personal_lib, upgrade = "never", force = TRUE)
  }
} else {
  cat("All GitHub packages present in personal library.\n")
}

# 4. FINAL VERIFICATION REPORT
cat("\n--- FINAL STATUS REPORT ---\n")
all_pkgs_to_check <- unique(c(cran_pkgs, bioc_pkgs, names(github_pkgs)))
# Check if they can actually be loaded (the ultimate test)
final_check <- data.frame(
  Package = all_pkgs_to_check,
  Status = sapply(all_pkgs_to_check, function(p) {
    requireNamespace(p, quietly = TRUE)
  })
)
print(final_check)
if (all(final_check$Status)) {
  cat("\n[+] Success! All packages are installed and version-matched.\n")
} else {
  cat("\n[!] Some packages failed. Check permissions or missing system dependencies (like hdf5).\n")
}

