# --- R Environment Setup for Single-Cell Analysis (v2) ---

# 0. INITIAL SETUP
# Install `remotes` if not already present, as it's needed for GitHub installations.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 1. SET THE PERSONAL LIBRARY PATH
# This is a best practice to keep project dependencies isolated.
personal_lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R_libs_scRNA")
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}
.libPaths(c(personal_lib, .libPaths())) # Prepend our new path to make it the default

# 2. DEFINE PACKAGE LISTS
# --- CRAN Packages ---
# CORRECTED: 'harmony' is now on CRAN and has been moved here.
cran_pkgs <- c(
  "Seurat", "devtools", "dplyr", "ggplot2", "Matrix", "ggpubr", "tidyr",
  "stringr", "tibble", "cowplot", "writexl", "readxl", "parallelly", "hdf5r",
  "harmony" # Added here from CRAN
)

# --- Bioconductor Packages ---
bioc_pkgs <- c(
  "BiocManager", "ComplexHeatmap", "Biobase", "BiocNeighbors", "BiocGenerics",
  "celda", "dittoSeq", "AUCell", "Gviz", "GenomicRanges", "rtracklayer",
  "BiocSingular", "SingleCellExperiment", "SummarizedExperiment"
)

# --- GitHub Packages ---
# Using a NAMED vector is a clean way to map package name to repo path.
# Format: "package_name" = "github_user/repo"
# CORRECTED: 'harmony' has been removed from this list.
github_pkgs <- c(
  "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
  "presto" = "immunogenomics/presto",
  "yaGST" = "miccec/yaGST",
  "ggtree" = "YuLab-SMU/ggtree",
  "fgsea" = "alserglab/fgsea",
  "SCEVAN" = "AntonioDeFalco/SCEVAN",
  "SplineDV" = "Xenon8778/SplineDV",
  "CellChat" = "jinworks/CellChat",
  "scSGS" = "Xenon8778/scSGS",
  "leidenbase" = "cole-trapnell-lab/leidenbase", # Monocle3 dependency
  "monocle3" = "cole-trapnell-lab/monocle3",
  "cicero" = "cole-trapnell-lab/cicero-release" # Installs as 'cicero'
)

# 3. AUTOMATED INSTALLATION LOGIC
cat("--- Starting Automated Installation ---\n")
cat("Packages will be installed in:", personal_lib, "\n\n")

# Check which packages are already installed in the target library
installed_in_lib <- installed.packages(lib.loc = personal_lib)[, "Package"]

# --- Install CRAN Packages ---
cat("--- Checking CRAN packages ---\n")
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed_in_lib)]
if (length(missing_cran) > 0) {
  cat("Installing from CRAN:", paste(missing_cran, collapse=", "), "\n")
  install.packages(missing_cran, lib = personal_lib)
} else {
  cat("All required CRAN packages are already installed.\n")
}

# --- Install Bioconductor Packages ---
cat("\n--- Checking Bioconductor packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = personal_lib)
}
# Make sure we check against all installed packages, not just those in our personal lib
all_installed <- installed.packages()[,"Package"]
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% all_installed)]
if (length(missing_bioc) > 0) {
  cat("Installing from Bioconductor:", paste(missing_bioc, collapse=", "), "\n")
  BiocManager::install(missing_bioc, lib = personal_lib, update = FALSE, ask = FALSE)
} else {
  cat("All required Bioconductor packages are already installed.\n")
}

# --- Install GitHub Packages ---
cat("\n--- Checking GitHub packages ---\n")
gh_package_names <- names(github_pkgs)
missing_gh <- gh_package_names[!(gh_package_names %in% all_installed)]

if (length(missing_gh) > 0) {
  repos_to_install <- github_pkgs[missing_gh]
  cat("Installing from GitHub:", paste(names(repos_to_install), collapse=", "), "\n")
  for (pkg_name in names(repos_to_install)) {
    repo <- repos_to_install[pkg_name]
    cat("Installing", pkg_name, "from", repo, "...\n")
    if (pkg_name == "cicero") {
      remotes::install_github(repo, ref = "monocle3", lib = personal_lib, upgrade = "never", force = TRUE)
    } else {
      remotes::install_github(repo, lib = personal_lib, upgrade = "never")
    }
  }
} else {
  cat("All required GitHub packages are already installed.\n")
}

# 4. FINAL VERIFICATION REPORT
cat("\n--- FINAL STATUS REPORT ---\n")
all_pkgs_to_check <- unique(c(cran_pkgs, bioc_pkgs, names(github_pkgs)))

final_check <- data.frame(
  Package = all_pkgs_to_check,
  Installed = sapply(all_pkgs_to_check, function(p) requireNamespace(p, quietly = TRUE))
)
rownames(final_check) <- NULL
print(final_check)

failed <- final_check$Package[!final_check$Installed]
if (length(failed) > 0) {
  cat("\n[!] The following packages failed to install or load:", paste(failed, collapse=", "), "\n")
  cat("[!] Please check the error messages above for details.\n")
} else {
  cat("\n[+] Success! Your R environment for single-cell analysis is ready.\n")
}
