# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 0: R ENVIRONMENT SETUP (FORCE MODE)
# Version: 3.0 (Red Hat 9 Optimized)
# =============================================================================

# 0. THE "FORCE" PROTOCOL - Locking the environment to official system paths
Sys.setenv(LD_LIBRARY_PATH = "/usr/lib64")
Sys.setenv(TMPDIR = paste0(Sys.getenv("HOME"), "/Rtemp"))
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
  "celda", "dittoSeq", "AUCell", "Gviz", "GenomicRanges", "rtracklayer",
  "BiocSingular", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder"
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

# 3. FINAL VERIFICATION REPORT
cat("\n--- FINAL STATUS REPORT ---\n")
all_pkgs_to_check <- unique(c(cran_pkgs, bioc_pkgs, names(github_pkgs)))
final_check <- data.frame(
  Package = all_pkgs_to_check,
  Status = sapply(all_pkgs_to_check, is_installed)
)
print(final_check)
