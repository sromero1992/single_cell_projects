# =============================================================================
# scRNA-seq PIPELINE - PACKAGE CHECKER
# Version: 1.0 (UNIFIED build)
#
# PURPOSE:
#   Fast, read-only audit of which pipeline dependencies are installed, grouped
#   by the script that needs them. Then, optionally, install ONLY what is
#   missing - nothing already present is touched or upgraded.
#
#   This is deliberately separate from 00_rlibs_installation.R:
#     00_rlibs_installation.R  - full setup, attempts installs, hard-stops
#     check_packages.R (this)  - reports in seconds, installs nothing by default
#
# HOW TO USE:
#   1. Run as-is. It only READS - nothing is installed, nothing changes.
#      Takes a few seconds. Read the report.
#   2. To install the missing ones, set INSTALL_MISSING <- TRUE below and
#      re-run. Only packages reported as missing are installed.
#
#   Safe to re-run at any time.
#
# NOTE ON PYTHON:
#   Script 10 (CellRank trajectories) needs reticulate plus a conda
#   environment. Since you are working on processing and annotation, that is
#   off by default. Set CHECK_PYTHON <- TRUE when you get to Script 10.
#   Installing the R package `reticulate` is harmless either way - it is the
#   conda environment that takes real setup, and this script never touches it.
# =============================================================================

# --- CONFIGURATION -----------------------------------------------------------

INSTALL_MISSING <- FALSE   # FALSE = report only (recommended first run)
                           # TRUE  = install missing packages, skip installed

CHECK_PYTHON    <- FALSE   # TRUE only when you start on Script 10

INCLUDE_OPTIONAL <- TRUE   # TRUE  = also report optional/feature-gated packages
                           # FALSE = only what is strictly required to run

options(timeout = 1200)    # GitHub tarballs routinely exceed the 60s default

# =============================================================================
# --- PACKAGE REGISTRY (do not edit unless adding a dependency) ---------------
# =============================================================================
# src: "cran" | "bioc" | "github"
# repo: for github, "user/repo" (optionally with subdir given separately)
# need: "required" | "optional"
# used: which script(s) consume it

pkgs <- list(
  # ---- Core: needed by almost everything -----------------------------------
  list(name="Seurat",         src="cran", used="01-10", need="required"),
  list(name="SeuratObject",   src="cran", used="01-10", need="required"),
  list(name="dplyr",          src="cran", used="01-10", need="required"),
  list(name="ggplot2",        src="cran", used="01-10", need="required"),
  list(name="patchwork",      src="cran", used="01-10", need="required"),
  list(name="Matrix",         src="cran", used="01-10", need="required"),
  list(name="tidyr",          src="cran", used="02-10", need="required"),
  list(name="tibble",         src="cran", used="02-10", need="required"),
  list(name="stringr",        src="cran", used="02-10", need="required"),
  list(name="cowplot",        src="cran", used="02-06", need="required"),
  list(name="ggpubr",         src="cran", used="02-07", need="required"),
  list(name="openxlsx",       src="cran", used="01",    need="required"),
  list(name="writexl",        src="cran", used="01-10", need="required"),
  list(name="readxl",         src="cran", used="01",    need="optional"),
  list(name="hdf5r",          src="cran", used="01",    need="required",
       note="Reads the 10x .h5 files. Needs system HDF5 dev headers."),
  list(name="harmony",        src="cran", used="01",    need="required",
       note="Batch integration."),
  list(name="R.utils",        src="cran", used="01-10", need="optional"),
  list(name="remotes",        src="cran", used="setup", need="required"),
  list(name="devtools",       src="cran", used="setup", need="required"),
  list(name="BiocManager",    src="cran", used="setup", need="required"),

  # ---- Script 01: processing ------------------------------------------------
  list(name="celda",          src="bioc", used="01", need="required",
       note="DecontX ambient RNA correction (RUN_DECONTX)."),
  list(name="DoubletFinder",  src="github", repo="chris-mcginnis-ucsf/DoubletFinder",
       used="01", need="required",
       note="DEFAULT doublet caller. Must be the post-2023 API (paramSweep, not paramSweep_v3)."),
  list(name="scDblFinder",    src="bioc", used="01", need="optional",
       note="Alternative doublet caller (DOUBLET_METHOD)."),
  list(name="SingleCellExperiment", src="bioc", used="01", need="optional",
       note="Required by scDblFinder."),
  list(name="SummarizedExperiment", src="bioc", used="01", need="optional"),
  list(name="BiocParallel",   src="bioc", used="01", need="optional",
       note="Required by scDblFinder."),
  list(name="xgboost",        src="cran", used="01", need="optional",
       note="Required by scDblFinder."),
  list(name="SCEVAN",         src="github", repo="AntonioDeFalco/SCEVAN",
       used="01", need="optional",
       note="CNA analysis (RUN_SCEVAN). Set FALSE if no tumour samples."),
  list(name="yaGST",          src="github", repo="miccec/yaGST",
       used="01", need="optional", note="SCEVAN dependency."),
  list(name="BPCells",        src="github", repo="bnprks/BPCells/r",
       used="01", need="optional",
       note="On-disk matrices (USE_BPCELLS). Needs system HDF5 dev headers. Set USE_BPCELLS <- FALSE to skip."),
  list(name="SeuratWrappers", src="github", repo="satijalab/seurat-wrappers",
       used="01-06", need="optional",
       note="NOTE capital S. Lowercase 'seuratWrappers' never matches an install."),

  # ---- Scripts 02-06: annotation -------------------------------------------
  list(name="presto",         src="github", repo="immunogenomics/presto",
       used="02-06", need="optional",
       note="Large speedup for FindMarkers / FindAllMarkers."),
  list(name="ComplexHeatmap", src="bioc", used="02-08", need="optional"),
  list(name="dittoSeq",       src="bioc", used="02-06", need="optional"),

  # ---- Script 07: DE --------------------------------------------------------
  list(name="MAST",           src="bioc", used="07", need="optional",
       note="MAST differential expression."),
  list(name="enrichR",        src="cran", used="07", need="optional",
       note="Enrichr ORA. Needs internet at run time."),
  list(name="fgsea",          src="github", repo="alserglab/fgsea",
       used="07", need="optional"),
  list(name="SplineDV",       src="github", repo="Xenon8778/SplineDV",
       used="07", need="optional", note="Differential variability."),
  list(name="broom",          src="cran", used="07", need="optional",
       note="ANOVA tidying."),
  list(name="AUCell",         src="bioc", used="07", need="optional"),

  # ---- Script 08: CellChat --------------------------------------------------
  list(name="CellChat",       src="github", repo="jinworks/CellChat",
       used="08", need="optional",
       note="Install ComplexHeatmap from Bioconductor FIRST or this fails confusingly."),
  list(name="BiocNeighbors",  src="bioc", used="08", need="optional"),

  # ---- Script 09: cell scores ----------------------------------------------
  list(name="CytoTRACE2",     src="github", repo="digitalcytometry/cytotrace2",
       subdir="cytotrace2_r", used="09", need="optional",
       note="Needs subdir='cytotrace2_r'. Without it you get 'does not appear to be a package'."),
  list(name="CytoTRACE",      src="github", repo="gunsagargulati/CytoTRACE",
       used="09", need="optional",
       note="Old package, stale deps. Skippable: set RUN_CYTOTRACE1 <- FALSE."),

  # ---- Script 10: trajectory (Python) --------------------------------------
  list(name="reticulate",     src="cran", used="10", need="python",
       note="R side only. The conda env is separate - see INSTALL_NOTES.md section 5.")
)

# =============================================================================
# --- AUDIT (read-only) -------------------------------------------------------
# =============================================================================

is_installed <- function(p) requireNamespace(p, quietly = TRUE)

# Filter by what the user asked to see
keep <- vapply(pkgs, function(p) {
  if (p$need == "python" && !CHECK_PYTHON) return(FALSE)
  if (p$need == "optional" && !INCLUDE_OPTIONAL) return(FALSE)
  TRUE
}, logical(1))
pkgs_use <- pkgs[keep]

status <- vapply(pkgs_use, function(p) is_installed(p$name), logical(1))

report <- data.frame(
  Package   = vapply(pkgs_use, function(p) p$name, character(1)),
  Installed = status,
  Source    = vapply(pkgs_use, function(p) p$src,  character(1)),
  Needed_By = vapply(pkgs_use, function(p) p$used, character(1)),
  Priority  = vapply(pkgs_use, function(p) p$need, character(1)),
  stringsAsFactors = FALSE
)

cat("\n")
cat("=====================================================================\n")
cat(" scRNA-seq PIPELINE - DEPENDENCY REPORT\n")
cat("=====================================================================\n")
cat(sprintf(" R version : %s\n", paste0(R.version$major, ".", R.version$minor)))
cat(sprintf(" Library   : %s\n", .libPaths()[1]))
cat(sprintf(" Checked   : %d packages\n\n", nrow(report)))

missing_req <- report[!report$Installed & report$Priority == "required", ]
missing_opt <- report[!report$Installed & report$Priority == "optional", ]
missing_py  <- report[!report$Installed & report$Priority == "python", ]

if (nrow(missing_req) == 0) {
  cat(" [OK] All REQUIRED packages are installed.\n")
} else {
  cat(sprintf(" [!!] %d REQUIRED package(s) MISSING:\n", nrow(missing_req)))
  for (i in seq_len(nrow(missing_req))) {
    cat(sprintf("        - %-18s (%s, needed by script %s)\n",
                missing_req$Package[i], missing_req$Source[i], missing_req$Needed_By[i]))
  }
}
cat("\n")

if (nrow(missing_opt) > 0) {
  cat(sprintf(" [--] %d OPTIONAL package(s) missing (features they gate will be skipped):\n",
              nrow(missing_opt)))
  for (i in seq_len(nrow(missing_opt))) {
    p <- pkgs_use[[which(vapply(pkgs_use, function(x) x$name, character(1)) == missing_opt$Package[i])[1]]]
    cat(sprintf("        - %-18s (%s, script %s)\n",
                p$name, p$src, p$used))
    if (!is.null(p$note)) cat(sprintf("            %s\n", p$note))
  }
  cat("\n")
}

if (CHECK_PYTHON && nrow(missing_py) > 0) {
  cat(" [--] Python-side (Script 10) missing:\n")
  for (i in seq_len(nrow(missing_py))) cat(sprintf("        - %s\n", missing_py$Package[i]))
  cat("\n")
}

# --- Readiness by script -----------------------------------------------------
# The practical question: what can I actually run right now?
cat("---------------------------------------------------------------------\n")
cat(" CAN I RUN IT?\n")
cat("---------------------------------------------------------------------\n")

ready <- function(names_needed) {
  miss <- names_needed[!vapply(names_needed, is_installed, logical(1))]
  if (length(miss) == 0) "READY" else paste0("BLOCKED (", paste(miss, collapse=", "), ")")
}

script_reqs <- list(
  "01 process_data"       = c("Seurat","dplyr","ggplot2","patchwork","Matrix",
                              "openxlsx","writexl","hdf5r","celda","harmony","DoubletFinder"),
  "02 global_annotation"  = c("Seurat","dplyr","ggplot2","patchwork","writexl"),
  "03-05 subannotation"   = c("Seurat","dplyr","ggplot2","patchwork","writexl"),
  "06 annotation_unifier" = c("Seurat","dplyr"),
  "07 DE + ANOVA"         = c("Seurat","dplyr","tidyr","tibble","writexl","MAST"),
  "08 cellchat"           = c("Seurat","CellChat"),
  "09 cell_scores"        = c("Seurat","dplyr","ggplot2","writexl"),
  "10 trajectory"         = c("Seurat","reticulate")
)
for (s in names(script_reqs)) {
  if (s == "10 trajectory" && !CHECK_PYTHON) {
    cat(sprintf("  %-24s %s\n", s, "SKIPPED (CHECK_PYTHON = FALSE)")); next
  }
  cat(sprintf("  %-24s %s\n", s, ready(script_reqs[[s]])))
}
cat("\n")
cat(" Note: Script 09 runs its entropy metrics with no extra packages.\n")
cat("       CytoTRACE2 / CytoTRACE only gate their own methods.\n\n")

# --- Functional checks -------------------------------------------------------
# Being installed is not the same as being the right version.
cat("---------------------------------------------------------------------\n")
cat(" VERSION / API CHECKS\n")
cat("---------------------------------------------------------------------\n")

if (is_installed("Seurat")) {
  v <- as.character(utils::packageVersion("Seurat"))
  cat(sprintf("  Seurat %s %s\n", v,
              if (utils::packageVersion("Seurat") >= "5.0.0") "(v5 - expected)" else "(v4 - pipeline assumes v5 layers)"))
}
if (is_installed("DoubletFinder")) {
  # The 2023 rewrite renamed paramSweep_v3 -> paramSweep. Script 01 uses the
  # NEW names, so an old install fails at run time with 'could not find function'.
  ns <- getNamespaceExports("DoubletFinder")
  if ("paramSweep" %in% ns) {
    cat("  DoubletFinder API: NEW (paramSweep) - correct for this pipeline\n")
  } else if ("paramSweep_v3" %in% ns) {
    cat("  DoubletFinder API: OLD (paramSweep_v3) - MUST UPGRADE:\n")
    cat("      remove.packages('DoubletFinder')\n")
    cat("      remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)\n")
  } else {
    cat("  DoubletFinder API: could not determine - inspect manually\n")
  }
}
if (is_installed("Matrix") && is_installed("Seurat")) {
  mv <- utils::packageVersion("Matrix")
  cat(sprintf("  Matrix %s\n", as.character(mv)))
}
cat("\n")

# =============================================================================
# --- INSTALL MISSING ONLY ----------------------------------------------------
# =============================================================================
to_install <- pkgs_use[!status]
# Never auto-install the Python-side package unless explicitly checking Python
to_install <- Filter(function(p) !(p$need == "python" && !CHECK_PYTHON), to_install)

if (!INSTALL_MISSING) {
  cat("=====================================================================\n")
  if (length(to_install) == 0) {
    cat(" Nothing missing. You are ready to go.\n")
  } else {
    cat(sprintf(" %d package(s) could be installed.\n", length(to_install)))
    cat(" This was a REPORT ONLY run - nothing was installed.\n")
    cat(" To install just those, set INSTALL_MISSING <- TRUE at the top and re-run.\n")
  }
  cat("=====================================================================\n\n")
} else {
  cat("=====================================================================\n")
  cat(sprintf(" INSTALLING %d MISSING PACKAGE(S)\n", length(to_install)))
  cat(" Already-installed packages are NOT touched or upgraded.\n")
  cat("=====================================================================\n\n")

  if (!is_installed("BiocManager")) install.packages("BiocManager", repos = "https://cran.rstudio.com")
  if (!is_installed("remotes"))     install.packages("remotes",     repos = "https://cran.rstudio.com")

  results <- list()
  for (p in to_install) {
    cat(sprintf("\n--- %s (%s) ---\n", p$name, p$src))
    ok <- tryCatch({
      if (p$src == "cran") {
        if (p$name == "hdf5r") {
          install.packages("hdf5r", type = "source",
                           configure.args = "--with-hdf5=/usr/bin/h5cc")
        } else {
          install.packages(p$name)
        }
      } else if (p$src == "bioc") {
        BiocManager::install(p$name, update = FALSE, ask = FALSE,
                             INSTALL_opts = "--no-test-load")
      } else if (p$src == "github") {
        if (!is.null(p$subdir)) {
          remotes::install_github(p$repo, subdir = p$subdir, upgrade = "never")
        } else {
          remotes::install_github(p$repo, upgrade = "never")
        }
      }
      is_installed(p$name)
    }, error = function(e) {
      cat(sprintf("  [FAILED] %s\n", conditionMessage(e)))
      FALSE
    })
    results[[p$name]] <- ok
    cat(sprintf("  -> %s\n", if (isTRUE(ok)) "OK" else "FAILED"))
  }

  cat("\n=====================================================================\n")
  cat(" INSTALL SUMMARY\n")
  cat("=====================================================================\n")
  okv <- unlist(results)
  if (length(okv)) {
    cat(sprintf("  succeeded: %d\n", sum(okv)))
    cat(sprintf("  failed   : %d\n", sum(!okv)))
    if (any(!okv)) {
      cat("\n  Failed packages:\n")
      for (n in names(okv)[!okv]) {
        pp <- Filter(function(x) x$name == n, pkgs_use)[[1]]
        cat(sprintf("    - %-18s (%s)\n", n, pp$src))
        if (!is.null(pp$note)) cat(sprintf("        %s\n", pp$note))
      }
      cat("\n  See INSTALL_NOTES.md for per-package troubleshooting.\n")
      cat("  Optional failures are survivable - the toggles above say how.\n")
    }
  }
  cat("\n  Re-run this script with INSTALL_MISSING <- FALSE to confirm.\n\n")
}

# --- Reminder about Python ---------------------------------------------------
if (!CHECK_PYTHON) {
  cat("---------------------------------------------------------------------\n")
  cat(" Python / Script 10 was NOT checked (CHECK_PYTHON = FALSE).\n")
  cat(" Nothing to do now - processing and annotation are pure R.\n")
  cat(" When you reach Script 10, set CHECK_PYTHON <- TRUE and follow\n")
  cat(" INSTALL_NOTES.md section 5 to build the conda environment.\n")
  cat("---------------------------------------------------------------------\n\n")
}
