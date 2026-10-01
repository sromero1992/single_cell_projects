# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 12a: EXPORT FOR pySCENIC (TF GRN)
# =============================================================================
# Exports raw-count expression matrices from a Seurat object into loom files that
# pySCENIC (grn -> ctx -> aucell) reads directly. No ATAC/peaks needed: pySCENIC
# infers regulons from scRNA alone using the pre-built cisTarget motif databases.
#
# Produces one loom PER cell type (Tumor epithelium / Stem cells / TA cells) and,
# optionally, one COMBINED loom of all three (recommended for GRN inference: more
# cells + the shared stem -> TA -> tumor regulatory context give cleaner networks).
# It also writes a run_pyscenic.sh with the exact grn/ctx/aucell commands wired to
# your database filenames.
#
# INPUT : a Seurat .rds carrying the CELLTYPE_COLUMN annotation (unified object).
# OUTPUT: <OUTPUT_DIR>/SCENIC/<label>/<label>.loom  + run_pyscenic.sh
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# =============================================================================
# --- PART 1: USER CONFIGURATION ----------------------------------------------
# =============================================================================
# ---- Shared portable config: nr4a1 defaults, override via env vars (see config.R).
# Run from the pipeline directory, or set NR4A1_CONFIG=/full/path/to/config.R. ----
.NR4A1_CFG <- Sys.getenv("NR4A1_CONFIG", "config.R")
if (!file.exists(.NR4A1_CFG)) stop("config.R not found at '", .NR4A1_CFG,
  "' - cd to the pipeline directory or set NR4A1_CONFIG.", call. = FALSE)
source(.NR4A1_CFG)

# PROJECT_NAME <- "Nr4a1_s17_ack"   # [portable] now set in config.R
# ROOT_PATH    <- "/home/ssromerogon/local_drive/optimus_drive/selim_working_dir/2026_nr4a1_ack/r_process"   # [portable] now set in config.R
# OUTPUT_DIR   <- file.path(ROOT_PATH, "seurat_output")   # [portable] now set in config.R

# Annotated object (the unified one carrying the cell-type column below).
RDS_PATH <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))

CELLTYPE_COLUMN <- "CellType"          # column holding the labels below
SAMPLE_COLUMN   <- "SampleID"
CONDITION_COLUMN<- "Genotype_sex"

# The populations to run SCENIC on.
# NOTE: adding CD8 T cells only -- the other 3 are already exported/run, so we
# don't re-export them. Set this back to all 4 if you ever rebuild from scratch.
# Use the EXACT label from the object (check: unique(data$CellType)).
CELL_TYPES <- c("CD8+ T cells")   # <- fix spelling to match your object if different
# CELL_TYPES <- c("Tumor epithelium", "Stem cells", "TA cells", "CD8+ T cells")

EXPORT_PER_TYPE  <- TRUE   # one loom per cell type (3 runs)
EXPORT_COMBINED  <- TRUE   # one loom of all three (recommended, single richer run)
COMBINED_LABEL   <- "TumorEpi_Stem_TA"

# SPLIT_BY: how the SCENIC networks are defined, beyond one-per-cell-type.
#   "none"      -> one network per cell type on ALL cells. Comparable AUCell across
#                  conditions; compare Nr4a1 activity by Genotype_sex afterwards.
#   "sex"       -> one network per cell type x SEX (2 groups). RECOMMENDED when sex
#                  strongly modifies the effect: enough cells for a stable network,
#                  and you compare genotype (WT/Polyp/KO) activity WITHIN each sex.
#   "condition" -> one network per cell type x Genotype_sex (up to 6). Exploratory
#                  REWIRING only; small-n and not directly comparable across groups.
# Genotype_sex is ALWAYS stored as a loom annotation regardless of the split, and
# any subset below MIN_CELLS_SCENIC is skipped.
SPLIT_BY <- "sex"           # "none" | "sex" | "condition"
MIN_CELLS_SCENIC <- 200     # skip a subset with fewer cells than this

# Where the exports go, and where your cisTarget databases live.
SCENIC_DIR    <- file.path(OUTPUT_DIR, "SCENIC")
# CISTARGET_DIR <- "/home/ssromerogon/cisTarget_databases"   # [portable] now set in config.R
TF_FILE       <- "allTFs_mm.txt"
DB_FEATHERS   <- c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
                   "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_TBL     <- "motifs-v9-nr.mgi-m0.001-o0.0.tbl"
N_WORKERS     <- 8         # cores for the printed pyscenic commands

# =============================================================================
# --- PART 2: LOAD + SANITY ---------------------------------------------------
# =============================================================================
if (!dir.exists(SCENIC_DIR)) dir.create(SCENIC_DIR, recursive = TRUE)
message("=== Loading Seurat object ===")
data <- readRDS(RDS_PATH)
DefaultAssay(data) <- "RNA"
# Seurat v5: collapse split layers so counts is a single matrix.
if (inherits(data[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(data[["RNA"]], search = "counts")) > 1) {
  data <- JoinLayers(data)
}
message("  Loaded: ", ncol(data), " cells x ", nrow(data), " genes")

if (!CELLTYPE_COLUMN %in% colnames(data@meta.data))
  stop("Column '", CELLTYPE_COLUMN, "' not found in metadata.", call. = FALSE)

present <- unique(as.character(data@meta.data[[CELLTYPE_COLUMN]]))
missing <- setdiff(CELL_TYPES, present)
if (length(missing) > 0)
  message("  [WARN] Not found in ", CELLTYPE_COLUMN, ": ",
          paste(missing, collapse = ", "),
          "\n         Present labels: ", paste(sort(present), collapse = " | "))
CELL_TYPES <- intersect(CELL_TYPES, present)
if (length(CELL_TYPES) == 0) stop("None of the requested cell types are present.", call. = FALSE)

# =============================================================================
# --- PART 3: LOOM EXPORT -----------------------------------------------------
# =============================================================================
# Writes raw COUNTS (genes x cells, symbol rownames) - what GRNBoost2 and AUCell
# expect. Primary: SCopeLoomR (pySCENIC-native). Fallback: gzipped cells x genes
# CSV (also readable by `pyscenic grn`), with a clear note.
have_loom <- requireNamespace("SCopeLoomR", quietly = TRUE)
if (!have_loom)
  message("\n  [NOTE] SCopeLoomR not installed - will write CSV.gz instead.\n",
          "         For native loom: remotes::install_github('aertslab/SCopeLoomR')\n")

export_one <- function(obj, label) {
  out_dir <- file.path(SCENIC_DIR, label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  counts <- as(counts, "CsparseMatrix")                       # genes x cells
  # Drop all-zero genes (SCENIC ignores them and they bloat the file).
  keep <- Matrix::rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  message(sprintf("  %-22s %d cells x %d genes", label, ncol(counts), nrow(counts)))

  meta <- obj@meta.data[, intersect(c(CELLTYPE_COLUMN, SAMPLE_COLUMN, CONDITION_COLUMN),
                                    colnames(obj@meta.data)), drop = FALSE]

  if (have_loom) {
    loom_path <- file.path(out_dir, paste0(label, ".loom"))
    if (file.exists(loom_path)) file.remove(loom_path)
    loom <- SCopeLoomR::build_loom(
      file.name       = loom_path,
      dgem            = counts,                                # genes x cells, raw counts
      title           = label,
      default.embedding = NULL
    )
    for (cn in colnames(meta))
      SCopeLoomR::add_col_attr(loom = loom, key = cn,
                               value = as.character(meta[[cn]]), as.annotation = TRUE)
    SCopeLoomR::close_loom(loom)
    message("    -> ", loom_path)
    return(loom_path)
  } else {
    # Fallback: cells x genes matrix (pyscenic grn reads csv/tsv).
    csv_path <- file.path(out_dir, paste0(label, "_expr_cells_x_genes.csv.gz"))
    m <- Matrix::t(counts)                                     # cells x genes
    gz <- gzfile(csv_path, "w")
    utils::write.csv(as.matrix(m), gz)                         # rownames = cells, header = genes
    close(gz)
    message("    -> ", csv_path, "  (CSV fallback)")
    return(csv_path)
  }
}

exports <- list()   # label -> file path

# Export one subset with a cell-count guard (GRNBoost2 needs enough cells).
export_guarded <- function(cells, label) {
  if (length(cells) < MIN_CELLS_SCENIC) {
    message(sprintf("    [SKIP] %-30s only %d cells (< MIN_CELLS_SCENIC=%d)",
                    label, length(cells), MIN_CELLS_SCENIC))
    return(invisible(NULL))
  }
  obj <- subset(data, cells = cells)
  exports[[label]] <<- export_one(obj, label)
  rm(obj); gc()
}

# Per-cell grouping vector for the chosen split ("sex" derives Female/Male from
# the Genotype_sex suffix; "condition" uses Genotype_sex itself).
split_grp <- rep(NA_character_, ncol(data))
if (SPLIT_BY != "none" && CONDITION_COLUMN %in% colnames(data@meta.data)) {
  gs <- as.character(data@meta.data[[CONDITION_COLUMN]])
  split_grp <- if (SPLIT_BY == "sex")
    ifelse(grepl("_Female$", gs), "Female",
    ifelse(grepl("_Male$",   gs), "Male", NA_character_))
  else gs
}

if (EXPORT_PER_TYPE) {
  message("\n=== Per-cell-type looms (split: ", SPLIT_BY, ") ===")
  ct_vec <- as.character(data@meta.data[[CELLTYPE_COLUMN]])
  for (ct in CELL_TYPES) {
    lab <- gsub("[^A-Za-z0-9]+", "_", ct)
    if (SPLIT_BY == "none") {
      export_guarded(colnames(data)[ct_vec == ct], lab)
    } else {
      for (lv in sort(unique(split_grp[!is.na(split_grp)]))) {
        cc <- colnames(data)[ct_vec == ct & !is.na(split_grp) & split_grp == lv]
        export_guarded(cc, paste0(lab, "__", gsub("[^A-Za-z0-9]+", "_", lv)))
      }
    }
  }
}

if (EXPORT_COMBINED && length(CELL_TYPES) > 1) {
  message("\n=== Combined loom (", COMBINED_LABEL, ") ===")
  cells <- colnames(data)[as.character(data@meta.data[[CELLTYPE_COLUMN]]) %in% CELL_TYPES]
  obj   <- subset(data, cells = cells)
  exports[[COMBINED_LABEL]] <- export_one(obj, COMBINED_LABEL)
  rm(obj); gc()
}

# =============================================================================
# --- PART 4: WRITE pySCENIC COMMANDS -----------------------------------------
# =============================================================================
# Emits a run_pyscenic.sh with grn -> ctx -> aucell for each export, wired to the
# cisTarget database filenames. Run it from the activated scanpy/pyscenic env.
sh_path <- file.path(SCENIC_DIR, "run_pyscenic.sh")
con <- file(sh_path, "w")
writeLines(c(
  "#!/usr/bin/env bash",
  "# Auto-generated by 12a_export_for_pyscenic.R",
  "# Run from the activated env:  conda activate scanpy_env_311 && bash run_pyscenic.sh",
  "set -euo pipefail",
  "",
  "# System libstdc++ lacks GLIBCXX_3.4.30 that scipy/arboreto need; preload the",
  "# env's newer one so the pyscenic CLI can import (same fix as Script 10).",
  "export LD_PRELOAD=\"${CONDA_PREFIX}/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}\"",
  "",
  paste0("DB_DIR=\"", CISTARGET_DIR, "\""),
  paste0("TFS=\"$DB_DIR/", TF_FILE, "\""),
  paste0("DB1=\"$DB_DIR/", DB_FEATHERS[1], "\""),
  paste0("DB2=\"$DB_DIR/", DB_FEATHERS[2], "\""),
  paste0("MOTIF=\"$DB_DIR/", MOTIF_TBL, "\""),
  paste0("NW=", N_WORKERS),
  ""
), con)
for (lab in names(exports)) {
  ex   <- exports[[lab]]
  d    <- dirname(ex)
  is_loom <- grepl("\\.loom$", ex)
  expr_flag <- if (is_loom) paste0("\"", ex, "\"") else paste0("\"", ex, "\"")
  writeLines(c(
    paste0("### ---- ", lab, " ----"),
    paste0("pyscenic grn ", expr_flag, " \"$TFS\" \\"),
    paste0("  -o \"", file.path(d, "adj.tsv"), "\" --num_workers $NW --method grnboost2"),
    "",
    paste0("pyscenic ctx \"", file.path(d, "adj.tsv"), "\" \"$DB1\" \"$DB2\" \\"),
    paste0("  --annotations_fname \"$MOTIF\" \\"),
    paste0("  --expression_mtx_fname ", expr_flag, " \\"),
    paste0("  -o \"", file.path(d, "regulons.csv"), "\" --mask_dropouts --num_workers $NW"),
    "",
    paste0("pyscenic aucell ", expr_flag, " \"", file.path(d, "regulons.csv"), "\" \\"),
    paste0("  -o \"", file.path(d, paste0(lab, "_pyscenic.loom")), "\" --num_workers $NW"),
    ""
  ), con)
}
close(con)
Sys.chmod(sh_path, mode = "0755")

message("\n=== Done ===")
message("  Looms/CSVs under: ", SCENIC_DIR)
message("  pySCENIC commands: ", sh_path)
message("\n  Next (in the activated env):")
message("    conda activate scanpy_env_311")
message("    bash ", sh_path)
message("\n  Then read <label>_pyscenic.loom back for regulon activity (AUCell),")
message("  and compare Nr4a1 regulon activity KO vs WT.")
