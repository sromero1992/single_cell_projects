# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 13: scTenifoldNet DIFFERENTIAL REGULATION
# =============================================================================
# Compares gene-regulatory networks between two conditions per cell type and
# reports differentially regulated (DR) genes, then enriches them (Enrichr).
#
# Method (per cell type, per contrast X vs Y):
#   1. counts -> keep genes expressed in >10% of cells in EACH group
#   2. pcNet on X and on Y  (principal-component-regression networks)
#   3. manifoldAlignment(X, Y)  -> shared low-dim embedding
#   4. dRegulation()            -> per-gene distance + p-value (DR genes)
#   5. Enrichr on the top significant DR genes
#
# Generalises the original 4 near-identical loops (broad / subtypes / colonocyte
# subtypes / T-cell subtypes) into one function driven by config.
#
# INPUT : annotated Seurat .rds with CELLTYPE_COLUMN + CONDITION_COLUMN.
# OUTPUT: <OUTPUT_DIR>/scTenifoldNet_results/<contrast>/<celltype>.xlsx
# =============================================================================
suppressPackageStartupMessages({
  library(scTenifoldNet)
  library(Seurat)
  library(Matrix)
  library(openxlsx)
  library(dplyr)
  library(enrichR)
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
RDS_PATH     <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))

# Grouping column that holds the labels in CELL_TYPES below.
CELLTYPE_COLUMN  <- "CellType"
CONDITION_COLUMN <- "Genotype_sex"

# Focus the run on specific cell types (must match CELLTYPE_COLUMN labels exactly).
# Set to NULL to run EVERY cell type in the column.
CELL_TYPES <- c("Tumor epithelium", "Stem cells", "TA cells")

# Optional: restrict to one lineage before splitting into its subtypes, mirroring
# the original "subcolonocytes"/"subtcells" loops. NULL = use all cells.
#   e.g. SUBSET_TO <- list(column = "CellType_broad", value = "Colonocytes")
SUBSET_TO <- NULL

# Two-group comparisons: each is c(X, Y). The network of X is aligned to Y and DR
# genes are reported; keep the perturbed/treatment group FIRST (as X = 'Nr4a1 KO'
# was in the original).
CONTRASTS_LIST <- list(
  Nr4a1_KO_polyp_vs_polyp_Female = c("Polyp_NR4a1_KO_Female", "Polyp_Female"),
  Nr4a1_KO_polyp_vs_polyp_Male   = c("Polyp_NR4a1_KO_Male",   "Polyp_Male"),
  Polyp_vs_WT_Female             = c("Polyp_Female",          "WT_Female"),
  Polyp_vs_WT_Male               = c("Polyp_Male",            "WT_Male")
)

# scTenifoldNet parameters
N_CORES        <- 8
MANIFOLD_D     <- 30       # manifold alignment dimensions
DROPOUT_FRAC   <- 0.10     # keep genes expressed in > this fraction of cells / group
MIN_CELLS      <- 30       # skip a cell type x group with fewer cells than this

# RESUME: skip any contrast x cell type whose output xlsx already exists AND holds
# a non-empty 'scTenifoldNet' result sheet (>0 gene rows). Lets you restart after a
# crash/reboot without recomputing the pairs that already finished. A partially
# written / corrupt file (e.g. killed mid-save) fails the content check and is
# recomputed. Set FALSE to force a full recompute (overwrites everything).
SKIP_IF_DONE   <- TRUE

# Enrichr
RUN_ENRICHR   <- TRUE
ENRICHR_DBS   <- c("GO_Biological_Process_2025",
                   "GO_Molecular_Function_2025")
DR_TOP_N      <- 250       # rank DR genes by p, take top N, keep p < 0.05

OUT_DIR <- file.path(OUTPUT_DIR, "scTenifoldNet_results")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
set.seed(42)

# =============================================================================
# --- PART 2: HELPERS ---------------------------------------------------------
# =============================================================================
# Drop mito / unnamed / ORF / antisense genes (as in the original workflow).
clean_genes <- function(obj) {
  g <- rownames(obj)
  keep <- !(startsWith(g, "MT-") | startsWith(g, "Mt-") | startsWith(g, "mt-") |
            startsWith(g, "ENSG0") |
            grepl("orf", g, ignore.case = TRUE) |
            grepl("-AS", g, ignore.case = TRUE))
  obj[keep, ]
}

# Excel sheet names cap at 31 chars; keep them short + unique.
safe_sheet <- function(x) substr(gsub("[^A-Za-z0-9_]", "_", x), 1, 31)

# Is an output xlsx a COMPLETE, usable result? TRUE only if the file exists, is
# non-empty, opens, carries the 'scTenifoldNet' sheet, and that sheet has >=1 gene
# row. A truncated/corrupt file (interrupted save) or an empty sheet returns FALSE
# so the pair is recomputed rather than wrongly skipped.
is_complete_output <- function(f, sheet = "scTenifoldNet", min_rows = 1) {
  if (!file.exists(f) || file.size(f) == 0) return(FALSE)
  isTRUE(tryCatch({
    if (!sheet %in% openxlsx::getSheetNames(f)) return(FALSE)
    df <- openxlsx::read.xlsx(f, sheet = sheet)
    !is.null(df) && nrow(df) >= min_rows && ncol(df) >= 1
  }, error = function(e) FALSE))
}

if (RUN_ENRICHR) {
  ok <- tryCatch({ setEnrichrSite("Enrichr"); TRUE },
                 error = function(e) { message("  [WARN] Enrichr offline; skipping."); FALSE })
  RUN_ENRICHR <- ok
}

# Core: one contrast (X vs Y) for one cell group -> xlsx.
run_sctenifold_pair <- function(counts_X, counts_Y, out_file, label) {
  # Dropout filter within each group.
  counts_X <- counts_X[Matrix::rowSums(counts_X > 0) / ncol(counts_X) > DROPOUT_FRAC, , drop = FALSE]
  counts_Y <- counts_Y[Matrix::rowSums(counts_Y > 0) / ncol(counts_Y) > DROPOUT_FRAC, , drop = FALSE]
  shared   <- intersect(rownames(counts_X), rownames(counts_Y))
  if (length(shared) < 50) { message("    [SKIP] ", label, ": only ", length(shared), " shared genes."); return(invisible(NULL)) }
  message("    ", label, ": ", length(shared), " shared genes | X=", ncol(counts_X), " Y=", ncol(counts_Y), " cells")
  X <- counts_X[shared, ]; Y <- counts_Y[shared, ]

  wb <- createWorkbook()
  tryCatch({
    set.seed(1); xNet <- pcNet(X, nCores = N_CORES); message("      net X done")
    set.seed(1); yNet <- pcNet(Y, nCores = N_CORES); message("      net Y done")
    set.seed(1)
    mA <- manifoldAlignment(xNet, yNet, d = MANIFOLD_D, nCores = N_CORES)
    rownames(mA) <- c(paste0("X_", shared), paste0("Y_", shared))
    dR <- dRegulation(manifoldOutput = mA)

    addWorksheet(wb, "scTenifoldNet", gridLines = FALSE)
    writeDataTable(wb, "scTenifoldNet", dR)

    if (RUN_ENRICHR) {
      glist <- dR %>% arrange(p.value) %>% head(DR_TOP_N) %>%
        filter(p.value < 0.05) %>% pull(gene)
      addWorksheet(wb, "DR_genes_for_Enrichr", gridLines = FALSE)
      writeDataTable(wb, "DR_genes_for_Enrichr", dR[dR$gene %in% glist, ])
      if (length(glist) >= 5) {
        er <- tryCatch(enrichr(glist, ENRICHR_DBS), error = function(e) NULL)
        if (!is.null(er)) for (db in names(er)) {
          df <- er[[db]]
          if (!is.null(df) && nrow(df) > 0) {
            df <- df[df$Adjusted.P.value < 0.05, , drop = FALSE]
            sn <- safe_sheet(paste0("enr_", db))
            addWorksheet(wb, sn, gridLines = FALSE)
            writeDataTable(wb, sn, df)
          }
        }
      } else message("      [note] <5 significant DR genes; Enrichr skipped.")
    }
    saveWorkbook(wb, out_file, overwrite = TRUE)
    message("      saved -> ", basename(out_file))
  }, error = function(e) cat("    ERROR", label, ":", conditionMessage(e), "\n"))
}

# =============================================================================
# --- PART 3: LOAD + CLEAN ----------------------------------------------------
# =============================================================================
message("=== Loading Seurat object ===")
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"
if (inherits(so[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(so[["RNA"]], search = "counts")) > 1) so <- JoinLayers(so)

if (!is.null(SUBSET_TO)) {
  keep <- colnames(so)[as.character(so@meta.data[[SUBSET_TO$column]]) == SUBSET_TO$value]
  so <- subset(so, cells = keep)
  message("  Restricted to ", SUBSET_TO$column, " == ", SUBSET_TO$value, ": ", ncol(so), " cells")
}
so <- clean_genes(so)
message("  ", ncol(so), " cells x ", nrow(so), " genes after cleaning")
gc()

for (col in c(CELLTYPE_COLUMN, CONDITION_COLUMN))
  if (!col %in% colnames(so@meta.data)) stop("Missing column: ", col, call. = FALSE)

# =============================================================================
# --- PART 4: RUN (contrast x cell type) --------------------------------------
# =============================================================================
ct_vec   <- as.character(so@meta.data[[CELLTYPE_COLUMN]])
cond_vec <- as.character(so@meta.data[[CONDITION_COLUMN]])
cell_types <- sort(unique(ct_vec[!is.na(ct_vec)]))
# Focus on the requested cell types (if any).
if (!is.null(CELL_TYPES)) {
  miss <- setdiff(CELL_TYPES, cell_types)
  if (length(miss) > 0)
    message("  [WARN] CELL_TYPES not found in ", CELLTYPE_COLUMN, ": ",
            paste(miss, collapse = ", "),
            "\n         Present: ", paste(cell_types, collapse = " | "))
  cell_types <- intersect(CELL_TYPES, cell_types)
  if (length(cell_types) == 0)
    stop("None of CELL_TYPES are present in ", CELLTYPE_COLUMN, ".", call. = FALSE)
}
message("  Cell types to run: ", paste(cell_types, collapse = ", "))

for (contrast in names(CONTRASTS_LIST)) {
  gX <- CONTRASTS_LIST[[contrast]][1]; gY <- CONTRASTS_LIST[[contrast]][2]
  cdir <- file.path(OUT_DIR, safe_sheet(contrast))
  if (!dir.exists(cdir)) dir.create(cdir, recursive = TRUE)
  message("\n=== Contrast: ", contrast, "  (X=", gX, " vs Y=", gY, ") ===")

  for (CT in cell_types) {
    cells_X <- colnames(so)[ct_vec == CT & cond_vec == gX]
    cells_Y <- colnames(so)[ct_vec == CT & cond_vec == gY]
    if (length(cells_X) < MIN_CELLS || length(cells_Y) < MIN_CELLS) {
      message("  [SKIP] ", CT, ": X=", length(cells_X), " Y=", length(cells_Y),
              " (< MIN_CELLS=", MIN_CELLS, ")")
      next
    }
    out_file <- file.path(cdir, paste0(gsub("[^A-Za-z0-9]+", "_", CT), "_sctenifold.xlsx"))
    if (SKIP_IF_DONE && is_complete_output(out_file)) {
      message("  [SKIP done] ", CT, " -> ", basename(out_file), " already complete.")
      next
    }
    message("  ", CT)
    cX <- GetAssayData(so, assay = "RNA", layer = "counts")[, cells_X, drop = FALSE]
    cY <- GetAssayData(so, assay = "RNA", layer = "counts")[, cells_Y, drop = FALSE]
    run_sctenifold_pair(cX, cY, out_file, paste0(CT, " | ", contrast))
    rm(cX, cY); gc()
  }
}

message("\n=== Script 13 complete ===")
message("  Results: ", OUT_DIR)
message("  Re-run with CELLTYPE_COLUMN <- \"CellType\" (or SUBSET_TO a lineage) for subtypes.")
