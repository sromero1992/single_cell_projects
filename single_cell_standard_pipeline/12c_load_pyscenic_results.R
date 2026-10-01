# =============================================================================
# scRNA-seq PIPELINE - SCRIPT 12c: LOAD pySCENIC RESULTS BACK INTO R
# =============================================================================
# Reads the per-label regulon activity produced by 12b (grn -> ctx -> aucell),
# attaches it to the Seurat object as an "AUC" assay, and makes the standard
# regulon-activity plots. Focused on the Nr4a1 regulon but works for any TF.
#
# INPUT  (per label dir under <OUTPUT_DIR>/SCENIC/<label>/):
#   aucell.csv   -- cells x regulons continuous AUC matrix   (PRIMARY source)
#   regulons.p   -- pickled regulon gene sets (optional, for membership export)
#   <label>_pyscenic.loom -- only if 12b was run with EXPORT_LOOM=True (optional)
# Because our 12b defaults to EXPORT_LOOM=False, aucell.csv is the source of truth.
#
# OUTPUT (per label, under <OUTPUT_DIR>/SCENIC/<label>/R_plots/):
#   FeaturePlot_<TF>.png           -- TF regulon AUC on the subset's UMAP
#   Violin_<TF>_by_condition.png   -- TF regulon AUC by CONDITION_COLUMN
#   Heatmap_meanAUC_topRegulons.png-- mean AUC per condition (top variable regulons)
#   <label>_<TF>_AUC_percell.csv   -- per-cell TF AUC + metadata (for stats)
#   <label>_meanAUC_by_condition.csv
#
# NOTES
#   * Regulon names carry a "(+)" suffix (e.g. "Nr4a1(+)"); we match the focus TF
#     by prefix so you don't have to type it.
#   * High regulon AUC = the TF's inferred TARGETS are coordinately expressed, NOT
#     literal TF on/off. Binarization is a separate optional step (see 12b docs).
#   * Cell-name alignment is the one real gotcha: AUC cells are intersected with
#     colnames(seurat) before attaching.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# --- Shared portable config: nr4a1 defaults, override via env vars (config.R).
# Run from the pipeline directory, or set NR4A1_CONFIG=/full/path/to/config.R. ----
.NR4A1_CFG <- Sys.getenv("NR4A1_CONFIG", "config.R")
if (!file.exists(.NR4A1_CFG)) stop("config.R not found at '", .NR4A1_CFG,
  "' - cd to the pipeline directory or set NR4A1_CONFIG.", call. = FALSE)
source(.NR4A1_CFG)

# =============================================================================
# --- CONFIG ------------------------------------------------------------------
# =============================================================================
SCENIC_DIR <- Sys.getenv("NR4A1_SCENIC_DIR", file.path(OUTPUT_DIR, "SCENIC"))
# Seurat object to attach AUC to (same one 12a exported from).
RDS_PATH   <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_with_cell_scores.rds"))

FOCUS_TF         <- "Nr4a1"        # regulon of interest (matched as "<TF>(+)")
CELLTYPE_COLUMN  <- "CellType"
CONDITION_COLUMN <- "Genotype_sex" # grouping for violins / mean-AUC heatmap
UMAP_REDUCTION   <- "umap_harmony" # falls back to any umap_* / umap present
HEATMAP_TOP_N    <- 30             # most-variable regulons shown in the heatmap
DPI_SETTING      <- 300
# Only run specific labels (subfolders of SCENIC_DIR); NULL = all with aucell.csv.
LABELS           <- NULL

# =============================================================================
# --- HELPERS -----------------------------------------------------------------
# =============================================================================
# Read a label's AUC matrix -> regulons x cells (cells as columns). aucell.csv is
# written by pySCENIC as cells x regulons with the cell id in the first column.
read_auc <- function(label_dir) {
  f <- file.path(label_dir, "aucell.csv")
  if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
  m <- tryCatch(utils::read.csv(f, row.names = 1, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) return(NULL)
  t(as.matrix(m))                       # -> regulons x cells
}

discover_labels <- function() {
  if (!dir.exists(SCENIC_DIR)) stop("SCENIC_DIR not found: ", SCENIC_DIR, call. = FALSE)
  d <- list.dirs(SCENIC_DIR, recursive = FALSE, full.names = FALSE)
  d[vapply(d, function(x) file.exists(file.path(SCENIC_DIR, x, "aucell.csv")), logical(1))]
}

pick_umap <- function(obj) {
  reds <- SeuratObject::Reductions(obj)
  if (UMAP_REDUCTION %in% reds) return(UMAP_REDUCTION)
  cand <- grep("^umap", reds, value = TRUE, ignore.case = TRUE)
  if (length(cand)) cand[1] else NA_character_
}

# =============================================================================
# --- LOAD OBJECT -------------------------------------------------------------
# =============================================================================
message("=== Loading Seurat object ===")
if (!file.exists(RDS_PATH)) stop("RDS not found: ", RDS_PATH, call. = FALSE)
so <- readRDS(RDS_PATH)
DefaultAssay(so) <- "RNA"
if (inherits(so[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(so[["RNA"]], search = "counts")) > 1) so <- JoinLayers(so)
message("  ", ncol(so), " cells in object")

labels <- if (is.null(LABELS)) discover_labels() else intersect(LABELS, discover_labels())
if (length(labels) == 0) stop("No labels with aucell.csv under ", SCENIC_DIR, call. = FALSE)
message("  Labels with AUC: ", paste(labels, collapse = ", "))

# =============================================================================
# --- PER-LABEL: attach AUC, plot, summarize ----------------------------------
# =============================================================================
for (label in labels) {
  ldir <- file.path(SCENIC_DIR, label)
  message("\n########## ", label, " ##########")
  auc <- read_auc(ldir)
  if (is.null(auc)) { message("  [SKIP] empty/unreadable aucell.csv"); next }

  # --- cell-name alignment (THE gotcha): keep only shared barcodes ------------
  shared <- intersect(colnames(auc), colnames(so))
  if (length(shared) < 20) {
    message("  [SKIP] only ", length(shared), " AUC cells match the object ",
            "(barcode mismatch?). AUC cols e.g.: ",
            paste(utils::head(colnames(auc), 3), collapse = ", "))
    next
  }
  message("  ", nrow(auc), " regulons x ", length(shared), " shared cells")
  sub <- subset(so, cells = shared)
  auc <- auc[, colnames(sub), drop = FALSE]        # reorder to match the object

  # Attach as a new assay (regulons x cells). data slot = AUC (already scaled 0-1ish).
  sub[["AUC"]] <- CreateAssayObject(data = auc)
  DefaultAssay(sub) <- "AUC"

  out_dir <- file.path(ldir, "R_plots")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # --- locate the focus-TF regulon by prefix ("Nr4a1(+)") ---------------------
  tf_hit <- grep(paste0("^", FOCUS_TF, "\\("), rownames(sub), value = TRUE)
  if (length(tf_hit) == 0)
    message("  [note] no ", FOCUS_TF, " regulon recovered for this label.")

  # --- Plot 1: FeaturePlot of the TF regulon on the subset's UMAP -------------
  umap <- pick_umap(sub)
  if (length(tf_hit) && !is.na(umap)) {
    p1 <- FeaturePlot(sub, features = tf_hit[1], reduction = umap, order = TRUE) +
      ggtitle(paste0(tf_hit[1], " regulon AUC — ", label))
    ggsave(file.path(out_dir, paste0("FeaturePlot_", FOCUS_TF, ".png")), p1,
           width = 7, height = 6, dpi = DPI_SETTING)
  } else if (length(tf_hit) && is.na(umap)) {
    message("  [note] no UMAP reduction on the subset; skipping FeaturePlot.")
  }

  # --- Plot 2: violin of TF regulon AUC by condition --------------------------
  if (length(tf_hit) && CONDITION_COLUMN %in% colnames(sub@meta.data)) {
    p2 <- VlnPlot(sub, features = tf_hit[1], group.by = CONDITION_COLUMN,
                  pt.size = 0) + ggtitle(paste0(tf_hit[1], " — ", label)) +
      theme(legend.position = "none")
    ggsave(file.path(out_dir, paste0("Violin_", FOCUS_TF, "_by_condition.png")), p2,
           width = 7, height = 5, dpi = DPI_SETTING)
  }

  # --- Plot 3: mean-AUC heatmap (top-variable regulons x condition) -----------
  if (CONDITION_COLUMN %in% colnames(sub@meta.data)) {
    grp <- as.character(sub@meta.data[[CONDITION_COLUMN]])
    # mean AUC per regulon per condition
    mean_by <- t(apply(auc, 1, function(v) tapply(v, grp, mean, na.rm = TRUE)))
    vars    <- apply(mean_by, 1, function(v) stats::var(v, na.rm = TRUE))
    top     <- names(sort(vars, decreasing = TRUE))[seq_len(min(HEATMAP_TOP_N, nrow(mean_by)))]
    hm_df   <- as.data.frame(as.table(mean_by[top, , drop = FALSE]))
    colnames(hm_df) <- c("regulon", "condition", "meanAUC")
    p3 <- ggplot(hm_df, aes(condition, regulon, fill = meanAUC)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "grey92", high = "#08306B", name = "mean\nAUC") +
      labs(title = paste0("Mean regulon AUC — ", label),
           subtitle = paste0("top ", length(top), " most-variable regulons across ",
                             CONDITION_COLUMN),
           x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
    ggsave(file.path(out_dir, "Heatmap_meanAUC_topRegulons.png"), p3,
           width = 8, height = max(6, length(top) * 0.22 + 2),
           dpi = DPI_SETTING, limitsize = FALSE)
    utils::write.csv(as.data.frame(mean_by) |> tibble::rownames_to_column("regulon"),
                     file.path(out_dir, paste0(label, "_meanAUC_by_condition.csv")),
                     row.names = FALSE)
  }

  # --- per-cell TF AUC + metadata (for downstream stats WT/Polyp/KO) ----------
  if (length(tf_hit)) {
    md_cols <- intersect(c(CELLTYPE_COLUMN, CONDITION_COLUMN, "SampleID"),
                         colnames(sub@meta.data))
    per_cell <- data.frame(cell = colnames(sub),
                           AUC  = auc[tf_hit[1], ],
                           sub@meta.data[, md_cols, drop = FALSE],
                           check.names = FALSE, row.names = NULL)
    utils::write.csv(per_cell,
                     file.path(out_dir, paste0(label, "_", FOCUS_TF, "_AUC_percell.csv")),
                     row.names = FALSE)
  }
  rm(sub); gc()
}

message("\n=== Script 12c complete ===")
message("  Per-label plots + tables under: ", SCENIC_DIR, "/<label>/R_plots/")
message("  Binarization (optional): pyscenic binarize on the python side, or ",
        "AUCell::AUCell_exploreThresholds(assignCells=TRUE) in R, then attach as a ",
        "second 0/1 assay and plot identically.")
