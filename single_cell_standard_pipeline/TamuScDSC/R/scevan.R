# =============================================================================
# scevan.R — SCEVAN copy-number / malignant-cell inference (standalone, optional)
# =============================================================================
# Separate, optional step. Call run_scevan() on a per-sample object or a sample
# list; results land in the `scevan_class` metadata column. Runs on raw counts.
# =============================================================================

#' Run SCEVAN on a per-sample object (or sample list)
#'
#' @param x A Seurat object or a sample list.
#' @param sample_id Label used for SCEVAN's output files.
#' @param par_cores Number of processor cores for SCEVAN (its `par_cores`).
#' @param organism "mouse" or "human".
#' @param subclones Passed to `pipelineCNA(SUBCLONES=)`.
#' @param plot_tree Passed to `pipelineCNA(plotTree=)`.
#' @param out_dir Directory SCEVAN writes its output into (per sample). SCEVAN
#'   writes an `output/` folder into the working directory, so we run it with the
#'   wd set here and restore it afterwards.
#' @param assay Assay providing the counts matrix (default "RNA").
#' @return `x` with a `scevan_class` column ("tumor"/"normal"/"filtered"/NA).
#' @export
run_scevan <- function(x, sample_id = NULL, par_cores = 1, organism = "mouse",
                       subclones = FALSE, plot_tree = FALSE, out_dir = ".",
                       assay = "RNA") {
  if (.is_sample_list(x)) {
    return(.map_samples(x, function(seu) {
      sid <- unique(as.character(seu@meta.data[[TamuScDSC_schema()$sample]]))[1] %||% "sample"
      run_scevan(seu, sample_id = sid, par_cores = par_cores, organism = organism,
                 subclones = subclones, plot_tree = plot_tree,
                 out_dir = file.path(out_dir, sid), assay = assay)
    }))
  }
  stopifnot(.is_seurat(x))
  if (!requireNamespace("SCEVAN", quietly = TRUE))
    stop("SCEVAN is not installed. install_github('AntonioDeFalco/SCEVAN').")

  s   <- TamuScDSC_schema()
  sid <- sample_id %||% "sample"
  cnts <- as.matrix(SeuratObject::GetAssayData(x, assay = assay, layer = "counts"))

  # SCEVAN scatters its output into the working directory and does not fully
  # honour output_dir, so run it from out_dir and restore the wd no matter what.
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  res <- tryCatch(
    SCEVAN::pipelineCNA(cnts, sample = sid, par_cores = par_cores,
                        organism = organism, SUBCLONES = subclones,
                        plotTree = plot_tree),
    error = function(e) { message("  [WARNING] SCEVAN failed on ", sid, ": ", e$message); NULL })

  cls <- stats::setNames(rep(NA_character_, ncol(x)), colnames(x))
  if (!is.null(res) && "class" %in% colnames(res)) {
    common <- intersect(rownames(res), colnames(x))
    cls[common] <- as.character(res[common, "class"])
  }
  x@meta.data[[s$scevan_class]] <- cls
  .stamp(x, "run_scevan", list(sample = sid, par_cores = par_cores, organism = organism))
}
