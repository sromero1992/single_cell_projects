# =============================================================================
# decontx.R — DecontX ambient-RNA correction (list- or single-object aware)
# =============================================================================
# Optionally uses the raw/droplet matrix as the ambient `background`, which the
# DecontX authors recommend: the ambient RNA distribution is then estimated from
# ALL barcodes (empty droplets, cells lost to QC/flow) rather than only the
# heuristic clusters within the filtered cells.
# =============================================================================

#' Run DecontX ambient-RNA correction on a Seurat object or list of them
#'
#' @param x A Seurat object or a list of Seurat objects.
#' @param assay Assay to decontaminate (default "RNA").
#' @param round_counts If TRUE, round corrected counts to integers and drop
#'   structural zeros.
#' @param mt_pattern Mitochondrial regex pattern (default "^MT-|^mt-").
#' @param background Optional raw/droplet counts matrix to use as the ambient
#'   background (passed to `decontX(background=)`). If supplied it is used for
#'   every object; usually you instead set `background_dir` to auto-load per sample.
#' @param background_dir Directory holding one sub-folder per sample with a raw
#'   Cell Ranger h5 (`<dir>/<SampleID>/<background_h5>`). NULL = no background.
#' @param background_h5 Raw h5 file name (default "sample_raw_feature_bc_matrix.h5").
#' @param sample_col Metadata column giving each object's SampleID (default "SampleID").
#' @return `x` with decontaminated counts in `assay`, a `decontX_contamination`
#'   column, and refreshed `nCount_RNA`/`nFeature_RNA`/`percent.mt`.
#' @export
run_decontx <- function(x, assay = "RNA", round_counts = TRUE,
                        mt_pattern = "^MT-|^mt-",
                        background = NULL, background_dir = NULL,
                        background_h5 = "sample_raw_feature_bc_matrix.h5",
                        sample_col = "SampleID") {

  if (is.list(x) && !inherits(x, "Seurat")) {
    message("   -> Processing list of ", length(x), " sample(s) with DecontX...")
    return(lapply(x, run_decontx, assay = assay, round_counts = round_counts,
                  mt_pattern = mt_pattern, background = background,
                  background_dir = background_dir, background_h5 = background_h5,
                  sample_col = sample_col))
  }
  if (!inherits(x, "Seurat"))
    stop("Input 'x' must be a Seurat object or a list of Seurat objects.")
  if (!requireNamespace("celda", quietly = TRUE))
    stop("celda (DecontX) is not installed. BiocManager::install('celda').")

  sid <- if (sample_col %in% colnames(x@meta.data))
    unique(as.character(x@meta.data[[sample_col]]))[1]
  else if ("orig.ident" %in% colnames(x@meta.data))
    unique(as.character(x$orig.ident))[1] else "Sample"
  message(paste("    -> Running DecontX on", sid, "..."))

  # ---- resolve the raw/droplet background matrix ----------------------------
  bg <- background
  if (is.null(bg) && !is.null(background_dir)) {
    raw_path <- file.path(background_dir, sid, background_h5)
    if (file.exists(raw_path)) {
      message("      loading raw/droplet background: ", raw_path)
      bg <- tryCatch({
        r <- Seurat::Read10X_h5(raw_path); if (is.list(r)) r[[1]] else r
      }, error = function(e) {
        message("      [WARN] could not read raw h5: ", e$message); NULL })
    } else {
      message("      [NOTE] raw background not found for ", sid,
              "; DecontX will run without a background.")
    }
  }

  tryCatch({
    counts_sparse <- SeuratObject::GetAssayData(x, assay = assay, layer = "counts")

    if (!is.null(bg)) {
      # DecontX needs the background genes to match. Use the intersection; any
      # custom feature absent from the raw matrix (e.g. a summed probe) keeps its
      # original counts and is stitched back afterwards.
      common <- intersect(rownames(counts_sparse), rownames(bg))
      if (length(common) < 2) stop("background shares too few genes with the sample")
      res <- celda::decontX(x = counts_sparse[common, , drop = FALSE],
                            background = bg[common, , drop = FALSE])
      dec <- res$decontXcounts
      extra <- setdiff(rownames(counts_sparse), common)
      corrected_counts <- if (length(extra) == 0) dec
        else rbind(dec, counts_sparse[extra, , drop = FALSE])[rownames(counts_sparse), ]
    } else {
      res <- celda::decontX(x = counts_sparse)
      corrected_counts <- res$decontXcounts
    }

    x$decontX_contamination <- res$estimates$all$contamination

    if (round_counts) corrected_counts <- round(corrected_counts)
    corrected_counts <- Matrix::drop0(corrected_counts)
    SeuratObject::LayerData(x, assay = assay, layer = "counts") <- corrected_counts

    x$nCount_RNA   <- Matrix::colSums(SeuratObject::GetAssayData(x, assay = assay, layer = "counts"))
    x$nFeature_RNA <- Matrix::colSums(SeuratObject::GetAssayData(x, assay = assay, layer = "counts") > 0)
    x$percent.mt   <- Seurat::PercentageFeatureSet(x, pattern = mt_pattern)

    message(paste("    -> DecontX complete for", sid,
                  if (!is.null(bg)) "(with raw background)." else "(no background)."))
    x
  }, error = function(e) {
    message(paste("    -> [WARNING] DecontX failed for", sid, "| Error:", e$message,
                  "| Proceeding with uncorrected counts."))
    x
  })
}
