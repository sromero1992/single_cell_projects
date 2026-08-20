# =============================================================================
# ingest.R — normalize all three input scenarios to one sample list
# =============================================================================
# as_sample_list() is the single boundary that decouples "where the data came
# from" from "what we do to it". Everything downstream operates on a named list
# of per-sample Seurat objects; merge_samples() collapses it back when ready.
# =============================================================================

#' Normalize any supported input into a named list of per-sample Seurat objects
#'
#' Handles the three scenarios uniformly:
#' * **character** — file paths (or a directory) of per-sample `.rds` checkpoints.
#' * **Seurat** — one (possibly merged) object, split by `sample_col`.
#' * **list** — an already-loaded list of Seurat objects.
#'
#' @param x Character paths/dir, a Seurat object, or a list of Seurat objects.
#' @param sample_col Metadata column identifying the sample (default "SampleID").
#' @param ... Passed to methods.
#' @return A named list of Seurat objects, one per sample.
#' @export
as_sample_list <- function(x, sample_col = "SampleID", ...) UseMethod("as_sample_list")

#' @export
as_sample_list.default <- function(x, sample_col = "SampleID", ...)
  stop("as_sample_list(): unsupported input of class '", class(x)[1],
       "'. Expected file paths, a Seurat object, or a list of Seurat objects.")

#' @export
as_sample_list.character <- function(x, sample_col = "SampleID", pattern = "\\.rds$", ...) {
  # A single existing directory -> list its .rds files; otherwise treat as paths.
  if (length(x) == 1 && dir.exists(x)) {
    files <- list.files(x, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  } else {
    files <- x
  }
  if (length(files) == 0) stop("as_sample_list(): no .rds files found.")
  missing <- files[!file.exists(files)]
  if (length(missing) > 0)
    stop("as_sample_list(): missing files:\n  ", paste(missing, collapse = "\n  "))

  message(sprintf("  as_sample_list: reading %d checkpoint file(s)...", length(files)))
  out <- lapply(files, readRDS)
  nm  <- vapply(seq_along(out), function(i) {
    seu <- out[[i]]
    sid <- if (sample_col %in% colnames(seu@meta.data))
      unique(as.character(seu@meta.data[[sample_col]]))[1] else NA
    sid %||% tools::file_path_sans_ext(basename(files[i]))
  }, character(1))
  nm[is.na(nm)] <- tools::file_path_sans_ext(basename(files[is.na(nm)]))
  stats::setNames(out, make.unique(nm))
}

#' @export
as_sample_list.Seurat <- function(x, sample_col = "SampleID", ...) {
  if (!sample_col %in% colnames(x@meta.data)) {
    message(sprintf("  as_sample_list: '%s' not found; treating the object as a single sample.",
                    sample_col))
    return(stats::setNames(list(x), "all"))
  }
  # JoinLayers first so the split produces clean per-sample objects (Seurat v5).
  if (inherits(x[["RNA"]], "Assay5") &&
      length(SeuratObject::Layers(x[["RNA"]], search = "counts")) > 1) {
    x <- SeuratObject::JoinLayers(x)
  }
  Seurat::SplitObject(x, split.by = sample_col)
}

#' @export
as_sample_list.list <- function(x, sample_col = "SampleID", ...) {
  if (!.is_sample_list(x))
    stop("as_sample_list(): list elements must all be Seurat objects.")
  if (is.null(names(x)) || any(names(x) == "")) {
    nm <- vapply(seq_along(x), function(i) {
      seu <- x[[i]]
      if (sample_col %in% colnames(seu@meta.data))
        unique(as.character(seu@meta.data[[sample_col]]))[1] else paste0("sample", i)
    }, character(1))
    x <- stats::setNames(x, make.unique(nm))
  }
  x
}

#' Read per-sample checkpoints from a directory (thin wrapper)
#'
#' @param dir Directory containing per-sample `.rds` files.
#' @param pattern File-name pattern.
#' @param sample_col Sample metadata column.
#' @return A named sample list.
#' @export
read_checkpoints <- function(dir, pattern = "\\.rds$", sample_col = "SampleID")
  as_sample_list(dir, sample_col = sample_col, pattern = pattern)

#' Load per-sample 10x/H5 matrices into a sample list (metadata-driven)
#'
#' Migrated from the STEP 2.1a loader in `01_process_data.R`. Loops over the rows
#' of a metadata table, reads each sample's Cell Ranger H5, builds a Seurat
#' object, and broadcasts that sample's metadata row onto its cells (skipping
#' names reserved by Seurat / the pipeline). Optional probe summing for Flex/KO
#' designs is supported via `probe`.
#'
#' @param metadata A data.frame with one row per sample and a `sample_col`.
#' @param h5_dir Directory containing one sub-folder per sample.
#' @param sample_col Column naming each sample (default "SampleID").
#' @param h5_name File name of the filtered matrix inside each sample folder.
#' @param min_cells Passed to [Seurat::CreateSeuratObject()].
#' @param probe Optional list to enable probe summing: `list(h5_name=, feature=,
#'   probes_for_sum=, mapping=)`. NULL (default) skips all probe handling.
#' @return A named list of per-sample Seurat objects.
#' @export
read_10x_samples <- function(metadata, h5_dir, sample_col = "SampleID",
                             h5_name = "sample_filtered_feature_bc_matrix.h5",
                             min_cells = 5, probe = NULL) {
  stopifnot(is.data.frame(metadata))
  colnames(metadata) <- clean_header(colnames(metadata))
  key <- resolve_sample_col(metadata, c(sample_col, "SampleID", "Sample.ID", "Sample"))
  if (is.null(key)) stop("read_10x_samples(): no sample column found in metadata.")

  reserved <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent_mt",
                "Doublet_Status", "DF_score", "DF_class", "scDblFinder_score",
                "scDblFinder_class", "seurat_clusters", "CellType", "sub_cell_types")

  out <- list()
  for (i in seq_len(nrow(metadata))) {
    row <- metadata[i, , drop = FALSE]
    sid <- as.character(row[[key]])
    h5  <- file.path(h5_dir, sid, h5_name)
    if (!file.exists(h5)) { warning("read_10x_samples(): H5 not found for ", sid); next }

    counts <- Seurat::Read10X_h5(h5)
    if (is.list(counts)) counts <- counts[[1]]           # gene-expression matrix

    # Optional probe summing (Flex / KO designs).
    if (!is.null(probe)) counts <- .add_probe_feature(counts, h5_dir, sid, probe)

    seu <- Seurat::CreateSeuratObject(counts = counts, project = sid, min.cells = min_cells)
    seu[[sample_col]] <- sid
    for (cn in setdiff(colnames(row), c(key, reserved))) seu[[cn]] <- row[[cn]][1]
    out[[sid]] <- seu
    rm(counts); gc()
  }
  if (length(out) == 0) stop("read_10x_samples(): no samples loaded.")
  message(sprintf("  read_10x_samples: loaded %d sample(s).", length(out)))
  out
}

# Internal: append a summed custom probe feature (e.g. Nr4a1_cust) to a matrix.
.add_probe_feature <- function(counts, h5_dir, sid, probe) {
  ph5 <- file.path(h5_dir, sid, probe$h5_name %||% "sample_raw_probe_bc_matrix.h5")
  if (!file.exists(ph5)) { warning("  probe H5 not found for ", sid); return(counts) }
  pm <- Seurat::Read10X_h5(ph5)
  if (is.list(pm) && "Probe Barcode" %in% names(pm)) pm <- pm[["Probe Barcode"]]
  common <- intersect(colnames(counts), colnames(pm))
  found  <- intersect(probe$probes_for_sum %||% character(0), rownames(pm))
  vec <- stats::setNames(numeric(length(common)), common)
  if (length(found) > 0) {
    s <- Matrix::colSums(pm[found, common, drop = FALSE]); vec[names(s)] <- s
  }
  crow <- Matrix::Matrix(0, 1, ncol(counts), sparse = TRUE)
  colnames(crow) <- colnames(counts)
  rownames(crow) <- probe$feature %||% "custom_probe"
  crow[1, names(vec)] <- vec
  rbind(counts, crow)
}

#' Merge a sample list into one object
#'
#' @param x A sample list (or, passed through, a Seurat object).
#' @param add_cell_ids Prefix cell names with the sample name to guarantee
#'   unique barcodes (default TRUE).
#' @param join Join the per-sample counts layers into one after merging
#'   (default TRUE). Seurat v5 keeps one layer per sample otherwise, which blocks
#'   any step that needs a single counts matrix (e.g. the global gene filter).
#'   Safe for both integration routes: RunHarmony wants the joined matrix, and
#'   IntegrateLayers re-splits internally. Set FALSE only if you specifically
#'   need the split layers preserved on the returned object.
#' @param project Project name for the merged object.
#' @return A single Seurat object (layers joined unless `join = FALSE`).
#' @export
merge_samples <- function(x, add_cell_ids = TRUE, join = TRUE, project = "TamuScDSC") {
  if (.is_seurat(x)) return(x)
  if (!.is_sample_list(x)) stop("merge_samples(): expected a list of Seurat objects.")
  if (length(x) == 1) return(x[[1]])
  ids <- if (add_cell_ids) names(x) else NULL
  message(sprintf("  merge_samples: merging %d samples...", length(x)))
  merged <- merge(x[[1]], y = x[-1], add.cell.ids = ids, project = project)
  if (join && inherits(merged[["RNA"]], "Assay5") &&
      length(SeuratObject::Layers(merged[["RNA"]], search = "counts")) > 1) {
    message("  merge_samples: joining per-sample layers into one...")
    merged <- SeuratObject::JoinLayers(merged)
  }
  .stamp(merged, "merge_samples", list(n_samples = length(x), joined = join))
}
