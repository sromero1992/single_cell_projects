# =============================================================================
# metadata.R — tolerant metadata cleaning + broadcast onto cells
# =============================================================================
# clean_header() and norm_key() are ported from 01_process_data_new.R. They make
# the SampleID join survive openxlsx dot-conversion, stray XML fragments, and
# casing/spacing differences between the metadata sheet and the object.
# =============================================================================

#' Strip XML fragments / whitespace artifacts from a header string
#' @param x Character vector of column headers.
#' @return Cleaned character vector.
#' @export
clean_header <- function(x) {
  x <- as.character(x)
  x <- gsub("<[^>]*>", "", x)          # drop XML/HTML fragments
  x <- gsub("[\r\n\t]", " ", x)
  x <- trimws(gsub("\\s+", " ", x))
  x
}

#' Normalize a key for tolerant matching (lowercase, alphanumeric only)
#' @param x Character vector.
#' @return Normalized character vector.
#' @export
norm_key <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(as.character(x))))

#' Resolve the sample-id column of an object or data.frame, tolerantly
#'
#' @param obj A Seurat object or data.frame.
#' @param candidates Preferred column names to try, in order.
#' @return The resolved column name, or NULL if none match.
#' @export
resolve_sample_col <- function(obj, candidates = c("SampleID", "Sample.ID",
                                                   "Sample_ID", "Sample", "orig.ident")) {
  cols <- if (.is_seurat(obj)) colnames(obj@meta.data) else colnames(obj)
  # exact first
  hit <- intersect(candidates, cols)
  if (length(hit) > 0) return(hit[1])
  # normalized fallback
  nk <- norm_key(cols)
  for (cand in candidates) {
    j <- which(nk == norm_key(cand))
    if (length(j) > 0) return(cols[j[1]])
  }
  NULL
}

#' Broadcast per-sample metadata rows onto cells
#'
#' Works on a sample list or a merged object. Joins the metadata sheet to the
#' object by a tolerantly-resolved sample column, and writes each metadata
#' field as a per-cell column.
#'
#' @param x A Seurat object or a sample list.
#' @param metadata A data.frame with one row per sample.
#' @param sample_col Column in `metadata` (and target) identifying the sample.
#' @param overwrite Overwrite existing columns of the same name (default FALSE).
#' @return `x` with the metadata columns attached.
#' @export
attach_metadata <- function(x, metadata, sample_col = "SampleID", overwrite = FALSE) {
  if (.is_sample_list(x))
    return(.map_samples(x, attach_metadata, metadata = metadata,
                        sample_col = sample_col, overwrite = overwrite))
  stopifnot(.is_seurat(x), is.data.frame(metadata))

  colnames(metadata) <- clean_header(colnames(metadata))
  meta_key <- resolve_sample_col(metadata, c(sample_col, "SampleID", "Sample.ID",
                                             "Sample_ID", "Sample"))
  obj_key  <- resolve_sample_col(x,        c(sample_col, "SampleID", "Sample.ID",
                                             "Sample_ID", "Sample", "orig.ident"))
  if (is.null(meta_key)) stop("attach_metadata(): no sample column found in metadata.")
  if (is.null(obj_key))  stop("attach_metadata(): no sample column found on the object.")

  # Build normalized lookup: sample-key -> metadata row.
  m <- metadata
  m$.k <- norm_key(m[[meta_key]])
  cell_k <- norm_key(as.character(x@meta.data[[obj_key]]))
  idx <- match(cell_k, m$.k)
  unmatched <- unique(as.character(x@meta.data[[obj_key]])[is.na(idx)])
  if (length(unmatched) > 0)
    warning("attach_metadata(): no metadata row for sample(s): ",
            paste(unmatched, collapse = ", "))

  add_cols <- setdiff(colnames(m), c(meta_key, ".k"))
  for (cn in add_cols) {
    if (cn %in% colnames(x@meta.data) && !overwrite) next
    x@meta.data[[cn]] <- m[[cn]][idx]
  }
  .stamp(x, "attach_metadata", list(fields = add_cols))
}
