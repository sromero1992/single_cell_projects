# =============================================================================
# subset.R — extract a sample, and graft recomputed columns back by barcode
# =============================================================================
# The "redo a step on one sample and put it back" workflow. Every join is keyed
# on cell barcode, never on row position, so a re-run cannot silently misalign.
# =============================================================================

#' Extract one or more samples from a merged object (or a sample list)
#'
#' @param x A Seurat object or a sample list.
#' @param ids Character vector of sample ids to keep.
#' @param sample_col Sample metadata column (default "SampleID").
#' @return The same type as `x`, restricted to `ids`.
#' @export
subset_samples <- function(x, ids, sample_col = "SampleID") {
  if (.is_sample_list(x)) {
    keep <- intersect(names(x), ids)
    if (length(keep) == 0) stop("subset_samples(): no matching sample names.")
    return(x[keep])
  }
  stopifnot(.is_seurat(x))
  col <- resolve_sample_col(x, c(sample_col, "SampleID", "orig.ident"))
  if (is.null(col)) stop("subset_samples(): no sample column on the object.")
  cells <- colnames(x)[as.character(x@meta.data[[col]]) %in% ids]
  if (length(cells) == 0) stop("subset_samples(): no cells match ids: ",
                               paste(ids, collapse = ", "))
  subset(x, cells = cells)
}

#' Graft metadata columns from a child object back into a parent, by barcode
#'
#' Use after re-running a step on a subset: writes the recomputed columns onto
#' the matching cells of the parent, leaving all other cells untouched.
#'
#' @param parent The full Seurat object to update.
#' @param child A subset/re-run object carrying the recomputed columns.
#' @param cols Character vector of metadata columns to copy.
#' @return `parent` with the columns updated for the child's cells.
#' @export
graft_meta <- function(parent, child, cols) {
  stopifnot(.is_seurat(parent), .is_seurat(child))
  missing <- setdiff(cols, colnames(child@meta.data))
  if (length(missing) > 0)
    stop("graft_meta(): child is missing column(s): ", paste(missing, collapse = ", "))

  child_cells <- colnames(child)
  idx <- match(child_cells, colnames(parent))
  if (anyNA(idx))
    warning("graft_meta(): ", sum(is.na(idx)),
            " child cell(s) not present in parent; skipped.")
  keep <- !is.na(idx)

  for (cn in cols) {
    if (!cn %in% colnames(parent@meta.data)) {
      # initialise the column with a compatible NA type
      template <- child@meta.data[[cn]]
      parent@meta.data[[cn]] <- if (is.numeric(template)) NA_real_
                                else if (is.factor(template)) NA_character_
                                else NA
    }
    parent@meta.data[[cn]][idx[keep]] <- child@meta.data[[cn]][keep]
  }
  .stamp(parent, "graft_meta", list(cols = cols, n_cells = sum(keep)))
}
