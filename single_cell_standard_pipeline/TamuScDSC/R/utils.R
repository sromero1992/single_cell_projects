# =============================================================================
# utils.R — small shared helpers
# =============================================================================

#' NULL-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Is this a Seurat object?
#' @keywords internal
#' @noRd
.is_seurat <- function(x) methods::is(x, "Seurat")

#' Is this a list of Seurat objects?
#' @keywords internal
#' @noRd
.is_sample_list <- function(x) {
  is.list(x) && !.is_seurat(x) && length(x) > 0 &&
    all(vapply(x, .is_seurat, logical(1)))
}

#' Stamp a provenance record into obj@misc$TamuScDSC_provenance
#'
#' Appends one record (step name, timestamp, parameters) so that after any
#' sequence of steps the full history is inspectable with provenance().
#' @keywords internal
#' @noRd
.stamp <- function(obj, step, params = list()) {
  if (!.is_seurat(obj)) return(obj)
  rec <- list(step = step, time = as.character(Sys.time()),
              n_cells = ncol(obj), params = params)
  obj@misc$TamuScDSC_provenance <- c(obj@misc$TamuScDSC_provenance %||% list(), list(rec))
  obj
}

#' Map a step over a sample list while preserving names and provenance
#' @keywords internal
#' @noRd
.map_samples <- function(x, fun, ...) {
  out <- lapply(x, fun, ...)
  names(out) <- names(x)
  out
}

#' Print the TamuScDSC provenance log of an object
#'
#' @param obj A Seurat object processed by TamuScDSC.
#' @return Invisibly, the provenance list; prints a compact summary.
#' @export
provenance <- function(obj) {
  p <- if (.is_seurat(obj)) obj@misc$TamuScDSC_provenance else NULL
  if (is.null(p) || length(p) == 0) {
    message("No TamuScDSC provenance recorded on this object.")
    return(invisible(NULL))
  }
  for (i in seq_along(p)) {
    r <- p[[i]]
    message(sprintf("  %2d. %-22s %s  (%d cells)",
                    i, r$step, r$time, r$n_cells %||% NA))
  }
  invisible(p)
}
