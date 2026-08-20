# =============================================================================
# pipeline.R — declarative runner for arbitrary step sequences
# =============================================================================
# When you want ordering freedom *and* an auditable record, describe the run as
# a list of steps and let run_pipeline() execute them, threading the object and
# stamping provenance. This is optional sugar over plain piping.
# =============================================================================

#' Run an ordered list of steps over an object
#'
#' Each step is either a function (called as `fn(x)`) or a list
#' `list(fn = <function>, args = <named list>)` called as
#' `do.call(fn, c(list(x), args))`. The object is threaded through in order.
#'
#' @param x The starting object (sample list or Seurat object).
#' @param steps A list of steps (see above).
#' @param verbose Announce each step (default TRUE).
#' @return The final object, with provenance recorded.
#' @examples
#' \dontrun{
#' data <- run_pipeline(paths, list(
#'   as_sample_list,
#'   function(x) apply_qc(x, p = list(min_features = 500), apply = FALSE),
#'   list(fn = detect_doublets, args = list(method = "both")),
#'   merge_samples,
#'   list(fn = integrate_data, args = list(method = "RunHarmony"))
#' ))
#' }
#' @export
run_pipeline <- function(x, steps, verbose = TRUE) {
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    if (is.function(step)) {
      fn <- step; args <- list()
    } else if (is.list(step) && is.function(step$fn)) {
      fn <- step$fn; args <- step$args %||% list()
    } else {
      stop("run_pipeline(): step ", i,
           " must be a function or list(fn=, args=).")
    }
    label <- names(steps)[i]
    if (is.null(label) || label == "") label <- paste0("step_", i)
    if (verbose) message(sprintf("[run_pipeline] %d/%d: %s", i, length(steps), label))
    x <- do.call(fn, c(list(x), args))
  }
  x
}
