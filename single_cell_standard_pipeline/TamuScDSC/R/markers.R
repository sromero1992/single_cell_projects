# =============================================================================
# markers.R — extensible marker registry
# =============================================================================
# The default marker set ships as inst/extdata/cell_type_markers.csv. Edit it
# through the API (add_markers / load_markers), not by touching annotation code.
# Scripts 02-05 should call get_markers() instead of carrying their own copy.
#
# The shipped CSV is the lab's own wide format:
#   tier, cell_type, parent_cell_type, markers, subcluster_resolution, description
# where `markers` is a pipe-delimited gene list (e.g. "Jchain|Mzb1|Igha").
# get_markers() also accepts a simple long format (cell_type, gene[, weight]).
# Genes are exploded to one row per gene internally so markers_as_list() and the
# annotation scripts get a clean cell_type -> gene mapping.
# =============================================================================

#' Path to the built-in marker CSV shipped with the package
#' @keywords internal
#' @noRd
.markers_default_path <- function() {
  p <- system.file("extdata", "cell_type_markers.csv", package = "TamuScDSC")
  if (identical(p, "")) NA_character_ else p
}

#' Get the marker table (built-in default, or extended)
#'
#' Returns one row per (cell_type, gene). Pipe-delimited gene lists in a
#' `markers` column are exploded automatically. `tier` and `parent_cell_type`
#' are carried through when present.
#'
#' @param path Optional CSV to read instead of the built-in default.
#' @return A data.frame with columns cell_type, gene, weight, and (if present)
#'   tier and parent_cell_type.
#' @export
get_markers <- function(path = NULL) {
  path <- path %||% .markers_default_path()
  if (is.na(path) || !file.exists(path)) {
    warning("get_markers(): marker CSV not found; returning empty table.")
    return(data.frame(cell_type = character(), gene = character(),
                      weight = numeric(), stringsAsFactors = FALSE))
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- clean_header(colnames(df))
  ct <- resolve_sample_col(df, c("cell_type", "celltype", "CellType", "type"))
  gn <- resolve_sample_col(df, c("markers", "gene", "Gene", "symbol", "marker"))
  if (is.null(ct) || is.null(gn))
    stop("get_markers(): CSV must have a cell-type column and a gene/markers column.")

  tier_col   <- resolve_sample_col(df, c("tier"))
  parent_col <- resolve_sample_col(df, c("parent_cell_type", "parent"))
  wt_col     <- resolve_sample_col(df, c("weight", "Weight", "w"))

  rows <- lapply(seq_len(nrow(df)), function(i) {
    genes <- unlist(strsplit(as.character(df[[gn]][i]), "\\|"))
    genes <- trimws(genes); genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) return(NULL)
    data.frame(
      tier      = if (!is.null(tier_col))   as.character(df[[tier_col]][i])   else NA_character_,
      cell_type = as.character(df[[ct]][i]),
      parent_cell_type = if (!is.null(parent_col)) as.character(df[[parent_col]][i]) else NA_character_,
      gene      = genes,
      weight    = if (!is.null(wt_col)) suppressWarnings(as.numeric(df[[wt_col]][i])) else 1,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$weight[is.na(out$weight)] <- 1
  out[!duplicated(out[, c("cell_type", "gene")]), , drop = FALSE]
}

#' Add markers for a cell type (returns an extended table)
#'
#' @param markers A marker data.frame from [get_markers()].
#' @param cell_type Character, the cell type to add genes to.
#' @param genes Character vector of gene symbols.
#' @param weight Weight applied to the added genes (default 1).
#' @return The extended marker data.frame (de-duplicated).
#' @export
add_markers <- function(markers, cell_type, genes, weight = 1) {
  stopifnot(is.data.frame(markers))
  add <- data.frame(cell_type = cell_type, gene = as.character(genes),
                    weight = weight, stringsAsFactors = FALSE)
  # Align to any extra columns the table carries (tier, parent_cell_type).
  for (extra in setdiff(colnames(markers), colnames(add))) add[[extra]] <- NA
  out <- rbind(markers, add[, colnames(markers), drop = FALSE])
  out[!duplicated(out[, c("cell_type", "gene")]), , drop = FALSE]
}

#' Load (and optionally merge) markers from a CSV
#'
#' @param path CSV path in the same long format as the built-in file.
#' @param base Optional existing marker table to merge onto.
#' @return A marker data.frame.
#' @export
load_markers <- function(path, base = NULL) {
  new <- get_markers(path)
  if (is.null(base)) return(new)
  common <- intersect(colnames(base), colnames(new))
  out <- rbind(base[, common, drop = FALSE], new[, common, drop = FALSE])
  out[!duplicated(out[, c("cell_type", "gene")]), , drop = FALSE]
}

#' Convert a marker table to a named list (cell_type -> gene vector)
#'
#' @param markers A marker data.frame.
#' @return A named list of character vectors, suitable for AddModuleScore etc.
#'   Cell types keep their first-appearance order (important for dotplot axes).
#' @export
markers_as_list <- function(markers) {
  ct <- factor(markers$cell_type, levels = unique(markers$cell_type))
  split(markers$gene, ct)
}
