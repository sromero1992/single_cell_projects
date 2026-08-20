# =============================================================================
# schema.R — the metadata-column contract shared by every step
# =============================================================================
# Steps never assume what ran before them; they read and write these agreed
# column names. Keeping them in one place means a downstream step (or a re-run)
# can look for "doublet_consensus" without caring which recipe produced it.
# =============================================================================

#' TamuScDSC metadata schema
#'
#' The canonical metadata column names written by TamuScDSC steps. Use these
#' constants rather than hard-coding strings, so a rename happens in one place.
#'
#' @return A named list of column-name constants.
#' @export
TamuScDSC_schema <- function() {
  list(
    sample        = "SampleID",
    # QC
    percent_mt        = "percent_mt",
    qc_pass_light     = "qc_pass_light",
    qc_pass_stringent = "qc_pass_stringent",
    # SCEVAN
    scevan_class  = "scevan_class",
    # Doublets (names chosen to match 01b_cluster_doublet_review.R verbatim)
    df_score      = "DF_score",
    df_class      = "DF_class",
    sc_score      = "scDblFinder_score",
    sc_class      = "scDblFinder_class",
    consensus     = "doublet_consensus",
    # Cluster review
    cluster_flag  = "cluster_dbl_flag"
  )
}

# Internal convenience: which of the schema columns are present on an object.
.schema_present <- function(obj) {
  s <- unlist(TamuScDSC_schema())
  cols <- colnames(obj@meta.data)
  s[s %in% cols]
}
