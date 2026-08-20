#' TamuScDSC: modular single-cell preprocessing, QC and doublet detection
#'
#' Atomic, order-free, re-runnable building blocks for the preprocessing stage
#' of a single-cell RNA-seq pipeline. See `vignette("three-scenarios")` and the
#' package `ARCHITECTURE.md` for the design.
#'
#' @keywords internal
"_PACKAGE"

# Quiet R CMD check about tidy-eval / NSE column names used in dplyr/ggplot2.
utils::globalVariables(c(
  ".data", "cluster", "N_Cells", "DF_mean", "scDbl_mean",
  "DF_frac_dbl", "scDbl_frac_dbl", "Method", "Mean_Score",
  "Freq", "DoubletFinder", "scDblFinder",
  # de.R (Enrichr / DE / DV / plotting NSE columns)
  "database", "Adjusted.P.value", "Overlap", "Genes", "Term", "group_key",
  "gene_count", "neg_log10_padj", "score", "mean_score", "mean_gene_count",
  "n_cell_types", "cell_type", "contrast", "direction", "gene", "Term_wrapped",
  "avg_log2FC", "p_val_adj", "member", "gene_vec", "fill_val", "hub_score",
  "n_pathways", "genes", "Gene", "Pval", "PValue", "pvalue", "pval"))
