#' Apply QC Filtering Bounds
#'
#' @param x A Seurat object or a list of Seurat objects.
#' @param p A list of QC parameters (e.g., qc$pre or qc$post).
#' @param apply If TRUE, physically subset cells; if FALSE, only write metadata flag.
#' @export
apply_qc <- function(x, p, apply = TRUE) {
  
  # 1. Map recursively if input is a sample list
  if (.is_sample_list(x)) {
    return(.map_samples(x, apply_qc, p = p, apply = apply))
  }
  
  stopifnot(.is_seurat(x))
  
  # Set default fallbacks if optional fields are missing
  p$min_features       <- p$min_features       %||% 200
  p$max_features       <- p$max_features       %||% Inf
  p$min_counts         <- p$min_counts         %||% 0
  p$max_counts         <- p$max_counts         %||% Inf
  p$max_mt             <- p$max_mt             %||% 20
  p$min_cells_per_gene <- p$min_cells_per_gene %||% 0
  p$mt_pattern         <- p$mt_pattern         %||% "^mt-"
  
  # 2. Compute % MT if not present
  s <- TamuScDSC_schema()
  if (!s$percent_mt %in% colnames(x@meta.data)) {
    x[[s$percent_mt]] <- Seurat::PercentageFeatureSet(x, pattern = p$mt_pattern)
  }
  
  # 3. Evaluate Pass/Fail Condition
  pass <- x$nFeature_RNA >= p$min_features &
    x$nFeature_RNA <= p$max_features &
    x$nCount_RNA   >= p$min_counts   &
    x$nCount_RNA   <= p$max_counts   &
    x@meta.data[[s$percent_mt]] <= p$max_mt
  
  x@meta.data[["qc_pass"]] <- pass
  
  # 4. Subset if requested
  if (apply) {
    n_before <- ncol(x)
    x <- subset(x, cells = colnames(x)[pass])
    
    # Filter low-expression genes across dataset if specified
    if (!is.null(p$min_cells_per_gene) && p$min_cells_per_gene > 0) {
      # A merged Seurat v5 object keeps one counts layer per sample, and
      # GetAssayData() cannot collapse multiple layers. Join them first so the
      # gene-detection count is computed across the whole dataset.
      if (inherits(x[["RNA"]], "Assay5") &&
          length(SeuratObject::Layers(x[["RNA"]], search = "counts")) > 1) {
        x <- SeuratObject::JoinLayers(x)
      }
      counts <- SeuratObject::GetAssayData(x, assay = "RNA", layer = "counts")
      gene_cell_counts <- Matrix::rowSums(counts > 0)
      keep_genes <- names(gene_cell_counts[gene_cell_counts >= p$min_cells_per_gene])
      x <- x[keep_genes, ]
    }
    
    message(sprintf("   [QC] Kept %d cells (removed %d).", ncol(x), n_before - ncol(x)))
  }
  
  return(x)
}