# =============================================================================
# process.R — reusable Seurat processing helpers
# =============================================================================
# These assume the data are already normalized (they do not call NormalizeData),
# matching the lab's existing usage. Use process_rna() on a whole object,
# process_and_extract_cell_types() to pull one or more cell types from any
# metadata column and re-embed them, and clustering() to re-cluster cheaply on
# an existing graph (no UMAP / neighbours recomputed).
# =============================================================================

#' Standard RNA processing: HVG -> scale -> PCA -> UMAP -> neighbours -> clusters
#'
#' Assumes the data are already normalized.
#'
#' @param data A Seurat object.
#' @param assay_name Assay to use (default "RNA").
#' @param num_hvg Number of variable features (default 2000).
#' @param dims_pca PCA dims used for UMAP and neighbours (default 50).
#' @param resolution Clustering resolution (default 1.0).
#' @return The processed Seurat object.
#' @export
process_rna <- function(data, assay_name = "RNA", num_hvg = 2000,
                        dims_pca = 50, resolution = 1.0) {
  Seurat::DefaultAssay(data) <- assay_name
  data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = num_hvg)
  data <- Seurat::ScaleData(data)
  data <- Seurat::RunPCA(data)
  data <- Seurat::RunUMAP(data, dims = 1:dims_pca, n.epochs = 500)
  data <- Seurat::FindNeighbors(data, dims = 1:dims_pca)
  data <- Seurat::FindClusters(data, resolution = resolution)
  data
}

#' Subset to one or more cell types (from any column) and re-embed them
#'
#' Accepts a vector of cell types and the metadata column to pull them from,
#' subsets, then re-runs HVG -> scale -> PCA -> UMAP -> neighbours -> clusters
#' on the subset. Assumes the data are already normalized.
#'
#' @param data A Seurat object.
#' @param cell_types Character vector of cell type(s) to keep.
#' @param cell_type_col Metadata column to pull them from (default "broad_cell_types").
#' @param assay_name Assay to use (default "RNA").
#' @param num_hvg Number of variable features (default 2000).
#' @param dims_pca PCA dims used for UMAP and neighbours (default 50).
#' @param resolution Clustering resolution (default 1.0).
#' @param min_dist UMAP min.dist (default 0.3).
#' @param kneigh Neighbours for UMAP (n.neighbors) and the graph (k.param) (default 15).
#' @return The subset, re-embedded Seurat object.
#' @export
process_and_extract_cell_types <- function(data, cell_types,
                                           cell_type_col = "broad_cell_types",
                                           assay_name = "RNA", num_hvg = 2000,
                                           dims_pca = 50, resolution = 1.0,
                                           min_dist = 0.3, kneigh = 15) {
  Seurat::DefaultAssay(data) <- assay_name
  if (!cell_type_col %in% colnames(data@meta.data))
    stop("process_and_extract_cell_types(): column '", cell_type_col, "' not found.")
  keep <- colnames(data)[data@meta.data[[cell_type_col]] %in% cell_types]
  if (length(keep) == 0)
    stop("process_and_extract_cell_types(): no cells match {",
         paste(cell_types, collapse = ", "), "} in column '", cell_type_col, "'.")
  data_sub <- subset(data, cells = keep)
  data_sub <- Seurat::FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = num_hvg)
  data_sub <- Seurat::ScaleData(data_sub)
  data_sub <- Seurat::RunPCA(data_sub)
  data_sub <- Seurat::RunUMAP(data_sub, dims = 1:dims_pca, n.epochs = 500,
                              min.dist = min_dist, n.neighbors = kneigh)
  data_sub <- Seurat::FindNeighbors(data_sub, dims = 1:dims_pca, k.param = kneigh)
  data_sub <- Seurat::FindClusters(data_sub, resolution = resolution)
  data_sub
}

#' Cluster only — re-run FindClusters on an existing neighbour graph
#'
#' Computes neither UMAP nor a neighbour graph; it re-clusters using a graph a
#' prior step already built (e.g. from [process_rna()] or [integrate_data()]).
#' Handy for trying resolutions cheaply.
#'
#' @param data A Seurat object that already carries a neighbour graph.
#' @param resolution Clustering resolution (default 1.0).
#' @param graph_name Graph to cluster on (default: Seurat's active graph).
#' @param cluster_name Name for the resulting cluster column (optional).
#' @return The Seurat object with the new clustering.
#' @export
clustering <- function(data, resolution = 1.0, graph_name = NULL, cluster_name = NULL) {
  if (length(data@graphs) == 0)
    stop("clustering(): no neighbour graph found. Run process_rna() / ",
         "FindNeighbors() first, or use process_and_extract_cell_types().")
  args <- list(object = data, resolution = resolution)
  if (!is.null(graph_name))   args$graph.name   <- graph_name
  if (!is.null(cluster_name)) args$cluster.name <- cluster_name
  do.call(Seurat::FindClusters, args)
}
