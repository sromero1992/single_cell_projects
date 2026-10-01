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
#' @param normalization "LogNormalize" (default, re-uses existing normalized data)
#'   or "SCT" (re-run SCTransform on the subset — its HVGs/variance differ from the
#'   global object). SCT writes a separate 'SCT' assay used only for this embedding;
#'   the RNA assay is restored as default before returning, so DE stays on RNA.
#' @param vars_to_regress Covariates for SCTransform's vars.to.regress (e.g.
#'   "percent_mt"). Ignored for LogNormalize. Do NOT pass nCount_RNA (SCT models
#'   depth itself).
#' @param harmonize_by Optional batch column (e.g. "SampleID"). When set, BOTH
#'   embeddings are produced so you can compare batch effects: the un-integrated
#'   PCA track (`umap_none` / `clusters_none`) AND the Harmony-corrected track
#'   (`umap_harmony` / `clusters_harmony`), with Idents set to `clusters_harmony`.
#'   NULL (default) produces only the un-integrated track and sets Idents to
#'   `clusters_none`. Naming matches the subannotation scripts' convention.
#' @return The subset, re-embedded Seurat object.
#' @export
process_and_extract_cell_types <- function(data, cell_types,
                                           cell_type_col = "broad_cell_types",
                                           assay_name = "RNA", num_hvg = 2000,
                                           dims_pca = 50, resolution = 1.0,
                                           min_dist = 0.3, kneigh = 15,
                                           normalization = c("LogNormalize", "SCT"),
                                           vars_to_regress = NULL,
                                           harmonize_by = NULL) {
  normalization <- match.arg(normalization)
  Seurat::DefaultAssay(data) <- assay_name
  if (!cell_type_col %in% colnames(data@meta.data))
    stop("process_and_extract_cell_types(): column '", cell_type_col, "' not found.")
  keep <- colnames(data)[data@meta.data[[cell_type_col]] %in% cell_types]
  if (length(keep) == 0)
    stop("process_and_extract_cell_types(): no cells match {",
         paste(cell_types, collapse = ", "), "} in column '", cell_type_col, "'.")
  data_sub <- subset(data, cells = keep)
  # Drop reductions/graphs inherited from the PARENT object — they were computed on
  # all cells and are meaningless for this subset; we recompute everything below.
  data_sub@reductions <- list()
  data_sub@graphs     <- list()

  if (normalization == "SCT") {
    if (!requireNamespace("sctransform", quietly = TRUE))
      stop("normalization='SCT' needs the 'sctransform' package installed.")
    vtr <- intersect(vars_to_regress, colnames(data_sub@meta.data))
    if (length(vtr) < length(vars_to_regress))
      warning("process_and_extract_cell_types(): vars_to_regress not in metadata, dropped: ",
              paste(setdiff(vars_to_regress, vtr), collapse = ", "))
    if (length(vtr)) {
      message("  [SCT] regressing out: ", paste(vtr, collapse = ", "),
              " (affects the SCT embedding/clustering only; RNA assay & DE unchanged).")
      if ("percent_mt" %in% vtr)
        message("  [SCT] NOTE: percent_mt removes the low-quality/dying-cell axis for ",
                "cleaner clusters, but if mito fraction genuinely differs by condition/",
                "cell type (plausible tumor/polyp vs WT) this can erase real signal and ",
                "shift cluster boundaries. Compare with normalization set to LogNormalize if unsure.")
    }
    data_sub <- Seurat::SCTransform(data_sub, assay = assay_name, new.assay = "SCT",
                                    variable.features.n = num_hvg,
                                    vars.to.regress = if (length(vtr)) vtr else NULL,
                                    vst.flavor = "v2", verbose = FALSE)
    Seurat::DefaultAssay(data_sub) <- "SCT"
    data_sub <- Seurat::RunPCA(data_sub, npcs = dims_pca, verbose = FALSE)
  } else {
    data_sub <- Seurat::FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = num_hvg)
    data_sub <- Seurat::ScaleData(data_sub)
    data_sub <- Seurat::RunPCA(data_sub)
  }

  # neighbours -> clusters -> UMAP for one reduction (mirrors integrate_data()).
  embed <- function(o, reduction, graph, clusters, umap) {
    n_dims <- min(dims_pca, ncol(SeuratObject::Embeddings(o, reduction)))
    o <- Seurat::FindNeighbors(o, reduction = reduction, dims = 1:n_dims,
                               k.param = kneigh, graph.name = graph, verbose = FALSE)
    o <- Seurat::FindClusters(o, resolution = resolution, graph.name = graph,
                              cluster.name = clusters, verbose = FALSE)
    Seurat::RunUMAP(o, reduction = reduction, dims = 1:n_dims, n.epochs = 500,
                    min.dist = min_dist, n.neighbors = kneigh,
                    reduction.name = umap, verbose = FALSE)
  }

  # Track A: un-integrated PCA — always produced (diagnostic for batch effects).
  data_sub <- embed(data_sub, "pca", "pca_nn", "clusters_none", "umap_none")
  active   <- "clusters_none"

  # Track B: Harmony-corrected — only when a batch column is supplied.
  if (!is.null(harmonize_by)) {
    if (!requireNamespace("harmony", quietly = TRUE))
      stop("process_and_extract_cell_types(): harmonize_by needs the 'harmony' package.")
    if (!harmonize_by %in% colnames(data_sub@meta.data))
      stop("process_and_extract_cell_types(): harmonize_by column '", harmonize_by, "' not found.")
    data_sub <- harmony::RunHarmony(data_sub, group.by.vars = harmonize_by,
                                    reduction.use = "pca", dims.use = 1:dims_pca,
                                    reduction.save = "harmony")
    data_sub <- embed(data_sub, "harmony", "harmony_nn", "clusters_harmony", "umap_harmony")
    active   <- "clusters_harmony"
  }

  Seurat::Idents(data_sub) <- active
  Seurat::DefaultAssay(data_sub) <- assay_name   # RNA back as default for downstream DE
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
