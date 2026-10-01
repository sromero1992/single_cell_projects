# =============================================================================
# integrate.R — normalize + HVG + PCA + integration + clustering + UMAP
# =============================================================================
# Two routes, matching Script 01's INTEGRATION_METHOD:
#   "RunHarmony"      — direct harmony call on the joined matrix, no split/join.
#   "IntegrateLayers" — Seurat v5 split-layer HarmonyIntegration path.
# Produces TWO embeddings (like the legacy Track A / Track B):
#   * unintegrated (diagnostic): reduction "pca"     -> clusters_unintegrated / umap_unintegrated
#   * harmony      (working)    : reduction "harmony" -> clusters_harmony / umap_harmony
# The harmony clustering is set as the active identity, so downstream always
# works on harmony unless you switch reductions yourself.
# =============================================================================

#' Normalize, reduce, integrate, cluster and embed a merged object
#'
#' @param obj A merged Seurat object.
#' @param method "RunHarmony" (default) or "IntegrateLayers".
#' @param group_by Batch/sample variable to correct on (default "SampleID").
#' @param n_pcs Number of PCs / harmony dims (default 30).
#' @param n_features HVG count (default 2000).
#' @param resolution Clustering resolution (default 1.0).
#' @param unintegrated Also compute the unintegrated pca embedding
#'   (`clusters_unintegrated` / `umap_unintegrated`) as a diagnostic (default TRUE).
#' @param normalization Embedding front-end: `"LogNormalize"` (default) or `"SCT"`
#'   (SCTransform v2 / glmGamPoi then PCA). Affects ONLY the SCT assay and the
#'   embedding/clustering; the RNA assay is left LogNormalized and reset as the
#'   default, so downstream DE is unchanged. `"SCT"` is supported with
#'   `method = "RunHarmony"` only.
#' @param vars_to_regress Character vector of covariates for SCTransform's
#'   `vars.to.regress` (e.g. `"percent_mt"`), used only when `normalization = "SCT"`.
#'   Do NOT pass `nCount_RNA` — SCT already models sequencing depth. `NULL` (default)
#'   regresses nothing.
#' @return The object with `pca` + `harmony` reductions, `umap_unintegrated` +
#'   `umap_harmony`, `clusters_unintegrated` + `clusters_harmony`, and Idents set
#'   to `clusters_harmony`.
#' @export
integrate_data <- function(obj, method = c("RunHarmony", "IntegrateLayers"),
                           group_by = "SampleID", n_pcs = 30, n_features = 2000,
                           resolution = 1.0, unintegrated = TRUE,
                           normalization = c("LogNormalize", "SCT"),
                           vars_to_regress = NULL) {
  method <- match.arg(method)
  normalization <- match.arg(normalization)
  stopifnot(.is_seurat(obj))
  Seurat::DefaultAssay(obj) <- "RNA"

  # normalize + HVG + scale + PCA front-end. SCTransform is an OPTIONAL swap for
  # the LogNormalize path used ONLY for the embedding/clustering (annotation): it
  # writes a separate 'SCT' assay and PCA, while the RNA assay is left untouched so
  # DE (Scripts 07/09) still reads LogNormalized RNA. SCT already models sequencing
  # depth internally, so do NOT put nCount_RNA in vars_to_regress; percent_mt is a
  # legitimate quality covariate to regress. SCT is wired to the RunHarmony route.
  .sct <- function(o, hvg) {
    if (!requireNamespace("sctransform", quietly = TRUE))
      stop("normalization='SCT' needs the 'sctransform' package installed.")
    vtr <- intersect(vars_to_regress, colnames(o@meta.data))
    if (length(vtr) < length(vars_to_regress))
      warning("integrate_data(): vars_to_regress not in metadata, dropped: ",
              paste(setdiff(vars_to_regress, vtr), collapse = ", "))
    if (length(vtr)) {
      message("  [SCT] regressing out: ", paste(vtr, collapse = ", "),
              " (affects the SCT embedding/clustering only; RNA assay & DE unchanged).")
      if ("percent_mt" %in% vtr)
        message("  [SCT] NOTE: percent_mt removes the low-quality/dying-cell axis for ",
                "cleaner clusters, but if mito fraction genuinely differs by condition/",
                "cell type (plausible tumor/polyp vs WT) this can erase real signal and ",
                "shift cluster boundaries. Compare with NORMALIZATION toggled off if unsure.")
    }
    o <- Seurat::SCTransform(o, assay = "RNA", new.assay = "SCT",
                             variable.features.n = hvg,
                             vars.to.regress = if (length(vtr)) vtr else NULL,
                             vst.flavor = "v2", verbose = FALSE)   # v2 uses glmGamPoi
    Seurat::DefaultAssay(o) <- "SCT"
    Seurat::RunPCA(o, npcs = n_pcs, verbose = FALSE)
  }

  # --- normalize + HVG + scale + PCA (+ harmony) -----------------------------
  if (method == "RunHarmony") {
    if (!requireNamespace("harmony", quietly = TRUE))
      stop("harmony is not installed. install.packages('harmony').")
    if (inherits(obj[["RNA"]], "Assay5") &&
        length(SeuratObject::Layers(obj[["RNA"]], search = "counts")) > 1)
      obj <- SeuratObject::JoinLayers(obj)
    if (normalization == "SCT") {
      obj <- .sct(obj, n_features)
    } else {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_features, verbose = FALSE)
      obj <- Seurat::ScaleData(obj, verbose = FALSE)
      obj <- Seurat::RunPCA(obj, npcs = n_pcs, verbose = FALSE)
    }
    obj <- harmony::RunHarmony(obj, group.by.vars = group_by,
                               reduction.use = "pca", dims.use = 1:n_pcs,
                               reduction.save = "harmony")
  } else {
    if (normalization == "SCT")
      stop("integrate_data(): normalization='SCT' is supported with method='RunHarmony' ",
           "only (SCTIntegration is not wired into the IntegrateLayers path). ",
           "Use method='RunHarmony'.")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[group_by]])
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_features, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = n_pcs, verbose = FALSE)
    obj <- Seurat::IntegrateLayers(obj, method = Seurat::HarmonyIntegration,
                                   orig.reduction = "pca", new.reduction = "harmony",
                                   verbose = FALSE)
    obj <- SeuratObject::JoinLayers(obj)
  }

  # --- helper: neighbors -> clusters -> umap for one reduction ---------------
  embed <- function(obj, reduction, graph, clusters, umap) {
    n_dims <- min(n_pcs, ncol(SeuratObject::Embeddings(obj, reduction)))
    obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = 1:n_dims,
                                 graph.name = graph, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = resolution, graph.name = graph,
                                cluster.name = clusters, verbose = FALSE)
    Seurat::RunUMAP(obj, reduction = reduction, dims = 1:n_dims,
                    reduction.name = umap, verbose = FALSE)
  }

  # Track A: unintegrated PCA (diagnostic)
  if (unintegrated)
    obj <- embed(obj, "pca", "pca_nn", "clusters_unintegrated", "umap_unintegrated")

  # Track B: harmony (working default)
  obj <- embed(obj, "harmony", "harmony_nn", "clusters_harmony", "umap_harmony")

  # Work on harmony by default.
  Seurat::Idents(obj) <- "clusters_harmony"

  # RNA is the source of truth for DE downstream; reset it as the default assay
  # (the SCT assay stays on the object for reference, but is not the default).
  Seurat::DefaultAssay(obj) <- "RNA"

  .stamp(obj, "integrate_data", list(method = method, group_by = group_by,
                                     n_pcs = n_pcs, resolution = resolution,
                                     unintegrated = unintegrated,
                                     normalization = normalization,
                                     vars_to_regress = paste(vars_to_regress, collapse = ",")))
}
