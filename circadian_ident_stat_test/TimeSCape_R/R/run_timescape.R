# ============================================================================
# TimeSCape v0.2 — R Main Analysis Pipeline
# Seurat-native version
# ============================================================================
#
# DIFFERENCES FROM MATLAB v0.2:
# 1. INPUT: Accepts Seurat objects directly. Cell type and ZT time are
#    specified as metadata column names (e.g. celltype_col="cell_type",
#    zt_col="ZT_str"). No manual tmeta table needed — ZT times are parsed
#    automatically from strings like "ZT00", "ZT03", "ZT12".
#    If your ZT column uses a different format, use build_tmeta()
#    to inspect and adjust the parsed numeric values before running.
# 2. NLS SOLVER: minpack.lm::nlsLM (Levenberg-Marquardt) instead of MATLAB
#    Trust-Region. Results are numerically nearly identical (<0.1% difference).
# 3. NORMALIZATION: manual lib-size norm (equivalent to scGEAtoolbox's
#    pkg.norm_libsize(X, 1e4) + log1p). Seurat's NormalizedData slot is also
#    supported if already computed.
# 4. PARALLELISM: future.apply::future_lapply() per gene block.
#    Set plan(multisession, workers=N) before calling run_timescape().
# 5. GUI: Shiny web app (app.R) instead of MATLAB uicontrol.
# 6. HEATMAP: pheatmap instead of MATLAB imagesc. Data/scaling identical.
# 7. PLOTS: ggplot2 instead of MATLAB plot/violinplot.
# ============================================================================


# ── Top-level block-fitting helper ────────────────────────────────────────────
# Defined OUTSIDE run_timescape() so its enclosing environment is the
# source-file global scope (where estimate_phaseR lives).  This is critical
# for parallel workers: future serialises .timescape_fit_block and its
# closure chain — if the function were defined inside run_timescape(), the
# entire function frame (including the full Seurat object and global expression
# matrix) would be captured in the closure and sent to every worker.
#
# By defining it here at file-scope, the closure only contains this file's
# globals (primarily estimate_phaseR), which is the only free variable it uses.
#
#' @keywords internal
.timescape_fit_block <- function(g_idx_block, X_sp, zt_v, zts, nz, p12,
                                  phase_fn = estimate_phaseR,
                                  norm_fac  = NULL) {
  nb        <- length(g_idx_block)

  # Densify only this block (≤200 genes × n_cells_ct doubles).
  # If norm_fac is provided (lib_size path), apply log1p-normalisation here
  # rather than pre-normalising the full sparse slice — avoids one extra
  # genes × cells copy in the caller before serialisation to workers.
  X_blk_raw <- as.matrix(X_sp[g_idx_block, , drop = FALSE])
  X_blk <- if (!is.null(norm_fac)) log1p(X_blk_raw * norm_fac[col(X_blk_raw)])
            else X_blk_raw
  rm(X_blk_raw)

  res <- list(
    amp   = numeric(nb),
    mesor = numeric(nb),
    acro  = numeric(nb),
    pval  = numeric(nb),
    rho   = numeric(nb),
    pmac  = numeric(nb),
    R0mat = matrix(NA_real_, nrow = nb, ncol = nz)
  )

  for (i in seq_len(nb)) {
    # Per-ZT expression vectors for gene i
    Xg_zts <- lapply(zts, function(z) {
      idx <- which(zt_v == z)
      if (length(idx) == 0L) return(numeric(0))
      as.numeric(X_blk[i, idx])
    })

    fit          <- phase_fn(Xg_zts, zts, p12, "Ftest")
    res$amp[i]   <- fit$amp
    res$mesor[i] <- fit$mesor
    res$acro[i]  <- fit$acrophase
    res$pval[i]  <- fit$p_value
    res$rho[i]   <- fit$rho
    res$pmac[i]  <- fit$p_value_macro

    # Per-ZT means (used later to build T2 / heatmap matrix)
    for (j in seq_along(zts)) {
      idx <- which(zt_v == zts[j])
      if (length(idx) > 0L)
        res$R0mat[i, j] <- mean(as.numeric(X_blk[i, idx]), na.rm = TRUE)
    }
  }

  res
}


#' Extract metadata from a Seurat or SingleCellExperiment object
#'
#' @param obj  A Seurat or SingleCellExperiment object.
#' @return A plain data.frame of cell-level metadata.
#' @keywords internal
.get_meta <- function(obj) {
  if (inherits(obj, "Seurat")) {
    return(obj@meta.data)
  } else if (inherits(obj, c("SingleCellExperiment", "SummarizedExperiment"))) {
    return(as.data.frame(SummarizedExperiment::colData(obj)))
  }
  stop("obj must be a Seurat or SingleCellExperiment object. Got: ", class(obj)[1])
}


#' Get expression matrix from Seurat or SingleCellExperiment (handles v4/v5 API)
#'
#' @param obj            A Seurat or SingleCellExperiment object.
#' @param use_normalized If TRUE, use pre-normalised slot (Seurat: "data";
#'                       SCE: "logcounts" then "normcounts").
#' @return A genes × cells sparse matrix.
#' @keywords internal
.get_matrix <- function(obj, use_normalized = FALSE) {
  if (inherits(obj, "Seurat")) {
    slot_name <- if (use_normalized) "data" else "counts"

    # ── Seurat v4: GetAssayData(slot = ...) ─────────────────────────────────
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = "RNA", slot = slot_name),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(mat) && prod(dim(mat)) > 0) return(mat)

    # ── Seurat v5: GetAssayData(layer = ...) ────────────────────────────────
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = "RNA", layer = slot_name),
      error = function(e) NULL
    )
    if (!is.null(mat) && prod(dim(mat)) > 0) return(mat)

    # ── Seurat v5: SeuratObject::LayerData (LayerData is in SeuratObject,
    #    not re-exported by Seurat — must call via SeuratObject::) ───────────
    mat <- tryCatch(
      SeuratObject::LayerData(obj, layer = slot_name),
      error = function(e) NULL
    )
    if (!is.null(mat) && prod(dim(mat)) > 0) return(mat)

    # ── Seurat v5: direct bracket accessor obj[["RNA"]][[layer]] ────────────
    mat <- tryCatch(
      obj[["RNA"]][[slot_name]],
      error = function(e) NULL
    )
    if (!is.null(mat) && prod(dim(mat)) > 0) return(mat)

    stop(
      "Cannot extract '", slot_name, "' matrix from Seurat object. ",
      "Make sure the RNA assay is present. ",
      "For Seurat v5 with split layers (multi-sample before integration), ",
      "run: obj <- JoinLayers(obj)  — then retry. ",
      "If your data only has a 'data' (log-normalised) slot and counts is empty, ",
      "use norm_str='logcounts' instead of norm_str='lib_size'."
    )

  } else {
    # ── SingleCellExperiment / SummarizedExperiment ──────────────────────────
    if (use_normalized) {
      tryCatch(
        SummarizedExperiment::assay(obj, "logcounts"),
        error = function(e) tryCatch(
          SummarizedExperiment::assay(obj, "normcounts"),
          error = function(e2) stop(
            "No 'logcounts' or 'normcounts' assay found. ",
            "Run scuttle::logNormCounts() first, or use norm_str='lib_size'.")))
    } else {
      tryCatch(
        SummarizedExperiment::assay(obj, "counts"),
        error = function(e) SummarizedExperiment::assay(obj, 1))   # first assay as fallback
    }
  }
}

# Keep old name as alias for backward compatibility
.get_seurat_matrix <- function(seurat_obj, use_normalized = FALSE)
  .get_matrix(seurat_obj, use_normalized)


#' Build tmeta from a Seurat or SingleCellExperiment object
#'
#' Extracts unique ZT time-point strings from a metadata column and parses
#' numeric ZT hours.  Supports both Seurat and SingleCellExperiment inputs.
#' Formats recognised: "ZT00", "ZT03", "zt12", "ZT_00", "ZT 0", plain
#' numeric strings "0", "3", etc.
#'
#' @param obj     A Seurat or SingleCellExperiment object.
#' @param zt_col  Name of the metadata column holding ZT strings.
#' @return A data.frame with columns: zt_str (original values), ZT_times (numeric hrs).
#'         Rows are sorted by ZT_times.  Inspect and adjust ZT_times manually
#'         if automatic parsing fails for your naming convention.
#' @export
build_tmeta <- function(obj, zt_col) {
  meta <- .get_meta(obj)
  if (!zt_col %in% colnames(meta))
    stop(sprintf("Column '%s' not found in metadata. Available: %s",
                 zt_col, paste(colnames(meta), collapse=", ")))

  zt_vals  <- sort(unique(as.character(meta[[zt_col]])))
  zt_clean <- gsub("(?i)^ZT[_ ]?", "", zt_vals, perl = TRUE)
  zt_num   <- suppressWarnings(as.numeric(zt_clean))

  if (any(is.na(zt_num)))
    warning(sprintf(
      "Could not parse numeric ZT from: %s\nSet ZT_times manually in the returned data.frame.",
      paste(zt_vals[is.na(zt_num)], collapse=", ")))

  tmeta <- data.frame(zt_str = zt_vals, ZT_times = zt_num, stringsAsFactors = FALSE)
  tmeta <- tmeta[order(tmeta$ZT_times), ]
  rownames(tmeta) <- NULL
  tmeta
}

#' @rdname build_tmeta
#' @export
build_tmeta_from_seurat <- function(seurat_obj, zt_col) build_tmeta(seurat_obj, zt_col)



#' Run Full TimeSCape Circadian Analysis Pipeline
#'
#' Loops over cell types, fits a cosinor model per gene, applies dual
#' significance criterion (F-test AND Pearson correlation), BH-corrects
#' p-values, writes 6 CSVs per cell type, and optionally generates heatmaps.
#'
#' Accepts either a \strong{Seurat} object or a \strong{SingleCellExperiment}
#' (Bioconductor) object — the type is detected automatically.
#'
#' @param obj            A Seurat \emph{or} SingleCellExperiment object.
#'                       Cell-level metadata is read from \code{obj@meta.data}
#'                       (Seurat) or \code{colData(obj)} (SCE).
#' @param celltype_col   Metadata column name for cell type annotations.
#'                       Can be any layer of annotation (e.g. "cell_type",
#'                       "subtype", "seurat_clusters").
#' @param zt_col         Metadata column name for ZT time strings.
#'                       Values like "ZT00","ZT03","ZT12" are parsed automatically.
#' @param tmeta          Optional. data.frame from build_tmeta().
#'                       If NULL it is built automatically from zt_col.
#'                       Provide explicitly if you need to override parsed ZT_times.
#' @param rm_low_conf    Logical. Write confident-only CSV subset (default TRUE).
#' @param period12       Logical. 12-hr period (TRUE) or 24-hr (FALSE, default).
#' @param custom_genelist Character vector of gene names to restrict analysis.
#'                       NULL (default) = all genes.
#' @param custom_celltype Character vector of cell types to process.
#'                       NULL (default) = all cell types found in celltype_col.
#' @param plot_heat      Logical. Generate heatmap PNG after each cell type (default TRUE).
#' @param norm_str       Normalisation strategy:
#'                       \describe{
#'                         \item{"lib_size"}{Library-size to 10 000 counts + log1p.
#'                           Requires a non-empty \code{counts} slot.
#'                           \strong{Recommended} when raw UMI counts are available.}
#'                         \item{"logcounts"}{Use the pre-computed normalised slot:
#'                           Seurat \code{data} slot (NormalizeData / SCTransform output)
#'                           or SCE \code{logcounts} assay.
#'                           Use this when \code{counts} is empty but \code{data} is not.}
#'                         \item{"none"}{Raw counts as-is. Use when data is already
#'                           normalised externally and stored in the counts slot.}
#'                       }
#' @param outdir         Root output directory. Per-cell-type subdirs are created inside.
#' @param group_col      Optional. Metadata column name for a second grouping variable
#'                       (e.g. cancer stage, treatment, replicate).  When provided,
#'                       the analysis is run independently for every
#'                       (cell_type × group) combination, each written to its own
#'                       subdirectory  outdir/CellType_Group/.
#'                       NULL (default) = no secondary grouping (single loop over
#'                       cell types only — equivalent to previous behaviour).
#' @param custom_group   Character vector of group values to restrict analysis.
#'                       NULL (default) = all unique values in group_col.
#' @param custom_zt      Character vector of ZT string values to keep,
#'                       e.g. \code{c("ZT00","ZT06","ZT12","ZT18")}.
#'                       NULL (default) = use all ZT time points found in \code{zt_col}.
#'                       Subsetting ZT reduces the matrix column count and therefore
#'                       peak RAM inside each cell-type loop.
#' @param n_workers      Integer. Number of parallel workers for gene-block fitting.
#'                       \describe{
#'                         \item{1L}{Sequential (plain \code{lapply}) — no data copies,
#'                           minimum RAM. Recommended for near-dense matrices (decontX,
#'                           MAGIC-imputed) or when total RAM is close to the matrix size.}
#'                         \item{>1}{Parallel via \code{future.apply::future_lapply}.
#'                           Only the per-cell-type sparse slice is sent to workers.
#'                           Useful for sparse scRNA-seq matrices with many cell types.
#'                           Set \code{plan(multisession, workers=N)} before calling, or
#'                           pass \code{n_workers > 1} here and the plan is managed here.}
#'                       }
#'
#' @return Named list of lists; one entry per cell type processed, each with:
#'   \describe{
#'     \item{$T1}{Stats data.frame (Genes, Amp, Abs_Amp, Mesor, Acrophase, Acrophase_24,
#'            Period, pvalue, pvalue_adj, Sine_corr, pvalue_corr, pvalue_adj_corr).}
#'     \item{$T2}{Per-ZT mean expression data.frame (Genes + one column per ZT label).}
#'   }
#'
#' @examples
#' \dontrun{
#'   library(Seurat); library(future)
#'
#'   seu <- readRDS("my_seurat.rds")
#'   tmeta <- build_tmeta(seu, zt_col = "ZT_str")
#'   print(tmeta)   # verify ZT_times parsed correctly; edit if needed
#'
#'   results <- run_timescape(seu,
#'                            celltype_col = "cell_type",
#'                            zt_col       = "ZT_str",
#'                            tmeta        = tmeta,
#'                            outdir       = "TimeSCape_output")
#' }
#' @export
run_timescape <- function(
  obj,
  celltype_col    = "cell_type",
  zt_col          = "ZT_str",
  tmeta           = NULL,
  rm_low_conf     = TRUE,
  period12        = FALSE,
  custom_genelist = NULL,
  custom_celltype = NULL,
  plot_heat       = TRUE,
  norm_str        = "lib_size",
  outdir          = getwd(),
  group_col       = NULL,
  custom_group    = NULL,
  custom_zt       = NULL,
  n_workers       = 1L
) {
  tic <- proc.time()["elapsed"]

  # ── Validate inputs ─────────────────────────────────────────────────────────
  if (!inherits(obj, c("Seurat", "SingleCellExperiment", "SummarizedExperiment")))
    stop("obj must be a Seurat or SingleCellExperiment object. Got: ", class(obj)[1])

  n_workers <- max(1L, as.integer(round(n_workers)))

  obj_type <- if (inherits(obj, "Seurat")) "Seurat" else "SCE"
  message("Input type detected: ", obj_type)

  # ── Parallel plan ─────────────────────────────────────────────────────────────
  # n_workers = 1 → plain lapply, no data serialisation, no worker copies.
  # n_workers > 1 → future_lapply; sends only the per-cell-type sparse slice.
  #
  # NOTE: .timescape_fit_block is defined at file scope (outside this function)
  # so its closure captures only estimate_phaseR and base R functions — it does
  # NOT capture the Seurat object or global expression matrix even in parallel.
  if (n_workers == 1L) {
    message("  Workers: 1 (sequential — no copies, recommended for dense matrices)")
    use_future <- FALSE
  } else {
    use_future <- TRUE
    cur_plan   <- future::plan()
    is_seq     <- inherits(cur_plan, "sequential") ||
                  identical(class(cur_plan)[1L], "sequential")
    if (is_seq) {
      message(sprintf("  Workers: %d (starting multisession plan)", n_workers))
      future::plan(future::multisession, workers = n_workers)
      on.exit(future::plan(future::sequential), add = TRUE)
    } else {
      cur_workers <- tryCatch(cur_plan$workers, error = function(e) "?")
      message(sprintf("  Workers: reusing existing plan (%s workers)", cur_workers))
    }
    old_maxsize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB
    on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)
    message("  future.globals.maxSize set to 4 GiB for this call")
  }

  # Normalise "logcounts" alias
  if (norm_str == "seurat") norm_str <- "logcounts"   # unify name

  meta <- .get_meta(obj)
  if (!celltype_col %in% colnames(meta))
    stop(sprintf("celltype_col='%s' not in metadata. Available: %s",
                 celltype_col, paste(colnames(meta), collapse=", ")))
  if (!zt_col %in% colnames(meta))
    stop(sprintf("zt_col='%s' not in metadata. Available: %s",
                 zt_col, paste(colnames(meta), collapse=", ")))

  # ── Optional second grouping column ──────────────────────────────────────────
  use_groups <- !is.null(group_col) &&
                length(trimws(group_col)) > 0 &&
                trimws(group_col) != "None"
  if (use_groups) {
    if (!group_col %in% colnames(meta))
      stop(sprintf("group_col='%s' not in metadata. Available: %s",
                   group_col, paste(colnames(meta), collapse=", ")))
    group_vec  <- as.character(meta[[group_col]])
    groups_all <- sort(unique(group_vec))
    if (!is.null(custom_group))
      groups_all <- groups_all[groups_all %in% custom_group]
    if (length(groups_all) == 0)
      stop("No group values remain after applying custom_group filter.")
    message("  Group column: '", group_col, "'  (",
            length(groups_all), " values: ",
            paste(groups_all, collapse=", "), ")")
  } else {
    group_vec  <- NULL
    groups_all <- NULL
  }

  # ── Build tmeta if not provided ─────────────────────────────────────────────
  if (is.null(tmeta)) {
    tmeta <- build_tmeta(obj, zt_col)
    message("Auto-built tmeta from '", zt_col, "':")
    print(tmeta)
  }
  if (any(is.na(tmeta$ZT_times)))
    stop("tmeta$ZT_times contains NA. Parse ZT times manually and pass via tmeta=.")

  per_label <- if (period12) "_period_12_" else "_period_24_"
  message("TimeSCape — period: ", if (period12) "12 hr" else "24 hr")
  message("  Normalisation: ", norm_str)

  if (norm_str == "none")
    message("  ⚠  norm_str='none': ensure counts slot holds raw or log-normalised ",
            "counts, NOT ScaleData output (mean=0, sd=1 per gene — removes MESOR).")

  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # ── Pull global expression matrix (SPARSE — never densify at full scale) ─────
  # "lib_size" / "none" → raw counts slot (sparse stays sparse)
  # "logcounts"         → pre-normalised data slot
  #
  # For decontX outputs (near-dense, 50-80% non-zero), keeping as sparse still
  # avoids the extra copies that as.matrix() would create globally.
  # Densification happens only per gene-block inside .timescape_fit_block().
  use_prebuilt_norm <- (norm_str == "logcounts")
  if (use_prebuilt_norm) {
    X_global          <- .get_matrix(obj, use_normalized = TRUE)   # sparse
    gene_names_global <- rownames(X_global)
    message(sprintf("  Global matrix (log-norm): %d genes × %d cells",
                    nrow(X_global), ncol(X_global)))
  } else {
    X_global_raw      <- .get_matrix(obj, use_normalized = FALSE)  # sparse
    gene_names_global <- rownames(X_global_raw)
    message(sprintf("  Global matrix (raw counts): %d genes × %d cells",
                    nrow(X_global_raw), ncol(X_global_raw)))
  }

  # Optional gene filter
  if (!is.null(custom_genelist)) {
    keep_genes        <- gene_names_global %in% custom_genelist
    gene_names_global <- gene_names_global[keep_genes]
    if (use_prebuilt_norm) {
      X_global     <- X_global[keep_genes, , drop = FALSE]
    } else {
      X_global_raw <- X_global_raw[keep_genes, , drop = FALSE]
    }
    message("  Gene filter applied: ", sum(keep_genes), " / ",
            length(keep_genes), " genes kept")
  }

  # ── Cell-type loop ───────────────────────────────────────────────────────────
  cell_type_vec  <- as.character(meta[[celltype_col]])
  zt_str_vec     <- as.character(meta[[zt_col]])
  cell_types_all <- sort(unique(cell_type_vec))

  if (!is.null(custom_celltype))
    cell_types_all <- cell_types_all[cell_types_all %in% custom_celltype]

  nexpected_zts <- nrow(tmeta)
  all_results   <- list()
  summary_rows  <- list()

  for (ct in cell_types_all) {

    iter_groups <- if (use_groups) groups_all else list(NULL)

    for (grp in iter_groups) {

      # ── Combo label for directory / file naming ──────────────────────────
      ct_safe    <- gsub("[^[:alnum:]_]", "_", trimws(ct))
      grp_safe   <- if (!is.null(grp))
                      gsub("[^[:alnum:]_]", "_", trimws(as.character(grp)))
                    else NULL
      combo_name <- if (!is.null(grp_safe)) paste0(ct_safe, "_", grp_safe)
                    else ct_safe

      if (use_groups) {
        message("\nProcessing: ", ct, "  [", group_col, " = ", grp, "]")
      } else {
        message("\nProcessing: ", ct)
      }

      ct_outdir <- file.path(outdir, combo_name)
      if (!dir.exists(ct_outdir)) dir.create(ct_outdir, showWarnings=FALSE, recursive=TRUE)
      message("  Output → ", ct_outdir)

      # ── Cell mask ─────────────────────────────────────────────────────────
      ct_mask <- cell_type_vec == ct
      if (!is.null(grp))       ct_mask <- ct_mask & (group_vec == grp)
      if (!is.null(custom_zt)) {
        ct_mask <- ct_mask & (zt_str_vec %in% custom_zt)
        message("  ZT filter: keeping ", sum(ct_mask),
                " cells at ZT: ", paste(custom_zt, collapse=", "))
      }

      zt_ct_str <- zt_str_vec[ct_mask]

      # ── Per-cell-type sparse slice ─────────────────────────────────────────
      # NEVER call as.matrix() here — that would densify genes × cells_ct.
      # Densification is deferred to per-block level inside .timescape_fit_block.
      if (use_prebuilt_norm) {
        X_ct <- X_global[, ct_mask, drop = FALSE]

      } else if (norm_str == "lib_size") {
        X_ct_raw <- X_global_raw[, ct_mask, drop = FALSE]
        cs       <- Matrix::colSums(X_ct_raw)
        cs[cs == 0] <- 1  # avoid divide-by-zero for empty cells
        # Store per-cell norm factors as a plain vector.  The raw sparse slice is
        # passed to workers AS-IS; each worker applies log1p( X_blk * nf ) on its
        # own small dense block.  This avoids creating a full genes × cells
        # normalised copy before serialisation — saves one large intermediate matrix.
        norm_fac <- pmax(1e4 / cs, 0)  # pmax guards against negative cs edge cases
        X_ct     <- X_ct_raw
        rm(X_ct_raw); gc(verbose = FALSE)

      } else {  # "none"
        X_ct <- X_global_raw[, ct_mask, drop = FALSE]
      }

      # norm_fac is NULL unless lib_size path was chosen
      if (!exists("norm_fac")) norm_fac <- NULL

      message(sprintf("  Slice: %d genes × %d cells  (%.1f MB sparse)",
                      nrow(X_ct), ncol(X_ct),
                      object.size(X_ct) / 1024^2))

      gene_names <- gene_names_global

      # ── ZT mapping ────────────────────────────────────────────────────────
      zt_num_ct <- tmeta$ZT_times[match(zt_ct_str, tmeta$zt_str)]
      zt_present <- sort(unique(zt_num_ct[!is.na(zt_num_ct)]))
      nzts       <- length(zt_present)
      message("  Time points: ", nzts, " / ", nexpected_zts, " expected")

      if (nzts < 4) {
        warning("  '", ct, "' has fewer than 4 time points — skipping.")
        rm(X_ct); gc(verbose = FALSE)
        next
      }
      if (nzts < nexpected_zts)
        message("  ⚠ Missing ", nexpected_zts - nzts,
                " time point(s) — fitting on available points only.")

      zt_labels <- tmeta$zt_str[match(zt_present, tmeta$ZT_times)]
      ngenes    <- nrow(X_ct)
      if (ngenes == 0) {
        warning("  No genes for '", ct, "' — skipping.")
        rm(X_ct); gc(verbose = FALSE)
        next
      }

      # ── Block cosinor fitting ──────────────────────────────────────────────
      # block_size controls peak dense RAM per iteration:
      #   e.g. 100 genes × 50k cells × 8 bytes = 40 MB per block
      # For decontX / MAGIC (near-dense), smaller blocks reduce RAM spikes.
      n_cells_ct <- ncol(X_ct)
      block_size <- if (n_cells_ct > 50000L) 50L else if (n_cells_ct > 20000L) 100L else 200L
      gene_blocks <- split(seq_len(ngenes), ceiling(seq_len(ngenes) / block_size))
      message(sprintf("  Fitting %d genes in %d blocks of ≤%d  [%d cells, %s]",
                      ngenes, length(gene_blocks), block_size, n_cells_ct,
                      if (use_future) sprintf("%d workers", n_workers) else "sequential"))

      if (!use_future) {
        # ── Sequential: no copies, no serialisation ──────────────────────────
        results_blocks <- lapply(
          gene_blocks, .timescape_fit_block,
          X_sp     = X_ct,
          zt_v     = zt_num_ct,
          zts      = zt_present,
          nz       = nzts,
          p12      = period12,
          norm_fac = norm_fac   # NULL for logcounts/none; vector for lib_size
        )
      } else {
        # ── Parallel: send sparse slice + norm_fac vector to workers ─────────
        # .timescape_fit_block is top-level so its closure captures only
        # estimate_phaseR — not the Seurat object or global expression matrix.
        # norm_fac (lib_size path) replaces sending a pre-normalised dense copy.
        results_blocks <- future.apply::future_lapply(
          gene_blocks,
          FUN      = .timescape_fit_block,
          X_sp     = X_ct,
          zt_v     = zt_num_ct,
          zts      = zt_present,
          nz       = nzts,
          p12      = period12,
          norm_fac = norm_fac,
          future.globals = list(
            .timescape_fit_block = .timescape_fit_block,
            estimate_phaseR      = estimate_phaseR
          ),
          future.packages = c("Matrix", "minpack.lm"),
          future.seed = TRUE
        )
      }

      # ── Free the per-cell-type matrix — no longer needed ──────────────────
      rm(X_ct); gc(verbose = FALSE)

      # ── Flatten blocks ────────────────────────────────────────────────────
      amp_v  <- unlist(lapply(results_blocks, `[[`, "amp"))
      mesor_v<- unlist(lapply(results_blocks, `[[`, "mesor"))
      acro_v <- unlist(lapply(results_blocks, `[[`, "acro"))
      pval_v <- unlist(lapply(results_blocks, `[[`, "pval"))
      rho_v  <- unlist(lapply(results_blocks, `[[`, "rho"))
      pmac_v <- unlist(lapply(results_blocks, `[[`, "pmac"))
      R0_mat <- do.call(rbind, lapply(results_blocks, `[[`, "R0mat"))
      rm(results_blocks); gc(verbose = FALSE)

      # ── Diagnostic: warn when all genes failed ─────────────────────────────
      n_na <- sum(is.na(pval_v))
      if (n_na == ngenes) {
        warning(
          "  ⚠  ALL ", ngenes, " genes returned NA for '", ct, "'.\n",
          "  Possible causes:\n",
          "  1. norm_str='lib_size' but counts slot is empty → use norm_str='logcounts'\n",
          "  2. estimate_phaseR unavailable in parallel workers → retry n_workers=1\n",
          "  3. Severe NaN/Inf in the expression matrix after normalisation\n",
          "  Retry with n_workers=1 to isolate the issue."
        )
      } else if (n_na > 0) {
        message(sprintf("  (%d / %d genes returned NA — likely zero/very low expression)",
                        n_na, ngenes))
      }

      # Wrap acrophase into [0, period)
      pd     <- if (period12) 12 else 24
      acro24 <- acro_v %% 24

      # BH correction — include ALL genes in the denominator, treating genes
      # that failed NLS fitting (NA p-value) as p = 1.  This matches MATLAB
      # behaviour where zero-expression genes contribute to the rank denominator
      # rather than being silently dropped before correction.  The NA rows are
      # filtered from the output table afterwards (line: valid <- !is.na(T1$pvalue)).
      pval_for_adj <- ifelse(is.na(pval_v), 1, pval_v)
      pmac_for_adj <- ifelse(is.na(pmac_v), 1, pmac_v)
      padj_v   <- stats::p.adjust(pval_for_adj, method = "BH")
      padj_mac <- stats::p.adjust(pmac_for_adj, method = "BH")

      # ── Build T1 (stats table) ────────────────────────────────────────────
      T1 <- data.frame(
        Genes           = gene_names,
        Amp             = amp_v,
        Abs_Amp         = abs(amp_v),
        Mesor           = mesor_v,
        Acrophase       = acro_v,
        Acrophase_24    = acro24,
        Period          = rep(pd, ngenes),
        pvalue          = pval_v,
        pvalue_adj      = padj_v,
        Sine_corr       = rho_v,
        pvalue_corr     = pmac_v,
        pvalue_adj_corr = padj_mac,
        stringsAsFactors = FALSE
      )

      # Remove rows where fitting failed
      valid  <- !is.na(T1$pvalue) & !is.na(T1$pvalue_corr)
      T1     <- T1[valid, ]
      R0_mat <- R0_mat[valid, , drop = FALSE]

      # Sort: pvalue_adj_corr ↑, pvalue_adj ↑, Acrophase_24 ↑, Abs_Amp ↓
      ord    <- order(T1$pvalue_adj_corr, T1$pvalue_adj, T1$Acrophase_24, -T1$Abs_Amp)
      T1     <- T1[ord, ]
      R0_mat <- R0_mat[ord, , drop = FALSE]
      rownames(T1) <- NULL

      # ── Build T2 (per-ZT mean expression) ────────────────────────────────
      T2 <- data.frame(Genes = T1$Genes, stringsAsFactors = FALSE)
      for (j in seq_along(zt_labels)) T2[[zt_labels[j]]] <- R0_mat[, j]

      # T3: ZT-normalised (divide by ZT00 column; fallback to next non-zero)
      T3      <- T2
      ref_col <- T3[, 2]
      zero_ref <- ref_col == 0 | is.na(ref_col)
      if (any(zero_ref) && ncol(T3) > 2) ref_col[zero_ref] <- T3[zero_ref, 3]
      ref_col[ref_col == 0] <- 1
      T3[, 2:ncol(T3)] <- T3[, 2:ncol(T3)] / ref_col

      # Confidence flags
      conf_both <- (T1$pvalue < 0.05) & (T1$pvalue_corr < 0.05)
      T1_conf   <- T1[conf_both, ]
      T2_conf   <- T2[conf_both, ]
      T3_conf   <- T3[conf_both, ]

      message("  Genes tested   : ", nrow(T1))
      message("  Confident (F+corr p<0.05): ", sum(conf_both))

      # ── Write CSVs ────────────────────────────────────────────────────────
      fbase <- paste0(combo_name, per_label)
      write.csv(T1, file.path(ct_outdir, paste0(fbase, "circadian_analysis_all.csv")),
                row.names = FALSE)
      write.csv(T2, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean.csv")),
                row.names = FALSE)
      write.csv(T3, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_normalized.csv")),
                row.names = FALSE)
      if (rm_low_conf) {
        write.csv(T1_conf, file.path(ct_outdir, paste0(fbase, "circadian_analysis_confident.csv")),
                  row.names = FALSE)
        write.csv(T2_conf, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_confident.csv")),
                  row.names = FALSE)
        write.csv(T3_conf, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_normalized_confident.csv")),
                  row.names = FALSE)
      }

      # ── Optional heatmap ─────────────────────────────────────────────────
      if (plot_heat && nrow(T1_conf) > 1) {
        tryCatch(
          generate_heatmap(combo_name, strict = TRUE, custom_name = "",
                           circ = FALSE, period12 = period12, outdir = ct_outdir),
          error = function(e) message("  Heatmap warning: ", e$message)
        )
      }

      all_results[[combo_name]] <- list(T1 = T1, T2 = T2)
      summary_rows[[combo_name]] <- data.frame(
        CellType          = ct,
        Group             = if (!is.null(grp)) grp else NA_character_,
        ComboName         = combo_name,
        NumCells          = sum(ct_mask),
        NumGenes          = ngenes,
        NumTested         = nrow(T1),
        NumConfident_Both = sum(conf_both),
        NumNonConfident   = nrow(T1) - sum(conf_both),
        stringsAsFactors  = FALSE
      )

      norm_fac <- NULL  # reset for next iteration
      rm(T1, T2, T3, T1_conf, T2_conf, T3_conf, R0_mat,
         amp_v, mesor_v, acro_v, pval_v, rho_v, pmac_v,
         padj_v, padj_mac, valid, conf_both, ord)
      gc(verbose = FALSE)

    }  # end group loop
  }  # end cell-type loop

  # ── Cross-cell-type summary CSV ──────────────────────────────────────────────
  if (length(summary_rows) > 0) {
    summary_tbl <- do.call(rbind, summary_rows)
    write.csv(summary_tbl,
              file.path(outdir, paste0("all_cell_types", per_label, "summary_results.csv")),
              row.names = FALSE)
  }

  elapsed <- proc.time()["elapsed"] - tic
  message(sprintf("\n=== TimeSCape complete. Total time: %.1f s ===", elapsed))
  invisible(all_results)
}
