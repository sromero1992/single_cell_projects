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
#    If your ZT column uses a different format, use build_tmeta_from_seurat()
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
    tryCatch(
      Seurat::GetAssayData(obj, assay = "RNA", slot = slot_name),
      error = function(e) tryCatch(
        Seurat::LayerData(obj, layer = slot_name),
        error = function(e2) stop(
          "Cannot extract matrix from Seurat object (RNA assay missing?). ",
          e2$message)))
  } else {
    # SingleCellExperiment / SummarizedExperiment
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


#' Get expression matrix from Seurat (handles v4 and v5 API)
#'
#' @param seurat_obj A Seurat object.
#' @param use_normalized If TRUE, use normalized data slot; else raw counts.
#' @return A genes × cells sparse or dense matrix.
#'
#' @keywords internal
.get_seurat_matrix <- function(seurat_obj, use_normalized = FALSE) {
  slot_name <- if (use_normalized) "data" else "counts"
  # Seurat v5 uses Layers API; v4 uses GetAssayData
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = slot_name),
    error = function(e) {
      tryCatch(
        Seurat::LayerData(seurat_obj, layer = slot_name),
        error = function(e2) {
          stop("Cannot extract expression matrix from Seurat object. ",
               "Make sure the RNA assay is present. Original error: ", e2$message)
        })
    })
}


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
#' @param tmeta          Optional. data.frame from build_tmeta_from_seurat().
#'                       If NULL it is built automatically from zt_col.
#'                       Provide explicitly if you need to override parsed ZT_times.
#' @param rm_low_conf    Logical. Write confident-only CSV subset (default TRUE).
#' @param period12       Logical. 12-hr period (TRUE) or 24-hr (FALSE, default).
#' @param custom_genelist Character vector of gene names to restrict analysis.
#'                       NULL (default) = all genes.
#' @param custom_celltype Character vector of cell types to process.
#'                       NULL (default) = all cell types found in celltype_col.
#' @param plot_heat      Logical. Generate heatmap PNG after each cell type (default TRUE).
#' @param norm_str       "lib_size" (default): library-size norm + log1p, identical
#'                       to MATLAB pkg.norm_libsize(X,1e4).
#'                       "seurat": use the Seurat normalized data slot directly.
#'                       "logcounts": use the SCE \code{logcounts} assay (or
#'                       \code{normcounts} as fallback) — SCE equivalent of
#'                       "seurat".  Seurat objects also accept "logcounts" as
#'                       an alias for "seurat".
#'                       "none": use raw counts as-is (pass pre-normalised data or
#'                       raw counts when no transformation is desired).
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
#'
#' @return Named list of lists; one entry per cell type processed, each with:
#'   $T1 — stats data.frame (Genes, Amp, Abs_Amp, Mesor, Acrophase, Acrophase_24,
#'           Period, pvalue, pvalue_adj, Sine_corr, pvalue_corr, pvalue_adj_corr)
#'   $T2 — per-ZT mean expression data.frame (Genes + one column per ZT label)
#'
#' @examples
#' \dontrun{
#'   library(Seurat); library(future)
#'   plan(multisession, workers = 4)   # parallel fitting
#'
#'   seu <- readRDS("my_seurat.rds")
#'
#'   # Inspect available metadata columns and ZT values
#'   head(seu@meta.data)
#'   # Suppose metadata has: cell_type, ZT_str = c("ZT00","ZT03","ZT06",...)
#'
#'   tmeta <- build_tmeta_from_seurat(seu, zt_col = "ZT_str")
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
  celltype_col   = "cell_type",
  zt_col         = "ZT_str",
  tmeta          = NULL,
  rm_low_conf    = TRUE,
  period12       = FALSE,
  custom_genelist = NULL,
  custom_celltype = NULL,
  plot_heat      = TRUE,
  norm_str       = "lib_size",
  outdir         = getwd(),
  group_col      = NULL,
  custom_group   = NULL
) {
  tic <- proc.time()["elapsed"]

  # ── Validate inputs ─────────────────────────────────────────────────────────
  if (!inherits(obj, c("Seurat", "SingleCellExperiment", "SummarizedExperiment")))
    stop("obj must be a Seurat or SingleCellExperiment object. Got: ", class(obj)[1])

  obj_type <- if (inherits(obj, "Seurat")) "Seurat" else "SCE"
  message("Input type detected: ", obj_type)

  # "logcounts" is accepted as an alias for "seurat" (pre-normalised slot)
  if (norm_str == "logcounts") norm_str <- if (obj_type == "Seurat") "seurat" else "logcounts"

  meta <- .get_meta(obj)
  if (!celltype_col %in% colnames(meta))
    stop(sprintf("celltype_col='%s' not in metadata. Available: %s",
                 celltype_col, paste(colnames(meta), collapse=", ")))
  if (!zt_col %in% colnames(meta))
    stop(sprintf("zt_col='%s' not in metadata. Available: %s",
                 zt_col, paste(colnames(meta), collapse=", ")))

  # ── Optional second grouping column ──────────────────────────────────────────
  # When group_col is provided, the analysis loops over every
  # (cell_type × group) combination independently.  Each combo is saved to its
  # own sub-directory outdir/CellType_Group/ and file prefix CellType_Group_.
  # Because lib_size normalisation is per-cell, splitting by group is always
  # mathematically equivalent to processing all cells together — no bias.
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
    groups_all <- NULL   # sentinel: single NULL means "no inner loop"
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

  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # ── Normalisation strategy ────────────────────────────────────────────────
  # norm_str options:
  #   "lib_size" — library-size to 10k + log1p, computed per-cell-type inside
  #                the loop below to minimise peak RAM.  IDENTICAL result to
  #                normalising the full matrix first because lib-size is a
  #                per-cell operation (each cell divided by its own sum).
  #                Cancer-stage mixing, replicate mixing, etc. do NOT affect
  #                the result.  Safe to use even when stages are in the same
  #                Seurat object.
  #   "seurat"   — use Seurat's pre-computed NormalizedData slot (RNA$data or
  #                GetAssayData slot="data").  Use this if you already ran
  #                NormalizeData() / SCTransform.
  #   "none"     — use raw counts as-is from the counts slot.  Use when you
  #                have already normalised and stored results in counts, or
  #                when you explicitly want no transformation.
  message("  Normalisation: ", norm_str)

  # For "seurat" we pull the full normalised matrix once (it is already in RAM
  # inside the Seurat object, so no extra cost).
  # For "lib_size" and "none" we defer to the cell-type loop.
  use_prebuilt_norm <- norm_str %in% c("seurat", "logcounts")
  if (use_prebuilt_norm) {
    X_global          <- as.matrix(.get_matrix(obj, use_normalized = TRUE))
    gene_names_global <- rownames(X_global)
  } else {
    X_global_raw      <- .get_matrix(obj, use_normalized = FALSE)
    gene_names_global <- rownames(X_global_raw)
  }

  # Optional gene filter applied to the global gene list
  if (!is.null(custom_genelist)) {
    keep_genes <- gene_names_global %in% custom_genelist
    gene_names_global <- gene_names_global[keep_genes]
    if (use_prebuilt_norm) {
      X_global <- X_global[keep_genes, , drop = FALSE]
    } else {
      X_global_raw <- X_global_raw[keep_genes, , drop = FALSE]
    }
  }

  # ── Cell-type loop ───────────────────────────────────────────────────────────
  cell_type_vec   <- as.character(meta[[celltype_col]])
  zt_str_vec      <- as.character(meta[[zt_col]])
  cell_types_all  <- sort(unique(cell_type_vec))

  if (!is.null(custom_celltype))
    cell_types_all <- cell_types_all[cell_types_all %in% custom_celltype]

  nexpected_zts   <- nrow(tmeta)
  all_results     <- list()
  summary_rows    <- list()

  for (ct in cell_types_all) {

    # When group_col is active we loop over group values; otherwise a single
    # NULL pass replicates the original single-loop behaviour.
    iter_groups <- if (use_groups) groups_all else list(NULL)

    for (grp in iter_groups) {

      # ── Combo label for directory / file naming ──────────────────────────
      ct_safe   <- gsub("[^[:alnum:]_]", "_", trimws(ct))
      grp_safe  <- if (!is.null(grp))
                     gsub("[^[:alnum:]_]", "_", trimws(as.character(grp)))
                   else NULL
      combo_name <- if (!is.null(grp_safe)) paste0(ct_safe, "_", grp_safe)
                    else ct_safe

      if (use_groups) {
        message("\nProcessing cell type: ", ct,
                "  [", group_col, " = ", grp, "]")
      } else {
        message("\nProcessing cell type: ", ct)
      }

      # Per-(cell-type × group) output directory
      ct_outdir <- file.path(outdir, combo_name)
      if (!dir.exists(ct_outdir)) dir.create(ct_outdir, showWarnings=FALSE, recursive=TRUE)
      message("  Output directory: ", ct_outdir)

      # Subset cells — cell type, and optionally group
      ct_mask   <- cell_type_vec == ct
      if (!is.null(grp)) ct_mask <- ct_mask & (group_vec == grp)
      zt_ct_str <- zt_str_vec[ct_mask]

    # ── Per-cell-type normalisation (low-RAM path) ─────────────────────────
    # Extracting and normalising only the cells for this cell type keeps peak
    # RAM to one dense cell-type matrix at a time instead of the full dataset.
    if (use_prebuilt_norm) {
      X_ct <- as.matrix(X_global[, ct_mask, drop = FALSE])
    } else if (norm_str == "lib_size") {
      X_ct_raw  <- as.matrix(X_global_raw[, ct_mask, drop = FALSE])
      csumsX    <- colSums(X_ct_raw)
      csumsX[csumsX == 0] <- 1
      X_ct      <- log1p(t(t(X_ct_raw) / csumsX * 1e4))
      rm(X_ct_raw)
    } else {  # "none"
      X_ct <- as.matrix(X_global_raw[, ct_mask, drop = FALSE])
    }

    gene_names <- gene_names_global

    # Map ZT strings → numeric times using tmeta
    zt_num_ct <- tmeta$ZT_times[match(zt_ct_str, tmeta$zt_str)]

    # Discover which ZT time points actually have cells
    zt_present <- sort(unique(zt_num_ct[!is.na(zt_num_ct)]))
    nzts       <- length(zt_present)
    message("  Time points found: ", nzts, " / ", nexpected_zts, " expected")

    if (nzts < 4) {
      warning("  Cell type '", ct, "' has fewer than 4 time points — skipping.")
      next
    }
    if (nzts < nexpected_zts)
      message("  ⚠ Missing ", nexpected_zts - nzts, " time point(s) — fitting on available points only.")

    # ZT label strings in order (for column names in output CSVs)
    zt_labels <- tmeta$zt_str[match(zt_present, tmeta$ZT_times)]

    ngenes    <- nrow(X_ct)
    if (ngenes == 0) { warning("No genes for cell type '", ct, "' — skipping."); next }

    # ── Block-parallel cosinor fitting ────────────────────────────────────────
    block_size   <- ifelse(ncol(X_ct) > 15000, 50, 500)
    gene_blocks  <- split(seq_len(ngenes), ceiling(seq_len(ngenes) / block_size))

    message("  Fitting ", ngenes, " genes in ", length(gene_blocks),
            " blocks of ≤", block_size, " ...")

    results_blocks <- future.apply::future_lapply(gene_blocks, function(g_idx_block) {
      nb <- length(g_idx_block)
      res <- list(
        amp   = numeric(nb), mesor = numeric(nb), acro = numeric(nb),
        pval  = numeric(nb), rho   = numeric(nb), pmac = numeric(nb),
        R0mat = matrix(NA_real_, nrow=nb, ncol=nzts)
      )
      for (i in seq_along(g_idx_block)) {
        gi <- g_idx_block[i]
        Xg_zts <- lapply(zt_present, function(z) {
          idx <- which(zt_num_ct == z)
          if (length(idx) == 0) return(numeric(0))
          as.numeric(X_ct[gi, idx])
        })
        fit <- estimate_phaseR(Xg_zts, zt_present, period12, "Ftest")
        res$amp[i]     <- fit$amp
        res$mesor[i]   <- fit$mesor
        res$acro[i]    <- fit$acrophase
        res$pval[i]    <- fit$p_value
        res$rho[i]     <- fit$rho
        res$pmac[i]    <- fit$p_value_macro
        # Per-ZT means
        for (j in seq_along(zt_present)) {
          idx <- which(zt_num_ct == zt_present[j])
          if (length(idx) > 0) res$R0mat[i, j] <- mean(as.numeric(X_ct[gi, idx]), na.rm=TRUE)
        }
      }
      res
    })

    # Flatten blocks
    amp_v  <- unlist(lapply(results_blocks, `[[`, "amp"))
    mesor_v<- unlist(lapply(results_blocks, `[[`, "mesor"))
    acro_v <- unlist(lapply(results_blocks, `[[`, "acro"))
    pval_v <- unlist(lapply(results_blocks, `[[`, "pval"))
    rho_v  <- unlist(lapply(results_blocks, `[[`, "rho"))
    pmac_v <- unlist(lapply(results_blocks, `[[`, "pmac"))
    R0_mat <- do.call(rbind, lapply(results_blocks, `[[`, "R0mat"))

    # Wrap acrophase into [0, period)
    pd     <- if (period12) 12 else 24
    acro24 <- acro_v %% 24

    # BH correction
    padj_v    <- stats::p.adjust(pval_v, method = "BH")
    padj_mac  <- stats::p.adjust(pmac_v, method = "BH")

    # ── Build T1 (stats table) ────────────────────────────────────────────────
    T1 <- data.frame(
      Genes          = gene_names,
      Amp            = amp_v,
      Abs_Amp        = abs(amp_v),
      Mesor          = mesor_v,
      Acrophase      = acro_v,
      Acrophase_24   = acro24,
      Period         = rep(pd, ngenes),
      pvalue         = pval_v,
      pvalue_adj     = padj_v,
      Sine_corr      = rho_v,
      pvalue_corr    = pmac_v,
      pvalue_adj_corr= padj_mac,
      stringsAsFactors = FALSE
    )

    # Remove rows where fitting failed
    valid <- !is.na(T1$pvalue) & !is.na(T1$pvalue_corr)
    T1    <- T1[valid, ]
    R0_mat<- R0_mat[valid, , drop=FALSE]

    # Sort: pvalue_adj_corr ↑, pvalue_adj ↑, Acrophase_24 ↑, Abs_Amp ↓
    ord  <- order(T1$pvalue_adj_corr, T1$pvalue_adj, T1$Acrophase_24, -T1$Abs_Amp)
    T1   <- T1[ord, ]
    R0_mat <- R0_mat[ord, , drop=FALSE]
    rownames(T1) <- NULL

    # ── Build T2 (per-ZT mean expression) ────────────────────────────────────
    T2 <- data.frame(Genes = T1$Genes, stringsAsFactors=FALSE)
    for (j in seq_along(zt_labels)) T2[[zt_labels[j]]] <- R0_mat[, j]

    # T3: ZT-normalised (divide by first non-zero column)
    T3       <- T2
    ref_col  <- T3[, 2]
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

    # ── Write CSVs ────────────────────────────────────────────────────────────
    fbase <- paste0(combo_name, per_label)
    write.csv(T1,      file.path(ct_outdir, paste0(fbase, "circadian_analysis_all.csv")),      row.names=FALSE)
    write.csv(T2,      file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean.csv")),          row.names=FALSE)
    write.csv(T3,      file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_normalized.csv")),row.names=FALSE)
    if (rm_low_conf) {
      write.csv(T1_conf, file.path(ct_outdir, paste0(fbase, "circadian_analysis_confident.csv")),      row.names=FALSE)
      write.csv(T2_conf, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_confident.csv")),      row.names=FALSE)
      write.csv(T3_conf, file.path(ct_outdir, paste0(fbase, "circadian_ZTs_mean_normalized_confident.csv")), row.names=FALSE)
    }

    # ── Optional heatmap ─────────────────────────────────────────────────────
    if (plot_heat && nrow(T1_conf) > 1) {
      tryCatch(
        generate_heatmap(combo_name, strict=TRUE, custom_name="", circ=FALSE,
                         period12=period12, outdir=ct_outdir),
        error = function(e) message("  Heatmap warning: ", e$message)
      )
    }

    all_results[[combo_name]] <- list(T1 = T1, T2 = T2)
    summary_rows[[combo_name]] <- data.frame(
      CellType            = ct,
      Group               = if (!is.null(grp)) grp else NA_character_,
      ComboName           = combo_name,
      NumCells            = sum(ct_mask),
      NumGenes            = ngenes,
      NumTested           = nrow(T1),
      NumConfident_Both   = sum(conf_both),
      NumNonConfident     = sum(!conf_both),
      stringsAsFactors    = FALSE
    )

    }  # end group loop
  }  # end cell-type loop

  # ── Cross-cell-type summary CSV ──────────────────────────────────────────────
  if (length(summary_rows) > 0) {
    summary_tbl <- do.call(rbind, summary_rows)
    write.csv(summary_tbl, file.path(outdir, paste0("all_cell_types", per_label, "summary_results.csv")),
              row.names=FALSE)
  }

  elapsed <- proc.time()["elapsed"] - tic
  message(sprintf("\n=== TimeSCape complete. Total time: %.1f s ===", elapsed))
  invisible(all_results)
}
