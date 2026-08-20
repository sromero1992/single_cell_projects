# =============================================================================
# doublets.R — dual-caller doublet engine (DoubletFinder + scDblFinder)
# =============================================================================
# Migrated verbatim in logic from 01_process_data_new.R. The only structural
# change: every value that used to be a global CONFIG constant is now read from
# a `p` parameter list built by doublet_params(), so each detector is callable
# in isolation and each sample can carry its own resolved rate.
#
# CONTRACT (unchanged): each detector takes a per-sample Seurat object and
# returns list(calls, score, info) with names = cell barcodes, so the caller
# never has to know which method ran.
# =============================================================================

#' Default doublet-detection parameters
#'
#' One place to set everything the engine needs. Pass a modified copy to
#' [detect_doublets()] / [run_doubletfinder()] / [run_scdblfinder()].
#'
#' @param ... Named overrides for any of the defaults below.
#' @return A named list of parameters.
#' @export
doublet_params <- function(...) {
  p <- list(
    # method + consensus
    method                    = "both",          # "DoubletFinder" | "scDblFinder" | "both"
    consensus_rule            = "union",          # union | intersect | DoubletFinder | scDblFinder
    # rate model (10x linear multiplet model)
    rate                      = "auto",           # "auto" or a fraction in (0,1)
    expected_cells_per_sample = NULL,             # NULL = use measured count
    multiplet_rate_per_1k     = 0.008,            # ~0.8% per 1000 cells (10x)
    rate_floor                = 0.01,
    rate_ceiling              = 0.16,
    rollback_threshold        = 0.30,
    # dim-reduction used by the throwaway DoubletFinder embedding
    n_pcs                     = 30,
    n_variable_features       = 2000,
    # DoubletFinder specifics
    df_pn                     = 0.25,
    df_pk_range_min           = 0.01,
    df_pk_range_max           = 0.30,
    df_pk_fallback            = 0.09,
    df_sct                    = FALSE,
    df_homotypic_adjust       = TRUE,
    df_cluster_res            = 0.8,
    # scDblFinder specifics
    scdbl_dbr_null            = FALSE,            # TRUE => let scDblFinder estimate dbr
    scdbl_score_threshold     = 0.50,
    # plotting
    dpi                       = 300
  )
  ov <- list(...)
  p[names(ov)] <- ov
  p
}

# -----------------------------------------------------------------------------
# Rate resolver — expected doublet rate for ONE sample (10x linear model).
# -----------------------------------------------------------------------------
#' Resolve the expected doublet rate for one sample
#'
#' @param n_cells Recovered cell count for this sample (post pre-QC / DecontX).
#' @param sample_id Character label, used in messages.
#' @param p A [doublet_params()] list.
#' @return list(rate, source, raw, basis_cells).
#' @export
resolve_doublet_rate <- function(n_cells, sample_id = "", p = doublet_params()) {

  if (!identical(p$rate, "auto")) {
    message(sprintf("    -> Doublet rate (fixed): %.2f%%  (expected ~%d doublets)",
                    p$rate * 100, round(n_cells * p$rate)))
    return(list(rate = p$rate, source = "FIXED", raw = p$rate, basis_cells = n_cells))
  }

  if (is.null(p$expected_cells_per_sample)) {
    basis_cells <- n_cells;                        basis_label <- "measured"
  } else {
    basis_cells <- p$expected_cells_per_sample;    basis_label <- "user-stated target"
  }

  raw_rate <- (basis_cells / 1000) * p$multiplet_rate_per_1k
  rate     <- min(max(raw_rate, p$rate_floor), p$rate_ceiling)
  source   <- if (raw_rate < p$rate_floor) "AUTO_FLOORED"
              else if (raw_rate > p$rate_ceiling) "AUTO_CAPPED" else "AUTO"

  message(sprintf("    -> Doublet rate (auto): %s cells (%s) -> %.2f%%%s  (expected ~%d doublets)",
                  basis_cells, basis_label, rate * 100,
                  if (source != "AUTO")
                    sprintf("  [%s; unclamped %.2f%%]", source, raw_rate * 100) else "",
                  round(n_cells * rate)))
  if (source == "AUTO_CAPPED") {
    warning(sprintf("  [RATE CAP] Sample '%s': %s cells implies %.1f%% multiplets; clamped to %.1f%%.",
                    sample_id, basis_cells, raw_rate * 100, p$rate_ceiling * 100))
  }
  list(rate = rate, source = source, raw = raw_rate, basis_cells = basis_cells)
}

# -----------------------------------------------------------------------------
# DETECTOR 1: DoubletFinder (paramSweep + constrained pK + homotypic adjust)
# -----------------------------------------------------------------------------
#' Run DoubletFinder on one per-sample Seurat object
#'
#' @param seu A per-sample Seurat object (post-DecontX / pre-QC).
#' @param sample_id Character label for plot filenames.
#' @param diag_dir Directory for diagnostic plots (NULL to skip plotting).
#' @param dbl_rate Expected doublet rate for this sample (see [resolve_doublet_rate()]).
#' @param p A [doublet_params()] list.
#' @return list(calls, score, info) with names = cell barcodes.
#' @export
run_doubletfinder <- function(seu, sample_id, diag_dir, dbl_rate, p = doublet_params()) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE))
    stop("DoubletFinder is not installed. install_github('chris-mcginnis-ucsf/DoubletFinder').")

  n_cells  <- ncol(seu)
  npcs_use <- min(p$n_pcs, n_cells - 1)
  if (npcs_use < 5) {
    warning(sprintf("  [SKIP] Sample '%s' has only %d cells - too few for DoubletFinder.",
                    sample_id, n_cells))
    return(list(
      calls = stats::setNames(rep("Singlet", n_cells), colnames(seu)),
      score = stats::setNames(rep(NA_real_, n_cells), colnames(seu)),
      info  = list(method = "DoubletFinder", pk = NA, pk_initial = NA, nExp = 0,
                   nExp_adj = 0, rate_used = dbl_rate, homotypic_prop = NA,
                   note = "TOO_FEW_CELLS")))
  }

  # A) throwaway embedding required by DoubletFinder
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                      nfeatures = p$n_variable_features, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = npcs_use, verbose = FALSE)

  # B) paramSweep -> BCmetric curve
  message(sprintf("    -> [DoubletFinder] paramSweep for %s (slow step)...", sample_id))
  sweep.res <- DoubletFinder::paramSweep(seu, PCs = 1:npcs_use, sct = p$df_sct)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))

  # C) constrained pK selection
  initial_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  final_pk   <- initial_pk
  pk_note    <- "GLOBAL_PEAK_IN_RANGE"
  if (final_pk < p$df_pk_range_min || final_pk > p$df_pk_range_max) {
    in_range <- bcmvn[bcmvn$pK >= p$df_pk_range_min & bcmvn$pK <= p$df_pk_range_max, ]
    if (nrow(in_range) > 0 && any(is.finite(in_range$BCmetric))) {
      final_pk <- in_range$pK[which.max(in_range$BCmetric)]; pk_note <- "CONSTRAINED_TO_RANGE"
    } else {
      final_pk <- p$df_pk_fallback; pk_note <- "HARDCODED_FALLBACK"
    }
  }

  # D) diagnostic pK plot
  if (!is.null(diag_dir)) tryCatch({
    sel <- bcmvn[bcmvn$pK == final_pk, ]
    pk_plot <- ggplot2::ggplot(bcmvn, ggplot2::aes(x = .data$pK, y = .data$BCmetric, group = 1)) +
      ggplot2::geom_line(color = "grey60") + ggplot2::geom_point(color = "grey60") +
      ggplot2::annotate("rect", xmin = p$df_pk_range_min, xmax = p$df_pk_range_max,
                        ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#4393C3") +
      ggplot2::geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
      { if (nrow(sel) > 0) ggplot2::geom_point(data = sel,
          ggplot2::aes(x = .data$pK, y = .data$BCmetric), color = "red", size = 4, shape = 18) } +
      ggplot2::ggtitle(paste0("pK Finder: ", sample_id),
                       subtitle = paste0("Selected pK = ", final_pk, " | Global peak = ",
                                         initial_pk, " | ", pk_note)) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(diag_dir, paste0(sample_id, "_DF_pk_plot.png")),
                    pk_plot, width = 7, height = 5, dpi = p$dpi, bg = "white")
  }, error = function(e) message("      -> [WARNING] pK plot failed: ", e$message))

  # E) nExp with optional homotypic adjustment
  nExp_raw <- round(n_cells * dbl_rate); nExp_use <- nExp_raw; homotypic_prop <- NA_real_
  if (p$df_homotypic_adjust) tryCatch({
    seu <- Seurat::FindNeighbors(seu, dims = 1:npcs_use, verbose = FALSE)
    seu <- Seurat::FindClusters(seu, resolution = p$df_cluster_res, verbose = FALSE)
    homotypic_prop <- DoubletFinder::modelHomotypic(seu$seurat_clusters)
    nExp_use <- round(nExp_raw * (1 - homotypic_prop))
  }, error = function(e) message("      -> [WARNING] Homotypic adjustment failed: ", e$message))
  nExp_use <- max(1, min(nExp_use, n_cells - 1))

  # F) run DoubletFinder
  seu <- DoubletFinder::doubletFinder(seu, PCs = 1:npcs_use, pN = p$df_pn,
                                      pK = final_pk, nExp = nExp_use, sct = p$df_sct)
  res_col  <- utils::tail(grep("^DF.classifications", colnames(seu@meta.data), value = TRUE), 1)
  pann_col <- utils::tail(grep("^pANN", colnames(seu@meta.data), value = TRUE), 1)
  calls <- stats::setNames(as.character(seu@meta.data[[res_col]]), rownames(seu@meta.data))
  score <- if (length(pann_col) == 1 && !is.na(pann_col))
    stats::setNames(as.numeric(seu@meta.data[[pann_col]]), rownames(seu@meta.data))
  else stats::setNames(rep(NA_real_, n_cells), rownames(seu@meta.data))

  list(calls = calls, score = score,
       info = list(method = "DoubletFinder", pk = final_pk, pk_initial = initial_pk,
                   nExp = nExp_raw, nExp_adj = nExp_use, rate_used = dbl_rate,
                   homotypic_prop = homotypic_prop, note = pk_note))
}

# -----------------------------------------------------------------------------
# DETECTOR 2: scDblFinder
# -----------------------------------------------------------------------------
#' Run scDblFinder on one per-sample Seurat object
#'
#' @inheritParams run_doubletfinder
#' @return list(calls, score, info) with names = cell barcodes.
#' @export
run_scdblfinder <- function(seu, sample_id, diag_dir, dbl_rate, p = doublet_params()) {
  if (!requireNamespace("scDblFinder", quietly = TRUE))
    stop("scDblFinder is not installed. BiocManager::install('scDblFinder').")

  n_cells <- ncol(seu)
  sce <- Seurat::as.SingleCellExperiment(seu)
  sce <- scDblFinder::scDblFinder(
    sce, dbr = if (p$scdbl_dbr_null) NULL else dbl_rate,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 123), verbose = FALSE)

  score <- stats::setNames(as.numeric(sce$scDblFinder.score), colnames(seu))
  calls <- stats::setNames(
    ifelse(score >= p$scdbl_score_threshold, "Doublet", "Singlet"), colnames(seu))

  if (!is.null(diag_dir)) tryCatch({
    dblt_frac <- sum(calls == "Doublet") / n_cells
    df <- data.frame(score = score, call = calls)
    pp <- ggplot2::ggplot(df, ggplot2::aes(x = .data$score, fill = .data$call)) +
      ggplot2::geom_histogram(bins = 60, color = "black", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = p$scdbl_score_threshold, linetype = "dashed",
                          color = "red", linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = c(Singlet = "#4393C3", Doublet = "#D6604D")) +
      ggplot2::labs(title = paste0("scDblFinder Score: ", sample_id),
                    subtitle = sprintf("%.1f%% doublets | threshold = %s",
                                       dblt_frac * 100, p$scdbl_score_threshold),
                    x = "Doublet Score", y = "Cell Count", fill = "Call") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(diag_dir, paste0(sample_id, "_scDbl_score_distribution.png")),
                    pp, width = 7, height = 5, dpi = p$dpi, bg = "white")
  }, error = function(e) message("      -> [WARNING] Score plot failed: ", e$message))

  rm(sce); gc()
  list(calls = calls, score = score,
       info = list(method = "scDblFinder", pk = NA, pk_initial = NA, nExp = NA,
                   nExp_adj = NA, rate_used = dbl_rate, homotypic_prop = NA,
                   note = paste0("threshold=", p$scdbl_score_threshold)))
}

# -----------------------------------------------------------------------------
# Concordance reporter (used when method == "both")
# -----------------------------------------------------------------------------
#' Confusion matrix + Cohen's kappa between the two callers
#'
#' @param df_calls,sc_calls Named character vectors ("Singlet"/"Doublet").
#' @param sample_id Character label.
#' @param diag_dir Directory for the plot (NULL to skip).
#' @return list(agreement, kappa, n_df_only, n_sc_only, n_both).
#' @export
report_concordance <- function(df_calls, sc_calls, sample_id, diag_dir = NULL) {
  common <- intersect(names(df_calls), names(sc_calls))
  a <- factor(df_calls[common], levels = c("Singlet", "Doublet"))
  b <- factor(sc_calls[common], levels = c("Singlet", "Doublet"))
  tab <- table(DoubletFinder = a, scDblFinder = b)
  n <- length(common); agree <- sum(diag(tab)) / n
  p_exp <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  kappa <- if (p_exp < 1) (agree - p_exp) / (1 - p_exp) else NA_real_

  if (!is.null(diag_dir)) tryCatch({
    tdf <- as.data.frame(tab)
    pc <- ggplot2::ggplot(tdf, ggplot2::aes(x = .data$scDblFinder, y = .data$DoubletFinder,
                                            fill = .data$Freq)) +
      ggplot2::geom_tile(color = "white", linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Freq), size = 6, fontface = "bold") +
      ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#4393C3") +
      ggplot2::labs(title = paste0("Doublet Call Concordance: ", sample_id),
                    subtitle = sprintf("Agreement = %.1f%% | kappa = %.3f | n = %d",
                                       agree * 100, kappa, n)) +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(diag_dir, paste0(sample_id, "_doublet_concordance.png")),
                    pc, width = 6, height = 5, dpi = 300, bg = "white")
  }, error = function(e) message("      -> [WARNING] Concordance plot failed: ", e$message))

  message(sprintf("    -> Concordance: %.1f%% agreement | kappa = %.3f", agree * 100, kappa))
  list(agreement = agree, kappa = kappa,
       n_df_only = tab["Doublet", "Singlet"], n_sc_only = tab["Singlet", "Doublet"],
       n_both = tab["Doublet", "Doublet"])
}

# -----------------------------------------------------------------------------
# Consensus rule: combine two binary call vectors into one.
# -----------------------------------------------------------------------------
.apply_consensus <- function(df_calls, sc_calls, rule) {
  cells <- union(names(df_calls), names(sc_calls))
  d <- df_calls[cells] == "Doublet"; d[is.na(d)] <- FALSE
  s <- sc_calls[cells] == "Doublet"; s[is.na(s)] <- FALSE
  is_dbl <- switch(rule,
    union        = d | s,
    intersect    = d & s,
    DoubletFinder = d,
    scDblFinder  = s,
    d | s)
  stats::setNames(ifelse(is_dbl, "Doublet", "Singlet"), cells)
}

# -----------------------------------------------------------------------------
# detect_doublets — the level-agnostic entry point.
# -----------------------------------------------------------------------------
#' Detect doublets on a sample list or a merged object
#'
#' Writes the schema columns (`DF_score`, `DF_class`, `scDblFinder_score`,
#' `scDblFinder_class`, `doublet_consensus`). Never removes cells — use
#' [apply_doublet_action()] for that.
#'
#' @param x A Seurat object or a list of per-sample Seurat objects.
#' @param method,rate Convenience overrides for `p$method` / `p$rate`.
#' @param by For a merged object: metadata column to split on before detecting
#'   (per-sample detection). `NULL` runs on the whole object at once.
#' @param p A [doublet_params()] list.
#' @param diag_dir Directory for per-sample diagnostic plots (NULL to skip).
#' @return The same type as `x`, with doublet columns added.
#' @export
detect_doublets <- function(x, method = NULL, rate = NULL, by = "SampleID",
                            p = doublet_params(), diag_dir = NULL) {
  if (!is.null(method)) p$method <- method
  if (!is.null(rate))   p$rate   <- rate

  # List input: map per sample (each object is one sample already).
  if (.is_sample_list(x)) {
    return(.map_samples(x, function(seu) {
      sid <- unique(as.character(seu@meta.data[[TamuScDSC_schema()$sample]]))[1] %||% "sample"
      .detect_one(seu, sid, p, diag_dir)
    }))
  }

  # Merged object: split by `by` (per-sample) or run global.
  if (.is_seurat(x)) {
    if (is.null(by) || !by %in% colnames(x@meta.data)) {
      return(.detect_one(x, "ALL", p, diag_dir))
    }
    ids <- unique(as.character(x@meta.data[[by]]))
    res_cols <- list()
    for (sid in ids) {
      cells <- colnames(x)[as.character(x@meta.data[[by]]) == sid]
      seu   <- subset(x, cells = cells)
      seu   <- .detect_one(seu, sid, p, diag_dir)
      res_cols[[sid]] <- seu@meta.data[, unlist(TamuScDSC_schema()[c(
        "df_score","df_class","sc_score","sc_class","consensus")]), drop = FALSE]
    }
    allrows <- do.call(rbind, res_cols)
    allrows <- allrows[match(colnames(x), rownames(allrows)), , drop = FALSE]
    for (cn in colnames(allrows)) x@meta.data[[cn]] <- allrows[[cn]]
    return(.stamp(x, "detect_doublets", list(method = p$method, by = by)))
  }
  stop("detect_doublets(): x must be a Seurat object or a list of them.")
}

# Run detection on ONE object (one sample) and write schema columns onto it.
.detect_one <- function(seu, sample_id, p, diag_dir) {
  s  <- TamuScDSC_schema()
  rr <- resolve_doublet_rate(ncol(seu), sample_id, p)$rate

  df_res <- NULL; sc_res <- NULL
  if (p$method %in% c("DoubletFinder", "both"))
    df_res <- run_doubletfinder(seu, sample_id, diag_dir, rr, p)
  if (p$method %in% c("scDblFinder", "both"))
    sc_res <- run_scdblfinder(seu, sample_id, diag_dir, rr, p)

  cells <- colnames(seu)
  na_v  <- stats::setNames(rep(NA_real_, length(cells)), cells)
  sing  <- stats::setNames(rep("Singlet", length(cells)), cells)

  df_score <- if (!is.null(df_res)) df_res$score[cells] else na_v
  df_class <- if (!is.null(df_res)) df_res$calls[cells] else sing
  sc_score <- if (!is.null(sc_res)) sc_res$score[cells] else na_v
  sc_class <- if (!is.null(sc_res)) sc_res$calls[cells] else sing

  seu@meta.data[[s$df_score]] <- as.numeric(df_score)
  seu@meta.data[[s$df_class]] <- as.character(df_class)
  seu@meta.data[[s$sc_score]] <- as.numeric(sc_score)
  seu@meta.data[[s$sc_class]] <- as.character(sc_class)

  if (!is.null(df_res) && !is.null(sc_res)) {
    report_concordance(df_res$calls, sc_res$calls, sample_id, diag_dir)
    cons <- .apply_consensus(df_res$calls, sc_res$calls, p$consensus_rule)[cells]
  } else if (!is.null(df_res)) {
    cons <- df_class
  } else {
    cons <- sc_class
  }
  seu@meta.data[[s$consensus]] <- as.character(cons)
  .stamp(seu, "detect_doublets", list(method = p$method, rate = rr,
                                       consensus_rule = p$consensus_rule))
}

# -----------------------------------------------------------------------------
# apply_doublet_action — the ONLY function that removes cells.
# -----------------------------------------------------------------------------
#' Label or remove cells based on a doublet decision column
#'
#' @param x A Seurat object or a sample list.
#' @param action "label" (default; no cells removed) or "remove".
#' @param target Which decision column to act on: "consensus" (default),
#'   "DoubletFinder", "scDblFinder", or "cluster_flag".
#' @return `x` with a `Doublet_Status` column, and (if action="remove") the
#'   flagged cells dropped.
#' @export
apply_doublet_action <- function(x, action = c("label", "remove"),
                                 target = c("consensus", "DoubletFinder",
                                            "scDblFinder", "cluster_flag")) {
  action <- match.arg(action); target <- match.arg(target)
  if (.is_sample_list(x))
    return(.map_samples(x, apply_doublet_action, action = action, target = target))

  s <- TamuScDSC_schema()
  col <- switch(target, consensus = s$consensus, DoubletFinder = s$df_class,
                scDblFinder = s$sc_class, cluster_flag = s$cluster_flag)
  if (!col %in% colnames(x@meta.data))
    stop(sprintf("apply_doublet_action(): column '%s' not found. Run the detection step first.", col))

  is_dbl <- if (target == "cluster_flag")
    x@meta.data[[col]] == "flagged_doublet_cluster"
  else x@meta.data[[col]] == "Doublet"
  x$Doublet_Status <- ifelse(is_dbl, "Doublet", "Singlet")

  if (action == "remove") {
    keep <- colnames(x)[!is_dbl]
    n0 <- ncol(x); x <- subset(x, cells = keep)
    message(sprintf("  apply_doublet_action: removed %d of %d cells (target=%s).",
                    n0 - ncol(x), n0, target))
    x <- .stamp(x, "apply_doublet_action", list(action = "remove", target = target))
  } else {
    x <- .stamp(x, "apply_doublet_action", list(action = "label", target = target))
  }
  x
}

# =============================================================================
# CLUSTER-LEVEL DOUBLET REVIEW (merged in from cluster_review.R)
# =============================================================================
# Kept in the same file as the doublet engine so the whole doublet workflow —
# per-sample detection (scores) then cluster-level review — lives in one place.

#' Default cluster-review parameters
#' @param ... Named overrides.
#' @return A named list.
#' @export
cluster_review_params <- function(...) {
  p <- list(
    recluster            = TRUE,
    recluster_resolution = 4.0,
    recluster_reduction  = "harmony",
    recluster_dims       = 50,
    cluster_col          = "clusters_review",
    flag_rule            = "hybrid",   # hybrid | either | both
    weak_k               = 1.5,
    strong_k             = 3.5,
    df_mean_floor        = 0.30,
    scdbl_mean_floor     = 0.30,
    fraction_threshold   = 0.60,
    annotation_resolution = 1.0,
    dpi                  = 300
  )
  ov <- list(...); p[names(ov)] <- ov; p
}

#' Cluster-level doublet review
#'
#' @param obj An integrated Seurat object carrying `DF_score` / `scDblFinder_score`
#'   (and their `_class` columns) from [detect_doublets()].
#' @param p A [cluster_review_params()] list.
#' @param action "flag" (default; add `cluster_dbl_flag`, remove nothing) or
#'   "remove" (drop flagged clusters and re-cluster the clean object).
#' @param out_dir Optional directory for review plots + summary table.
#' @return `obj` with a `cluster_dbl_flag` column (re-clustered if action="remove").
#' @export
cluster_doublet_review <- function(obj, p = cluster_review_params(),
                                   action = c("flag", "remove"), out_dir = NULL) {
  action <- match.arg(action)
  stopifnot(.is_seurat(obj))
  s <- TamuScDSC_schema()
  need <- c(s$df_score, s$sc_score, s$df_class, s$sc_class)
  miss <- setdiff(need, colnames(obj@meta.data))
  if (length(miss) > 0)
    stop("cluster_doublet_review(): missing column(s) ", paste(miss, collapse = ", "),
         ". Run detect_doublets(method = 'both') first.")
  if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # --- clustering to score ---------------------------------------------------
  cc <- p$cluster_col
  if (p$recluster) {
    n_dims <- min(p$recluster_dims,
                  ncol(SeuratObject::Embeddings(obj, p$recluster_reduction)))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n_dims, reduction = p$recluster_reduction,
                                 graph.name = "review_nn", verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = p$recluster_resolution,
                                graph.name = "review_nn", cluster.name = cc,
                                verbose = FALSE)
  }
  obj@meta.data[[cc]] <- factor(as.character(obj@meta.data[[cc]]))

  # --- per-cluster stats -----------------------------------------------------
  md <- obj@meta.data
  summ <- dplyr::summarise(
    dplyr::group_by(md, cluster = .data[[cc]]),
    N_Cells        = dplyr::n(),
    DF_mean        = mean(.data[[s$df_score]], na.rm = TRUE),
    scDbl_mean     = mean(.data[[s$sc_score]], na.rm = TRUE),
    DF_frac_dbl    = mean(.data[[s$df_class]] == "Doublet", na.rm = TRUE),
    scDbl_frac_dbl = mean(.data[[s$sc_class]] == "Doublet", na.rm = TRUE),
    .groups = "drop")

  # --- adaptive cutoffs: median + k*MAD, floored -----------------------------
  madcut <- function(x, k, fl) max(stats::median(x, na.rm = TRUE) +
                                    k * stats::mad(x, na.rm = TRUE), fl)
  DF_W <- madcut(summ$DF_mean,    p$weak_k,   p$df_mean_floor)
  SC_W <- madcut(summ$scDbl_mean, p$weak_k,   p$scdbl_mean_floor)
  DF_S <- madcut(summ$DF_mean,    p$strong_k, p$df_mean_floor)
  SC_S <- madcut(summ$scDbl_mean, p$strong_k, p$scdbl_mean_floor)
  message(sprintf("  Cutoffs: weak(DF>=%.2f, scDbl>=%.2f) strong(DF>=%.2f, scDbl>=%.2f) rule='%s'",
                  DF_W, SC_W, DF_S, SC_S, p$flag_rule))

  ft <- p$fraction_threshold
  summ <- dplyr::mutate(summ,
    df_weak   = .data$DF_mean    >= DF_W | .data$DF_frac_dbl    >= ft,
    sc_weak   = .data$scDbl_mean >= SC_W | .data$scDbl_frac_dbl >= ft,
    df_strong = .data$DF_mean    >= DF_S,
    sc_strong = .data$scDbl_mean >= SC_S,
    flagged   = dplyr::case_when(
      p$flag_rule == "both"   ~ .data$df_weak & .data$sc_weak,
      p$flag_rule == "either" ~ .data$df_weak | .data$sc_weak,
      TRUE                    ~ (.data$df_strong | .data$sc_strong) |
                                (.data$df_weak & .data$sc_weak)),
    reason = dplyr::case_when(
      !.data$flagged                       ~ "",
      .data$df_strong & !.data$sc_strong   ~ "DoubletFinder extreme",
      .data$sc_strong & !.data$df_strong   ~ "scDblFinder extreme",
      .data$df_strong & .data$sc_strong    ~ "both extreme",
      .data$df_weak & .data$sc_weak        ~ "both moderate",
      .data$df_weak                        ~ "DoubletFinder only",
      TRUE                                 ~ "scDblFinder only"))

  flagged_clusters <- as.character(summ$cluster[summ$flagged])
  obj@meta.data[[s$cluster_flag]] <- ifelse(
    as.character(obj@meta.data[[cc]]) %in% flagged_clusters,
    "flagged_doublet_cluster", "keep")
  n_flag <- sum(obj@meta.data[[s$cluster_flag]] == "flagged_doublet_cluster")
  message(sprintf("  rule='%s' -> %d cluster(s) flagged (%d cells, %.1f%%): %s",
                  p$flag_rule, length(flagged_clusters), n_flag,
                  100 * n_flag / ncol(obj), paste(flagged_clusters, collapse = ", ")))

  # --- optional plots + table ------------------------------------------------
  if (!is.null(out_dir)) tryCatch(
    .cluster_review_plots(obj, summ, cc, DF_W, SC_W, DF_S, SC_S,
                          flagged_clusters, p, out_dir),
    error = function(e) message("  [WARNING] Review plots failed: ", e$message))

  obj <- .stamp(obj, "cluster_doublet_review",
                list(flag_rule = p$flag_rule, flagged = flagged_clusters))

  # --- optional remove + re-cluster ------------------------------------------
  if (action == "remove" && length(flagged_clusters) > 0) {
    obj <- subset(obj, cells = colnames(obj)[
      obj@meta.data[[s$cluster_flag]] == "keep"])
    n_dims <- min(p$recluster_dims, ncol(SeuratObject::Embeddings(obj, "harmony")))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n_dims, reduction = "harmony",
                                 graph.name = "harmony_nn", verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = p$annotation_resolution,
                                graph.name = "harmony_nn",
                                cluster.name = "clusters_harmony", verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, dims = 1:n_dims, reduction = "harmony",
                           reduction.name = "umap_harmony", n.epochs = 500,
                           verbose = FALSE)
    message(sprintf("  Removed flagged clusters -> %d cells; re-clustered at %.1f.",
                    ncol(obj), p$annotation_resolution))
    obj <- .stamp(obj, "cluster_doublet_review:remove",
                  list(removed_clusters = flagged_clusters))
  }
  obj
}

# Internal: the three review figures + xlsx summary.
.cluster_review_plots <- function(obj, summ, cc, DF_W, SC_W, DF_S, SC_S,
                                  flagged_clusters, p, out_dir) {
  s <- TamuScDSC_schema()
  lev <- levels(summ$cluster); num <- suppressWarnings(as.numeric(lev))
  ord <- if (!any(is.na(num))) lev[order(num)] else lev
  bar_df <- tidyr::pivot_longer(
    dplyr::select(summ, cluster, DoubletFinder = .data$DF_mean,
                  scDblFinder = .data$scDbl_mean),
    -cluster, names_to = "Method", values_to = "Mean_Score")
  bar_df$cluster <- factor(as.character(bar_df$cluster), levels = ord)

  p_bar <- ggplot2::ggplot(bar_df, ggplot2::aes(.data$cluster, .data$Mean_Score,
                                                fill = .data$Method)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::geom_hline(yintercept = DF_W, linetype = "dashed", color = "#F8766D") +
    ggplot2::geom_hline(yintercept = SC_W, linetype = "dashed", color = "#00BFC4") +
    ggplot2::geom_hline(yintercept = DF_S, linetype = "dotted", color = "#F8766D") +
    ggplot2::geom_hline(yintercept = SC_S, linetype = "dotted", color = "#00BFC4") +
    ggplot2::scale_fill_manual(values = c(DoubletFinder = "#F8766D",
                                          scDblFinder = "#00BFC4")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Mean doublet score per cluster",
                  subtitle = paste0("Dashed = moderate (weak_k), dotted = extreme (strong_k). Flagged: ",
                                    if (length(flagged_clusters)) paste(flagged_clusters, collapse = ", ") else "none"),
                  x = cc, y = "Mean score") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(out_dir, "01_mean_doublet_score_per_cluster.png"),
                  p_bar, width = 13, height = 5, dpi = p$dpi, bg = "white")

  red <- if ("umap_harmony" %in% names(obj@reductions)) "umap_harmony" else "umap"

  # 02: continuous per-cell doublet scores on the UMAP (both callers side by side)
  p_scores <- Seurat::FeaturePlot(
    obj, features = c(s$df_score, s$sc_score), reduction = red,
    order = TRUE, pt.size = 0.2, ncol = 2, cols = c("grey85", "#B2182B"))
  ggplot2::ggsave(file.path(out_dir, "02_umap_doublet_scores.png"),
                  p_scores, width = 14, height = 6, dpi = p$dpi, bg = "white")

  # 03: which clusters got flagged
  p_um <- Seurat::DimPlot(obj, group.by = s$cluster_flag, reduction = red,
                          cols = c(keep = "grey85",
                                   flagged_doublet_cluster = "#B2182B"))
  ggplot2::ggsave(file.path(out_dir, "03_umap_flagged_clusters.png"),
                  p_um, width = 8, height = 7, dpi = p$dpi, bg = "white")

  if (requireNamespace("writexl", quietly = TRUE)) {
    out_tab <- dplyr::arrange(
      dplyr::transmute(summ, Cluster = as.character(cluster), .data$N_Cells,
                       DF_mean = round(.data$DF_mean, 4),
                       scDbl_mean = round(.data$scDbl_mean, 4),
                       DF_pct_doublet = round(100 * .data$DF_frac_dbl, 1),
                       scDbl_pct_doublet = round(100 * .data$scDbl_frac_dbl, 1),
                       Flagged = .data$flagged, Reason = .data$reason),
      dplyr::desc(.data$Flagged), dplyr::desc(pmax(.data$DF_mean, .data$scDbl_mean)))
    writexl::write_xlsx(list(Per_Cluster = out_tab),
                        file.path(out_dir, "cluster_doublet_summary.xlsx"))
  }
}
