#' Fit Cosine Model and Test Significance for Single Cell Circadian Rhythms
#'
#' Fits a cosine model to expression data from multiple time points using
#' nonlinear least squares and performs statistical significance testing
#' via F-test or likelihood ratio test.
#'
#' @param Xg_zts Named list of numeric vectors, each element contains expression
#'   values for cells measured at a specific ZT time point.
#' @param actual_times Numeric vector of true ZT hours (must match names of Xg_zts).
#' @param period12 Logical. If TRUE, use 12-hour period; if FALSE, use 24-hour period.
#' @param test_type Character: "Ftest" for F-test or "LRT" for likelihood ratio test.
#'
#' @return A list with elements:
#'   - acrophase: peak time of cosine (in hours, 0 to period)
#'   - amp: amplitude of cosine fit
#'   - period: period used (12 or 24)
#'   - mesor: mean level (midline estimating statistic of rhythm)
#'   - p_value: significance test p-value
#'   - rho: Pearson correlation between per-ZT means and fitted cosine
#'   - p_value_macro: p-value for the correlation
#'
#' @importFrom minpack.lm nlsLM
#' @importFrom stats cor.test pf pchisq
#' @export
#'
estimate_phaseR <- function(Xg_zts, actual_times, period12, test_type) {

  # Initialize return structure with NAs
  init_result <- function() {
    list(
      acrophase     = NA_real_,
      amp           = NA_real_,
      period        = NA_real_,
      mesor         = NA_real_,
      p_value       = NA_real_,
      rho           = NA_real_,
      p_value_macro = NA_real_
    )
  }

  tryCatch({

    # ── 1. Filter empty time points ─────────────────────────────────────────────
    non_empty    <- vapply(Xg_zts, function(x) !is.null(x) && length(x) > 0, logical(1))
    Xg_zts       <- Xg_zts[non_empty]
    actual_times <- actual_times[non_empty]
    nzts         <- length(Xg_zts)

    if (nzts < 4) return(init_result())

    # ── 2. Flatten expression vectors ──────────────────────────────────────────
    icells    <- vapply(Xg_zts, length, integer(1))
    num_cells <- sum(icells)

    R         <- numeric(num_cells)
    time_grid <- numeric(num_cells)
    R0        <- numeric(nzts)

    ic        <- 0L
    max_R0    <- -Inf
    max_peak_t <- actual_times[1]

    for (it in seq_len(nzts)) {
      n_it      <- icells[it]
      cell_vals <- as.numeric(Xg_zts[[it]])
      R[(ic + 1L):(ic + n_it)]         <- cell_vals
      time_grid[(ic + 1L):(ic + n_it)] <- actual_times[it]

      meanval <- mean(cell_vals, na.rm = TRUE)
      R0[it]  <- meanval
      ic      <- ic + n_it

      if (is.finite(meanval) && meanval > max_R0) {
        max_R0     <- meanval
        max_peak_t <- actual_times[it]
      }
    }

    # Guard: if all cells have NaN/NA expression, bail out early
    if (!any(is.finite(R))) return(init_result())

    # ── 3. Null model (flat mean) ───────────────────────────────────────────────
    mean_model <- mean(R, na.rm = TRUE)
    if (!is.finite(mean_model)) return(init_result())
    SSR_null   <- sum((R - mean_model)^2, na.rm = TRUE)

    # ── 4. Cosine fit via Levenberg-Marquardt ──────────────────────────────────
    period        <- if (period12) 12 else 24
    max_amp_guess <- max_R0 - mean_model

    start_params <- list(
      acro  = max_peak_t,
      amp   = max_amp_guess,
      mesor = mean_model
    )

    # Provide data explicitly so nlsLM finds R and time_grid regardless of
    # how the package resolves the parent frame across different R/minpack versions.
    fit_data <- data.frame(R = R, time_grid = time_grid)

    fit_obj <- minpack.lm::nlsLM(
      R ~ amp * cos(2 * pi * (time_grid - acro) / period) + mesor,
      data    = fit_data,
      start   = start_params,
      lower   = c(acro = 0,      amp = -Inf, mesor = 0),
      upper   = c(acro = period,  amp =  Inf, mesor = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 100)
    )

    # ── 5. Extract parameters ──────────────────────────────────────────────────
    cf        <- coef(fit_obj)
    acro_fit  <- cf["acro"]
    amp_fit   <- cf["amp"]
    mesor_fit <- cf["mesor"]

    # SSR via deviance() — reads directly from the nlsLM object without
    # calling predict(), saving one large vector allocation per gene.
    SSR_sine  <- deviance(fit_obj)

    # ── 6. Significance test ───────────────────────────────────────────────────
    if (test_type == "Ftest") {
      d1 <- 2L
      d2 <- num_cells - 3L
      p_value <- if (d2 > 0L && is.finite(SSR_sine) && SSR_sine > 0) {
        F_stat <- ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2)
        1 - stats::pf(F_stat, d1, d2)
      } else if (d2 > 0L && is.finite(SSR_sine) && SSR_sine == 0 && SSR_null == 0) {
        # Perfectly flat gene — no circadian rhythm; return p=1
        1
      } else {
        NA_real_
      }

    } else if (test_type == "LRT") {
      s2_null   <- SSR_null / num_cells
      logL_null <- -0.5 * num_cells * (log(2 * pi * s2_null) + 1)
      s2_sine   <- SSR_sine / num_cells
      logL_sine <- -0.5 * num_cells * (log(2 * pi * s2_sine) + 1)
      lrt_stat  <- -2 * (logL_null - logL_sine)
      p_value   <- 1 - stats::pchisq(lrt_stat, df = 2L)

    } else {
      p_value <- NA_real_
    }

    # ── 7. Macro-level Pearson correlation (per-ZT means vs fitted cosine) ─────
    fval      <- amp_fit * cos(2 * pi * (actual_times - acro_fit) / period) + mesor_fit
    cor_test  <- stats::cor.test(R0, fval, method = "pearson")

    list(
      acrophase     = acro_fit,
      amp           = amp_fit,
      period        = period,
      mesor         = mesor_fit,
      p_value       = p_value,
      rho           = cor_test$estimate,
      p_value_macro = cor_test$p.value
    )

  }, error = function(e) {
    # Return NA-filled result on any fitting error.
    # To diagnose what's failing, call diagnose_phaseR() on a few genes first.
    init_result()
  })
}


#' Diagnose estimate_phaseR failures for a cell type
#'
#' Runs estimate_phaseR on the first \code{n} genes for a given cell type and
#' prints whether each succeeded or failed, along with the error message if any.
#' Use this when run_timescape() reports "All genes returned NA" to find the
#' root cause.
#'
#' @param obj         Seurat or SCE object.
#' @param celltype_col Metadata column for cell type.
#' @param zt_col      Metadata column for ZT strings.
#' @param tmeta       data.frame from build_tmeta().
#' @param target_ct   Cell type to test (must match celltype_col values).
#' @param norm_str    Normalisation: "logcounts", "lib_size", or "none".
#' @param period12    Logical. 12-hr (TRUE) or 24-hr (FALSE).
#' @param n           Number of genes to test (default 5).
#'
#' @return Invisibly: list of fit results (or error strings) for each gene.
#' @export
diagnose_phaseR <- function(
    obj,
    celltype_col,
    zt_col,
    tmeta,
    target_ct,
    norm_str = "logcounts",
    period12 = FALSE,
    n        = 5L
) {
  meta      <- .get_meta(obj)
  ct_mask   <- as.character(meta[[celltype_col]]) == target_ct
  zt_strs   <- as.character(meta[[zt_col]])[ct_mask]
  zt_num    <- tmeta$ZT_times[match(zt_strs, tmeta$zt_str)]
  zt_present <- sort(unique(zt_num[!is.na(zt_num)]))

  use_norm  <- norm_str == "logcounts"
  X_raw     <- .get_matrix(obj, use_normalized = use_norm)

  if (!use_norm && norm_str == "lib_size") {
    X_ct_raw <- X_raw[, ct_mask, drop = FALSE]
    cs       <- Matrix::colSums(X_ct_raw)
    cs[cs == 0] <- 1
    X_ct     <- log1p(Matrix::t(Matrix::t(X_ct_raw) * (1e4 / cs)))
  } else {
    X_ct     <- X_raw[, ct_mask, drop = FALSE]
  }

  gene_names <- rownames(X_ct)
  n_test     <- min(n, nrow(X_ct))

  cat(sprintf("Diagnosing estimate_phaseR for '%s' — first %d genes\n",
              target_ct, n_test))
  cat(sprintf("  Cells: %d | ZT points: %d\n", sum(ct_mask), length(zt_present)))

  results <- vector("list", n_test)
  for (i in seq_len(n_test)) {
    Xg_zts <- lapply(zt_present, function(z) {
      idx <- which(zt_num == z)
      if (length(idx) == 0L) return(numeric(0))
      as.numeric(X_ct[i, idx])
    })

    # Run WITHOUT tryCatch to see the actual error
    fit <- tryCatch(
      estimate_phaseR(Xg_zts, zt_present, period12, "Ftest"),
      error = function(e) {
        cat(sprintf("  Gene %d (%s): ERROR — %s\n", i, gene_names[i], e$message))
        NULL
      }
    )

    if (!is.null(fit)) {
      status <- if (!is.na(fit$p_value)) "OK" else "NA (fitting failed)"
      cat(sprintf("  Gene %d (%s): %s  [amp=%.4f, mesor=%.4f, p=%.3g]\n",
                  i, gene_names[i], status,
                  ifelse(is.na(fit$amp), 0, fit$amp),
                  ifelse(is.na(fit$mesor), 0, fit$mesor),
                  ifelse(is.na(fit$p_value), NA, fit$p_value)))
    }
    results[[i]] <- fit
  }

  # Also try a gene that's actually expressed
  expr_row_sums <- Matrix::rowSums(X_ct > 0)
  top_gene_idx  <- which.max(expr_row_sums)
  cat(sprintf("\n  Most-expressed gene: %s (%d / %d cells non-zero)\n",
              gene_names[top_gene_idx], expr_row_sums[top_gene_idx], sum(ct_mask)))

  Xg_top <- lapply(zt_present, function(z) {
    idx <- which(zt_num == z)
    if (length(idx) == 0L) return(numeric(0))
    as.numeric(X_ct[top_gene_idx, idx])
  })

  fit_top <- tryCatch({
    # Run directly without the outer tryCatch — shows real error
    non_empty    <- vapply(Xg_top, function(x) !is.null(x) && length(x) > 0, logical(1))
    Xg_sub       <- Xg_top[non_empty]
    times_sub    <- zt_present[non_empty]
    icells       <- vapply(Xg_sub, length, integer(1))
    num_cells    <- sum(icells)
    R         <- numeric(num_cells); time_grid <- numeric(num_cells); ic <- 0L
    R0 <- numeric(length(times_sub)); max_R0 <- -Inf; max_peak_t <- times_sub[1]
    for (it in seq_along(times_sub)) {
      n_it <- icells[it]; vals <- as.numeric(Xg_sub[[it]])
      R[(ic+1L):(ic+n_it)] <- vals; time_grid[(ic+1L):(ic+n_it)] <- times_sub[it]
      mv <- mean(vals, na.rm=TRUE); R0[it] <- mv; ic <- ic + n_it
      if (is.finite(mv) && mv > max_R0) { max_R0 <- mv; max_peak_t <- times_sub[it] }
    }
    mean_model <- mean(R, na.rm=TRUE)
    period     <- if (period12) 12 else 24
    cat(sprintf("  R range: [%.4f, %.4f]  mean=%.4f  amp_guess=%.4f\n",
                min(R), max(R), mean_model, max_R0 - mean_model))
    fit_data   <- data.frame(R = R, time_grid = time_grid)
    sp         <- list(acro = max_peak_t, amp = max_R0 - mean_model, mesor = mean_model)
    fo         <- minpack.lm::nlsLM(
      R ~ amp * cos(2*pi*(time_grid - acro)/period) + mesor,
      data = fit_data, start = sp,
      lower = c(acro=0, amp=-Inf, mesor=0),
      upper = c(acro=period, amp=Inf, mesor=Inf),
      control = minpack.lm::nls.lm.control(maxiter=100)
    )
    cat(sprintf("  NLS fit OK: amp=%.4f acro=%.2f mesor=%.4f  deviance=%.6f\n",
                coef(fo)["amp"], coef(fo)["acro"], coef(fo)["mesor"], deviance(fo)))
    "OK"
  }, error = function(e) {
    cat(sprintf("  NLS fit ERROR: %s\n", e$message))
    e$message
  })

  invisible(results)
}
