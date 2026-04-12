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
#' @importFrom stats cor.test ecdf qf qchisq
#' @export
#'
estimate_phaseR <- function(Xg_zts, actual_times, period12, test_type) {

  # Initialize return structure with NAs
  init_result <- function() {
    list(
      acrophase = NA_real_,
      amp = NA_real_,
      period = NA_real_,
      mesor = NA_real_,
      p_value = NA_real_,
      rho = NA_real_,
      p_value_macro = NA_real_
    )
  }

  tryCatch({
    # Remove empty cells
    non_empty <- sapply(Xg_zts, function(x) {
      !is.null(x) && length(x) > 0
    })

    Xg_zts <- Xg_zts[non_empty]
    actual_times <- actual_times[non_empty]
    nzts <- length(Xg_zts)

    # Check minimum number of time points
    if (nzts < 4) {
      return(init_result())
    }

    # Flatten expression values and create time grid
    icells <- sapply(Xg_zts, length)
    num_cells <- sum(icells)

    # Build flattened vectors
    R <- numeric(num_cells)
    time_grid <- numeric(num_cells)
    R0 <- numeric(nzts)  # Per-ZT means

    ic <- 0
    max_R0 <- -Inf
    max_peak_t <- actual_times[1]

    for (it in seq_len(nzts)) {
      n_it <- icells[it]
      cell_vals <- as.numeric(Xg_zts[[it]])
      R[(ic + 1):(ic + n_it)] <- cell_vals

      meanval <- mean(cell_vals, na.rm = TRUE)
      R0[it] <- meanval
      time_grid[(ic + 1):(ic + n_it)] <- actual_times[it]

      ic <- ic + n_it

      if (meanval > max_R0) {
        max_R0 <- meanval
        max_peak_t <- actual_times[it]
      }
    }

    # Null model: mean
    mean_model <- mean(R, na.rm = TRUE)
    SSR_null <- sum((R - mean_model)^2, na.rm = TRUE)

    # Set period
    period <- if (period12) 12 else 24

    # Initial guess for amplitude
    max_amp_guess <- max_R0 - mean_model

    # Define cosine model
    # Formula: amp * cos(2*pi*(t - acro)/period) + mesor
    cosine_formula <- formula(
      R ~ amp * cos(2 * pi * (time_grid - acro) / period) + mesor
    )

    # Fit using nonlinear least squares (Levenberg-Marquardt via minpack.lm)
    fit_data <- data.frame(R = R, time_grid = time_grid)

    start_params <- list(
      acro = max_peak_t,
      amp = max_amp_guess,
      mesor = mean_model
    )

    fit_obj <- minpack.lm::nlsLM(
      cosine_formula,
      data = fit_data,
      start = start_params,
      lower = c(acro = 0, amp = -Inf, mesor = 0),
      upper = c(acro = period, amp = Inf, mesor = Inf),
      control = minpack.lm::nls.lm.control(maxiter = 100)
    )

    # Extract fitted parameters
    acro_fit <- coef(fit_obj)["acro"]
    amp_fit <- coef(fit_obj)["amp"]
    mesor_fit <- coef(fit_obj)["mesor"]

    # Compute SSR for sine model
    fitted_vals <- predict(fit_obj)
    SSR_sine <- sum((R - fitted_vals)^2, na.rm = TRUE)

    # Perform significance test
    if (test_type == "Ftest") {
      d1 <- 2
      d2 <- num_cells - 3

      if (d2 > 0) {
        F_stat <- ((SSR_null - SSR_sine) / d1) / (SSR_sine / d2)
        p_value <- 1 - stats::pf(F_stat, d1, d2)
      } else {
        p_value <- NA_real_
      }

    } else if (test_type == "LRT") {
      s2_null <- SSR_null / num_cells
      logL_null <- -0.5 * num_cells * (log(2 * pi * s2_null) + 1)

      s2_sine <- SSR_sine / num_cells
      logL_sine <- -0.5 * num_cells * (log(2 * pi * s2_sine) + 1)

      lrt_stat <- -2 * (logL_null - logL_sine)
      p_value <- 1 - stats::pchisq(lrt_stat, df = 2)

    } else {
      p_value <- NA_real_
    }

    # Compute macro-level correlation: per-ZT means vs fitted cosine at those times
    fval <- amp_fit * cos(2 * pi * (actual_times - acro_fit) / period) + mesor_fit

    cor_test <- stats::cor.test(R0, fval, method = "pearson")
    rho <- cor_test$estimate
    p_value_macro <- cor_test$p.value

    list(
      acrophase = acro_fit,
      amp = amp_fit,
      period = period,
      mesor = mesor_fit,
      p_value = p_value,
      rho = rho,
      p_value_macro = p_value_macro
    )

  }, error = function(e) {
    # Return NA-filled list on any error
    init_result()
  })
}
