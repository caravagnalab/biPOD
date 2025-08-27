
#' Plot growth model selection results
#'
#' Visualizes the model comparison (BIC or LOO) scores for a set of candidate growth models
#' and highlights the best-performing model.
#'
#' @param x A list containing:
#' \describe{
#'   \item{criterion}{Model selection criterion, e.g., `"bic"` or `"loo"`.}
#'   \item{model_table}{A data frame with model names as row names and a column for the criterion values.}
#' }
#'
#' @return A `ggplot` object showing model scores and highlighting the best model.
#'
#' @export
plot_growth_model_selection <- function(x) {

  # Validate inputs
  if (!is.list(x) || !all(c("criterion", "model_table") %in% names(x))) {
    stop("Input must be a list with 'criterion' and 'model_table' elements")
  }

  criterion <- x$criterion
  model_table <- x$model_table
  rownames(model_table) <- x$model_table$model

  if (!criterion %in% c("bic", "BIC", "loo", "LOO")) {
    stop("Criterion must be 'bic', 'BIC', 'loo', or 'LOO'")
  }

  if (criterion == "loo") {model_table$loo = model_table$looic}

  if (!is.data.frame(model_table)) model_table <- as.data.frame(model_table)

  criterion_col <- NULL
  for (col_name in names(model_table)) {
    if (toupper(col_name) == toupper(criterion)) {
      criterion_col <- col_name
      break
    }
  }

  if (is.null(criterion_col)) stop(paste("Could not find", toupper(criterion), "column in model_table"))

  df = data.frame(
    model = factor(rownames(model_table), levels = rownames(model_table)),
    score = model_table[[criterion_col]],
    QC = model_table$qc,
    is_best = model_table$model == x$best_model,
    stringsAsFactors = FALSE
  )

  y_label <- ifelse(tolower(criterion) == "bic", "BIC", "LOOIC")

  ggplot2::ggplot(df, ggplot2::aes(x = model, y = score, col = QC)) +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "grey60") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_point(
      data = dplyr::filter(df, is_best),
      ggplot2::aes(x = model, y = score),
      size = 5, shape = 1, color = "black", stroke = 1.2
    ) +
    ggplot2::labs(x = "Model", y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("FAIL"="indianred", "PASS"="forestgreen"))
}

#' Plot fitted growth model with uncertainty intervals
#'
#' Plots the fitted growth model along with uncertainty bands and observed data points.
#' Supports exponential, powerlaw, logistic, and Gompertz models.
#'
#' @param x A list containing:
#' \describe{
#'   \item{fit}{Posterior draws from the fitted model.}
#'   \item{breakpoints}{Numeric vector of breakpoints used in the fit.}
#'   \item{best_model}{Character string specifying the model name.}
#' }
#'
#' @param data Data frame containing observed data with time and count columns.
#' @param time_grid Optional numeric vector of time points at which to evaluate predictions.
#' @param time_col Column name for time (default `"time"`).
#' @param count_col Column name for counts (default `"count"`).
#' @param color Line and ribbon color (default `"dodgerblue"`).
#' @param alpha Transparency for the ribbon (default `0.3`).
#' @param CI Credible interval width (default `0.9`).
#'
#' @return A `ggplot` object showing the model fit, credible intervals, and observed data.
#'
#' @export
plot_growth_fit <- function(x, data,
                            time_grid = NULL, time_col = "time", count_col = "count",
                            color = "black", alpha = 0.3, CI = 0.95) {

  ribbon_df = get_data_for_growth_plot(x = x, data = data, time_grid = time_grid, CI = CI, time_col = time_col, count_col = count_col)
  breakpoints <- x$breakpoints
  model_name <- x$best_model

  p = ggplot2::ggplot(ribbon_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = color, alpha = alpha) +
    ggplot2::geom_line(ggplot2::aes(y = median), color = color) +
    ggplot2::geom_point(data = data, ggplot2::aes_string(x = time_col, y = count_col), color = "black", size = 2) +
    ggplot2::geom_vline(xintercept = breakpoints, linetype = "dashed", color = "gray") +
    ggplot2::labs(title = paste("Growth fit:", model_name), x = "Time", y = "Count") +
    ggplot2::theme_bw()

  if (length(breakpoints) != 0) {
    p = add_breakpoint_shadows(p, breakpoints)
  }

  p
}


old_get_data_for_growth_plot = function(x, data, time_grid = NULL, CI = 0.9,
                                    time_col = "time", count_col = "count") {
  fit <- x$fit
  breakpoints <- x$breakpoints
  draws <- posterior::as_draws_matrix(fit$draws)
  model_name <- x$best_model

  # normalize model name to be safe
  model_key <- tolower(model_name)
  if (model_key == "exp_quadratic") model_key <- "quadraticexp"

  time_obs <- data[[time_col]]
  count_obs <- data[[count_col]]
  if (is.null(time_grid)) time_grid <- seq(min(time_obs), max(time_obs), length.out = 200)

  G <- length(breakpoints) + 1
  segment_edges <- c(-Inf, breakpoints, Inf)

  rho_idx <- grep("^rho\\[", colnames(draws))
  rho_mat <- draws[, rho_idx, drop = FALSE]

  has_t0 <- "t0" %in% colnames(draws)
  has_n0 <- "n0" %in% colnames(draws)
  has_K  <- "K"  %in% colnames(draws)
  has_b  <- "b"  %in% colnames(draws)

  t0_vec <- if (has_t0) draws[, "t0"] else rep(min(time_obs), nrow(draws))
  n0_vec <- if (has_n0) draws[, "n0"] else rep(1, nrow(draws))
  K_vec  <- if (has_K)  draws[, "K"]  else rep(NA_real_, nrow(draws))
  b_vec  <- if (has_b)  draws[, "b"]  else rep(NA_real_, nrow(draws))

  # piecewise integral of rho from t0 to t
  integrate_r_matrix <- function(t_grid, t0_vec, rho_mat) {
    n_draws <- nrow(rho_mat)
    n_times <- length(t_grid)
    r_mat <- matrix(0, nrow = n_draws, ncol = n_times)

    for (i in seq_along(t_grid)) {
      t_current <- t_grid[i]
      for (g in seq_len(G)) {
        seg_start <- if (g == 1) segment_edges[1] else segment_edges[g]
        seg_end   <- segment_edges[g + 1]
        if (seg_start >= t_current) next

        int_start <- if (g == 1) pmax(t0_vec, seg_start) else rep(seg_start, n_draws)
        int_end   <- rep(pmin(t_current, seg_end), n_draws)

        valid_idx <- int_start < int_end
        if (!any(valid_idx)) next

        rho_vec <- rho_mat[, g]
        contrib <- if (model_key == "powerlaw") {
          # special case: ∫ rho(t)/t dt -> rho * log(t_end/t_start)
          ratio <- int_end / pmax(int_start, 1e-8)
          rho_vec * log(pmax(ratio, 1e-8))
        } else {
          # default: ∫ rho dt over segment
          rho_vec * (int_end - int_start)
        }

        contrib[!valid_idx] <- 0
        r_mat[, i] <- r_mat[, i] + contrib
      }
    }
    r_mat
  }

  rint_mat <- integrate_r_matrix(time_grid, t0_vec, rho_mat)

  # convenience matrices
  n_draws <- nrow(rho_mat)
  n_times <- length(time_grid)
  n0_mat  <- matrix(n0_vec, nrow = n_draws, ncol = n_times)
  K_mat   <- matrix(K_vec,  nrow = n_draws, ncol = n_times)
  b_mat   <- matrix(b_vec,  nrow = n_draws, ncol = n_times)

  # dt = t - t0 (per-draw)
  dt_mat <- sweep(matrix(rep(time_grid, each = n_draws), nrow = n_draws, ncol = n_times),
                  1, t0_vec, FUN = "-")

  # predictions by model
  mean_mat <- switch(model_key,
                     exponential = {
                       n0_mat * exp(rint_mat)
                     },
                     powerlaw = {
                       n0_mat * exp(rint_mat)
                     },
                     logistic = {
                       n0m <- pmax(n0_mat, 1e-6)
                       Km  <- pmax(K_mat, n0m + 1e-3)  # ensure K > n0
                       ratio  <- (Km - n0m) / n0m
                       denom  <- pmax(1 + ratio * exp(-rint_mat), 1e-6)
                       pred   <- Km / denom
                       pred[!is.finite(pred)] <- NA_real_
                       pred
                     },
                     gompertz = {
                       Km  <- pmax(K_mat, 1e-6)
                       n0m <- pmax(n0_mat, 1e-12) # avoid loglog issues
                       loglog_n0 <- log(-log(pmin(n0m / Km, 1 - 1e-12)))
                       exponent  <- loglog_n0 - rint_mat
                       Km * exp(-exp(exponent))
                     },
                     monomolecular = {
                       # y(t) = K - (K - n0) * exp(-∫rho)
                       if (!has_K) stop("Monomolecular model requires parameter K in draws.")
                       Km  <- pmax(K_mat, 1e-8)
                       n0m <- pmax(n0_mat, 1e-12)
                       Km - (Km - n0m) * exp(-rint_mat)
                     },
                     quadraticexp = ,  # alias: exp_quadratic
                     `exp_quadratic` = {
                       # y(t) = n0 * exp( ∫rho + b * (t - t0)^2 )
                       if (!has_b) stop("Exponential quadratic model requires parameter b in draws.")
                       n0_mat * exp(rint_mat + b_mat * (dt_mat ^ 2))
                     },
                     {
                       stop("Unknown model: ", model_name)
                     }
  )

  alpha <- (1 - CI) / 2
  ribbon_df <- apply(mean_mat, 2, function(x) {
    qs <- stats::quantile(x, probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE, names = FALSE)
    c(median = qs[2], lower = qs[1], upper = qs[3])
  }) |> t() |> as.data.frame()

  ribbon_df$time <- time_grid
  rownames(ribbon_df) <- NULL

  ribbon_df[, c("time", "median", "lower", "upper")]
}

get_data_for_growth_plot <- function(x, data, time_grid = NULL, CI = 0.95,
                                        time_col = "time", count_col = "count") {

  fit <- x$fit
  breakpoints <- x$breakpoints
  draws <- posterior::as_draws_matrix(fit$draws)
  model_name <- x$best_model
  model_key <- tolower(as.character(model_name))

  time_obs <- data[[time_col]]
  count_obs <- data[[count_col]]
  if (is.null(time_grid)) time_grid <- seq(min(time_obs), max(time_obs), length.out = 200)

  # Indices / presence of parameters in draws
  rho_idx <- grep("^rho\\[", colnames(draws))
  rho_mat <- draws[, rho_idx, drop = FALSE]

  has_t0 <- "t0" %in% colnames(draws)
  has_n0 <- "n0" %in% colnames(draws)
  has_K  <- "K"  %in% colnames(draws)
  has_A  <- "A"  %in% colnames(draws)  # optional, for monomolecular if available
  has_s  <- "sigma" %in% colnames(draws)

  # Fallbacks
  t0_vec <- if (has_t0) draws[, "t0"] else rep(min(time_obs), nrow(draws))
  n0_vec <- if (has_n0) draws[, "n0"] else rep(1, nrow(draws))
  K_vec  <- if (has_K)  draws[, "K"]  else rep(NA_real_, nrow(draws))
  A_vec  <- if (has_A)  draws[, "A"]  else if (has_K) draws[, "K"] else rep(NA_real_, nrow(draws)) # reuse K for A if needed
  s_vec  <- if (has_s)  draws[, "sigma"] else rep(NA_real_, nrow(draws))

  # Pick the appropriate mean function and a small adapter to call it
  # mean_fun_scalar <- switch(
  #   model_key,
  #   "exponential"    = function(t, t0, n0, rho, K, A) mean_piecewise_exponential(t, t0, n0, breakpoints, rho),
  #   "logistic"       = function(t, t0, n0, rho, K, A) mean_piecewise_logistic(t, t0, n0, breakpoints, rho, L = K),
  #   "gompertz"       = function(t, t0, n0, rho, K, A) mean_piecewise_gompertz(t, t0, n0, breakpoints, rho, K = K),
  #   "monomolecular"  = function(t, t0, n0, rho, K, A) mean_piecewise_monomolecular(t, t0, n0, breakpoints, rho, A = A),
  #   "quadraticexp"   = function(t, t0, n0, rho, K, A) mean_piecewise_quadraticexp(t, t0, n0, breakpoints, rho),
  #   # default fallback
  #   function(t, t0, n0, rho, K, A) mean_piecewise_exponential(t, t0, n0, breakpoints, rho)
  # )
  #
  # # Vectorize over time_grid safely (the mean_* functions are scalar in 't')
  # mean_fun_vec <- function(tt, t0, n0, rho, K, A) {
  #   unlist(lapply(tt, function(x) {mean_fun_scalar(x, t0, n0, rho, K, A)}))
  # }

  mean_fun_vec <- switch(
    model_key,
    "exponential"    = function(t, t0, n0, rho, K, A) mean_piecewise_exponential_vec(t, t0, n0, breakpoints, rho),
    "logistic"       = function(t, t0, n0, rho, K, A) mean_piecewise_logistic_vec(t, t0, n0, breakpoints, rho, L = K),
    "gompertz"       = function(t, t0, n0, rho, K, A) mean_piecewise_gompertz_vec(t, t0, n0, breakpoints, rho, K = K),
    "monomolecular"  = function(t, t0, n0, rho, K, A) mean_piecewise_monomolecular_vec(t, t0, n0, breakpoints, rho, A = A),
    "quadraticexp"   = function(t, t0, n0, rho, K, A) mean_piecewise_quadraticexp_vec(t, t0, n0, breakpoints, rho),
    # default fallback
    function(t, t0, n0, rho, K, A) mean_piecewise_exponential_vec(t, t0, n0, breakpoints, rho)
  )

  # Build predictive draws matrix (rows = posterior draws, cols = |time_grid|)
  n_draws <- nrow(draws)
  mean_mat <- matrix(0, nrow = n_draws, ncol = length(time_grid))

  for (i in seq_len(n_draws)) {
    if (has_t0) {
      mu <- mean_fun_vec(time_grid, t0 = t0_vec[i], n0 = 1, rho = rho_mat[i, ], K = K_vec[i], A = A_vec[i])
    } else {
      mu <- mean_fun_vec(time_grid, t0 = time_obs[1], n0 = n0_vec[i], rho = rho_mat[i, ], K = K_vec[i], A = A_vec[i])
    }

    if (has_s && !any(is.na(mu))) {
      mean_mat[i, ] <- rlnorm(length(mu), meanlog = log(pmax(mu, 1e-12)), sdlog = s_vec[i])
    } else {
      mean_mat[i, ] <- rpois(length(mu), lambda = pmax(mu, 0))
    }
  }

  # Credible ribbon
  alpha <- (1 - CI) / 2
  probs <- c(alpha, 0.5, 1 - alpha)

  # Use apply with pre-defined quantile function for efficiency
  ribbon_array <- apply(mean_mat, 2, quantile, probs = probs, na.rm = TRUE, names = FALSE)

  ribbon_df <- data.frame(
    time = time_grid,
    median = ribbon_array[2, ],
    lower = ribbon_array[1, ],
    upper = ribbon_array[3, ]
  )

  # IMPROVEMENT 5: Optional CI smoothing
  smooth_bandwidth = 0.1
  if (length(time_grid) > 10) {
    # Use loess smoothing for CI bounds to reduce noise
    time_range <- max(time_grid) - min(time_grid)
    span <- min(smooth_bandwidth, 1.0)  # Ensure span <= 1

    tryCatch({
      # Smooth the CI bounds but keep median as is (or smooth lightly)
      lower_smooth <- stats::loess(lower ~ time, data = ribbon_df, span = span)$fitted
      upper_smooth <- stats::loess(upper ~ time, data = ribbon_df, span = span)$fitted

      # Optional: light smoothing for median
      median_smooth <- stats::loess(median ~ time, data = ribbon_df, span = span * 0.5)$fitted

      ribbon_df$lower <- lower_smooth
      ribbon_df$upper <- upper_smooth
      ribbon_df$median <- median_smooth

    }, error = function(e) {
      warning("CI smoothing failed, using original values: ", e$message)
    })
  }

  return(ribbon_df[, c("time", "median", "lower", "upper")])
}
