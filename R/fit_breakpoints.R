#' Propose breakpoints using Differential Evolution
#'
#' Uses the DEoptim algorithm to propose optimal breakpoints for piecewise linear regression
#' on log-transformed counts over time.
#'
#' @param d A data frame with columns `time` and `count`.
#' @param K Integer. Number of segments.
#' @param avg_points_per_window Minimum average number of points per segment.
#' @param n_trials Number of optimization iterations.
#'
#' @return A list with:
#' \item{optimal_breakpoints}{Vector of optimal breakpoint positions.}
#' \item{segment_sizes}{Number of points in each resulting segment.}
propose_breakpoints_DE <- function(d, K, avg_points_per_window = 5, n_trials = 10) {
  # Input validation
  stopifnot(
    is.data.frame(d),
    all(c("time", "count") %in% colnames(d)),
    is.numeric(K) && K > 0,
    is.numeric(avg_points_per_window) && avg_points_per_window > 0,
    nrow(d) >= K * avg_points_per_window
  )

  # Prepare data
  x <- d$time
  y <- log(d$count)
  n <- length(x)

  if (K == 1) {
    cli::cli_alert_info("K=1: No breakpoints needed for single segment")
    return(list(optimal_breakpoints = numeric(0), segment_sizes = nrow(d)))
  }

  assign_segments <- function(time_vec, breakpoints) {
    if (length(breakpoints) == 0) return(rep(1L, length(time_vec)))
    bp_sorted <- sort(breakpoints)
    cuts <- c(-Inf, bp_sorted, Inf)
    cut(time_vec, breaks = cuts, labels = FALSE, right = TRUE, include.lowest = TRUE)
  }

  build_regression_matrix <- function(x_vec, breakpoints) {
    if (length(breakpoints) == 0) return(cbind(intercept = 1, slope = x_vec - min(x_vec)))
    bp_sorted <- sort(c(min(x_vec), breakpoints))
    n_segments <- length(bp_sorted)
    A <- matrix(0, nrow = length(x_vec), ncol = n_segments + 1)
    A[, 1] <- 1
    A[, 2] <- x_vec - bp_sorted[1]
    if (n_segments > 1) {
      for (i in 2:n_segments) A[, i + 1] <- pmax(0, x_vec - bp_sorted[i])
    }
    A
  }

  compute_sse <- function(A, y_vec) {
    tryCatch({
      qr_decomp <- qr(A)
      if (qr_decomp$rank < ncol(A)) return(Inf)
      beta <- qr.coef(qr_decomp, y_vec)
      if (any(is.na(beta))) return(Inf)
      y_hat <- A %*% beta
      sum((y_hat - y_vec)^2)
    }, error = function(e) Inf)
  }

  objective_function <- function(breakpoints, x_vec, y_vec, min_points) {
    segments <- assign_segments(x_vec, breakpoints)
    if (any(tabulate(segments) < min_points)) return(Inf)
    compute_sse(build_regression_matrix(x_vec, breakpoints), y_vec)
  }

  cli::cli_alert_info("Optimizing {K} segments with {K-1} breakpoints")

  n_breakpoints <- K - 1
  if (n_breakpoints == 0) return(numeric(0))
  time_range <- range(x)
  buffer <- diff(time_range) * 0.01
  lower_bounds <- rep(time_range[1] + buffer, n_breakpoints)
  upper_bounds <- rep(time_range[2] - buffer, n_breakpoints)

  control_params <- DEoptim::DEoptim.control(
    NP = max(50, 10 * n_breakpoints), itermax = n_trials,
    reltol = 1e-4, CR = 0.7, strategy = 2, F = 0.8,
    steptol = 50, trace = FALSE
  )

  result <- tryCatch({
    DEoptim::DEoptim(
      fn = function(bp) objective_function(bp, x, y, avg_points_per_window),
      lower = lower_bounds, upper = upper_bounds, control = control_params
    )
  }, error = function(e) {
    cli::cli_alert_danger("Optimization failed: {e$message}")
    return(NULL)
  })

  if (is.null(result)) return(NULL)

  optimal_breakpoints <- sort(result$optim$bestmem)
  final_sse <- result$optim$bestval
  segment_sizes <- tabulate(assign_segments(x, optimal_breakpoints))

  if (any(segment_sizes < avg_points_per_window)) {
    cli::cli_alert_warning("Final breakpoints violate constraints.")
    return(NULL)
  }

  cli::cli_alert_success("Found {length(optimal_breakpoints)} breakpoints with SSE = {round(final_sse, 4)}")
  list(optimal_breakpoints = optimal_breakpoints, segment_sizes = segment_sizes)
}

#' Segment fitting using growth models
#'
#' Iteratively evaluates candidate segmentations and selects the best one
#' using model comparison metrics (BIC or LOO).
#'
#' @inheritParams fit_breakpoints
#' @return A list with candidate evaluations and the best segmentation.
segment_fit <- function(data,
                        with_initiation,
                        noise_model,
                        max_segments = 4,
                        min_segment_size = 3,
                        comparison = c("bic", "loo"),
                        enforce_rho_separation = TRUE,
                        alpha_rho = 0.05,
                        models_to_fit = c("exponential"),
                        chains = 4,
                        iter = 2000,
                        seed = 123,
                        cores = 4) {

  comparison <- match.arg(comparison)
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]

  rhos_distinct <- function(draws, alpha = 0.05) {
    rho_idx <- grep("^rho\\[", colnames(draws))
    rho_mat <- draws[, rho_idx, drop = FALSE]
    G <- ncol(rho_mat)
    if (G <= 1) return(list(distinct = TRUE, merge_pairs = NULL))
    ci_bounds <- apply(rho_mat, 2, stats::quantile, probs = c(alpha/2, 1 - alpha/2))
    merge_pairs <- which(sapply(1:(G - 1), function(g) {
      !(ci_bounds[2, g] < ci_bounds[1, g + 1] || ci_bounds[2, g + 1] < ci_bounds[1, g])
    }))
    list(distinct = length(merge_pairs) == 0, merge_pairs = merge_pairs)
  }

  evaluate_model <- function(breaks) {
    fit_growth_models(
      data = data, breakpoints = breaks, with_initiation = with_initiation,
      chains = chains, iter = iter, seed = seed, cores = cores,
      comparison = comparison, models_to_fit = models_to_fit, noise_model = noise_model
    )
  }

  all_candidates <- list()
  all_evaluations <- list()

  for (k in seq(max_segments, 1)) {
    proposed_bp <- propose_breakpoints_DE(d = data, K = k, avg_points_per_window = min_segment_size, n_trials = 10)$optimal_breakpoints
    if (is.null(proposed_bp)) next

    segments <- sort(c(min(data$time), proposed_bp, max(data$time)))
    seg_lengths <- sapply(seq_along(segments[-1]), function(i) {
      sum(data$time >= segments[i] & data$time <= segments[i + 1])
    })
    if (any(seg_lengths < min_segment_size)) {
      all_candidates[[k]] <- list(valid = FALSE, reason = "min_segment_size")
      next
    }

    res <- evaluate_model(proposed_bp)
    all_evaluations <- append(all_evaluations, list(data.frame(
      num_segments = length(proposed_bp) + 1,
      breakpoints = if (length(proposed_bp) == 0) "none" else paste(round(proposed_bp, 3), collapse = ", "),
      criterion = comparison,
      score = if (comparison == "bic") min(res$model_table$BIC) else res$model_table[1, "looic"]
    )))

    if (enforce_rho_separation) {
      draws <- posterior::as_draws_matrix(res$fits[[1]]$draws())
      if (!rhos_distinct(draws, alpha = alpha_rho)$distinct) {
        all_candidates[[k]] <- list(valid = FALSE, reason = "overlapping rhos")
        next
      }
    }

    metric <- if (comparison == "bic") min(res$model_table$BIC) else res$model_table[1, "looic"]
    all_candidates[[k]] <- list(breakpoints = proposed_bp, result = res, metric = metric, valid = TRUE)
  }

  valid_candidates <- Filter(function(x) isTRUE(x$valid), all_candidates)
  if (length(valid_candidates) == 0) stop("No valid segmentations found.")
  best <- valid_candidates[[which.min(sapply(valid_candidates, function(x) x$metric))]]

  list(fit = parse_stan_fit(best$result$fits$exponential), criterion = comparison,
       best_breakpoints = best$breakpoints, evaluation_table = do.call(rbind, all_evaluations))
}

#' Fit a Bayesian growth model with breakpoints
#'
#' Fits the exponential growth model (with or without initiation) given
#' a set of fixed breakpoints.
#'
#' @inheritParams fit_breakpoints
#' @param breakpoints Proposed breakpoints
#' @return A list with the fitted object and summary statistics.
fit_growth_model_breakpoints <- function(data,
                                         breakpoints,
                                         noise_model,
                                         with_initiation = TRUE,
                                         chains = 4,
                                         iter = 4000,
                                         seed = 123,
                                         cores = 4,
                                         t_prior_sd = 1) {
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]
  G <- length(breakpoints) + 1
  stan_data <- list(S = nrow(data), G = G, N = data$count, T = data$time,
                    t_prior = as.vector(breakpoints), t_prior_sd = t_prior_sd, prior_only = 0)
  model <- if (with_initiation) get_model("exponential_with_init_bp", noise_model) else get_model("exponential_no_init_bp", noise_model)
  message("Fitting breakpoint growth model")
  fit <- suppressMessages(suppressWarnings(model$sample(
    data = stan_data, chains = chains, iter_warmup = iter / 2, iter_sampling = iter / 2,
    seed = seed, parallel_chains = cores, refresh = 0)))
  parsed_fit <- parse_stan_fit(fit)
  list(fit = parsed_fit, summary = fit$summary())
}

#' Fit breakpoints and final growth model
#'
#' This is the main exported function. It first proposes breakpoints via DE,
#' selects the best segmentation using model comparison, and finally refits
#' the model using the best breakpoints with full iterations.
#'
#' @param data Data frame with columns `time` and `count`.
#' @param with_initiation Logical; fit model with initiation parameter.
#' @param max_segments Maximum number of segments to test.
#' @param min_segment_size Minimum size for each segment.
#' @param comparison Criterion for model selection: `"bic"` or `"loo"`.
#' @param enforce_rho_separation Logical; enforce separation of rho parameters.
#' @param alpha_rho Confidence level for rho separation.
#' @param models_to_fit Models to fit, default `"exponential"`.
#' @param chains Number of MCMC chains.
#' @param iter Total iterations.
#' @param seed Random seed.
#' @param cores Parallel chains.
#' @param t_prior_sd Standard deviation for breakpoint prior.
#'
#' @return A list with:
#' \item{evaluation_table}{Model evaluation results.}
#' \item{first_breakpoints}{Breakpoints from segmentation step.}
#' \item{final_breakpoints}{Breakpoints estimated in final fit.}
#' \item{final_fit}{Final fitted model object.}
#' \item{final_summary}{Summary of the final fit.}
#' @export
fit_breakpoints <- function(data,
                            with_initiation,
                            max_segments = 4,
                            min_segment_size = 3,
                            comparison = c("bic", "loo"),
                            enforce_rho_separation = TRUE,
                            alpha_rho = 0.05,
                            models_to_fit = c("exponential"),
                            noise_model = c("lognormal", "poisson"),
                            chains = 4,
                            iter = 4000,
                            seed = 1234,
                            cores = 4,
                            t_prior_sd = 0.5) {

  comparison <- match.arg(comparison)
  noise_model <- match.arg(noise_model)

  if (floor(nrow(data) / min_segment_size) < max_segments) {
    message("Reducing max_segements due to low number of observations")
    max_segments = floor(nrow(data) / min_segment_size)
  }

  seg_res <- segment_fit(
    data = data,
    with_initiation = with_initiation,
    noise_model = noise_model,
    max_segments = max_segments,
    min_segment_size = min_segment_size,
    comparison = comparison,
    enforce_rho_separation = enforce_rho_separation,
    alpha_rho = alpha_rho,
    models_to_fit = models_to_fit,
    chains = chains,
    iter = as.integer(iter / 10),
    seed = seed,
    cores = cores
  )

  eval_table <- seg_res$evaluation_table
  first_bp <- seg_res$best_breakpoints

  final_fit <- fit_growth_model_breakpoints(
    data = data,
    breakpoints = first_bp,
    with_initiation = with_initiation,
    noise_model = noise_model,
    chains = chains,
    iter = iter,
    seed = seed,
    cores = cores,
    t_prior_sd = t_prior_sd
  )

  final_bps <- final_fit$summary %>%
    dplyr::filter(grepl("t_array", .data$variable)) %>%
    dplyr::pull(.data$median) %>%
    sort()

  list(evaluation_table = eval_table, first_breakpoints = first_bp,
       final_breakpoints = final_bps, final_fit = final_fit$fit,
       final_summary = final_fit$summary)
}
