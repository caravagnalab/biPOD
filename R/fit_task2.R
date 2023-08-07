#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param dt half-size of the interval to consider for each proposed breakpoints
#' @param lambda degree of smoothness of the spline to use
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#' @param desired_changepoints .
#'
#' @return the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
#' @export
breakpoints_inference <- function(
    x,
    dt = NULL,
    variational = FALSE,
    lambda = .5,
    desired_changepoints = NULL,
    factor_size = 1,
    chains = 4,
    iter = 5000,
    cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  # Prepare input data
  input_data <- prep_data_bp_inference(
    x = x,
    factor_size = factor_size,
    dt = dt,
    lambda = lambda,
    n_nodes = desired_changepoints
  )

  if (is.null(input_data)) {
    cli::cli_alert_danger("Not possible to infer the desired number of break points!")
    return(x)
  }

  if (length(input_data$changing_times_prior) == 0) {
    cli::cli_alert_danger("Not possible to infer break points. Zero probable breakpoints were found!")
    return(x)
  }

  # Get the model
  model_name <- "infer_changepoints"
  model <- get_model(model_name = model_name)

  # Fit the model with either MCMC or Variational
  if (variational) {
    sampling <- "variational"
    res <- variational_fit(model = model, data = input_data, iter = iter)
    fit_model <- res$fit_model
    elbo_d <- res$elbo_d
  } else {
    sampling <- "mcmc"
    tmp <- utils::capture.output(
      suppressMessages(
        fit_model <- model$sample(data = input_data, chains = chains, parallel_chains = cores, iter_warmup = iter, iter_sampling = iter, refresh = iter)
      )
    )
  }

  elbo_data <- c()
  if (variational) elbo_data <- elbo_d %>% stats::na.omit()
  fit <- fit_model

  # Add results to bipod object
  x$breakpoints_elbo <- elbo_data
  x$breakpoints_fit <- fit

  # Write fit info
  x$metadata$sampling <- sampling
  x$metadata$factor_size <- factor_size
  x$metadata$prior_K <- input_data$prior_K

  # Add median of breakpoints
  n_changepoints <- length(input_data$changing_times_prior)
  breakpoints_names <- lapply(1:n_changepoints, function(i) {
    paste0("changing_times[", i, "]")
  }) %>% unlist()

  breakpoints_sample <- get_parameters(x$breakpoints_fit, par_list = breakpoints_names)
  x$metadata$breakpoints <- breakpoints_sample %>%
    dplyr::group_by(.data$parameter) %>%
    dplyr::summarise(median = stats::median(.data$value)) %>%
    dplyr::pull(.data$median) %>%
    sort()

  x$counts$group <- bp_to_groups(x$counts, x$metadata$breakpoints)

  cli::cli_alert_success("Breakpoints have been inferred. Inspect the results using the {.field plot_breakpoints_posterior} function.")
  cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")

  x
}


## Utils for breakpoints inference

prep_data_bp_inference <- function(x, factor_size = 1, dt = NULL, lambda = .5, n_nodes = NULL) {
  # Get spline prior for changing times
  xs <- x$counts$time
  ys <- x$counts$count

  if (!(is.null(n_nodes))) {
    tuned_lambda <- tune_lambda(xs, ys, n_nodes = n_nodes)
    if (is.null(tuned_lambda)) {
      return(NULL)
    } # not found
    change_x <- spline_nodes(xs, ys, tuned_lambda)
  } else {
    change_x <- spline_nodes(xs, ys, lambda)
  }

  if (is.null(dt)) {
    if (length(change_x) >= 2) {
      dt <- min(diff(change_x)) / 2
    } else {
      dt <- (max(xs) - min(xs)) / 10
    }
  }

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = length(as.array(change_x)) + 1,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    changing_times_prior = as.array(change_x),
    dt = dt,
    control = list(max_treedepth = 15, adapt_delta = .8)
  )

  return(input_data)
}


tune_lambda <- function(xs, ys, n_nodes) {
  upper <- 1
  lower <- 0

  n_lower <- n_spline_nodes(xs, ys, lower)
  n_upper <- n_spline_nodes(xs, ys, upper)

  if ((n_lower > n_nodes) & (n_upper > n_nodes)) {
    return(NULL)
  }
  if ((n_lower < n_nodes) & (n_upper < n_nodes)) {
    return(NULL)
  }

  while (TRUE) {
    mid <- (lower + upper) / 2 # Calculate the midpoint of the interval
    if (n_spline_nodes(xs, ys, mid) == n_nodes) {
      return(mid)
    } else if (n_spline_nodes(xs, ys, mid) > n_nodes) {
      lower <- mid # Update the lower bound if f(mid) < n
    } else {
      upper <- mid # Update the upper bound otherwise
    }
  }
}

spline_nodes <- function(xs, ys, lambda, n_points = 1000) {
  spl <- stats::smooth.spline(x = xs, y = ys, spar = lambda)

  x_smooth <- seq(min(xs), max(xs), length = n_points)
  y_smooth <- stats::predict(spl, x = x_smooth)[2] %>% unlist()
  dy_smooth <- stats::predict(spl, x = x_smooth, deriv = 1)[2] %>% unlist()

  dy <- diff(y_smooth) / diff(x_smooth)
  dy <- c(dy[1], dy)

  sign_dy <- sign(dy)
  change_idx <- lapply(2:length(dy), function(i) {
    if (sign_dy[i - 1] != sign_dy[i]) {
      return(i)
    }
  }) %>%
    unlist() %>%
    stats::na.omit()

  x_smooth[change_idx]
}

n_spline_nodes <- function(xs, ys, lambda, n_points = 1000) {
  length(spline_nodes(xs, ys, lambda, n_points))
}
