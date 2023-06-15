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
breakpoints_inference <- function(x,
                                  dt = NULL,
                                  variational = FALSE,
                                  lambda = .5,
                                  desired_changepoints = NULL,
                                  factor_size = 1,
                                  chains = 4,
                                  iter = 4000,
                                  cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  # Prepare input data
  input_data <- prep_data_bp_inference(x = x, factor_size = factor_size, dt = dt, lambda = lambda, n_nodes = desired_changepoints)

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
    tmp <- utils::capture.output(suppressMessages(fit_model <- model$sample(data = input_data, chains = chains, parallel_chains = cores, iter_warmup = iter, iter_sampling = iter, refresh = iter)))
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
