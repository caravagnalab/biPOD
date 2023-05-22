
prep_data_bp_inference = function(x, factor_size = 1, prior_K = NULL, dt = NULL) {

  # Parameters check
  if (is.null(prior_K)) {
    prior_K = max(x$counts$count) / factor_size
  } else {
    prior_K = prior_K / factor_size
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Get spline prior for changing times
  spl <- stats::smooth.spline(x = x$counts$time, y = x$counts$count, spar = .4)

  x_smooth <- seq(min(x$counts$time), max(x$counts$time), length = 1000)
  y_smooth <- stats::predict(spl, x = x_smooth)[2] %>% unlist()
  dy_smooth <- stats::predict(spl, x = x_smooth, deriv = 1)[2] %>% unlist()

  dy <- diff(y_smooth) / diff(x_smooth)
  dy <- c(dy[1], dy)

  sign_dy <- sign(dy)
  change_idx <- lapply(2:length(dy), function(i) {
    if (sign_dy[i-1] != sign_dy[i]) {return(i)}
  }) %>% unlist() %>% stats::na.omit()

  change_x <- x_smooth[change_idx]

  if (is.null(dt)) { dt <- min(diff(change_x)) / 2}

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = length(as.array(change_x)) + 1,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    changing_times_prior = as.array(change_x),
    dt = dt,
    prior_K = prior_K,
    control = list(max_treedepth = 15, adapt_delta = .8)
  )

  return(input_data)
}

#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param dt half-size of the interval to consider for each proposed breakpoints
#' @param growth_type character string specifying the type of growth assumed,
#'  one of "exponential", "logistic".
#'
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param prior_K Prior mean for the carrying capacity.
#' @param model_selection Boolean, if TRUE the best model between exponential and logistic will be used
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'breakpoints_fit' slot containing the fitted model for the breakpoints
#' @export
breakpoints_inference = function(x, dt = NULL,
                        growth_type = "exponential",
                        variational = FALSE,
                        factor_size = 1, prior_K = NULL,
                        model_selection = FALSE,
                        chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(growth_type %in% c("exponential", "logistic"))) stop("growth_type must be one of 'exponential' and 'logistic'")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if(variational) "variational inference" else "MCMC sampling"

  # Prepare input data
  input_data <- biPOD:::prep_data_bp_inference(x=x, factor_size=factor_size, prior_K=prior_K, dt=dt)

  if (length(input_data$changing_times_prior) == 0) {
    cli::cli_alert_danger("Not possible to infer break points. Zero probable breakpoints were found!")
    return(x)
  }

  # Get the model
  if (model_selection) {
    cli::cli_abort("MODEL SELECTION FOR BP TO DO YET")
  } else {
    if (growth_type == "exponential") {
      model_name <- "exponential_with_changing_points"
    } else if (growth_type == "logistic") {
      model_name <- "logistic_with_changing_points"
    }
  }

  model <- get(model_name, stanmodels)

  # Fit the model with either MCMC or Variational
  if (variational) {
    sampling = "variational"
    res <- suppressWarnings(suppressMessages(biPOD:::iterative_variational(model, input_data, iter, warmup)))
    fit_model <- res$fit_model
    elbo_d <- res$elbo_d

  } else {
    sampling = "mcmc"
    tmp <- capture.output(fit_model <- rstan::sampling(
      model,
      data = input_data,
      chains = chains, iter = iter, warmup = warmup,
      cores = cores
    ))
  }

  elbo_data <- c()
  if (variational) elbo_data <- elbo_d
  fit <- fit_model

  # Add results to bipod object
  x$elbo_data <- elbo_data
  x$breakpoints_fit <- fit

  # Write fit info
  x$metadata$sampling <- sampling
  x$metadata$growth_type <- growth_type
  x$metadata$factor_size <- factor_size
  x$metadata$prior_K <- input_data$prior_K

  # Add median of breakpoints
  breakpoints_names <- names(x$breakpoints_fit)[grepl("changing_times\\[", names(x$breakpoints_fit))]
  breakpoints_sample <- biPOD:::extract_parameters(x$breakpoints_fit, par_list = breakpoints_names)
  x$metadata$breakpoints <- breakpoints_sample %>%
    dplyr::group_by(parameter) %>%
    dplyr::summarise(median = median(value)) %>%
    dplyr::pull(median) %>%
    sort()

  x$counts$group <- biPOD:::bp_to_groups(x$counts, x$metadata$breakpoints)

  cli::cli_alert_success("Breakpoints have been inferred. Inspect the results using the {.field plot_breakpoints_posterior} function.")
  cli::cli_alert_info("Median of the inferred breakpoints have been succesfully stored.")

  return(x)
}
