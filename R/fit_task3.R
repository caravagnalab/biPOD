#' Fit a Two-Population Growth Model to a bipod Object
#'
#' This function fits a growth model that accounts for two populations within a bipod object. The fitting can be done using either Variational Inference or Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @param x A `bipod` object.
#' @param variational A logical value indicating whether to use Variational Inference instead of MCMC sampling.
#'  If `TRUE`, the model will be fitted using Variational Inference; otherwise, MCMC sampling will be used. (default is FALSE)
#' @param factor_size A numeric value representing the factor by which to divide the counts in the bipod object.
#'  This value must be positive and appropriate for the data scale. (default is 1)
#' @param chains An integer specifying the number of chains to run in the MCMC algorithm.
#'  This parameter is ignored if `variational = TRUE`. (default is 4)
#' @param iter An integer specifying the number of iterations to run in the MCMC algorithm.
#'  This parameter is ignored if `variational = TRUE`. (default is 5000)
#' @param cores An integer specifying the number of cores to use for parallel processing during model fitting. (default is 4)
#'
#' @return The input `bipod` object with added slots:
#' - `'two_pop_fit'`: Contains the fitted two-population growth model.
#' - `'two_pop_fit_info'`: Contains information about the fitting process, including metadata such as sampling type and factor size.
#'
#' @export
fit_two_pop_model <- function(
    x,
    variational = FALSE,
    factor_size = 1,
    chains = 4,
    iter = 5000,
    cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  cli::cli_alert_info(paste("Fitting two population model using", sampling_type, "..."))
  cat("\n")

  res <- fit_two_pop_data(
    x = x,
    factor_size = factor_size,
    variational = variational,
    chains = chains,
    iter = iter,
    cores = cores
  )

  # Add results to bipod object
  x$two_pop_fit_elbo <- res$elbo_data
  x$two_pop_fit <- convert_mcmc_fit_to_biPOD(res$fit)

  if (sampling_type == "mcmc") {
    x$metadata$status <- diagnose_mcmc_fit(res$fit)
  } else {
    x$metadata$status <- diagnose_variational_fit(res$fit, res$elbo_data)
  }

  # Add metadata
  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size

  return(x)
}

# Utils fitting taks 3

prep_data_two_pop_fit <- function(x, factor_size) {
  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time)
  )
  return(input_data)
}

fit_two_pop_data <- function(x, factor_size, variational, chains, iter, cores) {
  input_data <- prep_data_two_pop_fit(x = x, factor_size = factor_size)
  model <- get_model(model_name = "two_pop")

  # Fit with either MCMC or Variational
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
  if (variational) elbo_data <- elbo_d
  fit <- fit_model

  # Write fit info
  fit_info <- list(
    sampling = sampling,
    factor_size = factor_size
  )

  res <- list(
    elbo_data = elbo_data,
    fit = fit,
    fit_info = fit_info
  )

  return(res)
}
