
#' Fit exponential growth model to bipod object
#'
#' @param x a bipod object
#' @param model_name character string specifying the model to fit, one of "gauss", "exact", or "poisson"
#'
#'  * `exact` Slowest model. Particularly good for low counts data.
#'  * `gauss` Use a gaussian approximation. Works well when counts data are far from 0.
#'  * `poisson` Fastest model. Can be applied both to high and low counts data.
#'
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param fix_rates logical indicating whether to fix rates at values specified in the bipod object
#'
#'  * `1` Birth rates are free. Death rate is fixed
#'  * `0` Both birth and death rate are fixed.
#'  * `-1` Death rates are free. Birth rate is fixed
#'
#' @param prior_par numeric parameter for the prior distribution on rates
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing

#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit

#' @importFrom rstan sampling
#' @export
fit_exp <- function(x, model_name, factor_size = 1, fix_rates = 0, prior_par = 2, chains = 4, iter = 4000, warmup = 2000, cores = 4) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that(model_name %in% c("gauss", "exact", "poisson"), msg = "model_name must be one of gauss, exact, or poisson")
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  assertthat::assert_that(fix_rates %in% c(-1, 0, 1), msg = "fix_rates must be one of -1, 0, or 1")
  assertthat::assert_that(prior_par > 0, msg = "prior_par must be positive")

  model_name <- paste0(model_name, "_exp")

  # prepare data
  Ns <- x$counts$count
  Ts <- x$counts$time

  data_model <- list(
    S = length(Ns) - 1,
    n0 = as.integer(Ns[1] / factor_size),
    N = as.integer(Ns[2:length(Ns)] / factor_size),
    T = Ts[2:length(Ts)], # / Ts[length(Ts)],
    k = fix_rates,
    prior_par = prior_par
  )

  # fit model
  model <- get(model_name, stanmodels)
  fit_model <- rstan::sampling(
    model,
    data = data_model,
    chains = chains, iter = iter, warmup = warmup,
    cores = cores
  )

  # write fit info
  fit_info <- list(
    factor_size = factor_size,
    model_name = model_name,
    fix_rates = fix_rates,
    prior_par = prior_par,
    growth_type = "exponential"
  )

  x$fit <- fit_model
  x$fit_info <- fit_info
  x
}

#' Fit logistic growth model to bipod object
#'
#' @param x a bipod object
#' @param model_name character string specifying the model to fit, currently only "gauss" is supported
#'
#'  * `gauss` Use a gaussian approximation. Works well when counts data are far from 0.
#'
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param fix_rates logical indicating whether to fix rates at values specified in the bipod object
#'
#'  * `1` Growth rate is not fixed at every time step.
#'  * `0` Growth rate is fixed throughout all time steps.
#'
#' @param prior_par numeric parameter for the prior distribution on rates
#' @param prior_K numeric parameter for the prior distribution for the carrying capacity,
#'  if NULL a normal prior centered on the highest count will be used
#'
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit
#'
#' @importFrom rstan sampling
#' @export
fit_log <- function(x, model_name, factor_size = 1, fix_rates = 0, prior_par = 2, prior_K = NULL, chains = 4, iter = 4000, warmup = 2000, cores = 4) {
  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  assertthat::assert_that(model_name %in% c("gauss", "exact", "poisson"), msg = "model_name must be one of gauss, exact, or poisson")
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  assertthat::assert_that(fix_rates %in% c(0, 1), msg = "fix_rates must be either 0 or 1")
  assertthat::assert_that(prior_par > 0, msg = "prior_par must be positive")

  model_name <- paste0(model_name, "_log")

  # prepare data
  Ns <- x$counts$count
  Ts <- x$counts$time

  if (is.null(prior_K)) {
    prior_K = max(x$counts$count)
  }

  data_model <- list(
    S = length(Ns) - 1,
    n0 = as.integer(Ns[1] / factor_size),
    N = as.integer(Ns[2:length(Ns)] / factor_size),
    T = Ts[2:length(Ts)], # / Ts[length(Ts)],
    k = fix_rates,
    prior_par = prior_par,
    max_K = prior_K
  )

  # fit model
  model <- get(model_name, stanmodels)
  fit_model <- rstan::sampling(
    model,
    data = data_model,
    chains = chains, iter = iter, warmup = warmup,
    cores = cores
  )

  # write fit info
  fit_info <- list(
    factor_size = factor_size,
    model_name = model_name,
    fix_rates = fix_rates,
    prior_par = prior_par,
    growth_type = "logistic"
  )

  x$fit <- fit_model
  x$fit_info <- fit_info
  x
}
