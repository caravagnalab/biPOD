
#' Fit exponential growth model to bipod object
#'
#' @param x a bipod object
#' @param model_type character string specifying the model to fit, one of "gauss", "exact".
#'
#'  * `exact` Slowest model. Particularly good for low counts data.
#'  * `gauss` Use a Gaussian approximation. Works well when counts data are far from 0.
#'
#' @param prior character string specifying the prior to use, "uniform" or "invagamma".
#'
#'  * `uniform` Uniform U(a,b), where at every iteration you sample with width w = (b-a) / g.
#'  * `invgamma` Inverse gamma prior invGamma(a,b)
#'
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param fix_rates logical indicating whether to fix rates at values specified in the bipod object
#'
#'  * `1` The birth rate is fixed whereas the death rates are variable.
#'  * `0` Both birth and death rate are fixed.
#'
#' @param a Minimum value for uniform prior, shape parameter for inverse gamma prior
#' @param b Maximum value for uniform prior, scale parameter for inverse gamma prior
#' @param g Extra parameters for uniform prior. Sample is done from a uniform centered on the previous
#'  sample and with width w = (b - a) / g.
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing

#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit

#' @importFrom rstan sampling
#' @export
fit_exp <- function(x, model_type = c("gauss", "exact"), prior = c("uniform", "invgamma"),
                    factor_size = 1, fix_rates = 0, a = 0, b = 1, g = 1,
                    chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  assertthat::assert_that(fix_rates %in% c(0,1), msg = "fix_rates must be either '0' or '1'")
  if (prior == "uniform") {
    assertthat::assert_that(a >= 0, msg = "'a' must be greater or equal to 0")
    assertthat::assert_that(b >= 0, msg = "'b' must be greater or equal to 0")
    assertthat::assert_that(b > a, msg = "'b' must be greater than 'a'")
    assertthat::assert_that(g >= 1, msg = "'g' must be greater or equal to 0" )
  } else {
    assertthat::assert_that(a > 0, msg = "'a' must be greater than 0")
    assertthat::assert_that(b > 0, msg = "'b' must be greater than 0")
  }

  # prepare data
  Ns <- x$counts$count
  Ts <- x$counts$time

  data_model <- list(
    S = length(Ns) - 1,
    n0 = as.integer(Ns[1] / factor_size),
    N = as.integer(Ns[2:length(Ns)] / factor_size),
    T = Ts[2:length(Ts)], # / Ts[length(Ts)],
    k = fix_rates,
    a = a,
    b = b,
    g = g
  )

  # fit model
  model_name <- paste("exponential", model_type, prior, sep = "_")
  model <- get(model_name, stanmodels)
  fit_model <- rstan::sampling(
    model,
    data = data_model,
    chains = chains, iter = iter, warmup = warmup,
    cores = cores
  )

  # write fit info
  fit_info <- list(
    growth_type = "exponential",
    model_type = model_type,
    prior = prior,
    factor_size = factor_size,
    fix_rates = fix_rates,
    a = a,
    b = b,
    g = g
  )

  x$fit <- fit_model
  x$fit_info <- fit_info
  x
}

#' Fit logistic growth model to bipod object
#'
#' @param x a bipod object
#' @param model_type character string specifying the model to fit, one of "gauss".
#'
#'  * `gauss` Use a Gaussian approximation. Works well when counts data are far from 0.
#'
#' @param prior character string specifying the prior to use, "uniform" or "invagamma".
#'
#'  * `uniform` Uniform U(a,b), where at every iteration you sample with width w = (b-a) / g.
#'  * `invgamma` Inverse gamma prior invGamma(a,b)
#'
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param fix_rates logical indicating whether to fix rates at values specified in the bipod object
#'
#'  * `1` The birth rate is fixed whereas the death rates are variable.
#'  * `0` Both birth and death rate are fixed.
#'
#' @param a Minimum value for uniform prior, shape parameter for inverse gamma prior
#' @param b Maximum value for uniform prior, scale parameter for inverse gamma prior
#' @param g Extra parameters for uniform prior. Sample is done from a uniform centered on the previous
#'  sample and with width w = (b - a) / g.
#' @param prior_K Prior mean for the carrying capacity.
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing

#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit

#' @importFrom rstan sampling
#' @export
fit_log <- function(x, model_type = c("gauss"), prior = c("uniform", "invgamma"),
                    factor_size = 1, fix_rates = 1, a = 0, b = 1, g = 1, prior_K = NULL,
                    chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  assertthat::assert_that(fix_rates %in% c(0,1), msg = "fix_rates must be either '0' or '1'")
  if (prior == "uniform") {
    assertthat::assert_that(a >= 0, msg = "'a' must be greater or equal to 0")
    assertthat::assert_that(b >= 0, msg = "'b' must be greater or equal to 0")
    assertthat::assert_that(b > a, msg = "'b' must be greater than 'a'")
    assertthat::assert_that(g >= 1, msg = "'g' must be greater or equal to 0" )
  } else {
    assertthat::assert_that(a > 0, msg = "'a' must be greater than 0")
    assertthat::assert_that(b > 0, msg = "'b' must be greater than 0")
  }

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
    prior_K = prior_K
  )

  # fit model
  model_name <- paste("logistic", model_type, prior, sep = "_")
  model <- get(model_name, stanmodels)
  fit_model <- rstan::sampling(
    model,
    data = data_model,
    chains = chains, iter = iter, warmup = warmup,
    cores = cores
  )

  # write fit info
  # write fit info
  fit_info <- list(
    growth_type = "logistic",
    model_type = model_type,
    prior = prior,
    factor_size = factor_size,
    fix_rates = fix_rates,
    a = a,
    b = b,
    g = g,
    prior_K = prior_K
  )

  x$fit <- fit_model
  x$fit_info <- fit_info
  x
}
