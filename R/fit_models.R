
#' Fit a Bayesian model for the evolution of a population.
#'
#' @param x A biPOD object of class `bipod`.
#' @param factor_size Integer. Counts will be divided by this factor size.
#' @param model_name String. Indicates the name of the model desired.
#'
#'  * `exact` Slowest model. Particularly good for low counts data.
#'  * `gauss` Use a gaussian approximation. Works well when counts data are far from 0.
#'  * `poisson` Fastest model. Can be applied both to high and low counts data.
#'
#' @param fix_rates An integer value. Indicates which rates must be fixed.
#'
#'  * `1` Birth rates are free. Death rate is fixed
#'  * `0` Both birth and death rate are fixed.
#'  * `-1` Death rates are free. Birth rate is fixed
#'
#' @param prior_par Parameter for inverse gamma prior. InvGamma(p,p).
#' @param chains Integer. Number of MCMC chains
#' @param iter Integer. Number of MCMC steps
#' @param warmup Integer. Number of MCMC warmup steps, must be smaller than 'chains'.
#' @param cores Integer. Number of cores to use.
#' @returns A plot. Represents the evolution of the population over time.
#' @importFrom rstan sampling
#' @export
fit_exp <- function(x, model_name, factor_size = 1, fix_rates = 0, prior_par = 2, chains = 4, iter = 4000, warmup = 2000, cores = 4) {
  stopifnot(inherits(x, "bipod"))
  stopifnot(model_name %in% c("gauss", "exact", "poisson"))

  # prepare data
  Ns <- x$counts$count
  Ts <- x$counts$time

  data_model <- list(
    S = length(Ns) - 1,
    n0 = as.integer(Ns[1] / factor_size),
    N = as.integer(Ns[2:length(Ns)] / factor_size),
    T = Ts[2:length(Ts)], #/ Ts[length(Ts)],
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
    prior_par = prior_par
  )

  x$fit <- fit_model
  x$fit_info <- fit_info
  x
}
