
#' Fit a Bayesian model for the evolution of a population.
#'
#' @param x A biPOD object of class `bipod`.
#' @param factor_size Integer. Counts will be divided by this factor size.
#' @param gaussian_approx Boolean. If True, a Gaussian approximation will be applied
#' @param rates_mode Boolean. If True, a Gaussian approximation will be applied
#' to the distributions used for the inference.
#' @param chains Integer. Number of MCMC chains
#' @param iter Integer. Number of MCMC steps
#' @param warmup Integer. Number of MCMC warmup steps, must be smaller than 'chains'.
#' @param cores Integer. Number of cores to use.
#' @returns A plot. Represents the evolution of the population over time.
#' @importFrom rstan stan_model sampling
#' @export
fit_exp <- function(x, factor_size = 1, gaussian_approx = F, rates_mode = "fixed_all", chains = 4, iter = 4000, warmup = 2000, cores = 4) {
  stopifnot(inherits(x, "bipod"))

  # prepare data
  Ns <- x$counts$count
  Ts <- x$counts$time

  data_model <- list(
    S = length(Ns) - 1,
    n0 = as.integer(Ns[1] / factor_size),
    N = as.integer(Ns[2:length(Ns)] / factor_size),
    T = Ts[2:length(Ts)] / Ts[length(Ts)]
  )

  # choose model name
  m0 <- if (gaussian_approx) "gauss" else "exact"
  model_name <- paste(m0, rates_mode, sep = "_")

  # fit model
  model <- get(model_name, stanmodels)
  fit_model <- sampling(
    model,
    data = data_model,
    chains = chains, iter = iter, warmup = warmup,
    cores = cores
  )

  x$fit <- fit_model
  x
}
