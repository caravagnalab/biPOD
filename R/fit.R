#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param growth_type character string specifying the type of growth assumed,
#'  one of "exponential", "logistic".
#'
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param t0_lower_bound lower bound of t0, which is the instant of time in which the population is born
#' @param prior_K Prior mean for the carrying capacity.
#' @param model_selection Boolean, if TRUE the best model between exponential and logistic will be used
#' @param model_selection_algo Algorithm to use for model selection, either 'bayes_factor' or 'mixture_model'
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit
#' @export
fit <- function(
    x,
    growth_type = "exponential",
    variational = FALSE,
    t0_lower_bound = -10,
    factor_size = 1, prior_K = NULL,
    model_selection = FALSE, model_selection_algo = "bayes_factor",
    chains = 4, iter = 4000, cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(growth_type %in% c("exponential", "logistic"))) stop("growth_type must be one of 'exponential' and 'logistic'")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  if (model_selection) {
    if (!(model_selection_algo %in% c("bayes_factor", "mixture_model"))) stop("model_selection_algo must be one of 'bayes_factor' and 'mixture_model'")
    cli::cli_alert_info(paste("Fitting with model selection."))
    cat("\n")

    if (model_selection_algo == "bayes_factor") {
      res <- fit_with_bayes_factor(
        x = x, factor_size = factor_size,
        variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
        chains = chains, iter = iter, cores = cores
      )
    } else {
      res <- fit_with_mixture_model(
        x = x, factor_size = factor_size, variational = variational, t0_lower_bound = t0_lower_bound,
        prior_K = prior_K, chains = chains, iter = iter, cores = cores
      )
    }
  } else {
    cli::cli_alert_info(paste("Fitting", growth_type, "growth using", sampling_type, "..."))
    cat("\n")

    res <- fit_data(
      x = x, growth_type = growth_type, factor_size = factor_size,
      variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
      chains = chains, iter = iter, cores = cores
    )
  }

  # Add results to bipod object
  x$fit_elbo <- res$elbo_data
  x$fit <- res$fit

  # Add metadata
  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size
  x$metadata$growth_type <- res$fit_info$growth_type
  x$metadata$t0_lower_bound <- res$fit_info$t0_lower_bound
  x$metadata$prior_K <- res$fit_info$prior_K

  x$metadata$best_growth <- res$fit_info$best_growth
  x$metadata$bayes_factor <- res$fit_info$bayes_factor
  x$metadata$evidence <- res$fit_info$evidence

  x$metadata$odd <- res$fit_info$odds
  x$metadata$model_selection_algo <- res$fit_info$model_selection_algo
  x$metadata$omega_mixture_model <- res$fit_info$omega_mixture_model

  return(x)
}
