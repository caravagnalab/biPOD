#' Fit a two population growth
#'
#' @param x a bipod object
#'
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'two_pop_fit' slot containing the fitted model and an added 'two_pop_fit_info' slot containing information about the fit
#' @export
fit_two_pop_model <- function(
    x,
    variational = FALSE,
    factor_size = 1,
    chains = 4,
    iter = 4000,
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
  x$two_pop_fit <- res$fit

  # Add metadata
  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size

  return(x)
}
