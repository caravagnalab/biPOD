
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
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
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
    model_selection = FALSE,
    chains = 4, iter = 4000, warmup = 2000, cores = 4){

  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(growth_type %in% c("exponential", "logistic"))) stop("growth_type must be one of 'exponential' and 'logistic'")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if(variational) "variational inference" else "MCMC sampling"

  if (model_selection) {
    cli::cli_alert_info(paste("Fitting with model selection."))
    cat("\n")

    res <- fit_with_model_selection(x=x, factor_size = factor_size,
                                    variational = variational, t0_lower_bound = t0_lower_bound, prior_K=prior_K,
                                    chains=chains, iter = iter, warmup = warmup, cores = cores)
  } else {
    cli::cli_alert_info(paste("Fitting", growth_type, "growth using", sampling_type, "..."))
    cat("\n")

    res <- fit_data(x=x, growth_type=growth_type, factor_size=factor_size,
                    variational=variational, t0_lower_bound=t0_lower_bound, prior_K=prior_K,
                    chains=chains, iter=iter, warmup=warmup, cores=cores)
  }

  # Add results to bipod object
  x$elbo_data <- res$elbo_data
  x$fit <- res$fit

  # Add metadata
  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size
  x$metadata$growth_type <- res$fit_info$growth_type
  x$metadata$t0_lower_bound <- res$fit_info$t0_lower_bound
  x$metadata$prior_K <- res$fit_info$prior_K

  return(x)
}

prep_data_fit = function(x, factor_size, prior_K, t0_lower_bound) {
  # Parameters check
  if (is.null(prior_K)) {
    prior_K = max(x$counts$count) / factor_size
  } else {
    prior_K = prior_K / factor_size
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Prepare data
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints = array(0, dim=c(0))
  } else {
    G = length(x$metadata$breakpoints) + 1
    breakpoints <- x$metadata$breakpoints
  }

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = G,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    t_array = as.array(breakpoints),
    t0_lower_bound = t0_lower_bound,
    prior_K = prior_K
  )

  return(input_data)
}

# Fit exponential growth model to bipod object
fit_data <- function(x,
                     growth_type = "exponential",
                     factor_size = 1,
                     variational = FALSE,
                     t0_lower_bound = -10,
                     prior_K = NULL,
                     chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  input_data <- prep_data_fit(x=x, factor_size=factor_size, prior_K=prior_K, t0_lower_bound=t0_lower_bound)

  # Get the model
  if (t0_lower_bound == x$counts$time[1]) {
    model_name <- paste0(growth_type, '_start_at_1')
  } else {
    model_name <- growth_type
  }
  model <- get(model_name, stanmodels)

  # Fit with either MCMC or Variational
  if (variational) {
    sampling = "variational"
    res <- suppressWarnings(suppressMessages(iterative_variational(model, input_data, iter, warmup)))
    fit_model <- res$fit_model
    elbo_d <- res$elbo_d

  } else {
    sampling = "mcmc"
    out <- utils::capture.output(fit_model <- rstan::sampling(
      model,
      data = input_data,
      chains = chains, iter = iter, warmup = warmup,
      cores = cores
    ))
  }

  elbo_data <- c()
  if (variational) elbo_data <- elbo_d
  fit <- fit_model

  # Write fit info
  fit_info <- list(
    sampling = sampling,
    growth_type = growth_type,
    factor_size = factor_size,
    t0_lower_bound = t0_lower_bound,
    prior_K = input_data$prior_K
  )

  res <- list(
    elbo_data = elbo_data,
    fit = fit,
    fit_info = fit_info
  )

  return(res)
}
