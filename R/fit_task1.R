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
    factor_size = 1,
    prior_K = NULL,
    model_selection = FALSE,
    model_selection_algo = "bayes_factor",
    chains = 4,
    iter = 5000,
    cores = 4) {
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
        x = x,
        factor_size = factor_size,
        variational = variational,
        t0_lower_bound = t0_lower_bound,
        prior_K = prior_K,
        chains = chains,
        iter = iter,
        cores = cores
      )
    } else {
      res <- fit_with_mixture_model(
        x = x,
        factor_size = factor_size,
        variational = variational,
        t0_lower_bound = t0_lower_bound,
        prior_K = prior_K,
        chains = chains,
        iter = iter,
        cores = cores
      )
    }
  } else {
    cli::cli_alert_info(paste("Fitting", growth_type, "growth using", sampling_type, "..."))
    cat("\n")

    res <- fit_data(
      x = x,
      growth_type = growth_type,
      factor_size = factor_size,
      variational = variational,
      t0_lower_bound = t0_lower_bound,
      prior_K = prior_K,
      chains = chains,
      iter = iter,
      cores = cores
    )
  }

  # Add results to bipod object
  x$fit_elbo <- res$elbo_data
  x$fit <- res$fit

  # Add metadata
  if (sampling_type == "mcmc") {
    x$metadata$status <- diagnose_mcmc_fit(res$fit)
  } else {
    x$metadata$status <- diagnose_variational_fit(res$fit, res$elbo_data)
  }

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

## Utils fot fitting

prep_data_fit <- function(x, factor_size, prior_K, t0_lower_bound) {
  # Parameters check
  if (is.null(prior_K)) {
    prior_K <- max(x$counts$count) / factor_size
  } else {
    prior_K <- prior_K / factor_size
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Prepare data
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints <- array(0, dim = c(0))
  } else {
    G <- length(x$metadata$breakpoints) + 1
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

fit_data <- function(
    x,
    growth_type = "exponential",
    factor_size = 1,
    variational = FALSE,
    t0_lower_bound = -10,
    prior_K = NULL,
    chains = 4,
    iter = 4000,
    cores = 4) {
  input_data <- prep_data_fit(
    x = x,
    factor_size = factor_size,
    prior_K = prior_K,
    t0_lower_bound = t0_lower_bound
  )

  # Get the model
  if (t0_lower_bound == x$counts$time[1]) {
    model_name <- paste0(growth_type, "_no_t0")
  } else {
    model_name <- growth_type
  }
  model <- get_model(model_name = model_name)

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

## Bayes factors utils

fit_with_bayes_factor <- function(x,
                                  factor_size = 1,
                                  variational = FALSE,
                                  t0_lower_bound = -10,
                                  prior_K = NULL,
                                  chains = 4, iter = 4000, cores = 4) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size, prior_K = prior_K, t0_lower_bound = t0_lower_bound)

  res_exp <- fit_data(
    x = x, growth_type = "exponential", factor_size = factor_size,
    variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
    chains = chains, iter = iter, cores = cores
  )

  res_log <- fit_data(
    x = x, growth_type = "logistic", factor_size = factor_size,
    variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
    chains = chains, iter = iter, cores = cores
  )

  fits <- list(res_exp$fit, res_log$fit)

  # Compute Marginal Likelihoods for every model

  # Approximated version
  marginal_likelihoods <- lapply(fits, function(f) {
    # b <- bridgesampling::bridge_sampler(f, silent = TRUE)
    f$lp() %>% stats::median()
  })

  # Exact version
  # marginal_likelihoods <- lapply(fits, function(f) {
  #   b <- bridgesampling::bridge_sampler(f, silent = TRUE)
  #   unname(b$logml)
  # })

  marginal_likelihoods <- marginal_likelihoods %>% unlist()
  # Compute pairwise Bayes factor
  bayes_factors <- outer(marginal_likelihoods, marginal_likelihoods, function(x, y) {
    return(exp(x - y))
  })
  bayes_factors <- dplyr::as_tibble(bayes_factors, .name_repair = "minimal")

  if (bayes_factors[1, 2] > 1) {
    K <- bayes_factors[1, 2]
    best_growth <- "Exponential"
    evidence <- bayes_factor_evidence(K)
    res <- res_exp
  } else {
    K <- bayes_factors[2, 1]
    best_growth <- "Logistic"
    evidence <- bayes_factor_evidence(K)
    res <- res_log
  }

  res$fit_info$best_growth <- best_growth
  res$fit_info$bayes_factor <- K
  res$fit_info$evidence <- evidence
  res$fit_info$model_selection_algo <- "bayes_factor"

  cli::cli_alert_info("Model selection finished!")
  cli::cli_alert_info("Model with {.val {best_growth}} growth deemed better with {.val {evidence}} evidence. (BF = {.val {K}})")

  return(res)
}

bayes_factor_evidence <- function(K) {
  if (K < 1) {
    evidence <- "Negative (supports alternative model)"
  } else if (K >= 1 && K < 10^(1 / 2)) {
    evidence <- "Barely worth mentioning"
  } else if (K >= 10^(1 / 2) && K < 10^(1)) {
    evidence <- "Substantial"
  } else if (K >= 10^1 && K < 10^(3 / 2)) {
    evidence <- "Strong"
  } else if (K >= 10^(3 / 2) && K < 10^(2)) {
    evidence <- "Very strong"
  } else {
    evidence <- "Decisive"
  }
  evidence
}

## Mixture model utils

fit_with_mixture_model <- function(
    x,
    factor_size = 1,
    variational = FALSE,
    t0_lower_bound = -10,
    prior_K = NULL,
    chains = 4,
    iter = 4000,
    cores = 4) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size, prior_K = prior_K, t0_lower_bound = t0_lower_bound)

  model <- get_model(model_name = "exp_log_mixture")

  tmp <- utils::capture.output(suppressMessages(fit <- model$sample(input_data, chains = chains, iter_warmup = iter, iter_sampling = iter, parallel_chains = cores)))

  # extract omega draws and compute probability of log and of exp
  # omega >= 0.5 suggests Exponential, and vice versa
  n_exponential <- (fit$draws("omega", format = "matrix") >= .5) %>% sum()
  n_logistic <- (fit$draws("omega", format = "matrix") < .5) %>% sum()

  p_exp <- n_exponential / (n_exponential + n_logistic)
  p_log <- 1 - p_exp

  if (p_exp >= p_log) {
    best_growth <- "Exponential"
    odds <- p_exp / p_log
    res <- fit_data(
      x = x,
      growth_type = "exponential",
      factor_size = factor_size,
      variational = variational,
      t0_lower_bound = t0_lower_bound,
      prior_K = prior_K,
      chains = chains,
      iter = iter,
      cores = cores
    )
  } else {
    best_growth <- "Logistic"
    odds <- p_log / p_exp
    res <- fit_data(
      x = x,
      growth_type = "logistic",
      factor_size = factor_size,
      variational = variational,
      t0_lower_bound = t0_lower_bound,
      prior_K = prior_K,
      chains = chains,
      iter = iter,
      cores = cores
    )
  }

  res$fit_info$best_growth <- best_growth
  res$fit_info$odds <- odds
  res$fit_info$model_selection_algo <- "mixture_model"
  res$fit_info$omega_mixture_model <- fit$draws("omega", format = "matrix")

  return(res)
}
