#' Fit a Growth Model to a Bipod Object
#'
#' This function fits a specified growth model to a `bipod` object, allowing for exponential
#' or logistic growth fitting. If both models are considered, an automatic model selection
#' process determines the best fit.
#'
#' @param x A `bipod` object containing population count data over time.
#' @param growth_type Character string specifying the growth model to fit. Options are:
#'   - `"exponential"` – Fits an exponential growth model.
#'   - `"logistic"` – Fits a logistic growth model.
#'   - `"both"` – Performs model selection between exponential and logistic growth.
#'   Default is `"exponential"`.
#' @param factor_size Numeric value used to scale population counts in the `bipod` object.
#'   Must be positive and no larger than the minimum count value. Default is `1`.
#' @param infer_t0 Logical value indicating whether to infer the initial time of population origin (`t0`).
#'   If `TRUE`, `t0` is estimated as part of the model fitting. Default is `TRUE`.
#' @param model_selection_algo Character string specifying the model selection algorithm when `growth_type = "both"`.
#'   Options are:
#'   - `"bayes_factor"` – Compares models using Bayes factors.
#'   - `"mixture_model"` – Uses a mixture modeling approach to estimate probabilities.
#'   Default is `"bayes_factor"`.
#' @param variational Logical value indicating whether to use variational inference instead of Markov Chain Monte Carlo (MCMC) sampling.
#'   If `TRUE`, variational inference is applied; otherwise, MCMC is used. Default is `FALSE`.
#' @param chains Integer specifying the number of MCMC chains. Ignored if `variational = TRUE`. Default is `4`.
#' @param cores Integer specifying the number of CPU cores to use for parallel processing. Default is `4`.
#' @param iter Integer specifying the number of MCMC iterations. Ignored if `variational = TRUE`. Default is `5000`.
#'
#' @return Returns the input `bipod` object with additional attributes:
#'   * `fit` – The fitted growth model.
#'   * `fit_info` – Metadata about the fitting process, including:
#'     * Sampling method (MCMC or variational inference)
#'     * Factor size used for scaling
#'     * Selected growth model
#'     * Model selection details (if applicable)
#'
#' @examples
#' # Create a bipod object with your data
#' data = biPOD::sim_stochastic_exponential(100, 1, 0, 10, .25)
#' x = biPOD::init(data, "sample")
#' x <- fit(x, growth_type = "both", model_selection_algo = "bayes_factor")
#' biPOD::plot_fit(x, CI = .8)
#'
#' @export
fit <- function(
    x,
    growth_type = "exponential",
    infer_t0 = TRUE,
    variational = FALSE,
    factor_size = 1,
    model_selection_algo = "bayes_factor",
    chains = 4,
    iter = 5000,
    cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(growth_type %in% c("exponential", "logistic", "both"))) stop("growth_type must be one of 'exponential' and 'logistic'")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  if (factor_size > min(x$counts$count)) stop(paste0("given the input, factor_size must be smaller or equal to ", min(x$counts$count)))
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  if (growth_type == "both") {
    if (!(model_selection_algo %in% c("bayes_factor", "mixture_model"))) stop("model_selection_algo must be one of 'bayes_factor' and 'mixture_model'")
    cli::cli_alert_info(paste("Fitting with model selection."))
    cat("\n")

    if (model_selection_algo == "bayes_factor") {
      func <- fit_with_bayes_factor
    } else {
      func <- fit_with_mixture_model
    }

    res <- func(
      x = x,
      factor_size = factor_size,
      infer_t0 = infer_t0,
      variational = variational,
      chains = chains,
      cores = cores,
      iter = iter
    )
  } else {
    cli::cli_alert_info(paste("Fitting", growth_type, "growth using", sampling_type, "..."))
    cat("\n")

    res <- fit_data(
      x = x,
      infer_t0 = infer_t0,
      growth_type = growth_type,
      factor_size = factor_size,
      variational = variational,
      chains = chains,
      iter = iter,
      cores = cores
    )
  }

  # Add results to bipod object
  x$fit_elbo <- res$elbo_data
  x$fit <- convert_mcmc_fit_to_biPOD(res$fit, variational = variational)

  # Add metadata
  if (sampling_type == "mcmc") {
    x$metadata$status <- diagnose_mcmc_fit(res$fit)
  } else {
    x$metadata$status <- diagnose_variational_fit(res$fit, res$elbo_data)
  }

  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size
  x$metadata$growth_type <- res$fit_info$growth_type
  x$metadata$t0_inferred <- res$fit_info$t0_inferred

  x$metadata$best_growth <- res$fit_info$best_growth
  x$metadata$bayes_factor <- res$fit_info$bayes_factor
  x$metadata$evidence <- res$fit_info$evidence

  x$metadata$odd <- res$fit_info$odds
  x$metadata$model_selection_algo <- res$fit_info$model_selection_algo
  x$metadata$omega_mixture_model <- res$fit_info$omega_mixture_model

  return(x)
}

## Utils fot fitting

prep_data_fit <- function(x, factor_size) {
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
    prior_K = .85 * max(x$counts$count) / factor_size
  )

  return(input_data)
}

fit_data <- function(x, growth_type, factor_size, infer_t0, variational, chains, cores, iter) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size)

  # Get the model
  if (infer_t0) {
    model_name <- growth_type
  } else {
    model_name <- paste0(growth_type, "_no_t0")
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
    t0_inferred = infer_t0
  )

  res <- list(
    elbo_data = elbo_data,
    fit = fit,
    fit_info = fit_info
  )

  return(res)
}

## Bayes factors utils

fit_with_bayes_factor <- function(x, factor_size, infer_t0, variational, chains, cores, iter) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size)

  results <- lapply(c("exponential", "logistic"), function(g) {
    fit_data(
      x = x,
      growth_type = g,
      factor_size = factor_size,
      infer_t0 = infer_t0,
      variational = variational,
      chains = chains,
      cores = cores,
      iter = iter
    )
  })

  fits <- lapply(results, function(r) {
    r$fit
  })
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
    res <- results[[1]]
  } else {
    K <- bayes_factors[2, 1]
    best_growth <- "Logistic"
    evidence <- bayes_factor_evidence(K)
    res <- results[[2]]
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
fit_with_mixture_model <- function(x, factor_size, infer_t0, variational, chains, cores, iter) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size)

  if (infer_t0) {
    model <- get_model(model_name = "exp_log_mixture")
  } else {
    stop("Mixture model not ready without the inference of t0!")
  }

  if (variational) cli::cli_alert_warning("The first step of the inference will be performed using MCMC...")
  tmp <- utils::capture.output(suppressMessages(fit <- model$sample(input_data, chains = chains, iter_warmup = iter, iter_sampling = iter, parallel_chains = cores)))

  # extract omega draws and compute probability of log and of exp
  # omega >= 0.5 suggests Exponential, and vice versa
  n_exponential <- (fit$draws("omega", format = "matrix") >= .5) %>% sum()
  n_logistic <- (fit$draws("omega", format = "matrix") < .5) %>% sum()

  p_exp <- n_exponential / (n_exponential + n_logistic)
  p_log <- 1 - p_exp

  if (p_exp >= p_log) {
    best_growth <- "exponential"
    odds <- p_exp / p_log
  } else {
    best_growth <- "logistic"
    odds <- p_log / p_exp
  }

  res <- fit_data(
    x = x,
    growth_type = best_growth,
    factor_size = factor_size,
    infer_t0 = infer_t0,
    variational = variational,
    chains = chains,
    cores = cores,
    iter = iter
  )

  res$fit_info$best_growth <- best_growth
  res$fit_info$odds <- odds
  res$fit_info$model_selection_algo <- "mixture_model"
  res$fit_info$omega_mixture_model <- fit$draws("omega", format = "matrix")

  return(res)
}
