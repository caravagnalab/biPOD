
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
  x$fit_info <- res$fit_info
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
  if (is.null(x$breakpoints)) {
    G <- 1
    breakpoints = array(0, dim=c(0))
  } else {
    G = length(x$breakpoints) + 1
    breakpoints <- x$breakpoints
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
  print(model_name)

  model <- get(model_name, stanmodels)

  print(model)

  # Fit with either MCMC or Variational
  if (variational) {
    sampling = "variational"
    res <- suppressWarnings(suppressMessages(iterative_variational(model, input_data, iter, warmup)))
    fit_model <- res$fit_model
    elbo_d <- res$elbo_d

  } else {
    sampling = "mcmc"
    fit_model <- rstan::sampling(
      model,
      data = input_data,
      chains = chains, iter = iter, warmup = warmup,
      cores = cores
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

bayes_factor_evidence <- function(K) {
  if (K < 1) {
    evidence = "Negative (supports alternative model)"
  } else if (K >= 1 && K < 10^(1/2)) {
    evidence = "Barely worth mentioning"
  } else if (K >= 10^(1/2) && K < 10^(1)) {
    evidence = "Substantial"
  } else if (K >= 10^1 && K < 10^(3/2)) {
    evidence = "Strong"
  } else if (K >= 10^(3/2) && K < 10^(2)) {
    evidence = "Very strong"
  } else {
    evidence = "Decisive"
  }
  evidence
}

fit_with_model_selection <- function(x,
                                     factor_size = 1,
                                     variational = FALSE,
                                     t0_lower_bound = -10,
                                     prior_K = NULL,
                                     chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  input_data <- prep_data_fit(x=x, factor_size=factor_size, prior_K=prior_K, t0_lower_bound=t0_lower_bound)

  res_exp <- fit_data(x=x, growth_type="exponential", factor_size=factor_size,
                  variational=variational, t0_lower_bound=t0_lower_bound, prior_K=prior_K,
                  chains=chains, iter=iter, warmup=warmup, cores=cores)
  res_log <- fit_data(x=x, growth_type="logistic", factor_size=factor_size,
                  variational=variational, t0_lower_bound=t0_lower_bound, prior_K=prior_K,
                  chains=chains, iter=iter, warmup=warmup, cores=cores)

  fits <- list(res_exp$fit, res_log$fit)

  # Compute Marginal Likelihoods for every model
  marginal_likelihoods <- lapply(fits, function(f) {
    b <- bridgesampling::bridge_sampler(f, silent = TRUE)
    unname(b$logml)
  })

  marginal_likelihoods <- marginal_likelihoods %>% unlist()
  # Compute pairwise Bayes factor
  bayes_factors <- outer(marginal_likelihoods, marginal_likelihoods, function(x,y) {
    return(exp(x-y))
  })
  bayes_factors <- dplyr::as_tibble(bayes_factors, .name_repair = "minimal")

  if (bayes_factors[1,2] > 1) {
    K <- bayes_factors[1,2]
    best_growth = "Exponential"
    evidence = bayes_factor_evidence(K)
    res <- res_exp
  } else {
    K <- bayes_factors[2,1]
    best_growth = "Logistic"
    evidence = bayes_factor_evidence(K)
    res <- res_log
  }

  cli::cli_alert_info("Model selection finished!")
  cli::cli_alert_info("Model with {.val {best_growth}} growth deemed better with {.val {evidence}} evidence. (BF = {.val {K}})")

  return(res)
}

iterative_variational = function(model, data_model, iter, warmpup) {
  # Iteratively fit with rstan::vb
  # For the first 3/4 of iterations, stop if pareto k value is lower than 0.5
  # For the remaining iterations, it stops if pareto k value is lower than 1
  # If this does not happen, you obtain a fit which is not credible and robus

  N = 20
  for (i in 1:N) {
    out <- utils::capture.output({
      warnings <- NULL
      fit_model <- withCallingHandlers({
        rstan::vb(
          model, data_model, iter = 100000,  eval_elbo = 50,
          output_samples = iter - warmpup
        )
      }, warning = function(w) { warnings <<- c(warnings, w$message) })
      warnings
    })

    pareto_k <- parse_pareto_warning(w = warnings)
    if (i <= N*3/4) {
      if (pareto_k == 0) break
    } else {
      if (pareto_k <= 1) break
    }
  }

  elbo_d <- parse_variational_output(out = out) %>%
    dplyr::mutate(pareto_k = pareto_k)

  return(list(fit_model = fit_model, elbo_d = elbo_d))
}

parse_variational_output = function(out) {
  # Parse the output of rstan::vb in order to obtain
  # the values of ELBO and delta_ELBO_mean during the samples
  limits <- c()
  for (i in 1:length(out)) {
    if (length(grep("delta_ELBO_mean", out[i], value=TRUE))) limits <- c(limits, i + 1)
    if (length(grep("Drawing a sample of size", out[i], value=TRUE))) limits <- c(limits, i - 2)
  }

  elbo_lines = out[c(limits[1]:limits[2])]
  has_converged = FALSE

  ELBO <- c()
  delta_ELBO_mean <- c()
  delta_ELBO_med <- c()
  iter <- c()
  for (elbo_line in elbo_lines) {
    split_string <- strsplit(elbo_line, " +")[[1]]

    iter <- c(iter, as.numeric(split_string[3]))
    ELBO <- c(ELBO, as.numeric(split_string[4]))
    delta_ELBO_mean <- c(delta_ELBO_mean, as.numeric(split_string[5]))
    delta_ELBO_med <- c(delta_ELBO_med, as.numeric(split_string[6]))
    if (length(grep("MEDIAN ELBO CONVERGED", elbo_line))) has_converged = TRUE
  }

  elbo_data <- dplyr::tibble(iter=iter, ELBO=ELBO, delta_ELBO_mean=delta_ELBO_mean, delta_ELBO_med=delta_ELBO_med) %>%
    dplyr::mutate(convergence = has_converged)

  elbo_data
}

parse_pareto_warning = function(w) {
  # Extract the value of the Pareto k diagnostic
  # from the warning or rstan::vb
  if (is.null(w)) return(0)
  w <- gsub(". Resampling", " ", w)
  w <- strsplit(w, " ")[[1]]
  pareto_k <- as.numeric(w[6])
  pareto_k
}

