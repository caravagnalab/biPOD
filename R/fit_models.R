
#' Fit growth model to bipod object
#'
#' @param x a bipod object
#' @param growth_type character string specifying the type of growth assumed,
#'  one of "exponential", "logistic".
#' @param model_type character string specifying the model to fit, one of "gauss", "exact".
#'
#'  * `gauss` Use a Gaussian approximation. Works well when counts data are far from 0.
#'
#' @param prior character string specifying the prior to use, "uniform" or "invagamma".
#'
#'  * `uniform` Uniform U(a,b), where at every iteration you sample with width w = (b-a) / g.
#'  * `invgamma` Inverse gamma prior invGamma(a,b)
#'
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param a Minimum value for uniform prior, shape parameter for inverse gamma prior
#' @param b Maximum value for uniform prior, scale parameter for inverse gamma prior
#' @param g Extra parameters for uniform prior. Sample is done from a uniform centered on the previous
#'  sample and with width w = (b - a) / g.
#' @param prior_K Prior mean for the carrying capacity.
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit
#' @export
fit <- function(
    x,
    growth_type = c("exponential", "logistic"),
    model_type = c("gauss", "exact"),
    prior = c("uniform", "invgamma"),
    variational = FALSE,
    factor_size = 1, a = 0, b = 1, g = 1, prior_K = NULL,
    chains = 4, iter = 4000, warmup = 2000, cores = 4){

  growth_type <- match.arg(growth_type)
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  sampling_type <- if(variational) "variational inference" else "MCMC sampling"

  cli::cli_alert_info(paste("Fitting", growth_type, "growth with", model_type, "model and", prior, "prior using", sampling_type, "..."))
  cat("\n")

  if (growth_type == "exponential") {
    res <- fit_exp(x, model_type, prior, factor_size, a, b, g, variational, chains, iter, warmup, cores)
  } else if (growth_type == "logistic") {
    res <- fit_log(x, model_type, prior, variational, factor_size, a, b, g, prior_K, chains, iter, warmup, cores)
  } else {
    stop("'growth_type' should be either 'exponential' or 'logistic'")
  }

  # Add results to bipod object
  x$elbo_data <- res$elbo_data
  x$fits <- res$fits
  x$fit_info <- res$fit_info
  return(x)
}

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
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param a Minimum value for uniform prior, shape parameter for inverse gamma prior
#' @param b Maximum value for uniform prior, scale parameter for inverse gamma prior
#' @param g Extra parameters for uniform prior. Sample is done from a uniform centered on the previous
#'  sample and with width w = (b - a) / g.
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit
#' @export
fit_exp <- function(x, model_type = c("gauss", "exact"), prior = c("uniform", "invgamma"),
                    factor_size = 1, a = 0, b = 1, g = 1, variational = FALSE,
                    chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  # Parameters check
  if (!inherits(x, "bipod")) stop("Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  if (!(factor_size > 0)) stop("'factor_size' must be positive")
  if (prior == "uniform") {
    if (!(a >= 0)) stop("with an uniform prior, 'a' must be >= 0")
    if (!(b > 0)) stop("with an uniform prior, 'b' must be > 0")
    if (!(b > a)) stop("with an uniform prior, 'b' must be greater than 'a'")
    if (!(g >= 1)) stop("with an uniform prior, 'g' must be >= 1")
  } else {
    if (!(a > 0)) stop("with an invgamma prior, 'a' must be positive")
    if (!(b > 0)) stop("with an invgamma prior, 'b' must be positive")
  }

  # Initialize list of fits and ELBO results
  elbo_data <- list()
  fits <- list()
  groups <- unique(x$counts$group)

  for (i in 2:length(groups)) {
    # Obtain info of previous group
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    t0 <- previous$time[nrow(previous)]
    n0 <- previous$count[nrow(previous)]

    # Obtain info of current group
    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    # Prepare the data
    data_model <- list(
      S = nrow(current),
      t0 = t0,
      T = as.array(current$time), # / Ts[length(Ts)],
      n0 = as.integer(n0 / factor_size),
      N = as.array(as.integer(current$count / factor_size)),
      k = 0,
      a = a,
      b = b,
      g = g
    )

    # Get the model
    model_name <- paste("exponential", model_type, prior, sep = "_")
    model <- get(model_name, stanmodels)

    # FIt with either MCMC or Variational
    if (variational) {
      sampling = "variational"
      res <- suppressWarnings(suppressMessages(iterative_variational(model, data_model, iter, warmup)))
      fit_model <- res$fit_model
      elbo_d <- res$elbo_d

    } else {
      sampling = "mcmc"
      fit_model <- rstan::sampling(
        model,
        data = data_model,
        chains = chains, iter = iter, warmup = warmup,
        cores = cores
      )
    }

    if (variational) {
      elbo_data[[paste0("elbo", groups[i])]] <- elbo_d
    }
    fits[[paste0("fit", groups[i])]] <- fit_model
  }

  # Write fit info
  fit_info <- list(
    sampling = sampling,
    growth_type = "exponential",
    model_type = model_type,
    prior = prior,
    factor_size = factor_size,
    a = a,
    b = b,
    g = g
  )

  res <- list(
    elbo_data = elbo_data,
    fits = fits,
    fit_info = fit_info
  )

  return(res)
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
#' @param variational Boolean specifying whether using variational as opposed to mcmc sampling
#' @param factor_size numeric factor by which to divide counts in the bipod object
#' @param a Minimum value for uniform prior, shape parameter for inverse gamma prior
#' @param b Maximum value for uniform prior, scale parameter for inverse gamma prior
#' @param g Extra parameters for uniform prior. Sample is done from a uniform centered on the previous
#'  sample and with width w = (b - a) / g.
#' @param prior_K Prior mean for the carrying capacity.
#' @param chains integer number of chains to run in the Markov Chain Monte Carlo (MCMC) algorithm
#' @param iter integer number of iterations to run in the MCMC algorithm
#' @param warmup integer number of warmup iterations to run in the MCMC algorithm
#' @param cores integer number of cores to use in parallel processing
#'
#' @return the input bipod object with an added 'fit' slot containing the fitted model and an added 'fit_info' slot containing information about the fit
#' @export
fit_log <- function(x, model_type = c("gauss"), prior = c("uniform", "invgamma"), variational = FALSE,
                    factor_size = 1, a = 0, b = 1, g = 1, prior_K = NULL,
                    chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  # Parameters check
  if (!inherits(x, "bipod")) stop("Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  if (!(factor_size > 0)) stop("'factor_size' must be positive")
  if (prior == "uniform") {
    if (!(a >= 0)) stop("with an uniform prior, 'a' must be >= 0")
    if (!(b > 0)) stop("with an uniform prior, 'b' must be > 0")
    if (!(b > a)) stop("with an uniform prior, 'b' must be greater than 'a'")
    if (!(g >= 1)) stop("with an uniform prior, 'g' must be >= 1")
  } else {
    if (!(a > 0)) stop("with an invgamma prior, 'a' must be positive")
    if (!(b > 0)) stop("with an invgamma prior, 'b' must be positive")
  }
  if (is.null(prior_K)) {
    prior_K = max(x$counts$count)
  } else {
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Initialize list of fits and ELBO results
  elbo_data <- list()
  fits <- list()
  groups <- unique(x$counts$group)
  for (i in 2:length(groups)) {
    # Obtain info of previous group
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    t0 <- previous$time[nrow(previous)]
    n0 <- previous$count[nrow(previous)]

    # Obtain info of current group
    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    # Prepare the data
    data_model <- list(
      S = nrow(current),
      t0 = t0,
      T = as.array(current$time), # / Ts[length(Ts)],
      n0 = as.integer(n0 / factor_size),
      N = as.array(as.integer(current$count / factor_size)),
      a = a,
      b = b,
      prior_K = prior_K
    )

    # Get the model
    model_name <- paste("logistic", model_type, prior, sep = "_")
    model <- get(model_name, stanmodels)

    # Fit the model with either MCMC or Variational
    if (variational) {
      sampling = "variational"
      res <- suppressWarnings(suppressMessages(iterative_variational(model, data_model, iter, warmup)))
      fit_model <- res$fit_model
      elbo_d <- res$elbo_d

    } else {
      sampling = "mcmc"
      fit_model <- rstan::sampling(
        model,
        data = data_model,
        chains = chains, iter = iter, warmup = warmup,
        cores = cores
      )
    }

    if (variational) {
      elbo_data[[paste0("elbo", groups[i])]] <- elbo_d
    }
    fits[[paste0("fit", groups[i])]] <- fit_model
  }

  # Write fit info
  fit_info <- list(
    sampling = sampling,
    growth_type = "logistic",
    model_type = model_type,
    prior = prior,
    factor_size = factor_size,
    a = a,
    b = b,
    g = g,
    prior_K = prior_K
  )

  res = list(
    elbo_data = elbo_data,
    fits = fits,
    fit_info = fit_info
  )

  return(res)
}

iterative_variational = function(model, data_model, iter, warmpup) {
  # Iteratively fit with rstan::vb
  # For the first 3/4 of iterations, stop if pareto k value is lower than 0.5
  # For the remaining iterations, it stops if pareto k value is lower than 1
  # If this does not happen, you obtain a fit which is not credible and robus

  N = 20
  for (i in 1:N) {
    out <- capture.output({
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

  elbo_data <- data.frame(iter=iter, ELBO=ELBO, delta_ELBO_mean=delta_ELBO_mean, delta_ELBO_med=delta_ELBO_med) %>%
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
