
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

  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  if (prior == "uniform") {
    assertthat::assert_that(a >= 0, msg = "'a' must be greater or equal to 0")
    assertthat::assert_that(b >= 0, msg = "'b' must be greater or equal to 0")
    assertthat::assert_that(b > a, msg = "'b' must be greater than 'a'")
    assertthat::assert_that(g >= 1, msg = "'g' must be greater or equal to 0" )
  } else {
    assertthat::assert_that(a > 0, msg = "'a' must be greater than 0")
    assertthat::assert_that(b > 0, msg = "'b' must be greater than 0")
  }

  elbo_data <- list()
  fits <- list()
  groups <- unique(x$counts$group)
  for (i in 2:length(groups)) {
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    t0 <- previous$time[nrow(previous)]
    n0 <- previous$count[nrow(previous)]

    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    # Prepare the data
    data_model <- list(
      S = nrow(current),
      t0 = t0,
      T = current$time, # / Ts[length(Ts)],
      n0 = as.integer(n0 / factor_size),
      N = c(as.integer(current$count / factor_size)),
      k = 0,
      a = a,
      b = b,
      g = g
    )

    # fit model
    model_name <- paste("exponential", model_type, prior, sep = "_")
    model <- get(model_name, stanmodels)

    # MCMC or Variational
    if (variational) {
      sampling = "variational"
      out <- capture.output(
        fit_model <- rstan::vb(
          model, data_model, iter = 100000,  eval_elbo = 50,
          output_samples = iter - warmup
        )
      )
      elbo_d <- parse_variational_output(out = out)

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

  # write fit info
  fit_info <- list(
    growth_type = "exponential",
    model_type = model_type,
    prior = prior,
    factor_size = factor_size,
    a = a,
    b = b,
    g = g
  )

  x$elbo_data <- elbo_data
  x$fits <- fits
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
fit_log <- function(x, model_type = c("gauss", "exact"), prior = c("uniform", "invgamma"), variational = FALSE,
                    factor_size = 1, a = 0, b = 1, g = 1, prior_K = NULL,
                    chains = 4, iter = 4000, warmup = 2000, cores = 4) {

  assertthat::assert_that(inherits(x, "bipod"), msg = "Input must be a bipod object")
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  assertthat::assert_that(factor_size > 0, msg = "factor_size must be positive")
  if (prior == "uniform") {
    assertthat::assert_that(a >= 0, msg = "'a' must be greater or equal to 0")
    assertthat::assert_that(b >= 0, msg = "'b' must be greater or equal to 0")
    assertthat::assert_that(b > a, msg = "'b' must be greater than 'a'")
    assertthat::assert_that(g >= 1, msg = "'g' must be greater or equal to 0" )
  } else {
    assertthat::assert_that(a > 0, msg = "'a' must be greater than 0")
    assertthat::assert_that(b > 0, msg = "'b' must be greater than 0")
  }

  if (is.null(prior_K)) {
    prior_K = max(x$counts$count)
  }

  elbo_data <- list()
  fits <- list()
  groups <- unique(x$counts$group)
  for (i in 2:length(groups)) {
    previous <- x$counts %>%
      dplyr::filter(.data$group == groups[i-1])

    t0 <- previous$time[nrow(previous)]
    n0 <- previous$count[nrow(previous)]

    current <- x$counts %>%
      dplyr::filter(.data$group == groups[i])

    data_model <- list(
      S = nrow(current),
      t0 = t0,
      T = current$time,
      n0 = as.integer(n0 / factor_size),
      N = as.integer(current$count / factor_size),
      a = a,
      b = b,
      prior_K = prior_K
    )

    # fit model
    model_name <- paste("logistic", model_type, prior, sep = "_")
    model <- get(model_name, stanmodels)

    # MCMC or Variational
    if (variational) {
      sampling = "variational"
      out <- capture.output(
        fit_model <- rstan::vb(
          model, data_model, iter = 100000,
          importance_resampling = TRUE, output_samples = iter - warmup
        )
      )
      elbo_d <- parse_variational_output(out = out)

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

  # write fit info
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

  x$elbo_data <- elbo_data
  x$fits <- fits
  x$fit_info <- fit_info
  x
}

parse_variational_output = function(out) {
  limits <- c()
  for (i in 1:length(out)) {
    if (length(grep("delta_ELBO_mean", out[i], value=TRUE))) limits <- c(limits, i + 1)
    if (length(grep("Drawing a sample of size", out[i], value=TRUE))) limits <- c(limits, i - 2)
  }

  elbo_lines = out[c(limits[1]:limits[2])]

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
  }

  elbo_data <- data.frame(iter=iter, ELBO=ELBO, delta_ELBO_mean=delta_ELBO_mean, delta_ELBO_med=delta_ELBO_med)
  elbo_data
}
