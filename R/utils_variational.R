iterative_variational <- function(model, data, iter, max_iterations = 100) {
  # Iteratively fit with rstan::vb
  # For the first 3/4 of iterations, stop if pareto k value is lower than 0.5
  # For the remaining iterations, it stops if pareto k value is lower than 1
  # If this does not happen, you obtain a fit which is not credible and robus

  for (i in 1:max_iterations) {
    out <- utils::capture.output({
      warnings <- NULL
      fit_model <- withCallingHandlers(
        {
          model$variational(data = data, output_samples = iter, refresh = iter, adapt_engaged = TRUE)
        },
        warning = function(w) {
          warnings <<- c(warnings, w$message)
        }
      )
      warnings
    })

    pareto_k <- biPOD:::parse_pareto_warning(w = warnings)
    if (i <= max_iterations * 3 / 4) {
      if (pareto_k == 0) break
    } else {
      if (pareto_k <= 1) break
    }
  }

  elbo_d <- suppressWarnings(parse_variational_output(out = out)) %>%
    dplyr::mutate(pareto_k = pareto_k)

  return(list(fit_model = fit_model, elbo_d = elbo_d))
}

parse_variational_output <- function(out) {
  # Parse the output of rstan::vb in order to obtain
  # the values of ELBO and delta_ELBO_mean during the samples
  limits <- c()
  for (i in 1:length(out)) {
    if (length(grep("delta_ELBO_mean", out[i], value = TRUE))) limits <- c(limits, i + 1)
    if (length(grep("Drawing a sample of size", out[i], value = TRUE))) limits <- c(limits, i - 1)
  }

  elbo_lines <- out[c(limits[1]:limits[2])]
  has_converged <- FALSE

  ELBO <- c()
  delta_ELBO_mean <- c()
  delta_ELBO_med <- c()
  iter <- c()
  for (elbo_line in elbo_lines) {
    split_string <- strsplit(elbo_line, " +")[[1]]

    iter <- c(iter, as.numeric(split_string[2]))
    ELBO <- c(ELBO, as.numeric(split_string[3]))
    delta_ELBO_mean <- c(delta_ELBO_mean, as.numeric(split_string[4]))
    delta_ELBO_med <- c(delta_ELBO_med, as.numeric(split_string[5]))

    if (length(grep("MEDIAN ELBO CONVERGED", elbo_line))) has_converged <- TRUE
  }

  elbo_data <- dplyr::tibble(iter = iter, ELBO = ELBO, delta_ELBO_mean = delta_ELBO_mean, delta_ELBO_med = delta_ELBO_med) %>%
    dplyr::mutate(convergence = has_converged)

  elbo_data
}

parse_pareto_warning <- function(w) {
  # Extract the value of the Pareto k diagnostic
  # from the warning or rstan::vb
  if (is.null(w)) {
    return(0)
  }
  w <- gsub(". Resampling", " ", w)
  w <- strsplit(w, " ")[[1]]
  pareto_k <- as.numeric(w[6])
  pareto_k
}
