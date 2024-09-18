group_contiguous <- function(x) {
  rle_x <- rle(x)
  return(rep(seq_along(rle_x$lengths) - 1, times = rle_x$lengths))
}

exp_growth <- function(t, t0, t_array, rho_array, n0) {
  if (length(t_array) == 0) {
    return(n0 * exp(rho_array[1] * (t - t0)))
  }

  if (t <= t_array[1]) {
    return(n0 * exp(rho_array[1] * (t - t0)))
  }

  res <- n0 * exp(rho_array[1] * (t_array[1] - t0))

  if (length(t_array) >= 2) {
    for (i in 2:length(t_array)) {
      if (t <= t_array[i]) {
        return(res * exp(rho_array[i] * (t - t_array[i - 1])))
      } else {
        res <- res * exp(rho_array[i] * (t_array[i] - t_array[i - 1]))
      }
    }
  }

  res <- res * exp(rho_array[length(rho_array)] * (t - t_array[length(t_array)]))
  return(res)
}

log_growth <- function(t, n0, rho, K) {
  num <- K * n0
  den <- n0 + (K - n0) * exp(-rho * t)
  return(num / den)
}

log_growth_multiple <- function(t, t0, t_array, rho_array, K, n0) {
  current_n0 <- n0
  if (length(t_array) == 0) {
    return(log_growth(t - t0, current_n0, rho_array[1], K))
  }

  if (t <= t_array[1]) {
    return(log_growth(t - t0, current_n0, rho_array[1], K))
  }

  dt <- t_array[1] - t0
  current_n0 <- log_growth(dt, current_n0, rho_array[1], K)
  if (length(t_array) >= 2) {
    for (i in 2:length(t_array)) {
      if (t <= t_array[i]) {
        dt <- t - t_array[i - 1]
        return(log_growth(dt, current_n0, rho_array[i], K))
      } else {
        dt <- t_array[i] - t_array[i - 1]
        current_n0 <- log_growth(dt, current_n0, rho_array[i], K)
      }
    }
  }

  dt <- t - t_array[length(t_array)]
  return(log_growth(dt, current_n0, rho_array[length(rho_array)], K))
}


two_pops_evo <- function(t, ns, t0_r, rho_s, rho_r) {
  s_pop <- ns * exp(rho_s * (t))
  if (t > t0_r) {
    r_pop <- 1 * exp(rho_r * (t - t0_r))
  } else {
    r_pop <- 0
  }

  list(
    r_pop = r_pop,
    s_pop = s_pop
  )
}


diagnose_mcmc_fit <- function(fit) {
  all_pars <- colnames(fit$draws(format = "matrix"))
  n_chains <- ncol(fit$draws("lp__"))

  status <- "PASS"
  for (par in all_pars) {
    rhat <- posterior::rhat(fit$draws(par))
    if (is.na(rhat)) {
      status <- "PASS"
    } else if (rhat > 1.01) {
      status <- "FAIL"
      break()
    }
  }

  status
}

diagnose_variational_fit <- function(fit, elbo_data) {
  elbo_converged <- all(elbo_data$convergence)

  if (elbo_converged) {
    status <- "PASS"
  } else {
    status <- "FAIL"
  }

  status
}


variational_fit <- function(model, data, iter, max_iterations = 10) {
  fitted <- FALSE
  attempts <- 1

  while (attempts <= max_iterations) {
    tryCatch(
      expr = {
        out <- utils::capture.output(
          fit_model <- model$variational(data = data, output_samples = iter, iter = iter * 10, eta = .1, adapt_engaged = T, adapt_iter = 200, algorithm = "meanfield")
        )
        fitted <- TRUE
      },
      error = function(err) {
        fitted <- FALSE
      },
      warning = function(warn) {
        fitted <- FALSE
      }
    )

    if (fitted) {
      elbo_d <- suppressWarnings(parse_variational_output(out = out))
      return(list(fit_model = fit_model, elbo_d = elbo_d))
    } else {
      attempts <- attempts + 1
    }
  }

  return(NULL)
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


convert_mcmc_fit_to_biPOD <- function(mcmc_fit) {
  # convert mcmc fit into dplyr::tibble
  parameters <- mcmc_fit$metadata()$model_params
  draws <- mcmc_fit$draws() %>% dplyr::as_tibble() %>% dplyr::mutate_all(as.numeric)

  rhats <- lapply(parameters, function(p) { compute_rhat(draws, p) })
  names(rhats) <- parameters
  rhats

  list(parameters = parameters, draws = draws, rhat = rhats)
}

compute_rhat <- function(draws, par_name) {
  X <- draws[grepl(par_name, colnames(draws), fixed = T)] %>% as.matrix()

  N <- nrow(X)
  M <- ncol(X)

  sm_squared <- function(x) {
    x.mean = mean(x)
    s <- sum((x - x.mean)**2)
    return(s / (length(x) - 1))
  }

  sms <- lapply(1:ncol(X), function(m) {
    sm_squared(X[,m])
  }) %>% unlist()

  W <- sum(sms) / M
  B <- sum((colMeans(X) - mean(X))**2) * N / (M - 1)

  V <- W * (N - 1) / N + B * (M + 1) / (N * M)

  sqrt(V / W)
}
