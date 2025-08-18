#' Fit growth models and select the best one
#'
#' Fits multiple candidate growth models (exponential, logistic, Gompertz) to time-series
#' count data with optional breakpoints, using either MCMC sampling or variational inference (VI).
#' The best model is selected based on the specified criterion (LOO, BIC, or ELBO).
#'
#' @param data Data frame with columns `time` and `count`.
#' @param breakpoints Numeric vector of breakpoints.
#' @param with_initiation Logical; whether to include an initiation parameter in the models.
#' @param chains Number of MCMC chains.
#' @param iter Number of iterations (or output samples for VI).
#' @param seed Random seed for reproducibility.
#' @param cores Number of CPU cores for parallelization.
#' @param comparison Criterion for model selection, either `"loo"` or `"bic"`.
#' @param models_to_fit Models to fit, default `c("exponential", "logistic", "gompertz")`.
#' @param method Fitting method, `"sampling"` or `"vi"`.
#' @param use_elbo Logical; if `TRUE` and `method="vi"`, use the ELBO for model comparison.
#'
#' @return A list containing:
#'   \item{best_model}{The name of the best model.}
#'   \item{fit}{Parsed Stan fit object for the best model.}
#'   \item{model_table}{Model comparison table.}
#'   \item{criterion}{The model selection criterion used.}
#'   \item{method}{The fitting method used.}
#'   \item{breakpoints}{Breakpoints used in the fit.}
#'
#' @export
fit_growth <- function(data,
                       breakpoints = numeric(0),
                       with_initiation = TRUE,
                       chains = 4,
                       iter = 2000,
                       seed = 123,
                       cores = 4,
                       comparison = c("loo", "bic"),
                       models_to_fit = c("exponential", "logistic", "gompertz"),
                       method = c("sampling", "vi"),
                       use_elbo = FALSE) {

  comparison <- match.arg(comparison)
  method <- match.arg(method)

  res <- if (method == "sampling") {
    fit_growth_models(
      data = data, breakpoints = breakpoints, with_initiation = with_initiation,
      chains = chains, iter = iter, seed = seed, cores = cores,
      comparison = comparison, models_to_fit = models_to_fit
    )
  } else {
    fit_growth_models_VI(
      data = data, breakpoints = breakpoints, with_initiation = with_initiation,
      chains = chains, iter = iter, seed = seed, cores = cores,
      comparison = comparison, models_to_fit = models_to_fit,
      method = "vi", use_elbo = use_elbo
    )
  }

  res$model_table$model = rownames(res$model_table)

  if (res$criterion %in% c("bic", "elbo")) {
    best_model <- rownames(res$model_table)[which.min(res$model_table[[1]])]
  } else if (res$criterion == "loo") {
    best_model <- as.character(res$model_table$model[1])
  } else {
    stop("Unknown criterion: cannot select best model")
  }

  best_fit <- res$fits[[best_model]]

  list(
    best_model = best_model,
    fit = parse_stan_fit(best_fit),
    model_table = res$model_table,
    criterion = res$criterion,
    method = method,
    breakpoints = breakpoints
  )
}

#' Fit growth models via MCMC sampling
#'
#' Internal function that fits multiple growth models to the data using MCMC sampling.
#'
#' @inheritParams fit_growth
#' @return A list containing model fits, comparisons, and model selection metrics.
fit_growth_models <- function(data, breakpoints, with_initiation = TRUE,
                              chains = 4, iter = 2000, seed = 123, cores = 4,
                              comparison = c("loo", "bic"),
                              models_to_fit = c("exponential", "logistic", "gompertz")) {
  comparison <- match.arg(comparison)
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]

  t_array <- if (length(breakpoints) == 0) array(0, dim = 0) else as.vector(breakpoints)
  G <- length(breakpoints) + 1

  stan_data <- list(S = nrow(data), G = G, N = data$count, T = data$time, t_array = t_array)
  model_files <- paste0(models_to_fit, if (with_initiation) "_with_init" else "_no_init")

  fits <- list()
  comparisons <- list()
  info_criteria <- numeric(length(models_to_fit))
  names(info_criteria) <- models_to_fit

  for (i in seq_along(models_to_fit)) {
    model_name <- models_to_fit[i]
    model <- get_model(model_files[i])

    message(sprintf("Fitting model: %s", model_name))
    fit <- suppressMessages(suppressWarnings(model$sample(
      data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
      seed = seed, parallel_chains = cores, refresh = 0
    )))

    fits[[model_name]] <- fit
    draws <- fit$draws(format = "draws_matrix")
    log_lik <- draws[, grep("^log_lik", colnames(draws)), drop = FALSE]

    if (comparison == "loo") {
      comparisons[[model_name]] <- loo::loo(log_lik)
      info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
    } else {
      log_mean_lik <- apply(log_lik, 2, mean)
      log_lik_sum <- sum(log_mean_lik)
      param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
      n_params <- length(param_cols)
      N <- length(data$count)
      info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
    }
  }

  comp_table <- if (comparison == "bic") {
    data.frame(BIC = sort(info_criteria))
  } else {
    tbl <- loo::loo_compare(comparisons)
    tbl <- as.data.frame(tbl)
    tbl$model <- rownames(tbl)
    tbl[, c("model", setdiff(names(tbl), "model"))]
  }
  comp_table$model = models_to_fit

  list(fits = fits, comparisons = comparisons, model_table = comp_table, criterion = comparison)
}

#' Fit growth models using variational inference (VI)
#'
#' Internal function that fits models using variational inference and optionally uses ELBO for model selection.
#'
#' @inheritParams fit_growth
#' @return A list containing model fits, comparisons, and model selection metrics.
fit_growth_models_VI <- function(data, breakpoints, with_initiation = TRUE,
                                 chains = 4, iter = 2000, seed = 123, cores = 4,
                                 comparison = c("loo", "bic"),
                                 models_to_fit = c("exponential", "logistic", "gompertz"),
                                 method = c("sampling", "vi"),
                                 use_elbo = FALSE) {

  comparison <- match.arg(comparison)
  method <- match.arg(method)

  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]

  t_array <- if (length(breakpoints) == 0) array(0, dim = 0) else as.vector(breakpoints)
  G <- length(breakpoints) + 1

  stan_data <- list(
    S = nrow(data),
    G = G,
    N = data$count,
    T = data$time,
    t_array = t_array
  )

  model_files <- paste0(models_to_fit, if (with_initiation) "_with_init" else "_no_init")

  fits <- list()
  comparisons <- list()
  info_criteria <- numeric(length(models_to_fit))
  names(info_criteria) <- models_to_fit

  for (i in seq_along(models_to_fit)) {
    model_name <- models_to_fit[i]
    model = get_model(model_files[i])
    # model_path <- file.path(stan_dir, model_files[i])
    # model <- cmdstanr::cmdstan_model(model_path)

    message(sprintf("Fitting model: %s (%s)", model_name, method))

    if (method == "sampling") {
      fit <- suppressMessages(
        suppressWarnings(
          model$sample(
            data = stan_data,
            chains = chains,
            iter_warmup = iter,
            iter_sampling = iter,
            seed = seed,
            parallel_chains = cores,
            refresh = 0
          )
        )
      )

      draws <- fit$draws(format = "draws_matrix")
      log_lik <- draws[, grep("^log_lik", colnames(draws)), drop = FALSE]

      if (comparison == "loo") {
        comparisons[[model_name]] <- loo::loo(log_lik)
        info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
      } else if (comparison == "bic") {
        log_mean_lik <- apply(log_lik, 2, mean)
        log_lik_sum <- sum(log_mean_lik)
        param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
        n_params <- length(param_cols)
        N <- length(data$count)
        info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
      }

    } else if (method == "vi") {
      fit <- suppressMessages(
        suppressWarnings(
          model$variational(
            data = stan_data,
            seed = seed,
            output_samples = iter
          )
        )
      )

      draws <- fit$draws(format = "draws_matrix")
      log_lik_idx <- grep("^log_lik", colnames(draws))

      if (use_elbo) {
        elbo <- fit$metadata()$elbo
        info_criteria[i] <- -elbo
        comparisons[[model_name]] <- list(elbo = elbo)
      } else if (length(log_lik_idx) > 0) {
        log_lik <- draws[, log_lik_idx, drop = FALSE]

        if (comparison == "loo") {
          comparisons[[model_name]] <- loo::loo(log_lik)
          info_criteria[i] <- comparisons[[model_name]]$estimates["looic", "Estimate"]
        } else if (comparison == "bic") {
          log_mean_lik <- apply(log_lik, 2, mean)
          log_lik_sum <- sum(log_mean_lik)
          param_cols <- grep("^rho\\[|^n0$|^t0$|^K$", colnames(draws), value = TRUE)
          n_params <- length(param_cols)
          N <- length(data$count)
          info_criteria[i] <- -2 * log_lik_sum + log(N) * n_params
        }
      } else {
        warning(paste0("No log_lik found in draws for model ", model_name," - cannot compute ", comparison))
        info_criteria[i] <- NA
        comparisons[[model_name]] <- NULL
      }
    }

    fits[[model_name]] <- fit
  }

  comp_table <- if (use_elbo) {
    data.frame(ELBO = sort(info_criteria, na.last = TRUE))
  } else if (comparison == "bic") {
    data.frame(BIC = sort(info_criteria, na.last = TRUE))
  } else if (comparison == "loo") {
    valid_comparisons <- comparisons[!vapply(comparisons, is.null, logical(1))]
    tbl <- loo::loo_compare(valid_comparisons)
    tbl <- as.data.frame(tbl)
    tbl$model <- rownames(tbl)
    tbl <- tbl[, c("model", setdiff(names(tbl), "model"))]
    tbl
  } else {
    NULL
  }

  list(
    fits = fits,
    comparisons = comparisons,
    model_table = comp_table,
    criterion = if (use_elbo) "elbo" else comparison,
    method = method
  )
}

#' Fit and compare tumor recovery models
#'
#' Fits and compares multiple recovery models (two-population de-novo, pre-existing, and single-population)
#' to determine the best one using either LOO or BIC as selection criterion.
#'
#' @param data Data frame with `time` and `count` columns.
#' @param chains Number of MCMC chains.
#' @param iter Number of iterations.
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param comparison Criterion for model selection: `"loo"` or `"bic"`.
#'
#' @return A list containing:
#'   \item{best_model}{Name of the best recovery model.}
#'   \item{best_fit}{Parsed Stan fit object for the best model.}
#'   \item{all_fits}{List of all fitted models.}
#'   \item{model_table}{Comparison table with IC values.}
#'   \item{criterion}{Criterion used for model selection.}
#'
#' @examples
#' \dontrun{
#'   fit_best_recovery_model(my_data)
#' }
#' @export
fit_best_recovery_model <- function(data,
                                    chains = 4,
                                    iter = 4000,
                                    seed = 123,
                                    cores = 4,
                                    comparison = c("loo", "bic")) {
  comparison <- match.arg(comparison)
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]

  stan_data <- list(S = nrow(data), N = data$count, T = data$time)
  model_files <- c("two_pop_both", "two_pop_single")

  fits <- list()
  ic_values <- numeric(length(model_files))

  for (i in seq_along(model_files)) {
    mod <- biPOD:::get_model(model_files[i])
    fit <- suppressMessages(suppressWarnings(mod$sample(
      data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
      seed = seed, parallel_chains = cores, refresh = 0
    )))
    fits[[i]] <- fit

    log_lik <- fit$draws("log_lik") %>% posterior::as_draws_matrix()
    if (comparison == "loo") {
      ic_values[i] <- loo::loo(log_lik)$estimates["elpd_loo", "Estimate"]
    } else {
      log_lik_sum <- sum(apply(log_lik, 1, mean))
      k <- length(fit$metadata()$parameters)
      n <- nrow(data)
      ic_values[i] <- -2 * log_lik_sum + k * log(n)
    }
  }

  best_idx <- if (comparison == "loo") which.max(ic_values) else which.min(ic_values)
  #all_fits <- lapply(fits, biPOD:::parse_stan_fit)

  if (model_files[best_idx] == "two_pop_both") {
    draws = fits[[best_idx]]$draws(format = "draws_list")
    t0_draws = unlist(lapply(draws, function(c) {c[["t0_r"]]}))

    if (stats::median(t0_draws) <= min(data$time)) {
      mod = biPOD:::get_model("two_pop_preexisting")
      best_model = "pre-existing"
    } else {
      mod = biPOD:::get_model("two_pop_denovo")
      best_model = "de-novo"
    }

    best_fit = suppressMessages(suppressWarnings(mod$sample(
      data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
      seed = seed, parallel_chains = cores, refresh = 0
    )))
  } else {
    best_fit = fits[[best_idx]]
    best_model = "single-pop"
  }

  best_fit <- biPOD:::parse_stan_fit(best_fit)

  model_table <- if (comparison == "loo") {
    dplyr::tibble(model = model_files, LOO = ic_values)
  } else {
    dplyr::tibble(model = model_files, BIC = ic_values)
  }

  list(best_model = best_model, first_fit = biPOD:::parse_stan_fit(fits[[best_idx]]),
       final_fit = best_fit, model_table = model_table, criterion = comparison)
}
