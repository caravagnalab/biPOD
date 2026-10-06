
#' Fit growth models and select the best one
#'
#' Fits multiple candidate growth models (exponential, logistic, Gompertz) to time-series
#' count data with optional breakpoints, using either MCMC sampling or variational inference (VI).
#' The best model is selected based on the specified criterion (LOO, BIC, or ELBO).
#'
#' @param x biPOD object
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
fit_growth <- function(x,
                       with_initiation = TRUE,
                       chains = 4,
                       iter = 2000,
                       seed = 123,
                       cores = 4,
                       comparison = c("bic", "loo"),
                       models_to_fit = c("exponential", "logistic", "gompertz", "monomolecular", "quadraticexp"),
                       method = c("sampling", "vi"),
                       noise_model = c("lognormal", "poisson", "negbinomial"),
                       use_elbo = FALSE) {

  data = x$counts
  breakpoints = x$metadata$breakpoints

  comparison <- match.arg(comparison)
  noise_model <- match.arg(noise_model)
  method <- match.arg(method)

  res <- if (method == "sampling") {
    fit_growth_models(
      data = data, breakpoints = breakpoints, with_initiation = with_initiation,
      chains = chains, iter = iter, seed = seed, cores = cores,
      comparison = comparison, models_to_fit = models_to_fit, noise_model = noise_model
    )
  } else {
    fit_growth_models_VI(
      data = data, breakpoints = breakpoints, with_initiation = with_initiation,
      chains = chains, iter = iter, seed = seed, cores = cores,
      comparison = comparison, models_to_fit = models_to_fit, noise_model = noise_model,
      method = "vi", use_elbo = use_elbo
    )
  }

  res$model_table$qc = unlist(lapply(res$model_table$model, function(n) res$fits_qc[[n]]$verdict))

  # Candidate rows: all PASS if any; otherwise all rows (all FAIL case)
  pass_idx <- which(res$model_table$qc == "PASS")
  cand_idx <- if (length(pass_idx) > 0) pass_idx else seq_len(nrow(res$model_table))

  # Pick best among candidates according to criterion
  if (res$criterion %in% c("bic", "elbo")) {
    # Use a named metric column if available; otherwise fall back to first column (as in original)
    metric_values <- if (res$criterion %in% names(res$model_table)) {
      res$model_table[[res$criterion]]
    } else {
      res$model_table[[1]]
    }
    best_idx   <- cand_idx[which.min(metric_values[cand_idx])]
    best_model <- rownames(res$model_table)[best_idx]

  } else if (res$criterion == "loo") {
    # Assume model_table already sorted by LOO (best first), keep first among candidates
    best_model <- as.character(res$model_table$model[cand_idx][1])

  } else {
    stop("Unknown criterion: cannot select best model")
  }

  best_fit <- res$fits[[best_model]]

  x$growth_fit = list(
    best_model = best_model,
    fit = parse_stan_fit(best_fit),
    qc = res$fits_qc[[best_model]],
    model_table = res$model_table,
    criterion = res$criterion,
    method = method,
    breakpoints = breakpoints
  )
  x
}

#' Fit growth models via MCMC sampling
#'
#' Internal function that fits multiple growth models to the data using MCMC sampling.
#'
#' @inheritParams fit_growth
#' @return A list containing model fits, comparisons, and model selection metrics.
fit_growth_models <- function(data, breakpoints, with_initiation = TRUE,
                              chains = 4, iter = 2000, seed = 123, cores = 4,
                              comparison = c("bic", "loo"),
                              models_to_fit = c("exponential", "logistic", "gompertz"),
                              noise_model = c("lognormal", "poisson", "negbinomial")) {
  comparison <- match.arg(comparison)
  noise_model <- match.arg(noise_model)
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]

  t_array <- if (length(breakpoints) == 0) array(0, dim = 0) else as.vector(breakpoints)
  G <- length(breakpoints) + 1

  stan_data <- list(S = nrow(data), G = G, N = data$count, T = data$time, t_array = t_array, prior_only = 0)
  model_files <- paste0(models_to_fit, if (with_initiation) "_with_init" else "_no_init")

  fits <- list()
  fits_qc <- list()
  comparisons <- list()
  info_criteria <- numeric(length(models_to_fit))
  names(info_criteria) <- models_to_fit

  for (i in seq_along(models_to_fit)) {
    model_name <- models_to_fit[i]
    model <- get_model(model_files[i], noise_model)

    stan_data$prior_only <- 0
    message(sprintf("Fitting model: %s", model_name))
    fit <- suppressMessages(suppressWarnings(model$sample(
      data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
      seed = seed, parallel_chains = cores, refresh = 0
    )))

    qc_result <- stan_qc(model, fit, stan_data, require_no_tdhit = F, require_no_div = F)

    fits[[model_name]] <- fit
    fits_qc[[model_name]] <- qc_result
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
    tbl = data.frame(BIC = info_criteria)
    tbl$model = models_to_fit
    tbl
  } else {
    tbl <- loo::loo_compare(comparisons)
    tbl <- as.data.frame(tbl)
    tbl$model <- rownames(tbl)
    tbl[, c("model", setdiff(names(tbl), "model"))]
  }

  list(fits = fits, fits_qc = fits_qc, comparisons = comparisons, model_table = comp_table, criterion = comparison)
}

#' Fit growth models using variational inference (VI)
#'
#' Internal function that fits models using variational inference and optionally uses ELBO for model selection.
#'
#' @inheritParams fit_growth
#' @return A list containing model fits, comparisons, and model selection metrics.
fit_growth_models_VI <- function(data, breakpoints, with_initiation = TRUE,
                                 chains = 4, iter = 2000, seed = 123, cores = 4,
                                 comparison = c("bic", "loo"),
                                 models_to_fit = c("exponential", "logistic", "gompertz"),
                                 noise_model = c("lognormal", "poisson", "negbinomial"),
                                 method = c("sampling", "vi"),
                                 use_elbo = FALSE) {

  comparison <- match.arg(comparison)
  noise_model <- match.arg(noise_model)
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
  fits_qc <- list()
  comparisons <- list()
  info_criteria <- numeric(length(models_to_fit))
  names(info_criteria) <- models_to_fit

  for (i in seq_along(models_to_fit)) {
    model_name <- models_to_fit[i]
    model = get_model(model_files[i], noise_model)
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
    fits_qc[[model_name]] <- list(verdict = "SKIPPED", fail_reasons = "QC diagnostics not available for VI fits")
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
    fits_qc = fits_qc,
    comparisons = comparisons,
    model_table = comp_table,
    criterion = if (use_elbo) "elbo" else comparison,
    method = method
  )
}

#' Initial values for the recovery models
#'
#' Rates are taken from the decline to the minimum (sensitive) and the rise from
#' the minimum to the last observation (resistant); `t_end` and `t0_r` are set so
#' that the curves pass through the first and last observations, then clipped
#' to each model's constraints. Chains after the first are jittered.
#'
#' @param data Data frame with `time` and `count`, sorted by time.
#' @param model_name One of the recovery model names.
#' @param noise_model `"poisson"` or `"negbinomial"`.
#' @param chain_id Chain index.
#' @return A named list of initial values for one chain.
recovery_inits <- function(data, model_name, noise_model, chain_id = 1) {
  t <- data$time
  y <- pmax(data$count, 1)
  S <- length(t)
  i_min <- which.min(y)
  jit <- if (chain_id == 1) 1 else exp(stats::rnorm(1, 0, 0.3))
  shift <- if (chain_id == 1) 0 else stats::rnorm(1, 0, 0.01 * (t[S] - t[1]))

  rho_s <- max(1e-3, log(y[1] / y[i_min]) / max(t[i_min] - t[1], 1)) * jit
  rho_r <- max(1e-3, log(y[S] / y[i_min]) / max(t[S] - t[i_min], 1)) * jit
  t_end <- t[1] + max(log(y[1]), 0.5) / rho_s + abs(shift)
  t0_r <- t[S] - max(log(y[S]), 0.5) / rho_r + shift
  rho_single <- max(1e-3, log(y[S] / y[1]) / max(t[S] - t[1], 1)) * jit
  rho_decay <- max(1e-3, log(y[1]) / max(t[S] - t[1], 1)) * jit

  init <- switch(model_name,
    two_pop_both        = list(rho_s = rho_s, rho_r = rho_r, t0_r = t0_r, t_end = t_end),
    two_pop_preexisting = list(rho_s = rho_s, rho_r = rho_r, t0_r = min(t0_r, t[1] - 1), t_end = t_end),
    two_pop_denovo      = list(rho_s = rho_s, rho_r = rho_r, t0_r = max(t0_r, t[1] + 1), t_end = t_end),
    two_pop_single      = list(rho_r = rho_single, t0_r = t[1] - 1),
    single_pop_decay    = list(rho_s = rho_decay, t_end = t[1] + max(log(y[1]), 0.5) / rho_decay)
  )
  if (noise_model == "negbinomial") init$phi <- 10
  init
}

# Sampled parameters of each recovery model (the BIC penalty counts these only)
recovery_model_params <- function(model_name, noise_model) {
  p <- switch(model_name,
    two_pop_both = , two_pop_preexisting = , two_pop_denovo = c("rho_s", "rho_r", "t0_r", "t_end"),
    two_pop_single = c("rho_r", "t0_r"),
    single_pop_decay = c("rho_s", "t_end")
  )
  if (noise_model == "negbinomial") p <- c(p, "phi")
  p
}

sample_recovery_model <- function(model_name, stan_data, data, noise_model, chains,
                                  iter, seed, cores) {
  mod <- get_model(model_name, noise_model)
  run <- function(seed) {
    set.seed(seed)
    inits <- lapply(seq_len(chains), function(i) recovery_inits(data, model_name, noise_model, i))
    # cmdstanr names output CSVs by model and minute; a unique directory per
    # run stops a garbage-collected earlier fit from deleting these files
    out_dir <- tempfile("biPOD_recovery_")
    dir.create(out_dir)
    fit <- suppressMessages(suppressWarnings(mod$sample(
      data = stan_data, chains = chains, iter_warmup = iter, iter_sampling = iter,
      seed = seed, parallel_chains = cores, refresh = 0, init = inits,
      output_dir = out_dir, show_messages = FALSE
    )))
    fit$draws("lp__")
    fit
  }
  # one retry with another seed if a chain fails
  tryCatch(run(seed), error = function(e) run(seed + 1000))
}

# Fits one candidate. A fit that has not converged is refitted once with twice
# the iterations and another seed, and the better-converged fit is kept.
fit_recovery_candidate <- function(model_name, stan_data, data, noise_model, chains,
                                   iter, seed, cores, max_rhat) {
  params <- recovery_model_params(model_name, noise_model)
  fit <- sample_recovery_model(model_name, stan_data, data, noise_model, chains, iter, seed, cores)
  rhat <- max(fit$summary(params)$rhat, na.rm = TRUE)
  refit <- rhat >= max_rhat
  if (refit) {
    fit2 <- sample_recovery_model(model_name, stan_data, data, noise_model, chains,
                                  2 * iter, seed + 1, cores)
    rhat2 <- max(fit2$summary(params)$rhat, na.rm = TRUE)
    if (rhat2 < rhat) {
      fit <- fit2
      rhat <- rhat2
    }
  }
  list(fit = fit, rhat = rhat, refit = refit)
}

#' Fit and compare tumor recovery models
#'
#' Fits and compares three candidate recovery models: a two-population mixture
#' (`two_pop_both`, a decaying "sensitive" population plus a growing "resistant"
#' population, i.e. relapse), a single growing population (`two_pop_single`, no
#' mixture detected), and a single shrinking-only population (`single_pop_decay`,
#' modeling a tumor that never relapses). The best of the three is selected using
#' either LOO or BIC. If the two-population mixture wins, a further posterior-median
#' check on the resistant clone's birth time (`t0_r`) refits with either the
#' de-novo or pre-existing resistant-clone model.
#'
#' Every chain is initialised from data-driven values (see `recovery_inits`).
#' BIC is computed from the maximum log-likelihood over draws and the number of
#' sampled parameters, so a chain stuck in a poor mode cannot decide the model
#' choice. A candidate whose fit has not converged (Rhat at or above
#' `max_rhat`) is refitted once with twice the iterations and another seed.
#' Selection still uses the criterion over all candidates: their BIC rests on
#' the best draw, which non-convergence barely affects. Non-convergence makes
#' the parameters of the chosen model unreliable instead, so `converged`
#' requires both the chosen candidate and the final fit to have converged.
#'
#' @param data Data frame with `time` and `count` columns.
#' @param noise_model `"poisson"` or `"negbinomial"`.
#' @param chains Number of MCMC chains.
#' @param iter Number of warmup and of sampling iterations per chain.
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param comparison Criterion for model selection: `"bic"` or `"loo"`.
#' @param max_rhat Rhat threshold for refitting a candidate and for the
#'   `converged` flag.
#'
#' @return A list containing:
#'   \item{best_model}{Name of the best recovery model: `"pre-existing"`, `"de-novo"`,
#'     `"single-pop-growing"`, or `"single-pop-shrinking"`.}
#'   \item{first_fit}{Parsed Stan fit of the selected candidate model.}
#'   \item{final_fit}{Parsed Stan fit of the final model.}
#'   \item{model_table}{Comparison table of all candidates: IC value, Rhat,
#'     and whether the candidate was refitted (`refit`).}
#'   \item{criterion}{Criterion used for model selection.}
#'   \item{max_rhat}{Maximum Rhat of the final fit's parameters.}
#'   \item{winner_rhat}{Maximum Rhat of the chosen candidate's fit.}
#'   \item{n_candidates_converged}{Number of candidates with Rhat below the
#'     threshold.}
#'   \item{converged}{Whether both `winner_rhat` and `max_rhat` are below the
#'     threshold.}
#'
#' @examples
#' \dontrun{
#'   fit_best_recovery_model(my_data)
#' }
#' @export
fit_best_recovery_model <- function(data,
                                    noise_model = c("poisson", "negbinomial"),
                                    chains = 4,
                                    iter = 2000,
                                    seed = 123,
                                    cores = 4,
                                    comparison = c("bic", "loo"),
                                    max_rhat = 1.05) {
  comparison <- match.arg(comparison)
  noise_model <- match.arg(noise_model)
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]
  n <- nrow(data)

  stan_data <- list(S = n, N = as.integer(round(data$count)), T = data$time, prior_only = 0)
  model_files <- c("two_pop_both", "two_pop_single", "single_pop_decay")

  fits <- list()
  model_table <- dplyr::bind_rows(lapply(model_files, function(m) {
    cand <- fit_recovery_candidate(m, stan_data, data, noise_model, chains, iter, seed,
                                   cores, max_rhat)
    fit <- cand$fit
    fits[[m]] <<- fit
    log_lik <- posterior::as_draws_matrix(fit$draws("log_lik"))
    params <- recovery_model_params(m, noise_model)
    ic <- if (comparison == "loo") {
      suppressWarnings(loo::loo(log_lik)$estimates["elpd_loo", "Estimate"])
    } else {
      -2 * max(rowSums(log_lik)) + length(params) * log(n)
    }
    row <- dplyr::tibble(model = m, ic = ic, rhat = cand$rhat, refit = cand$refit)
    names(row)[2] <- toupper(comparison)
    row
  }))

  ic_values <- model_table[[toupper(comparison)]]
  best <- model_files[if (comparison == "loo") which.max(ic_values) else which.min(ic_values)]

  best_fit <- fits[[best]]
  final_name <- best
  best_model <- switch(best, two_pop_single = "single-pop-growing",
                       single_pop_decay = "single-pop-shrinking", NA_character_)
  if (best == "two_pop_both") {
    t0 <- stats::median(posterior::as_draws_matrix(best_fit$draws("t0_r")))
    final_name <- if (t0 <= min(data$time)) "two_pop_preexisting" else "two_pop_denovo"
    best_model <- if (final_name == "two_pop_preexisting") "pre-existing" else "de-novo"
    best_fit <- sample_recovery_model(final_name, stan_data, data, noise_model, chains,
                                      2 * iter, seed, cores)
  }

  winner_rhat <- model_table$rhat[model_table$model == best]
  if (winner_rhat >= max_rhat) {
    cli::cli_alert_warning("Selected candidate {best} did not converge (max Rhat {round(winner_rhat, 2)}).")
  }
  final_rhat <- max(best_fit$summary(recovery_model_params(final_name, noise_model))$rhat, na.rm = TRUE)
  if (final_rhat >= max_rhat) {
    cli::cli_alert_warning("Recovery model {final_name} did not converge (max Rhat {round(final_rhat, 2)}).")
  }

  list(best_model = best_model, first_fit = parse_stan_fit(fits[[best]]),
       final_fit = parse_stan_fit(best_fit), model_table = model_table, criterion = comparison,
       max_rhat = final_rhat, winner_rhat = winner_rhat,
       n_candidates_converged = sum(model_table$rhat < max_rhat),
       converged = winner_rhat < max_rhat && final_rhat < max_rhat)
}
