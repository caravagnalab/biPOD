#' Fit the held-phase model
#'
#' Model for a continuous marker that has a non-tumour background, after an
#' intervention that lowers it. The observed value is `b + T(t)`, a background `b`
#' plus a tumour-derived part `T(t)` that declines at rate `kill` until `t_n`,
#' then either regrows at rate `r` straight away (`held = FALSE`, continuous
#' regrowth) or stays flat until an exit time `t_e` and regrows at `r` after it
#' (`held = TRUE`). Noise is Student-t (5 degrees of freedom) on the log scale.
#'
#' The decline is allowed only until `t_decline_max` (for example the end
#' of the intervention plus a lag), since the decline needs it. Without this
#' bound, a deep dip below the background can mimic any flat stretch.
#'
#' Read as continuous growth, a held phase of length `H` followed by regrowth at
#' rate `r` is equivalent to regrowth from a population that was only
#' `exp(-r * H)` of the residual disease. This seed fraction is returned as a
#' derived quantity; the data cannot tell a true hold from a small seed.
#'
#' Priors: `kill ~ lognormal(log 0.03, 1)`, `r ~ lognormal(log 0.01, 1)` (per
#' time unit; suited to days), `log b ~ normal(log(background_centre * min(y)),
#' background_sd)`, `sigma ~ half-normal(0, 0.3)`, exit time uniform over the
#' remaining window.
#'
#' @param data Data frame with `time` and `count` (the marker value, > 0).
#' @param t_decline_max Latest time at which the decline may end.
#' @param held Fit the held-phase model (`TRUE`) or continuous regrowth (`FALSE`).
#' @param background_centre Prior centre of the background, as a multiple of the
#'   lowest observed value.
#' @param background_sd Prior sd of the log background.
#' @param weights Likelihood weights per observation (0 leaves one out); used by
#'   [compare_held_phase()] for exact leave-one-out refits.
#' @param chains Number of MCMC chains.
#' @param iter Number of sampling iterations per chain (warmup is 1.5 times longer).
#' @param seed Random seed.
#' @param cores Number of CPU cores.
#' @param adapt_delta Target acceptance rate for NUTS.
#' @param max_rhat Rhat threshold used for the `converged` flag.
#'
#' @return A list containing:
#'   \item{parameters}{Tibble with the posterior median and 90% interval of
#'     `kill`, `r`, `t_n`, `t_e`, `held_days`, `background` and, for
#'     `held = TRUE`, `log10_seed_fraction`.}
#'   \item{fit}{Parsed Stan fit.}
#'   \item{log_lik}{Draws by observations matrix of pointwise log-likelihoods.}
#'   \item{max_rhat}{Maximum Rhat over the sampled parameters.}
#'   \item{converged}{Whether `max_rhat` is below the threshold.}
#'
#' @examples
#' \dontrun{
#'   fit_held_phase(marker, t_decline_max = 170)
#' }
#' @export
fit_held_phase <- function(data,
                           t_decline_max,
                           held = TRUE,
                           background_centre = 1,
                           background_sd = 0.5,
                           weights = NULL,
                           chains = 4,
                           iter = 1000,
                           seed = 123,
                           cores = 4,
                           adapt_delta = 0.97,
                           max_rhat = 1.05) {
  stopifnot(all(c("time", "count") %in% colnames(data)))
  data <- data[order(data$time), ]
  if (any(data$count <= 0)) stop("'count' must be positive for the held-phase model")
  if (anyDuplicated(data$time)) stop("Times must be unique")
  t <- data$time
  y <- data$count
  S <- length(t)
  if (S < 5) stop("At least five observations are needed")
  if (t_decline_max <= t[1]) stop("'t_decline_max' must be after the first observation")
  if (is.null(weights)) weights <- rep(1, S)

  stan_data <- list(
    S = S, t = t, y = y, t_n_max = t_decline_max, held = as.integer(held), w = weights,
    log_b_prior = log(background_centre * min(y)), log_b_sd = background_sd,
    log_T0_prior = log(max(y[1] - 0.5 * min(y), 1))
  )
  inits <- held_phase_inits(t, y, t_decline_max, held, chains, stan_data$log_T0_prior)

  mod <- get_model("held_phase", "lognormal")
  out_dir <- tempfile("biPOD_held_")
  dir.create(out_dir)
  fit <- suppressMessages(suppressWarnings(mod$sample(
    data = stan_data, chains = chains, iter_warmup = round(1.5 * iter), iter_sampling = iter,
    seed = seed, parallel_chains = cores, refresh = 0, init = inits,
    adapt_delta = adapt_delta, max_treedepth = 12, output_dir = out_dir,
    show_messages = FALSE, show_exceptions = FALSE
  )))

  dr <- posterior::as_draws_df(fit$draws(c("kill", "r", "t_n", "t_e", "held_days", "log_b")))
  vals <- list(kill = dr$kill, r = dr$r, t_n = dr$t_n, t_e = dr$t_e,
               held_days = dr$held_days, background = exp(dr$log_b))
  if (held) vals$log10_seed_fraction <- -dr$r * dr$held_days / log(10)
  params <- dplyr::bind_rows(lapply(names(vals), function(v) {
    qv <- stats::quantile(vals[[v]], c(0.5, 0.05, 0.95), names = FALSE)
    dplyr::tibble(parameter = v, median = qv[1], q5 = qv[2], q95 = qv[3])
  }))

  sampled <- c("log_T0", "kill", "r", "t_n", "log_b", "sigma", if (held) "u")
  rhat <- max(fit$summary(sampled)$rhat, na.rm = TRUE)
  if (rhat >= max_rhat) {
    cli::cli_alert_warning("Held-phase model did not converge (max Rhat {round(rhat, 2)}).")
  }

  list(parameters = params, fit = parse_stan_fit(fit),
       log_lik = posterior::as_draws_matrix(fit$draws("log_lik")),
       max_rhat = rhat, converged = rhat < max_rhat)
}

#' Compare continuous regrowth with a held phase
#'
#' Fits [fit_held_phase()] with and without a held phase and compares them by
#' PSIS leave-one-out cross-validation on the pointwise log-likelihood.
#' Observations with Pareto k above `exact_k` are refitted exactly with that
#' observation left out, so unreliable importance weights cannot decide the
#' comparison.
#'
#' Evidence for a held phase is `"needed"` when the elpd difference (held minus
#' continuous) exceeds 2 standard errors, `"suggestive"` when it is positive and
#' the 5% quantile of the held phase exceeds `min_held`, and `"not needed"`
#' otherwise.
#'
#' @inheritParams fit_held_phase
#' @param exact_k Pareto k above which an observation is refitted exactly.
#' @param min_held Minimum held phase (5% quantile) for `"suggestive"` evidence.
#' @param ... Further arguments passed to [fit_held_phase()].
#'
#' @return A list containing:
#'   \item{elpd_diff}{elpd of the held-phase model minus continuous regrowth.}
#'   \item{se}{Standard error of `elpd_diff`.}
#'   \item{evidence}{`"needed"`, `"suggestive"` or `"not needed"`.}
#'   \item{held}{Output of [fit_held_phase()] with `held = TRUE`.}
#'   \item{continuous}{Output of [fit_held_phase()] with `held = FALSE`.}
#'   \item{n_exact_refits}{Number of exact leave-one-out refits.}
#'
#' @examples
#' \dontrun{
#'   compare_held_phase(marker, t_decline_max = 170)
#' }
#' @export
compare_held_phase <- function(data, t_decline_max, exact_k = 0.7, min_held = 60, ...) {
  data <- data[order(data$time), ]
  m0 <- fit_held_phase(data, t_decline_max, held = FALSE, ...)
  m1 <- fit_held_phase(data, t_decline_max, held = TRUE, ...)
  l0 <- loo_with_exact_refits(m0, data, t_decline_max, held = FALSE, exact_k = exact_k, ...)
  l1 <- loo_with_exact_refits(m1, data, t_decline_max, held = TRUE, exact_k = exact_k, ...)

  d_i <- l1$elpd - l0$elpd
  diff <- sum(d_i)
  se <- sqrt(length(d_i) * stats::var(d_i))
  held_q5 <- m1$parameters$q5[m1$parameters$parameter == "held_days"]
  evidence <- if (diff > 2 * se) {
    "needed"
  } else if (diff > 0 && held_q5 > min_held) {
    "suggestive"
  } else {
    "not needed"
  }

  list(elpd_diff = diff, se = se, evidence = evidence, held = m1, continuous = m0,
       n_exact_refits = l0$n_refit + l1$n_refit)
}

# Pointwise PSIS-LOO elpd, with exact refits where Pareto k > exact_k
loo_with_exact_refits <- function(m, data, t_decline_max, held, exact_k, ...) {
  ll <- m$log_lik
  n_chains <- posterior::nchains(posterior::as_draws_array(m$fit$draws))
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(seq_len(n_chains), each = nrow(ll) / n_chains))
  lo <- suppressWarnings(loo::loo(ll, r_eff = r_eff))
  elpd <- lo$pointwise[, "elpd_loo"]
  refit <- which(lo$diagnostics$pareto_k > exact_k)
  for (i in refit) {
    w <- rep(1, nrow(data)); w[i] <- 0
    mi <- fit_held_phase(data, t_decline_max, held = held, weights = w, ...)
    lli <- mi$log_lik[, i]
    elpd[i] <- max(lli) + log(mean(exp(lli - max(lli))))
  }
  list(elpd = elpd, n_refit = length(refit))
}

# Initial values: nadir at the decline bound, exit at the last sample near the minimum
held_phase_inits <- function(t, y, t_decline_max, held, chains, log_T0) {
  S <- length(t)
  near_min <- which(y <= 1.5 * min(y))
  t_e0 <- t[max(near_min)]
  t_n0 <- min(t_decline_max - 1e-3, max(t[1] + (t_decline_max - t[1]) / 2, t_decline_max - 30))
  u0 <- if (held) min(0.99, max(0.01, (t_e0 - t_n0) / (t[S] - t_n0))) else 0.5
  lapply(seq_len(chains), function(i) {
    list(log_T0 = log_T0, kill = 0.03, r = 0.01, t_n = t_n0, u = u0,
         log_b = log(min(y)), sigma = 0.15)
  })
}
