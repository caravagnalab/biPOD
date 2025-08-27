
stan_qc <- function(model, fit, stan_data,
                    y_name      = "count",
                    loglik_name = "log_lik",
                    # ----- Thresholds (tweak to taste) -----
                    thr_rhat         = 1.01,
                    thr_ess_min      = 200,      # bulk ESS minimum
                    require_no_div   = TRUE,
                    require_no_tdhit = TRUE,
                    thr_pareto_k     = 1.0,
                    thr_info_gain    = 1.25,
                    max_pareto_k     = 1
) {
  suppressPackageStartupMessages({
    library(cmdstanr); library(posterior); library(bayesplot); library(loo)
  })

  out <- list()
  draws_df <- as_draws_df(fit$draws())

  # ---------------- Convergence / NUTS ----------------
  summ <- summarise_draws(draws_df)
  out$summary <- summ[, c("variable","mean","sd","rhat","ess_bulk","ess_tail")]

  nuts <- fit$sampler_diagnostics()
  div_per_chain   <- colSums(nuts[,,"divergent__"])
  td_max          <- max(nuts[,,"treedepth__"])
  td_hits_per_ch  <- colSums(nuts[,,"treedepth__"] >= td_max)
  out$divergences <- div_per_chain
  out$treedepth_hits <- td_hits_per_ch

  # helper: extract parameter names to evaluate (exclude rep/pred/lik/lp__)
  model_params <- fit$metadata()$model_params
  pars <- grep("(rep|pred|lik)", model_params, invert = TRUE, value = TRUE)
  pars <- setdiff(pars, "lp__")
  # keep only those actually present in draws (defensive)
  pars <- intersect(pars, colnames(draws_df))

  # ---------------- Goodness of fit (LOO) -------------
  out$loo <- NULL
  loo_k <- NULL
  ll_mat <- posterior::as_draws_matrix(fit$draws(loglik_name))
  loo_fit <- loo::loo(ll_mat)
  out$loo <- loo_fit
  loo_k <- as.numeric(loo_fit$diagnostics$pareto_k)

  # ---------------- Identifiability -------------------
  # Prior-only run
  fit$summary()

  stan_data$prior_only = 1
  fit_prior <- model$sample(
    data = stan_data, iter_sampling = 2000, chains = 1, refresh = 0
  )

  prior_df <- as_draws_df(fit_prior$draws())

  # prior–posterior contraction (sd_prior / sd_post)
  pp <- lapply(pars, function(p) {
    sd_prior <- tryCatch(sd(prior_df[[p]], na.rm = TRUE), error = function(e) NA_real_)
    sd_post  <- tryCatch(sd(draws_df[[p]], na.rm = TRUE), error = function(e) NA_real_)
    data.frame(parameter = p, sd_prior = sd_prior, sd_post = sd_post,
               info_gain = sd_prior / sd_post)
  })
  out$prior_posterior <- do.call(rbind, pp)

  # ---------------- Verdict rules ---------------------
  fails <- c()

  # Convergence rules
  bad_rhat <- out$summary$rhat[match(pars, out$summary$variable)]
  if (any(!is.na(bad_rhat) & bad_rhat > thr_rhat)) {
    fails <- c(fails, sprintf("Rhat > %.3f for: %s",
                              thr_rhat,
                              paste(pars[which(bad_rhat > thr_rhat)], collapse=", ")))
  }
  bad_ess <- out$summary$ess_bulk[match(pars, out$summary$variable)]
  if (any(!is.na(bad_ess) & bad_ess < thr_ess_min)) {
    fails <- c(fails, sprintf("bulk ESS < %d for: %s",
                              thr_ess_min,
                              paste(pars[which(bad_ess < thr_ess_min)], collapse=", ")))
  }
  if (require_no_div && any(div_per_chain > 0)) {
    fails <- c(fails, sprintf("NUTS divergences present (per-chain: %s)",
                              paste(div_per_chain, collapse = ",")))
  }
  if (require_no_tdhit && any(td_hits_per_ch > 0)) {
    fails <- c(fails, sprintf("Tree depth saturations present (per-chain: %s)",
                              paste(td_hits_per_ch, collapse=",")))
  }

  # GoF rule (Pareto-k)
  if (!is.null(loo_k) && sum(loo_k > thr_pareto_k, na.rm = TRUE) > max_pareto_k) {
    n_bad <- sum(loo_k > thr_pareto_k, na.rm = TRUE)
    fails <- c(fails, sprintf("LOO Pareto-k > %.2f for %d point(s)", thr_pareto_k, n_bad))
  }

  # Identifiability rules
  if (!is.null(out$prior_posterior) && nrow(out$prior_posterior) > 0) {
    ig <- out$prior_posterior
    bad_ig_idx <- which(!is.na(ig$info_gain) & ig$info_gain < thr_info_gain)
    if (length(bad_ig_idx) > 0) {
      fails <- c(fails, sprintf("Weak prior→posterior contraction (info_gain < %.2f) for: %s",
                                thr_info_gain,
                                paste(ig$parameter[bad_ig_idx], collapse=", ")))
    }
  }

  out$verdict <- if (length(fails) == 0) "PASS" else "FAIL"
  out$fail_reasons <- fails
  out$checked_parameters <- pars
  class(out) <- c("stan_qc_report", class(out))
  return(out)
}
#
# print.stan_qc_report <- function(x, ...) {
#   cat("\n=== Stan QC Report ===\n")
#   cat("Verdict:", x$verdict, "\n")
#   if (length(x$fail_reasons)) {
#     cat("Reasons:\n -", paste(x$fail_reasons, collapse = "\n - "), "\n")
#   }
#   cat("\nConvergence:\n")
#   cat("  Divergences per chain:", paste(x$divergences, collapse=", "), "\n")
#   cat("  Treedepth hits per chain:", paste(x$treedepth_hits, collapse=", "), "\n")
#
#   if (!is.null(x$loo)) {
#     cat("\nLOO:\n")
#     print(x$loo)
#   }
#
#   if (!is.null(x$prior_posterior)) {
#     cat("\nPrior→Posterior contraction (first 10):\n")
#     print(utils::head(x$prior_posterior[order(x$prior_posterior$info_gain), ], 10))
#   }
#
#   if (!is.null(x$posterior_PCs)) {
#     cat("\nPosterior principal directions (variance fractions):\n")
#     print(utils::head(x$posterior_PCs, 5))
#   }
#
#   if (!is.null(x$sensitivity)) {
#     cat("\nSensitivity (PRCC to", deparse(substitute(x$target_name)), "):\n")
#     print(x$sensitivity)
#   }
#   invisible(x)
# }

