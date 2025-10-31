
stan_qc <- function(model, fit, stan_data,
                    y_name      = "count",
                    loglik_name = "log_lik",
                    thr_rhat         = 1.01,
                    thr_ess_min      = 200,      # bulk ESS minimum
                    require_no_div   = TRUE,
                    require_no_tdhit = TRUE,
                    thr_pareto_k     = 1.0,
                    thr_info_gain    = 1.1,
                    thr_z_shift_abs  = 0.25,
                    max_pareto_k     = 1
) {
  out <- list()

  draws_df <- posterior::as_draws_df(fit$draws())

  # ---------------- Convergence / NUTS ----------------
  summ <- posterior::summarise_draws(draws_df)
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

  # ---- Prior-only run (unchanged) ----
  stan_data$prior_only = 1
  fit_prior <- model$sample(data = stan_data, iter_sampling = 2000, chains = 1, refresh = 0)
  prior_df <- posterior::as_draws_df(fit_prior$draws())

  # ---- Prior–posterior table: info_gain + z_shift only ----
  pp <- lapply(pars, function(p) {
    pr <- prior_df[[p]]; po <- draws_df[[p]]
    mu_prior <- tryCatch(mean(pr, na.rm = TRUE), error = function(e) NA_real_)
    mu_post  <- tryCatch(mean(po, na.rm = TRUE), error = function(e) NA_real_)
    sd_prior <- tryCatch(sd(pr,   na.rm = TRUE), error = function(e) NA_real_)
    sd_post  <- tryCatch(sd(po,   na.rm = TRUE), error = function(e) NA_real_)
    info_gain <- sd_prior / sd_post
    z_shift <- if (is.finite(sd_prior) && sd_prior > 0) (mu_post - mu_prior) / sd_prior else NA_real_
    data.frame(parameter = p,
               mu_prior = mu_prior, sd_prior = sd_prior,
               mu_post  = mu_post,  sd_post  = sd_post,
               info_gain = info_gain, z_shift = z_shift)
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

  # Identifiability rules (info_gain + z_shift)
  if (!is.null(out$prior_posterior) && nrow(out$prior_posterior) > 0) {
    ig <- out$prior_posterior
    weak_contr_idx <- which(!is.na(ig$info_gain) & ig$info_gain < thr_info_gain)
    weak_shift_idx <- which(!is.na(ig$z_shift) & abs(ig$z_shift) < thr_z_shift_abs)
    to_fail <- intersect(weak_contr_idx, weak_shift_idx)
    if (length(to_fail) > 0) {
      fails <- c(fails, sprintf("Uninformative posterior (info_gain < %.2f & |z_shift| < %.2f) for: %s",
                                thr_info_gain, thr_z_shift_abs,
                                paste(ig$parameter[to_fail], collapse = ", ")))
    }
  }

  out$verdict <- if (length(fails) == 0) "PASS" else "FAIL"
  out$fail_reasons <- fails
  out$checked_parameters <- pars
  class(out) <- c("stan_qc_report", class(out))
  return(out)
}
