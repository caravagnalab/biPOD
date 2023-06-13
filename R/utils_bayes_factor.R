fit_with_bayes_factor <- function(x,
                                 factor_size = 1,
                                 variational = FALSE,
                                 t0_lower_bound = -10,
                                 prior_K = NULL,
                                 chains = 4, iter = 4000, cores = 4) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size, prior_K = prior_K, t0_lower_bound = t0_lower_bound)

  res_exp <- fit_data(
    x = x, growth_type = "exponential", factor_size = factor_size,
    variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
    chains = chains, iter = iter, cores = cores
  )

  res_log <- fit_data(
    x = x, growth_type = "logistic", factor_size = factor_size,
    variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
    chains = chains, iter = iter, cores = cores
  )

  fits <- list(res_exp$fit, res_log$fit)

  # Compute Marginal Likelihoods for every model

  # Approximated version
  marginal_likelihoods <- lapply(fits, function(f) {
    # b <- bridgesampling::bridge_sampler(f, silent = TRUE)
    f$lp() %>% stats::median()
  })

  # Exact version
  # marginal_likelihoods <- lapply(fits, function(f) {
  #   b <- bridgesampling::bridge_sampler(f, silent = TRUE)
  #   unname(b$logml)
  # })

  marginal_likelihoods <- marginal_likelihoods %>% unlist()
  # Compute pairwise Bayes factor
  bayes_factors <- outer(marginal_likelihoods, marginal_likelihoods, function(x, y) {
    return(exp(x - y))
  })
  bayes_factors <- dplyr::as_tibble(bayes_factors, .name_repair = "minimal")

  if (bayes_factors[1, 2] > 1) {
    K <- bayes_factors[1, 2]
    best_growth <- "Exponential"
    evidence <- bayes_factor_evidence(K)
    res <- res_exp
  } else {
    K <- bayes_factors[2, 1]
    best_growth <- "Logistic"
    evidence <- bayes_factor_evidence(K)
    res <- res_log
  }

  res$fit_info$best_growth <- best_growth
  res$fit_info$bayes_factor <- K
  res$fit_info$evidence <- evidence
  res$fit_info$model_selection_algo <- "bayes_factor"

  cli::cli_alert_info("Model selection finished!")
  cli::cli_alert_info("Model with {.val {best_growth}} growth deemed better with {.val {evidence}} evidence. (BF = {.val {K}})")

  return(res)
}

bayes_factor_evidence <- function(K) {
  if (K < 1) {
    evidence <- "Negative (supports alternative model)"
  } else if (K >= 1 && K < 10^(1 / 2)) {
    evidence <- "Barely worth mentioning"
  } else if (K >= 10^(1 / 2) && K < 10^(1)) {
    evidence <- "Substantial"
  } else if (K >= 10^1 && K < 10^(3 / 2)) {
    evidence <- "Strong"
  } else if (K >= 10^(3 / 2) && K < 10^(2)) {
    evidence <- "Very strong"
  } else {
    evidence <- "Decisive"
  }
  evidence
}
