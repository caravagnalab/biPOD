# BIC of the growth models: best draw's log-likelihood, every sampled parameter

test_that("growth-model BIC uses the best draw and counts every sampled parameter", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
  d <- data.frame(time = c(0, 20, 40, 60, 80, 100, 120, 140), count = c(800, 300, 90, 40, 60, 150, 400, 900))
  res <- fit_growth_models(d, breakpoints = 70, with_initiation = FALSE, chains = 2, iter = 300,
                           cores = 2, models_to_fit = "exponential", noise_model = "negbinomial",
                           comparison = "bic")
  fit <- res$fits$exponential
  ll <- posterior::as_draws_matrix(fit$draws("log_lik"))
  # n0, rho[1], rho[2], phi
  expect_equal(res$model_table$BIC, -2 * max(rowSums(ll)) + log(nrow(d)) * 4)
})
