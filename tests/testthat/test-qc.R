# stan_qc() must survive a non-finite log-likelihood: model choice by BIC
# does not need loo.

slow <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
}

zero_heavy <- data.frame(
  time = c(0, 109, 178, 253, 262, 351, 369, 513, 600, 691, 1160, 1189, 1458, 1483),
  count = c(129569, 538, 25, 0, 0, 0, 0, 0, 0, 0, 869, 558, 18622, 46952)
)

growth_stan_data <- function(d) {
  list(S = nrow(d), G = 1L, N = as.integer(d$count), T = d$time,
       t_array = array(0, dim = 0), prior_only = 0)
}

# A fit whose log_lik has one NaN, standing in for any non-finite likelihood
with_nan_loglik <- function(fit) {
  list(
    draws = function(variables = NULL, ..., format = "draws_array") {
      d <- fit$draws(variables, ..., format = format)
      if ("log_lik[1]" %in% posterior::variables(d)) {
        if (posterior::is_draws_array(d)) d[1, 1, "log_lik[1]"] <- NaN else d[1, "log_lik[1]"] <- NaN
      }
      d
    },
    sampler_diagnostics = function(...) fit$sampler_diagnostics(...),
    metadata = function(...) fit$metadata(...)
  )
}

test_that("stan_qc reports a non-finite log_lik instead of stopping in loo", {
  slow()
  stan_data <- growth_stan_data(zero_heavy)
  model <- biPOD:::get_model("exponential_no_init", "negbinomial")
  fit <- model$sample(data = stan_data, chains = 2, iter_warmup = 200, iter_sampling = 200,
                      seed = 1, refresh = 0, show_messages = FALSE)
  qc <- biPOD:::stan_qc(model, with_nan_loglik(fit), stan_data)
  expect_null(qc$loo)
  expect_equal(qc$n_nonfinite_loglik, 1)
  expect_match(qc$loo_error, "non-finite")
})

test_that("fit_growth_models with BIC survives a non-finite log_lik", {
  slow()
  model <- biPOD:::get_model("exponential_no_init", "negbinomial")
  nan_model <- list(sample = function(...) with_nan_loglik(model$sample(...)))
  local_mocked_bindings(get_model = function(...) nan_model)
  res <- fit_growth_models(zero_heavy, breakpoints = NULL, with_initiation = FALSE, chains = 2,
                           iter = 200, cores = 2, models_to_fit = "exponential",
                           noise_model = "negbinomial", comparison = "bic")
  expect_null(res$fits_qc$exponential$loo)
  expect_equal(res$fits_qc$exponential$n_nonfinite_loglik, 1)
})
