# Sparse, zero-heavy series on long windows: the combination that produced
# non-finite log-likelihoods and crashed fit_breakpoints().

slow <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
}

sparse <- data.frame(time = c(0, 60, 120, 180, 260, 340, 420), count = c(1500, 0, 0, 0, 2, 1, 40))

# 50% zeros, 1483-day window, counts up to ~1.3e5
zero_heavy <- data.frame(
  time = c(0, 109, 178, 253, 262, 351, 369, 513, 600, 691, 1160, 1189, 1458, 1483),
  count = c(129569, 538, 25, 0, 0, 0, 0, 0, 0, 0, 869, 558, 18622, 46952)
)

growth_stan_data <- function(d) {
  list(S = nrow(d), G = 1L, N = as.integer(d$count), T = d$time,
       t_array = array(0, dim = 0), prior_only = 0)
}

test_that("recovery models start from a rate whose exp() would overflow", {
  slow()
  stan_data <- list(S = nrow(sparse), N = as.integer(sparse$count), T = sparse$time, prior_only = 0)
  # rho_s * (t_end - T[1]) = 2000 > 709.78. Only the negative binomial can
  # start here: the Poisson log-likelihood of such a mean is -inf in any
  # parameterisation.
  init <- list(rho_s = 5, rho_r = 0.02, t0_r = 300, t_end = 400, phi = 10)
  fit <- biPOD:::get_model("two_pop_both", "negbinomial")$sample(
    data = stan_data, chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 1,
    init = list(init, init), refresh = 0, show_messages = FALSE
  )
  expect_equal(fit$num_chains(), 2)
  expect_true(all(is.finite(posterior::as_draws_matrix(fit$draws("log_lik")))))

  for (noise in c("poisson", "negbinomial")) {
    rec <- fit_best_recovery_model(sparse, noise_model = noise, iter = 500, cores = 4)
    expect_true(all(is.finite(rec$first_fit$summary$median)))
  }
})

test_that("growth models keep log_lik finite when the predicted mean is huge", {
  slow()
  # with this seed, neg_binomial_2_rng / poisson_rng in generated quantities
  # used to throw on means above 2^30, which set every draw's log_lik to NaN
  for (noise in c("poisson", "negbinomial")) {
    fit <- suppressWarnings(biPOD:::get_model("exponential_with_init", noise)$sample(
      data = growth_stan_data(zero_heavy), chains = 4, iter_warmup = 200, iter_sampling = 200,
      seed = 1234, refresh = 0, show_messages = FALSE
    ))
    expect_true(all(is.finite(posterior::as_draws_matrix(fit$draws("log_lik")))))
  }
})

test_that("series with 40%+ zeros fit without filtering", {
  slow()
  expect_gte(mean(zero_heavy$count == 0), 0.4)
  for (noise in c("poisson", "negbinomial")) {
    rec <- fit_best_recovery_model(zero_heavy, noise_model = noise, iter = 500, cores = 4)
    expect_true(rec$best_model %in% c("pre-existing", "de-novo", "single-pop-growing", "single-pop-shrinking"))
  }
  x <- suppressMessages(init(transform(zero_heavy, id = "z"), sample = "z"))
  x <- fit_breakpoints(x, with_initiation = FALSE, noise_model = "negbinomial", iter = 1000, cores = 4)
  expect_true(all(is.finite(x$breakpoints_fit$final_breakpoints)))
})
