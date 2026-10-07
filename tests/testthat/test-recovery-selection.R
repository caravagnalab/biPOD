# Selection and convergence in fit_best_recovery_model(), with sampling
# replaced by fixture fits of known Rhat and log-likelihood

d <- data.frame(time = c(0, 50, 100, 150, 200), count = c(500, 20, 3, 40, 900))

# A stand-in for a CmdStanMCMC fit: every draw has log-likelihood ll
fixture_fit <- function(rhat, ll) {
  ll_mat <- posterior::as_draws_matrix(matrix(ll / nrow(d), nrow = 10, ncol = nrow(d),
                                              dimnames = list(NULL, sprintf("log_lik[%d]", seq_len(nrow(d))))))
  list(
    summary = function(variables = NULL, ...) {
      dplyr::tibble(variable = if (is.null(variables)) "rho_r" else variables, median = 0, rhat = rhat)
    },
    draws = function(variables = NULL, format = "draws_matrix", ...) {
      posterior::as_draws(ll_mat) |> posterior::as_draws_list()
    }
  )
}

# fixtures[[model]] gives the Rhat of successive fits of that model
mock_sampling <- function(fixtures, ll) {
  calls <- list()
  fake <- function(model_name, stan_data, data, noise_model, chains, iter, seed, cores) {
    calls[[length(calls) + 1]] <<- list(model = model_name, iter = iter, seed = seed)
    k <- sum(vapply(calls, function(x) x$model == model_name, logical(1)))
    fixture_fit(fixtures[[model_name]][k], ll[[model_name]])
  }
  list(fake = fake, calls = function() calls)
}

# two_pop_single has the best BIC
ll <- list(two_pop_both = -40, two_pop_single = -20, single_pop_decay = -60)

test_that("a non-converged best candidate is refitted, still selected, and flagged", {
  m <- mock_sampling(list(two_pop_both = 1.00, two_pop_single = c(1.50, 1.40), single_pop_decay = 1.00), ll)
  local_mocked_bindings(sample_recovery_model = m$fake)
  res <- suppressMessages(fit_best_recovery_model(d, iter = 100))

  expect_equal(res$best_model, "single-pop-growing")
  expect_setequal(res$model_table$model, c("two_pop_both", "two_pop_single", "single_pop_decay"))
  expect_equal(res$model_table$refit, c(FALSE, TRUE, FALSE))
  expect_equal(res$model_table$rhat[res$model_table$model == "two_pop_single"], 1.40)
  expect_equal(res$winner_rhat, 1.40)
  expect_equal(res$n_candidates_converged, 2)
  expect_false(res$converged)

  refit <- Filter(function(x) x$model == "two_pop_single", m$calls())
  expect_length(refit, 2)
  expect_equal(refit[[2]]$iter, 200)
  expect_false(refit[[2]]$seed == refit[[1]]$seed)
})

test_that("a refit that converges replaces the first fit", {
  m <- mock_sampling(list(two_pop_both = 1.00, two_pop_single = c(1.50, 1.01), single_pop_decay = 1.00), ll)
  local_mocked_bindings(sample_recovery_model = m$fake)
  res <- suppressMessages(fit_best_recovery_model(d, iter = 100))

  expect_equal(res$winner_rhat, 1.01)
  expect_equal(res$n_candidates_converged, 3)
  expect_true(res$converged)
})

test_that("a converged candidate is not refitted", {
  m <- mock_sampling(list(two_pop_both = 1.00, two_pop_single = 1.00, single_pop_decay = 1.00), ll)
  local_mocked_bindings(sample_recovery_model = m$fake)
  res <- suppressMessages(fit_best_recovery_model(d, iter = 100))

  expect_length(m$calls(), 3)
  expect_false(any(res$model_table$refit))
})

test_that("the rate priors of the recovery models follow rho_r_prior_sd and rho_s_prior_sd", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
  stan_data <- list(S = nrow(d), N = as.integer(d$count), T = d$time, prior_only = 1,
                    rho_r_prior_sd = 0.1, rho_s_prior_sd = 0.5,
                    fit_background = 1L, background_fixed = 0, background_prior_mean = 1)
  fit <- biPOD:::get_model("two_pop_both", "negbinomial")$sample(
    data = stan_data, chains = 2, iter_warmup = 1000, iter_sampling = 2000, seed = 1,
    refresh = 0, show_messages = FALSE
  )
  med <- apply(posterior::as_draws_matrix(fit$draws(c("rho_r", "rho_s"))), 2, stats::median)
  # median of a half-normal is 0.674 sd
  expect_equal(unname(med), c(0.0674, 0.337), tolerance = 0.1)
})

test_that("a sampled background is counted in the BIC of every candidate", {
  converged <- list(two_pop_both = 1.00, two_pop_single = 1.00, single_pop_decay = 1.00)
  local_mocked_bindings(sample_recovery_model = mock_sampling(converged, ll)$fake)
  sampled <- suppressMessages(fit_best_recovery_model(d, iter = 100))
  local_mocked_bindings(sample_recovery_model = mock_sampling(converged, ll)$fake)
  fixed <- suppressMessages(fit_best_recovery_model(d, iter = 100, background = 2))

  expect_equal(sampled$model_table$BIC - fixed$model_table$BIC, rep(log(nrow(d)), 3))
  expect_equal(sampled$model_table$BIC[sampled$model_table$model == "two_pop_single"],
               -2 * ll$two_pop_single + 3 * log(nrow(d)))
  expect_true("background" %in% names(sampled))
  expect_true("background" %in% names(sampled$model_table))
})

test_that("the background setting reaches the Stan data and the inits", {
  seen <- NULL
  fake <- function(model_name, stan_data, data, noise_model, chains, iter, seed, cores) {
    seen <<- stan_data
    fixture_fit(1.00, ll[[model_name]])
  }
  local_mocked_bindings(sample_recovery_model = fake)

  suppressMessages(fit_best_recovery_model(d, iter = 100, background_prior_mean = 3))
  expect_equal(seen$fit_background, 1L)
  expect_equal(seen$background_prior_mean, 3)

  suppressMessages(fit_best_recovery_model(d, iter = 100, background = 0.5))
  expect_equal(seen$fit_background, 0L)
  expect_equal(seen$background_fixed, 0.5)

  expect_error(fit_best_recovery_model(d, background = 0))
  expect_error(fit_best_recovery_model(d, background = c(1, 2)))
  expect_error(fit_best_recovery_model(d, background_prior_mean = 0))

  expect_equal(dim(recovery_inits(d, "single_pop_decay", "poisson", 1, background_prior_mean = 2)$background_par), 1)
  expect_null(recovery_inits(d, "single_pop_decay", "poisson", 1)$background_par)
})
