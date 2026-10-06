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
