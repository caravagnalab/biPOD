# Input checks run everywhere; model fits need CmdStan and run only when
# BIPOD_SLOW_TESTS=true (each takes about a minute).

slow <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
}

test_that("fit_held_phase rejects invalid input", {
  d <- data.frame(time = c(0, 10, 20, 30, 40, 50), count = c(100, 50, 20, 10, 10, 30))
  expect_error(fit_held_phase(d[, "time", drop = FALSE], t_decline_max = 20), "count")
  expect_error(fit_held_phase(transform(d, count = c(0, d$count[-1])), t_decline_max = 20), "positive")
  expect_error(fit_held_phase(d[1:4, ], t_decline_max = 20), "five")
  expect_error(fit_held_phase(d, t_decline_max = -1), "after the first")
  expect_error(fit_held_phase(rbind(d, d[6, ]), t_decline_max = 20), "unique")
})

test_that("fit_clone_selection rejects invalid input", {
  d <- data.frame(time = rep(c(0, 30, 60), 2), clone = rep(c("A", "B"), each = 3), count = c(10, 1, 0, 5, 1, 4))
  expect_error(fit_clone_selection(d[, c("time", "count")]), "clone")
  expect_error(fit_clone_selection(transform(d, count = c(-1, d$count[-1]))), "non-negative")
  expect_error(fit_clone_selection(d[d$clone == "A", ]), "two clones")
  expect_error(fit_clone_selection(d[-2, ]), "one row per time point")
})

test_that("compare_held_phase finds a simulated held phase and not a false one", {
  slow()
  set.seed(1)
  t <- c(0, 20, 40, 60, 80, 100, 130, 160, seq(200, 1000, by = 50))
  sim <- function(held_days) {
    t_n <- 120; log_T <- log(300) - 0.04 * (pmin(t, t_n)) + 0.02 * pmax(t - (t_n + held_days), 0)
    exp(log(8 + exp(log_T)) + 0.1 * stats::rt(length(t), 5))
  }
  with_hold <- compare_held_phase(data.frame(time = t, count = sim(600)), t_decline_max = 150, cores = 2)
  expect_equal(with_hold$evidence, "needed")
  h <- with_hold$held$parameters
  expect_true(h$q5[h$parameter == "held_days"] < 600 && h$q95[h$parameter == "held_days"] > 450)
  without <- compare_held_phase(data.frame(time = t, count = sim(0)), t_decline_max = 150, cores = 2)
  expect_equal(without$evidence, "not needed")
})

test_that("fit_clone_selection detects the selected clone", {
  slow()
  set.seed(2)
  times <- c(0, 30, 60, 120, 250, 400, 550)
  mu_a <- 300 * exp(-0.08 * pmin(times, 100))
  mu_b <- 30 * exp(-0.05 * pmin(times, 100) + 0.03 * pmax(times - 100, 0))
  d <- data.frame(time = rep(times, 2), clone = rep(c("A", "B"), each = length(times)),
                  count = c(stats::rnbinom(length(times), mu = mu_a, size = 20),
                            stats::rnbinom(length(times), mu = mu_b, size = 20)), depth = 1000)
  fit <- fit_clone_selection(d, cores = 2, iter = 1000)
  expect_gt(fit$clones$p_selected[fit$clones$clone == "B"], 0.95)
  expect_lt(fit$clones$p_selected[fit$clones$clone == "A"], 0.05)
})
