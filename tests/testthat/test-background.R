# Assay background in the recovery models: a cleared series with a small
# count after clearance must not be forced into a regrowth model.

slow <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(identical(Sys.getenv("BIPOD_SLOW_TESTS"), "true"), "set BIPOD_SLOW_TESTS=true to run")
}

cleared <- data.frame(time = c(0, 10, 20, 30, 50, 70, 90, 110, 130, 150, 170),
                      count = c(4000, 2000, 500, 100, 0, 0, 0, 0, 0, 0, 0))
with_final <- function(...) {
  v <- c(...)
  rbind(cleared, data.frame(time = 170 + 10 * seq_along(v), count = v))
}
regrowth <- c("pre-existing", "de-novo")

fit <- function(d, ...) {
  suppressMessages(fit_best_recovery_model(d, noise_model = "negbinomial", seed = 1, iter = 1000, ...))
}

test_that("a single count of one after clearance does not select regrowth", {
  slow()
  # without a background this selected de-novo by a BIC margin of ~17
  expect_equal(fit(with_final(1))$best_model, "single-pop-shrinking")
})

test_that("a known background keeps small final counts on the decay model", {
  slow()
  for (v in c(0, 1, 3, 5)) {
    expect_equal(fit(with_final(v), background = 1)$best_model, "single-pop-shrinking")
  }
  expect_equal(fit(with_final(3, 0, 0), background = 1)$best_model, "single-pop-shrinking")
})

test_that("a sustained rise after clearance is still regrowth", {
  slow()
  d <- with_final(5, 20, 80)
  expect_true(fit(d)$best_model %in% regrowth)
  expect_true(fit(d, background = 1)$best_model %in% regrowth)
})

test_that("a series without counts after clearance has a background near zero", {
  slow()
  res <- fit(cleared)
  expect_equal(res$best_model, "single-pop-shrinking")
  expect_lt(res$background, 0.5)
})

test_that("fixing the background at its fitted value gives the same fit", {
  slow()
  d <- with_final(1)
  sampled <- fit(d)
  fixed <- fit(d, background = sampled$background)
  expect_equal(fixed$best_model, sampled$best_model)
  rho <- function(r) r$final_fit$summary$median[r$final_fit$summary$variable == "rho_s"]
  expect_equal(rho(fixed), rho(sampled), tolerance = 0.1)
})
