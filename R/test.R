
test = function() {

  # Monomolecular
  monomolecular_f = function(xs, K, n0, alpha) {
    K - ((K - n0) * exp(- alpha * (xs - 0)))
  }

  quadexp_f = function(xs, n0, rho, beta) {
    n0 * exp(rho * (xs - 0) + beta * (xs - 0)**2)
  }

  xs = c(0:10)
  ys = monomolecular_f(xs, 1000, 10, .25)
  #ys = quadexp_f(xs, 10, .25, .02)
  ys = rpois(length(ys), ys)

  d = dplyr::tibble(time = xs, count = ys)

  d$count = as.integer(d$count)
  x = biPOD:::fit_growth(data = d,
                         with_initiation = T,
                         noise_model = "lognormal",
                         models_to_fit = c("monomolecular", "exponential", "gompertz", "logistic", "quadraticexp"),
                         chains = 4,
                         iter = 4000,
                         seed = 1234,
                         comparison = "loo",
                         cores = 4,
                         method = "sampling",
                         use_elbo = F,
                         breakpoints = NULL)

  biPOD:::plot_ribbon(x$fit, d) +
    ggplot2::scale_y_continuous(transform = "log10")


  biPOD::plot_growth_model_selection(x)
  biPOD:::plot_growth_fit(x=x, data = d, CI = 0.95)

  library(tidyverse)
  all_data = readRDS("/Users/jovoni/Dropbox/Zenodo_biPOD_v2/Zenodo/method_validation/sauer/data/processed_data.rds")
  data = all_data %>%
    dplyr::filter(mouse == 544) %>%
    dplyr::filter(time <= 20, time >= -10)

  # Fit changepoints
  bp_fit = biPOD:::fit_breakpoints(data = data,
                                   with_initiation = F,
                                   max_segments = 3,
                                   min_segment_size = 2,
                                   comparison = "bic",
                                   chains = 4, iter = 4000, cores = 4, t_prior_sd = 0.5, noise_model = "poisson",
                                   enforce_rho_separation = T, alpha_rho = .05, models_to_fit = c("exponential"),
                                   seed = 1234)


  biPOD:::plot_breakpoint_model_selection(bp_fit)

  p = biPOD:::plot_ribbon(bp_fit$final_fit, data = data, ci = .9, shadow_breakpoints = bp_fit$final_breakpoints)
  biPOD:::plot_breakpoint_posterior(bp_fit = bp_fit, data = data, colors = NULL)



  # Extract best breakpoints
  bps = bp_fit$final_breakpoints

  res = biPOD:::fit_growth(data = data, breakpoints = bps, with_initiation = T,
            chains = 4, iter = 4000,
            seed = 123, cores = 4, comparison = "bic", method = "sampling", use_elbo = F)

  biPOD:::plot_growth_model_selection(res)


  biPOD::plot_growth_fit(x = res, data = data, color = "black")

  biPOD:::plot_ribbon(fit = res$fit, data = data, ci = .5, shadow_breakpoints = bps) +
    geom_vline(xintercept = bps, linetype = "dashed")

  biPOD:::plot_breakpoint_posterior(bp_fit = bp_fit, data = data, colors = NULL)


  biPOD:::plot_parameter_posteriors(res$fit$draws, params = c("rho[1]"), faceted = F)
  biPOD:::plot_parameter_posteriors(fit_draws = res$fit$draws, params = c("rho\\["), faceted = F)

  # Fit U-shape
  data_u = data %>% dplyr::filter(time >= 0)
  u_results = biPOD:::fit_best_recovery_model(data = data_u,
                                              chains = 4,
                                              noise_model = "poisson",
                                              iter = 4000,
                                              seed = 123,
                                              cores = 4,
                                              comparison = c("bic"))

  x = u_results
  #biPOD:::plot_growth_model_selection(u_results)
  biPOD::plot_u_ribbon(fit = u_results, data = data_u, ci = .5)

}

