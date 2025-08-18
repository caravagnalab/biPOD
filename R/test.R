
test = function() {

  library(tidyverse)
  all_data = readRDS("/Users/jovoni/Dropbox/Zenodo_biPOD_v2/Zenodo/method_validation/sauer/data/processed_data.rds")
  data = all_data %>%
    dplyr::filter(mouse == 544) %>%
    dplyr::filter(time <= 20, time >= -10)

  # Fit changepoints
  bp_fit = biPOD:::fit_breakpoints(data = data, with_initiation = F, max_segments = 3, min_segment_size = 2, comparison = "bic",
                           chains = 4, iter = 4000, cores = 4, t_prior_sd = 1,
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

  res = biPOD:::fit_growth(data = data, breakpoints = bps, with_initiation = T,
                           chains = 4, iter = 4000,
                           seed = 123, cores = 4, comparison = "bic", method = "vi", use_elbo = F)

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
                                              iter = 4000,
                                              seed = 123,
                                              cores = 4,
                                              comparison = c("bic"))

  x = u_results
  biPOD:::plot_growth_model_selection(u_results)
  biPOD::plot_u_ribbon(fit = u_results, data = data_u, ci = .5)

}

