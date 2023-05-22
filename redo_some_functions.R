

d <- biPOD::sim_stochastic_exponential(100, 1, c(0,0,0,1.3,1.3,1.3, .2, .2, .2), steps = 9, delta_t = .5)

x <- biPOD::init(d, "")

biPOD::plot_input(x, add_highlights = T)

x <- biPOD::breakpoints_inference(x, dt = .5, warmup = 5000, iter = 10000, variational = F)
biPOD::plot_breakpoints_posterior(x)

biPOD::plot_input(x, add_highlights = T)

biPOD::plot_elbo(x)


x <- biPOD::fit(x, t0_lower_bound = -10, model_selection = T)

biPOD::plot_fit(x, CI = .5)
biPOD::plot_simple_fit(x, CI = 0)

biPOD::plot_traces(x, x$fit, pars = c("rho[1]", "K", "rho[2]"), diagnose = T)


biPOD::plot_normalized_growth_rate_posteriors(x)

biPOD::get_growth_rate_posteriors(x)
