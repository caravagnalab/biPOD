fit_with_mixture_model <- function(x,
                                   factor_size = 1,
                                   variational = FALSE,
                                   t0_lower_bound = -10,
                                   prior_K = NULL,
                                   chains = 4, iter = 4000, cores = 4) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size, prior_K = prior_K, t0_lower_bound = t0_lower_bound)

  model <- get_model(model_name = "exp_log_mixture")

  tmp <- utils::capture.output(suppressMessages(fit <- model$sample(input_data, chains = chains, iter_warmup = iter, iter_sampling = iter, parallel_chains = cores)))

  # extract omega draws and compute probability of log and of exp
  # omega >= 0.5 suggests Exponential, and vice versa
  n_exponential <- (fit$draws("omega", format = "matrix") >= .5) %>% sum()
  n_logistic <- (fit$draws("omega", format = "matrix") < .5) %>% sum()

  p_exp <- n_exponential / (n_exponential + n_logistic)
  p_log <- 1 - p_exp

  if (p_exp >= p_log) {
    best_growth <- "Exponential"
    odds <- p_exp / p_log
    res <- fit_data(
      x = x, growth_type = "exponential", factor_size = factor_size,
      variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
      chains = chains, iter = iter, cores = cores
    )
  } else {
    best_growth <- "Logistic"
    odds <- p_log / p_exp
    res <- fit_data(
      x = x, growth_type = "logistic", factor_size = factor_size,
      variational = variational, t0_lower_bound = t0_lower_bound, prior_K = prior_K,
      chains = chains, iter = iter, cores = cores
    )
  }

  res$fit_info$best_growth <- best_growth
  res$fit_info$odds <- odds
  res$fit_info$model_selection_algo <- "mixture_model"
  res$fit_info$omega_mixture_model <- fit$draws("omega", format = "matrix")

  return(res)
}
