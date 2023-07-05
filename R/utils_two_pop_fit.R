prep_data_two_pop_fit <- function(x, factor_size) {
  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time)
  )
  return(input_data)
}

fit_two_pop_data <- function(x, factor_size, variational, chains, iter, cores) {
  input_data <- prep_data_two_pop_fit(x = x, factor_size = factor_size)
  model <- get_model(model_name = "two_pop")

  # Fit with either MCMC or Variational
  if (variational) {
    sampling <- "variational"
    res <- variational_fit(model = model, data = input_data, iter = iter)
    fit_model <- res$fit_model
    elbo_d <- res$elbo_d
  } else {
    sampling <- "mcmc"
    tmp <- utils::capture.output(suppressMessages(fit_model <- model$sample(data = input_data, chains = chains, parallel_chains = cores, iter_warmup = iter, iter_sampling = iter, refresh = iter)))
  }

  elbo_data <- c()
  if (variational) elbo_data <- elbo_d
  fit <- fit_model

  # Write fit info
  fit_info <- list(
    sampling = sampling,
    factor_size = factor_size
  )

  res <- list(
    elbo_data = elbo_data,
    fit = fit,
    fit_info = fit_info
  )

  return(res)
}
