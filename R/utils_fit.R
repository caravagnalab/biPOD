prep_data_fit <- function(x, factor_size, prior_K, t0_lower_bound) {
  # Parameters check
  if (is.null(prior_K)) {
    prior_K <- max(x$counts$count) / factor_size
  } else {
    prior_K <- prior_K / factor_size
    if (prior_K <= 0) stop("'prior_K' should eiter be NULL or positive")
  }

  # Prepare data
  if (is.null(x$metadata$breakpoints)) {
    G <- 1
    breakpoints <- array(0, dim = c(0))
  } else {
    G <- length(x$metadata$breakpoints) + 1
    breakpoints <- x$metadata$breakpoints
  }

  # Prepare input data list
  input_data <- list(
    S = nrow(x$counts),
    G = G,
    N = as.array(as.integer(x$counts$count / factor_size)),
    T = as.array(x$counts$time),
    t_array = as.array(breakpoints),
    t0_lower_bound = t0_lower_bound,
    prior_K = prior_K
  )

  return(input_data)
}

# Fit exponential growth model to bipod object
fit_data <- function(x,
                     growth_type = "exponential",
                     factor_size = 1,
                     variational = FALSE,
                     t0_lower_bound = -10,
                     prior_K = NULL,
                     chains = 4, iter = 4000, cores = 4) {
  input_data <- prep_data_fit(x = x, factor_size = factor_size, prior_K = prior_K, t0_lower_bound = t0_lower_bound)

  # Get the model
  if (t0_lower_bound == x$counts$time[1]) {
    model_name <- paste0(growth_type, "_no_t0")
  } else {
    model_name <- growth_type
  }
  model <- get_model(model_name = model_name)

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
    growth_type = growth_type,
    factor_size = factor_size,
    t0_lower_bound = t0_lower_bound,
    prior_K = input_data$prior_K
  )

  res <- list(
    elbo_data = elbo_data,
    fit = fit,
    fit_info = fit_info
  )

  return(res)
}
