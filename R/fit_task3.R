#' Fit a Two-Population Growth Model to a bipod Object
#'
#' This function fits a growth model that accounts for two populations within a bipod object. The fitting can be done using either Variational Inference or Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @param x A `bipod` object.
#' @param variational A logical value indicating whether to use Variational Inference instead of MCMC sampling.
#'  If `TRUE`, the model will be fitted using Variational Inference; otherwise, MCMC sampling will be used. (default is FALSE)
#' @param factor_size A numeric value representing the factor by which to divide the counts in the bipod object.
#'  This value must be positive and appropriate for the data scale. (default is 1)
#' @param chains An integer specifying the number of chains to run in the MCMC algorithm.
#'  This parameter is ignored if `variational = TRUE`. (default is 4)
#' @param iter An integer specifying the number of iterations to run in the MCMC algorithm.
#'  This parameter is ignored if `variational = TRUE`. (default is 5000)
#' @param cores An integer specifying the number of cores to use for parallel processing during model fitting. (default is 4)
#'
#' @return The input `bipod` object with added slots:
#' - `'two_pop_fit'`: Contains the fitted two-population growth model.
#' - `'two_pop_fit_info'`: Contains information about the fitting process, including metadata such as sampling type and factor size.
#'
#' @export
fit_two_pop_model <- function(
    x,
    variational = FALSE,
    factor_size = 1,
    chains = 4,
    iter = 5000,
    cores = 4) {
  # Check input
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  if (!(factor_size > 0)) stop("factor_size must be positive")
  sampling_type <- if (variational) "variational inference" else "MCMC sampling"

  cli::cli_alert_info(paste("Fitting two population model using", sampling_type, "..."))
  cat("\n")

  res <- fit_two_pop_data(
    x = x,
    factor_size = factor_size,
    variational = variational,
    chains = chains,
    iter = iter,
    cores = cores
  )

  # Add results to bipod object
  x$two_pop_fit_elbo <- res$elbo_data
  x$metadata$chains = chains

  if (sampling_type == "MCMC sampling") {
    x$metadata$status <- diagnose_mcmc_fit(res$fit)
  } else {
    x$metadata$status <- diagnose_variational_fit(res$fit, res$elbo_data)
  }

  # Add metadata
  x$metadata$sampling <- res$fit_info$sampling
  x$metadata$factor_size <- res$fit_info$factor_size

  # Produce plots ####
  ## Produce evo plot ####
  best_fit <- res$fit
  draws <- best_fit$draws(format = "df", variables = "yrep")

  mu <- draws %>%
    dplyr::as_tibble() %>%
    dplyr::select(!c(.data$.chain, .data$.iteration, .data$.draw)) %>%
    #select(starts_with("mu")) %>%
    apply(2, stats::quantile, c(00.05, 0.5, 0.95)) %>%
    t() %>%
    data.frame(x = x$counts$time) %>%
    tidyr::gather(pct, y, -x)
  mu$y <- mu$y * factor_size

  evo_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=.data$time, y=.data$count), data = x$counts, size = 1) +
    ggplot2::geom_line(ggplot2::aes(x=.data$x, y=.data$y, linetype = .data$pct), data = mu, color = 'darkgreen') +
    ggplot2::scale_linetype_manual(values = c(2,1,2)) +
    ggplot2::guides(linetype = "none") +
    ggplot2::theme_bw()

  ## Produce evo separated plot ####
  draws <- best_fit$draws(format = "df", variables = "ns")
  ns <- draws %>%
    dplyr::as_tibble() %>%
    dplyr::select(!c(.data$.chain, .data$.iteration, .data$.draw)) %>%
    apply(2, stats::quantile, c(00.05, 0.5, 0.95)) %>%
    t() %>%
    data.frame(x = x$counts$time) %>%
    tidyr::gather(pct, y, -x) %>%
    dplyr::mutate(y = .data$y * factor_size)

  draws <- best_fit$draws(format = "df", variables = "nr")
  nr <- draws %>%
    dplyr::as_tibble() %>%
    dplyr::select(!c(.data$.chain, .data$.iteration, .data$.draw)) %>%
    apply(2, stats::quantile, c(00.05, 0.5, 0.95)) %>%
    t() %>%
    data.frame(x = x$counts$time) %>%
    tidyr::gather(pct, y, -x) %>%
    dplyr::mutate(y = .data$y * factor_size)

  evo_plot_separated <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=.data$time, y=.data$count), data = x$counts, size = 1) +
    ggplot2::geom_line(ggplot2::aes(x=.data$x, y=.data$y, linetype = .data$pct), data = nr, color = 'steelblue3') +
    ggplot2::geom_line(ggplot2::aes(x=.data$x, y=.data$y, linetype = .data$pct), data = ns, color = 'indianred') +
    ggplot2::scale_linetype_manual(values = c(2,1,2)) +
    ggplot2::guides(linetype = "none") +
    ggplot2::theme_bw()

  ## Produce times plot ####
  variables <- best_fit$summary()$variable

  times <- dplyr::tibble(x = best_fit$draws(format = "df", variables = "t0_r")[,1] %>% unlist() %>% as.numeric(), par = 't0_r')
  t0_r_draws <- times$x
  if ("t_end" %in% variables) {
    times <- dplyr::bind_rows(
      times,
      dplyr::tibble(x = best_fit$draws(format = "df", variables = "t_end")[,1] %>% unlist() %>% as.numeric(), par = 't_end')
    )
  }

  times_plot <- times %>%
    dplyr::group_by(.data$par) %>%
    dplyr::mutate(q_low = stats::quantile(x, .05), q_high = stats::quantile(x, .95)) %>%
    dplyr::filter(x >= .data$q_low, x <= .data$q_high) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$x, col=.data$par)) +
    ggplot2::geom_density() +
    ggplot2::labs(x = "Time", y="Density") +
    ggplot2::geom_vline(xintercept = x$counts$time[1], linetype = "dashed") +
    ggplot2::theme_bw()

  ## Produce plot of F_R ####
  nr_first_obs <- best_fit$draws(format = "df", variables = "nr")[,1] %>% unlist() %>% as.numeric()
  ns_first_obs <- best_fit$draws(format = "df", variables = "ns")[,1] %>% unlist() %>% as.numeric()

  nr_first_obs <- nr_first_obs[t0_r_draws <= 0]
  ns_first_obs <- ns_first_obs[t0_r_draws <= 0]

  f_r = nr_first_obs / (nr_first_obs + ns_first_obs)
  fr_plot <- dplyr::tibble(x = f_r) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=x)) +
    ggplot2::geom_density() +
    ggplot2::labs(x = "Resistant population fraction", y="Density") +
    ggplot2::lims(x = c(0,1)) +
    ggplot2::theme_bw()

  ## Rates plot ####
  rho_r_draws <- best_fit$draws(format = "df", variables = "rho_r")[,1] %>% unlist() %>% as.numeric()
  rates <- dplyr::tibble(x=best_fit$draws(format = "df", variables = "rho_r")[,1] %>% unlist() %>% as.numeric(), par = "rho_r")
  if ("rho_s" %in% variables) {
    rates <- dplyr::bind_rows(
      rates,
      dplyr::tibble(x=-(best_fit$draws(format = "df", variables = "rho_s")[,1] %>% unlist() %>% as.numeric()), par = "rho_s")
    )
  }

  rates_plot <- rates %>%
    dplyr::group_by(.data$par) %>%
    dplyr::mutate(q_low = stats::quantile(x, .05), q_high = stats::quantile(x, .95)) %>%
    dplyr::filter(x >= .data$q_low, x <= .data$q_high) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$x, col=.data$par)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "Growth rate", y="Density") +
    ggplot2::theme_bw()

  x$two_pop_plots <- list(
    evo_plot = evo_plot,
    evo_plot_separated = evo_plot_separated,
    rates_plot = rates_plot,
    times_plot = times_plot,
    fr_plot = fr_plot
  )

  x$two_pop_fit <- convert_mcmc_fit_to_biPOD(res$fit, variational=variational)
  #x$two_pop_fit <- res$fit
  return(x)
}

# Utils fitting taks 3

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
  #model <- get_model(model_name = "two_pop")

  m1 <- get_model("two_pop_single")
  m2 <- get_model("two_pop_both")

  tmp <- utils::capture.output(f1 <- m1$sample(input_data, parallel_chains = chains, iter_warmup = iter, iter_sampling = iter, chains = chains))
  tmp <- utils::capture.output(f2 <- m2$sample(input_data, parallel_chains = chains, iter_warmup = iter, iter_sampling = iter, chains = chains))

  fits <- list(f1, f2)

  loo1 <- f1$loo()
  loo2 <- f2$loo()

  loos <- loo::loo_compare(list(loo1, loo2))
  fit_model <- fits[[as.numeric(stringr::str_replace(rownames(loos)[1], "model", ""))]]

  sampling <- "mcmc"
  elbo_d <- NULL
  # Fit with either MCMC or Variational
  # if (variational) {
  #   sampling <- "variational"
  #   res <- variational_fit(model = model, data = input_data, iter = iter)
  #   fit_model <- res$fit_model
  #   elbo_d <- res$elbo_d
  # } else {
  #   sampling <- "mcmc"
  #   tmp <- utils::capture.output(
  #     suppressMessages(
  #       fit_model <- model$sample(data = input_data, chains = chains, parallel_chains = cores, iter_warmup = iter, iter_sampling = iter, refresh = iter)
  #     )
  #   )
  # }

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
