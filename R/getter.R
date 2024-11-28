get_model <- function(model_name) {
  all_paths <- list(
    "exponential" = "exponential.stan",
    "exponential_no_t0" = "exponential_no_t0.stan",
    "logistic" = "logistic.stan",
    "logistic_no_t0" = "logistic_no_t0.stan",
    "exp_log_mixture" = "exp_log_mixture.stan",
    #"infer_changepoints" = "infer_changepoints.stan",
    "fit_breakpoints" = "fit_breakpoints.stan",
    "piecewise_changepoints" = "piecewise_linear_regression.stan",
    #"two_pop" = "two_population.stan",
    "two_pop_both" = 'two_pop_both_v2.stan',
    "two_pop_single" = 'two_pop_single.stan',
    "two_pop_pre" = 'two_pop_pre.stan',
    "two_pop_post" = 'two_pop_post.stan'
    #"pw_lin_fixed_b" = "pw_linear_b_fixed.stan",
  )

  if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")

  model_path <- system.file("cmdstan", all_paths[[model_name]], package = "biPOD", mustWork = T)

  tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
  model
}

get_parameter <- function(fit, par_name, variational) {
  if (variational) {
    v <- fit$draws[,which(par_name==fit$parameters)] %>% as.vector() %>% unlist() %>% unname()
  } else {
    v <- fit$draws[,,which(par_name==fit$parameters)] %>% as.vector() %>% unlist() %>% unname()
  }
  #v <- fit$draws[,grepl(par_name, colnames(fit$draws), fixed = T)] %>% unlist() %>% unname()

  dplyr::tibble(value = v, parameter = par_name)

  # fit$draws[[par_name]]
  #   dplyr::as_tibble() %>%
  #   dplyr::mutate(par_name = par_name) %>%
  #   dplyr::rename("value" = 1, "parameter" = 2) %>%
  #   dplyr::mutate(value = as.numeric(.data$value))
}

get_parameters <- function(fit, par_list, variational) {
  d <- dplyr::tibble()
  for (p in par_list) {
    d <- dplyr::bind_rows(d, get_parameter(fit, p, variational = variational))
  }
  d
}
