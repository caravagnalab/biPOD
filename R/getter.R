get_model <- function(model_name) {
  all_paths <- list(
    "exponential" = "exponential.stan",
    "exponential_no_t0" = "exponential_no_t0.stan",
    "logistic" = "logistic.stan",
    "logistic_no_t0" = "logistic_no_t0.stan",
    "infer_changepoints" = "infer_changepoints.stan",
    "exp_log_mixture" = "exp_log_mixture.stan",
    "two_pop" = "two_population.stan",
    "piecewise_changepoints" = "piecewise_linear_regression.stan"
  )

  if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")

  model_path <- system.file("cmdstan", all_paths[[model_name]], package = "biPOD", mustWork = T)

  tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
  model
}

get_parameter <- function(fit, par_name) {
  v <- fit$draws[,grepl(par_name, colnames(fit$draws), fixed = T)] %>% unlist() %>% unname()
  dplyr::tibble(value = v, parameter = par_name)

  # fit$draws[[par_name]]
  #   dplyr::as_tibble() %>%
  #   dplyr::mutate(par_name = par_name) %>%
  #   dplyr::rename("value" = 1, "parameter" = 2) %>%
  #   dplyr::mutate(value = as.numeric(.data$value))
}

get_parameters <- function(fit, par_list) {
  d <- dplyr::tibble()
  for (p in par_list) {
    d <- dplyr::bind_rows(d, get_parameter(fit, p))
  }
  d
}
