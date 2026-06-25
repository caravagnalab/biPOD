
get_model <- function(model_name, noise_model) {
  # Standard ODE models (poisson / lognormal / negbinomial)
  ode_models <- list(
    "exponential_no_init_bp"    = "exponential_no_init_bp.stan",
    "exponential_no_init"       = "exponential_no_init.stan",
    "exponential_with_init_bp"  = "exponential_with_init_bp.stan",
    "exponential_with_init"     = "exponential_with_init.stan",
    "gompertz_no_init"          = "gompertz_no_init.stan",
    "gompertz_with_init"        = "gompertz_with_init.stan",
    "logistic_no_init"          = "logistic_no_init.stan",
    "logistic_with_init"        = "logistic_with_init.stan",
    "monomolecular_with_init"   = "monomolecular_with_init.stan",
    "quadraticexp_with_init"    = "quadraticexp_with_init.stan",
    "monomolecular_no_init"     = "monomolecular_no_init.stan",
    "quadraticexp_no_init"      = "quadraticexp_no_init.stan",
    "two_pop_both"              = "two_pop_both.stan",
    "two_pop_denovo"            = "two_pop_denovo.stan",
    "two_pop_preexisting"       = "two_pop_preexisting.stan",
    "two_pop_single"            = "two_pop_single.stan"
  )

  # GBM models live under inst/cmdstan/gbm/{noise_model}/
  gbm_models <- list(
    "gbm_exponential"   = "exponential.stan",
    "gbm_logistic"      = "logistic.stan",
    "gbm_gompertz"      = "gompertz.stan",
    "gbm_monomolecular" = "monomolecular.stan",
    "gbm_quadraticexp"  = "quadraticexp.stan"
  )

  if (model_name %in% names(ode_models)) {
    model_file <- ode_models[[model_name]]
    model_path <- system.file("cmdstan", noise_model, model_file, package = "biPOD", mustWork = TRUE)
  } else if (model_name %in% names(gbm_models)) {
    model_file <- gbm_models[[model_name]]
    model_path <- system.file("cmdstan", "gbm", noise_model, model_file, package = "biPOD", mustWork = TRUE)
  } else {
    all_names <- c(names(ode_models), names(gbm_models))
    stop(sprintf("model_name '%s' not recognized. Available models: %s",
                 model_name, paste(all_names, collapse = ", ")))
  }

  if (!file.exists(model_path)) stop(paste0(model_path, " model does not exist"))

  suppressMessages(suppressWarnings(cmdstanr::cmdstan_model(model_path)))
}
