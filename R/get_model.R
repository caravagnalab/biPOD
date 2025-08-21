
get_model <- function(model_name) {
  # Define available models and their file names
  all_paths <- list(
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

  # Check if model_name is valid
  if (!(model_name %in% names(all_paths))) {
    stop(sprintf("model_name '%s' not recognized. Available models: %s",
                 model_name, paste(names(all_paths), collapse = ", ")))
  }

  # Get the model file path from package
  model_file <- all_paths[[model_name]]
  model_path <- system.file("cmdstan", model_file, package = "biPOD", mustWork = TRUE)

  # Compile model quietly
  suppressMessages({
    suppressWarnings({
      model <- cmdstanr::cmdstan_model(model_path)
    })
  })

  return(model)
}
