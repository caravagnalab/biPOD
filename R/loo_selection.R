
#' Compute the pairwise LOO between models
#'
#' @param models A list of biPOD object of class `bipod` that have been fitted with the same prior.
#'
#' @returns A tibble with pairwise LOOs. The column names are the sample names of the models.
#' @export
loo_selection = function(models) {
  check_model_selection_input(models = models)

  loos <- lapply(models, function(m) {
    log_lik <- loo::extract_log_lik(m$fit, merge_chains = F)
    r_eff <- loo::relative_eff(exp(log_lik), cores=4)
    loo <- loo::loo(log_lik, r_eff = r_eff, cores=4)
    return(loo)
  })

  loos_comparison <- loo::loo_compare(loos)
  rownames(loos_comparison) <- sapply(models, function(m) {return(m$sample)})
  loos_comparison
}
