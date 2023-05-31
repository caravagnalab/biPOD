
#' Compute the pairwise Bayes Factors between models
#'
#' @param models A list of biPOD object of class `bipod` that have been fitted with the same prior.
#'
#' @returns A tibble with pairwise Bayes Factors. The column names are the sample names of the models.
#' @export
bf_selection = function(models) {
  check_model_selection_input(models = models)

  # Compute Marginal Likelihoods for every model
  marginal_likelihoods <- lapply(models, function(m) {
    b <- bridgesampling::bridge_sampler(m$fit, silent = TRUE)
    unname(b$logml)
  })

  marginal_likelihoods <- marginal_likelihoods %>% unlist()
  # Compute pairwise Bayes factor
  bayes_factors <- outer(marginal_likelihoods, marginal_likelihoods, FUN = function(x,y) {
    return(exp(x-y))
  })
  bayes_factors <- dplyr::as_tibble(bayes_factors)
  colnames(bayes_factors) <- sapply(models, function(m) {return(m$sample)})
  bayes_factors
}
