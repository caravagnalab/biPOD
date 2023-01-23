
#' Compute the pairwise LOO between models
#'
#' @param models A list of biPOD object of class `bipod` that have been fitted with the same prior.
#'
#' @returns A tibble with pairwise LOOs. The column names are the sample names of the models.
#' @export
loo_selection = function(models) {
  check_model_selection_input(models = models)

  # Compute Marginal LOO for every model
  loos <- rep(0, length(models))
  for (i in 1:length(models)) {
    m = models[[i]]

    loo <- sapply(m$fits, function(fit) {
      log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
      r_eff <- loo::relative_eff(exp(log_lik), cores = 4)
      loo <- loo::loo(log_lik, r_eff = r_eff, cores = 4)
      return(loo$estimate)
    })

    loos[i] <- sum(loo[1,])
  }

  # Compute Bayes factor for every couple of model
  loos_differences <- outer(loos, loos, FUN = function(x,y) {
    return(x-y)
  })

  loos_differences <- dplyr::as_tibble(loos_differences)
  colnames(loos_differences) <- sapply(models, function(m) {return(m$sample)})

}

#' Compute the pairwise Bayes Factors between models
#'
#' @param models A list of biPOD object of class `bipod` that have been fitted with the same prior.
#'
#' @returns A tibble with pairwise Bayes Factors. The column names are the sample names of the models.
#' @export
bf_selection = function(models) {
  check_model_selection_input(models = models)

  # Compute Marginal Likelihoods for every model
  marginal_likelihoods <- rep(0, length(models))
  for (i in 1:length(models)) {
    m = models[[i]]

    bridge <- sapply(m$fits, function(fit) {
      b <- bridgesampling::bridge_sampler(fit, silent = TRUE)
      unname(b$logml)
    })

    marginal_likelihoods[i] <- sum(bridge)
  }

  # Compute pairwise Bayes factor
  bayes_factors <- outer(marginal_likelihoods, marginal_likelihoods, FUN = function(x,y) {
    return(exp(x-y))
  })
  bayes_factors <- dplyr::as_tibble(bayes_factors)
  colnames(bayes_factors) <- sapply(models, function(m) {return(m$sample)})
  bayes_factors
}

check_model_selection_input = function(models) {
  # Check similarities in Priors
  for (i in 1:length(models)) {
    if (!(inherits(models[[i]], "bipod"))) stop("Models must be a vector of bipod objects")
  }

  for (i in 2:length(models)) {
    info1 = models[[i-1]]$fit_info
    info2 = models[[i]]$fit_info

    if (!(info1$sampling == info2$sampling)) stop("Models must have the samesampling method!")
    if (!(info1$model_type == info2$model_type)) stop("Models must have the same model_type!")
    if (!(info1$prior == info2$prior)) stop("Models must have the same type of prior!")
    if (!(info1$a == info2$a)) stop("Models must have the same parameter 'a' for the prior!")
    if (!(info1$b == info2$b)) stop("Models must have the same parameter 'b' for the prior!")
    if (!(info1$g == info2$g)) stop("Models must have the same parameter 'g' for the prior!")
  }
}
