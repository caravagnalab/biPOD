#' Overlay Posterior Predictive Density
#'
#' Creates a density plot overlaying the posterior predictive density over the observed data.
#' This plot allows you to compare the distribution of observed data with the distribution predicted by the model, using posterior predictive draws.
#' It is useful for evaluating how well the model captures the data distribution.
#'
#' @param x A `bipod` object containing a 'counts' field with the observed data. This data is typically a time series or count data.
#' @param x_fit A fitted model object that contains posterior draws of the replicated data (`yrep`). This model should have been fitted to the data.
#' @param n An integer specifying the number of posterior predictive draws to use in the plot. (default is 500). The more draws used, the smoother the density overlay will be.
#'
#' @return A density plot generated using `bayesplot::ppc_dens_overlay`. The plot overlays the posterior predictive distribution with the observed data.
#' @export
biPOD_ppc_dens_overlay <- function(x, x_fit, n=500) {
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  y <- x$counts$count
  if (x$metadata$sampling == 'variational') {
    yrep <- x_fit$draws[,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  } else {
    yrep <- x_fit$draws[,,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  }
  yrep <- yrep[1:n,]
  if ("q" %in% x_fit$parameters) y = log(y)
  if (!(is.null(x_fit$normalization_pars))) y = (y - mean(y)) / stats::sd(y)
  bayesplot::ppc_dens_overlay(y, yrep)
}

#' Posterior Predictive Intervals Plot
#'
#' Creates a plot displaying posterior predictive intervals for the observed data, providing a visual assessment of how well the model captures data variation over time.
#' The plot shows the credible intervals (uncertainty) of the posterior predictive distribution at each time point.
#'
#' @param x A `bipod` object containing a 'counts' field with the observed data. This data typically represents counts or time series values.
#' @param x_fit A fitted model object containing posterior draws of the replicated data (`yrep`). This object should include posterior samples of the model predictions.
#' @param n An integer specifying the number of posterior predictive draws to use in the plot. (default is 500). The higher the number, the more accurate the representation of the predictive intervals.
#'
#' @return A plot displaying the posterior predictive intervals generated using `bayesplot::ppc_intervals`. The plot shows how well the model accounts for the observed data.
#' @export
biPOD_ppc_intervals <- function(x, x_fit, n=500) {
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  y <- x$counts$count
  if (x$metadata$sampling == 'variational') {
    yrep <- x_fit$draws[,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  } else {
    yrep <- x_fit$draws[,,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  }
  yrep <- yrep[1:n,]
  x <- x$counts$time
  if ("q" %in% x_fit$parameters) y = log(y)
  if (!(is.null(x_fit$normalization_pars))) y = (y - mean(y)) / stats::sd(y)
  bayesplot::ppc_intervals(y, yrep, x=x)
}

#' Posterior Predictive Ribbon Plot
#'
#' Creates a ribbon plot showing central credible intervals (CI) and outer intervals of the posterior predictive distribution.
#' This plot helps visualize how well the model fits the observed data across a range of probability intervals.
#' It can be used to assess the model's ability to capture the central tendency and variation in the data.
#'
#' @param x A `bipod` object containing a 'counts' field with the observed data. This is typically a time series or count data.
#' @param x_fit A fitted model object containing posterior draws of the replicated data (`yrep`). This model should be fitted using Bayesian methods.
#' @param n An integer specifying the number of posterior predictive draws to use in the plot. (default is 500). A higher value results in a more stable representation of the posterior distribution.
#' @param prob A numeric value between 0 and 1 representing the probability mass for the inner interval (e.g., the 50% CI). The default is 0.5.
#' @param prob_outer A numeric value between 0 and 1 representing the probability mass for the outer interval (e.g., the 90% CI). The default is 0.9.
#'
#' @return A ribbon plot created using `bayesplot::ppc_ribbon`. The plot shows the posterior predictive distribution with shaded regions representing credible intervals at different probability levels.
#' @export
biPOD_ppc_ribbon <- function(x, x_fit, n=500, prob = 0.5, prob_outer = 0.9) {
  if (!(inherits(x, "bipod"))) stop("Input must be a bipod object")
  y <- x$counts$count
  if (x$metadata$sampling == 'variational') {
    yrep <- x_fit$draws[,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  } else {
    yrep <- x_fit$draws[,,which(grepl("yrep", x_fit$parameters))] %>% posterior::as_draws_matrix()
  }
  yrep <- yrep[1:n,]
  x <- x$counts$time
  if ("q" %in% x_fit$parameters) y = log(y)
  if (!(is.null(x_fit$normalization_pars))) y = (y - mean(y)) / stats::sd(y)
  bayesplot::ppc_ribbon(y, yrep, x=x, prob = prob, prob_outer = prob_outer)
}
