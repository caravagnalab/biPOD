#' Overlay Posterior Predictive Density
#'
#' Creates an overlay of the posterior predictive density over the observed data.
#' This function uses the posterior predictive draws to compare the observed data distribution with the model-predicted distribution.
#'
#' @param x A `bipod` object that contains a 'counts' field with the observed data.
#' @param x_fit A fitted model object that contains posterior draws of the replicated data (yrep).
#' @param n Number of posterior predictive draws to use in the plot. (default is 500)
#'
#' @return A density overlay plot using `bayesplot::ppc_dens_overlay`.
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
#' Creates a plot displaying the posterior predictive intervals for the observed data.
#' This helps visualize how well the model captures the variation in the data over time.
#'
#' @param x A `bipod` object that contains a 'counts' field with the observed data.
#' @param x_fit A fitted model object that contains posterior draws of the replicated data (yrep).
#' @param n Number of posterior predictive draws to use in the plot. (default is 500)
#'
#' @return A plot showing posterior predictive intervals using `bayesplot::ppc_intervals`.
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
#' Creates a ribbon plot showing the central credible intervals (CI) and outer intervals of the posterior predictive distribution.
#' The plot visualizes how well the model fits the data across a range of probability intervals.
#'
#' @param x A `bipod` object that contains a 'counts' field with the observed data.
#' @param x_fit A fitted model object that contains posterior draws of the replicated data (yrep).
#' @param n Number of posterior predictive draws to use in the plot. (default is 500)
#' @param prob Numeric value between 0 and 1 representing the probability mass for the inner interval. (default is 0.5)
#' @param prob_outer Numeric value between 0 and 1 representing the probability mass for the outer interval. (default is 0.9)
#'
#' @return A ribbon plot using `bayesplot::ppc_ribbon`.
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
